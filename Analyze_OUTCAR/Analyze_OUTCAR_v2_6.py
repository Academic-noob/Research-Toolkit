#!/usr/bin/env python3
"""
Analyze_OUTCAR v2.6
---------------

A practical offline snapshot audit tool for VASP structure optimizations.

Main improvements over the original script:
1. Reads actual VASP settings (EDIFFG, EDIFF, NELM, NSW, IBRION, ISIF).
2. Separates official VASP convergence from heuristic warnings.
3. Treats total drift as a warning, not a formal convergence criterion.
4. Exports CONTCAR preferentially; falls back to OUTCAR only when necessary.
5. Preserves lattice vectors and Selective Dynamics flags where possible.
6. Writes CSV and JSON reports for batch workflows.
7. Warns explicitly when OUTCAR and OSZICAR ionic-step counts disagree.
8. Uses a robust floating-point parser supporting scientific notation.
9. Estimates remaining ionic steps and approximate wall-clock time from a downloaded snapshot when defensible.
10. Ignores benign POTCAR "kinetic energy error for atom" metadata and prints context for real error hits.
11. Does not infer whether a scheduler job is running from local file timestamps.
12. Locks each ionic-step energy to the last OUTCAR ``free energy TOTEN`` value
    immediately preceding its completed ``POSITION TOTAL-FORCE`` block.
13. Uses OSZICAR ionic ``N F=`` summaries as an independent cross-check rather
    than silently replacing the OUTCAR energy convention.

Typical usage:
    python Analyze_OUTCAR_v2_6.py
    python Analyze_OUTCAR_v2_6.py OUTCAR OSZICAR POSCAR
    python Analyze_OUTCAR_v2_6.py --incar INCAR --output-dir vasp_doctor_report
    python Analyze_OUTCAR_v2_6.py --no-plots --no-export

The default output directory is:
    vasp_doctor_report/
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import sys
import statistics
from datetime import datetime
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import numpy as np

FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
FLOAT_RE = re.compile(rf"^{FLOAT}$")

# Fatal patterns change the final audit status to ERROR_KEYWORD_DETECTED.
# A bare word "ERROR" is intentionally *not* fatal: normal VASP OUTCAR files can
# contain benign POTCAR metadata such as "kinetic energy error for atom=...".
FATAL_ERROR_PATTERNS = {
    "BRMIX": re.compile(r"\bBRMIX\b", re.I),
    "ZBRENT": re.compile(r"\bZBRENT\b", re.I),
    "EDDDAV": re.compile(r"\bEDDDAV\b", re.I),
    "internal error": re.compile(r"internal\s+error", re.I),
    "segmentation fault": re.compile(r"segmentation\s+fault", re.I),
    "LAPACK": re.compile(r"\bLAPACK\b", re.I),
    "very serious problems": re.compile(r"very\s+serious\s+problems", re.I),
}

# These are shown for manual inspection but do not automatically mark the run as
# failed. This avoids false alarms while still exposing unusual messages.
SUSPICIOUS_PATTERNS = {
    "generic ERROR text": re.compile(r"\bERROR\b", re.I),
}

# Known benign phrases that should not be shown as failures or suspicious hits.
BENIGN_ERROR_PATTERNS = {
    "POTCAR kinetic-energy metadata": re.compile(
        r"kinetic\s+energy\s+error\s+for\s+atom\s*=", re.I
    ),
}

SETTING_KEYS_FLOAT = ("EDIFFG", "EDIFF")
SETTING_KEYS_INT = ("NELM", "NSW", "IBRION", "ISIF", "ISYM")


@dataclass
class ScalarEvent:
    line_no: int
    value: float
    vector: Optional[list[float]] = None


@dataclass
class LatticeEvent:
    line_no: int
    lattice: list[list[float]]


@dataclass
class PositionEvent:
    header_line: int
    end_line: int
    positions: list[list[float]]
    forces: list[list[float]]
    max_force: float
    max_atom_idx: int


@dataclass
class IonicStep:
    step: int
    header_line: int
    end_line: int
    energy: Optional[float]
    energy_source: Optional[str]
    delta_energy: Optional[float]
    max_force: float
    max_atom_idx: int
    drift: Optional[float]
    drift_vec: Optional[list[float]]
    pressure: Optional[float]
    scf_steps: Optional[int]
    lattice: Optional[list[list[float]]] = None
    positions: list[list[float]] = field(default_factory=list)
    forces: list[list[float]] = field(default_factory=list)


def is_float_token(token: str) -> bool:
    return bool(FLOAT_RE.match(token))


def safe_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def format_optional(value: Optional[float], fmt: str = ".6f") -> str:
    if value is None or not math.isfinite(value):
        return "N/A"
    return format(value, fmt)


def json_safe(value: Any) -> Any:
    """Convert dataclass/numpy/NaN-rich objects into standards-compliant JSON."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def parse_key_value_text(text: str, keys: Iterable[str]) -> dict[str, str]:
    """
    Extract VASP-style KEY = VALUE settings. It supports multiple assignments
    on a single line and ignores comments after # or !.
    """
    cleaned = text.split("#", 1)[0].split("!", 1)[0]
    result: dict[str, str] = {}
    for key in keys:
        match = re.search(rf"\b{re.escape(key)}\s*=\s*([^\s;]+)", cleaned, re.I)
        if match:
            result[key.upper()] = match.group(1)
    return result


def parse_incar(path: Path) -> dict[str, Any]:
    settings: dict[str, Any] = {}
    if not path.is_file():
        return settings

    all_keys = SETTING_KEYS_FLOAT + SETTING_KEYS_INT
    for line in path.read_text(errors="replace").splitlines():
        extracted = parse_key_value_text(line, all_keys)
        for key, raw in extracted.items():
            try:
                settings[key] = float(raw) if key in SETTING_KEYS_FLOAT else int(float(raw))
            except ValueError:
                settings[key] = raw
    return settings


@dataclass
class PoscarData:
    symbols: list[str]
    counts: list[int]
    lattice: list[list[float]]
    positions_cart: list[list[float]]
    selective_flags: Optional[list[list[str]]]
    coordinate_mode: str

    @property
    def nions(self) -> int:
        return sum(self.counts)

    @property
    def atom_labels(self) -> list[str]:
        labels: list[str] = []
        for symbol, count in zip(self.symbols, self.counts):
            labels.extend([symbol] * count)
        return labels


def parse_poscar(path: Path) -> Optional[PoscarData]:
    if not path.is_file():
        return None

    lines = path.read_text(errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError(f"{path}: POSCAR/CONTCAR 内容过短。")

    try:
        scale_raw = float(lines[1].split()[0])
        lattice_raw = np.array(
            [[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)],
            dtype=float,
        )
    except Exception as exc:
        raise ValueError(f"{path}: 无法解析缩放因子或晶格矢量。") from exc

    if scale_raw > 0:
        scale = scale_raw
    elif scale_raw < 0:
        raw_volume = abs(float(np.linalg.det(lattice_raw)))
        if raw_volume <= 0:
            raise ValueError(f"{path}: 晶格体积无效，无法应用负缩放因子。")
        scale = (abs(scale_raw) / raw_volume) ** (1.0 / 3.0)
    else:
        raise ValueError(f"{path}: POSCAR 第二行缩放因子不能为 0。")

    lattice = lattice_raw * scale

    line_idx = 5
    tokens = lines[line_idx].split()
    if not tokens:
        raise ValueError(f"{path}: 元素或原子数量行为空。")

    if all(token.lstrip("+-").isdigit() for token in tokens):
        counts = [int(token) for token in tokens]
        symbols = [f"X{i + 1}" for i in range(len(counts))]
        line_idx += 1
    else:
        symbols = tokens
        line_idx += 1
        count_tokens = lines[line_idx].split()
        if not count_tokens or not all(token.lstrip("+-").isdigit() for token in count_tokens):
            raise ValueError(f"{path}: 无法解析原子数量。")
        counts = [int(token) for token in count_tokens]
        line_idx += 1

    nions = sum(counts)
    if nions <= 0:
        raise ValueError(f"{path}: 原子总数必须大于 0。")

    selective = False
    if line_idx < len(lines) and lines[line_idx].strip().lower().startswith("s"):
        selective = True
        line_idx += 1

    if line_idx >= len(lines):
        raise ValueError(f"{path}: 缺少 Direct/Cartesian 坐标模式。")
    mode_line = lines[line_idx].strip().lower()
    line_idx += 1
    direct_mode = mode_line.startswith("d")
    cart_mode = mode_line.startswith("c") or mode_line.startswith("k")
    if not (direct_mode or cart_mode):
        raise ValueError(f"{path}: 无法识别坐标模式: {mode_line!r}")

    coords: list[list[float]] = []
    flags: list[list[str]] = []
    for atom_i in range(nions):
        if line_idx + atom_i >= len(lines):
            raise ValueError(f"{path}: 坐标行不足，预期 {nions} 个原子。")
        parts = lines[line_idx + atom_i].split()
        if len(parts) < 3:
            raise ValueError(f"{path}: 第 {atom_i + 1} 个原子坐标不足三列。")
        try:
            coords.append([float(parts[0]), float(parts[1]), float(parts[2])])
        except ValueError as exc:
            raise ValueError(f"{path}: 第 {atom_i + 1} 个原子坐标无法解析。") from exc
        if selective:
            atom_flags = parts[3:6] if len(parts) >= 6 else ["T", "T", "T"]
            flags.append([flag.upper()[0] for flag in atom_flags])

    coord_array = np.array(coords, dtype=float)
    if direct_mode:
        positions_cart = coord_array @ lattice
        coordinate_mode = "Direct"
    else:
        positions_cart = coord_array * scale
        coordinate_mode = "Cartesian"

    return PoscarData(
        symbols=symbols,
        counts=counts,
        lattice=lattice.tolist(),
        positions_cart=positions_cart.tolist(),
        selective_flags=flags if selective else None,
        coordinate_mode=coordinate_mode,
    )


class VaspDoctor:
    def __init__(
        self,
        outcar_path: Path,
        oszicar_path: Path,
        poscar_path: Path,
        incar_path: Path,
        contcar_path: Path,
        output_dir: Path,
        fallback_force_threshold: float = 0.02,
        drift_warning_threshold: float = 0.005,
        prefer_contcar: bool = True,
        trend_window: int = 8,
    ) -> None:
        self.outcar_path = outcar_path
        self.oszicar_path = oszicar_path
        self.poscar_path = poscar_path
        self.incar_path = incar_path
        self.contcar_path = contcar_path
        self.output_dir = output_dir
        self.fallback_force_threshold = fallback_force_threshold
        self.drift_warning_threshold = drift_warning_threshold
        self.prefer_contcar = prefer_contcar
        self.trend_window = max(3, int(trend_window))

        self.poscar: Optional[PoscarData] = None
        self.settings: dict[str, Any] = {}
        self.nions: Optional[int] = None

        self.position_events: list[PositionEvent] = []
        self.energy_events: list[ScalarEvent] = []
        self.drift_events: list[ScalarEvent] = []
        self.pressure_events: list[ScalarEvent] = []
        self.lattice_events: list[LatticeEvent] = []
        self.loop_plus_real_times: list[float] = []

        self.steps: list[IonicStep] = []
        self.nelm_list: list[int] = []
        self.oszicar_ionic_step_numbers: list[int] = []
        self.oszicar_ionic_energies: list[float] = []
        self.oszicar_incomplete_scf_steps: int = 0
        self.reliably_paired_energy_steps: int = 0
        self.force_only_analysis_steps: int = 0
        self.outcar_energy_fallback_steps: int = 0
        self.last_energy_crosscheck: Optional[dict[str, Any]] = None
        self.energy_crosscheck_mismatches: list[dict[str, Any]] = []

        self.official_converged = False
        self.error_hits: dict[str, list[int]] = {}
        self.error_context: dict[str, list[dict[str, Any]]] = {}
        self.suspicious_hits: dict[str, list[int]] = {}
        self.suspicious_context: dict[str, list[dict[str, Any]]] = {}
        self.benign_error_hits: dict[str, list[int]] = {}
        self.info_messages: list[str] = []
        self.warnings: list[str] = []
        self.status = "UNKNOWN"

    def read_inputs(self) -> None:
        if not self.outcar_path.is_file():
            raise FileNotFoundError(f"找不到 OUTCAR: {self.outcar_path}")

        self.settings.update(parse_incar(self.incar_path))
        if self.poscar_path.is_file():
            self.poscar = parse_poscar(self.poscar_path)
            self.nions = self.poscar.nions

    @staticmethod
    def _compact_line(line: str, limit: int = 260) -> str:
        text = line.rstrip("\n\r")
        return text if len(text) <= limit else text[: limit - 3] + "..."

    def _record_error_line(self, line_no: int, line: str) -> None:
        compact = self._compact_line(line)
        for name, pattern in BENIGN_ERROR_PATTERNS.items():
            if pattern.search(line):
                self.benign_error_hits.setdefault(name, []).append(line_no)
                return

        for name, pattern in FATAL_ERROR_PATTERNS.items():
            if pattern.search(line):
                self.error_hits.setdefault(name, []).append(line_no)
                self.error_context.setdefault(name, []).append(
                    {"line": line_no, "text": compact}
                )

        for name, pattern in SUSPICIOUS_PATTERNS.items():
            if pattern.search(line):
                self.suspicious_hits.setdefault(name, []).append(line_no)
                self.suspicious_context.setdefault(name, []).append(
                    {"line": line_no, "text": compact}
                )

    def parse_outcar(self) -> None:
        re_nions = re.compile(r"\bNIONS\s*=\s*(\d+)", re.I)
        re_energy = re.compile(rf"free\s+energy\s+TOTEN\s*=\s*({FLOAT})\s+eV", re.I)
        re_drift = re.compile(rf"total\s+drift\s*:\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})", re.I)
        re_pressure = re.compile(rf"external\s+pressure\s*=\s*({FLOAT})\s+kB", re.I)
        re_loop_plus = re.compile(
            rf"LOOP\+\s*:\s*cpu\s+time\s+{FLOAT}\s*:\s*real\s+time\s+({FLOAT})",
            re.I,
        )
        re_pos_head = re.compile(r"^\s*POSITION\s+TOTAL-FORCE", re.I)
        re_lattice_head = re.compile(r"^\s*direct\s+lattice\s+vectors", re.I)
        all_setting_keys = SETTING_KEYS_FLOAT + SETTING_KEYS_INT

        reading_positions = False
        position_header_line = 0
        positions: list[list[float]] = []
        forces: list[list[float]] = []

        reading_lattice = False
        lattice_header_line = 0
        lattice_rows: list[list[float]] = []

        def finalize_position(end_line: int) -> None:
            nonlocal positions, forces, reading_positions
            if not positions:
                reading_positions = False
                return
            force_norms = np.linalg.norm(np.array(forces, dtype=float), axis=1)
            max_idx = int(np.argmax(force_norms))
            self.position_events.append(
                PositionEvent(
                    header_line=position_header_line,
                    end_line=end_line,
                    positions=positions,
                    forces=forces,
                    max_force=float(force_norms[max_idx]),
                    max_atom_idx=max_idx,
                )
            )
            positions = []
            forces = []
            reading_positions = False

        with self.outcar_path.open("r", errors="replace") as handle:
            for line_no, line in enumerate(handle, start=1):
                self._record_error_line(line_no, line)

                if "reached required accuracy" in line.lower():
                    self.official_converged = True

                if self.nions is None:
                    nions_match = re_nions.search(line)
                    if nions_match:
                        self.nions = int(nions_match.group(1))

                extracted = parse_key_value_text(line, all_setting_keys)
                for key, raw in extracted.items():
                    try:
                        self.settings[key] = float(raw) if key in SETTING_KEYS_FLOAT else int(float(raw))
                    except ValueError:
                        self.settings.setdefault(key, raw)

                if reading_lattice:
                    parts = line.split()
                    if len(parts) >= 3 and all(is_float_token(token) for token in parts[:3]):
                        lattice_rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                        if len(lattice_rows) == 3:
                            self.lattice_events.append(
                                LatticeEvent(line_no=lattice_header_line, lattice=lattice_rows)
                            )
                            lattice_rows = []
                            reading_lattice = False
                        continue
                    lattice_rows = []
                    reading_lattice = False

                if re_lattice_head.search(line):
                    reading_lattice = True
                    lattice_header_line = line_no
                    lattice_rows = []
                    continue

                if reading_positions:
                    parts = line.split()
                    if len(parts) >= 6 and all(is_float_token(token) for token in parts[:6]):
                        positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
                        forces.append([float(parts[3]), float(parts[4]), float(parts[5])])
                        if self.nions is not None and len(positions) >= self.nions:
                            finalize_position(line_no)
                        continue

                    if positions:
                        finalize_position(line_no - 1)
                    # Continue processing the current line after finalization.

                if re_pos_head.search(line):
                    reading_positions = True
                    position_header_line = line_no
                    positions = []
                    forces = []
                    continue

                match = re_energy.search(line)
                if match:
                    self.energy_events.append(ScalarEvent(line_no=line_no, value=float(match.group(1))))

                match = re_drift.search(line)
                if match:
                    vector = [float(match.group(i)) for i in range(1, 4)]
                    self.drift_events.append(
                        ScalarEvent(line_no=line_no, value=float(np.linalg.norm(vector)), vector=vector)
                    )

                match = re_pressure.search(line)
                if match:
                    self.pressure_events.append(
                        ScalarEvent(line_no=line_no, value=float(match.group(1)))
                    )

                match = re_loop_plus.search(line)
                if match:
                    value = safe_float(match.group(1))
                    if value is not None and value > 0:
                        self.loop_plus_real_times.append(value)

        if reading_positions and positions:
            finalize_position(position_header_line + len(positions))

        if self.nions is None and self.position_events:
            self.nions = len(self.position_events[0].positions)

    def parse_oszicar(self) -> None:
        if not self.oszicar_path.is_file():
            self.warnings.append(f"找不到 OSZICAR: {self.oszicar_path}；不会输出电子步统计。")
            return

        electronic_re = re.compile(r"^\s*(?:DAV|RMM|CG|BROYDEN|DMP)\s*:", re.I)
        ionic_re = re.compile(rf"^\s*(\d+)\s+F\s*=\s*({FLOAT})", re.I)
        current_nelm = 0

        with self.oszicar_path.open("r", errors="replace") as handle:
            for line in handle:
                if electronic_re.search(line):
                    current_nelm += 1
                else:
                    ionic_match = ionic_re.search(line)
                    if ionic_match:
                        self.nelm_list.append(current_nelm if current_nelm > 0 else 1)
                        self.oszicar_ionic_step_numbers.append(int(ionic_match.group(1)))
                        self.oszicar_ionic_energies.append(float(ionic_match.group(2)))
                        current_nelm = 0

        self.oszicar_incomplete_scf_steps = current_nelm
        if current_nelm:
            self.warnings.append(
                f"OSZICAR 末尾存在 {current_nelm} 个尚未归入已完成离子步的电子步：下载快照可能截取于该 SCF 中，也可能对应一次中途停止的计算。"
            )
        if len(self.oszicar_ionic_energies) != len(self.nelm_list):
            self.warnings.append(
                "OSZICAR 离子步能量数量与电子步统计数量不一致；离子步能量将仅在可可靠解析时使用。"
            )

    @staticmethod
    def _events_in_window(
        events: list[ScalarEvent],
        lower: int,
        upper: float,
    ) -> list[ScalarEvent]:
        return [event for event in events if lower < event.line_no < upper]

    @staticmethod
    def _lattice_before(events: list[LatticeEvent], line_no: int) -> Optional[list[list[float]]]:
        candidates = [event for event in events if event.line_no < line_no]
        return candidates[-1].lattice if candidates else None

    @staticmethod
    def _choose_after_or_before(
        events: list[ScalarEvent],
        start: int,
        end: int,
        lower: int,
        upper: float,
        prefer: str = "after",
    ) -> Optional[ScalarEvent]:
        candidates = [event for event in events if lower < event.line_no < upper]
        after = [event for event in candidates if event.line_no > end]
        before = [event for event in candidates if event.line_no < start]

        if prefer == "before":
            if before:
                return before[-1]
            if after:
                return after[0]
        else:
            if after:
                return after[0]
            if before:
                return before[-1]

        return min(candidates, key=lambda event: abs(event.line_no - start)) if candidates else None

    def build_steps(self) -> None:
        """Build ionic-step records by anchoring OUTCAR TOTEN values to force blocks.

        For a completed ionic step, VASP writes one or more ``free energy TOTEN``
        records during electronic relaxation and then a ``POSITION TOTAL-FORCE``
        block. The relevant energy for that force block is the *last* TOTEN found
        after the previous completed force block and before the current force
        block. This is the robust form of the state-machine logic used by the
        original script.

        OSZICAR ``N F=`` values are retained as an independent cross-check. They
        are not silently substituted for the OUTCAR-locked energy. If no locked
        OUTCAR energy is available for a force block, OSZICAR is used only as an
        explicit fallback and a warning is emitted.
        """
        self.steps = []
        self.reliably_paired_energy_steps = 0
        self.force_only_analysis_steps = 0
        self.outcar_energy_fallback_steps = 0
        self.last_energy_crosscheck = None
        self.energy_crosscheck_mismatches = []

        if not self.position_events:
            return

        outcar_force_count = len(self.position_events)
        oszicar_count = len(self.oszicar_ionic_energies)
        oszicar_available = self.oszicar_path.is_file()

        for idx, position in enumerate(self.position_events):
            prev_end = self.position_events[idx - 1].end_line if idx > 0 else 0
            next_header: float = (
                self.position_events[idx + 1].header_line
                if idx + 1 < outcar_force_count
                else float("inf")
            )

            # Primary rule: last TOTEN before the completed force block.
            before_candidates = [
                event
                for event in self.energy_events
                if prev_end < event.line_no < position.header_line
            ]
            energy_event = before_candidates[-1] if before_candidates else None
            energy_source: Optional[str] = None

            # Rare fallback for non-standard output ordering.
            if energy_event is None:
                after_candidates = [
                    event
                    for event in self.energy_events
                    if position.end_line < event.line_no < next_header
                ]
                if after_candidates:
                    energy_event = after_candidates[0]
                    energy_source = (
                        "OUTCAR TOTEN fallback: first value after completed force block"
                    )
                    self.outcar_energy_fallback_steps += 1
                    self.warnings.append(
                        f"离子步 {idx + 1}: 未找到力块之前的 TOTEN，"
                        "回退到力块之后、下一个力块之前的首个 TOTEN。请人工复核。"
                    )

            if energy_event is not None:
                ionic_energy = energy_event.value
                if energy_source is None:
                    energy_source = (
                        "OUTCAR TOTEN locked: last value before completed POSITION TOTAL-FORCE block"
                    )
                self.reliably_paired_energy_steps += 1
            elif idx < oszicar_count:
                ionic_energy = self.oszicar_ionic_energies[idx]
                energy_source = "OSZICAR ionic F= fallback: locked OUTCAR TOTEN unavailable"
                self.outcar_energy_fallback_steps += 1
                self.warnings.append(
                    f"离子步 {idx + 1}: 无法从 OUTCAR 锁定 TOTEN，"
                    "临时回退到 OSZICAR F=。"
                )
            else:
                ionic_energy = None
                energy_source = "N/A: no locked OUTCAR TOTEN and no OSZICAR fallback"
                self.force_only_analysis_steps += 1

            step_number = (
                self.oszicar_ionic_step_numbers[idx]
                if idx < len(self.oszicar_ionic_step_numbers)
                else idx + 1
            )

            # Independent OSZICAR audit. Do not replace OUTCAR values silently.
            if energy_event is not None and idx < oszicar_count:
                oszicar_energy = self.oszicar_ionic_energies[idx]
                delta = energy_event.value - oszicar_energy
                crosscheck = {
                    "step": step_number,
                    "outcar_locked_toten_eV": energy_event.value,
                    "outcar_toten_line": energy_event.line_no,
                    "oszicar_F_eV": oszicar_energy,
                    "outcar_minus_oszicar_eV": delta,
                }
                self.last_energy_crosscheck = crosscheck
                if abs(delta) > 1.0e-4:
                    self.energy_crosscheck_mismatches.append(crosscheck)

            drift_event = self._choose_after_or_before(
                self.drift_events,
                start=position.header_line,
                end=position.end_line,
                lower=prev_end,
                upper=next_header,
                prefer="after",
            )
            pressure_event = self._choose_after_or_before(
                self.pressure_events,
                start=position.header_line,
                end=position.end_line,
                lower=prev_end,
                upper=next_header,
                prefer="before",
            )

            scf_steps = self.nelm_list[idx] if idx < len(self.nelm_list) else None
            self.steps.append(
                IonicStep(
                    step=step_number,
                    header_line=position.header_line,
                    end_line=position.end_line,
                    energy=ionic_energy,
                    energy_source=energy_source,
                    delta_energy=None,
                    max_force=position.max_force,
                    max_atom_idx=position.max_atom_idx,
                    drift=drift_event.value if drift_event else None,
                    drift_vec=drift_event.vector if drift_event else None,
                    pressure=pressure_event.value if pressure_event else None,
                    scf_steps=scf_steps,
                    lattice=self._lattice_before(self.lattice_events, position.header_line),
                    positions=position.positions,
                    forces=position.forces,
                )
            )

        previous_energy: Optional[float] = None
        for step in self.steps:
            if step.energy is not None and previous_energy is not None:
                step.delta_energy = abs(step.energy - previous_energy)
            if step.energy is not None:
                previous_energy = step.energy

        if oszicar_available and outcar_force_count != oszicar_count:
            self.warnings.append(
                "OUTCAR 力块数量与 OSZICAR 已完成离子步数量不一致："
                f"OUTCAR force blocks={outcar_force_count}，OSZICAR ionic F= rows={oszicar_count}。"
                "OUTCAR 锁定能量仍按已完成力块分析；OSZICAR 仅用于交叉核对。"
            )

        if self.energy_crosscheck_mismatches:
            last = self.energy_crosscheck_mismatches[-1]
            self.warnings.append(
                "OUTCAR 锁定 TOTEN 与 OSZICAR F= 存在差异。"
                f"最近一次差异位于 step {last['step']}: "
                f"OUTCAR={last['outcar_locked_toten_eV']:.8f} eV, "
                f"OSZICAR={last['oszicar_F_eV']:.8f} eV, "
                f"OUTCAR-OSZICAR={last['outcar_minus_oszicar_eV']:.8f} eV。"
                "脚本不会自动替换任一能量，请结合原始文件判断原因。"
            )

        if not self.energy_events:
            self.warnings.append(
                "OUTCAR 中未解析到 free energy TOTEN；离子步能量可能需要 OSZICAR 回退。"
            )


    @staticmethod
    def format_duration(seconds: Optional[float]) -> str:
        if seconds is None or not math.isfinite(seconds) or seconds < 0:
            return "N/A"
        total = int(round(seconds))
        days, rem = divmod(total, 86400)
        hours, rem = divmod(rem, 3600)
        minutes, secs = divmod(rem, 60)
        chunks: list[str] = []
        if days:
            chunks.append(f"{days} d")
        if hours or days:
            chunks.append(f"{hours} h")
        if minutes or hours or days:
            chunks.append(f"{minutes} min")
        chunks.append(f"{secs} s")
        return " ".join(chunks)

    def estimate_progress(self) -> dict[str, Any]:
        """
        Estimate remaining ionic steps from the recent logarithmic decline of the active
        convergence metric. This is deliberately conservative: noisy, flat, or rising
        trends are reported as not estimable instead of producing a fake ETA.
        """
        ediffg = safe_float(self.settings.get("EDIFFG"))
        nsw = self.settings.get("NSW")
        completed_steps = len(self.steps)
        remaining_nsw_budget = (
            max(0, nsw - completed_steps)
            if isinstance(nsw, int) and nsw > 0
            else None
        )

        base: dict[str, Any] = {
            "available": False,
            "metric": None,
            "metric_unit": None,
            "target": None,
            "current": None,
            "current_to_target_ratio": None,
            "trend_window_requested": self.trend_window,
            "trend_points_used": 0,
            "estimated_remaining_ionic_steps": None,
            "estimated_completion_ionic_step": None,
            "median_loop_plus_walltime_seconds": None,
            "estimated_remaining_walltime_seconds": None,
            "estimated_remaining_walltime_human": "N/A",
            "confidence": "UNAVAILABLE",
            "reason": None,
            "remaining_nsw_budget": remaining_nsw_budget,
            "estimated_to_fit_within_current_nsw": None,
            "warning": "This is a trend-based estimate, not a guaranteed completion time.",
        }

        if not self.steps:
            base["reason"] = "No completed ionic steps are available yet."
            return base

        ibrion = self.settings.get("IBRION")
        if ibrion == -1 or nsw == 0:
            base["reason"] = "Static calculation: no ionic-optimization ETA is applicable."
            return base

        if self.official_converged:
            base.update(
                {
                    "available": True,
                    "estimated_remaining_ionic_steps": 0,
                    "estimated_completion_ionic_step": completed_steps,
                    "estimated_remaining_walltime_seconds": 0.0,
                    "estimated_remaining_walltime_human": "0 s",
                    "confidence": "OFFICIAL_CONVERGENCE_CONFIRMED",
                    "reason": "VASP official convergence stop message was detected.",
                    "estimated_to_fit_within_current_nsw": True,
                }
            )
            return base

        if ediffg is not None and ediffg < 0:
            metric = "max_force"
            unit = "eV/Å"
            target = abs(ediffg)
            series = [(step.step, step.max_force) for step in self.steps]
        elif ediffg is not None and ediffg > 0:
            metric = "delta_energy"
            unit = "eV"
            target = ediffg
            series = [
                (step.step, step.delta_energy)
                for step in self.steps
                if step.delta_energy is not None
            ]
        else:
            metric = "max_force_heuristic"
            unit = "eV/Å"
            target = self.fallback_force_threshold
            series = [(step.step, step.max_force) for step in self.steps]
            base["warning"] = (
                "EDIFFG is unavailable. ETA uses the fallback force threshold and is heuristic only."
            )

        valid = [
            (int(step), float(value))
            for step, value in series
            if value is not None and math.isfinite(float(value)) and float(value) > 0
        ]
        base.update({"metric": metric, "metric_unit": unit, "target": target})
        if not valid:
            base["reason"] = "The active convergence metric could not be parsed."
            return base

        current_step, current = valid[-1]
        base["current"] = current
        base["current_to_target_ratio"] = current / target if target > 0 else None

        loop_times = [value for value in self.loop_plus_real_times if math.isfinite(value) and value > 0]
        if loop_times:
            recent_loop_times = loop_times[-self.trend_window :]
            base["median_loop_plus_walltime_seconds"] = statistics.median(recent_loop_times)

        if current <= target:
            base.update(
                {
                    "available": True,
                    "estimated_remaining_ionic_steps": 0,
                    "estimated_completion_ionic_step": current_step,
                    "estimated_remaining_walltime_seconds": 0.0,
                    "estimated_remaining_walltime_human": "0 s",
                    "confidence": "THRESHOLD_ALREADY_MET",
                    "reason": "The parsed convergence metric already satisfies the target.",
                    "estimated_to_fit_within_current_nsw": True,
                }
            )
            return base

        recent = valid[-self.trend_window :]
        base["trend_points_used"] = len(recent)
        if len(recent) < 3:
            base["reason"] = "At least 3 completed metric points are required for a trend estimate."
            return base

        x = np.array([step for step, _ in recent], dtype=float)
        y = np.log(np.array([value for _, value in recent], dtype=float))
        slope, intercept = np.polyfit(x, y, 1)
        predicted = slope * x + intercept
        ss_res = float(np.sum((y - predicted) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
        consecutive_changes = np.diff(np.array([value for _, value in recent], dtype=float))
        non_increasing_fraction = (
            float(np.mean(consecutive_changes <= 0)) if len(consecutive_changes) else None
        )

        base.update(
            {
                "trend_log_slope_per_ionic_step": float(slope),
                "trend_multiplicative_factor_per_ionic_step": float(math.exp(slope)),
                "trend_r_squared": r_squared,
                "recent_non_increasing_fraction": non_increasing_fraction,
            }
        )

        if not math.isfinite(slope) or slope >= -1e-10:
            base["reason"] = "Recent convergence metric is flat or rising; a defensible ETA cannot be calculated."
            return base

        estimated_steps_float = math.log(target / current) / slope
        if not math.isfinite(estimated_steps_float) or estimated_steps_float < 0:
            base["reason"] = "Trend extrapolation produced an invalid remaining-step estimate."
            return base

        estimated_steps = int(math.ceil(estimated_steps_float))
        if estimated_steps > 100000:
            base["reason"] = "Trend extrapolation exceeds 100000 ionic steps and is not useful."
            return base

        if len(recent) >= 6 and r_squared >= 0.80 and (non_increasing_fraction or 0.0) >= 0.60:
            confidence = "HIGH"
        elif len(recent) >= 4 and r_squared >= 0.45:
            confidence = "MEDIUM"
        else:
            confidence = "LOW"

        median_step_seconds = base["median_loop_plus_walltime_seconds"]
        eta_seconds = estimated_steps * median_step_seconds if median_step_seconds is not None else None
        fits_nsw = (
            estimated_steps <= remaining_nsw_budget
            if remaining_nsw_budget is not None
            else None
        )
        base.update(
            {
                "available": True,
                "estimated_remaining_ionic_steps": estimated_steps,
                "estimated_completion_ionic_step": current_step + estimated_steps,
                "estimated_remaining_walltime_seconds": eta_seconds,
                "estimated_remaining_walltime_human": self.format_duration(eta_seconds),
                "confidence": confidence,
                "reason": "Estimated from the recent logarithmic decline of the active convergence metric.",
                "estimated_to_fit_within_current_nsw": fits_nsw,
            }
        )
        return base

    def assess(self) -> dict[str, Any]:
        if not self.steps:
            self.status = (
                "ERROR_KEYWORD_DETECTED"
                if self.error_hits
                else "NO_COMPLETED_IONIC_STEPS_IN_SNAPSHOT"
            )
            return {
                "status": self.status,
                "criterion_met_by_parsed_values": None,
                "convergence_basis": None,
                "official_vasp_stop_message": self.official_converged,
                "formal_convergence_confirmed": False,
                "heuristic_force_ok": None,
                "progress_estimate": self.estimate_progress(),
            }

        last = self.steps[-1]
        ediffg = safe_float(self.settings.get("EDIFFG"))
        nelm = self.settings.get("NELM")
        nsw = self.settings.get("NSW")
        ibrion = self.settings.get("IBRION")

        formal_ok: Optional[bool] = None
        convergence_basis = "VASP official stop message"

        if self.official_converged:
            formal_ok = True
        elif ediffg is not None and ediffg < 0:
            formal_ok = last.max_force <= abs(ediffg)
            convergence_basis = f"max force <= |EDIFFG| ({abs(ediffg):.6g} eV/Å)"
        elif ediffg is not None and ediffg > 0:
            formal_ok = (
                last.delta_energy is not None and last.delta_energy <= ediffg
            )
            convergence_basis = f"|ΔE| <= EDIFFG ({ediffg:.6g} eV)"
        elif ibrion == -1 or nsw == 0:
            formal_ok = None
            convergence_basis = "static calculation; no ionic optimization criterion"
        else:
            convergence_basis = "EDIFFG unavailable; only heuristic force threshold is available"

        heuristic_force_ok = last.max_force <= self.fallback_force_threshold

        saturated_steps: list[int] = []
        if isinstance(nelm, int) and nelm > 0:
            saturated_steps = [
                step.step
                for step in self.steps
                if step.scf_steps is not None and step.scf_steps >= nelm
            ]

        nsw_reached = isinstance(nsw, int) and nsw > 0 and len(self.steps) >= nsw
        progress = self.estimate_progress()

        if self.error_hits:
            self.status = "ERROR_KEYWORD_DETECTED"
        elif ibrion == -1 or nsw == 0:
            self.status = "STATIC_CALCULATION"
        elif self.official_converged:
            self.status = "CONVERGED_OFFICIAL"
        elif formal_ok is True:
            self.status = "CRITERION_MET_NOT_CONFIRMED"
        elif nsw_reached:
            self.status = "NSW_REACHED_NOT_CONVERGED"
        elif formal_ok is False:
            self.status = "NOT_CONVERGED_IN_SNAPSHOT"
        else:
            self.status = "INCOMPLETE_OR_UNKNOWN_SNAPSHOT"

        if last.drift is not None and last.drift > self.drift_warning_threshold:
            self.warnings.append(
                f"最后一步 total drift={last.drift:.6g} eV/Å，"
                f"高于辅助警戒值 {self.drift_warning_threshold:.6g} eV/Å。"
                "该指标仅用于数值质量检查，不是 VASP 正式离子收敛条件。"
            )

        if saturated_steps:
            preview = ", ".join(str(step) for step in saturated_steps[:10])
            suffix = " ..." if len(saturated_steps) > 10 else ""
            self.warnings.append(
                f"{len(saturated_steps)} 个离子步达到或超过 NELM={nelm}：{preview}{suffix}。"
                "电子自洽可能不稳定。"
            )

        if nsw_reached and not self.official_converged:
            self.warnings.append(
                f"离子步数已达到 NSW={nsw}，但未检测到 VASP 正式收敛停止信息。"
            )

        return {
            "status": self.status,
            "criterion_met_by_parsed_values": formal_ok,
            "convergence_basis": convergence_basis,
            "official_vasp_stop_message": self.official_converged,
            "formal_convergence_confirmed": self.official_converged,
            "heuristic_force_threshold_eV_per_A": self.fallback_force_threshold,
            "heuristic_force_ok": heuristic_force_ok,
            "nelm_saturated_steps": saturated_steps,
            "nsw_reached": nsw_reached,
            "progress_estimate": progress,
        }

    def atom_label(self, atom_idx: int) -> str:
        if self.poscar and atom_idx < len(self.poscar.atom_labels):
            return self.poscar.atom_labels[atom_idx]
        return "?"

    def bottleneck_atoms(self, top_n: int = 5) -> list[dict[str, Any]]:
        if not self.steps:
            return []

        counter = Counter(step.max_atom_idx for step in self.steps)
        output: list[dict[str, Any]] = []
        for atom_idx, count in counter.most_common(top_n):
            selected = [step for step in self.steps if step.max_atom_idx == atom_idx]
            max_step = max(selected, key=lambda step: step.max_force)
            force_vec = max_step.forces[atom_idx]
            output.append(
                {
                    "atom_index_1based": atom_idx + 1,
                    "element": self.atom_label(atom_idx),
                    "bottleneck_count": count,
                    "fraction": count / len(self.steps),
                    "historical_max_force_eV_per_A": max_step.max_force,
                    "historical_max_force_step": max_step.step,
                    "force_vector_at_historical_max_eV_per_A": force_vec,
                }
            )
        return output

    def export_latest_structure(self, output_path: Path) -> str:
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if self.prefer_contcar and self.contcar_path.is_file() and self.contcar_path.stat().st_size > 0:
            try:
                parse_poscar(self.contcar_path)
            except Exception as exc:
                self.warnings.append(
                    f"CONTCAR 存在但无法解析，将回退到 OUTCAR 导出：{exc}"
                )
            else:
                shutil.copy2(self.contcar_path, output_path)
                return f"copied from {self.contcar_path}"

        if not self.steps:
            raise RuntimeError("OUTCAR 中没有可导出的离子步。")

        last = self.steps[-1]
        lattice = last.lattice
        if lattice is None and self.poscar is not None:
            lattice = self.poscar.lattice
            self.warnings.append(
                "OUTCAR 中未解析到最终晶格，导出结构使用 POSCAR 初始晶格。"
                "该回退仅适合固定晶格任务；若 ISIF 允许晶格变化，请不要直接使用该文件。"
            )

        if lattice is None:
            raise RuntimeError(
                "无法导出结构：既没有可用 CONTCAR，也无法从 OUTCAR/POSCAR 获得晶格矢量。"
            )

        nions = len(last.positions)
        if self.poscar and self.poscar.nions == nions:
            symbols = self.poscar.symbols
            counts = self.poscar.counts
            flags = self.poscar.selective_flags
        else:
            symbols = ["X"]
            counts = [nions]
            flags = None
            self.warnings.append(
                "POSCAR 元素信息不可用或原子总数不一致；导出文件将使用占位元素 X。"
            )

        with output_path.open("w") as handle:
            handle.write("Generated by VaspDoctor v2.6 from OUTCAR\n")
            handle.write("1.0\n")
            for vector in lattice:
                handle.write(" ".join(f"{value:18.10f}" for value in vector) + "\n")
            handle.write(" " + " ".join(symbols) + "\n")
            handle.write(" " + " ".join(str(count) for count in counts) + "\n")
            if flags and len(flags) == nions:
                handle.write("Selective dynamics\n")
            handle.write("Cartesian\n")
            for atom_idx, coords in enumerate(last.positions):
                line = " ".join(f"{value:18.10f}" for value in coords)
                if flags and len(flags) == nions:
                    line += "   " + " ".join(flags[atom_idx])
                handle.write(line + "\n")

        return "rebuilt from OUTCAR"

    def write_csv(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = [
            "step",
            "energy_eV",
            "energy_source",
            "delta_energy_eV",
            "max_force_eV_per_A",
            "max_atom_index_1based",
            "max_atom_element",
            "drift_eV_per_A",
            "drift_x",
            "drift_y",
            "drift_z",
            "pressure_kB",
            "scf_steps",
            "outcar_header_line",
        ]
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for step in self.steps:
                drift_vec = step.drift_vec or [None, None, None]
                writer.writerow(
                    {
                        "step": step.step,
                        "energy_eV": step.energy,
                        "energy_source": step.energy_source,
                        "delta_energy_eV": step.delta_energy,
                        "max_force_eV_per_A": step.max_force,
                        "max_atom_index_1based": step.max_atom_idx + 1,
                        "max_atom_element": self.atom_label(step.max_atom_idx),
                        "drift_eV_per_A": step.drift,
                        "drift_x": drift_vec[0],
                        "drift_y": drift_vec[1],
                        "drift_z": drift_vec[2],
                        "pressure_kB": step.pressure,
                        "scf_steps": step.scf_steps,
                        "outcar_header_line": step.header_line,
                    }
                )

    def write_json(self, path: Path, assessment: dict[str, Any], export_source: Optional[str]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        last = self.steps[-1] if self.steps else None
        last_step_summary = None
        if last is not None:
            max_force_vector = last.forces[last.max_atom_idx] if last.forces else None
            last_step_summary = {
                "step": last.step,
                "header_line": last.header_line,
                "end_line": last.end_line,
                "energy_eV": last.energy,
                "energy_source": last.energy_source,
                "delta_energy_eV": last.delta_energy,
                "max_force_eV_per_A": last.max_force,
                "max_atom_index_1based": last.max_atom_idx + 1,
                "max_atom_element": self.atom_label(last.max_atom_idx),
                "max_force_vector_eV_per_A": max_force_vector,
                "drift_eV_per_A": last.drift,
                "drift_vector_eV_per_A": last.drift_vec,
                "pressure_kB": last.pressure,
                "scf_steps": last.scf_steps,
                "lattice_A": last.lattice,
            }

        report = {
            "tool": "VaspDoctor v2.6",
            "inputs": {
                "OUTCAR": str(self.outcar_path),
                "OSZICAR": str(self.oszicar_path),
                "POSCAR": str(self.poscar_path),
                "INCAR": str(self.incar_path),
                "CONTCAR": str(self.contcar_path),
            },
            "settings": self.settings,
            "nions": self.nions,
            "outcar_force_blocks": len(self.position_events),
            "oszicar_completed_ionic_steps": len(self.oszicar_ionic_energies),
            "analysis_records": len(self.steps),
            "reliably_paired_energy_steps": self.reliably_paired_energy_steps,
            "force_only_analysis_steps": self.force_only_analysis_steps,
            "oszicar_incomplete_scf_steps": self.oszicar_incomplete_scf_steps,
            "loop_plus_walltime_samples_seconds": self.loop_plus_real_times,
            "latest_raw_outcar_toten_eV_diagnostic_only": self.energy_events[-1].value if self.energy_events else None,
            "latest_raw_outcar_toten_line": self.energy_events[-1].line_no if self.energy_events else None,
            "latest_raw_outcar_toten_note": "Raw tail reference only. Ionic-step energies are separately locked to completed POSITION TOTAL-FORCE blocks.",
            "outcar_energy_fallback_steps": self.outcar_energy_fallback_steps,
            "last_energy_crosscheck": self.last_energy_crosscheck,
            "energy_crosscheck_mismatches": self.energy_crosscheck_mismatches,
            "assessment": assessment,
            "last_step": last_step_summary,
            "bottleneck_atoms": self.bottleneck_atoms(),
            "fatal_error_hits": self.error_hits,
            "fatal_error_context": self.error_context,
            "suspicious_hits": self.suspicious_hits,
            "suspicious_context": self.suspicious_context,
            "ignored_benign_error_hits": self.benign_error_hits,
            "info_messages": self.info_messages,
            "warnings": self.warnings,
            "export_source": export_source,
        }
        path.write_text(
            json.dumps(json_safe(report), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )


    def write_snapshot_status(self, path: Path, assessment: dict[str, Any]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        progress = assessment.get("progress_estimate", {})
        last = self.steps[-1] if self.steps else None
        lines = [
            "VaspDoctor v2.6 offline snapshot status",
            f"generated_at = {datetime.now().astimezone().isoformat(timespec='seconds')}",
            f"status = {assessment.get('status')}",
            f"outcar_force_blocks = {len(self.position_events)}",
            f"oszicar_completed_ionic_steps = {len(self.oszicar_ionic_energies)}",
            f"analysis_records = {len(self.steps)}",
            f"reliably_paired_energy_steps = {self.reliably_paired_energy_steps}",
            f"force_only_analysis_steps = {self.force_only_analysis_steps}",
            f"outcar_energy_fallback_steps = {self.outcar_energy_fallback_steps}",
            f"incomplete_scf_steps_at_snapshot_end = {self.oszicar_incomplete_scf_steps}",
        ]
        if last is not None:
            lines.extend(
                [
                    f"last_completed_ionic_energy_eV = {format_optional(last.energy, '.10g')}",
                    f"last_completed_ionic_energy_source = {last.energy_source or 'N/A'}",
                    f"last_max_force_eV_per_A = {last.max_force:.8g}",
                    f"last_delta_energy_eV = {format_optional(last.delta_energy, '.8g')}",
                    f"last_scf_steps = {last.scf_steps if last.scf_steps is not None else 'N/A'}",
                ]
            )
        lines.extend(
            [
                f"eta_available = {progress.get('available')}",
                f"eta_metric = {progress.get('metric')}",
                f"eta_target = {progress.get('target')}",
                f"eta_current = {progress.get('current')}",
                f"eta_remaining_ionic_steps = {progress.get('estimated_remaining_ionic_steps')}",
                f"eta_remaining_walltime = {progress.get('estimated_remaining_walltime_human')}",
                f"eta_confidence = {progress.get('confidence')}",
                f"eta_reason = {progress.get('reason')}",
            ]
        )
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def plot(self) -> None:
        if not self.steps:
            return

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            self.warnings.append("未安装 matplotlib；跳过绘图。")
            return

        self.output_dir.mkdir(parents=True, exist_ok=True)
        x = np.array([step.step for step in self.steps], dtype=int)

        def arr(attr: str) -> np.ndarray:
            values = []
            for step in self.steps:
                value = getattr(step, attr)
                values.append(np.nan if value is None else value)
            return np.array(values, dtype=float)

        energy = arr("energy")
        delta_energy = arr("delta_energy")
        force = arr("max_force")
        drift = arr("drift")
        pressure = arr("pressure")
        scf = arr("scf_steps")

        ediffg = safe_float(self.settings.get("EDIFFG"))
        force_threshold = abs(ediffg) if ediffg is not None and ediffg < 0 else self.fallback_force_threshold
        energy_threshold = ediffg if ediffg is not None and ediffg > 0 else None
        nelm = self.settings.get("NELM")

        def add_hline(ax: Any, value: Optional[float], label: str) -> None:
            if value is not None and math.isfinite(value):
                ax.axhline(value, linestyle="--", label=label)
                ax.legend()

        fig, axes = plt.subplots(3, 2, figsize=(14, 12), sharex=True)
        fig.suptitle(f"VASP Optimization Audit: {len(self.steps)} ionic steps", fontsize=15)

        ax = axes[0, 0]
        ax.plot(x, energy, marker="o", markersize=3)
        ax.set_title("Total Energy")
        ax.set_ylabel("Energy (eV)")
        ax.grid(True)

        ax = axes[0, 1]
        ax.plot(x, delta_energy, marker="o", markersize=3)
        positive = delta_energy[np.isfinite(delta_energy) & (delta_energy > 0)]
        if len(positive):
            ax.set_yscale("log")
        add_hline(ax, energy_threshold, f"EDIFFG={energy_threshold:g} eV" if energy_threshold is not None else "EDIFFG")
        ax.set_title("|ΔE| between ionic steps")
        ax.set_ylabel("|ΔE| (eV)")
        ax.grid(True)

        ax = axes[1, 0]
        ax.plot(x, force, marker="o", markersize=3)
        positive = force[np.isfinite(force) & (force > 0)]
        if len(positive):
            ax.set_yscale("log")
        add_hline(ax, force_threshold, f"force threshold={force_threshold:g} eV/Å")
        ax.set_title("Maximum Atomic Force")
        ax.set_ylabel("Fmax (eV/Å)")
        ax.grid(True)

        ax = axes[1, 1]
        ax.plot(x, drift, marker="o", markersize=3)
        add_hline(ax, self.drift_warning_threshold, f"drift warning={self.drift_warning_threshold:g}")
        ax.set_title("Total Drift (warning only)")
        ax.set_ylabel("Drift (eV/Å)")
        ax.grid(True)

        ax = axes[2, 0]
        ax.bar(x, scf)
        add_hline(ax, float(nelm) if isinstance(nelm, int) else None, f"NELM={nelm}")
        ax.set_title("Electronic Steps per Ionic Step")
        ax.set_xlabel("Ionic Step")
        ax.set_ylabel("SCF Steps")
        ax.grid(True, axis="y")

        ax = axes[2, 1]
        ax.plot(x, pressure, marker="o", markersize=3)
        ax.set_title("External Pressure (diagnostic only)")
        ax.set_xlabel("Ionic Step")
        ax.set_ylabel("Pressure (kB)")
        ax.grid(True)

        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(self.output_dir / "plot_summary.png", dpi=180)
        plt.close(fig)

        series = [
            ("plot_energy.png", energy, "Total Energy", "Energy (eV)", None, False, "line"),
            ("plot_delta_energy.png", delta_energy, "|ΔE| between Ionic Steps", "|ΔE| (eV)", energy_threshold, True, "line"),
            ("plot_force.png", force, "Maximum Atomic Force", "Fmax (eV/Å)", force_threshold, True, "line"),
            ("plot_drift.png", drift, "Total Drift (warning only)", "Drift (eV/Å)", self.drift_warning_threshold, False, "line"),
            ("plot_nelm.png", scf, "Electronic Steps per Ionic Step", "SCF Steps", float(nelm) if isinstance(nelm, int) else None, False, "bar"),
            ("plot_pressure.png", pressure, "External Pressure (diagnostic only)", "Pressure (kB)", None, False, "line"),
        ]

        for filename, y, title, ylabel, threshold, use_log, kind in series:
            fig, ax = plt.subplots(figsize=(8, 5.5))
            if kind == "bar":
                ax.bar(x, y)
            else:
                ax.plot(x, y, marker="o", markersize=3)
            if use_log:
                positive = y[np.isfinite(y) & (y > 0)]
                if len(positive):
                    ax.set_yscale("log")
            add_hline(ax, threshold, f"threshold={threshold:g}" if threshold is not None else "threshold")
            ax.set_title(title)
            ax.set_xlabel("Ionic Step")
            ax.set_ylabel(ylabel)
            ax.grid(True)
            fig.tight_layout()
            fig.savefig(self.output_dir / filename, dpi=180)
            plt.close(fig)

    def print_report(self, assessment: dict[str, Any], export_source: Optional[str]) -> None:
        print("=" * 72)
        print("VaspDoctor v2.6 — VASP 结构优化审核")
        print("=" * 72)
        print(f"OUTCAR 力块数量 : {len(self.position_events)}")
        print(f"OSZICAR 完成步数: {len(self.oszicar_ionic_energies)}")
        print(f"分析记录数      : {len(self.steps)}")
        print(f"可靠能量配对步数: {self.reliably_paired_energy_steps}")
        print(f"仅含力记录数    : {self.force_only_analysis_steps}")
        print(f"原子总数        : {self.nions if self.nions is not None else 'N/A'}")
        print(f"状态            : {assessment['status']}")
        print(f"按参数推断满足  : {assessment.get('criterion_met_by_parsed_values')}")
        print(f"判断依据        : {assessment.get('convergence_basis')}")
        print(f"正式收敛确认    : {assessment.get('formal_convergence_confirmed')}")
        print(f"VASP 停止信息   : {assessment.get('official_vasp_stop_message')}")
        if export_source:
            print(f"结构导出来源    : {export_source}")

        if self.settings:
            selected = ["EDIFFG", "EDIFF", "NELM", "NSW", "IBRION", "ISIF", "ISYM"]
            text = ", ".join(
                f"{key}={self.settings[key]}" for key in selected if key in self.settings
            )
            if text:
                print(f"关键参数        : {text}")

        print("离子步能量来源  : OUTCAR TOTEN locked to completed POSITION TOTAL-FORCE blocks")
        print("OSZICAR F= 用途 : independent cross-check; fallback only when OUTCAR locking fails")

        if self.steps:
            last = self.steps[-1]
            print("\n最后一个离子步")
            print("-" * 72)
            print(f"Step            : {last.step}")
            print(f"Energy          : {format_optional(last.energy)} eV")
            print(f"Energy source   : {last.energy_source or 'N/A'}")
            if self.energy_events:
                raw = self.energy_events[-1]
                print(
                    f"OUTCAR tail TOTEN: {raw.value:.8f} eV "
                    f"(line {raw.line_no}; raw tail reference, not blindly used)"
                )
            if self.last_energy_crosscheck is not None:
                audit = self.last_energy_crosscheck
                print(
                    "Energy crosscheck : "
                    f"OUTCAR locked={audit['outcar_locked_toten_eV']:.8f} eV; "
                    f"OSZICAR F={audit['oszicar_F_eV']:.8f} eV; "
                    f"Δ={audit['outcar_minus_oszicar_eV']:.8f} eV"
                )
            print(f"|ΔE|            : {format_optional(last.delta_energy)} eV")
            print(
                f"Fmax            : {last.max_force:.6f} eV/Å "
                f"(atom #{last.max_atom_idx + 1}, {self.atom_label(last.max_atom_idx)})"
            )
            print(f"Drift           : {format_optional(last.drift)} eV/Å")
            print(f"Pressure        : {format_optional(last.pressure)} kB")
            print(f"SCF steps       : {last.scf_steps if last.scf_steps is not None else 'N/A'}")
            print(f"Incomplete SCF  : {self.oszicar_incomplete_scf_steps}")

            print("\n快照收敛估计")
            print("-" * 72)
            progress = assessment.get("progress_estimate", {})
            print(f"Metric          : {progress.get('metric')}")
            print(f"Current / target: {format_optional(progress.get('current'), '.6g')} / {format_optional(progress.get('target'), '.6g')} {progress.get('metric_unit') or ''}")
            print(f"Estimate usable : {progress.get('available')}")
            print(f"Remaining steps : {progress.get('estimated_remaining_ionic_steps')}")
            print(f"Expected step   : {progress.get('estimated_completion_ionic_step')}")
            print(f"Median LOOP+    : {format_optional(progress.get('median_loop_plus_walltime_seconds'), '.1f')} s / ionic step")
            print(f"Remaining time  : {progress.get('estimated_remaining_walltime_human')}")
            print(f"Confidence      : {progress.get('confidence')}")
            print(f"Fit in NSW      : {progress.get('estimated_to_fit_within_current_nsw')}")
            print(f"Reason          : {progress.get('reason')}")

            print("\n高频瓶颈原子")
            print("-" * 72)
            for item in self.bottleneck_atoms(top_n=5):
                vector = item["force_vector_at_historical_max_eV_per_A"]
                print(
                    f"#{item['atom_index_1based']} ({item['element']}): "
                    f"{item['bottleneck_count']}/{len(self.steps)} 步成为最大受力原子；"
                    f"历史最大力={item['historical_max_force_eV_per_A']:.6f} eV/Å "
                    f"(step {item['historical_max_force_step']}); "
                    f"F=({vector[0]:.6f}, {vector[1]:.6f}, {vector[2]:.6f})"
                )

        if not self.steps:
            print("\n快照收敛估计")
            print("-" * 72)
            print(f"Incomplete SCF  : {self.oszicar_incomplete_scf_steps}")
            print("Estimate usable : False")
            print("Remaining time  : N/A")
            print("Reason          : 尚无已完成离子步；至少需要 3 个完整离子步才能做趋势外推。")

        if self.error_hits:
            print("\n检测到高风险错误关键词")
            print("-" * 72)
            for name, lines in self.error_hits.items():
                preview = ", ".join(str(value) for value in lines[:10])
                suffix = " ..." if len(lines) > 10 else ""
                print(f"{name}: lines {preview}{suffix}")
                for item in self.error_context.get(name, [])[:3]:
                    print(f"  L{item['line']}: {item['text']}")

        if self.suspicious_hits:
            print("\n检测到需人工确认的文本（不自动判定失败）")
            print("-" * 72)
            for name, lines in self.suspicious_hits.items():
                preview = ", ".join(str(value) for value in lines[:10])
                suffix = " ..." if len(lines) > 10 else ""
                print(f"{name}: lines {preview}{suffix}")
                for item in self.suspicious_context.get(name, [])[:3]:
                    print(f"  L{item['line']}: {item['text']}")

        if self.benign_error_hits:
            total = sum(len(lines) for lines in self.benign_error_hits.values())
            print(f"\n已忽略已知无害的 ERROR 文本: {total} 行（例如 POTCAR kinetic energy error 元数据）")

        if self.info_messages:
            print("\n提示")
            print("-" * 72)
            for message in self.info_messages:
                print(f"- {message}")

        if self.warnings:
            print("\n警告")
            print("-" * 72)
            for warning in self.warnings:
                print(f"- {warning}")

        print("\n输出目录")
        print("-" * 72)
        print(self.output_dir.resolve())

    def run(self, make_plots: bool, export_structure: bool) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.read_inputs()
        self.parse_outcar()
        self.parse_oszicar()
        self.build_steps()
        assessment = self.assess()

        export_source: Optional[str] = None
        if export_structure:
            export_path = self.output_dir / "POSCAR_latest.vasp"
            has_nonempty_contcar = self.contcar_path.is_file() and self.contcar_path.stat().st_size > 0
            if self.steps or has_nonempty_contcar:
                try:
                    export_source = self.export_latest_structure(export_path)
                except Exception as exc:
                    self.warnings.append(f"结构导出失败：{exc}")
            else:
                self.info_messages.append(
                    "尚无已完成离子步且 CONTCAR 不可用；暂不导出 POSCAR_latest.vasp。"
                )

        self.write_csv(self.output_dir / "ionic_steps.csv")
        if make_plots:
            self.plot()
        self.write_json(self.output_dir / "summary.json", assessment, export_source)
        self.write_snapshot_status(self.output_dir / "snapshot_status.txt", assessment)
        self.print_report(assessment, export_source)
        return assessment


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit VASP optimization convergence and export structured reports.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("outcar", nargs="?", default="OUTCAR", help="OUTCAR file")
    parser.add_argument("oszicar", nargs="?", default="OSZICAR", help="OSZICAR file")
    parser.add_argument("poscar", nargs="?", default="POSCAR", help="POSCAR file")
    parser.add_argument("--incar", default="INCAR", help="INCAR file")
    parser.add_argument("--contcar", default="CONTCAR", help="CONTCAR file")
    parser.add_argument("--output-dir", default="vasp_doctor_report", help="Report directory")
    parser.add_argument(
        "--force-threshold",
        type=float,
        default=0.02,
        help="Fallback heuristic Fmax threshold when EDIFFG is unavailable (eV/Å)",
    )
    parser.add_argument(
        "--drift-warning",
        type=float,
        default=0.005,
        help="Auxiliary total-drift warning threshold (eV/Å)",
    )
    parser.add_argument(
        "--no-prefer-contcar",
        action="store_true",
        help="Rebuild latest POSCAR from OUTCAR even when CONTCAR exists",
    )
    parser.add_argument(
        "--trend-window",
        type=int,
        default=8,
        help="Number of recent completed ionic steps used for trend-based ETA",
    )
    parser.add_argument("--no-plots", action="store_true", help="Skip PNG plots")
    parser.add_argument("--no-export", action="store_true", help="Skip POSCAR_latest.vasp export")
    return parser


def make_doctor(args: argparse.Namespace) -> VaspDoctor:
    return VaspDoctor(
        outcar_path=Path(args.outcar),
        oszicar_path=Path(args.oszicar),
        poscar_path=Path(args.poscar),
        incar_path=Path(args.incar),
        contcar_path=Path(args.contcar),
        output_dir=Path(args.output_dir),
        fallback_force_threshold=args.force_threshold,
        drift_warning_threshold=args.drift_warning,
        prefer_contcar=not args.no_prefer_contcar,
        trend_window=args.trend_window,
    )


def main() -> int:
    args = build_parser().parse_args()
    doctor = make_doctor(args)
    try:
        doctor.run(make_plots=not args.no_plots, export_structure=not args.no_export)
    except Exception as exc:
        print(f"FATAL: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
