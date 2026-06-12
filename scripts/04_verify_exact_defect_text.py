#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_verify_exact_defect_text.py

严格审计 doped Stage-1 缺陷 POSCAR：
- pristine bulk 必须与 authoritative bulk 逐字节一致；
- vacancy：除删除一个目标坐标行、对应元素计数减 1 外，所有坐标文本必须原样保留；
- substitution：坐标文本集合必须与 bulk 完全一致，只允许元素归属和计数变化；
- interstitial：bulk 坐标文本必须全部原样保留，只允许新增一个坐标行；
- 所有缺陷 POSCAR 的比例因子、三行晶格矢量、Selective dynamics 行和坐标模式行
  必须与 authoritative bulk 逐字符一致。

脚本不调用 pymatgen 或 ASE，不进行浮点重写。它只比较文本。

典型用法：
    python .\scripts\04_verify_exact_defect_text.py `
        --bulk .\bulk_POSCAR `
        --workflow-root .\doped_snb_workflow

退出码：
    0：全部通过；
    2：至少一个文件不满足严格要求。
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any


def fail(message: str) -> "NoReturn":  # type: ignore[name-defined]
    raise SystemExit(f"ERROR: {message}")


@dataclass(frozen=True)
class RawPoscar:
    path: Path
    lines: list[str]
    comment: str
    scale_line: str
    lattice_lines: list[str]
    symbols_line: str
    counts_line: str
    symbols: list[str]
    counts: list[int]
    selective_line: str | None
    coordinate_mode_line: str
    coord_lines: list[str]
    trailing_lines: list[str]

    @classmethod
    def from_file(cls, path: Path) -> "RawPoscar":
        if not path.is_file():
            fail(f"POSCAR not found: {path}")
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
        if len(lines) < 8:
            fail(f"POSCAR too short: {path}")

        symbols_line = lines[5]
        counts_line = lines[6]
        symbols = symbols_line.split()
        counts_text = counts_line.split()

        if not symbols or not counts_text:
            fail(f"Missing symbols/counts rows: {path}")
        if all(token.lstrip("+-").isdigit() for token in symbols):
            fail(f"VASP4-style POSCAR is not supported: {path}")
        if not all(token.lstrip("+-").isdigit() for token in counts_text):
            fail(f"Could not parse counts row in {path}: {counts_line!r}")

        counts = [int(token) for token in counts_text]
        if len(symbols) != len(counts):
            fail(f"Symbols/counts length mismatch: {path}")

        cursor = 7
        selective_line: str | None = None
        if cursor < len(lines) and lines[cursor].strip().lower().startswith("s"):
            selective_line = lines[cursor]
            cursor += 1

        if cursor >= len(lines):
            fail(f"Missing coordinate mode row: {path}")
        coordinate_mode_line = lines[cursor]
        if not coordinate_mode_line.strip().lower().startswith(("d", "c", "k")):
            fail(f"Unknown coordinate mode row in {path}: {coordinate_mode_line!r}")
        cursor += 1

        natoms = sum(counts)
        coord_lines = lines[cursor: cursor + natoms]
        if len(coord_lines) != natoms:
            fail(
                f"Expected {natoms} coordinate rows but found {len(coord_lines)}: {path}"
            )

        trailing = lines[cursor + natoms:]
        return cls(
            path=path,
            lines=lines,
            comment=lines[0],
            scale_line=lines[1],
            lattice_lines=lines[2:5],
            symbols_line=symbols_line,
            counts_line=counts_line,
            symbols=symbols,
            counts=counts,
            selective_line=selective_line,
            coordinate_mode_line=coordinate_mode_line,
            coord_lines=coord_lines,
            trailing_lines=trailing,
        )

    @property
    def atom_count(self) -> int:
        return sum(self.counts)

    @property
    def lattice_block_sha256(self) -> str:
        text = "\n".join([self.scale_line, *self.lattice_lines]) + "\n"
        return hashlib.sha256(text.encode("utf-8")).hexdigest()


def compare_common_structure_text(bulk: RawPoscar, defect: RawPoscar) -> list[str]:
    errors: list[str] = []
    if defect.scale_line != bulk.scale_line:
        errors.append("scale line changed")
    if defect.lattice_lines != bulk.lattice_lines:
        errors.append("one or more lattice-vector lines changed")
    if defect.selective_line != bulk.selective_line:
        errors.append("Selective dynamics line changed")
    if defect.coordinate_mode_line != bulk.coordinate_mode_line:
        errors.append("coordinate mode line changed")
    if any(line.strip() for line in defect.trailing_lines):
        errors.append("defect POSCAR contains non-empty trailing rows")
    return errors


def counter_difference(left: Counter[str], right: Counter[str]) -> Counter[str]:
    """Return positive multiplicities in left - right."""
    return left - right


def verify_defect(
    *,
    bulk: RawPoscar,
    defect: RawPoscar,
    defect_type: str,
) -> tuple[list[str], dict[str, Any]]:
    errors = compare_common_structure_text(bulk, defect)
    bulk_rows = Counter(bulk.coord_lines)
    defect_rows = Counter(defect.coord_lines)

    removed_rows = counter_difference(bulk_rows, defect_rows)
    added_rows = counter_difference(defect_rows, bulk_rows)
    delta_atoms = defect.atom_count - bulk.atom_count

    normalized_type = defect_type.lower().strip()
    if normalized_type == "vacancy":
        if delta_atoms != -1:
            errors.append(f"vacancy atom-count delta must be -1, observed {delta_atoms}")
        if sum(removed_rows.values()) != 1:
            errors.append(
                "vacancy must remove exactly one original coordinate row, "
                f"observed {sum(removed_rows.values())}"
            )
        if sum(added_rows.values()) != 0:
            errors.append(
                "vacancy unexpectedly adds or modifies coordinate rows: "
                f"{sum(added_rows.values())} added rows"
            )

    elif normalized_type == "substitution":
        if delta_atoms != 0:
            errors.append(f"substitution atom-count delta must be 0, observed {delta_atoms}")
        if removed_rows or added_rows:
            errors.append(
                "substitution must preserve every coordinate row verbatim; "
                f"removed={sum(removed_rows.values())}, added={sum(added_rows.values())}"
            )

    elif normalized_type == "interstitial":
        if delta_atoms != 1:
            errors.append(f"interstitial atom-count delta must be +1, observed {delta_atoms}")
        if sum(removed_rows.values()) != 0:
            errors.append(
                "interstitial unexpectedly removes or modifies original coordinate rows: "
                f"{sum(removed_rows.values())} removed rows"
            )
        if sum(added_rows.values()) != 1:
            errors.append(
                "interstitial must add exactly one new coordinate row, "
                f"observed {sum(added_rows.values())}"
            )

    else:
        errors.append(f"unsupported or unknown defect type: {defect_type!r}")

    details = {
        "bulk_atom_count": bulk.atom_count,
        "defect_atom_count": defect.atom_count,
        "delta_atoms": delta_atoms,
        "removed_coordinate_rows": sum(removed_rows.values()),
        "added_coordinate_rows": sum(added_rows.values()),
        "lattice_text_identical": defect.lattice_block_sha256 == bulk.lattice_block_sha256,
    }
    return errors, details


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit exact Stage-1 defect POSCAR text preservation."
    )
    parser.add_argument("--bulk", required=True, help="Authoritative bulk POSCAR")
    parser.add_argument(
        "--workflow-root",
        default="doped_snb_workflow",
        help="Workflow output root; default: doped_snb_workflow",
    )
    args = parser.parse_args()

    bulk_path = Path(args.bulk).expanduser().resolve()
    root = Path(args.workflow_root).expanduser().resolve()
    stage_doped = root / "01_doped_defects"
    by_charge = stage_doped / "by_charge"
    reports = root / "00_reports"

    bulk = RawPoscar.from_file(bulk_path)

    if any(line.strip() for line in bulk.trailing_lines):
        fail(
            "Authoritative bulk POSCAR contains non-empty trailing rows. "
            "Remove MD velocities or auxiliary tail rows before strict export."
        )

    exported_bulk = stage_doped / "bulk_supercell_POSCAR"
    if not exported_bulk.is_file():
        fail(f"Missing exported pristine bulk: {exported_bulk}")

    bulk_byte_identical = bulk_path.read_bytes() == exported_bulk.read_bytes()

    if not by_charge.is_dir():
        fail(f"Missing Stage-1 defect directory: {by_charge}")

    rows: list[dict[str, Any]] = []
    failures = 0

    for defect_dir in sorted(path for path in by_charge.iterdir() if path.is_dir()):
        poscar_path = defect_dir / "POSCAR"
        meta_path = defect_dir / "DEFECT_ENTRY_META.json"
        if not poscar_path.is_file() or not meta_path.is_file():
            failures += 1
            rows.append(
                {
                    "defect_directory": str(defect_dir),
                    "status": "ERROR",
                    "errors": "missing POSCAR or DEFECT_ENTRY_META.json",
                }
            )
            continue

        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        defect_type = str(meta.get("defect_type", "unknown"))
        defect = RawPoscar.from_file(poscar_path)
        errors, details = verify_defect(
            bulk=bulk,
            defect=defect,
            defect_type=defect_type,
        )
        status = "OK" if not errors else "ERROR"
        if errors:
            failures += 1

        row: dict[str, Any] = {
            "defect_directory": str(defect_dir),
            "entry_name": meta.get("entry_name", defect_dir.name),
            "defect_type": defect_type,
            "charge_state": meta.get("charge_state", ""),
            "status": status,
            "errors": " | ".join(errors),
        }
        row.update(details)
        rows.append(row)

    report = reports / "strict_stage1_text_audit.csv"
    write_csv(report, rows)

    print("=" * 88)
    print("Strict Stage-1 POSCAR text audit")
    print("=" * 88)
    print(f"authoritative_bulk        = {bulk_path}")
    print(f"exported_bulk             = {exported_bulk}")
    print(f"bulk_byte_identical       = {bulk_byte_identical}")
    print(f"defect_directories        = {len(rows)}")
    print(f"failed_defect_directories = {failures}")
    print(f"audit                     = {report}")
    print("=" * 88)

    if not bulk_byte_identical:
        print("ERROR: Exported pristine bulk differs byte-for-byte from authoritative bulk.")
        raise SystemExit(2)
    if failures:
        print("ERROR: At least one Stage-1 defect POSCAR violated strict text preservation.")
        raise SystemExit(2)

    print("PASS: Every Stage-1 defect POSCAR satisfies strict text preservation.")
    raise SystemExit(0)


if __name__ == "__main__":
    main()
