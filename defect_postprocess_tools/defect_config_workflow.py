# -*- coding: utf-8 -*-
"""
General defect Spinney + KO semi-automatic configuration workflow.

Python: 3.11+
Dependencies: standard library only.

Expected project layout by default:

project_root/
├─ supercell/
│  ├─ POSCAR
│  ├─ OUTCAR
│  └─ vasprun.xml
├─ defect_A/
│  ├─ POSCAR                         # unrelaxed initial defect structure
│  ├─ charge_state_0/scf/INCAR
│  ├─ charge_state_0/scf/OUTCAR
│  └─ charge_state_q/scf/OUTCAR
└─ _defect_postprocess/
   ├─ defect_components.toml
   ├─ ko_center_reports/
   ├─ generated_defect_db.py
   ├─ logs/
   └─ backups/

The script never modifies POSCAR, INCAR, OUTCAR, vasprun.xml, the finder,
or the Spinney script.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import pprint
import py_compile
import re
import shutil
import socket
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

# -----------------------------------------------------------------------------
# Constants and regular expressions
# -----------------------------------------------------------------------------

WORKFLOW_VERSION = 1
CHARGE_DIR_RE = re.compile(r"^charge_state_([+-]?\d+)$")
NELECT_RE = re.compile(
    r"(?im)^\s*NELECT\s*=\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
)
INTEGER_LINE_RE = re.compile(r"^\s*(?:\d+\s+)*\d+\s*$")

COMPONENT_PATTERNS = {
    "substitution": re.compile(
        r"^(?:(?P<multiplicity>\d+))?sub_"
        r"(?P<new>[A-Z][A-Za-z0-9]*)_(?P<host>[A-Z][A-Za-z0-9]*)$"
    ),
    "vacancy": re.compile(
        r"^(?:(?P<multiplicity>\d+))?vac_(?P<species>[A-Z][A-Za-z0-9]*)$"
    ),
    "interstitial": re.compile(
        r"^(?:(?P<multiplicity>\d+))?int_(?P<species>[A-Z][A-Za-z0-9]*)$"
    ),
}

NORMAL_END_MARKERS = (
    "General timing and accounting informations for this job",
    "General timing and accounting information for this job",
)
TOTEN_MARKER = "free  energy   TOTEN"
CORE_POTENTIAL_MARKERS = (
    "electrostatic potential at core",
    "potential at core",
)

DEFAULT_IGNORES = {
    "supercell",
    "_defect_postprocess",
    "__pycache__",
    ".git",
    ".agents",
    "tests",
    "plot",
    "unused",
    "不用",
}

# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class ProjectPaths:
    root: Path
    tool_script: Path
    finder: Path
    bulk_poscar: Path
    bulk_outcar: Path
    bulk_vasprun: Path
    defect_poscar_name: str
    neutral_incar_rel: Path
    charge_outcar_template: str
    postprocess: Path
    components_toml: Path
    reports_dir: Path
    logs_dir: Path
    backups_dir: Path
    generated_python: Path
    generated_report: Path
    generated_validation: Path
    scan_report: Path
    check_report: Path
    audit_dir: Path
    audit_summary: Path
    species_groups_toml: Path


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def timestamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_json(data: Any) -> str:
    payload = json.dumps(
        data,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f"{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_json(path: Path, data: Any) -> None:
    atomic_write_text(
        path,
        json.dumps(data, ensure_ascii=False, indent=2) + "\n",
    )


def backup_file(path: Path, backup_dir: Path) -> Path | None:
    if not path.exists():
        return None
    backup_dir.mkdir(parents=True, exist_ok=True)
    candidate = backup_dir / f"{path.name}.{timestamp()}.bak"
    counter = 1
    while candidate.exists():
        candidate = backup_dir / f"{path.name}.{timestamp()}.{counter}.bak"
        counter += 1
    shutil.copy2(path, candidate)
    return candidate


def finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def sort_charges(charges: Iterable[int]) -> list[int]:
    def key(q: int) -> tuple[int, int]:
        if q == 0:
            return (0, 0)
        if q < 0:
            return (1, abs(q))
        return (2, q)

    return sorted(set(charges), key=key)


def safe_relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


# -----------------------------------------------------------------------------
# TOML read/write
# -----------------------------------------------------------------------------


def load_toml(path: Path) -> dict[str, Any]:
    import tomllib

    with path.open("rb") as handle:
        data = tomllib.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"TOML root must be a table: {path}")
    return data


def toml_quote(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return toml_quote(value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, list):
        return "[" + ", ".join(toml_value(item) for item in value) + "]"
    raise TypeError(f"Unsupported TOML value: {value!r}")


def dump_components_toml(data: dict[str, Any]) -> str:
    lines = [
        "# Defect components maintained manually.",
        "# Edit enabled, components, label and optional correction_center_mode.",
        "# component forms: sub_Zn_Cu, 2sub_Zn_Cu, vac_Cu, int_O",
        "# correction_center_mode: auto (default) or skip",
        "",
    ]
    for defect_name in sorted(data):
        config = data[defect_name]
        if not isinstance(config, dict):
            raise TypeError(f"Configuration for {defect_name!r} must be a table.")
        lines.append(f"[{toml_quote(defect_name)}]")
        lines.append(f"enabled = {toml_value(config.get('enabled', True))}")
        lines.append(f"components = {toml_value(config.get('components', []))}")
        lines.append(f"label = {toml_value(config.get('label', ''))}")
        # Preserve simple additional scalar/list fields for forward compatibility.
        for key in sorted(set(config) - {"enabled", "components", "label"}):
            lines.append(f"{key} = {toml_value(config[key])}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


# -----------------------------------------------------------------------------
# POSCAR and component parsing
# -----------------------------------------------------------------------------


def parse_poscar_composition(path: Path) -> dict[str, int]:
    """
    Parse VASP 5/6 POSCAR element names and counts.

    VASP 4 POSCAR files without an element-symbol line are rejected because
    component validation requires explicit chemical symbols.
    """
    if not path.is_file():
        raise FileNotFoundError(f"POSCAR not found: {path}")
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 7:
        raise ValueError(f"POSCAR is too short: {path}")

    symbols = lines[5].split()
    if not symbols:
        raise ValueError(f"Empty element-symbol line in POSCAR: {path}")
    if all(token.lstrip("+-").isdigit() for token in symbols):
        raise ValueError(
            f"VASP 4 style POSCAR without element symbols is unsupported: {path}"
        )

    count_tokens = lines[6].split()
    if len(symbols) != len(count_tokens):
        raise ValueError(
            f"Element/count length mismatch in POSCAR {path}: "
            f"{symbols!r} versus {count_tokens!r}"
        )
    try:
        counts = [int(token) for token in count_tokens]
    except ValueError as exc:
        raise ValueError(f"Invalid atom-count line in POSCAR {path}") from exc

    if any(value < 0 for value in counts):
        raise ValueError(f"Negative atom count in POSCAR: {path}")
    if len(set(symbols)) != len(symbols):
        raise ValueError(f"Duplicate element symbols in POSCAR: {path}")

    return dict(zip(symbols, counts))


def composition_delta(
    bulk: dict[str, int],
    defect: dict[str, int],
) -> dict[str, int]:
    """Return defect composition minus bulk composition."""
    species = sorted(set(bulk) | set(defect))
    return {
        element: defect.get(element, 0) - bulk.get(element, 0)
        for element in species
        if defect.get(element, 0) - bulk.get(element, 0) != 0
    }


def parse_component(component: str) -> dict[str, Any]:
    for kind, pattern in COMPONENT_PATTERNS.items():
        match = pattern.fullmatch(component)
        if not match:
            continue
        multiplicity = int(match.group("multiplicity") or 1)
        if multiplicity <= 0:
            raise ValueError(f"Component multiplicity must be positive: {component}")
        if kind == "substitution":
            return {
                "raw": component,
                "kind": kind,
                "multiplicity": multiplicity,
                "new_species": match.group("new"),
                "host_species": match.group("host"),
            }
        return {
            "raw": component,
            "kind": kind,
            "multiplicity": multiplicity,
            "species": match.group("species"),
        }
    raise ValueError(
        f"Unsupported component {component!r}. Valid forms include "
        "sub_Zn_Cu, 2sub_Zn_Cu, vac_Cu, 2vac_Cu, int_O and 2int_O."
    )


def expected_delta_from_components(components: list[str]) -> dict[str, int]:
    delta: Counter[str] = Counter()
    for raw in components:
        component = parse_component(raw)
        multiplicity = component["multiplicity"]
        if component["kind"] == "vacancy":
            delta[component["species"]] -= multiplicity
        elif component["kind"] == "interstitial":
            delta[component["species"]] += multiplicity
        elif component["kind"] == "substitution":
            delta[component["host_species"]] -= multiplicity
            delta[component["new_species"]] += multiplicity
        else:
            raise AssertionError(component)
    return dict(sorted((key, value) for key, value in delta.items() if value != 0))


def total_component_multiplicity(components: list[str]) -> int:
    return sum(parse_component(component)["multiplicity"] for component in components)



# -----------------------------------------------------------------------------
# Independent structure audit
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Structure:
    comment: str
    lattice: list[list[float]]
    symbols: list[str]
    frac_coords: list[list[float]]


@dataclass(frozen=True)
class AuditSite:
    species: str
    frac_coords: list[float]
    atom_indices: list[int]
    kind: str  # "atom" or "group"


def determinant3(matrix: list[list[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def inverse3(matrix: list[list[float]]) -> list[list[float]]:
    det = determinant3(matrix)
    if abs(det) < 1e-14:
        raise ValueError("Singular lattice matrix.")
    a, b, c = matrix
    cof = [
        [b[1] * c[2] - b[2] * c[1], -(b[0] * c[2] - b[2] * c[0]), b[0] * c[1] - b[1] * c[0]],
        [-(a[1] * c[2] - a[2] * c[1]), a[0] * c[2] - a[2] * c[0], -(a[0] * c[1] - a[1] * c[0])],
        [a[1] * b[2] - a[2] * b[1], -(a[0] * b[2] - a[2] * b[0]), a[0] * b[1] - a[1] * b[0]],
    ]
    # inverse = transpose(cofactor) / det
    return [[cof[j][i] / det for j in range(3)] for i in range(3)]


def row_vector_times_matrix(vector: list[float], matrix: list[list[float]]) -> list[float]:
    return [
        sum(vector[k] * matrix[k][j] for k in range(3))
        for j in range(3)
    ]


def normalize_frac(coords: list[float]) -> list[float]:
    return [value - math.floor(value) for value in coords]


def parse_poscar_structure(path: Path) -> Structure:
    if not path.is_file():
        raise FileNotFoundError(f"POSCAR not found: {path}")
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError(f"POSCAR is too short: {path}")

    comment = lines[0].strip()
    try:
        scale = float(lines[1].split()[0])
        raw_lattice = [
            [float(value) for value in lines[index].split()[:3]]
            for index in range(2, 5)
        ]
    except Exception as exc:
        raise ValueError(f"Invalid scale/lattice in POSCAR: {path}") from exc

    if scale == 0:
        raise ValueError(f"POSCAR scale factor cannot be zero: {path}")
    if scale > 0:
        lattice_factor = scale
    else:
        raw_volume = abs(determinant3(raw_lattice))
        if raw_volume <= 0:
            raise ValueError(f"Invalid lattice volume in POSCAR: {path}")
        lattice_factor = (abs(scale) / raw_volume) ** (1.0 / 3.0)
    lattice = [[lattice_factor * value for value in row] for row in raw_lattice]

    symbols = lines[5].split()
    if not symbols or all(token.lstrip("+-").isdigit() for token in symbols):
        raise ValueError(
            f"VASP 4 POSCAR without an explicit element-symbol line is unsupported: {path}"
        )
    try:
        counts = [int(token) for token in lines[6].split()]
    except ValueError as exc:
        raise ValueError(f"Invalid atom-count line in POSCAR: {path}") from exc
    if len(symbols) != len(counts):
        raise ValueError(f"Element/count length mismatch in POSCAR: {path}")

    cursor = 7
    if cursor < len(lines) and lines[cursor].strip().lower().startswith("s"):
        cursor += 1
    if cursor >= len(lines):
        raise ValueError(f"Missing coordinate mode in POSCAR: {path}")
    mode = lines[cursor].strip().lower()
    cursor += 1
    direct = mode.startswith("d")
    cartesian = mode.startswith("c") or mode.startswith("k")
    if not direct and not cartesian:
        raise ValueError(f"Unknown POSCAR coordinate mode {lines[cursor - 1]!r}: {path}")

    total = sum(counts)
    if len(lines) < cursor + total:
        raise ValueError(f"POSCAR has fewer coordinate rows than declared atoms: {path}")
    expanded_symbols: list[str] = []
    for symbol, count in zip(symbols, counts):
        expanded_symbols.extend([symbol] * count)

    inverse_lattice = inverse3(lattice)
    frac_coords: list[list[float]] = []
    for offset in range(total):
        tokens = lines[cursor + offset].split()
        if len(tokens) < 3:
            raise ValueError(f"Invalid coordinate row {cursor + offset + 1}: {path}")
        values = [float(tokens[index]) for index in range(3)]
        if direct:
            frac = values
        else:
            cart = [lattice_factor * value for value in values]
            frac = row_vector_times_matrix(cart, inverse_lattice)
        frac_coords.append(normalize_frac(frac))

    return Structure(
        comment=comment,
        lattice=lattice,
        symbols=expanded_symbols,
        frac_coords=frac_coords,
    )


def periodic_distance(
    frac_a: list[float],
    frac_b: list[float],
    lattice: list[list[float]],
) -> float:
    delta = [frac_a[index] - frac_b[index] for index in range(3)]
    delta = [value - round(value) for value in delta]
    cart = row_vector_times_matrix(delta, lattice)
    return math.sqrt(sum(value * value for value in cart))


def lattice_max_difference(
    lattice_a: list[list[float]],
    lattice_b: list[list[float]],
) -> float:
    return max(
        abs(lattice_a[i][j] - lattice_b[i][j])
        for i in range(3)
        for j in range(3)
    )


def minimum_pair_distance(structure: Structure) -> tuple[float | None, tuple[int, int] | None]:
    best_distance: float | None = None
    best_pair: tuple[int, int] | None = None
    for first in range(len(structure.symbols)):
        for second in range(first + 1, len(structure.symbols)):
            distance = periodic_distance(
                structure.frac_coords[first],
                structure.frac_coords[second],
                structure.lattice,
            )
            if best_distance is None or distance < best_distance:
                best_distance = distance
                best_pair = (first, second)
    return best_distance, best_pair


def load_species_groups(path: Path) -> dict[str, dict[str, Any]]:
    if not path.is_file():
        return {}
    raw = load_toml(path)
    groups: dict[str, dict[str, Any]] = {}
    for name, config in raw.items():
        if not isinstance(config, dict):
            raise ValueError(f"Species-group {name!r} must be a TOML table.")
        composition = config.get("composition")
        anchor = config.get("anchor_species")
        radius = config.get("radius", 2.5)
        if not isinstance(composition, dict) or not composition:
            raise ValueError(f"Species-group {name!r} needs a non-empty composition table.")
        normalized: dict[str, int] = {}
        for species, count in composition.items():
            if not isinstance(species, str) or not isinstance(count, int) or count <= 0:
                raise ValueError(f"Invalid composition in species-group {name!r}.")
            normalized[species] = count
        if not isinstance(anchor, str) or normalized.get(anchor) != 1:
            raise ValueError(
                f"Species-group {name!r} requires anchor_species with composition count 1."
            )
        if not finite_number(radius) or float(radius) <= 0:
            raise ValueError(f"Invalid radius for species-group {name!r}.")
        groups[name] = {
            "composition": normalized,
            "anchor_species": anchor,
            "radius": float(radius),
        }
    return groups


def build_audit_sites(
    structure: Structure,
    groups: dict[str, dict[str, Any]],
    group_ambiguity_tolerance: float,
) -> tuple[list[AuditSite], list[str]]:
    consumed: set[int] = set()
    sites: list[AuditSite] = []
    warnings: list[str] = []

    for group_name in sorted(groups):
        definition = groups[group_name]
        composition = definition["composition"]
        anchor_species = definition["anchor_species"]
        radius = definition["radius"]
        anchors = [
            index
            for index, symbol in enumerate(structure.symbols)
            if symbol == anchor_species and index not in consumed
        ]

        for anchor_index in anchors:
            members = [anchor_index]
            selected_now: set[int] = {anchor_index}
            complete = True
            ambiguous = False
            for species, count in sorted(composition.items()):
                remaining = count - (1 if species == anchor_species else 0)
                if remaining <= 0:
                    continue
                candidates: list[tuple[float, int]] = []
                for atom_index, symbol in enumerate(structure.symbols):
                    if symbol != species or atom_index in consumed or atom_index in selected_now:
                        continue
                    distance = periodic_distance(
                        structure.frac_coords[anchor_index],
                        structure.frac_coords[atom_index],
                        structure.lattice,
                    )
                    if distance <= radius:
                        candidates.append((distance, atom_index))
                candidates.sort()
                if len(candidates) < remaining:
                    complete = False
                    break
                chosen = candidates[:remaining]
                if len(candidates) > remaining:
                    boundary_gap = candidates[remaining][0] - candidates[remaining - 1][0]
                    if boundary_gap <= group_ambiguity_tolerance:
                        ambiguous = True
                for _, atom_index in chosen:
                    members.append(atom_index)
                    selected_now.add(atom_index)

            if not complete:
                continue
            if ambiguous:
                warnings.append(
                    f"Ambiguous membership while constructing group {group_name} "
                    f"around anchor atom {anchor_index}."
                )
            consumed.update(members)
            sites.append(
                AuditSite(
                    species=group_name,
                    frac_coords=structure.frac_coords[anchor_index],
                    atom_indices=sorted(members),
                    kind="group",
                )
            )

    for atom_index, symbol in enumerate(structure.symbols):
        if atom_index in consumed:
            continue
        sites.append(
            AuditSite(
                species=symbol,
                frac_coords=structure.frac_coords[atom_index],
                atom_indices=[atom_index],
                kind="atom",
            )
        )
    return sites, warnings


def hungarian_assignment(cost: list[list[float]]) -> list[int]:
    """Minimum-cost assignment for a square matrix; returns row -> column."""
    n = len(cost)
    if n == 0:
        return []
    if any(len(row) != n for row in cost):
        raise ValueError("Hungarian assignment requires a square matrix.")
    u = [0.0] * (n + 1)
    v = [0.0] * (n + 1)
    p = [0] * (n + 1)
    way = [0] * (n + 1)

    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [float("inf")] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = float("inf")
            j1 = 0
            for j in range(1, n + 1):
                if used[j]:
                    continue
                current = cost[i0 - 1][j - 1] - u[i0] - v[j]
                if current < minv[j]:
                    minv[j] = current
                    way[j] = j0
                if minv[j] < delta:
                    delta = minv[j]
                    j1 = j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break

    assignment = [-1] * n
    for column in range(1, n + 1):
        if p[column] != 0:
            assignment[p[column] - 1] = column - 1
    return assignment


def match_audit_sites(
    bulk_sites: list[AuditSite],
    defect_sites: list[AuditSite],
    lattice: list[list[float]],
    site_tolerance: float,
    ambiguity_tolerance: float,
) -> dict[str, Any]:
    nb = len(bulk_sites)
    nd = len(defect_sites)
    size = nb + nd
    if size == 0:
        return {"matches": [], "unmatched_bulk": [], "unmatched_defect": [], "ambiguities": []}

    big = 1.0e9
    unmatched_cost = site_tolerance
    distances = [
        [
            periodic_distance(bulk.frac_coords, defect.frac_coords, lattice)
            for defect in defect_sites
        ]
        for bulk in bulk_sites
    ]
    cost = [[big] * size for _ in range(size)]

    for i in range(nb):
        for j in range(nd):
            if distances[i][j] <= site_tolerance:
                cost[i][j] = distances[i][j]
        cost[i][nd + i] = unmatched_cost

    for j in range(nd):
        cost[nb + j][j] = unmatched_cost
        for dummy_column in range(nd, size):
            cost[nb + j][dummy_column] = 0.0

    assignment = hungarian_assignment(cost)
    matches: list[dict[str, Any]] = []
    unmatched_bulk: list[int] = []
    matched_defect: set[int] = set()

    for i in range(nb):
        column = assignment[i]
        if 0 <= column < nd and cost[i][column] < big / 2:
            matches.append({"bulk_index": i, "defect_index": column, "distance": distances[i][column]})
            matched_defect.add(column)
        else:
            unmatched_bulk.append(i)

    unmatched_defect = [j for j in range(nd) if j not in matched_defect]
    ambiguities: list[dict[str, Any]] = []
    for match in matches:
        i = match["bulk_index"]
        j = match["defect_index"]
        assigned = match["distance"]
        bulk_candidates = sorted(
            (distances[i][candidate], candidate)
            for candidate in range(nd)
            if distances[i][candidate] <= site_tolerance
        )
        defect_candidates = sorted(
            (distances[candidate][j], candidate)
            for candidate in range(nb)
            if distances[candidate][j] <= site_tolerance
        )
        if len(bulk_candidates) > 1 and bulk_candidates[1][0] - bulk_candidates[0][0] <= ambiguity_tolerance:
            ambiguities.append(
                {
                    "kind": "bulk_site_multiple_defect_candidates",
                    "bulk_site": i,
                    "assigned_defect_site": j,
                    "candidates": bulk_candidates[:5],
                }
            )
        if len(defect_candidates) > 1 and defect_candidates[1][0] - defect_candidates[0][0] <= ambiguity_tolerance:
            ambiguities.append(
                {
                    "kind": "defect_site_multiple_bulk_candidates",
                    "defect_site": j,
                    "assigned_bulk_site": i,
                    "candidates": defect_candidates[:5],
                }
            )

    return {
        "matches": matches,
        "unmatched_bulk": unmatched_bulk,
        "unmatched_defect": unmatched_defect,
        "ambiguities": ambiguities,
    }


def canonical_component(kind: str, *, host: str | None = None, new: str | None = None, species: str | None = None) -> str:
    if kind == "substitution":
        if not host or not new:
            raise ValueError("Substitution needs host and new species.")
        return f"sub_{new}_{host}"
    if kind == "vacancy":
        if not species:
            raise ValueError("Vacancy needs species.")
        return f"vac_{species}"
    if kind == "interstitial":
        if not species:
            raise ValueError("Interstitial needs species.")
        return f"int_{species}"
    raise ValueError(f"Unsupported event kind: {kind}")


def expand_components(components: list[str]) -> list[str]:
    expanded: list[str] = []
    for raw in components:
        parsed = parse_component(raw)
        if parsed["kind"] == "substitution":
            canonical = canonical_component(
                "substitution",
                host=parsed["host_species"],
                new=parsed["new_species"],
            )
        elif parsed["kind"] == "vacancy":
            canonical = canonical_component("vacancy", species=parsed["species"])
        else:
            canonical = canonical_component("interstitial", species=parsed["species"])
        expanded.extend([canonical] * parsed["multiplicity"])
    return sorted(expanded)


def expected_atomic_delta(
    components: list[str],
    groups: dict[str, dict[str, Any]],
) -> dict[str, int]:
    delta: Counter[str] = Counter()

    def add_species(species: str, amount: int) -> None:
        if species in groups:
            for element, count in groups[species]["composition"].items():
                delta[element] += amount * count
        else:
            delta[species] += amount

    for raw in components:
        parsed = parse_component(raw)
        multiplicity = parsed["multiplicity"]
        if parsed["kind"] == "vacancy":
            add_species(parsed["species"], -multiplicity)
        elif parsed["kind"] == "interstitial":
            add_species(parsed["species"], multiplicity)
        else:
            add_species(parsed["host_species"], -multiplicity)
            add_species(parsed["new_species"], multiplicity)
    return dict(sorted((species, amount) for species, amount in delta.items() if amount))


def independently_detect_events(
    bulk_structure: Structure,
    defect_structure: Structure,
    bulk_sites: list[AuditSite],
    defect_sites: list[AuditSite],
    match_result: dict[str, Any],
) -> tuple[list[str], list[dict[str, Any]]]:
    components: list[str] = []
    events: list[dict[str, Any]] = []

    for match in match_result["matches"]:
        bulk_site = bulk_sites[match["bulk_index"]]
        defect_site = defect_sites[match["defect_index"]]
        if bulk_site.species == defect_site.species:
            continue
        component = canonical_component(
            "substitution",
            host=bulk_site.species,
            new=defect_site.species,
        )
        components.append(component)
        events.append(
            {
                "kind": "substitution",
                "component": component,
                "host_species": bulk_site.species,
                "new_species": defect_site.species,
                "bulk_site_index": match["bulk_index"],
                "defect_site_index": match["defect_index"],
                "bulk_atom_indices": bulk_site.atom_indices,
                "defect_atom_indices": defect_site.atom_indices,
                "bulk_frac_coords": bulk_site.frac_coords,
                "defect_frac_coords": defect_site.frac_coords,
                "mapping_distance_angstrom": match["distance"],
            }
        )

    for index in match_result["unmatched_bulk"]:
        site = bulk_sites[index]
        component = canonical_component("vacancy", species=site.species)
        components.append(component)
        events.append(
            {
                "kind": "vacancy",
                "component": component,
                "species": site.species,
                "bulk_site_index": index,
                "bulk_atom_indices": site.atom_indices,
                "bulk_frac_coords": site.frac_coords,
            }
        )

    for index in match_result["unmatched_defect"]:
        site = defect_sites[index]
        component = canonical_component("interstitial", species=site.species)
        components.append(component)
        events.append(
            {
                "kind": "interstitial",
                "component": component,
                "species": site.species,
                "defect_site_index": index,
                "defect_atom_indices": site.atom_indices,
                "defect_frac_coords": site.frac_coords,
            }
        )

    return sorted(components), events


def compare_component_multisets(expected: list[str], detected: list[str]) -> dict[str, Any]:
    expected_counter = Counter(expected)
    detected_counter = Counter(detected)
    missing_counter = expected_counter - detected_counter
    unexpected_counter = detected_counter - expected_counter
    missing = sorted(component for component, count in missing_counter.items() for _ in range(count))
    unexpected = sorted(component for component, count in unexpected_counter.items() for _ in range(count))
    return {
        "match": not missing and not unexpected,
        "missing": missing,
        "unexpected": unexpected,
    }


def current_audit_metadata(
    paths: ProjectPaths,
    defect_dir: Path,
    expected_components: list[str],
) -> dict[str, Any]:
    defect_poscar = defect_dir / paths.defect_poscar_name
    return {
        "workflow_version": WORKFLOW_VERSION,
        "defect_name": defect_dir.name,
        "expected_components": expected_components,
        "expected_components_sha256": sha256_json(expected_components),
        "bulk_poscar_sha256": sha256_file(paths.bulk_poscar),
        "defect_poscar_sha256": sha256_file(defect_poscar),
        "species_groups_sha256": sha256_file(paths.species_groups_toml),
    }


def audit_report_is_current_pass(
    paths: ProjectPaths,
    defect_dir: Path,
    expected_components: list[str],
) -> tuple[bool, list[str]]:
    report_path = paths.audit_dir / f"{defect_dir.name}.json"
    if not report_path.is_file():
        return False, ["Structure-audit report is missing."]
    try:
        report = load_json_object(report_path)
        if report.get("status") != "PASS":
            return False, [f"Structure-audit status is {report.get('status')!r}, not PASS."]
        metadata = report.get("_workflow_metadata")
        expected = current_audit_metadata(paths, defect_dir, expected_components)
        if not isinstance(metadata, dict):
            return False, ["Structure-audit report lacks metadata."]
        differences = [
            f"{key}: report={metadata.get(key)!r}, current={expected.get(key)!r}"
            for key in expected
            if metadata.get(key) != expected.get(key)
        ]
        return not differences, differences
    except Exception as exc:
        return False, [f"Structure-audit report is unreadable: {exc}"]


def discover_structure_defects(paths: ProjectPaths) -> list[Path]:
    ignores = set(DEFAULT_IGNORES)
    ignores.add(paths.postprocess.name)
    ignores.add(Path(paths.bulk_poscar).parts[0])
    defects: list[Path] = []
    for child in sorted(paths.root.iterdir(), key=lambda item: item.name):
        if not child.is_dir() or child.name.startswith(".") or child.name in ignores:
            continue
        if (child / paths.defect_poscar_name).is_file():
            defects.append(child)
    return defects


# -----------------------------------------------------------------------------
# INCAR, charge-state and OUTCAR parsing
# -----------------------------------------------------------------------------


def strip_incar_comment(line: str) -> str:
    positions = [position for position in (line.find("!"), line.find("#")) if position >= 0]
    if positions:
        return line[: min(positions)]
    return line


def parse_nelect(path: Path) -> tuple[int | float | None, list[str], list[str]]:
    if not path.is_file():
        return None, [], [f"INCAR not found: {path}"]

    values: list[float] = []
    raw_lines: list[str] = []
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        active = strip_incar_comment(raw_line)
        match = NELECT_RE.search(active)
        if match:
            values.append(float(match.group(1)))
            raw_lines.append(raw_line)

    warnings: list[str] = []
    if not values:
        return None, raw_lines, [f"No explicit NELECT found in {path}"]
    if len(values) > 1:
        warnings.append(
            f"Multiple active NELECT assignments found in {path}; the last value is used."
        )
    value = values[-1]
    normalized: int | float = int(value) if value.is_integer() else value
    return normalized, raw_lines, warnings


def charge_directories(defect_dir: Path) -> dict[int, Path]:
    found: dict[int, Path] = {}
    for child in defect_dir.iterdir():
        if not child.is_dir():
            continue
        match = CHARGE_DIR_RE.fullmatch(child.name)
        if not match:
            continue
        charge = int(match.group(1))
        if charge in found:
            raise ValueError(f"Duplicate charge state {charge}: {defect_dir}")
        found[charge] = child
    return found


def outcar_path_for_charge(
    defect_dir: Path,
    charge: int,
    template: str,
) -> Path:
    try:
        relative = template.format(q=charge)
    except Exception as exc:
        raise ValueError(
            f"Invalid --charge-outcar template {template!r}; it must support {{q}}."
        ) from exc
    return defect_dir / Path(relative)


def validate_outcar(path: Path, strict: bool) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    if not path.is_file():
        return [f"OUTCAR missing: {path}"], warnings
    if path.stat().st_size == 0:
        return [f"OUTCAR is empty: {path}"], warnings

    text = path.read_text(encoding="utf-8", errors="replace")
    lower = text.lower()

    if TOTEN_MARKER.lower() not in lower:
        errors.append(f"OUTCAR lacks a TOTEN record: {path}")
    if not any(marker.lower() in lower for marker in NORMAL_END_MARKERS):
        message = f"OUTCAR lacks the standard VASP normal-end marker: {path}"
        (errors if strict else warnings).append(message)
    if not any(marker.lower() in lower for marker in CORE_POTENTIAL_MARKERS):
        message = f"OUTCAR lacks a recognizable core-potential block: {path}"
        (errors if strict else warnings).append(message)

    return errors, warnings


# -----------------------------------------------------------------------------
# Project discovery and paths
# -----------------------------------------------------------------------------


def build_paths(args: argparse.Namespace) -> ProjectPaths:
    root = (args.root if args.root is not None else Path.cwd()).resolve()
    tool_script = Path(__file__).resolve()
    finder = (
        args.finder.resolve()
        if args.finder
        else tool_script.with_name("find_cluster_ko_center_V1.py")
    )

    # v2.2 fixed layout:
    #   active configuration files are stored at the project root;
    #   all generated files are stored under _defect_postprocess.
    postprocess = root / args.postprocess_dir

    if getattr(args, "workspace", None) is not None:
        print(
            "WARNING: --workspace is deprecated and ignored in v2.2. "
            "Configuration is read from the project root and generated output "
            f"is written to {postprocess}.",
            file=sys.stderr,
        )

    return ProjectPaths(
        root=root,
        tool_script=tool_script,
        finder=finder,
        bulk_poscar=root / args.bulk_poscar,
        bulk_outcar=root / args.bulk_outcar,
        bulk_vasprun=root / args.bulk_vasprun,
        defect_poscar_name=args.defect_poscar_name,
        neutral_incar_rel=Path(args.neutral_incar),
        charge_outcar_template=args.charge_outcar,
        postprocess=postprocess,
        components_toml=root / args.components_file,
        reports_dir=postprocess / "ko_center_reports",
        logs_dir=postprocess / "logs",
        backups_dir=postprocess / "backups",
        generated_python=postprocess / "generated_defect_db.py",
        generated_report=postprocess / "generated_defect_db_report.json",
        generated_validation=postprocess / "generated_defect_db_validation.json",
        scan_report=postprocess / "scan_report.json",
        check_report=postprocess / "check_report.json",
        audit_dir=postprocess / "structure_audit",
        audit_summary=postprocess / "structure_audit_summary.json",
        species_groups_toml=root / args.species_groups_file,
    )


def initialize_directories(paths: ProjectPaths) -> None:
    paths.postprocess.mkdir(parents=True, exist_ok=True)
    paths.reports_dir.mkdir(parents=True, exist_ok=True)
    paths.logs_dir.mkdir(parents=True, exist_ok=True)
    paths.backups_dir.mkdir(parents=True, exist_ok=True)
    paths.audit_dir.mkdir(parents=True, exist_ok=True)
    if not paths.components_toml.exists():
        atomic_write_text(paths.components_toml, dump_components_toml({}))
    if not paths.species_groups_toml.exists():
        atomic_write_text(
            paths.species_groups_toml,
            "# Optional molecular/pseudo-species definitions.\n"
            "# Example:\n"
            "# [FA]\n"
            "# composition = { C = 1, H = 5, N = 2 }\n"
            "# anchor_species = \"C\"\n"
            "# radius = 2.2\n",
        )


def discover_defects(paths: ProjectPaths) -> list[Path]:
    ignores = set(DEFAULT_IGNORES)
    ignores.add(paths.postprocess.name)
    ignores.add(Path(paths.bulk_poscar).parts[0])

    defects: list[Path] = []
    for child in sorted(paths.root.iterdir(), key=lambda item: item.name):
        if not child.is_dir():
            continue
        if child.name.startswith(".") or child.name in ignores:
            continue
        try:
            charges = charge_directories(child)
        except ValueError:
            defects.append(child)
            continue
        if not charges or 0 not in charges:
            continue
        if not (child / paths.defect_poscar_name).is_file():
            continue
        defects.append(child)
    return defects


def sync_components_template(paths: ProjectPaths, defects: list[Path]) -> dict[str, Any]:
    if paths.components_toml.exists():
        config = load_toml(paths.components_toml)
    else:
        config = {}

    added: list[str] = []
    for defect_dir in defects:
        if defect_dir.name not in config:
            config[defect_dir.name] = {
                "enabled": True,
                "components": [],
                "label": "",
            }
            added.append(defect_dir.name)

    if added:
        backup_file(paths.components_toml, paths.backups_dir)
        atomic_write_text(paths.components_toml, dump_components_toml(config))

    return {
        "config": config,
        "added": added,
    }


# -----------------------------------------------------------------------------
# KO report parsing and freshness
# -----------------------------------------------------------------------------


def recursive_values(obj: Any, keys: set[str]) -> list[Any]:
    values: list[Any] = []
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key in keys:
                values.append(value)
            values.extend(recursive_values(value, keys))
    elif isinstance(obj, list):
        for value in obj:
            values.extend(recursive_values(value, keys))
    return values


def numeric_triplet(value: Any) -> list[float] | None:
    if isinstance(value, dict):
        for key in (
            "coords",
            "fractional_coordinates",
            "frac_coords",
            "center_frac",
            "correction_center_frac",
        ):
            if key in value:
                return numeric_triplet(value[key])
        return None

    if not isinstance(value, (list, tuple)) or len(value) != 3:
        return None
    if not all(finite_number(item) for item in value):
        return None
    return [float(item) for item in value]


def extract_correction_center(report: dict[str, Any]) -> list[float]:
    preferred_keys = {
        "recommended_correction_center",
        "correction_center",
        "recommended_center",
        "recommended_center_frac",
        "recommended_correction_center_frac",
        "centroid_frac",
    }
    candidates: list[list[float]] = []
    for value in recursive_values(report, preferred_keys):
        triplet = numeric_triplet(value)
        if triplet is not None:
            candidates.append(triplet)

    unique: list[list[float]] = []
    for candidate in candidates:
        if not any(
            all(abs(a - b) <= 1e-12 for a, b in zip(candidate, existing))
            for existing in unique
        ):
            unique.append(candidate)

    if len(unique) != 1:
        raise ValueError(
            f"Expected exactly one unique correction center; found {len(unique)}: {unique}"
        )
    return unique[0]


def extract_events(report: dict[str, Any]) -> list[Any]:
    candidates: list[list[Any]] = []
    for value in recursive_values(
        report,
        {"defect_events", "detected_events", "events", "identified_events"},
    ):
        if isinstance(value, list):
            candidates.append(value)
    nonempty = [candidate for candidate in candidates if candidate]
    return max(nonempty, key=len) if nonempty else []


def extract_explicit_defect_count(report: dict[str, Any]) -> int | None:
    candidates: set[int] = set()
    for value in recursive_values(
        report,
        {"defect_count", "event_count", "n_defects"},
    ):
        if isinstance(value, bool):
            continue
        if isinstance(value, int):
            candidates.add(value)
        elif isinstance(value, float) and value.is_integer():
            candidates.add(int(value))
    if not candidates:
        return None
    if len(candidates) != 1:
        raise ValueError(f"Conflicting explicit defect counts: {sorted(candidates)}")
    return next(iter(candidates))


def current_ko_metadata(
    paths: ProjectPaths,
    defect_dir: Path,
    components: list[str],
) -> dict[str, Any]:
    defect_poscar = defect_dir / paths.defect_poscar_name
    required = [paths.bulk_poscar, defect_poscar, paths.finder]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing KO input(s): " + ", ".join(missing))

    return {
        "workflow_version": WORKFLOW_VERSION,
        "defect_name": defect_dir.name,
        "components": components,
        "components_sha256": sha256_json(components),
        "bulk_poscar": str(paths.bulk_poscar.resolve()),
        "bulk_poscar_sha256": sha256_file(paths.bulk_poscar),
        "defect_poscar": str(defect_poscar.resolve()),
        "defect_poscar_sha256": sha256_file(defect_poscar),
        "finder": str(paths.finder.resolve()),
        "finder_sha256": sha256_file(paths.finder),
    }


def compare_ko_metadata(
    report: dict[str, Any],
    expected: dict[str, Any],
) -> list[str]:
    metadata = report.get("_workflow_metadata")
    if not isinstance(metadata, dict):
        return ["KO report has no _workflow_metadata and is treated as stale."]

    keys = (
        "workflow_version",
        "defect_name",
        "components",
        "components_sha256",
        "bulk_poscar_sha256",
        "defect_poscar_sha256",
        "finder_sha256",
    )
    differences: list[str] = []
    for key in keys:
        if metadata.get(key) != expected.get(key):
            differences.append(
                f"{key}: report={metadata.get(key)!r}, current={expected.get(key)!r}"
            )
    return differences


def load_json_object(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"JSON root must be an object: {path}")
    return data


# -----------------------------------------------------------------------------
# Defect analysis and generated entries
# -----------------------------------------------------------------------------



def correction_center_mode_from_config(config: dict[str, Any]) -> str:
    """
    Return the per-defect correction-centre mode.

    Supported values:
      auto   : run Finder and require a valid KO/FNV centre report.
      skip   : do not run Finder; still generate the defect entry, but omit correction_center.

    Legacy alias:
      skip_correction_center = true
    """
    if config.get("skip_correction_center", False) is True:
        return "skip"

    mode = config.get("correction_center_mode", "auto")
    if not isinstance(mode, str):
        raise ValueError("correction_center_mode must be 'auto' or 'skip'.")
    mode = mode.strip().lower()
    if mode not in {"auto", "skip"}:
        raise ValueError(
            f"Unsupported correction_center_mode={mode!r}; use 'auto' or 'skip'."
        )
    return mode


def analyze_defect(
    paths: ProjectPaths,
    defect_dir: Path,
    config: dict[str, Any],
    strict_outcar: bool,
    require_current_ko: bool,
    require_audit: bool,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    name = defect_dir.name
    item: dict[str, Any] = {
        "defect_name": name,
        "defect_directory": safe_relative(defect_dir, paths.root),
        "status": "pending",
        "warnings": [],
        "unresolved_reasons": [],
    }

    cfg = config.get(name)
    if not isinstance(cfg, dict):
        item["status"] = "unresolved"
        item["unresolved_reasons"].append(
            "Missing defect_components.toml table. Run scan first."
        )
        return None, item

    if cfg.get("enabled", True) is False:
        item["status"] = "disabled"
        return None, item

    components = cfg.get("expected_components", cfg.get("components"))
    label = cfg.get("label", "")
    if not isinstance(label, str):
        item["status"] = "unresolved"
        item["unresolved_reasons"].append("label must be a string.")
        return None, item
    if not isinstance(components, list) or not components or not all(
        isinstance(component, str) and component.strip() for component in components
    ):
        item["status"] = "unresolved"
        item["unresolved_reasons"].append(
            "components is empty or invalid in defect_components.toml."
        )
        return None, item
    components = [component.strip() for component in components]
    item["components"] = components
    item["components_source"] = "defect_components.toml"
    item["label"] = label

    try:
        center_mode = correction_center_mode_from_config(cfg)
    except Exception as exc:
        item["status"] = "unresolved"
        item["unresolved_reasons"].append(str(exc))
        return None, item
    item["correction_center_mode"] = center_mode

    if require_audit:
        audit_ok, audit_reasons = audit_report_is_current_pass(
            paths, defect_dir, components
        )
        item["audit_report"] = safe_relative(
            paths.audit_dir / f"{name}.json", paths.root
        )
        item["audit_reasons"] = audit_reasons
        if not audit_ok:
            item["status"] = "unresolved"
            item["unresolved_reasons"].append(
                "A current PASS structure audit is required: " + " | ".join(audit_reasons)
            )
            return None, item

    try:
        parsed_components = [parse_component(component) for component in components]
        groups = load_species_groups(paths.species_groups_toml)
        expected_delta = expected_atomic_delta(components, groups)
        multiplicity = total_component_multiplicity(components)
        item["parsed_components"] = parsed_components
        item["component_multiplicity"] = multiplicity
        item["expected_composition_delta_defect_minus_bulk"] = expected_delta
    except Exception as exc:
        item["status"] = "unresolved"
        item["unresolved_reasons"].append(str(exc))
        return None, item

    defect_poscar = defect_dir / paths.defect_poscar_name
    item["bulk_poscar"] = safe_relative(paths.bulk_poscar, paths.root)
    item["defect_poscar"] = safe_relative(defect_poscar, paths.root)

    try:
        bulk_composition = parse_poscar_composition(paths.bulk_poscar)
        defect_composition = parse_poscar_composition(defect_poscar)
        actual_delta = composition_delta(bulk_composition, defect_composition)
        item["bulk_composition"] = bulk_composition
        item["defect_composition"] = defect_composition
        item["actual_composition_delta_defect_minus_bulk"] = actual_delta
        if actual_delta != expected_delta:
            raise ValueError(
                "Configured components do not match POSCAR composition change: "
                f"expected {expected_delta}, actual {actual_delta}"
            )
    except Exception as exc:
        item["status"] = "unresolved"
        item["unresolved_reasons"].append(str(exc))
        return None, item

    try:
        charge_map = charge_directories(defect_dir)
        if 0 not in charge_map:
            raise ValueError("Neutral charge_state_0 directory is missing.")
        charges = sort_charges(charge_map)
        item["charges"] = charges
        item["charge_directories"] = {
            str(charge): safe_relative(charge_map[charge], paths.root)
            for charge in charges
        }

        outcar_hashes: dict[str, str] = {}
        for charge in charges:
            outcar = outcar_path_for_charge(
                defect_dir,
                charge,
                paths.charge_outcar_template,
            )
            errors, warnings = validate_outcar(outcar, strict_outcar)
            item["warnings"].extend(warnings)
            if errors:
                item["unresolved_reasons"].extend(errors)
            elif outcar.is_file():
                outcar_hashes[str(charge)] = sha256_file(outcar)
        item["outcar_sha256_by_charge"] = outcar_hashes
        if item["unresolved_reasons"]:
            item["status"] = "unresolved"
            return None, item
    except Exception as exc:
        item["status"] = "unresolved"
        item["unresolved_reasons"].append(f"Charge-state validation failed: {exc}")
        return None, item

    neutral_incar = defect_dir / paths.neutral_incar_rel
    item["neutral_incar"] = safe_relative(neutral_incar, paths.root)
    nelect, raw_nelect_lines, nelect_warnings = parse_nelect(neutral_incar)
    item["raw_nelect_lines"] = raw_nelect_lines
    item["neutral_nelect"] = nelect
    item["warnings"].extend(nelect_warnings)
    if nelect is None:
        item["status"] = "unresolved"
        item["unresolved_reasons"].append("No explicit neutral-state NELECT.")
        return None, item

    if center_mode == "skip":
        center: list[float] | None = None
        defect_count = multiplicity
        item["ko_report"] = None
        item["correction_center"] = None
        item["detected_events"] = []
        item["defect_count_sources"] = {
            "components": multiplicity,
            "events": None,
            "explicit": None,
        }
        item["defect_count"] = defect_count
        item["warnings"].append(
            "KO/FNV correction-centre search skipped by defect configuration."
        )
    else:
        ko_report_path = paths.reports_dir / f"{name}.json"
        item["ko_report"] = safe_relative(ko_report_path, paths.root)
        if not ko_report_path.is_file():
            item["status"] = "unresolved"
            item["unresolved_reasons"].append(
                "KO report is missing. Run the centers command."
            )
            return None, item

        try:
            ko_report = load_json_object(ko_report_path)
            expected_metadata = current_ko_metadata(paths, defect_dir, components)
            stale_reasons = compare_ko_metadata(ko_report, expected_metadata)
            item["ko_report_stale_reasons"] = stale_reasons
            if stale_reasons and require_current_ko:
                raise ValueError("KO report is stale: " + " | ".join(stale_reasons))

            center = extract_correction_center(ko_report)
            events = extract_events(ko_report)
            explicit_count = extract_explicit_defect_count(ko_report)
            event_count = len(events) if events else None
            count_sources = {
                "components": multiplicity,
                "events": event_count,
                "explicit": explicit_count,
            }
            counts = {value for value in count_sources.values() if value is not None}
            if len(counts) != 1:
                raise ValueError(f"Inconsistent defect counts: {count_sources}")
            defect_count = next(iter(counts))
            if defect_count <= 0:
                raise ValueError(f"Invalid defect count: {defect_count}")

            item["correction_center"] = center
            item["detected_events"] = events
            item["defect_count_sources"] = count_sources
            item["defect_count"] = defect_count
        except Exception as exc:
            item["status"] = "unresolved"
            item["unresolved_reasons"].append(str(exc))
            return None, item

    entry = build_defect_entry(
        components=components,
        label=label,
        charges=charges,
        center=center,
        defect_count=defect_count,
        nelect=nelect,
        correction_center_mode=center_mode,
    )
    item["generated_entry"] = entry
    item["defect_type"] = entry["type"]
    item["status"] = "generated"
    return entry, item


def build_defect_entry(
    *,
    components: list[str],
    label: str,
    charges: list[int],
    center: list[float] | None,
    defect_count: int,
    nelect: int | float,
    correction_center_mode: str = "auto",
) -> dict[str, Any]:
    parsed = [parse_component(component) for component in components]
    simple = len(parsed) == 1 and parsed[0]["multiplicity"] == 1

    correction_center = None
    if correction_center_mode == "auto":
        if center is None:
            raise ValueError("Automatic correction-centre mode requires a centre.")
        correction_center = {
            "kind": "fractional_coordinates",
            "coords": center,
        }

    concentration = {
        "defect": defect_count,
        "electron": nelect,
        "hole": nelect,
    }

    def finalize(entry: dict[str, Any]) -> dict[str, Any]:
        # The skip mode is a workflow/configuration instruction only.
        # Skipped defects are still generated normally, but the final database
        # entry omits the correction_center key entirely.
        if correction_center_mode == "auto":
            entry["correction_center"] = correction_center
        entry["conc"] = concentration
        return entry

    if simple:
        component = parsed[0]
        if component["kind"] == "substitution":
            return finalize({
                "type": "substitution",
                "charges": charges,
                "target": component["host_species"],
                "sub_species": component["new_species"],
            })
        if component["kind"] == "vacancy":
            return finalize({
                "type": "vacancy",
                "charges": charges,
                "target": component["species"],
            })
        if component["kind"] == "interstitial":
            return finalize({
                "type": "interstitial",
                "charges": charges,
                "target": component["species"],
            })

    return finalize({
        "type": "cluster_legacy",
        "label": label.strip() or " + ".join(components),
        "components": components,
        "charges": charges,
    })


def analyze_project(
    paths: ProjectPaths,
    strict_outcar: bool,
    require_audit: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    if not paths.components_toml.is_file():
        raise FileNotFoundError(
            f"Missing {paths.components_toml}. Run init and scan first."
        )
    if not paths.finder.is_file():
        raise FileNotFoundError(f"Finder script not found: {paths.finder}")

    config = load_toml(paths.components_toml)
    defects = discover_defects(paths)
    generated: dict[str, Any] = {}
    report: dict[str, Any] = {
        "workflow_version": WORKFLOW_VERSION,
        "root": str(paths.root),
        "generated_at": now_iso(),
        "defects": [],
    }

    for defect_dir in defects:
        entry, item = analyze_defect(
            paths,
            defect_dir,
            config,
            strict_outcar,
            require_current_ko=True,
            require_audit=require_audit,
        )
        report["defects"].append(item)
        if entry is not None:
            generated[defect_dir.name] = entry

    statuses = [item["status"] for item in report["defects"]]
    report.update(
        {
            "discovered_count": len(defects),
            "generated_count": statuses.count("generated"),
            "unresolved_count": statuses.count("unresolved"),
            "disabled_count": statuses.count("disabled"),
        }
    )
    return generated, report


def render_generated_python(generated: dict[str, Any]) -> str:
    lines = [
        "# -*- coding: utf-8 -*-",
        "# Auto-generated defect configurations.",
        "# Review before copying or importing into a Spinney script.",
        "",
        "GENERATED_DEFECT_DB = {",
    ]
    for defect_name in sorted(generated):
        entry_text = pprint.pformat(
            generated[defect_name],
            sort_dicts=False,
            width=100,
            compact=False,
        )
        entry_lines = entry_text.splitlines()
        lines.append(f"    {defect_name!r}: {entry_lines[0]}")
        lines.extend(f"    {line}" for line in entry_lines[1:])
        lines[-1] += ","
        lines.append("")
    lines.append("}")
    lines.append("")
    return "\n".join(lines)


def import_generated_database(path: Path) -> dict[str, Any]:
    spec = importlib.util.spec_from_file_location(
        "_generated_defect_database_validation",
        path,
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot create import spec for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    database = getattr(module, "GENERATED_DEFECT_DB", None)
    if not isinstance(database, dict):
        raise TypeError("GENERATED_DEFECT_DB is missing or is not a dict.")
    return database


# -----------------------------------------------------------------------------
# Commands
# -----------------------------------------------------------------------------


def command_init(paths: ProjectPaths, args: argparse.Namespace) -> int:
    initialize_directories(paths)
    result = {
        "workflow_version": WORKFLOW_VERSION,
        "root": str(paths.root),
        "output_directory": str(paths.postprocess),
        "created_at": now_iso(),
    }
    atomic_write_json(paths.postprocess / "init_report.json", result)
    print(f"initialized={paths.postprocess}")
    return 0


def command_scan(paths: ProjectPaths, args: argparse.Namespace) -> int:
    initialize_directories(paths)

    core_errors: list[str] = []
    if not paths.bulk_poscar.is_file():
        core_errors.append(f"Bulk POSCAR missing: {paths.bulk_poscar}")

    defects = discover_defects(paths)
    sync_result = sync_components_template(paths, defects)

    report: dict[str, Any] = {
        "workflow_version": WORKFLOW_VERSION,
        "root": str(paths.root),
        "generated_at": now_iso(),
        "core_errors": core_errors,
        "added_to_components_toml": sync_result["added"],
        "defects": [],
    }

    for defect_dir in defects:
        item: dict[str, Any] = {
            "defect_name": defect_dir.name,
            "defect_poscar": safe_relative(
                defect_dir / paths.defect_poscar_name,
                paths.root,
            ),
            "errors": [],
            "warnings": [],
        }
        try:
            charges = sort_charges(charge_directories(defect_dir))
            item["charges"] = charges
            if 0 not in charges:
                item["errors"].append("charge_state_0 is missing.")
            for charge in charges:
                outcar = outcar_path_for_charge(
                    defect_dir,
                    charge,
                    paths.charge_outcar_template,
                )
                errors, warnings = validate_outcar(outcar, args.strict_outcar)
                item["errors"].extend(errors)
                item["warnings"].extend(warnings)

            neutral_incar = defect_dir / paths.neutral_incar_rel
            nelect, raw_lines, nelect_warnings = parse_nelect(neutral_incar)
            item["neutral_incar"] = safe_relative(neutral_incar, paths.root)
            item["neutral_nelect"] = nelect
            item["raw_nelect_lines"] = raw_lines
            item["warnings"].extend(nelect_warnings)
            if nelect is None:
                item["errors"].append("No explicit neutral-state NELECT.")
        except Exception as exc:
            item["errors"].append(str(exc))

        item["status"] = "ready" if not item["errors"] else "unresolved"
        report["defects"].append(item)

    report["discovered_count"] = len(defects)
    report["unresolved_count"] = sum(
        item["status"] == "unresolved" for item in report["defects"]
    )
    atomic_write_json(paths.scan_report, report)

    print(
        f"discovered={report['discovered_count']} "
        f"unresolved={report['unresolved_count']} "
        f"added={len(sync_result['added'])}"
    )
    print(f"components={paths.components_toml}")
    print(f"report={paths.scan_report}")

    if core_errors:
        return 3
    return 2 if report["unresolved_count"] else 0



def command_audit(paths: ProjectPaths, args: argparse.Namespace) -> int:
    """Independently identify constructed defects before VASP calculations."""
    initialize_directories(paths)
    if not paths.bulk_poscar.is_file():
        print(f"ERROR: bulk POSCAR missing: {paths.bulk_poscar}", file=sys.stderr)
        return 3

    structure_defects = discover_structure_defects(paths)
    sync_result = sync_components_template(paths, structure_defects)
    config = load_toml(paths.components_toml)
    groups = load_species_groups(paths.species_groups_toml)

    selected = {path.name: path for path in structure_defects}
    selected_names = sorted(selected)
    if args.only:
        if args.only not in selected:
            print(
                f"ERROR: defect initial POSCAR not discovered at project root: {args.only}",
                file=sys.stderr,
            )
            return 3
        selected_names = [args.only]

    try:
        bulk_structure = parse_poscar_structure(paths.bulk_poscar)
        bulk_composition = parse_poscar_composition(paths.bulk_poscar)
        bulk_sites, bulk_group_warnings = build_audit_sites(
            bulk_structure,
            groups,
            args.group_ambiguity_tolerance,
        )
    except Exception as exc:
        print(f"ERROR: failed to prepare bulk reference: {exc}", file=sys.stderr)
        return 4

    summary: dict[str, Any] = {
        "workflow_version": WORKFLOW_VERSION,
        "root": str(paths.root),
        "output_directory": str(paths.postprocess),
        "generated_at": now_iso(),
        "bulk_poscar": safe_relative(paths.bulk_poscar, paths.root),
        "bulk_composition": bulk_composition,
        "species_groups": groups,
        "new_template_entries": sync_result["added"],
        "defects": [],
    }

    pass_count = 0
    fail_count = 0
    disabled_count = 0
    unresolved_count = 0

    for name in selected_names:
        defect_dir = selected[name]
        defect_poscar = defect_dir / paths.defect_poscar_name
        report_path = paths.audit_dir / f"{name}.json"
        item: dict[str, Any] = {
            "defect_name": name,
            "bulk_poscar": safe_relative(paths.bulk_poscar, paths.root),
            "defect_poscar": safe_relative(defect_poscar, paths.root),
            "status": "PENDING",
            "warnings": list(bulk_group_warnings),
            "errors": [],
        }

        cfg = config.get(name)
        if not isinstance(cfg, dict):
            item["status"] = "UNRESOLVED_CONFIGURATION"
            item["errors"].append("Missing defect configuration table.")
            unresolved_count += 1
            atomic_write_json(report_path, item)
            summary["defects"].append(item)
            continue
        if cfg.get("enabled", True) is False:
            item["status"] = "DISABLED"
            disabled_count += 1
            atomic_write_json(report_path, item)
            summary["defects"].append(item)
            continue

        expected_raw = cfg.get("expected_components", cfg.get("components"))
        if not isinstance(expected_raw, list) or not expected_raw or not all(
            isinstance(component, str) and component.strip()
            for component in expected_raw
        ):
            item["status"] = "UNRESOLVED_CONFIGURATION"
            item["errors"].append(
                "expected_components/components is empty or invalid."
            )
            unresolved_count += 1
            atomic_write_json(report_path, item)
            summary["defects"].append(item)
            continue

        expected_components = [component.strip() for component in expected_raw]
        item["expected_components"] = expected_components
        item["expected_components_source"] = (
            "expected_components"
            if "expected_components" in cfg
            else "components"
        )

        try:
            expected_expanded = expand_components(expected_components)
            expected_delta = expected_atomic_delta(expected_components, groups)
            defect_structure = parse_poscar_structure(defect_poscar)
            defect_composition = parse_poscar_composition(defect_poscar)
            actual_delta = composition_delta(bulk_composition, defect_composition)
            lattice_difference = lattice_max_difference(
                bulk_structure.lattice,
                defect_structure.lattice,
            )
            minimum_distance, minimum_pair = minimum_pair_distance(defect_structure)

            defect_sites, defect_group_warnings = build_audit_sites(
                defect_structure,
                groups,
                args.group_ambiguity_tolerance,
            )
            item["warnings"].extend(defect_group_warnings)

            match_result = match_audit_sites(
                bulk_sites,
                defect_sites,
                bulk_structure.lattice,
                args.site_tolerance,
                args.ambiguity_tolerance,
            )
            detected_components, events = independently_detect_events(
                bulk_structure,
                defect_structure,
                bulk_sites,
                defect_sites,
                match_result,
            )
            comparison = compare_component_multisets(
                expected_expanded,
                detected_components,
            )

            item.update(
                {
                    "label": cfg.get("label", ""),
                    "lattice_max_abs_difference_angstrom": lattice_difference,
                    "lattice_tolerance_angstrom": args.lattice_tolerance,
                    "bulk_site_count_after_grouping": len(bulk_sites),
                    "defect_site_count_after_grouping": len(defect_sites),
                    "bulk_composition": bulk_composition,
                    "defect_composition": defect_composition,
                    "expected_atomic_composition_delta_defect_minus_bulk": expected_delta,
                    "actual_atomic_composition_delta_defect_minus_bulk": actual_delta,
                    "expected_components_expanded": expected_expanded,
                    "detected_components": detected_components,
                    "events": events,
                    "missing_events": comparison["missing"],
                    "unexpected_events": comparison["unexpected"],
                    "mapping_ambiguities": match_result["ambiguities"],
                    "minimum_interatomic_distance_angstrom": minimum_distance,
                    "minimum_distance_atom_pair_0based": minimum_pair,
                    "site_tolerance_angstrom": args.site_tolerance,
                    "ambiguity_tolerance_angstrom": args.ambiguity_tolerance,
                    "minimum_allowed_distance_angstrom": args.min_distance,
                }
            )

            if lattice_difference > args.lattice_tolerance:
                status = "FAIL_LATTICE_MISMATCH"
            elif minimum_distance is not None and minimum_distance < args.min_distance:
                status = "FAIL_ATOM_OVERLAP"
            elif actual_delta != expected_delta:
                status = "FAIL_COMPOSITION_MISMATCH"
            elif match_result["ambiguities"] and not args.allow_ambiguous:
                status = "AMBIGUOUS_MAPPING"
            elif comparison["match"]:
                status = "PASS"
            else:
                expected_counter = Counter(expected_expanded)
                detected_counter = Counter(detected_components)
                if set(expected_counter) == set(detected_counter):
                    status = "FAIL_COUNT_MISMATCH"
                elif comparison["missing"] and comparison["unexpected"]:
                    status = "FAIL_TYPE_MISMATCH"
                elif comparison["missing"]:
                    status = "FAIL_MISSING_DEFECT"
                else:
                    status = "FAIL_UNEXPECTED_DEFECT"

            if args.strict_audit and item["warnings"] and status == "PASS":
                status = "AMBIGUOUS_MAPPING"
                item["errors"].append(
                    "Strict audit converted grouping/mapping warnings into a failure."
                )

            item["status"] = status
            item["_workflow_metadata"] = current_audit_metadata(
                paths,
                defect_dir,
                expected_components,
            )
            item["_workflow_metadata"].update(
                {
                    "generated_at": now_iso(),
                    "site_tolerance": args.site_tolerance,
                    "ambiguity_tolerance": args.ambiguity_tolerance,
                    "group_ambiguity_tolerance": args.group_ambiguity_tolerance,
                    "lattice_tolerance": args.lattice_tolerance,
                    "min_distance": args.min_distance,
                }
            )
        except Exception as exc:
            item["status"] = "AUDIT_ERROR"
            item["errors"].append(str(exc))

        if item["status"] == "PASS":
            pass_count += 1
        elif item["status"] in {"UNRESOLVED_CONFIGURATION", "AUDIT_ERROR"}:
            unresolved_count += 1
        else:
            fail_count += 1

        if report_path.exists():
            backup_file(report_path, paths.backups_dir)
        atomic_write_json(report_path, item)
        summary["defects"].append(item)

    summary.update(
        {
            "selected_count": len(selected_names),
            "pass_count": pass_count,
            "fail_count": fail_count,
            "unresolved_count": unresolved_count,
            "disabled_count": disabled_count,
        }
    )
    atomic_write_json(paths.audit_summary, summary)
    print(
        f"selected={len(selected_names)} pass={pass_count} fail={fail_count} "
        f"unresolved={unresolved_count} disabled={disabled_count}"
    )
    print(f"summary={paths.audit_summary}")
    if unresolved_count:
        return 2
    if fail_count:
        return 2
    return 0


def command_centers(paths: ProjectPaths, args: argparse.Namespace) -> int:
    initialize_directories(paths)

    if not paths.finder.is_file():
        print(f"ERROR: finder missing: {paths.finder}", file=sys.stderr)
        return 3
    if not paths.bulk_poscar.is_file():
        print(f"ERROR: bulk POSCAR missing: {paths.bulk_poscar}", file=sys.stderr)
        return 3
    if not paths.components_toml.is_file():
        print(
            f"ERROR: components TOML missing: {paths.components_toml}",
            file=sys.stderr,
        )
        return 3

    config = load_toml(paths.components_toml)
    defects = {path.name: path for path in discover_defects(paths)}
    selected_names = sorted(defects)
    if args.only:
        if args.only not in defects:
            print(
                f"ERROR: defect not discovered at project root: {args.only}",
                file=sys.stderr,
            )
            return 3
        selected_names = [args.only]

    summary: dict[str, Any] = {
        "workflow_version": WORKFLOW_VERSION,
        "root": str(paths.root),
        "generated_at": now_iso(),
        "check_only": args.check_only,
        "force": args.force,
        "defects": [],
    }

    generated_count = 0
    current_count = 0
    skipped_count = 0
    unresolved_count = 0
    failed_count = 0

    for name in selected_names:
        defect_dir = defects[name]
        item: dict[str, Any] = {
            "defect_name": name,
            "status": "pending",
            "warnings": [],
            "errors": [],
        }
        summary["defects"].append(item)

        cfg = config.get(name)
        if not isinstance(cfg, dict):
            item["status"] = "unresolved"
            item["errors"].append("Missing defect_components.toml table.")
            unresolved_count += 1
            continue
        if cfg.get("enabled", True) is False:
            item["status"] = "disabled"
            continue

        components = cfg.get("expected_components", cfg.get("components"))
        if not isinstance(components, list) or not components or not all(
            isinstance(component, str) and component.strip()
            for component in components
        ):
            item["status"] = "unresolved"
            item["errors"].append("components is empty or invalid.")
            unresolved_count += 1
            continue
        components = [component.strip() for component in components]
        item["components"] = components

        try:
            center_mode = correction_center_mode_from_config(cfg)
        except Exception as exc:
            item["status"] = "unresolved"
            item["errors"].append(str(exc))
            unresolved_count += 1
            continue
        item["correction_center_mode"] = center_mode

        if args.require_audit:
            audit_ok, audit_reasons = audit_report_is_current_pass(
                paths, defect_dir, components
            )
            item["audit_report"] = safe_relative(
                paths.audit_dir / f"{name}.json", paths.root
            )
            item["audit_reasons"] = audit_reasons
            if not audit_ok:
                item["status"] = "unresolved"
                item["errors"].append(
                    "A current PASS structure audit is required: "
                    + " | ".join(audit_reasons)
                )
                unresolved_count += 1
                continue

        # Validate composition before calling the finder.
        try:
            groups = load_species_groups(paths.species_groups_toml)
            expected_delta = expected_atomic_delta(components, groups)
            bulk_composition = parse_poscar_composition(paths.bulk_poscar)
            defect_composition = parse_poscar_composition(
                defect_dir / paths.defect_poscar_name
            )
            actual_delta = composition_delta(bulk_composition, defect_composition)
            item["expected_composition_delta"] = expected_delta
            item["actual_composition_delta"] = actual_delta
            if expected_delta != actual_delta:
                raise ValueError(
                    f"components/POSCAR mismatch: expected {expected_delta}, "
                    f"actual {actual_delta}"
                )
            expected_metadata = current_ko_metadata(paths, defect_dir, components)
        except Exception as exc:
            item["status"] = "unresolved"
            item["errors"].append(str(exc))
            unresolved_count += 1
            continue

        if center_mode == "skip":
            item["status"] = "skipped"
            item["report"] = None
            item["warnings"].append(
                "KO/FNV correction-centre search skipped by defect configuration."
            )
            skipped_count += 1
            continue

        report_path = paths.reports_dir / f"{name}.json"
        item["report"] = safe_relative(report_path, paths.root)

        stale_reasons: list[str]
        if report_path.is_file():
            try:
                report = load_json_object(report_path)
                stale_reasons = compare_ko_metadata(report, expected_metadata)
            except Exception as exc:
                stale_reasons = [f"Existing report is unreadable: {exc}"]
        else:
            stale_reasons = ["KO report does not exist."]
        item["stale_reasons"] = stale_reasons

        if not stale_reasons and not args.force:
            item["status"] = "current"
            current_count += 1
            continue

        if args.check_only:
            item["status"] = "stale" if report_path.exists() else "missing"
            unresolved_count += 1
            continue

        command = [
            sys.executable,
            "-X",
            "utf8",
            str(paths.finder),
            "--bulk",
            str(paths.bulk_poscar),
            "--defect",
            str(defect_dir / paths.defect_poscar_name),
        ]
        for component in components:
            command.extend(["--component", component])
        command.extend(["--report", str(report_path)])
        item["command"] = command

        stdout_path = paths.logs_dir / f"ko_{name}.stdout.log"
        stderr_path = paths.logs_dir / f"ko_{name}.stderr.log"
        item["stdout_log"] = safe_relative(stdout_path, paths.root)
        item["stderr_log"] = safe_relative(stderr_path, paths.root)

        if report_path.exists():
            backup = backup_file(report_path, paths.backups_dir)
            item["report_backup"] = (
                safe_relative(backup, paths.root) if backup else None
            )

        finder_env = os.environ.copy()
        finder_env["PYTHONUTF8"] = "1"
        finder_env["PYTHONIOENCODING"] = "utf-8"

        process = subprocess.run(
            command,
            cwd=paths.root,
            text=True,
            encoding="utf-8",
            errors="replace",
            capture_output=True,
            shell=False,
            env=finder_env,
        )
        atomic_write_text(stdout_path, process.stdout)
        atomic_write_text(stderr_path, process.stderr)
        item["return_code"] = process.returncode

        if process.returncode != 0:
            item["status"] = "failed"
            item["errors"].append(
                f"Finder returned {process.returncode}; inspect the stderr log."
            )
            failed_count += 1
            continue
        if not report_path.is_file():
            item["status"] = "failed"
            item["errors"].append(
                "Finder returned 0 but did not create the requested JSON report."
            )
            failed_count += 1
            continue

        try:
            report = load_json_object(report_path)
            metadata = dict(expected_metadata)
            metadata.update(
                {
                    "generated_at": now_iso(),
                    "command": command,
                    "return_code": process.returncode,
                    "hostname": socket.gethostname(),
                }
            )
            report["_workflow_metadata"] = metadata
            atomic_write_json(report_path, report)
            item["status"] = "generated"
            generated_count += 1
        except Exception as exc:
            item["status"] = "failed"
            item["errors"].append(f"Could not finalize KO report: {exc}")
            failed_count += 1

    summary.update(
        {
            "selected_count": len(selected_names),
            "generated_count": generated_count,
            "current_count": current_count,
            "skipped_count": skipped_count,
            "unresolved_count": unresolved_count,
            "failed_count": failed_count,
        }
    )
    summary_path = paths.reports_dir / "ko_center_summary.json"
    atomic_write_json(summary_path, summary)

    print(
        f"selected={len(selected_names)} generated={generated_count} "
        f"current={current_count} skipped={skipped_count} "
        f"unresolved={unresolved_count} failed={failed_count}"
    )
    print(f"summary={summary_path}")

    if failed_count:
        return 4
    if unresolved_count:
        return 2
    return 0


def command_check(paths: ProjectPaths, args: argparse.Namespace) -> int:
    initialize_directories(paths)
    try:
        generated, report = analyze_project(paths, args.strict_outcar, args.require_audit)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    except Exception as exc:
        print(f"ERROR: project analysis failed: {exc}", file=sys.stderr)
        return 4

    report["mode"] = "check"
    atomic_write_json(paths.check_report, report)
    print(
        f"discovered={report['discovered_count']} "
        f"generated={report['generated_count']} "
        f"unresolved={report['unresolved_count']} "
        f"disabled={report['disabled_count']}"
    )
    print(f"report={paths.check_report}")
    return 2 if report["unresolved_count"] else 0


def command_build(paths: ProjectPaths, args: argparse.Namespace) -> int:
    initialize_directories(paths)
    try:
        generated, report = analyze_project(paths, args.strict_outcar, args.require_audit)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    except Exception as exc:
        print(f"ERROR: project analysis failed: {exc}", file=sys.stderr)
        return 4

    report["mode"] = "build"
    if report["unresolved_count"]:
        report["py_compile_return_code"] = None
        atomic_write_json(paths.generated_report, report)
        print(
            f"ERROR: unresolved={report['unresolved_count']}; "
            "generated_defect_db.py was not updated.",
            file=sys.stderr,
        )
        return 2

    output_backup = backup_file(paths.generated_python, paths.backups_dir)
    if output_backup:
        report["output_backup"] = safe_relative(output_backup, paths.root)

    atomic_write_text(paths.generated_python, render_generated_python(generated))

    try:
        py_compile.compile(str(paths.generated_python), doraise=True)
        report["py_compile_return_code"] = 0
        imported = import_generated_database(paths.generated_python)
        if imported != generated:
            raise ValueError(
                "Imported GENERATED_DEFECT_DB differs from the in-memory database."
            )
    except Exception as exc:
        report["py_compile_return_code"] = 4
        report["build_error"] = str(exc)
        atomic_write_json(paths.generated_report, report)
        print(f"ERROR: generated Python validation failed: {exc}", file=sys.stderr)
        return 4

    report["generated_python"] = safe_relative(paths.generated_python, paths.root)
    atomic_write_json(paths.generated_report, report)
    print(
        f"generated={report['generated_count']} "
        f"python={paths.generated_python}"
    )
    print(f"report={paths.generated_report}")
    return 0


def command_validate(paths: ProjectPaths, args: argparse.Namespace) -> int:
    initialize_directories(paths)
    if not paths.generated_python.is_file():
        print(
            f"ERROR: generated database does not exist: {paths.generated_python}",
            file=sys.stderr,
        )
        return 3

    try:
        expected, report = analyze_project(paths, args.strict_outcar, args.require_audit)
        actual = import_generated_database(paths.generated_python)
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    except Exception as exc:
        print(f"ERROR: validation setup failed: {exc}", file=sys.stderr)
        return 4

    differences: list[str] = []
    if report["unresolved_count"]:
        differences.append(
            f"Current project has {report['unresolved_count']} unresolved defect(s)."
        )
    if set(actual) != set(expected):
        differences.append(
            f"Key set mismatch: actual={sorted(actual)}, expected={sorted(expected)}"
        )
    for name in sorted(set(actual) & set(expected)):
        if actual[name] != expected[name]:
            differences.append(f"Entry differs from current inputs: {name}")

    validation = {
        "workflow_version": WORKFLOW_VERSION,
        "validated_at": now_iso(),
        "root": str(paths.root),
        "generated_python": str(paths.generated_python),
        "valid": not differences,
        "differences": differences,
        "expected_count": len(expected),
        "actual_count": len(actual),
        "current_analysis": report,
    }
    atomic_write_json(paths.generated_validation, validation)

    if differences:
        print("generated_defect_db.py is stale or invalid:", file=sys.stderr)
        for difference in differences:
            print(f"  - {difference}", file=sys.stderr)
        print(f"report={paths.generated_validation}")
        return 5

    print("generated_defect_db.py is current and valid.")
    print(f"report={paths.generated_validation}")
    return 0


def command_all(paths: ProjectPaths, args: argparse.Namespace) -> int:
    sequence = (
        ("scan", command_scan),
        ("audit", command_audit),
        ("centers", command_centers),
        ("check", command_check),
        ("build", command_build),
        ("validate", command_validate),
    )
    for name, function in sequence:
        print(f"\n===== {name} =====")
        return_code = function(paths, args)
        if return_code != 0:
            print(
                f"Workflow stopped at {name} with return code {return_code}.",
                file=sys.stderr,
            )
            return return_code
    return 0


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="General defect Spinney + KO configuration workflow."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help=(
            "Project root containing supercell/ and defect directories. "
            "Default: current working directory."
        ),
    )
    parser.add_argument(
        "--finder",
        type=Path,
        help=(
            "Path to find_cluster_ko_center_V1.py. "
            "Default: next to this workflow script."
        ),
    )
    parser.add_argument(
        "--bulk-poscar",
        default="supercell/POSCAR",
        help="Bulk POSCAR path relative to project root.",
    )
    parser.add_argument(
        "--bulk-outcar",
        default="supercell/OUTCAR",
        help="Bulk OUTCAR path relative to project root.",
    )
    parser.add_argument(
        "--bulk-vasprun",
        default="supercell/vasprun.xml",
        help="Bulk vasprun.xml path relative to project root.",
    )
    parser.add_argument(
        "--defect-poscar-name",
        default="POSCAR",
        help="Unrelaxed defect POSCAR name/path relative to each defect directory.",
    )
    parser.add_argument(
        "--neutral-incar",
        default="charge_state_0/scf/INCAR",
        help="Neutral INCAR path relative to each defect directory.",
    )
    parser.add_argument(
        "--charge-outcar",
        default="charge_state_{q}/scf/OUTCAR",
        help="Charge OUTCAR template relative to each defect directory.",
    )
    parser.add_argument(
        "--workspace",
        type=Path,
        help=(
            "Deprecated compatibility option. Ignored in v2.2; active configuration "
            "is read from the project root and generated output is written under "
            "--postprocess-dir."
        ),
    )
    parser.add_argument(
        "--components-file",
        default="defect_components.toml",
        help="Active defect-component TOML file relative to the project root.",
    )
    parser.add_argument(
        "--species-groups-file",
        default="species_groups.toml",
        help="Active molecular/pseudo-species TOML file relative to the project root.",
    )
    parser.add_argument(
        "--postprocess-dir",
        default="_defect_postprocess",
        help="Generated-output directory name relative to project root.",
    )
    parser.add_argument(
        "--require-audit",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Require a current PASS structure-audit report before KO/result processing.",
    )
    parser.add_argument(
        "--strict-outcar",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Require TOTEN, normal-end marker and core-potential block in each OUTCAR. "
            "Use --no-strict-outcar to downgrade missing end/core markers to warnings."
        ),
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("init", help="Create the postprocess directory structure.")
    subparsers.add_parser(
        "scan",
        help="Discover result directories, validate basic inputs and sync the TOML template.",
    )

    audit = subparsers.add_parser(
        "audit",
        help="Independently audit unrelaxed defect POSCARs before VASP calculations.",
    )
    audit.add_argument("--only", help="Audit one exact defect directory name.")
    audit.add_argument(
        "--site-tolerance",
        type=float,
        default=0.35,
        help="Maximum periodic site-matching distance in angstrom.",
    )
    audit.add_argument(
        "--ambiguity-tolerance",
        type=float,
        default=0.02,
        help="Distance-gap threshold for ambiguous site mappings in angstrom.",
    )
    audit.add_argument(
        "--group-ambiguity-tolerance",
        type=float,
        default=0.10,
        help="Boundary-gap threshold for molecular-group membership ambiguity.",
    )
    audit.add_argument(
        "--lattice-tolerance",
        type=float,
        default=1.0e-6,
        help="Maximum absolute lattice-vector difference in angstrom.",
    )
    audit.add_argument(
        "--min-distance",
        type=float,
        default=0.50,
        help="Minimum allowed interatomic distance in defect POSCAR, angstrom.",
    )
    audit.add_argument(
        "--allow-ambiguous",
        action="store_true",
        help="Do not fail solely because mapping alternatives are nearly degenerate.",
    )
    audit.add_argument(
        "--strict",
        dest="strict_audit",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Treat molecular-group construction warnings as audit failures.",
    )

    centers = subparsers.add_parser(
        "centers",
        help="Generate or refresh KO-center reports.",
    )
    centers.add_argument("--only", help="Process one exact defect directory name.")
    centers.add_argument(
        "--force",
        action="store_true",
        help="Regenerate current KO reports.",
    )
    centers.add_argument(
        "--check-only",
        action="store_true",
        help="Check KO report freshness without running the finder.",
    )

    subparsers.add_parser(
        "check",
        help="Validate components, compositions, charges, NELECT and KO reports.",
    )
    subparsers.add_parser(
        "build",
        help="Generate _defect_postprocess/generated_defect_db.py.",
    )
    subparsers.add_parser(
        "validate",
        help="Validate the generated database against all current inputs.",
    )

    all_parser = subparsers.add_parser(
        "all",
        help="Run scan -> audit -> centers -> check -> build -> validate.",
    )
    # Keep a compatible attribute set for command_centers when invoked by all.
    all_parser.set_defaults(
        only=None,
        force=False,
        check_only=False,
        site_tolerance=0.35,
        ambiguity_tolerance=0.02,
        group_ambiguity_tolerance=0.10,
        lattice_tolerance=1.0e-6,
        min_distance=0.50,
        allow_ambiguous=False,
        strict_audit=False,
    )

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    paths = build_paths(args)

    if not paths.root.is_dir():
        print(f"ERROR: project root does not exist: {paths.root}", file=sys.stderr)
        return 3

    commands = {
        "init": command_init,
        "scan": command_scan,
        "audit": command_audit,
        "centers": command_centers,
        "check": command_check,
        "build": command_build,
        "validate": command_validate,
        "all": command_all,
    }

    try:
        return commands[args.command](paths, args)
    except KeyboardInterrupt:
        print("Interrupted by user.", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"UNHANDLED ERROR: {exc}", file=sys.stderr)
        return 4


if __name__ == "__main__":
    raise SystemExit(main())
