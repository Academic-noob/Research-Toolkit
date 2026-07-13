#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
find_ko_center_single_or_cluster.py
======================================

为单一缺陷或由 vacancy / interstitial / substitution 任意组合构成的复合缺陷，
自动识别 KO/FNV correction centre，并输出可直接写入 Spinney 配置的：

    'correction_center': {
        'kind': 'fractional_coordinates',
        'coords': [x, y, z],
    },

适用的 component 写法
----------------------
    vac_Cu, 2vac_Cu
    int_O, 2int_Cu
    sub_Zn_Cu, 2sub_Zn_Cu

物理中心定义
------------
对于每个 component，取一个固定的参考点：

1. vacancy V_X
   用“缺失原子在 bulk 中的原始格点坐标”。

2. substitution A_X
   用“被 A 占据的 bulk 宿主 X 格点坐标”。

3. interstitial A_i
   用“代表性已弛豫结构（推荐 q=0）中 A_i 的实际坐标”。

当只识别到一个缺陷参考点时，直接使用该参考点本身；当识别到多个参考点时，
对全部参考点作 minimum-image 周期性平均。此坐标可作为同一名义缺陷（单点或复合）
所有电荷态的固定 KO/FNV correction centre。

单缺陷使用方式
------------
本脚本不要求多个 component。只传入一个 `--component` 即可：

    --component vac_Cu
    --component sub_Zn_Sn
    --component int_O

此时推荐中心就是一个唯一参考点，不做多点几何平均。

重要限制
--------
- 输入 defect 结构应为同一缺陷簇的一个充分弛豫代表构型，推荐 q=0 的
  relax/CONTCAR 或已由其复制得到的 scf/POSCAR。
- 若间隙原子在不同电荷态发生远距离迁移、形成 split-interstitial，或簇在
  电荷态间发生明显重构，则“一个固定点中心”仅是近似。需检查 KO plateau，
  并比较不同合理中心下 correction 的敏感性。
- 本脚本只确定 Spinney 的 correction_center，不改变 Spinney 用实际 OUTCAR
  组分差计算的化学势项，也不设置 VASP 的 NELECT。

依赖
----
    numpy, scipy, pymatgen

示例（PowerShell）
------------------
# 单一替位：Zn_Sn
python .\find_ko_center_single_or_cluster.py `
  --bulk .\supercell\scf\POSCAR `
  --defect .\sub_Zn_Sn\charge_state_0\scf\POSCAR `
  --component sub_Zn_Sn `
  --report .\ko_center_reports\sub_Zn_Sn.json

# 替位 + 空位：V_Cu + Zn_Cu
python .\find_ko_center_single_or_cluster.py `
  --bulk .\supercell\scf\POSCAR `
  --defect .\sub_Zn_Cu_Vac_Cu\charge_state_0\scf\POSCAR `
  --component vac_Cu `
  --component sub_Zn_Cu `
  --report .\ko_center_reports\sub_Zn_Cu_Vac_Cu.json

# 空位 + 间隙：V_S + S_i
python .\find_ko_center_single_or_cluster.py `
  --bulk .\supercell\scf\POSCAR `
  --defect .\vac_S_plus_int_S\charge_state_0\scf\POSCAR `
  --component vac_S `
  --component int_S `
  --report .\ko_center_reports\vac_S_plus_int_S.json

# 两个空位 + 一个间隙：2V_Cu + Zn_i
python .\find_ko_center_single_or_cluster.py `
  --bulk .\supercell\scf\POSCAR `
  --defect .\2vac_Cu_plus_int_Zn\charge_state_0\scf\POSCAR `
  --component 2vac_Cu `
  --component int_Zn `
  --report .\ko_center_reports\2vac_Cu_plus_int_Zn.json
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
from pymatgen.core import Structure
from scipy.optimize import linear_sum_assignment


# component 形式： vac_Cu / 2vac_Cu / int_O / sub_Zn_Cu / 2sub_Zn_Cu
COMPONENT_PATTERN = re.compile(
    r"(?:(\d+))?"
    r"(vac|int|sub)_"
    r"([A-Z][a-z]?)"
    r"(?:_([A-Z][a-z]?))?"
    r"$"
)

BIG_COST = 1.0e8


@dataclass(frozen=True)
class Component:
    """One parsed cluster component."""

    kind: str                 # vac / int / sub
    species: str              # vac/int species, or substitution new species
    host: str | None          # substitution host species
    multiplicity: int
    source_text: str


def parse_components(items: Sequence[str]) -> List[Component]:
    """Parse --component labels into Component objects."""
    parsed: List[Component] = []

    for item in items:
        match = COMPONENT_PATTERN.fullmatch(item)
        if match is None:
            raise ValueError(
                f"Cannot parse component {item!r}.\n"
                "Supported forms: vac_Cu, 2vac_Cu, int_Zn, 2int_O, "
                "sub_Zn_Cu, 2sub_Zn_Cu."
            )

        multiplicity = int(match.group(1) or 1)
        kind = match.group(2)
        species = match.group(3)
        host = match.group(4)

        if kind == "sub" and host is None:
            raise ValueError(f"Substitution component has no host: {item!r}")
        if kind != "sub" and host is not None:
            raise ValueError(f"Invalid non-substitution component: {item!r}")

        parsed.append(
            Component(
                kind=kind,
                species=species,
                host=host,
                multiplicity=multiplicity,
                source_text=item,
            )
        )

    if not parsed:
        raise ValueError("At least one --component is required.")

    return parsed


def flatten_expected_events(
    components: Sequence[Component],
) -> Tuple[Counter[str], Counter[str], Counter[Tuple[str, str]]]:
    """
    Return expected vacancy / interstitial / substitution event counters.

    Returns
    -------
    vacancies : Counter[element]
    interstitials : Counter[element]
    substitutions : Counter[(new_species, host_species)]
    """
    vacancies: Counter[str] = Counter()
    interstitials: Counter[str] = Counter()
    substitutions: Counter[Tuple[str, str]] = Counter()

    for component in components:
        if component.kind == "vac":
            vacancies[component.species] += component.multiplicity
        elif component.kind == "int":
            interstitials[component.species] += component.multiplicity
        elif component.kind == "sub":
            substitutions[(component.species, str(component.host))] += component.multiplicity
        else:
            raise RuntimeError(f"Unexpected component kind: {component.kind}")

    return vacancies, interstitials, substitutions


def element_counts(structure: Structure) -> Counter[str]:
    return Counter(site.specie.symbol for site in structure)


def expected_composition_delta(components: Sequence[Component]) -> Dict[str, int]:
    """Return N_bulk - N_defect expected from component specification."""
    delta: Counter[str] = Counter()

    for component in components:
        if component.kind == "vac":
            delta[component.species] += component.multiplicity
        elif component.kind == "int":
            delta[component.species] -= component.multiplicity
        elif component.kind == "sub":
            delta[component.species] -= component.multiplicity
            delta[str(component.host)] += component.multiplicity

    return {
        element: int(value)
        for element, value in sorted(delta.items())
        if value != 0
    }


def actual_composition_delta(bulk: Structure, defect: Structure) -> Dict[str, int]:
    """Return actual N_bulk - N_defect from structure composition."""
    bulk_counts = element_counts(bulk)
    defect_counts = element_counts(defect)

    return {
        element: int(bulk_counts[element] - defect_counts[element])
        for element in sorted(set(bulk_counts) | set(defect_counts))
        if bulk_counts[element] - defect_counts[element] != 0
    }


def pbc_fractional_delta(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Minimum-image fractional displacement a - b."""
    delta = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return delta - np.round(delta)


def pbc_distance_matrix(bulk: Structure, defect: Structure) -> np.ndarray:
    """
    Return PBC distance matrix [defect_atom, bulk_atom] in Å.

    Distances use the bulk lattice because the requested correction centre is
    ultimately expressed in the bulk reference frame.
    """
    bulk_frac = np.asarray(bulk.frac_coords, dtype=float)
    defect_frac = np.asarray(defect.frac_coords, dtype=float)

    delta = defect_frac[:, None, :] - bulk_frac[None, :, :]
    delta -= np.round(delta)
    cart = np.einsum("...i,ij->...j", delta, bulk.lattice.matrix)
    return np.linalg.norm(cart, axis=2)


def periodic_distance_from_coord(
    frac_coord: Sequence[float],
    structure: Structure,
) -> np.ndarray:
    """PBC distance from frac_coord to every atom in structure, in Å."""
    coord = np.asarray(frac_coord, dtype=float)
    all_frac = np.asarray(structure.frac_coords, dtype=float)
    delta = pbc_fractional_delta(all_frac, coord)
    cart = delta @ structure.lattice.matrix
    return np.linalg.norm(cart, axis=1)


def flatten_counter(counter: Counter[str]) -> List[str]:
    """Expand Counter({'Cu': 2, 'S': 1}) -> ['Cu', 'Cu', 'S']."""
    result: List[str] = []
    for species in sorted(counter):
        result.extend([species] * int(counter[species]))
    return result


def build_constrained_assignment(
    bulk: Structure,
    defect: Structure,
    vacancies: Counter[str],
    interstitials: Counter[str],
    substitutions: Counter[Tuple[str, str]],
    unmatched_cost_ang: float,
) -> Tuple[np.ndarray, np.ndarray, List[str], List[str], np.ndarray]:
    """
    Construct and solve a constrained global assignment with dummy nodes.

    Rows:
      [all actual defect atoms] + [one dummy row for each expected vacancy]

    Columns:
      [all actual bulk atoms] + [one dummy column for each expected interstitial]

    Allowed mappings:
      defect atom -> bulk atom:
          same species (unchanged) OR declared substitution new_species -> host
      defect atom -> interstitial dummy:
          only if its species matches the declared interstitial species
      vacancy dummy -> bulk atom:
          only if its species matches declared vacancy species

    Dummy-dummy mapping is forbidden. Thus the solver must identify each
    declared vacancy and interstitial explicitly.
    """
    vacancy_slots = flatten_counter(vacancies)
    interstitial_slots = flatten_counter(interstitials)

    n_bulk = len(bulk)
    n_defect = len(defect)
    n_vac = len(vacancy_slots)
    n_int = len(interstitial_slots)

    if n_defect + n_vac != n_bulk + n_int:
        raise ValueError(
            "Atom-count relation is inconsistent with components.\n"
            f"bulk atoms={n_bulk}, defect atoms={n_defect}, "
            f"vacancies={n_vac}, interstitials={n_int}.\n"
            "Required: N_defect + N_vac = N_bulk + N_int."
        )

    if unmatched_cost_ang <= 0:
        raise ValueError("--unmatched-cost-ang must be positive.")

    distance_matrix = pbc_distance_matrix(bulk, defect)

    n_rows = n_defect + n_vac
    n_cols = n_bulk + n_int
    cost = np.full((n_rows, n_cols), BIG_COST, dtype=float)

    # Actual defect atom rows -> actual bulk atom columns.
    for defect_index, defect_site in enumerate(defect):
        defect_species = defect_site.specie.symbol

        for bulk_index, bulk_site in enumerate(bulk):
            bulk_species = bulk_site.specie.symbol
            is_unchanged = defect_species == bulk_species
            is_declared_substitution = (defect_species, bulk_species) in substitutions

            if is_unchanged or is_declared_substitution:
                cost[defect_index, bulk_index] = distance_matrix[defect_index, bulk_index]

        # Actual defect atom -> interstitial dummy column.
        for slot_index, int_species in enumerate(interstitial_slots):
            if defect_species == int_species:
                cost[defect_index, n_bulk + slot_index] = unmatched_cost_ang

    # Vacancy dummy rows -> actual bulk atom columns.
    for slot_index, vac_species in enumerate(vacancy_slots):
        row = n_defect + slot_index
        for bulk_index, bulk_site in enumerate(bulk):
            if bulk_site.specie.symbol == vac_species:
                cost[row, bulk_index] = unmatched_cost_ang

    rows, cols = linear_sum_assignment(cost)

    if len(rows) != n_rows or len(cols) != n_cols:
        raise RuntimeError("Assignment did not cover the full square matrix.")

    selected = cost[rows, cols]
    if np.any(selected >= BIG_COST / 2):
        raise ValueError(
            "No chemically allowed global mapping could be constructed.\n"
            "Check --component labels, bulk/defect structures, and whether the "
            "defect structure actually contains the requested vacancy/interstitial/substitution cluster."
        )

    row_to_col = np.full(n_rows, -1, dtype=int)
    row_to_col[rows] = cols

    return row_to_col, distance_matrix, vacancy_slots, interstitial_slots, selected


def collect_events(
    bulk: Structure,
    defect: Structure,
    row_to_col: np.ndarray,
    distance_matrix: np.ndarray,
    vacancy_slots: Sequence[str],
    interstitial_slots: Sequence[str],
    expected_vacancies: Counter[str],
    expected_interstitials: Counter[str],
    expected_substitutions: Counter[Tuple[str, str]],
) -> List[Dict[str, object]]:
    """Decode constrained assignment into vacancy/interstitial/substitution events."""
    n_bulk = len(bulk)
    n_defect = len(defect)

    found_vacancies: Counter[str] = Counter()
    found_interstitials: Counter[str] = Counter()
    found_substitutions: Counter[Tuple[str, str]] = Counter()
    records: List[Dict[str, object]] = []

    # Actual defect atoms.
    for defect_index in range(n_defect):
        col = int(row_to_col[defect_index])
        defect_species = defect[defect_index].specie.symbol

        if col < n_bulk:
            bulk_index = col
            bulk_species = bulk[bulk_index].specie.symbol

            if defect_species == bulk_species:
                continue

            pair = (defect_species, bulk_species)
            if pair not in expected_substitutions:
                raise RuntimeError(
                    "Internal mapping error: unexpected substitution "
                    f"{defect_species}_{bulk_species}."
                )

            found_substitutions[pair] += 1
            records.append(
                {
                    "event_type": "substitution",
                    "event_label": f"{defect_species}_{bulk_species}",
                    "new_species": defect_species,
                    "host_species": bulk_species,
                    "defect_index_0based": int(defect_index),
                    "bulk_host_index_0based": int(bulk_index),
                    "bulk_reference_frac_coords": [
                        float(x) for x in bulk[bulk_index].frac_coords
                    ],
                    "defect_relaxed_frac_coords": [
                        float(x) for x in defect[defect_index].frac_coords
                    ],
                    "mapping_displacement_ang": float(
                        distance_matrix[defect_index, bulk_index]
                    ),
                    "centre_point_source": "bulk_host_site",
                    "centre_point_frac_coords": [
                        float(x) for x in bulk[bulk_index].frac_coords
                    ],
                }
            )

        else:
            int_slot_index = col - n_bulk
            expected_species = interstitial_slots[int_slot_index]

            if defect_species != expected_species:
                raise RuntimeError(
                    "Internal mapping error: interstitial dummy species mismatch."
                )

            found_interstitials[defect_species] += 1
            nearest_distances = periodic_distance_from_coord(
                defect[defect_index].frac_coords,
                bulk,
            )
            nearest_index = int(np.argmin(nearest_distances))

            records.append(
                {
                    "event_type": "interstitial",
                    "event_label": f"{defect_species}_i",
                    "species": defect_species,
                    "defect_index_0based": int(defect_index),
                    "bulk_host_index_0based": None,
                    "bulk_reference_frac_coords": None,
                    "defect_relaxed_frac_coords": [
                        float(x) for x in defect[defect_index].frac_coords
                    ],
                    "nearest_bulk_index_0based": nearest_index,
                    "nearest_bulk_species": bulk[nearest_index].specie.symbol,
                    "nearest_bulk_distance_ang": float(nearest_distances[nearest_index]),
                    "centre_point_source": "representative_relaxed_interstitial_site",
                    "centre_point_frac_coords": [
                        float(x) for x in defect[defect_index].frac_coords
                    ],
                }
            )

    # Vacancy dummy rows.
    for slot_index, vac_species in enumerate(vacancy_slots):
        row = n_defect + slot_index
        col = int(row_to_col[row])

        if col >= n_bulk:
            raise RuntimeError("Internal mapping error: vacancy dummy mapped to non-bulk column.")

        bulk_index = col
        bulk_species = bulk[bulk_index].specie.symbol

        if bulk_species != vac_species:
            raise RuntimeError(
                "Internal mapping error: vacancy dummy species mismatch."
            )

        found_vacancies[vac_species] += 1
        records.append(
            {
                "event_type": "vacancy",
                "event_label": f"V_{vac_species}",
                "species": vac_species,
                "defect_index_0based": None,
                "bulk_host_index_0based": int(bulk_index),
                "bulk_reference_frac_coords": [
                    float(x) for x in bulk[bulk_index].frac_coords
                ],
                "defect_relaxed_frac_coords": None,
                "centre_point_source": "bulk_vacancy_site",
                "centre_point_frac_coords": [
                    float(x) for x in bulk[bulk_index].frac_coords
                ],
            }
        )

    if found_vacancies != expected_vacancies:
        raise ValueError(
            f"Vacancy identification mismatch. Expected={dict(expected_vacancies)}, "
            f"found={dict(found_vacancies)}"
        )
    if found_interstitials != expected_interstitials:
        raise ValueError(
            f"Interstitial identification mismatch. Expected={dict(expected_interstitials)}, "
            f"found={dict(found_interstitials)}"
        )
    if found_substitutions != expected_substitutions:
        expected_text = {
            f"{new}_{host}": int(n)
            for (new, host), n in sorted(expected_substitutions.items())
        }
        found_text = {
            f"{new}_{host}": int(n)
            for (new, host), n in sorted(found_substitutions.items())
        }
        raise ValueError(
            "Substitution identification mismatch.\n"
            f"Expected={expected_text}\nFound={found_text}"
        )

    # Sort output consistently: vacancies, substitutions, then interstitials;
    # within each category by bulk/defect index.
    order = {"vacancy": 0, "substitution": 1, "interstitial": 2}
    records.sort(
        key=lambda r: (
            order[r["event_type"]],
            -1 if r["bulk_host_index_0based"] is None else r["bulk_host_index_0based"],
            -1 if r["defect_index_0based"] is None else r["defect_index_0based"],
        )
    )

    return records


def periodic_centroid(
    frac_coords: Sequence[Sequence[float]],
    lattice_matrix: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute equal-weight minimum-image periodic centroid and point radii in Å.

    A normal arithmetic mean is wrong when coordinates straddle a periodic
    boundary, e.g. x = 0.98 and 0.02. This routine unwraps around the first
    point before averaging, then wraps the result back into [0, 1).
    """
    points = np.asarray(frac_coords, dtype=float)

    if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0:
        raise ValueError("Need at least one N×3 fractional-coordinate array.")

    reference = points[0]
    unwrapped = reference + pbc_fractional_delta(points, reference)
    centre = np.mean(unwrapped, axis=0) % 1.0

    delta = pbc_fractional_delta(points, centre)
    radii = np.linalg.norm(delta @ lattice_matrix, axis=1)
    return centre, radii


def print_event(record: Dict[str, object]) -> None:
    """Pretty-print a decoded defect event."""
    event_type = str(record["event_type"])

    if event_type == "vacancy":
        print(
            f"vacancy      {record['event_label']:<8s}  "
            f"bulk_site={record['bulk_host_index_0based']:>3d}  "
            f"bulk_ref={np.array(record['bulk_reference_frac_coords'])}"
        )

    elif event_type == "substitution":
        print(
            f"substitution {record['event_label']:<8s}  "
            f"bulk_site={record['bulk_host_index_0based']:>3d}  "
            f"defect_site={record['defect_index_0based']:>3d}  "
            f"bulk_ref={np.array(record['bulk_reference_frac_coords'])}  "
            f"relaxed={np.array(record['defect_relaxed_frac_coords'])}  "
            f"disp={record['mapping_displacement_ang']:.4f} Å"
        )

    elif event_type == "interstitial":
        print(
            f"interstitial {record['event_label']:<8s}  "
            f"defect_site={record['defect_index_0based']:>3d}  "
            f"relaxed={np.array(record['defect_relaxed_frac_coords'])}  "
            f"nearest_bulk={record['nearest_bulk_species']}[{record['nearest_bulk_index_0based']}] "
            f"at {record['nearest_bulk_distance_ang']:.4f} Å"
        )

    else:
        raise RuntimeError(f"Unknown event type: {event_type}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Find a fixed KO/FNV correction centre for a single defect or a "
            "vacancy-interstitial-substitution defect cluster."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("--bulk", required=True, help="Bulk POSCAR or CONTCAR.")
    parser.add_argument(
        "--defect",
        required=True,
        help=(
            "A relaxed representative defect structure, preferably q=0 relax/CONTCAR "
            "or a scf/POSCAR copied from it."
        ),
    )
    parser.add_argument(
        "--component",
        action="append",
        required=True,
        help=(
            "Provide one --component for a single defect, or repeat for a cluster. "
            "Examples: vac_Cu, int_O, sub_Zn_Cu, 2vac_Cu, 2sub_Zn_Cu."
        ),
    )
    parser.add_argument(
        "--report",
        default="cluster_ko_center_report.json",
        help="Output JSON audit report. Parent directories are created automatically.",
    )
    parser.add_argument(
        "--unmatched-cost-ang",
        type=float,
        default=0.75,
        help=(
            "Hungarian-assignment penalty used to declare a vacancy/interstitial. "
            "Increase only if normal relaxation is unusually large; lower if a "
            "wrong atom is being treated as unchanged."
        ),
    )
    parser.add_argument(
        "--mapping-warning-ang",
        type=float,
        default=1.50,
        help="Warn when an identified substitution displacement exceeds this value.",
    )
    parser.add_argument(
        "--radius-warning-ang",
        type=float,
        default=5.00,
        help="Warn when maximum cluster-centre radius exceeds this value.",
    )

    return parser


def main() -> None:
    args = build_parser().parse_args()

    bulk_path = Path(args.bulk)
    defect_path = Path(args.defect)
    report_path = Path(args.report)

    if not bulk_path.is_file():
        raise FileNotFoundError(f"Bulk structure not found: {bulk_path}")
    if not defect_path.is_file():
        raise FileNotFoundError(f"Defect structure not found: {defect_path}")

    components = parse_components(args.component)
    expected_vacancies, expected_interstitials, expected_substitutions = (
        flatten_expected_events(components)
    )

    bulk = Structure.from_file(str(bulk_path))
    defect = Structure.from_file(str(defect_path))

    expected_delta = expected_composition_delta(components)
    actual_delta = actual_composition_delta(bulk, defect)

    print("=" * 110)
    print("SINGLE-DEFECT / CLUSTER KO/FNV CORRECTION-CENTRE FINDER")
    print("=" * 110)
    print(f"bulk       : {bulk_path}")
    print(f"defect     : {defect_path}")
    print(f"components : {args.component}")
    print(f"bulk atoms / defect atoms: {len(bulk)} / {len(defect)}")
    print(f"expected N_bulk - N_defect: {expected_delta}")
    print(f"actual   N_bulk - N_defect: {actual_delta}")

    if expected_delta != actual_delta:
        raise ValueError(
            "Composition mismatch. The supplied bulk/defect structures do not "
            "match the requested components.\n"
            f"Expected: {expected_delta}\n"
            f"Actual:   {actual_delta}"
        )

    lattice_rel_diff = float(
        np.linalg.norm(defect.lattice.matrix - bulk.lattice.matrix)
        / np.linalg.norm(bulk.lattice.matrix)
    )
    print(f"bulk/defect lattice relative difference: {lattice_rel_diff:.3e}")
    if lattice_rel_diff > 1.0e-4:
        print(
            "WARNING: bulk and defect lattices differ. KO/FNV point-defect "
            "workflows normally use fixed-cell defect calculations."
        )

    row_to_col, distance_matrix, vacancy_slots, interstitial_slots, selected_costs = (
        build_constrained_assignment(
            bulk=bulk,
            defect=defect,
            vacancies=expected_vacancies,
            interstitials=expected_interstitials,
            substitutions=expected_substitutions,
            unmatched_cost_ang=args.unmatched_cost_ang,
        )
    )

    events = collect_events(
        bulk=bulk,
        defect=defect,
        row_to_col=row_to_col,
        distance_matrix=distance_matrix,
        vacancy_slots=vacancy_slots,
        interstitial_slots=interstitial_slots,
        expected_vacancies=expected_vacancies,
        expected_interstitials=expected_interstitials,
        expected_substitutions=expected_substitutions,
    )

    n_reference_points = len(events)
    defect_scope = "single_point_defect" if n_reference_points == 1 else "defect_cluster"
    print(
        "defect classification: "
        + ("single point defect (one correction reference point)"
           if defect_scope == "single_point_defect"
           else f"defect cluster ({n_reference_points} correction reference points)")
    )

    print("\nIdentified defect events:")
    for event in events:
        print_event(event)

    centre_points = [event["centre_point_frac_coords"] for event in events]
    centre, radii = periodic_centroid(centre_points, bulk.lattice.matrix)

    # Diagnostic centre: use relaxed positions for substitution + interstitial,
    # but retain bulk coordinate for vacancy because no actual atom exists there.
    diagnostic_points: List[Sequence[float]] = []
    for event in events:
        if event["event_type"] == "vacancy":
            diagnostic_points.append(event["bulk_reference_frac_coords"])
        else:
            diagnostic_points.append(event["defect_relaxed_frac_coords"])
    diagnostic_centre, diagnostic_radii = periodic_centroid(
        diagnostic_points,
        bulk.lattice.matrix,
    )

    print("\n" + "=" * 110)
    print("RECOMMENDED FIXED KO/FNV CORRECTION CENTRE")
    print("=" * 110)
    print("Centre-point convention:")
    print("  vacancy       -> original bulk vacancy site")
    print("  substitution  -> original bulk host site")
    print("  interstitial  -> representative relaxed interstitial site")
    print("Use this one centre for every charge state of the same nominal defect.")
    if defect_scope == "single_point_defect":
        print("single-defect reference frac = "
              f"[{centre[0]:.10f}, {centre[1]:.10f}, {centre[2]:.10f}]")
        print("single-defect radius (Å) = 0.0000")
    else:
        print(
            "mixed-reference centroid frac = "
            f"[{centre[0]:.10f}, {centre[1]:.10f}, {centre[2]:.10f}]"
        )
        print(
            "event-point distances from centroid (Å) = "
            + ", ".join(f"{radius:.4f}" for radius in radii)
        )
        print(f"max cluster radius (Å) = {np.max(radii):.4f}")
    print("\nPaste into LEGACY_CLUSTER_DB:")
    print(
        "'correction_center': {'kind': 'fractional_coordinates', "
        f"'coords': [{centre[0]:.10f}, {centre[1]:.10f}, {centre[2]:.10f}]}},"
    )

    print("\nDiagnostic only:")
    if defect_scope == "single_point_defect":
        print(
            "relaxed/reference point frac = "
            f"[{diagnostic_centre[0]:.10f}, {diagnostic_centre[1]:.10f}, {diagnostic_centre[2]:.10f}]"
        )
    else:
        print(
            "relaxed/mixed centroid frac = "
            f"[{diagnostic_centre[0]:.10f}, {diagnostic_centre[1]:.10f}, {diagnostic_centre[2]:.10f}]"
        )
        print(
            "diagnostic event-point distances from centroid (Å) = "
            + ", ".join(f"{radius:.4f}" for radius in diagnostic_radii)
        )

    substitution_displacements = [
        float(event["mapping_displacement_ang"])
        for event in events
        if event["event_type"] == "substitution"
    ]

    if substitution_displacements and max(substitution_displacements) > args.mapping_warning_ang:
        print(
            "\nWARNING: maximum substitution mapping displacement = "
            f"{max(substitution_displacements):.4f} Å > "
            f"{args.mapping_warning_ang:.4f} Å."
        )
        print("Inspect the event table and visually verify the mapping in VESTA.")

    if defect_scope == "defect_cluster" and np.max(radii) > args.radius_warning_ang:
        print(
            "\nWARNING: the cluster is spatially extended; maximum radius = "
            f"{np.max(radii):.4f} Å > {args.radius_warning_ang:.4f} Å."
        )
        print(
            "A single point-charge KO centre may be approximate. Compare KO "
            "alignment plateaus and formation energies for several reasonable centres."
        )

    report = {
        "bulk": str(bulk_path.resolve()),
        "defect": str(defect_path.resolve()),
        "components": list(args.component),
        "defect_classification": defect_scope,
        "n_correction_reference_points": int(n_reference_points),
        "settings": {
            "unmatched_cost_ang": float(args.unmatched_cost_ang),
            "mapping_warning_ang": float(args.mapping_warning_ang),
            "radius_warning_ang": float(args.radius_warning_ang),
        },
        "expected_composition_delta_Nbulk_minus_Ndefect": expected_delta,
        "actual_composition_delta_Nbulk_minus_Ndefect": actual_delta,
        "lattice_relative_difference": lattice_rel_diff,
        "identified_events": events,
        "recommended_correction_center": {
            "kind": "fractional_coordinates",
            "coords": [float(value) for value in centre],
            "definition": (
                "single reference point when one event is identified; otherwise periodic centroid: "
                "bulk vacancy sites + bulk substitution host sites + representative relaxed interstitial sites"
            ),
            "event_point_distances_ang": [float(value) for value in radii],
        },
        "diagnostic_relaxed_mixed_centroid": {
            "coords": [float(value) for value in diagnostic_centre],
            "event_point_distances_ang": [
                float(value) for value in diagnostic_radii
            ],
        },
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"\nJSON report written: {report_path.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        print(f"\nERROR: {error}", file=sys.stderr)
        sys.exit(1)
