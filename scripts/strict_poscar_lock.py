#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Strict coordinate lock for doped + ShakeNBreak workflows.

Purpose
-------
`doped` may reconstruct an input supercell into an equivalent coordinate set.
`pymatgen.io.vasp.inputs.Poscar.write_file()` also serialises floating-point
values again.  Both behaviours are physically acceptable, but they do not
preserve the literal decimal text of the user's input POSCAR.

This module adds an opt-in strict layer:
  * the input POSCAR is treated as the authoritative text template;
  * doped is still used to determine defect topology, symmetry and charges;
  * unperturbed vacancy/substitution POSCARs are rebuilt from the original raw
    POSCAR rows, so every untouched lattice/coordinate line is copied verbatim;
  * ShakeNBreak distorted structures keep their intentional coordinate changes,
    while their scale and lattice-vector lines are restored verbatim from the
    authoritative POSCAR;
  * any unexpected mapping failure aborts the workflow before VASP jobs exist.

The current strict text exporter fully supports VASP5-style POSCAR files with
vacancies and substitutions. Interstitials are supported when Direct
coordinates are used; the newly inserted interstitial row is necessarily new.

Patch version
-------------
v4:
  * Stage-1 unperturbed defect POSCAR export is literal text surgery;
  * all untouched coordinate rows are copied byte-for-byte as strings;
  * title, scale factor, lattice rows, Selective dynamics row and coordinate-mode
    row are preserved verbatim;
  * only the logically necessary header tokens and defect coordinate rows may
    change;
  * internal numerical tolerances are used only for defect-site identification,
    never for writing authoritative Stage-1 coordinates;
  * ShakeNBreak handling remains isolated for later discussion.
"""

from __future__ import annotations

import csv
import hashlib
import math
import re
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Poscar

try:
    from scipy.optimize import linear_sum_assignment
except Exception:  # pragma: no cover - fallback only
    linear_sum_assignment = None


class StrictLockError(RuntimeError):
    """Raised when exact-coordinate locking cannot be guaranteed."""


def _site_symbol(site: Any) -> str:
    specie = getattr(site, "specie", None)
    symbol = getattr(specie, "symbol", None)
    if symbol is not None:
        return str(symbol)
    return str(specie)


def _wrap(frac: Iterable[float]) -> np.ndarray:
    wrapped = np.mod(np.asarray(list(frac), dtype=float), 1.0)
    wrapped[np.isclose(wrapped, 1.0, atol=1e-12)] = 0.0
    wrapped[np.isclose(wrapped, 0.0, atol=1e-12)] = 0.0
    return wrapped


def _pbc_cart_distances(lattice_matrix: np.ndarray, frac: np.ndarray, targets: np.ndarray) -> np.ndarray:
    diff = np.asarray(targets, dtype=float) - np.asarray(frac, dtype=float)
    diff -= np.round(diff)
    cart = diff @ np.asarray(lattice_matrix, dtype=float)
    return np.linalg.norm(cart, axis=1)


def _nearest_site_index(structure: Structure, frac: Iterable[float], symbol: str | None = None) -> tuple[int, float]:
    candidates = [
        idx for idx, site in enumerate(structure)
        if symbol is None or _site_symbol(site) == symbol
    ]
    if not candidates:
        raise StrictLockError(f"No source sites found for element {symbol!r}.")
    target_coords = np.asarray([structure[idx].frac_coords for idx in candidates], dtype=float)
    dists = _pbc_cart_distances(structure.lattice.matrix, _wrap(frac), target_coords)
    local = int(np.argmin(dists))
    return candidates[local], float(dists[local])


@dataclass(frozen=True)
class RawPoscarTemplate:
    """Raw VASP5 POSCAR representation with per-site text rows."""

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
    direct_mode: bool
    coord_lines: list[str]
    site_symbols: list[str]
    trailing_lines: list[str]

    @classmethod
    def from_file(cls, path: Path) -> "RawPoscarTemplate":
        lines = path.read_text(encoding="utf-8", errors="strict").splitlines()
        if len(lines) < 8:
            raise StrictLockError(f"POSCAR is too short: {path}")

        symbols = lines[5].split()
        counts_text = lines[6].split()
        if not symbols or not counts_text:
            raise StrictLockError(f"Missing element/count rows in {path}")
        if all(token.lstrip("+-").isdigit() for token in symbols):
            raise StrictLockError("Strict lock requires a VASP5-style POSCAR with element symbols on line 6.")
        if not all(token.lstrip("+-").isdigit() for token in counts_text):
            raise StrictLockError(f"Could not parse atom counts in {path}: {lines[6]!r}")
        counts = [int(token) for token in counts_text]
        if len(symbols) != len(counts):
            raise StrictLockError("Element/count length mismatch in authoritative POSCAR.")

        cursor = 7
        selective_line: str | None = None
        if lines[cursor].strip().lower().startswith("s"):
            selective_line = lines[cursor]
            cursor += 1
        if cursor >= len(lines):
            raise StrictLockError("Missing coordinate mode row in authoritative POSCAR.")
        coordinate_mode_line = lines[cursor]
        direct_mode = coordinate_mode_line.strip().lower().startswith("d")
        if not direct_mode and not coordinate_mode_line.strip().lower().startswith(("c", "k")):
            raise StrictLockError(f"Unknown POSCAR coordinate mode: {coordinate_mode_line!r}")
        cursor += 1

        natoms = sum(counts)
        coord_lines = lines[cursor : cursor + natoms]
        if len(coord_lines) != natoms:
            raise StrictLockError(f"Expected {natoms} coordinate rows but found {len(coord_lines)} in {path}")
        trailing = lines[cursor + natoms :]
        if any(line.strip() for line in trailing):
            raise StrictLockError(
                "Authoritative POSCAR contains non-empty trailing rows, such as MD velocities. "
                "Strict defect text export refuses to guess how they should be edited."
            )

        expanded_symbols: list[str] = []
        for symbol, count in zip(symbols, counts, strict=True):
            expanded_symbols.extend([symbol] * count)
        return cls(
            path=path,
            lines=lines,
            comment=lines[0],
            scale_line=lines[1],
            lattice_lines=lines[2:5],
            symbols_line=lines[5],
            counts_line=lines[6],
            symbols=symbols,
            counts=counts,
            selective_line=selective_line,
            coordinate_mode_line=coordinate_mode_line,
            direct_mode=direct_mode,
            coord_lines=coord_lines,
            site_symbols=expanded_symbols,
            trailing_lines=trailing,
        )

    @property
    def lattice_block_sha256(self) -> str:
        payload = "\n".join([self.scale_line, *self.lattice_lines]) + "\n"
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()

    @staticmethod
    def _replace_tokens_preserving_spacing(line: str, new_tokens: list[str]) -> str:
        """Replace non-whitespace tokens while preserving existing separators.

        If the number of tokens changes, retain the original leading whitespace
        and use two spaces between tokens. This path is required only when a
        substitution introduces or removes an element species.
        """
        matches = list(re.finditer(r"\S+", line))
        if len(matches) == len(new_tokens):
            parts: list[str] = []
            cursor = 0
            for match, token in zip(matches, new_tokens, strict=True):
                parts.append(line[cursor:match.start()])
                parts.append(str(token))
                cursor = match.end()
            parts.append(line[cursor:])
            return "".join(parts)

        leading = re.match(r"^\s*", line).group(0)  # type: ignore[union-attr]
        trailing = re.search(r"\s*$", line).group(0)  # type: ignore[union-attr]
        return leading + "  ".join(str(token) for token in new_tokens) + trailing

    def render_rows(self, rows: list[tuple[str, str]], *, comment: str | None = None) -> str:
        """Render a Stage-1 defect POSCAR by literal text surgery.

        ``comment`` is intentionally ignored in strict Stage-1 mode. Metadata
        belongs in JSON/CSV reports, not in the authoritative POSCAR text.

        The only permitted edits are:
          * element-count token changes required by a defect;
          * element-symbol token changes if a new substitution species appears;
          * deletion, relocation or insertion of the defect coordinate row.

        Every untouched coordinate row is copied from the original POSCAR string
        without floating-point parsing or formatting.
        """
        del comment  # preserve the authoritative title line exactly

        output_symbols = list(self.symbols)
        for symbol, _ in rows:
            if symbol not in output_symbols:
                output_symbols.append(symbol)

        grouped: dict[str, list[str]] = defaultdict(list)
        for symbol, row in rows:
            grouped[symbol].append(row)

        # Remove species only if the defect makes their count zero. For common
        # vacancies this does not occur; for substitutions it can be necessary.
        output_symbols = [symbol for symbol in output_symbols if grouped.get(symbol)]
        output_counts = [len(grouped[symbol]) for symbol in output_symbols]

        if output_symbols == self.symbols:
            symbols_line = self.symbols_line
        else:
            symbols_line = self._replace_tokens_preserving_spacing(
                self.symbols_line,
                output_symbols,
            )

        counts_line = self._replace_tokens_preserving_spacing(
            self.counts_line,
            [str(count) for count in output_counts],
        )

        output = [
            self.comment,
            self.scale_line,
            *self.lattice_lines,
            symbols_line,
            counts_line,
        ]
        if self.selective_line is not None:
            output.append(self.selective_line)
        output.append(self.coordinate_mode_line)

        for symbol in output_symbols:
            output.extend(grouped[symbol])

        # Preserve any blank trailing rows from the original file. Non-empty
        # trailing rows were already rejected during parsing.
        output.extend(self.trailing_lines)
        return "\n".join(output) + "\n"


@dataclass
class LockedEntry:
    entry_name: str
    defect_species: str
    charge_state: int
    defect_type: str
    operation: str
    source_site_index: int | None
    source_site_symbol: str | None
    inserted_symbol: str | None
    strict_defect_frac_coords: list[float]
    poscar_text: str
    strict_structure: Structure


@dataclass
class StrictLockContext:
    enabled: bool
    bulk_path: Path
    raw_bulk: RawPoscarTemplate | None
    bulk_structure: Structure | None
    internal_lattice_matrix: np.ndarray | None
    generated_to_input_basis_transform: np.ndarray | None
    generated_to_input_shift: np.ndarray | None
    lattice_transform_max_residual: float | None
    mapping_max_distance_angstrom: float | None
    entries_by_name: dict[str, LockedEntry]
    entries_by_species_charge: dict[tuple[str, int], LockedEntry]

    @classmethod
    def disabled(cls, bulk_path: Path) -> "StrictLockContext":
        return cls(
            enabled=False,
            bulk_path=bulk_path,
            raw_bulk=None,
            bulk_structure=None,
            internal_lattice_matrix=None,
            generated_to_input_basis_transform=None,
            generated_to_input_shift=None,
            lattice_transform_max_residual=None,
            mapping_max_distance_angstrom=None,
            entries_by_name={},
            entries_by_species_charge={},
        )

    def write_bulk_exact(self, destination: Path) -> None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if self.enabled:
            shutil.copy2(self.bulk_path, destination)
        else:
            raise StrictLockError("write_bulk_exact() called while strict lock is disabled.")

    def write_entry_exact(self, entry_name: str, destination: Path) -> bool:
        if not self.enabled:
            return False
        record = self.entries_by_name.get(entry_name)
        if record is None:
            raise StrictLockError(f"Missing strict-lock record for {entry_name}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(record.poscar_text, encoding="utf-8")
        return True

    def write_snb_poscar(
        self,
        destination: Path,
        structure: Structure,
        *,
        comment: str,
        defect_species: str,
        charge_state: int,
        distortion_label: str,
    ) -> str:
        """Write SnB POSCAR in the authoritative input-cell basis.

        The unperturbed structure is written from the literal POSCAR text
        template.  Distorted structures are intentionally not text-identical:
        their fractional coordinates are transformed from doped's internal,
        physically equivalent basis back to the supplied POSCAR basis before
        writing.  The raw lattice text block is then restored verbatim.
        """
        destination.parent.mkdir(parents=True, exist_ok=True)
        record = self.entries_by_species_charge.get((str(defect_species), int(charge_state)))
        if self.enabled and str(distortion_label).lower() == "unperturbed":
            if record is None:
                raise StrictLockError(
                    "Strict lock could not find the exact unperturbed POSCAR template for "
                    f"defect={defect_species!r}, charge={charge_state}. "
                    "Refusing to fall back to pymatgen serialisation."
                )
            destination.write_text(record.poscar_text, encoding="utf-8")
            return "STRICT_EXACT_UNPERTURBED"

        if self.enabled:
            converted = self.convert_internal_snb_structure_to_authoritative_basis(structure)
            Poscar(converted, comment=comment).write_file(destination)
            self.restore_exact_lattice_block(destination)
            return "STRICT_AUTHORITATIVE_BASIS_INTENTIONAL_SNB_COORDINATES"

        Poscar(structure, comment=comment).write_file(destination)
        return "STANDARD_PYMATGEN_SERIALISATION"

    def convert_internal_snb_structure_to_authoritative_basis(
        self,
        structure: Structure,
    ) -> Structure:
        """Convert a distorted SnB structure into the supplied POSCAR basis.

        ShakeNBreak must distort the internally consistent doped entry rather
        than a text-reordered strict export, because substitution site indices
        belong to doped's internal structure.  This conversion is applied only
        after distortion generation.
        """
        if (
            not self.enabled
            or self.bulk_structure is None
            or self.internal_lattice_matrix is None
            or self.generated_to_input_basis_transform is None
            or self.generated_to_input_shift is None
        ):
            raise StrictLockError(
                "Strict SnB basis conversion requested without a complete strict-lock context."
            )

        snb_lattice = np.asarray(structure.lattice.matrix, dtype=float)
        if not np.allclose(
            snb_lattice,
            np.asarray(self.internal_lattice_matrix, dtype=float),
            atol=1e-5,
            rtol=0,
        ):
            raise StrictLockError(
                "ShakeNBreak structure lattice changed unexpectedly before export. "
                "Strict mode only supports fixed-cell SnB structure generation.\n"
                f"Expected internal lattice:\n{self.internal_lattice_matrix}\n"
                f"Observed SnB lattice:\n{snb_lattice}"
            )

        converted_frac = np.asarray(structure.frac_coords, dtype=float) @ np.asarray(
            self.generated_to_input_basis_transform,
            dtype=float,
        )
        converted_frac = np.asarray(
            [_wrap(coords + self.generated_to_input_shift) for coords in converted_frac],
            dtype=float,
        )

        return Structure(
            lattice=self.bulk_structure.lattice,
            species=[site.species for site in structure],
            coords=converted_frac,
            coords_are_cartesian=False,
            site_properties=structure.site_properties,
        )

    def restore_exact_lattice_block(self, poscar_path: Path) -> None:
        if not self.enabled or self.raw_bulk is None:
            return
        lines = poscar_path.read_text(encoding="utf-8", errors="strict").splitlines()
        if len(lines) < 5:
            raise StrictLockError(f"Generated POSCAR is too short: {poscar_path}")
        lines[1] = self.raw_bulk.scale_line
        lines[2:5] = self.raw_bulk.lattice_lines
        poscar_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _linear_assignment(cost: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if linear_sum_assignment is not None:
        return linear_sum_assignment(cost)

    # Conservative fallback. This is adequate only when each nearest mapping is
    # unambiguous, which is the expected strict-lock use case.
    pairs = sorted(
        (float(cost[i, j]), i, j)
        for i in range(cost.shape[0])
        for j in range(cost.shape[1])
    )
    used_i: set[int] = set()
    used_j: set[int] = set()
    rows: list[int] = []
    cols: list[int] = []
    for _, i, j in pairs:
        if i in used_i or j in used_j:
            continue
        used_i.add(i)
        used_j.add(j)
        rows.append(i)
        cols.append(j)
        if len(rows) == cost.shape[0]:
            break
    if len(rows) != cost.shape[0]:
        raise StrictLockError("Could not construct a one-to-one site mapping without SciPy.")
    return np.asarray(rows, dtype=int), np.asarray(cols, dtype=int)


def _generated_frac_to_source_basis(
    frac: Iterable[float],
    generated_to_source_basis: np.ndarray,
) -> np.ndarray:
    """Convert generated-cell fractional coordinates into the source-cell basis.

    pymatgen stores lattice vectors as rows, so:
        cart = frac_generated @ lattice_generated
             = frac_generated @ T @ lattice_source

    where:
        T = lattice_generated @ inv(lattice_source)

    For physically equivalent cells with the same Cartesian orientation, T must
    be an integer unimodular matrix, apart from small floating-point noise.
    """
    return np.asarray(list(frac), dtype=float) @ np.asarray(generated_to_source_basis, dtype=float)


def _assignment_for_shift(
    source: Structure,
    generated: Structure,
    shift: np.ndarray,
    generated_to_source_basis: np.ndarray,
) -> tuple[dict[int, int], float, float]:
    lattice = np.asarray(source.lattice.matrix, dtype=float)
    source_by_symbol: dict[str, list[int]] = defaultdict(list)
    generated_by_symbol: dict[str, list[int]] = defaultdict(list)
    for idx, site in enumerate(source):
        source_by_symbol[_site_symbol(site)].append(idx)
    for idx, site in enumerate(generated):
        generated_by_symbol[_site_symbol(site)].append(idx)

    if {k: len(v) for k, v in source_by_symbol.items()} != {k: len(v) for k, v in generated_by_symbol.items()}:
        raise StrictLockError("Input bulk and doped bulk have different element counts.")

    mapping: dict[int, int] = {}
    all_distances: list[float] = []
    for symbol, gen_indices in generated_by_symbol.items():
        src_indices = source_by_symbol[symbol]
        cost = np.zeros((len(gen_indices), len(src_indices)), dtype=float)
        targets = np.asarray([source[src_idx].frac_coords for src_idx in src_indices], dtype=float)
        for i, gen_idx in enumerate(gen_indices):
            converted = _generated_frac_to_source_basis(
                generated[gen_idx].frac_coords,
                generated_to_source_basis,
            )
            shifted = _wrap(converted + shift)
            cost[i, :] = _pbc_cart_distances(lattice, shifted, targets)
        rows, cols = _linear_assignment(cost)
        for row, col in zip(rows, cols, strict=True):
            mapping[gen_indices[int(row)]] = src_indices[int(col)]
            all_distances.append(float(cost[int(row), int(col)]))

    return mapping, max(all_distances, default=0.0), sum(all_distances)


def _find_bulk_alignment(
    source: Structure,
    generated: Structure,
    tolerance_angstrom: float,
    lattice_transform_tolerance: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Map doped's equivalent bulk representation back to the supplied cell.

    The supplied POSCAR remains authoritative. We permit doped's internal bulk
    to use an equivalent lattice basis, including:
      * small floating-point reconstruction noise;
      * lattice-vector permutations;
      * lattice-vector sign flips;
      * integer unimodular basis changes.

    We reject any true supercell-size or strain change.
    """
    if len(source) != len(generated):
        raise StrictLockError(
            f"Strict lock requires input bulk and doped bulk to contain the same number of atoms: "
            f"input={len(source)}, doped={len(generated)}"
        )
    if Counter(_site_symbol(site) for site in source) != Counter(_site_symbol(site) for site in generated):
        raise StrictLockError("Strict lock requires identical input/doped bulk compositions.")

    source_lattice = np.asarray(source.lattice.matrix, dtype=float)
    generated_lattice = np.asarray(generated.lattice.matrix, dtype=float)
    try:
        generated_to_source_basis = generated_lattice @ np.linalg.inv(source_lattice)
    except np.linalg.LinAlgError as exc:
        raise StrictLockError("Supplied POSCAR has a singular lattice matrix.") from exc

    nearest_integer_transform = np.rint(generated_to_source_basis).astype(int)
    lattice_transform_max_residual = float(
        np.max(np.abs(generated_to_source_basis - nearest_integer_transform))
    )
    determinant = int(round(float(np.linalg.det(nearest_integer_transform))))

    if abs(determinant) != 1 or lattice_transform_max_residual > lattice_transform_tolerance:
        raise StrictLockError(
            "Strict lock rejected doped's internal bulk lattice because it is not an equivalent "
            "integer-unimodular basis representation of the supplied supercell.\n"
            f"Maximum transform residual = {lattice_transform_max_residual:.6e}; "
            f"allowed = {lattice_transform_tolerance:.6e}\n"
            f"Rounded basis transform determinant = {determinant}\n"
            f"Generated-to-input basis transform:\n{generated_to_source_basis}\n"
            f"Nearest integer transform:\n{nearest_integer_transform}\n"
            "This usually means that the input structure was not preserved as the same physical "
            "supercell. Check [supercell] generate_supercell=false and use a clean output directory."
        )

    basis_transform = nearest_integer_transform.astype(float)

    source_symbols = [_site_symbol(site) for site in source]
    generated_symbols = [_site_symbol(site) for site in generated]
    counts = Counter(generated_symbols)
    anchor_symbol = min(counts, key=lambda key: (counts[key], key))
    gen_anchor = generated_symbols.index(anchor_symbol)
    converted_anchor = _generated_frac_to_source_basis(
        generated[gen_anchor].frac_coords,
        basis_transform,
    )
    candidates = [
        _wrap(np.asarray(source[src_idx].frac_coords) - converted_anchor)
        for src_idx, symbol in enumerate(source_symbols)
        if symbol == anchor_symbol
    ]

    best_shift: np.ndarray | None = None
    best_max = math.inf
    best_sum = math.inf
    for shift in candidates:
        _, max_dist, total_dist = _assignment_for_shift(
            source,
            generated,
            shift,
            basis_transform,
        )
        if (max_dist, total_dist) < (best_max, best_sum):
            best_shift = shift
            best_max = max_dist
            best_sum = total_dist

    if best_shift is None or best_max > tolerance_angstrom:
        raise StrictLockError(
            "Could not map doped's equivalent bulk coordinates back to the supplied POSCAR within the "
            f"strict tolerance. Best maximum deviation = {best_max:.6e} Å; "
            f"tolerance = {tolerance_angstrom:.6e} Å.\n"
            f"Accepted integer basis transform:\n{basis_transform}"
        )
    return basis_transform, best_shift, float(best_max), lattice_transform_max_residual

def _generated_site_symbol_at_defect(entry: Any, defect_frac_coords: np.ndarray, tolerance_angstrom: float) -> str:
    structure = entry.defect_supercell
    idx, distance = _nearest_site_index(structure, defect_frac_coords)
    if distance > tolerance_angstrom:
        raise StrictLockError(
            f"Could not identify inserted/replacement atom at the defect site: closest distance = {distance:.6e} Å"
        )
    return _site_symbol(structure[idx])


def _format_added_direct_row(frac: Iterable[float], selective: bool) -> str:
    coords = _wrap(frac)
    row = "  " + "  ".join(f"{value:.16f}" for value in coords)
    if selective:
        row += "  T  T  T"
    return row


def _preserve_internal_entry_for_snb(entry: Any) -> None:
    """Validate that the original doped entry remains usable by ShakeNBreak.

    Do not mutate ``entry.sc_entry.structure``. In current pymatgen versions,
    ``ComputedStructureEntry.structure`` is a read-only property backed by
    ``_structure``. More importantly, mutating the entry here can break the site
    index convention used by ShakeNBreak for substitution defects. Strict text
    locking is therefore an export-layer operation only.
    """
    sc_entry = getattr(entry, "sc_entry", None)
    if sc_entry is None or not hasattr(sc_entry, "structure"):
        raise StrictLockError("DefectEntry has no usable sc_entry.structure attribute.")


def apply_strict_coordinate_lock(
    *,
    cfg: dict[str, Any],
    bulk_path: Path,
    defect_gen: Any,
    selected: dict[str, Any],
    report_path: Path,
) -> StrictLockContext:
    """Rebuild selected unperturbed defect structures from the literal input POSCAR."""
    strict_cfg = cfg.get("strict_lock", {})
    if strict_cfg is None:
        strict_cfg = {}
    if not isinstance(strict_cfg, dict):
        raise StrictLockError("TOML section [strict_lock] must be a table.")
    enabled = bool(strict_cfg.get("enabled", False))
    if not enabled:
        return StrictLockContext.disabled(bulk_path)

    supercell_cfg = cfg.get("supercell", {})
    if bool(supercell_cfg.get("generate_supercell", True)):
        raise StrictLockError(
            "[strict_lock] enabled=true is incompatible with [supercell] generate_supercell=true. "
            "A supplied target supercell must be locked with generate_supercell=false."
        )

    tolerance = float(strict_cfg.get("mapping_tolerance_angstrom", 1e-4))
    if tolerance <= 0:
        raise StrictLockError("strict_lock.mapping_tolerance_angstrom must be positive.")

    lattice_transform_tolerance = float(
        strict_cfg.get("lattice_transform_tolerance", 1e-5)
    )
    if lattice_transform_tolerance <= 0:
        raise StrictLockError("strict_lock.lattice_transform_tolerance must be positive.")

    raw = RawPoscarTemplate.from_file(bulk_path)
    exact_bulk = Structure.from_file(str(bulk_path))
    if len(exact_bulk) != len(raw.coord_lines):
        raise StrictLockError("pymatgen atom count does not match authoritative POSCAR coordinate-row count.")

    basis_transform, shift, mapping_max_distance, lattice_transform_max_residual = _find_bulk_alignment(
        exact_bulk,
        defect_gen.bulk_supercell,
        tolerance,
        lattice_transform_tolerance,
    )
    source_rows = list(zip(raw.site_symbols, raw.coord_lines, strict=True))
    entries_by_name: dict[str, LockedEntry] = {}
    entries_by_species_charge: dict[tuple[str, int], LockedEntry] = {}
    report_rows: list[dict[str, Any]] = []

    for entry_name, entry in selected.items():
        defect_species = entry_name.rsplit("_", 1)[0] if entry_name.rsplit("_", 1)[-1].lstrip("+-").isdigit() else entry_name
        charge_state = int(entry.charge_state)
        defect_type_obj = getattr(getattr(entry, "defect", None), "defect_type", None)
        defect_type = getattr(defect_type_obj, "name", str(defect_type_obj)).lower()
        if "vacancy" in defect_type:
            defect_type = "vacancy"
        elif "substitution" in defect_type:
            defect_type = "substitution"
        elif "interstitial" in defect_type:
            defect_type = "interstitial"
        else:
            raise StrictLockError(f"Unsupported strict-lock defect type for {entry_name}: {defect_type}")

        generated_defect_frac = _wrap(entry.sc_defect_frac_coords)
        strict_frac = _wrap(
            _generated_frac_to_source_basis(generated_defect_frac, basis_transform) + shift
        )
        strict_structure = exact_bulk.copy()
        rows = list(source_rows)
        source_site_index: int | None = None
        source_site_symbol: str | None = None
        inserted_symbol: str | None = None

        if defect_type in {"vacancy", "substitution"}:
            source_site_index, distance = _nearest_site_index(exact_bulk, strict_frac)
            if distance > tolerance:
                raise StrictLockError(
                    f"{entry_name}: mapped defect site is not an input bulk site. "
                    f"Distance={distance:.6e} Å; tolerance={tolerance:.6e} Å"
                )
            source_site_symbol = _site_symbol(exact_bulk[source_site_index])

        if defect_type == "vacancy":
            strict_structure.remove_sites([source_site_index])  # type: ignore[list-item]
            rows.pop(source_site_index)  # type: ignore[arg-type]
            operation = f"remove source site #{source_site_index} ({source_site_symbol})"

        elif defect_type == "substitution":
            inserted_symbol = _generated_site_symbol_at_defect(entry, generated_defect_frac, tolerance)
            strict_structure.replace(source_site_index, inserted_symbol)  # type: ignore[arg-type]
            _, raw_row = rows[source_site_index]  # type: ignore[index]
            rows[source_site_index] = (inserted_symbol, raw_row)  # type: ignore[index]
            operation = f"replace source site #{source_site_index} ({source_site_symbol} -> {inserted_symbol})"

        else:  # interstitial
            if not raw.direct_mode:
                raise StrictLockError(
                    f"{entry_name}: strict interstitial text export currently requires Direct coordinates."
                )
            inserted_symbol = _generated_site_symbol_at_defect(entry, generated_defect_frac, tolerance)
            strict_structure.append(inserted_symbol, strict_frac, coords_are_cartesian=False)
            rows.append((inserted_symbol, _format_added_direct_row(strict_frac, raw.selective_line is not None)))
            operation = f"append interstitial ({inserted_symbol})"

        _preserve_internal_entry_for_snb(entry)
        poscar_text = raw.render_rows(rows)
        record = LockedEntry(
            entry_name=entry_name,
            defect_species=defect_species,
            charge_state=charge_state,
            defect_type=defect_type,
            operation=operation,
            source_site_index=source_site_index,
            source_site_symbol=source_site_symbol,
            inserted_symbol=inserted_symbol,
            strict_defect_frac_coords=[float(x) for x in strict_frac],
            poscar_text=poscar_text,
            strict_structure=strict_structure,
        )
        entries_by_name[entry_name] = record
        key = (defect_species, charge_state)
        existing = entries_by_species_charge.get(key)
        if existing is not None and existing.poscar_text != poscar_text:
            raise StrictLockError(f"Duplicate strict records disagree for defect/charge key {key!r}")
        entries_by_species_charge[key] = record
        report_rows.append(
            {
                "entry_name": entry_name,
                "defect_species": defect_species,
                "charge_state": charge_state,
                "defect_type": defect_type,
                "operation": operation,
                "source_site_index": source_site_index,
                "source_site_symbol": source_site_symbol,
                "inserted_symbol": inserted_symbol,
                "strict_defect_frac_coords": [float(x) for x in strict_frac],
                "bulk_mapping_max_distance_angstrom": mapping_max_distance,
                "lattice_transform_max_residual": lattice_transform_max_residual,
                "generated_to_input_basis_transform": basis_transform.tolist(),
                "generated_to_input_shift": shift.tolist(),
                "lattice_block_sha256": raw.lattice_block_sha256,
                "stage1_export_mode": "LITERAL_TEXT_SURGERY",
            }
        )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in report_rows for key in row})
    with report_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report_rows)

    return StrictLockContext(
        enabled=True,
        bulk_path=bulk_path,
        raw_bulk=raw,
        bulk_structure=exact_bulk,
        internal_lattice_matrix=np.asarray(defect_gen.bulk_supercell.lattice.matrix, dtype=float),
        generated_to_input_basis_transform=basis_transform,
        generated_to_input_shift=shift,
        lattice_transform_max_residual=lattice_transform_max_residual,
        mapping_max_distance_angstrom=mapping_max_distance,
        entries_by_name=entries_by_name,
        entries_by_species_charge=entries_by_species_charge,
    )
