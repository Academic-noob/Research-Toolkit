#!/usr/bin/env python3
"""
01_generate_doped_snb_structures_v3.py

Config-driven, POTCAR-free structure generation for doped + ShakeNBreak.

Stages performed automatically:
  1. Read a relaxed pristine bulk structure.
  2. Build an efficient doped defect supercell, or preserve an exact input cell.
  3. Generate intrinsic defects and optional extrinsic defects with doped.
  4. Apply manual charge-state edits and filters.
  5. Save inventories and representative defect POSCARs.
  6. Generate ShakeNBreak distorted structures with Distortions.apply_distortions().
  7. Create a ready-to-upload directory tree with POSCAR, metadata and optional
     INCAR / KPOINTS / job templates, without requiring local POTCAR files.

Typical use in PowerShell:
  Copy-Item .\\00_workflow_config.toml.example .\\00_workflow_config.toml
  python .\\01_generate_doped_snb_structures_v3.py --config .\\00_workflow_config.toml
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 fallback
    import tomli as tomllib  # type: ignore

from doped.generation import DefectsGenerator
from monty.json import MontyDecoder
from monty.serialization import dumpfn, loadfn
from pymatgen.io.vasp.inputs import Poscar
from shakenbreak.input import Distortions


@dataclass(frozen=True)
class WorkflowPaths:
    config: Path
    root: Path
    stage_doped: Path
    stage_snb: Path
    stage_reports: Path
    stage_templates: Path


def fail(message: str) -> "NoReturn":  # type: ignore[name-defined]
    raise SystemExit(f"ERROR: {message}")


def resolve(base: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def load_toml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        fail(f"Missing TOML config: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def cfg_table(cfg: dict[str, Any], key: str) -> dict[str, Any]:
    value = cfg.get(key, {})
    if not isinstance(value, dict):
        fail(f"TOML section [{key}] must be a table")
    return value


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.+%-]+", "_", str(value).strip())
    return cleaned.strip("_") or "unnamed"


def charge_label(charge: int | float) -> str:
    integer = int(round(float(charge)))
    if integer > 0:
        return f"+{integer}"
    return str(integer)


def strip_charge_suffix(entry_name: str) -> str:
    return re.sub(r"_[+-]?\d+$", "", entry_name)


def infer_defect_type(entry: Any) -> str:
    defect_type = getattr(getattr(entry, "defect", None), "defect_type", None)
    text = getattr(defect_type, "name", str(defect_type)).lower()
    for name in ("vacancy", "substitution", "interstitial", "other"):
        if name in text:
            return name
    # Fallback for older object representations.
    class_name = getattr(getattr(entry, "defect", None), "__class__", type(None)).__name__.lower()
    for name in ("vacancy", "substitution", "interstitial"):
        if name in class_name:
            return name
    return "unknown"


def frac_coords(entry: Any) -> list[float] | None:
    value = getattr(entry, "sc_defect_frac_coords", None)
    if value is None:
        return None
    return [float(x) for x in value]


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def copy_optional(source: Path | None, destination: Path) -> bool:
    if source is None:
        return False
    if not source.is_file():
        fail(f"Template file not found: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return True


def optional_template(config_dir: Path, text: str | None) -> Path | None:
    if text is None or not str(text).strip():
        return None
    return resolve(config_dir, str(text))


def prepare_paths(cfg: dict[str, Any], config_path: Path) -> tuple[WorkflowPaths, dict[str, Path | None]]:
    paths_cfg = cfg_table(cfg, "paths")
    workflow_cfg = cfg_table(cfg, "workflow")
    config_dir = config_path.parent

    root = resolve(config_dir, str(paths_cfg.get("workflow_root", "doped_snb_workflow")))
    overwrite = bool(workflow_cfg.get("overwrite", False))
    if root.exists() and any(root.iterdir()):
        if not overwrite:
            fail(f"Workflow root exists and is not empty: {root}. Set workflow.overwrite=true only when deletion is intended.")
        shutil.rmtree(root)

    root.mkdir(parents=True, exist_ok=True)
    paths = WorkflowPaths(
        config=config_path,
        root=root,
        stage_doped=root / "01_doped_defects",
        stage_snb=root / "02_snb_structures",
        stage_reports=root / "00_reports",
        stage_templates=root / "00_templates_snapshot",
    )
    for directory in (paths.stage_doped, paths.stage_snb, paths.stage_reports, paths.stage_templates):
        directory.mkdir(parents=True, exist_ok=True)

    templates = {
        "INCAR": optional_template(config_dir, paths_cfg.get("incar_template")),
        "KPOINTS": optional_template(config_dir, paths_cfg.get("kpoints_template")),
        "job.sh": optional_template(config_dir, paths_cfg.get("job_template")),
    }
    for target_name, source in templates.items():
        if source is not None:
            copy_optional(source, paths.stage_templates / f"{target_name}.template")
    shutil.copy2(config_path, paths.root / "00_workflow_config.used.toml")
    return paths, templates


def build_defects_generator(cfg: dict[str, Any], bulk_path: Path) -> DefectsGenerator:
    workflow_cfg = cfg_table(cfg, "workflow")
    supercell_cfg = cfg_table(cfg, "supercell")
    defects_cfg = cfg_table(cfg, "defects")
    interstitial_cfg = cfg_table(cfg, "interstitials")

    generate_supercell = bool(supercell_cfg.get("generate_supercell", True))
    kwargs: dict[str, Any] = {
        "generate_supercell": generate_supercell,
        "processes": int(workflow_cfg.get("processes", 1)),
        "charge_state_gen_kwargs": {
            "padding": int(defects_cfg.get("charge_padding", 1)),
            "probability_threshold": float(defects_cfg.get("charge_probability_threshold", 0.0075)),
        },
    }

    target_frac_coords = supercell_cfg.get("target_frac_coords", [])
    if target_frac_coords:
        if not isinstance(target_frac_coords, list) or len(target_frac_coords) != 3:
            fail("supercell.target_frac_coords must be [] or a three-number array")
        kwargs["target_frac_coords"] = [float(x) for x in target_frac_coords]

    if generate_supercell:
        kwargs["supercell_gen_kwargs"] = {
            "min_image_distance": float(supercell_cfg.get("min_image_distance", 10.0)),
            "min_atoms": int(supercell_cfg.get("min_atoms", 50)),
            "ideal_threshold": float(supercell_cfg.get("ideal_threshold", 0.10)),
            "force_cubic": bool(supercell_cfg.get("force_cubic", False)),
            "force_diagonal": bool(supercell_cfg.get("force_diagonal", False)),
        }

    extrinsic = defects_cfg.get("extrinsic", [])
    if extrinsic:
        kwargs["extrinsic"] = extrinsic

    generate_interstitials = bool(defects_cfg.get("generate_interstitials", True))
    if not generate_interstitials:
        kwargs["interstitial_gen_kwargs"] = False
    else:
        generator_kwargs = cfg_table(interstitial_cfg, "generator_kwargs") if "generator_kwargs" in interstitial_cfg else {}
        if generator_kwargs:
            kwargs["interstitial_gen_kwargs"] = generator_kwargs
        coords = interstitial_cfg.get("coords", [])
        if coords:
            kwargs["interstitial_coords"] = coords

    print("=" * 88)
    print("Stage 1: doped defect generation")
    print("=" * 88)
    print(f"bulk = {bulk_path}")
    print("DefectsGenerator kwargs =")
    print(json.dumps(kwargs, indent=2, ensure_ascii=False, default=str))
    print()
    return DefectsGenerator(str(bulk_path), **kwargs)


def edit_charge_states(defect_gen: DefectsGenerator, cfg: dict[str, Any]) -> None:
    defects_cfg = cfg_table(cfg, "defects")
    add_map = defects_cfg.get("add_charge_states", {})
    remove_map = defects_cfg.get("remove_charge_states", {})
    if not isinstance(add_map, dict) or not isinstance(remove_map, dict):
        fail("[defects.add_charge_states] and [defects.remove_charge_states] must be TOML tables")
    for species, charges in add_map.items():
        defect_gen.add_charge_states(defect_entry_name=str(species), charge_states=[int(x) for x in charges])
    for species, charges in remove_map.items():
        defect_gen.remove_charge_states(defect_entry_name=str(species), charge_states=[int(x) for x in charges])


def regex_match_any(patterns: Iterable[str], text: str) -> bool:
    return any(re.search(pattern, text) for pattern in patterns)


def select_entries(defect_gen: DefectsGenerator, cfg: dict[str, Any]) -> dict[str, Any]:
    defects_cfg = cfg_table(cfg, "defects")
    include_types = {str(x).lower() for x in defects_cfg.get("include_types", ["vacancy", "substitution", "interstitial"])}
    include_patterns = [str(x) for x in defects_cfg.get("include_name_regex", [])]
    exclude_patterns = [str(x) for x in defects_cfg.get("exclude_name_regex", [])]
    charge_filter = {int(x) for x in defects_cfg.get("include_charges", [])}

    selected: dict[str, Any] = {}
    for entry_name, entry in defect_gen.defect_entries.items():
        species = strip_charge_suffix(entry_name)
        defect_type = infer_defect_type(entry)
        charge = int(entry.charge_state)
        if defect_type not in include_types:
            continue
        if include_patterns and not regex_match_any(include_patterns, species) and not regex_match_any(include_patterns, entry_name):
            continue
        if exclude_patterns and (regex_match_any(exclude_patterns, species) or regex_match_any(exclude_patterns, entry_name)):
            continue
        if charge_filter and charge not in charge_filter:
            continue
        selected[entry_name] = entry

    if not selected:
        fail("No defect entries survived filtering. Check [defects] filters and charge-state settings.")
    return selected


def export_doped_stage(paths: WorkflowPaths, defect_gen: DefectsGenerator, selected: dict[str, Any], cfg: dict[str, Any]) -> None:
    workflow_cfg = cfg_table(cfg, "workflow")
    Poscar(defect_gen.bulk_supercell).write_file(paths.stage_doped / "bulk_supercell_POSCAR")
    defect_gen.to_json(str(paths.stage_doped / "defect_gen_all.json.gz"))
    dumpfn(selected, paths.stage_doped / "selected_defect_entries.json.gz")

    rows: list[dict[str, Any]] = []
    representative_written: set[str] = set()
    representative_dir = paths.stage_doped / "representative_defects"
    by_charge_dir = paths.stage_doped / "by_charge"

    for entry_name, entry in selected.items():
        species = strip_charge_suffix(entry_name)
        row = {
            "entry_name": entry_name,
            "defect_species": species,
            "defect_type": infer_defect_type(entry),
            "charge_state": int(entry.charge_state),
            "atom_count": len(entry.defect_supercell),
            "sc_defect_frac_coords": frac_coords(entry),
        }
        rows.append(row)
        entry_dir = by_charge_dir / safe_name(entry_name)
        entry_dir.mkdir(parents=True, exist_ok=True)
        Poscar(entry.defect_supercell).write_file(entry_dir / "POSCAR")
        write_json(entry_dir / "DEFECT_ENTRY_META.json", row)

        if bool(workflow_cfg.get("write_representative_defects", True)) and species not in representative_written:
            species_dir = representative_dir / safe_name(species)
            species_dir.mkdir(parents=True, exist_ok=True)
            Poscar(entry.defect_supercell).write_file(species_dir / "defect_POSCAR")
            write_json(species_dir / "DEFECT_SPECIES_META.json", row)
            representative_written.add(species)

    write_csv(
        paths.stage_reports / "doped_selected_defect_inventory.csv",
        rows,
        ["entry_name", "defect_species", "defect_type", "charge_state", "atom_count", "sc_defect_frac_coords"],
    )


def make_distortions(selected: dict[str, Any], cfg: dict[str, Any]) -> Distortions:
    snb_cfg = cfg_table(cfg, "shakenbreak")
    bond_distortions = list(snb_cfg.get("bond_distortions", []))
    include_dimer = bool(snb_cfg.get("include_dimer", True))
    if bond_distortions:
        bond_distortions = [x for x in bond_distortions if str(x).lower() != "dimer"]
        if include_dimer:
            bond_distortions.append("Dimer")
    elif not include_dimer:
        # Explicit default numeric list without Dimer.
        increment = float(snb_cfg.get("distortion_increment", 0.10))
        if increment <= 0:
            fail("shakenbreak.distortion_increment must be positive")
        count = int(round(0.6 / increment))
        bond_distortions = [round(i * increment, 10) for i in range(-count, count + 1)]

    kwargs: dict[str, Any] = {
        "defect_entries": list(selected.values()),
        "oxidation_states": dict(snb_cfg.get("oxidation_states", {})) or None,
        "distortion_increment": float(snb_cfg.get("distortion_increment", 0.10)),
        "bond_distortions": bond_distortions or None,
        "local_rattle": bool(snb_cfg.get("local_rattle", False)),
        "distorted_elements": dict(snb_cfg.get("distorted_elements", {})) or None,
        "distorted_atoms": dict(snb_cfg.get("distorted_atoms", {})) or None,
    }
    stdev = float(snb_cfg.get("stdev", 0.0))
    d_min = float(snb_cfg.get("d_min", 0.0))
    if stdev > 0:
        kwargs["stdev"] = stdev
    if d_min > 0:
        kwargs["d_min"] = d_min

    print()
    print("=" * 88)
    print("Stage 2: ShakeNBreak structure generation without local POTCAR")
    print("=" * 88)
    printable = {key: value for key, value in kwargs.items() if key != "defect_entries"}
    printable["selected_entry_count"] = len(selected)
    print(json.dumps(printable, indent=2, ensure_ascii=False, default=str))
    print()
    return Distortions(**kwargs)


def extract_snb_structure(
    payload: Any,
    *,
    defect_species: str,
    charge: int,
    distortion_label: str,
) -> tuple[Any, str | None]:
    """Extract a pymatgen Structure from a ShakeNBreak distortion payload.

    ``Distortions.apply_distortions()`` returns a nested dictionary.  Each
    distortion label maps to a metadata dictionary containing at least
    ``"Defect Structure"`` rather than mapping directly to a Structure object.

    This helper also tolerates already-decoded Structure objects and serialized
    Monty dictionaries, so the same exporter works both immediately after
    generation and during recovery from ``snb_defects_dict.json.gz``.
    """
    if hasattr(payload, "is_ordered"):
        return payload, None

    if not isinstance(payload, dict):
        raise TypeError(
            "Unexpected ShakeNBreak payload type for "
            f"{defect_species}, charge={charge}, distortion={distortion_label}: "
            f"{type(payload).__name__}"
        )

    comment = payload.get("POSCAR Comment")
    candidate_keys = (
        "Defect Structure",
        "defect_structure",
        "Structure",
        "structure",
    )

    for key in candidate_keys:
        if key not in payload:
            continue
        candidate = payload[key]

        if hasattr(candidate, "is_ordered"):
            return candidate, str(comment) if comment is not None else None

        if isinstance(candidate, dict):
            decoded = MontyDecoder().process_decoded(candidate)
            if hasattr(decoded, "is_ordered"):
                return decoded, str(comment) if comment is not None else None

    raise TypeError(
        "Could not extract a pymatgen Structure from ShakeNBreak payload for "
        f"{defect_species}, charge={charge}, distortion={distortion_label}. "
        f"Available keys: {sorted(payload.keys())}"
    )


def existing_paths_and_templates(
    cfg: dict[str, Any],
    config_path: Path,
) -> tuple[WorkflowPaths, dict[str, Path | None]]:
    """Resolve an already-generated workflow tree without deleting it."""
    paths_cfg = cfg_table(cfg, "paths")
    config_dir = config_path.parent
    root = resolve(config_dir, str(paths_cfg.get("workflow_root", "doped_snb_workflow")))

    paths = WorkflowPaths(
        config=config_path,
        root=root,
        stage_doped=root / "01_doped_defects",
        stage_snb=root / "02_snb_structures",
        stage_reports=root / "00_reports",
        stage_templates=root / "00_templates_snapshot",
    )
    for directory in (
        paths.root,
        paths.stage_doped,
        paths.stage_snb,
        paths.stage_reports,
        paths.stage_templates,
    ):
        directory.mkdir(parents=True, exist_ok=True)

    templates = {
        "INCAR": optional_template(config_dir, paths_cfg.get("incar_template")),
        "KPOINTS": optional_template(config_dir, paths_cfg.get("kpoints_template")),
        "job.sh": optional_template(config_dir, paths_cfg.get("job_template")),
    }
    return paths, templates


def iter_snb_structures(
    structures_payload: Any,
    *,
    defect_species: str,
    charge: int,
) -> Iterable[tuple[str, Any]]:
    """Yield ``(distortion_label, structure_payload)`` pairs from ShakeNBreak output.

    Current ShakeNBreak ``Distortions.apply_distortions()`` output uses:

    .. code-block:: python

        {
            "Unperturbed": Structure,
            "distortions": {
                "Bond_Distortion_-20.0%": Structure,
                "Bond_Distortion_0.0%": Structure,
                "Bond_Distortion_20.0%": Structure,
                "Dimer": Structure,
                ...
            },
        }

    Older/custom wrappers may already expose a flat dictionary. This iterator
    supports both layouts while refusing ambiguous nested dictionaries.
    """
    if not isinstance(structures_payload, dict):
        raise TypeError(
            "Unexpected ShakeNBreak structures payload type for "
            f"{defect_species}, charge={charge}: "
            f"{type(structures_payload).__name__}"
        )

    yielded: set[str] = set()

    if "Unperturbed" in structures_payload:
        yielded.add("Unperturbed")
        yield "Unperturbed", structures_payload["Unperturbed"]

    nested_distortions = structures_payload.get("distortions")
    if nested_distortions is not None:
        if not isinstance(nested_distortions, dict):
            raise TypeError(
                "ShakeNBreak key 'distortions' must map to a dictionary for "
                f"{defect_species}, charge={charge}; got "
                f"{type(nested_distortions).__name__}"
            )
        for distortion_label, structure_payload in nested_distortions.items():
            label = str(distortion_label)
            yielded.add(label)
            yield label, structure_payload

    # Compatibility path for already-flat dictionaries or custom extensions.
    for distortion_label, structure_payload in structures_payload.items():
        label = str(distortion_label)
        if label in {"Unperturbed", "distortions"} or label in yielded:
            continue
        yield label, structure_payload


def export_snb_stage(
    paths: WorkflowPaths,
    defects_dict: dict[str, Any],
    distortion_metadata: dict[str, Any],
    templates: dict[str, Path | None],
) -> None:
    dumpfn(defects_dict, paths.stage_snb / "snb_defects_dict.json.gz")
    dumpfn(distortion_metadata, paths.stage_snb / "distortion_metadata.json.gz")

    rows: list[dict[str, Any]] = []
    for defect_species, defect_payload in defects_dict.items():
        charges = defect_payload.get("charges", {})
        for charge, charge_payload in charges.items():
            integer_charge = int(charge)
            defect_charge_name = f"{safe_name(defect_species)}_{charge_label(integer_charge)}"
            structures = charge_payload.get("structures", {})
            for distortion_label, structure_payload in iter_snb_structures(
                structures,
                defect_species=str(defect_species),
                charge=integer_charge,
            ):
                calc_dir = paths.stage_snb / defect_charge_name / safe_name(str(distortion_label))
                calc_dir.mkdir(parents=True, exist_ok=True)

                structure, poscar_comment = extract_snb_structure(
                    structure_payload,
                    defect_species=str(defect_species),
                    charge=integer_charge,
                    distortion_label=str(distortion_label),
                )
                if poscar_comment is None:
                    poscar_comment = (
                        f"{defect_species}_{charge_label(integer_charge)} "
                        f"{distortion_label}"
                    )
                Poscar(structure, comment=poscar_comment).write_file(calc_dir / "POSCAR")

                copied = {
                    name: copy_optional(source, calc_dir / name)
                    for name, source in templates.items()
                }
                metadata = {
                    "defect_species": defect_species,
                    "defect_charge_name": defect_charge_name,
                    "charge_state": integer_charge,
                    "distortion_label": str(distortion_label),
                    "calc_dir": str(calc_dir),
                    "poscar_comment": poscar_comment,
                    "snb_payload_keys": (
                        sorted(structure_payload.keys())
                        if isinstance(structure_payload, dict)
                        else []
                    ),
                    "templates_copied": copied,
                    "potcar_status": "PENDING: populate later on cluster with 02_finalize_vasp_inputs.py",
                    "nelect_status": "PENDING: requires final POTCAR ZVAL audit",
                }
                write_json(calc_dir / "CALC_META.json", metadata)
                (calc_dir / "POTCAR_PENDING.txt").write_text(
                    "POTCAR intentionally not generated locally.\n"
                    "Run 02_finalize_vasp_inputs.py on the cluster after uploading this workflow.\n",
                    encoding="utf-8",
                )
                rows.append(metadata)

    write_csv(
        paths.stage_reports / "snb_calculation_inventory.csv",
        rows,
        [
            "defect_species",
            "defect_charge_name",
            "charge_state",
            "distortion_label",
            "calc_dir",
            "poscar_comment",
            "snb_payload_keys",
            "templates_copied",
            "potcar_status",
            "nelect_status",
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate doped defects and ShakeNBreak structures without local POTCAR.")
    parser.add_argument("--config", required=True, help="Path to 00_workflow_config.toml")
    parser.add_argument(
        "--recover-export-only",
        action="store_true",
        help=(
            "Skip doped and ShakeNBreak generation. Re-export POSCAR directories "
            "from an existing 02_snb_structures/snb_defects_dict.json.gz file."
        ),
    )
    args = parser.parse_args()

    config_path = Path(args.config).expanduser().resolve()
    cfg = load_toml(config_path)

    if args.recover_export_only:
        paths, templates = existing_paths_and_templates(cfg, config_path)
        defects_file = paths.stage_snb / "snb_defects_dict.json.gz"
        metadata_file = paths.stage_snb / "distortion_metadata.json.gz"
        if not defects_file.is_file():
            fail(f"Recovery file not found: {defects_file}")
        if not metadata_file.is_file():
            fail(f"Recovery metadata file not found: {metadata_file}")

        print("=" * 88)
        print("Recovery mode: export existing ShakeNBreak structures")
        print("=" * 88)
        print(f"defects_dict = {defects_file}")
        print(f"metadata     = {metadata_file}")
        print()

        defects_dict = loadfn(defects_file)
        distortion_metadata = loadfn(metadata_file)
        export_snb_stage(paths, defects_dict, distortion_metadata, templates)

        summary = {
            "workflow_root": str(paths.root),
            "recovery_mode": True,
            "source": str(defects_file),
            "next_step": (
                "Upload workflow root to cluster and run "
                "02_finalize_vasp_inputs.py --config 00_workflow_config.used.toml"
            ),
        }
        write_json(paths.stage_reports / "generation_summary_recovered.json", summary)
        print("Recovery export complete.")
        print(f"Exported tree: {paths.stage_snb}")
        return

    paths, templates = prepare_paths(cfg, config_path)
    bulk_path = resolve(config_path.parent, str(cfg_table(cfg, "paths").get("bulk", "bulk_POSCAR")))
    if not bulk_path.is_file():
        fail(f"Bulk structure not found: {bulk_path}")

    defect_gen = build_defects_generator(cfg, bulk_path)
    edit_charge_states(defect_gen, cfg)
    selected = select_entries(defect_gen, cfg)
    export_doped_stage(paths, defect_gen, selected, cfg)

    snb_cfg = cfg_table(cfg, "shakenbreak")
    if bool(snb_cfg.get("enabled", True)):
        distortions = make_distortions(selected, cfg)
        defects_dict, distortion_metadata = distortions.apply_distortions()
        export_snb_stage(paths, defects_dict, distortion_metadata, templates)
    else:
        print("ShakeNBreak disabled by configuration; exported doped defects only.")

    summary = {
        "workflow_root": str(paths.root),
        "bulk_input": str(bulk_path),
        "bulk_supercell": str(paths.stage_doped / "bulk_supercell_POSCAR"),
        "selected_defect_entries": len(selected),
        "snb_enabled": bool(snb_cfg.get("enabled", True)),
        "next_step": "Upload workflow root to cluster and run 02_finalize_vasp_inputs.py --config 00_workflow_config.used.toml",
    }
    write_json(paths.stage_reports / "generation_summary.json", summary)
    print()
    print("=" * 88)
    print("Generation complete")
    print("=" * 88)
    for key, value in summary.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
