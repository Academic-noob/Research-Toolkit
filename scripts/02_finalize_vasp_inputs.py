#!/usr/bin/env python3
"""
02_finalize_vasp_inputs.py

Cluster-side finalisation for a POTCAR-free doped + ShakeNBreak structure tree.

This script traverses every calculation directory containing CALC_META.json and:
  - optionally copies cluster INCAR / KPOINTS / job templates;
  - populates POTCAR by copying one supplied file or running a user command;
  - parses POTCAR ZVAL values and POSCAR species counts;
  - writes NELECT = neutral_valence_electrons - defect_charge into INCAR;
  - optionally writes parity-based NUPDOWN as a starting guess;
  - writes a CSV audit report.

Recommended mode when you already have a custom POTCAR assembly script:
  potcar_mode = "command"
  potcar_command = "python /path/to/potcar_5.4.4.py"

Typical use on Linux cluster:
  python 02_finalize_vasp_inputs.py --config doped_snb_workflow/00_workflow_config.used.toml
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 fallback
    import tomli as tomllib  # type: ignore


FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"


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


def parse_poscar_species_counts(path: Path) -> tuple[list[str], list[int]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 7:
        raise ValueError(f"POSCAR too short: {path}")
    symbols = lines[5].split()
    counts = lines[6].split()
    if not symbols or not counts:
        raise ValueError(f"Missing species/count lines in {path}")
    if all(token.lstrip("+-").isdigit() for token in symbols):
        raise ValueError(f"VASP4-style POSCAR without element symbols is not supported: {path}")
    if not all(token.lstrip("+-").isdigit() for token in counts):
        raise ValueError(f"Could not parse POSCAR counts line in {path}: {lines[6]}")
    integer_counts = [int(token) for token in counts]
    if len(symbols) != len(integer_counts):
        raise ValueError(f"Species/count length mismatch in {path}")
    return symbols, integer_counts


def parse_potcar_blocks(path: Path) -> list[dict[str, Any]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    titles = re.findall(r"^\s*TITEL\s*=\s*(.+?)\s*$", text, flags=re.MULTILINE)
    zvals = [float(value) for value in re.findall(rf"^\s*ZVAL\s*=\s*({FLOAT})", text, flags=re.MULTILINE)]
    if not titles or not zvals:
        raise ValueError(f"Could not parse TITEL/ZVAL records from POTCAR: {path}")
    if len(titles) != len(zvals):
        raise ValueError(f"POTCAR TITEL/ZVAL record count mismatch: {path}")
    blocks: list[dict[str, Any]] = []
    for title, zval in zip(titles, zvals):
        # Typical title: PAW_PBE Pb_d 06Sep2000
        tokens = title.split()
        potential_token = tokens[1] if len(tokens) >= 2 else tokens[0]
        base_symbol = potential_token.split("_")[0]
        blocks.append({"title": title, "potential": potential_token, "base_symbol": base_symbol, "zval": zval})
    return blocks


def set_incar_value(path: Path, key: str, value: str | int | float) -> None:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines() if path.is_file() else []
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=", flags=re.IGNORECASE)
    replacement = f"{key} = {value}"
    found = False
    output: list[str] = []
    for line in lines:
        if pattern.search(line) and not line.lstrip().startswith(("#", "!")):
            if not found:
                output.append(replacement)
                found = True
            else:
                output.append(f"# DUPLICATE_DISABLED_BY_FINALIZER: {line}")
        else:
            output.append(line)
    if not found:
        if output and output[-1].strip():
            output.append("")
        output.append(replacement)
    path.write_text("\n".join(output) + "\n", encoding="utf-8")


def copy_template(source_text: str, config_dir: Path, destination: Path) -> bool:
    if not str(source_text).strip():
        return False
    source = resolve(config_dir, source_text)
    if not source.is_file():
        fail(f"Template not found: {source}")
    shutil.copy2(source, destination)
    return True


def run_potcar_action(calc_dir: Path, finalize: dict[str, Any], config_dir: Path) -> str:
    mode = str(finalize.get("potcar_mode", "skip")).strip().lower()
    potcar = calc_dir / "POTCAR"
    protect = bool(finalize.get("protect_existing_potcar", True))
    if potcar.is_file() and protect:
        return "existing POTCAR preserved"
    if mode == "skip":
        return "POTCAR skipped"
    if mode == "copy":
        source_text = str(finalize.get("potcar_file", "")).strip()
        if not source_text:
            raise ValueError("finalize.potcar_file is required when potcar_mode='copy'")
        source = resolve(config_dir, source_text)
        if not source.is_file():
            raise FileNotFoundError(f"Supplied POTCAR not found: {source}")
        shutil.copy2(source, potcar)
        return f"POTCAR copied from {source}"
    if mode == "command":
        command = str(finalize.get("potcar_command", "")).strip()
        if not command:
            raise ValueError("finalize.potcar_command is required when potcar_mode='command'")
        subprocess.run(command, cwd=calc_dir, shell=True, check=True)
        if not potcar.is_file():
            raise FileNotFoundError(f"Custom POTCAR command completed but did not create {potcar}")
        return f"POTCAR generated by command: {command}"
    raise ValueError(f"Unknown finalize.potcar_mode: {mode}")


def update_nelect(calc_dir: Path, charge: int, finalize: dict[str, Any]) -> dict[str, Any]:
    poscar = calc_dir / "POSCAR"
    potcar = calc_dir / "POTCAR"
    incar = calc_dir / "INCAR"
    if not potcar.is_file():
        return {"status": "POTCAR missing; NELECT not written"}
    symbols, counts = parse_poscar_species_counts(poscar)
    blocks = parse_potcar_blocks(potcar)
    if len(symbols) != len(blocks):
        raise ValueError(
            f"POTCAR block count ({len(blocks)}) does not match POSCAR species count ({len(symbols)}) in {calc_dir}"
        )
    mismatches = []
    for poscar_symbol, block in zip(symbols, blocks):
        if poscar_symbol != block["base_symbol"]:
            mismatches.append(f"POSCAR {poscar_symbol} vs POTCAR {block['potential']}")
    neutral = sum(count * block["zval"] for count, block in zip(counts, blocks))
    target = neutral - charge
    if bool(finalize.get("write_nelect", True)):
        set_incar_value(incar, "NELECT", f"{target:.8f}".rstrip("0").rstrip("."))
    suggested_nupdown = int(round(target)) % 2
    mode = str(finalize.get("nupdown_mode", "audit_only")).strip().lower()
    if mode == "parity":
        set_incar_value(incar, "NUPDOWN", suggested_nupdown)
    elif mode not in {"audit_only", "preserve"}:
        raise ValueError(f"Unknown finalize.nupdown_mode: {mode}")
    return {
        "status": "OK",
        "poscar_symbols": symbols,
        "poscar_counts": counts,
        "potcar_blocks": blocks,
        "species_order_mismatches": mismatches,
        "neutral_valence_electrons": neutral,
        "charge_state": charge,
        "target_nelect": target,
        "suggested_parity_nupdown": suggested_nupdown,
        "nupdown_mode": mode,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Populate POTCAR / NELECT / templates for SnB calculation directories.")
    parser.add_argument("--config", required=True, help="Path to 00_workflow_config.used.toml or equivalent")
    args = parser.parse_args()

    config_path = Path(args.config).expanduser().resolve()
    cfg = load_toml(config_path)
    config_dir = config_path.parent
    paths_cfg = cfg_table(cfg, "paths")
    finalize = cfg_table(cfg, "finalize")
    workflow_root = resolve(config_dir, str(paths_cfg.get("workflow_root", ".")))
    # When config was copied into workflow_root, workflow_root may resolve to a nested path.
    # Prefer the config directory itself if it already contains 02_snb_structures.
    if (config_dir / "02_snb_structures").is_dir():
        workflow_root = config_dir
    snb_root = workflow_root / "02_snb_structures"
    reports = workflow_root / "00_reports"
    reports.mkdir(parents=True, exist_ok=True)
    if not snb_root.is_dir():
        fail(f"Missing SnB structure directory: {snb_root}")

    rows: list[dict[str, Any]] = []
    calc_dirs = sorted(path.parent for path in snb_root.rglob("CALC_META.json"))
    if not calc_dirs:
        fail(f"No CALC_META.json files found under {snb_root}")

    for calc_dir in calc_dirs:
        meta = json.loads((calc_dir / "CALC_META.json").read_text(encoding="utf-8"))
        charge = int(meta["charge_state"])
        row: dict[str, Any] = {
            "calc_dir": str(calc_dir),
            "defect_charge_name": meta.get("defect_charge_name", ""),
            "distortion_label": meta.get("distortion_label", ""),
            "charge_state": charge,
            "status": "OK",
        }
        try:
            copied = {
                "INCAR": copy_template(str(finalize.get("cluster_incar_template", "")), config_dir, calc_dir / "INCAR"),
                "KPOINTS": copy_template(str(finalize.get("cluster_kpoints_template", "")), config_dir, calc_dir / "KPOINTS"),
                "job.sh": copy_template(str(finalize.get("cluster_job_template", "")), config_dir, calc_dir / "job.sh"),
            }
            row["cluster_templates_copied"] = copied
            row["potcar_action"] = run_potcar_action(calc_dir, finalize, config_dir)
            audit = update_nelect(calc_dir, charge, finalize)
            row.update(audit)
            (calc_dir / "POTCAR_PENDING.txt").unlink(missing_ok=True)
            (calc_dir / "FINALIZE_AUDIT.json").write_text(json.dumps(row, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
        except Exception as exc:  # keep processing other directories
            row["status"] = f"ERROR: {exc}"
        rows.append(row)

    report = reports / "finalize_vasp_inputs_audit.csv"
    fieldnames = sorted({key for row in rows for key in row})
    with report.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    errors = [row for row in rows if not str(row.get("status", "")).startswith("OK")]
    print("=" * 88)
    print("VASP finalisation complete")
    print("=" * 88)
    print(f"workflow_root = {workflow_root}")
    print(f"calculation_directories = {len(rows)}")
    print(f"errors = {len(errors)}")
    print(f"audit = {report}")
    if errors:
        print("Some directories failed. Inspect the audit CSV before submission.")
        raise SystemExit(2)


if __name__ == "__main__":
    main()
