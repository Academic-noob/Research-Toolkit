#!/usr/bin/env python
"""Extract PES data with either CLI arguments or an internal config block."""
from pathlib import Path
from glob import glob
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pymatgen.core import Structure
from nonrad.ccd import get_dQ, get_PES_from_vaspruns, get_omega_from_PES

# ============================================================================
# USER-EDITABLE INTERNAL MODE
# ============================================================================
CONFIG_MODE = "cli"  # "cli" or "internal"
INTERNAL_CONFIG = {
    "case_root": "/absolute/path/to/nonrad_case",
    "dE": None,       # eV; fill this before using internal mode
    "ref_index": 4,
}
# ============================================================================

CLI_DEFAULTS = {
    "case_root": "nonrad_case",
    "ref_index": 4,
}


def pick(args, key, default=None):
    cli_value = getattr(args, key)
    if cli_value is not None:
        return cli_value
    if args.mode == "internal":
        return INTERNAL_CONFIG.get(key, default)
    return default


def require(parser, values, keys):
    missing = [key for key in keys if values.get(key) in (None, "")]
    if missing:
        parser.error(
            "Missing required setting(s): " + ", ".join(missing)
            + ". Supply CLI options or fill INTERNAL_CONFIG and use --mode internal."
        )


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract PES and effective phonon frequencies from CCD static vasprun.xml files."
    )
    p.add_argument("--mode", choices=["cli", "internal"], default=CONFIG_MODE)
    p.add_argument("--case-root", default=None)
    p.add_argument("--dE", type=float, default=None,
                   help="Energy offset between excited and ground PES in eV.")
    p.add_argument("--ref-index", type=int, default=None)
    return p, p.parse_args()


def main():
    parser, args = parse_args()
    cfg = {
        "case_root": pick(args, "case_root", CLI_DEFAULTS["case_root"]),
        "dE": pick(args, "dE"),
        "ref_index": pick(args, "ref_index", CLI_DEFAULTS["ref_index"]),
    }
    require(parser, cfg, ["case_root", "dE"])

    case = Path(cfg["case_root"]).expanduser()
    results = case / "05_results"
    results.mkdir(parents=True, exist_ok=True)

    ground_contcar = case / "01_relax" / "ground" / "CONTCAR"
    excited_contcar = case / "01_relax" / "excited" / "CONTCAR"
    for path in (ground_contcar, excited_contcar):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty endpoint CONTCAR: {path}")

    ground_struct = Structure.from_file(ground_contcar)
    excited_struct = Structure.from_file(excited_contcar)
    dQ = get_dQ(ground_struct, excited_struct)

    ground_vaspruns = glob(str(case / "03_static" / "ground" / "*" / "vasprun.xml"))
    excited_vaspruns = glob(str(case / "03_static" / "excited" / "*" / "vasprun.xml"))

    if not ground_vaspruns:
        raise RuntimeError("No ground vasprun.xml files found.")
    if not excited_vaspruns:
        raise RuntimeError("No excited vasprun.xml files found.")

    Q_ground, E_ground = get_PES_from_vaspruns(ground_struct, excited_struct, ground_vaspruns)
    Q_excited, E_excited = get_PES_from_vaspruns(ground_struct, excited_struct, excited_vaspruns)
    E_excited = cfg["dE"] + E_excited

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(Q_ground, E_ground, s=14, label="ground")
    ax.scatter(Q_excited, E_excited, s=14, label="excited")

    q = np.linspace(
        min(min(Q_ground), min(Q_excited)) - 1,
        max(max(Q_ground), max(Q_excited)) + 1.5,
        300,
    )
    ground_omega = get_omega_from_PES(Q_ground, E_ground, ax=ax, q=q)
    excited_omega = get_omega_from_PES(Q_excited, E_excited, ax=ax, q=q)

    ax.set_xlim(-3, 4)
    ax.set_ylim(-0.2, 3.5)
    ax.set_xlabel(r"$Q$ [amu$^{1/2}$ $\AA$]")
    ax.set_ylabel("E [eV]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(results / "ccd_pes.png", dpi=600)
    plt.close(fig)

    out = results / "pes_summary.txt"
    with out.open("w", encoding="utf-8") as f:
        f.write(f"mode = {args.mode}\n")
        f.write(f"case_root = {case}\n")
        f.write(f"ref_index = {cfg['ref_index']}\n")
        f.write(f"dQ = {dQ:.8f} amu^(1/2) Angstrom\n")
        f.write(f"dE = {cfg['dE']:.8f} eV\n")
        f.write(f"ground_omega = {ground_omega:.8f} eV\n")
        f.write(f"excited_omega = {excited_omega:.8f} eV\n")
        f.write("\n===== Ground PES =====\n")
        for qv, ev in sorted(zip(Q_ground, E_ground)):
            f.write(f"Q = {qv:14.8f}   E = {ev:16.10f}\n")
        f.write("\n===== Excited PES =====\n")
        for qv, ev in sorted(zip(Q_excited, E_excited)):
            f.write(f"Q = {qv:14.8f}   E = {ev:16.10f}\n")

    print(out.read_text(encoding="utf-8"))
    print(f"Wrote figure: {results / 'ccd_pes.png'}")


if __name__ == "__main__":
    main()
