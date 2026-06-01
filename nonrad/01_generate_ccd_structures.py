#!/usr/bin/env python
"""Generate CCD structures with either CLI arguments or an internal config block."""
from pathlib import Path
import argparse
import numpy as np
from pymatgen.core import Structure
from nonrad.ccd import get_dQ, get_cc_structures

# ============================================================================
# USER-EDITABLE INTERNAL MODE
# Set CONFIG_MODE = "internal" to run without writing paths in the command line.
# Keep CONFIG_MODE = "cli" to preserve the original command-line usage.
# Any explicitly supplied CLI option overrides the value below.
# ============================================================================
CONFIG_MODE = "cli"  # "cli" or "internal"
INTERNAL_CONFIG = {
    "ground_contcar": "/absolute/path/to/nonrad_case/01_relax/ground/CONTCAR",
    "excited_contcar": "/absolute/path/to/nonrad_case/01_relax/excited/CONTCAR",
    "outdir": "/absolute/path/to/nonrad_case/02_ccd_structures",
    "npoints": 9,
    "qmin": -0.5,
    "qmax": 0.5,
}
# ============================================================================

CLI_DEFAULTS = {
    "npoints": 9,
    "qmin": -0.5,
    "qmax": 0.5,
}


def pick(args, key, default=None):
    """Use explicit CLI value first, then internal config or CLI default."""
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
        description="Generate nonrad CCD structures from two relaxed endpoint CONTCAR files."
    )
    p.add_argument("--mode", choices=["cli", "internal"], default=CONFIG_MODE,
                   help="Input mode. Default is controlled by CONFIG_MODE at the top of the script.")
    p.add_argument("--ground-contcar", default=None)
    p.add_argument("--excited-contcar", default=None)
    p.add_argument("--outdir", default=None)
    p.add_argument("--npoints", type=int, default=None)
    p.add_argument("--qmin", type=float, default=None)
    p.add_argument("--qmax", type=float, default=None)
    return p, p.parse_args()


def main():
    parser, args = parse_args()
    cfg = {
        "ground_contcar": pick(args, "ground_contcar"),
        "excited_contcar": pick(args, "excited_contcar"),
        "outdir": pick(args, "outdir"),
        "npoints": pick(args, "npoints", CLI_DEFAULTS["npoints"]),
        "qmin": pick(args, "qmin", CLI_DEFAULTS["qmin"]),
        "qmax": pick(args, "qmax", CLI_DEFAULTS["qmax"]),
    }
    require(parser, cfg, ["ground_contcar", "excited_contcar", "outdir"])

    ground_contcar = Path(cfg["ground_contcar"]).expanduser()
    excited_contcar = Path(cfg["excited_contcar"]).expanduser()
    outdir = Path(cfg["outdir"]).expanduser()

    for path in (ground_contcar, excited_contcar):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty structure file: {path}")

    outdir.mkdir(parents=True, exist_ok=True)
    ground_struct = Structure.from_file(ground_contcar)
    excited_struct = Structure.from_file(excited_contcar)

    dQ = get_dQ(ground_struct, excited_struct)
    displacements = np.linspace(cfg["qmin"], cfg["qmax"], cfg["npoints"])

    ground_structs, excited_structs = get_cc_structures(
        ground_struct,
        excited_struct,
        displacements,
        remove_zero=False,
    )

    meta = outdir / "ccd_points.tsv"
    with meta.open("w", encoding="utf-8") as f:
        f.write("index\tdisplacement\tQ_approx\n")
        for i, disp in enumerate(displacements):
            f.write(f"{i}\t{disp:.12f}\t{disp*dQ:.12f}\n")

    for branch, structs in [("ground", ground_structs), ("excited", excited_structs)]:
        for i, struct in enumerate(structs):
            d = outdir / branch / str(i)
            d.mkdir(parents=True, exist_ok=True)
            struct.to(filename=str(d / "POSCAR"), fmt="poscar")

    summary = outdir / "ccd_summary.txt"
    summary.write_text(
        f"mode = {args.mode}\n"
        f"dQ = {dQ:.12f} amu^(1/2) Angstrom\n"
        f"npoints = {cfg['npoints']}\n"
        f"qmin = {cfg['qmin']}\n"
        f"qmax = {cfg['qmax']}\n"
        f"ground_contcar = {ground_contcar}\n"
        f"excited_contcar = {excited_contcar}\n",
        encoding="utf-8",
    )

    print(f"mode = {args.mode}")
    print(f"dQ = {dQ:.12f} amu^(1/2) Angstrom")
    print(f"Wrote CCD structures to: {outdir}")
    print(f"Wrote metadata: {meta}")


if __name__ == "__main__":
    main()
