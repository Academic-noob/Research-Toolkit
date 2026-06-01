#!/usr/bin/env python
"""Extract Wif with either CLI arguments or an internal config block."""
from pathlib import Path
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from nonrad.elphon import get_Wif_from_WSWQ

# ============================================================================
# USER-EDITABLE INTERNAL MODE
# Note: def_index and bulk_index use 1-based band numbering. spin uses
# 0=up and 1=down. kpoint=1 selects the first k-point. Keep these consistent
# with your tested command.
# ============================================================================
CONFIG_MODE = "cli"  # "cli" or "internal"
INTERNAL_CONFIG = {
    "case_root": "/absolute/path/to/nonrad_case",
    "branch": "ground",       # "ground" or "excited"
    "ref_index": 4,
    "def_index": None,         # example: 192
    "bulk_index": None,        # example: [189, 190, 191]
    "spin": 1,
    "kpoint": 1,
    "include_zero": False,
}
# ============================================================================

CLI_DEFAULTS = {
    "case_root": "nonrad_case",
    "branch": "ground",
    "ref_index": 4,
    "spin": 1,
    "kpoint": 1,
    "include_zero": False,
}


def pick(args, key, default=None):
    cli_value = getattr(args, key)
    if cli_value is not None:
        return cli_value
    if args.mode == "internal":
        return INTERNAL_CONFIG.get(key, default)
    return default


def require(parser, values, keys):
    missing = [key for key in keys if values.get(key) in (None, "", [])]
    if missing:
        parser.error(
            "Missing required setting(s): " + ", ".join(missing)
            + ". Supply CLI options or fill INTERNAL_CONFIG and use --mode internal."
        )


def read_ccd_points(path: Path):
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty CCD metadata file: {path}")
    data = {}
    with path.open("r", encoding="utf-8") as f:
        next(f)
        for line in f:
            if not line.strip():
                continue
            idx, _disp, q = line.split()[:3]
            data[int(idx)] = float(q)
    return data


def parse_args():
    p = argparse.ArgumentParser(description="Extract Wif from WSWQ files using nonrad.")
    p.add_argument("--mode", choices=["cli", "internal"], default=CONFIG_MODE)
    p.add_argument("--case-root", default=None)
    p.add_argument("--branch", choices=["ground", "excited"], default=None)
    p.add_argument("--ref-index", type=int, default=None)
    p.add_argument("--def-index", type=int, default=None)
    p.add_argument("--bulk-index", type=int, nargs="+", default=None)
    p.add_argument("--spin", type=int, default=None)
    p.add_argument("--kpoint", type=int, default=None)
    p.add_argument("--include-zero", action=argparse.BooleanOptionalAction, default=None,
                   help="Include Q=0 WSWQ point if present.")
    return p, p.parse_args()


def main():
    parser, args = parse_args()
    cfg = {
        "case_root": pick(args, "case_root", CLI_DEFAULTS["case_root"]),
        "branch": pick(args, "branch", CLI_DEFAULTS["branch"]),
        "ref_index": pick(args, "ref_index", CLI_DEFAULTS["ref_index"]),
        "def_index": pick(args, "def_index"),
        "bulk_index": pick(args, "bulk_index"),
        "spin": pick(args, "spin", CLI_DEFAULTS["spin"]),
        "kpoint": pick(args, "kpoint", CLI_DEFAULTS["kpoint"]),
        "include_zero": pick(args, "include_zero", CLI_DEFAULTS["include_zero"]),
    }
    require(parser, cfg, ["case_root", "def_index", "bulk_index"])

    case = Path(cfg["case_root"]).expanduser()
    results = case / "05_results"
    results.mkdir(parents=True, exist_ok=True)

    ccd_points = read_ccd_points(case / "02_ccd_structures" / "ccd_points.tsv")
    wswq_dir = case / "04_wswq" / cfg["branch"]
    if not wswq_dir.is_dir():
        raise FileNotFoundError(f"Missing WSWQ directory: {wswq_dir}")

    wswqs = []
    numbered_dirs = [x for x in wswq_dir.iterdir() if x.is_dir() and x.name.isdigit()]
    for d in sorted(numbered_dirs, key=lambda x: int(x.name)):
        idx = int(d.name)
        if idx not in ccd_points:
            raise KeyError(f"Index {idx} from {d} is absent from ccd_points.tsv")
        q = ccd_points[idx]
        if (not cfg["include_zero"]) and abs(q) < 1e-10:
            continue
        w = d / "WSWQ"
        if w.exists() and w.stat().st_size > 0:
            wswqs.append((q, str(w)))

    if not wswqs:
        raise RuntimeError(f"No WSWQ files found in {wswq_dir}")

    initial_vasprun = case / "03_static" / cfg["branch"] / str(cfg["ref_index"]) / "vasprun.xml"
    if not initial_vasprun.is_file() or initial_vasprun.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty initial_vasprun: {initial_vasprun}")

    fig = plt.figure(figsize=(9, 5))
    Wifs_list = get_Wif_from_WSWQ(
        wswqs=wswqs,
        initial_vasprun=str(initial_vasprun),
        def_index=cfg["def_index"],
        bulk_index=cfg["bulk_index"],
        spin=cfg["spin"],
        kpoint=cfg["kpoint"],
        fig=fig,
    )

    Wif_rms = float(np.sqrt(np.mean([x[1] ** 2 for x in Wifs_list])))
    tag = f"{cfg['branch']}_spin{cfg['spin']}_k{cfg['kpoint']}"
    fig.tight_layout()
    fig.savefig(results / f"Wif_{tag}.png", dpi=600)
    plt.close(fig)

    out = results / f"Wif_{tag}.txt"
    with out.open("w", encoding="utf-8") as f:
        f.write(f"mode = {args.mode}\n")
        f.write(f"case_root = {case}\n")
        f.write(f"branch = {cfg['branch']}\n")
        f.write(f"initial_vasprun = {initial_vasprun}\n")
        f.write(f"def_index = {cfg['def_index']}\n")
        f.write(f"bulk_index = {cfg['bulk_index']}\n")
        f.write(f"spin = {cfg['spin']}\n")
        f.write(f"kpoint = {cfg['kpoint']}\n")
        f.write(f"include_zero = {cfg['include_zero']}\n")
        f.write("\n===== Wif channels =====\n")
        for band, val in Wifs_list:
            f.write(f"{band}\t{val:.12f}\n")
        f.write(f"\nWif_rms = {Wif_rms:.12f}\n")

    print(out.read_text(encoding="utf-8"))
    print(f"Wrote figure: {results / f'Wif_{tag}.png'}")


if __name__ == "__main__":
    main()
