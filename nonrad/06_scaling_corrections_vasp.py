#!/usr/bin/env python
"""
Calculate optional Sommerfeld and charged-supercell scaling factors.
Supports both CLI mode and an editable internal config block.
"""
from pathlib import Path
import argparse
import sys
import numpy as np

# Optional: force use of a local nonrad source containing
# charged_supercell_scaling_VASP. Fill this only when needed.
NONRAD_SRC = None  # example: "/home/student/soft/nonrad-1.2.0"
if NONRAD_SRC:
    sys.path.insert(0, NONRAD_SRC)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from nonrad.scaling import sommerfeld_parameter, charged_supercell_scaling_VASP

# ============================================================================
# USER-EDITABLE INTERNAL MODE
# Sommerfeld and charged-supercell scaling remain independently optional.
# Fill the parameters for the part(s) you want to calculate.
# ============================================================================
CONFIG_MODE = "cli"  # "cli" or "internal"
INTERNAL_CONFIG = {
    "case_root": "/absolute/path/to/nonrad_case",
    "outdir": None,  # None -> <case_root>/05_results/scaling
    "Z": None,
    "m_eff": None,
    "eps_static": None,
    "T": 300.0,
    "Tmin": 25.0,
    "Tmax": 800.0,
    "nT": 1000,
    "wavecar": None,       # example: "/absolute/path/to/WAVECAR"
    "band_index": None,    # bulk-like state, example: 189
    "def_index": None,     # defect state, example: 192
    "spin": 0,
    "kpoint": 1,
    "cutoff": 0.02,
    "limit": 5.0,
    "full_range": False,
    "tag": "scaling",
}
# ============================================================================

CLI_DEFAULTS = {
    "case_root": "nonrad_case",
    "outdir": None,
    "T": 300.0,
    "Tmin": 25.0,
    "Tmax": 800.0,
    "nT": 1000,
    "spin": 0,
    "kpoint": 1,
    "cutoff": 0.02,
    "limit": 5.0,
    "full_range": False,
    "tag": "scaling",
}


def pick(args, key, default=None):
    cli_value = getattr(args, key)
    if cli_value is not None:
        return cli_value
    if args.mode == "internal":
        return INTERNAL_CONFIG.get(key, default)
    return default


def parse_args():
    p = argparse.ArgumentParser(
        description="Calculate Sommerfeld and charged-supercell scaling factors for nonrad."
    )
    p.add_argument("--mode", choices=["cli", "internal"], default=CONFIG_MODE)
    p.add_argument("--case-root", default=None)
    p.add_argument("--outdir", default=None)
    p.add_argument("--Z", type=float, default=None)
    p.add_argument("--m-eff", dest="m_eff", type=float, default=None)
    p.add_argument("--eps-static", dest="eps_static", type=float, default=None)
    p.add_argument("--T", type=float, default=None)
    p.add_argument("--Tmin", type=float, default=None)
    p.add_argument("--Tmax", type=float, default=None)
    p.add_argument("--nT", type=int, default=None)
    p.add_argument("--wavecar", default=None)
    p.add_argument("--band-index", dest="band_index", type=int, default=None)
    p.add_argument("--def-index", dest="def_index", type=int, default=None)
    p.add_argument("--spin", type=int, default=None)
    p.add_argument("--kpoint", type=int, default=None)
    p.add_argument("--cutoff", type=float, default=None)
    p.add_argument("--limit", type=float, default=None)
    p.add_argument("--full-range", dest="full_range", action=argparse.BooleanOptionalAction, default=None)
    p.add_argument("--tag", default=None)
    return p.parse_args()


def main():
    args = parse_args()
    cfg = {
        key: pick(args, key, CLI_DEFAULTS.get(key))
        for key in INTERNAL_CONFIG
    }

    case_root = Path(cfg["case_root"]).expanduser()
    outdir = Path(cfg["outdir"]).expanduser() if cfg["outdir"] else case_root / "05_results" / "scaling"
    outdir.mkdir(parents=True, exist_ok=True)

    report = [
        "===== nonrad scaling corrections report =====",
        f"mode = {args.mode}",
        f"case_root = {case_root}",
        f"outdir = {outdir}",
        "",
    ]

    do_sommerfeld = (
        cfg["Z"] is not None
        and cfg["m_eff"] is not None
        and cfg["eps_static"] is not None
    )
    if do_sommerfeld:
        if cfg["nT"] < 2:
            raise ValueError("nT must be at least 2.")
        if cfg["Tmax"] < cfg["Tmin"]:
            raise ValueError("Tmax must be greater than or equal to Tmin.")
        if cfg["m_eff"] <= 0 or cfg["eps_static"] <= 0:
            raise ValueError("m_eff and eps_static must be positive.")

        f_single = sommerfeld_parameter(cfg["T"], cfg["Z"], cfg["m_eff"], cfg["eps_static"])
        T_grid = np.linspace(cfg["Tmin"], cfg["Tmax"], cfg["nT"])
        f_grid = sommerfeld_parameter(T_grid, cfg["Z"], cfg["m_eff"], cfg["eps_static"])

        report.extend([
            "===== Sommerfeld factor =====",
            f"Z = {cfg['Z']}",
            f"m_eff = {cfg['m_eff']} m0",
            f"eps_static = {cfg['eps_static']}",
            f"Sommerfeld factor @ {cfg['T']:.2f} K = {float(f_single):.8f}",
        ])
        if cfg["Z"] < 0:
            report.append("interpretation = attractive center; capture coefficient is enhanced by f(T).")
        elif cfg["Z"] > 0:
            report.append("interpretation = repulsive center; capture coefficient is suppressed by f(T).")
        else:
            report.append("interpretation = neutral long-range interaction; f(T) should be close to 1.")

        csv_path = outdir / f"{cfg['tag']}_sommerfeld_vs_T.csv"
        np.savetxt(csv_path, np.column_stack([T_grid, f_grid]), delimiter=",",
                   header="T_K,sommerfeld_factor", comments="")

        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(T_grid, f_grid, lw=2.0)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Sommerfeld factor")
        ax.grid(True, ls="--", alpha=0.5)
        fig.tight_layout()
        fig_path = outdir / f"{cfg['tag']}_sommerfeld_vs_T.png"
        fig.savefig(fig_path, dpi=600)
        plt.close(fig)
        report.extend([f"wrote_csv = {csv_path}", f"wrote_fig = {fig_path}", ""])
    else:
        report.extend([
            "===== Sommerfeld factor =====",
            "Skipped. Required: Z, m_eff, eps_static",
            "",
        ])

    do_charged_scaling = (
        cfg["wavecar"] is not None
        and cfg["band_index"] is not None
        and cfg["def_index"] is not None
    )
    if do_charged_scaling:
        wavecar = Path(cfg["wavecar"]).expanduser()
        if not wavecar.is_file() or wavecar.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty WAVECAR: {wavecar}")

        report.extend([
            "===== Charged-supercell scaling =====",
            f"wavecar = {wavecar}",
            f"band_index = {cfg['band_index']}",
            f"def_index = {cfg['def_index']}",
            f"spin = {cfg['spin']}",
            f"kpoint = {cfg['kpoint']}",
        ])
        fig = plt.figure(figsize=(12, 5))
        factor = charged_supercell_scaling_VASP(
            str(wavecar),
            cfg["band_index"],
            def_index=cfg["def_index"],
            spin=cfg["spin"],
            kpoint=cfg["kpoint"],
            cutoff=cfg["cutoff"],
            limit=cfg["limit"],
            full_range=cfg["full_range"],
            fig=fig,
        )
        plt.tight_layout()
        fig_path = outdir / f"{cfg['tag']}_charged_supercell_scaling.png"
        fig.savefig(fig_path, dpi=600)
        plt.close(fig)

        scaling = 1.0 / factor
        report.extend([
            f"raw_factor_from_charged_supercell_scaling = {float(factor):.12f}",
            f"recommended_Wif_scaling_1_over_factor = {float(scaling):.12f}",
            f"wrote_fig = {fig_path}",
            "",
            "Interpretation:",
            "  If Wif was computed in a charged supercell and this correction is needed,",
            "  use Wif_corrected = Wif_raw * recommended_Wif_scaling_1_over_factor.",
            "  If Wif was computed in a neutral charge state, this correction is usually not needed.",
            "",
        ])
    else:
        report.extend([
            "===== Charged-supercell scaling =====",
            "Skipped. Required: wavecar, band_index, def_index",
            "",
        ])

    report_path = outdir / f"{cfg['tag']}_scaling_report.txt"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(report_path.read_text(encoding="utf-8"))
    print(f"Wrote report: {report_path}")


if __name__ == "__main__":
    main()
