#!/usr/bin/env python
"""Calculate C(T) with either CLI arguments or an internal config block."""
from pathlib import Path
import argparse
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from nonrad import get_C

# ============================================================================
# USER-EDITABLE INTERNAL MODE
# Fill the physical inputs below and set CONFIG_MODE = "internal" to run with:
#   python scripts/04_calculate_capture.py
# ============================================================================
CONFIG_MODE = "cli"  # "cli" or "internal"
INTERNAL_CONFIG = {
    "case_root": "/absolute/path/to/nonrad_case",
    "dQ": None,
    "dE": None,
    "wi": None,
    "wf": None,
    "Wif": None,
    "volume": None,
    "g": 1.0,
    "Tmin": 0.0,
    "Tmax": 800.0,
    "nT": 100,
    "m_eff": 0.2,
}
# ============================================================================

CLI_DEFAULTS = {
    "case_root": "nonrad_case",
    "g": 1.0,
    "Tmin": 0.0,
    "Tmax": 800.0,
    "nT": 100,
    "m_eff": 0.2,
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
        description="Calculate nonradiative capture coefficient C(T) using nonrad.get_C."
    )
    p.add_argument("--mode", choices=["cli", "internal"], default=CONFIG_MODE)
    p.add_argument("--case-root", default=None)
    p.add_argument("--dQ", type=float, default=None)
    p.add_argument("--dE", type=float, default=None)
    p.add_argument("--wi", type=float, default=None)
    p.add_argument("--wf", type=float, default=None)
    p.add_argument("--Wif", type=float, default=None)
    p.add_argument("--volume", type=float, default=None)
    p.add_argument("--g", type=float, default=None)
    p.add_argument("--Tmin", type=float, default=None)
    p.add_argument("--Tmax", type=float, default=None)
    p.add_argument("--nT", type=int, default=None)
    p.add_argument("--m-eff", dest="m_eff", type=float, default=None,
                   help="Effective mass in m0 for capture cross-section estimate.")
    return p, p.parse_args()


def main():
    parser, args = parse_args()
    cfg = {
        "case_root": pick(args, "case_root", CLI_DEFAULTS["case_root"]),
        "dQ": pick(args, "dQ"),
        "dE": pick(args, "dE"),
        "wi": pick(args, "wi"),
        "wf": pick(args, "wf"),
        "Wif": pick(args, "Wif"),
        "volume": pick(args, "volume"),
        "g": pick(args, "g", CLI_DEFAULTS["g"]),
        "Tmin": pick(args, "Tmin", CLI_DEFAULTS["Tmin"]),
        "Tmax": pick(args, "Tmax", CLI_DEFAULTS["Tmax"]),
        "nT": pick(args, "nT", CLI_DEFAULTS["nT"]),
        "m_eff": pick(args, "m_eff", CLI_DEFAULTS["m_eff"]),
    }
    require(parser, cfg, ["case_root", "dQ", "dE", "wi", "wf", "Wif", "volume"])

    if cfg["nT"] < 2:
        parser.error("nT must be at least 2.")
    if cfg["Tmax"] < cfg["Tmin"]:
        parser.error("Tmax must be greater than or equal to Tmin.")
    if cfg["m_eff"] <= 0:
        parser.error("m_eff must be positive.")

    case = Path(cfg["case_root"]).expanduser()
    results = case / "05_results"
    results.mkdir(parents=True, exist_ok=True)

    T = np.linspace(cfg["Tmin"], cfg["Tmax"], cfg["nT"])
    c_of_T = get_C(
        dQ=cfg["dQ"],
        dE=cfg["dE"],
        wi=cfg["wi"],
        wf=cfg["wf"],
        Wif=cfg["Wif"],
        volume=cfg["volume"],
        g=cfg["g"],
        T=T,
    )

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(T, c_of_T, lw=2.2, label="capture coefficient")
    ax.set_yscale("log")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"Capture coefficient $C$ (cm$^3$ s$^{-1}$)")
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(results / "capture_coefficient_vs_T.png", dpi=600)
    plt.close(fig)

    kb = 1.380649e-23
    m0 = 9.10938356e-31
    idx_300 = int(np.abs(T - 300.0).argmin())
    T_actual = T[idx_300]
    v_th_cm_s = np.sqrt((8.0 * kb * T_actual) / (np.pi * cfg["m_eff"] * m0)) * 100.0
    C_300 = c_of_T[idx_300]
    sigma_300 = C_300 / v_th_cm_s

    out = results / "capture_summary.txt"
    with out.open("w", encoding="utf-8") as f:
        f.write(f"mode = {args.mode}\n")
        f.write(f"case_root = {case}\n")
        for key in ("dQ", "dE", "wi", "wf", "Wif", "volume", "g", "m_eff"):
            f.write(f"{key} = {cfg[key]}\n")
        f.write(f"T_300_used = {T_actual:.6f}\n")
        f.write(f"C_300K = {C_300:.8e} cm^3 s^-1\n")
        f.write(f"v_th_300K = {v_th_cm_s:.8e} cm s^-1\n")
        f.write(f"sigma_300K = {sigma_300:.8e} cm^2\n")

    np.savetxt(
        results / "capture_coefficient_vs_T.csv",
        np.column_stack([T, c_of_T]),
        delimiter=",",
        header="T_K,C_cm3_s",
        comments="",
    )

    print(out.read_text(encoding="utf-8"))
    print(f"Wrote figure: {results / 'capture_coefficient_vs_T.png'}")


if __name__ == "__main__":
    main()
