"""Module 2 demo: the metrology engine.

    python3 examples/run_module2_demo.py

Figures (outputs/figures/):
    20_spectrum_fit.png          raw spectrum, model fit, residuals
    21_thickness_landscape.png   why fitting has many wrong answers
    22_uncertainty_vs_thickness  when refractive index can and cannot be measured
    23_gauge_study.png           repeatability and precision-to-tolerance
    24_wafer_measurement.png     a Module 1 wafer measured site by site
    25_wrong_stack_model.png     what a wrong stack model does to the answer
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from filmtwin.chamber import run_scenario  # noqa: E402
from filmtwin.metrology import Reflectometer, gauge_study, measure_wafer  # noqa: E402
from filmtwin.optics import cauchy_n  # noqa: E402
from filmtwin.scenarios import SCENARIOS  # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "outputs", "figures")


def fig_spectrum(meter):
    R = meter.measure(500.0, 1.478)
    f = meter.fit(R)
    model = meter.model(f.thickness_nm, f.n633 - 0.0035 / 0.633 ** 2, f.scale).ravel()
    fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax[0].plot(meter.wl, R, ".", ms=3, label="measured (simulated)")
    ax[0].plot(meter.wl, model, "-", lw=1.5, label="model fit")
    ax[0].set_ylabel("Reflectance")
    ax[0].set_title(f"PECVD SiO2 on Si. True 500.0 nm, n=1.4780 | Fit {f.thickness_nm:.2f} "
                    f"+/- {f.sigma_thickness_nm:.2f} nm, n={f.n633:.4f}, GOF {f.gof:.4f}", fontsize=9)
    ax[0].legend(fontsize=8)
    ax[1].plot(meter.wl, R - model, ".", ms=2, color="0.3")
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_ylabel("Residual")
    ax[1].set_xlabel("Wavelength (nm)")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "20_spectrum_fit.png"), dpi=140)
    plt.close(fig)
    return f


def fig_landscape(meter):
    R = meter.measure(500.0, 1.478)
    grid = np.arange(10, 2000, 1.0)
    sse = meter.thickness_landscape(R, grid, 1.47)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.semilogy(grid, sse, lw=1)
    ax.axvline(500, color="crimson", ls="--", lw=1, label="true thickness (500 nm)")
    ax.set_xlabel("Trial thickness (nm)")
    ax.set_ylabel("Sum of squared error")
    ax.set_title("Fit error vs thickness: many local minima. A local optimizer started in the "
                 "wrong valley returns the wrong answer, so search globally first.", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "21_thickness_landscape.png"), dpi=140)
    plt.close(fig)


def fig_uncertainty(meter):
    ds = np.array([15, 20, 30, 40, 60, 80, 100, 150, 200, 300, 500, 800, 1200, 1800])
    sn, bias = [], []
    for d in ds:
        R = meter.measure(float(d), 1.478)
        sn.append(meter.fit(R, fit_index=True).sigma_n633)
        # Price of fixing n: assume 1.470 when the film is really 1.478
        bias.append((meter.fit(R, fit_index=False, n633_guess=1.470).thickness_nm - d) / d * 100)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].loglog(ds, sn, "-o", ms=4)
    ax[0].axhline(0.005, color="crimson", ls="--", lw=1, label="auto mode fixes n above this")
    ax[0].set_xlabel("Film thickness (nm)")
    ax[0].set_ylabel("1-sigma uncertainty in n (633 nm)")
    ax[0].set_title("Thin films: too few fringes, so n and d cannot be separated", fontsize=9)
    ax[0].legend(fontsize=8)
    ax[1].semilogx(ds, bias, "-o", ms=4, color="tab:orange")
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_xlabel("Film thickness (nm)")
    ax[1].set_ylabel("Thickness error (%)")
    ax[1].set_title("Price of fixing n wrong (assumed 1.470, true 1.478)", fontsize=9)
    for a in ax:
        a.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "22_uncertainty_vs_thickness.png"), dpi=140)
    plt.close(fig)


def fig_gauge(meter):
    g = gauge_study(meter, 500.0, 1.478, repeats=40, tolerance_nm=30.0)  # spec 500 +/- 3 %
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].hist(g["thickness_values"], bins=12, color="tab:blue", alpha=0.8)
    ax[0].axvline(500, color="crimson", ls="--")
    ax[0].set_xlabel("Measured thickness (nm)")
    ax[0].set_title(f"Repeatability {g['repeatability_nm']:.2f} nm, bias {g['bias_nm']:+.2f} nm, "
                    f"P/T {g['p_to_t']:.3f}", fontsize=9)
    ax[1].hist(g["n633_values"], bins=12, color="tab:green", alpha=0.8)
    ax[1].axvline(1.478, color="crimson", ls="--")
    ax[1].set_xlabel("Measured n (633 nm)")
    ax[1].set_title(f"Repeatability {g['repeatability_n633']:.4f}, bias {g['bias_n633']:+.4f}", fontsize=9)
    fig.suptitle("Gauge study: same site measured 40 times", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "23_gauge_study.png"), dpi=140)
    plt.close(fig)
    return g


def fig_wafer(meter):
    _, wafers, sites = run_scenario(SCENARIOS["showerhead_clog"], n_wafers=560, seed=11)
    w = sites[sites.wafer_index == 559]
    fits = measure_wafer(meter, w.thickness_nm.values, w.refractive_index.values)
    meas = np.array([f.thickness_nm for f in fits])
    err = meas - w.thickness_nm.values
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    ax[0].plot(w.thickness_nm, meas, "o", ms=4)
    lo, hi = w.thickness_nm.min() - 2, w.thickness_nm.max() + 2
    ax[0].plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax[0].set_xlabel("True thickness at site (nm)")
    ax[0].set_ylabel("Measured by Module 2 (nm)")
    ax[0].set_title("Module 1 wafer (late showerhead clog) measured site by site", fontsize=9)
    sc = ax[1].scatter(w.x_mm, w.y_mm, c=err, cmap="RdBu_r", vmin=-1, vmax=1, s=60, edgecolors="k",
                       linewidths=0.3)
    ax[1].set_aspect("equal")
    ax[1].set_title(f"Measurement error per site: mean {err.mean():+.2f} nm, sd {err.std(ddof=1):.2f} nm",
                    fontsize=9)
    plt.colorbar(sc, ax=ax[1], label="nm")
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "24_wafer_measurement.png"), dpi=140)
    plt.close(fig)
    return err


def fig_wrong_stack(meter):
    """True stack has 80 nm silicon nitride under the oxide; a single-layer model ignores it."""
    N_sin = cauchy_n(meter.wl, 1.95, B=0.018).astype(complex)
    true_under = [(N_sin, 80.0)]
    R = meter.measure(500.0, 1.478, under_layers=true_under)
    wrong = meter.fit(R)  # model assumes oxide directly on Si
    right_meter = Reflectometer(under_layers=true_under, seed=99)
    right = right_meter.fit(R)
    m_wrong = meter.model(wrong.thickness_nm, wrong.n633 - 0.0035 / 0.633 ** 2, wrong.scale).ravel()
    m_right = right_meter.model(right.thickness_nm, right.n633 - 0.0035 / 0.633 ** 2, right.scale).ravel()
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(meter.wl, R, ".", ms=3, color="0.3", label="measured: 500 nm SiO2 / 80 nm SiN / Si")
    ax.plot(meter.wl, m_wrong, "-", label=f"single-layer model: {wrong.thickness_nm:.1f} nm, "
                                          f"n {wrong.n633:.3f}, GOF {wrong.gof:.3f}")
    ax.plot(meter.wl, m_right, "-", label=f"correct 2-layer model: {right.thickness_nm:.1f} nm, "
                                          f"n {right.n633:.4f}, GOF {right.gof:.4f}")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Reflectance")
    ax.set_title("A wrong stack model returns a wrong thickness. The low GOF is the warning: always check it, not just the number.",
                 fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "25_wrong_stack_model.png"), dpi=140)
    plt.close(fig)
    return wrong, right


def main():
    os.makedirs(FIG, exist_ok=True)
    meter = Reflectometer(seed=2026)
    f = fig_spectrum(meter)
    print(f"[1/6] Spectrum fit: {f.thickness_nm:.2f} nm, n {f.n633:.4f}, GOF {f.gof:.4f}")
    fig_landscape(meter)
    print("[2/6] Thickness landscape saved")
    fig_uncertainty(meter)
    print("[3/6] Uncertainty vs thickness saved")
    g = fig_gauge(meter)
    print(f"[4/6] Gauge study: repeatability {g['repeatability_nm']:.3f} nm, "
          f"n {g['repeatability_n633']:.5f}, P/T {g['p_to_t']:.3f}")
    err = fig_wafer(meter)
    print(f"[5/6] Wafer measured: error mean {err.mean():+.3f} nm, sd {err.std(ddof=1):.3f} nm")
    wrong, right = fig_wrong_stack(meter)
    print(f"[6/6] Wrong stack: single-layer {wrong.thickness_nm:.1f} nm (GOF {wrong.gof:.3f}, ok={wrong.ok}) "
          f"vs 2-layer {right.thickness_nm:.1f} nm (GOF {right.gof:.4f})")
    print(f"\nFigures in {os.path.abspath(FIG)}")


if __name__ == "__main__":
    main()
