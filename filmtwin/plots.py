"""Plotting helpers: wafer maps, run trends, fault gallery, regime plot."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # works from Terminal with no display
import matplotlib.pyplot as plt
import numpy as np

from .metrics import uniformity
from .physics import Condition, Hardware, Recipe, deposit
from .wafer import Wafer, dense_grid, sites_49pt


def _wafer_outline(ax, wafer: Wafer):
    th = np.linspace(0, 2 * np.pi, 361)
    R = wafer.radius_mm
    ax.plot(R * np.cos(th), R * np.sin(th), color="0.3", lw=1.0)
    # notch at the bottom
    ax.plot([0], [-R], marker="^", color="0.3", ms=5)
    ax.set_aspect("equal")
    ax.set_xlim(-R * 1.05, R * 1.05)
    ax.set_ylim(-R * 1.05, R * 1.05)
    ax.set_xticks([])
    ax.set_yticks([])


def wafer_map(ax, values_grid, X, Y, mask, wafer=Wafer(), sites=None, site_values=None,
              title="", cmap="viridis", vmin=None, vmax=None, label="nm"):
    """Contour map from a dense grid, with measurement sites on top."""
    Z = np.where(mask, values_grid, np.nan)
    levels = np.linspace(vmin, vmax, 25) if (vmin is not None and vmax is not None) else 24
    cs = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, extend="neither")
    if sites is not None:
        ax.scatter(sites[:, 0], sites[:, 1], c="white", s=9, edgecolors="black", linewidths=0.4)
    _wafer_outline(ax, wafer)
    ax.set_title(title, fontsize=9)
    cb = plt.colorbar(cs, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return cs


def ideal_dense_map(recipe: Recipe, hw: Hardware, cond: Condition, wafer=Wafer(), n=121,
                    quantity="thickness_nm"):
    X, Y, mask = dense_grid(wafer, n)
    xy = np.column_stack([X.ravel(), Y.ravel()])
    r = deposit(xy, recipe, hw, cond, noise=False, usable_radius_mm=wafer.usable_radius_mm)
    return X, Y, mask, getattr(r, quantity).reshape(X.shape)


def plot_baseline(recipe: Recipe, hw: Hardware, path: str):
    """Thickness, refractive index, WERR and temperature maps for one recipe."""
    wafer = Wafer()
    sites = sites_49pt(wafer)
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.8))
    specs = [("thickness_nm", "Thickness", "nm", "viridis"),
             ("temperature_C", "Wafer temperature", "deg C", "inferno"),
             ("werr", "Wet etch rate ratio", "x thermal ox", "magma_r"),
             ("damkohler", "Damkohler number k_s/h_g", "-", "cividis")]
    for ax, (q, name, lab, cmap) in zip(axes, specs):
        X, Y, mask, Z = ideal_dense_map(recipe, hw, Condition(), wafer, quantity=q)
        title = name
        if q == "thickness_nm":
            r = deposit(sites, recipe, hw, Condition(), noise=False)
            u = uniformity(r.thickness_nm)
            title = f"{name}: mean {u['mean']:.1f} nm, NU {u['nu_1sigma_pct']:.2f}% (1s)"
        wafer_map(ax, Z, X, Y, mask, wafer, sites if q == "thickness_nm" else None,
                  title=title, cmap=cmap, label=lab)
    fig.suptitle("FilmTwin baseline: healthy chamber, default PECVD SiO2 recipe", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_fault_gallery(recipe: Recipe, hw: Hardware, conditions: dict, path: str):
    """Baseline map plus the *difference* map each fault creates."""
    wafer = Wafer()
    X, Y, mask, base = ideal_dense_map(recipe, hw, Condition(), wafer)
    n = len(conditions) + 1
    cols = 4
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.7 * rows))
    axes = np.atleast_1d(axes).ravel()
    wafer_map(axes[0], base, X, Y, mask, wafer, title="Healthy baseline (thickness)", label="nm")
    for ax, (name, cond) in zip(axes[1:], conditions.items()):
        _, _, _, Z = ideal_dense_map(recipe, hw, cond, wafer)
        d = (Z - base) / base * 100.0
        # Each panel gets its own symmetric scale so small signatures stay visible.
        lim = max(float(np.nanmax(np.abs(np.where(mask, d, np.nan)))), 0.05)
        wafer_map(ax, d, X, Y, mask, wafer, title=f"{name}: change vs baseline",
                  cmap="RdBu_r", vmin=-lim, vmax=lim, label="% thickness")
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Fault signature gallery (each fault at full severity)", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_trend(wafers, scenario, path: str):
    """Run chart: thickness, NU, RI/WERR and the most useful sensors."""
    i = wafers["wafer_index"].values
    fig, axes = plt.subplots(5, 1, figsize=(11, 11), sharex=True)
    pm = wafers.loc[wafers["pm_event"], "wafer_index"].values
    starts = [f.start_wafer for f in scenario.faults]

    def deco(ax):
        for p in pm:
            ax.axvline(p, color="0.5", ls="--", lw=0.8)
        for s in starts:
            ax.axvline(s, color="crimson", ls=":", lw=1.2)
        ax.grid(alpha=0.3)

    axes[0].plot(i, wafers["m_thickness_mean_nm"], ".", ms=3)
    axes[0].set_ylabel("Mean thk (nm)")
    axes[1].plot(i, wafers["m_nu_1sigma_pct"], ".", ms=3, color="tab:orange")
    axes[1].set_ylabel("NU 1s (%)")
    ax2 = axes[2]
    ax2.plot(i, wafers["m_ri_mean"], ".", ms=3, color="tab:green")
    ax2.set_ylabel("RI", color="tab:green")
    ax2b = ax2.twinx()
    ax2b.plot(i, wafers["m_werr_mean"], ".", ms=3, color="tab:purple")
    ax2b.set_ylabel("WERR", color="tab:purple")
    axes[3].plot(i, wafers["s_rf_reflected_W"], ".", ms=3, color="tab:red", label="RF reflected (W)")
    axes[3].plot(i, wafers["s_vpp_V"] - 250, ".", ms=3, color="tab:brown", label="Vpp - 250 (V)")
    axes[3].legend(fontsize=7, loc="upper left")
    axes[3].set_ylabel("RF sensors")
    axes[4].plot(i, wafers["s_heater_edge_pct"], ".", ms=3, color="tab:blue", label="Edge heater %")
    axes[4].plot(i, wafers["s_heater_center_pct"], ".", ms=3, color="tab:cyan", label="Centre heater %")
    axes[4].legend(fontsize=7, loc="lower left")
    axes[4].set_ylabel("Heater power")
    axes[4].set_xlabel("Wafer count (grey dashed = chamber clean, red dotted = fault start)")
    for ax in axes:
        deco(ax)
    fig.suptitle(f"Scenario '{scenario.name}': {scenario.description}", fontsize=10, wrap=True)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def plot_regime(pecvd: dict, thermal: dict, path: str):
    """Arrhenius plot of rate vs 1000/T for PECVD and a thermal CVD contrast."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for d, lab in ((pecvd, "PECVD (Ea 8 kJ/mol)"), (thermal, "Thermal CVD (Ea 150 kJ/mol)")):
        invT = 1000.0 / (d["T_C"] + 273.15)
        axes[0].semilogy(invT, d["rate_nm_min"], "-o", ms=3, label=lab)
        axes[1].plot(d["T_C"], d["Da"], "-o", ms=3, label=lab)
    axes[0].set_xlabel("1000 / T (1/K)")
    axes[0].set_ylabel("Centre deposition rate (nm/min)")
    axes[0].set_title("Grove curve: steep = reaction-limited, flat = transport-limited", fontsize=9)
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3, which="both")
    axes[1].axhline(1.0, color="k", lw=0.8, ls="--")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Temperature (C)")
    axes[1].set_ylabel("Damkohler number k_s / h_g")
    axes[1].set_title("Da < 1 reaction-limited | Da > 1 transport-limited", fontsize=9)
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
