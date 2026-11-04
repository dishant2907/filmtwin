"""Module 4 demo: DOE, SPC with an executable OCAP, and run-to-run control.

    python3 examples/run_module4_demo.py

Outputs:
    outputs/figures/40-46_*.png
    outputs/doe/*.csv               design + results tables (open them in JMP or Excel)
    outputs/doe/doe_report.md       the experiment written up
    outputs/reports/ocap_*.md       OCAP execution log for a real excursion
"""
from __future__ import annotations

import dataclasses
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from filmtwin import doe  # noqa: E402
from filmtwin.chamber import ChamberSimulator, run_scenario  # noqa: E402
from filmtwin.diagnosis import Diagnoser  # noqa: E402
from filmtwin.features import wafer_features  # noqa: E402
from filmtwin.metrics import uniformity  # noqa: E402
from filmtwin.physics import Recipe, deposit  # noqa: E402
from filmtwin.r2r import EWMAController, estimate_post_pm_boost, run_closed_loop, thickness_error_stats  # noqa: E402
from filmtwin.scenarios import SCENARIOS  # noqa: E402
from filmtwin.spc import (OCAP_STEPS, ContextLimits, Limits, capability, count_events, cusum_chart,  # noqa: E402
                          ewma_chart, first_signal, run_ocap, western_electric)
from filmtwin.stats import t_ppf  # noqa: E402
from filmtwin.wafer import sites_49pt  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
FIG = os.path.join(ROOT, "outputs", "figures")
DOE_DIR = os.path.join(ROOT, "outputs", "doe")
REP = os.path.join(ROOT, "outputs", "reports")
TARGET, HALF_SPEC = 500.0, 15.0  # spec 500 +/- 15 nm (+/- 3 %)
RSM_FACTORS = {"edge_temp_offset_C": (0.0, 30.0), "pressure_torr": (2.0, 3.0), "precursor_sccm": (120.0, 180.0)}


# ===========================================================================
# Part A: DOE
# ===========================================================================
def part_a_doe():
    os.makedirs(DOE_DIR, exist_ok=True)
    print("[A1] Screening: 2^(6-2) fractional factorial + 4 centre points (20 wafers)")
    design = doe.fractional_factorial_2_6_2(n_center=4, seed=1)
    scr = doe.run_on_twin(design, seed=2)
    scr.to_csv(os.path.join(DOE_DIR, "screening_results.csv"), index=False)
    tables = {}
    for resp in ["rate_nm_min", "nu_pct", "ri", "werr_mean"]:
        fit, tab = doe.analyse_screening(scr, resp)
        tables[resp] = (fit, tab)
    fig_pareto(tables)

    # Vital few: the strongest drivers of uniformity and of rate
    nu_tab = tables["nu_pct"][1]
    rate_tab = tables["rate_nm_min"][1]
    top = lambda tab: tab[tab.term.isin(doe.FACTORS)].assign(a=lambda d: d.t.abs()).nlargest(2, "a").term.tolist()
    chosen = sorted(set(top(nu_tab) + top(rate_tab)), key=list(RSM_FACTORS).index)
    print(f"      vital few for uniformity and rate: {chosen}")

    print("[A2] Response surface: face-centred CCD on 3 factors (20 wafers)")
    factors = list(RSM_FACTORS)
    ccd = doe.ccd_face_centered(factors, n_center=6, seed=3)
    rsm = doe.run_on_twin(ccd, seed=4, factors=RSM_FACTORS)
    rsm.to_csv(os.path.join(DOE_DIR, "ccd_results.csv"), index=False)
    models = {r: doe.fit_response_surface(rsm, r, factors) for r in doe.RESPONSES}

    print("[A3] Desirability optimisation")
    best, pred, D = doe.optimise({k: models[k] for k in doe.DEFAULT_GOALS}, factors)
    best_real = doe.to_real(best.to_frame().T, RSM_FACTORS).iloc[0]
    recipe_opt = dataclasses.replace(Recipe(), **{k: float(v) for k, v in best_real.items()})
    fig_contours(models, factors, best)

    print("[A4] Confirmation: 5 wafers at the chosen recipe")
    X, _ = doe.model_matrix(best.to_frame().T, factors, "quadratic")
    rng = np.random.default_rng(9)
    conf = []
    for _ in range(5):
        r = deposit(sites_49pt(), recipe_opt, rng=rng)
        u = uniformity(r.thickness_nm)
        t = r.thickness_nm
        conf.append({"rate_nm_min": u["mean"] / recipe_opt.dep_time_s * 60, "nu_pct": u["nu_1sigma_pct"],
                     "profile_pct": (t[25:].mean() - t[:9].mean()) / t.mean() * 100,
                     "ri": r.refractive_index.mean(), "werr_mean": r.werr.mean()})
    conf = pd.DataFrame(conf)
    pis = {resp: models[resp].predict(X, 0.95, n_new=len(conf))
           for resp in ["profile_pct", "rate_nm_min", "ri", "werr_mean", "nu_pct"]}
    fig_confirmation(conf, pis)

    # baseline for comparison
    base = deposit(sites_49pt(), Recipe(), noise=False)
    ub = uniformity(base.thickness_nm)
    summary = {
        "baseline_nu": ub["nu_1sigma_pct"], "baseline_rate": ub["mean"] / Recipe().dep_time_s * 60,
        "opt_nu": conf.nu_pct.mean(), "opt_rate": conf.rate_nm_min.mean(),
        "opt_recipe": best_real.to_dict(), "desirability": D,
    }
    write_doe_report(tables, models, best_real, pred, pis, conf, summary, chosen)
    print(f"      uniformity {summary['baseline_nu']:.2f}% -> {summary['opt_nu']:.2f}%, "
          f"rate {summary['baseline_rate']:.0f} -> {summary['opt_rate']:.0f} nm/min")
    return summary


def fig_pareto(tables):
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
    for ax, (resp, (fit, tab)) in zip(axes, tables.items()):
        t = tab[tab.term.isin(doe.FACTORS)].copy()
        t["abs_t"] = t.t.abs()
        t = t.sort_values("abs_t")
        crit = t_ppf(0.975, fit.df_resid)
        ax.barh(t.term, t.abs_t, color=np.where(t.abs_t > crit, "tab:blue", "0.75"))
        ax.axvline(crit, color="crimson", ls="--", lw=1, label=f"p = 0.05 (|t| = {crit:.2f})")
        ax.set_xscale("log")
        ax.set_title(f"{resp}  (R2 {fit.r2:.3f})", fontsize=9)
        ax.tick_params(labelsize=7)
        ax.set_xlabel("|t| of effect (log scale)", fontsize=8)
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("Screening DOE, 20 wafers: which knobs matter for each response (blue = significant)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "40_doe_screening_pareto.png"), dpi=140)
    plt.close(fig)


def fig_contours(models, factors, best):
    a, b, c = factors  # edge offset (x), pressure (y), precursor held at optimum
    g = np.linspace(-1, 1, 81)
    A, B = np.meshgrid(g, g)
    coded = pd.DataFrame({a: A.ravel(), b: B.ravel(), c: np.full(A.size, best[c])})
    X, _ = doe.model_matrix(coded, factors, "quadratic")
    P = {r: (X @ m.coef).reshape(A.shape) for r, m in models.items()}
    real_a = doe.to_real(pd.DataFrame({a: g}), RSM_FACTORS)[a].to_numpy()
    real_b = doe.to_real(pd.DataFrame({b: g}), RSM_FACTORS)[b].to_numpy()
    RA, RB = np.meshgrid(real_a, real_b)
    ok = (np.abs(P["profile_pct"]) < 1.0) & (P["rate_nm_min"] > 300) & (P["werr_mean"] < 1.6)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.contourf(RA, RB, ok.astype(float), levels=[0.5, 1.5], colors=["#c8e6c9"])
    cs1 = ax.contour(RA, RB, P["profile_pct"], levels=[-3, -2, -1, 0, 1, 2, 3], colors="tab:blue", linewidths=1)
    ax.clabel(cs1, fmt="profile %.0f%%", fontsize=7)
    cs2 = ax.contour(RA, RB, P["rate_nm_min"], levels=8, colors="tab:orange", linewidths=1, linestyles="--")
    ax.clabel(cs2, fmt="%.0f nm/min", fontsize=7)
    br = doe.to_real(best.to_frame().T, RSM_FACTORS).iloc[0]
    ax.plot(br[a], br[b], "r*", ms=16, label="desirability optimum")
    ax.plot(0.0, 2.5, "ks", ms=9, label="original recipe", clip_on=False, zorder=5)
    ax.set_xlabel("Edge heater trim (C)")
    ax.set_ylabel("Pressure (Torr)")
    ax.set_title(f"Process window at precursor {br[c]:.0f} sccm. Green = flat profile (|edge-centre| < 1 %), "
                 "rate > 300 nm/min, WERR < 1.6", fontsize=8)
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "41_doe_process_window.png"), dpi=140)
    plt.close(fig)


def fig_confirmation(conf, pis):
    fig, axes = plt.subplots(1, len(pis), figsize=(15, 3.6))
    for ax, (resp, (yhat, half)) in zip(axes, pis.items()):
        obs = conf[resp]
        ax.errorbar([0], yhat, yerr=half, fmt="o", color="tab:blue", capsize=8,
                    label=f"prediction +/- 95 % PI for the mean of {len(obs)} runs")
        ax.plot(np.full(len(obs), 0.3), obs, "x", color="0.6", ms=6, label="individual wafers")
        ax.plot([0.3], [obs.mean()], "rD", ms=7, label="confirmation mean")
        inside = abs(obs.mean() - yhat[0]) <= half[0]
        ax.set_title(f"{resp}: {'confirmed' if inside else 'NOT confirmed'}", fontsize=9,
                     color="black" if inside else "crimson")
        ax.set_xlim(-0.5, 0.8)
        ax.set_xticks([])
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=6)
    fig.suptitle("Confirmation: mean of 5 runs at the optimum vs the model's prediction interval", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "42_doe_confirmation.png"), dpi=140)
    plt.close(fig)


def write_doe_report(tables, models, best_real, pred, pis, conf, summary, chosen):
    L = ["# DOE report: PECVD oxide recipe optimisation on FilmTwin", "",
         "## 1. Screening (2^(6-2) resolution IV, 16 runs + 4 centre points, randomised)", ""]
    for resp, (fit, tab) in tables.items():
        L.append(f"**{resp}** (R² {fit.r2:.3f}, model p {fit.f_p:.2g})")
        L.append("")
        L.append("| term | effect | t | p |")
        L.append("|---|---|---|---|")
        for _, r in tab[tab.term != "intercept"].iterrows():
            eff = "" if pd.isna(r.effect) else f"{r.effect:+.4g}"
            L.append(f"| {r.term} | {eff} | {r.t:+.2f} | {r.p_value:.3g} |")
        L.append("")
    L += [f"Vital few chosen for the response surface: {', '.join(chosen)}.", "",
          "## 2. Response surface (face-centred CCD, 20 runs)", "",
          "| response | R² | adj R² | RMSE |", "|---|---|---|---|"]
    for r, m in models.items():
        L.append(f"| {r} | {m.r2:.4f} | {m.r2_adj:.4f} | {m.rmse:.4g} |")
    L += ["", "Uniformity is modelled through the signed edge-minus-centre profile, not 1-sigma NU. "
          "NU is a standard deviation and folds into a V at its minimum, which a quadratic cannot fit "
          "(compare the R² values).", "",
          "## 3. Optimum (Derringer desirability)", "",
          "| factor | value |", "|---|---|"]
    for k, v in best_real.items():
        L.append(f"| {k} | {v:.3g} |")
    L += [f"| overall desirability | {summary['desirability']:.3f} |", "",
          "Pressure sits at the edge of the explored range. The next experiment should extend it, "
          "checking constraints this model does not include (particles from gas-phase nucleation at "
          "high pressure, film stress).", "",
          "## 4. Confirmation (5 wafers)", "",
          "Verdict compares the mean of the confirmation runs with the 95 % prediction interval for a "
          "mean of 5 new runs. The last column puts the miss in engineering terms: the difference as a "
          "share of the window that matters for that response.", "",
          "| response | predicted | 95 % PI | observed mean | verdict | miss as % of window |",
          "|---|---|---|---|---|---|"]
    windows = {"profile_pct": 3.0, "rate_nm_min": 10.0, "ri": 0.008, "werr_mean": 0.1, "nu_pct": 0.5}
    for resp, (yhat, half) in pis.items():
        obs = conf[resp]
        ok = abs(obs.mean() - yhat[0]) <= half[0]
        miss = abs(obs.mean() - yhat[0]) / windows[resp] * 100
        L.append(f"| {resp} | {yhat[0]:.6g} | ±{half[0]:.3g} | {obs.mean():.6g} | "
                 f"{'confirmed' if ok else 'outside PI'} | {miss:.1f} % (±{windows[resp]:g}) |")
    L += ["", "**Reading this table like an engineer.** A response can fail the statistical test and still "
          "not matter. Refractive index is measured so precisely that its prediction interval is tiny, so a "
          "model error far too small to affect the film becomes statistically detectable. Also, with five "
          "responses each tested at 95 %, there is roughly a 23 % chance that at least one fails by chance "
          "alone (1 - 0.95^5). The decision rule: a statistical miss that is also a small fraction of the "
          "engineering window is noted, not acted on; a miss that is a large fraction of the window means "
          "the model is wrong in that region and needs more experiments."]
    L += ["", f"**Result:** 1-sigma non-uniformity {summary['baseline_nu']:.2f} % → {summary['opt_nu']:.2f} %, "
          f"deposition rate {summary['baseline_rate']:.0f} → {summary['opt_rate']:.0f} nm/min, "
          "refractive index held on target.", "",
          "Raw tables (`screening_results.csv`, `ccd_results.csv`) can be re-analysed in JMP: "
          "Analyze → Fit Model, with the coded factor columns as effects."]
    with open(os.path.join(DOE_DIR, "doe_report.md"), "w") as fh:
        fh.write("\n".join(L))


# ===========================================================================
# Part B: SPC + OCAP (open loop)
# ===========================================================================
def part_b_spc(diag):
    print("[B1] Phase I limits from two healthy runs (context-corrected for wafers since clean)")
    p1 = pd.concat([run_scenario(SCENARIOS["baseline"], seed=s)[1] for s in (100, 101)], ignore_index=True)
    ctx = ContextLimits.from_phase1(p1.m_thickness_mean_nm, p1.wafers_since_clean)
    lim = Limits.from_phase1(ctx.residual(p1.m_thickness_mean_nm, p1.wafers_since_clean))

    _, h, _ = run_scenario(SCENARIOS["baseline"], seed=205)
    hr = ctx.residual(h.m_thickness_mean_nm, h.wafers_since_clean)
    print(f"      false-alarm episodes on a fresh healthy run (600 wafers): "
          f"WE1 only {count_events(western_electric(hr, lim).we1)}, all four WE rules "
          f"{count_events(western_electric(hr, lim)['any'])}, EWMA {count_events(ewma_chart(hr, lim).signal)}, "
          f"CUSUM {count_events(cusum_chart(hr, lim).signal)}")
    print("[B2] Phase II: monitoring the MFC drift excursion")
    sc = SCENARIOS["mfc_drift"]
    _, w, s = run_scenario(sc, seed=203)
    resid = ctx.residual(w.m_thickness_mean_nm, w.wafers_since_clean)
    we = western_electric(resid, lim)
    ew = ewma_chart(resid, lim, lam=0.2)
    cu = cusum_chart(resid, lim)
    st = sc.faults[0].start_wafer
    det = {"I chart (WE rules)": first_signal(we["any"], w.wafer_index, st),
           "EWMA": first_signal(ew.signal, w.wafer_index, st),
           "CUSUM": first_signal(cu.signal, w.wafer_index, st)}
    print("      first signal after fault start: " + ", ".join(f"{k} {v}" for k, v in det.items()))
    fig_spc(w, resid, lim, we, ew, cu, st)

    trig = det["I chart (WE rules)"]
    feats = wafer_features(w, s)
    sensor_limits = {c: (p1[c].mean() - 4 * p1[c].std(), p1[c].mean() + 4 * p1[c].std())
                     for c in ["s_rf_reflected_W", "s_heater_edge_pct", "s_heater_center_pct", "s_throttle_valve_pct"]}
    events = []
    for idx in [trig, 260, 450]:
        rule = "we1" if we.we1.iloc[idx] else "trend rule"
        expected_here = float(w.m_thickness_mean_nm.iloc[idx] - resid[idx])  # context-expected thickness
        events.append(run_ocap(idx, rule, w, Limits(expected_here, lim.sigma, lim.mr_bar),
                               (TARGET - HALF_SPEC, TARGET + HALF_SPEC), sensor_limits, diag, feats))
    write_ocap(events)
    for e in events:
        print(f"      OCAP at wafer {e.wafer_index}: {e.disposition}")
    return p1


def fig_spc(w, resid, lim, we, ew, cu, start):
    i = w.wafer_index
    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    ax[0].plot(i, resid, ".", ms=3, color="0.4")
    for k, c in (("we1", "red"), ("we2", "orange"), ("we3", "gold"), ("we4", "purple")):
        m = we[k].to_numpy()
        ax[0].plot(i[m], resid[m], "o", ms=4, mfc="none", color=c, label=k.upper())
    for v, ls in ((lim.center, "-"), (lim.ucl, "--"), (lim.lcl, "--")):
        ax[0].axhline(v, color="k", ls=ls, lw=0.8)
    ax[0].set_ylabel("Thickness residual (nm)")
    ax[0].legend(fontsize=7, ncol=4, loc="lower left")
    ax[0].set_title("I chart on context-corrected thickness, Western Electric rules marked", fontsize=9)
    ax[1].plot(i, ew.ewma, "-", lw=1)
    ax[1].plot(i, ew.ucl, "k--", lw=0.8)
    ax[1].plot(i, ew.lcl, "k--", lw=0.8)
    ax[1].set_ylabel("EWMA (nm)")
    ax[2].plot(i, cu.cusum_hi, label="CUSUM high")
    ax[2].plot(i, cu.cusum_lo, label="CUSUM low")
    ax[2].axhline(cu.h.iloc[0], color="k", ls="--", lw=0.8)
    ax[2].set_ylabel("CUSUM (sigmas)")
    ax[2].legend(fontsize=7)
    ax[2].set_xlabel("Wafer (red dotted = fault starts)")
    for a in ax:
        a.axvline(start, color="crimson", ls=":", lw=1.2)
        a.grid(alpha=0.3)
    fig.suptitle("SPC on the MFC drift excursion: three chart types, same data", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "43_spc_charts.png"), dpi=130)
    plt.close(fig)


def write_ocap(events):
    os.makedirs(REP, exist_ok=True)
    L = ["# OCAP execution log: MFC drift excursion", "", "## The plan", "", "```", OCAP_STEPS.strip(), "```", ""]
    for e in events:
        L += [f"## Trigger at wafer {e.wafer_index} ({e.trigger})", ""]
        L += [f"- {line}" for line in e.log]
        L += ["", f"**Disposition:** {e.disposition}", ""]
    with open(os.path.join(REP, "ocap_mfc_drift.md"), "w") as fh:
        fh.write("\n".join(L))


# ===========================================================================
# Part C: run-to-run control
# ===========================================================================
def part_c_r2r(diag, p1_open):
    boost = estimate_post_pm_boost(p1_open)
    print(f"[C1] Feedforward learned from history: first wafers after a clean need "
          f"+{boost[0] * 100:.1f} % and +{boost[1] * 100:.1f} % time")
    sc = SCENARIOS["mfc_drift"]
    runs = {}
    for label, ctrl in [("open loop", None), ("EWMA feedback", EWMAController(lam=0.6)),
                        ("EWMA feedback + feedforward", EWMAController(lam=0.6, post_pm_boost=boost))]:
        w, s = run_closed_loop(ChamberSimulator(faults=list(sc.faults), seed=7), 600, ctrl)
        runs[label] = (w, s)
        st = thickness_error_stats(w, TARGET, HALF_SPEC)
        cap = capability(w.m_thickness_mean_nm, TARGET - HALF_SPEC, TARGET + HALF_SPEC)
        print(f"      {label:30s} RMSE {st['rmse_nm']:5.2f} nm, worst {st['max_abs_error_nm']:5.1f} nm, "
              f"out of spec {st['pct_out_of_spec']:4.1f} %, Ppk {cap['Ppk']:.2f}")

    print("[C2] Monitoring the controller itself (engineering adjustment limit)")
    p1c = pd.concat([run_closed_loop(ChamberSimulator(seed=sd), 600,
                                     EWMAController(lam=0.6, post_pm_boost=boost))[0] for sd in (6, 8)],
                    ignore_index=True)
    ctx_t = ContextLimits.from_phase1(p1c.r2r_time_s, p1c.wafers_since_clean)
    w_ctrl, s_ctrl = runs["EWMA feedback + feedforward"]
    alarm, dev = ctx_t.adjustment_alarm(w_ctrl.r2r_time_s, w_ctrl.wafers_since_clean, 1.0)
    adj_det = first_signal(alarm, w_ctrl.wafer_index, sc.faults[0].start_wafer)
    f = wafer_features(w_ctrl, s_ctrl)
    r = diag.run(f)
    t2_det = first_signal(r.alarm, r.wafer_index, sc.faults[0].start_wafer)
    acc = (r.loc[r.out_of_limit & (f.gt_label != "none").values, "diagnosis"] == "mfc_drift").mean()
    print(f"      thickness held in spec, but the controller adjustment alarm fires at wafer {adj_det}; "
          f"Module 3 detects at {t2_det} and diagnoses mfc_drift on {acc:.0%} of flagged wafers")
    fig_r2r(runs, dev, alarm, sc.faults[0].start_wafer)

    print("[C3] Tuning lambda")
    lams = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    curves = {}
    for lag in (0, 2, 6):
        curves[lag] = [thickness_error_stats(run_closed_loop(ChamberSimulator(faults=list(sc.faults), seed=7),
                                                             600, EWMAController(lam=lam, post_pm_boost=boost),
                                                             metrology_lag=lag)[0], TARGET, HALF_SPEC)["rmse_nm"]
                       for lam in lams]
    fig_lambda(lams, curves)
    fig_capability(runs)
    return {"adj_detect": adj_det, "t2_detect": t2_det, "acc": acc}


def fig_r2r(runs, dev, alarm, start):
    fig, ax = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)
    for label, (w, _) in runs.items():
        ax[0].plot(w.wafer_index, w.m_thickness_mean_nm, ".", ms=2.5, label=label)
        if label != "open loop":
            ax[1].plot(w.wafer_index, w.r2r_time_s, "-", lw=1, label=label)
    ax[0].axhspan(TARGET - HALF_SPEC, TARGET + HALF_SPEC, color="green", alpha=0.07, label="spec")
    ax[0].axhline(TARGET, color="k", lw=0.6)
    ax[0].set_ylabel("Mean thickness (nm)")
    ax[0].legend(fontsize=7, loc="lower left")
    ax[0].set_title("MFC drift: open loop drifts out of spec; the controller holds thickness on target", fontsize=9)
    ax[1].set_ylabel("Deposition time (s)")
    ax[1].legend(fontsize=7)
    ax[1].set_title("...by stretching deposition time: the drift moved from thickness into the controller", fontsize=9)
    w = runs["EWMA feedback + feedforward"][0]
    ax[2].plot(w.wafer_index, dev, "-", lw=1, color="tab:purple")
    ax[2].axhspan(-1, 1, color="green", alpha=0.07)
    ax[2].plot(w.wafer_index[alarm], dev[alarm], "r.", ms=3, label="adjustment alarm")
    ax[2].set_ylabel("Time vs expected (%)")
    ax[2].set_xlabel("Wafer (red dotted = fault starts)")
    ax[2].legend(fontsize=7)
    ax[2].set_title("Engineering limit (+/-1 %) on the controller output catches what the thickness chart can't see",
                    fontsize=9)
    for a in ax:
        a.axvline(start, color="crimson", ls=":", lw=1.2)
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "44_r2r_control.png"), dpi=130)
    plt.close(fig)


def fig_lambda(lams, curves):
    fig, ax = plt.subplots(figsize=(7, 4))
    for lag, rm in curves.items():
        ax.plot(lams, rm, "-o", ms=4, label=f"metrology delay {lag} wafers")
    ax.set_xlabel("EWMA weight lambda (0 = ignore new data, 1 = trust only the latest)")
    ax.set_ylabel("Thickness RMSE vs target (nm)")
    ax.set_title("Controller tuning: too low lags the drift, too high chases noise", fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "45_r2r_lambda_tuning.png"), dpi=140)
    plt.close(fig)


def fig_capability(runs):
    fig, ax = plt.subplots(figsize=(8, 4))
    bins = np.linspace(460, 525, 66)
    for label in ("open loop", "EWMA feedback + feedforward"):
        w = runs[label][0]
        cap = capability(w.m_thickness_mean_nm, TARGET - HALF_SPEC, TARGET + HALF_SPEC)
        ax.hist(w.m_thickness_mean_nm, bins=bins, alpha=0.6, label=f"{label}: Ppk {cap['Ppk']:.2f}")
    for v in (TARGET - HALF_SPEC, TARGET + HALF_SPEC):
        ax.axvline(v, color="crimson", ls="--")
    ax.set_xlabel("Mean thickness (nm)")
    ax.set_ylabel("Wafers")
    ax.set_title("Process capability during the MFC drift (spec 500 +/- 15 nm)", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "46_capability.png"), dpi=140)
    plt.close(fig)


# ===========================================================================
def main():
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(REP, exist_ok=True)
    doe_summary = part_a_doe()
    print("[B0] Training the Module 3 diagnoser for the OCAP")
    train = pd.concat([wafer_features(*run_scenario(SCENARIOS["baseline"], seed=s)[1:]) for s in (100, 101)],
                      ignore_index=True)
    diag = Diagnoser.train(train)
    p1 = part_b_spc(diag)
    part_c_r2r(diag, p1)
    print(f"\nFigures in {os.path.abspath(FIG)}\nDOE tables and report in {os.path.abspath(DOE_DIR)}")


if __name__ == "__main__":
    main()
