"""Module 3 demo: detect and diagnose chamber faults, scored against ground truth.

    python3 examples/run_module3_demo.py

Trains on two healthy runs, then tests on every scenario (different random
seeds), on three chambers the signatures were NOT built from, and on a fault
that is not in the signature library.

Outputs:
    outputs/figures/30-36_*.png      figures
    outputs/module3_results.csv      scorecard
    outputs/reports/*.md             auto-generated incident reports
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

from filmtwin.chamber import ChamberSimulator, run_scenario  # noqa: E402
from filmtwin.diagnosis import (Diagnoser, ReferenceModel, build_signatures, score,  # noqa: E402
                                shewhart_detection)
from filmtwin.features import FEATURES, wafer_features  # noqa: E402
from filmtwin.physics import Hardware, Recipe  # noqa: E402
from filmtwin.plots import _wafer_outline  # noqa: E402
from filmtwin.scenarios import SCENARIOS  # noqa: E402
from filmtwin.wafer import Wafer, dense_grid  # noqa: E402
from filmtwin.zernike import zernike  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
FIG = os.path.join(ROOT, "outputs", "figures")
REP = os.path.join(ROOT, "outputs", "reports")
FAULTS = ["edge_tc_drift", "showerhead_clog", "mfc_drift", "rf_match_degradation",
          "edge_ring_erosion", "leveling_tilt"]
UNITS = {"edge_tc_drift": "C colder", "showerhead_clog": "fraction of centre gas lost",
         "mfc_drift": "fractional flow error", "rf_match_degradation": "fraction of RF power lost",
         "edge_ring_erosion": "extra edge plasma", "leveling_tilt": "tilt per radius"}


def features_for(scenario, seed, hardware=None):
    _, w, s = run_scenario(scenario, seed=seed, hardware=hardware)
    return wafer_features(w, s)


def fault_start(sc):
    return min(f.start_wafer for f in sc.faults) if sc.faults else None


# ---------------------------------------------------------------------------
def fig_zernike_basis():
    wafer = Wafer()
    X, Y, mask = dense_grid(wafer, 101)
    rho = np.hypot(X, Y) / wafer.usable_radius_mm
    th = np.arctan2(Y, X)
    shapes = [(1, 1, "tilt_x"), (1, -1, "tilt_y"), (2, 0, "defocus (bowl)"),
              (2, 2, "astig_0 (saddle)"), (2, -2, "astig_45"), (4, 0, "spherical (ring)")]
    fig, axes = plt.subplots(1, 6, figsize=(16, 3))
    for ax, (n, m, name) in zip(axes, shapes):
        Z = np.where(mask, zernike(n, m, rho, th), np.nan)
        ax.contourf(X, Y, Z, levels=21, cmap="RdBu_r")
        _wafer_outline(ax, wafer)
        ax.set_title(name, fontsize=9)
    fig.suptitle("Zernike building blocks: every wafer map is a weighted sum of shapes like these", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "30_zernike_basis.png"), dpi=130)
    plt.close(fig)


def fig_signatures(diag):
    sd = np.sqrt(np.diag(diag.reference.cov))
    S = diag.signatures.copy()
    unit = {k: v[1] for k, v in __import__("filmtwin.diagnosis", fromlist=["SIGNATURES"]).SIGNATURES.items()}
    Z = S.apply(lambda col: col * unit[col.name]).div(sd, axis=0)  # z-units per typical unit fault
    fig, ax = plt.subplots(figsize=(9, 6.5))
    lim = 10
    im = ax.imshow(np.clip(Z.values, -lim, lim), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(Z.shape[1]))
    ax.set_xticklabels(Z.columns, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(Z.shape[0]))
    ax.set_yticklabels(Z.index, fontsize=8)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            v = Z.values[i, j]
            if abs(v) >= 0.5:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(v) > 6 else "black")
    plt.colorbar(im, ax=ax, label="change in healthy-noise sigmas (clipped at +/-10)")
    ax.set_title("Fault fingerprint table: what one small unit of each fault does to each feature",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "31_signature_table.png"), dpi=140)
    plt.close(fig)


def fig_detection(all_results):
    names = list(all_results)
    fig, axes = plt.subplots(len(names), 1, figsize=(11, 2.0 * len(names)), sharex=True)
    for ax, name in zip(axes, names):
        r, sc, s = all_results[name]
        ax.semilogy(r.wafer_index, r.t2, ".", ms=2, color="tab:blue")
        ax.axhline(r.t2_limit.iloc[0], color="k", ls="--", lw=0.8)
        st = fault_start(sc)
        if st is not None:
            ax.axvline(st, color="crimson", ls=":", lw=1.2)
        if s.get("detect_wafer") is not None:
            ax.axvline(s["detect_wafer"], color="tab:green", lw=1.5)
        ax.set_ylabel(name, fontsize=7, rotation=0, ha="right", va="center")
        ax.grid(alpha=0.3, which="both")
    axes[-1].set_xlabel("Wafer (red dotted = fault starts, green = alarm raised, dashed = T2 limit)")
    fig.suptitle("Multivariate detection (Hotelling T-squared) on every scenario", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "32_detection.png"), dpi=130)
    plt.close(fig)


def fig_severity(all_results, feats_by_name):
    panels = [n for n in all_results if SCENARIOS[n].faults]
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, 2.1 * len(panels)), sharex=True)
    for ax, name in zip(axes, panels):
        r, sc, _ = all_results[name]
        f = feats_by_name[name]
        for kind in sorted({fe.kind for fe in sc.faults}):
            truth = f[f"gt_{kind}"].abs()
            est = r[f"est_{kind}"].abs()
            scale = truth.max() if truth.max() > 0 else 1
            ax.plot(r.wafer_index, est / scale, ".", ms=2, label=f"{kind} estimated")
            ax.plot(f.wafer_index, truth / scale, "k-", lw=1)
        ax.set_ylim(-0.1, 1.3)
        ax.set_ylabel(name, fontsize=7, rotation=0, ha="right", va="center")
        ax.legend(fontsize=6, loc="upper left")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("Wafer (black line = true severity, dots = estimated; both scaled to the fault's max)")
    fig.suptitle("Severity tracking: the decomposition estimates how bad each fault is, wafer by wafer",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "33_severity_tracking.png"), dpi=130)
    plt.close(fig)


def fig_comparison(table):
    t = table[table.scenario != "baseline"].copy()
    x = np.arange(len(t))
    fig, ax = plt.subplots(figsize=(11, 4))
    big = 650
    for k, (col, lab) in enumerate([("t2_delay", "Multivariate T2 (this work)"),
                                    ("nu_delay", "Uniformity (NU) chart"),
                                    ("mean_delay", "Mean thickness chart")]):
        vals = t[col].fillna(big).to_numpy(dtype=float)
        bars = ax.bar(x + (k - 1) * 0.27, vals, 0.27, label=lab)
        for b, v in zip(bars, t[col]):
            if pd.isna(v):
                ax.text(b.get_x() + b.get_width() / 2, big / 2, "MISSED", ha="center", va="center",
                        fontsize=7, rotation=90, color="white", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(t.scenario, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Wafers from fault start to alarm")
    ax.set_title("Detection delay (lower is better). Only the T2 method also says WHAT is wrong.", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "34_method_comparison.png"), dpi=140)
    plt.close(fig)


def fig_case(diag, feats, results, sites, wafer_index, title, path):
    r = results[results.wafer_index == wafer_index].iloc[0]
    s = sites[sites.wafer_index == wafer_index].sort_values("site")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    sc = ax[0].scatter(s.x_mm, s.y_mm, c=s.thickness_nm, cmap="viridis", s=80, edgecolors="k", linewidths=0.3)
    _wafer_outline(ax[0], Wafer())
    plt.colorbar(sc, ax=ax[0], label="nm")
    ax[0].set_title(f"Wafer {wafer_index}: measured 49-site map", fontsize=9)
    z = evidence(diag, feats, wafer_index)
    top = z.abs().sort_values(ascending=False).head(8).index
    ax[1].barh(range(len(top)), z[top].values, color=["tab:red" if v > 0 else "tab:blue" for v in z[top]])
    ax[1].set_yticks(range(len(top)))
    ax[1].set_yticklabels(top, fontsize=8)
    ax[1].invert_yaxis()
    ax[1].set_xlabel("Deviation from healthy (sigmas)")
    ax[1].set_title("Evidence: which features moved", fontsize=9)
    contrib = {g: r[f"contrib_{g}"] for g in FAULTS}
    ax[2].bar(range(len(contrib)), list(contrib.values()), color="tab:purple")
    ax[2].axhline(diag.min_significance, color="k", ls="--", lw=0.8)
    ax[2].set_xticks(range(len(contrib)))
    ax[2].set_xticklabels(list(contrib), rotation=30, ha="right", fontsize=8)
    ax[2].set_ylabel("Contribution (sigmas)")
    ax[2].set_title(f"Diagnosis: {r.diagnosis} (explains {r.explained:.0%})", fontsize=9)
    for a in ax[1:]:
        a.grid(alpha=0.3)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def evidence(diag, feats, wafer_index):
    row = feats[feats.wafer_index == wafer_index]
    R = diag.reference.residuals(row)[0]
    return pd.Series(R / np.sqrt(np.diag(diag.reference.cov)), index=diag.reference.features)


def incident_report(name, diag, feats, results, sc):
    """Short escalation-style write-up, generated from the diagnosis."""
    alarms = results[results.alarm]
    lines = [f"# Incident report: chamber CH-A, scenario `{name}`", ""]
    if alarms.empty:
        lines += ["No alarm raised. Chamber within its healthy reference."]
        return "\n".join(lines)
    first = int(alarms.wafer_index.iloc[0])
    last = results.iloc[-1]
    z = evidence(diag, feats, int(last.wafer_index))
    top = z.abs().sort_values(ascending=False).head(4).index
    lines += [
        f"**Detection.** Multivariate alarm first raised at wafer {first} "
        f"(lot {feats.loc[feats.wafer_index == first, 'lot_id'].iloc[0]}).",
        "",
        f"**Current diagnosis (wafer {int(last.wafer_index)}).** `{last.diagnosis}`, "
        f"explaining {last.explained:.0%} of the deviation.",
        "",
        "**Estimated severity:**",
    ]
    for g in FAULTS:
        v = last[f"est_{g}"]
        if last[f"contrib_{g}"] >= diag.min_significance:
            lines.append(f"- {g}: {v:+.3g} ({UNITS[g]})")
    lines += ["", "**Evidence (largest deviations from healthy, in sigmas):**"]
    lines += [f"- {k}: {z[k]:+.1f}" for k in top]
    # Playbook comes from the DIAGNOSED fault(s), not from the scenario's known answer.
    diagnosed = [] if last.diagnosis in ("none", "unknown") else last.diagnosis.split("+")
    if not diagnosed:
        lines += ["", "**Playbook.** No known signature explains this pattern. Escalate: check recent "
                  "maintenance, recipe or hardware changes, and anything not yet in the fault library."]
    for g in diagnosed:
        play = SCENARIOS[g]
        lines += ["", f"**Playbook for `{g}`.**",
                  f"- Expected map pattern: {play.map_signature}",
                  f"- Sensor check: {play.sensor_clue}",
                  f"- Action: {play.typical_fix}"]
    lines += ["", "_Ground truth for scoring (not available in a real fab): " +
              f"{feats.gt_label.iloc[-1]}._"]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
def main():
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(REP, exist_ok=True)

    print("[1/6] Training on two healthy runs")
    train = pd.concat([features_for(SCENARIOS["baseline"], s) for s in (100, 101)], ignore_index=True)
    diag = Diagnoser.train(train)
    fig_zernike_basis()
    fig_signatures(diag)

    print("[2/6] Testing on every scenario (new random seeds)")
    rows, all_results, feats_by_name, sites_by_name = [], {}, {}, {}
    for i, (name, sc) in enumerate(SCENARIOS.items()):
        _, w, s = run_scenario(sc, seed=300 + i)
        f = wafer_features(w, s)
        r = diag.run(f)
        st = fault_start(sc)
        sres = score(r, f.gt_label, st)
        ref_start = st if st is not None else 0
        nu_det, nu_fa = shewhart_detection(train.m_nu_1sigma_pct, f.m_nu_1sigma_pct, f.wafer_index, ref_start)
        mn_det, mn_fa = shewhart_detection(train.mean_nm, f.mean_nm, f.wafer_index, ref_start)
        rows.append({
            "scenario": name,
            "fault_start": st,
            "t2_detect": sres.get("detect_wafer"),
            "t2_delay": sres.get("delay_wafers"),
            "nu_delay": (nu_det - st) if (st is not None and nu_det is not None) else None,
            "mean_delay": (mn_det - st) if (st is not None and mn_det is not None) else None,
            "false_alarm_rate": sres["false_alarm_rate"],
            "diagnosis_accuracy": sres["diagnosis_accuracy"],
            "baseline_nu_chart_alarm": nu_det if st is None else None,
            "baseline_mean_chart_alarm": mn_det if st is None else None,
        })
        all_results[name] = (r, sc, sres)
        feats_by_name[name] = f
        sites_by_name[name] = s
        with open(os.path.join(REP, f"{name}_incident.md"), "w") as fh:
            fh.write(incident_report(name, diag, f, r, sc))
    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(ROOT, "outputs", "module3_results.csv"), index=False)

    print("[3/6] Figures")
    fig_detection(all_results)
    fig_severity(all_results, feats_by_name)
    fig_comparison(table)
    fig_case(diag, feats_by_name["mixed_clog_and_match"], all_results["mixed_clog_and_match"][0],
             sites_by_name["mixed_clog_and_match"], 480,
             "Case study: two faults at once (showerhead clog + RF match), wafer 480",
             os.path.join(FIG, "35_case_mixed.png"))

    print("[4/6] Robustness: three chambers the signatures were not built from")
    robust = []
    for k in range(3):
        hw = Hardware.sample(np.random.default_rng(50 + k), name=f"CH-{k + 1}")
        tr = pd.concat([features_for(SCENARIOS["baseline"], s, hw) for s in (100, 101)], ignore_index=True)
        d2 = Diagnoser(ReferenceModel().fit(tr), build_signatures(hw=Hardware()))  # nominal signatures
        accs = []
        for i, (name, sc) in enumerate(SCENARIOS.items()):
            if not sc.faults:
                continue
            f = features_for(sc, 400 + 10 * k + i, hw)
            accs.append(score(d2.run(f), f.gt_label, fault_start(sc))["diagnosis_accuracy"])
        robust.append((hw.name, float(np.mean(accs))))
        print(f"      {hw.name}: mean diagnosis accuracy {np.mean(accs):.1%}")

    print("[5/6] Unseen fault: heater height drift (gap 10 -> 11.5 mm), not in the library")
    changes = {i: dataclasses.replace(Recipe(), spacing_mm=10 + 1.5 * min(1, (i - 150) / 300))
               for i in range(150, 600)}
    w, s = ChamberSimulator(seed=77).run(600, recipe_changes=changes)
    w["gt_label"] = np.where(w.wafer_index >= 150, "spacing_drift", "none")
    f = wafer_features(w, s)
    r = diag.run(f)
    flagged = r[(r.wafer_index >= 150) & r.out_of_limit]
    share_unknown = (flagged.diagnosis == "unknown").mean()
    print(f"      flagged {len(flagged)} wafers; labelled 'unknown' {share_unknown:.0%} of the time")
    fig_case(diag, f, r, s, 560, "Unseen fault (heater height drift): the method says 'unknown' "
             "instead of blaming the wrong part", os.path.join(FIG, "36_case_unseen_fault.png"))

    print("[6/6] Scorecard")
    show = table[["scenario", "fault_start", "t2_delay", "nu_delay", "mean_delay", "false_alarm_rate",
                  "diagnosis_accuracy"]]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print(f"\nFigures in {os.path.abspath(FIG)}\nReports in {os.path.abspath(REP)}")


if __name__ == "__main__":
    main()
