"""Module 4b: Statistical Process Control, capability and an executable OCAP.

Charts
    I-MR    individuals and moving range (one value per wafer)
    EWMA    exponentially weighted moving average: catches small sustained shifts
    CUSUM   cumulative sum: also built for small shifts, and shows when one began
Rules
    Western Electric rules 1-4 on the individuals chart
Capability
    Cp, Cpk (within-process sigma from the moving range) and Ppk (overall sigma)
OCAP
    An Out-of-Control Action Plan written as code: when a chart fires, it
    walks the same decision steps an engineer would, including calling the
    Module 3 diagnoser, and writes down what it did and why.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

D2_N2 = 1.128  # control-chart constant d2 for moving ranges of 2
D4_N2 = 3.267


# ---------------------------------------------------------------------------
# Limits from a Phase I (known good) period
# ---------------------------------------------------------------------------
@dataclass
class Limits:
    center: float
    sigma: float  # within-process sigma estimated from the average moving range
    mr_bar: float

    @classmethod
    def from_phase1(cls, x):
        x = np.asarray(x, float)
        mr = np.abs(np.diff(x))
        mr_bar = float(mr.mean())
        return cls(float(x.mean()), mr_bar / D2_N2, mr_bar)

    @property
    def ucl(self):
        return self.center + 3 * self.sigma

    @property
    def lcl(self):
        return self.center - 3 * self.sigma


@dataclass
class ContextLimits:
    """Limits for a signal that moves normally with chamber context.

    Example: a run-to-run controller's deposition time ramps up between
    cleans (compensating wall buildup) and drops back after each clean. A
    plain chart would alarm on that every cycle. Here the expected value is
    modelled against wafers-since-clean from Phase I data, and the chart runs
    on the residual. Sigma is the overall residual sigma, because a stepwise
    controller output has near-zero moving ranges.
    """

    coef: np.ndarray
    sigma: float

    @classmethod
    def from_phase1(cls, x, wafers_since_clean):
        from .diagnosis import _context_matrix
        X = _context_matrix(wafers_since_clean)
        coef, *_ = np.linalg.lstsq(X, np.asarray(x, float), rcond=None)
        r = np.asarray(x, float) - X @ coef
        return cls(coef, float(r.std(ddof=1)))

    def residual(self, x, wafers_since_clean):
        from .diagnosis import _context_matrix
        return np.asarray(x, float) - _context_matrix(wafers_since_clean) @ self.coef

    def as_limits(self):
        return Limits(0.0, self.sigma, self.sigma * D2_N2)

    def adjustment_alarm(self, x, wafers_since_clean, limit_pct: float = 1.0, run_rule=(3, 5)):
        """Engineering limit on a controller's output.

        Statistical limits fail here: a controller's output is autocorrelated
        by design (it is itself a smoothed average), so a 3-sigma chart
        assuming independent points raises false alarms constantly. Fabs
        instead set an engineering limit on how far the controller may move
        from its expected value before someone must look. Returns (alarm,
        deviation_pct).
        """
        x = np.asarray(x, float)
        resid = self.residual(x, wafers_since_clean)
        expected = x - resid
        dev = resid / expected * 100.0
        raw = (np.abs(dev) > limit_pct).astype(int)
        k, n = run_rule
        alarm = pd.Series(raw).rolling(n, min_periods=1).sum().to_numpy() >= k
        return alarm, dev


# ---------------------------------------------------------------------------
# Western Electric rules
# ---------------------------------------------------------------------------
def western_electric(x, lim: Limits) -> pd.DataFrame:
    """Flag each point against WE rules 1-4. Returns a boolean frame."""
    z = (np.asarray(x, float) - lim.center) / lim.sigma
    n = len(z)
    r1 = np.abs(z) > 3

    def window_rule(k, m, thresh):
        hit = np.zeros(n, bool)
        for side in (1, -1):
            beyond = (side * z) > thresh
            cnt = pd.Series(beyond.astype(int)).rolling(m, min_periods=m).sum().to_numpy()
            hit |= (cnt >= k) & beyond
        return hit

    r2 = window_rule(2, 3, 2.0)  # 2 of 3 beyond 2 sigma, same side
    r3 = window_rule(4, 5, 1.0)  # 4 of 5 beyond 1 sigma, same side
    r4 = window_rule(8, 8, 0.0)  # 8 in a row on one side
    return pd.DataFrame({"we1": r1, "we2": r2, "we3": r3, "we4": r4, "any": r1 | r2 | r3 | r4})


# ---------------------------------------------------------------------------
# EWMA and CUSUM
# ---------------------------------------------------------------------------
def ewma_chart(x, lim: Limits, lam: float = 0.2, L: float = 3.0):
    x = np.asarray(x, float)
    z = np.empty_like(x)
    prev = lim.center
    for i, v in enumerate(x):
        prev = lam * v + (1 - lam) * prev
        z[i] = prev
    i = np.arange(1, len(x) + 1)
    width = L * lim.sigma * np.sqrt(lam / (2 - lam) * (1 - (1 - lam) ** (2 * i)))
    return pd.DataFrame({"ewma": z, "ucl": lim.center + width, "lcl": lim.center - width,
                         "signal": (z > lim.center + width) | (z < lim.center - width)})


def cusum_chart(x, lim: Limits, k: float = 0.5, h: float = 5.0):
    """Tabular CUSUM with reference value k and decision interval h (in sigmas)."""
    z = (np.asarray(x, float) - lim.center) / lim.sigma
    hi = np.zeros_like(z)
    lo = np.zeros_like(z)
    for i in range(len(z)):
        hi[i] = max(0.0, (hi[i - 1] if i else 0.0) + z[i] - k)
        lo[i] = max(0.0, (lo[i - 1] if i else 0.0) - z[i] - k)
    return pd.DataFrame({"cusum_hi": hi, "cusum_lo": lo, "h": h, "signal": (hi > h) | (lo > h)})


def count_events(signal) -> int:
    """Number of separate alarm episodes (a run of consecutive signals counts once)."""
    s = np.asarray(signal, bool).astype(int)
    return int(np.sum(np.diff(np.concatenate([[0], s])) == 1))


def first_signal(signal, wafer_index, start=0):
    idx = np.asarray(wafer_index)
    s = np.asarray(signal, bool) & (idx >= start)
    return int(idx[s][0]) if s.any() else None


# ---------------------------------------------------------------------------
# Capability
# ---------------------------------------------------------------------------
def capability(x, lsl: float, usl: float) -> dict:
    """Cp/Cpk use within-process sigma (moving range); Pp/Ppk use overall sigma."""
    x = np.asarray(x, float)
    mu = float(x.mean())
    s_within = float(np.mean(np.abs(np.diff(x))) / D2_N2)
    s_overall = float(x.std(ddof=1))
    return {
        "mean": mu,
        "sigma_within": s_within,
        "sigma_overall": s_overall,
        "Cp": (usl - lsl) / (6 * s_within),
        "Cpk": min(usl - mu, mu - lsl) / (3 * s_within),
        "Pp": (usl - lsl) / (6 * s_overall),
        "Ppk": min(usl - mu, mu - lsl) / (3 * s_overall),
        "pct_out_of_spec": float(np.mean((x < lsl) | (x > usl)) * 100),
    }


# ---------------------------------------------------------------------------
# OCAP
# ---------------------------------------------------------------------------
OCAP_STEPS = """
OCAP: PECVD oxide mean thickness, chamber CH-A

Trigger: any Western Electric rule on the mean-thickness I chart, or an EWMA signal.

1. VERIFY THE MEASUREMENT. Re-measure the wafer. If the re-measurement is
   back inside the limits and the fit GOF is good, log a metrology false
   alarm and release the lot.
2. CHECK CONTEXT. Is this the first or second wafer after a chamber clean?
   If so, a low reading is the known seasoning effect: log it, continue.
3. CHECK TOOL SENSORS (FDC). Reflected RF power, edge and centre heater power,
   throttle valve. Any sensor outside its own limits points at hardware.
4. RUN THE DIAGNOSER (Module 3) on the last wafers. Record the probable
   cause, its estimated severity, and the share of the deviation explained.
5. DISPOSITION.
   - Out of spec: HOLD the lot, take the chamber down, escalate to the
     process owner. Restart only after a qualification wafer passes.
   - In spec, nothing confirms the chart (sensors normal, no multivariate
     anomaly): continue under watch, sample every wafer, re-run if it persists.
     Charts with several run rules do raise occasional false alarms.
   - In spec, cause identified: continue running, open a maintenance ticket
     for the diagnosed part, tighten sampling to every wafer.
   - In spec but an anomaly nobody can explain: HOLD and escalate.
"""


@dataclass
class OCAPResult:
    wafer_index: int
    trigger: str
    log: list = field(default_factory=list)
    disposition: str = ""


def run_ocap(wafer_index: int, trigger: str, wafers: pd.DataFrame, lim: Limits, spec: tuple,
             sensor_limits: dict, diagnoser=None, features: pd.DataFrame = None,
             remeasure_sigma_nm: float = 0.4, rng=None) -> OCAPResult:
    """Execute the OCAP for one triggering wafer and record every decision."""
    rng = rng or np.random.default_rng(0)
    res = OCAPResult(wafer_index, trigger)
    row = wafers.loc[wafers.wafer_index == wafer_index].iloc[0]
    y = float(row["m_thickness_mean_nm"])

    # 1. Verify measurement
    y2 = y + rng.normal(0, remeasure_sigma_nm)
    inside = lim.lcl <= y2 <= lim.ucl
    res.log.append(f"Step 1, verify: re-measured {y2:.1f} nm (first reading {y:.1f} nm). "
                   + ("Back inside limits: metrology false alarm." if inside and trigger == "we1"
                      else "Confirms the signal."))
    if inside and trigger == "we1":
        res.disposition = "RELEASE: metrology false alarm"
        return res

    # 2. Context
    wsc = int(row["wafers_since_clean"])
    if wsc <= 1 and y < lim.center:
        res.log.append(f"Step 2, context: wafer {wsc} after a clean, low reading is the known "
                       "seasoning effect.")
        res.disposition = "CONTINUE: expected first-wafer effect"
        return res
    res.log.append(f"Step 2, context: {wsc} wafers since clean, not a seasoning effect.")

    # 3. Sensors
    hits = []
    for col, (lo, hi) in sensor_limits.items():
        v = float(row[col])
        if v < lo or v > hi:
            hits.append(f"{col} = {v:.2f} (limits {lo:.2f} to {hi:.2f})")
    res.log.append("Step 3, sensors: " + ("; ".join(hits) if hits else "all inside their limits."))

    # 4. Diagnose
    cause, explained = "not run", np.nan
    if diagnoser is not None and features is not None:
        window = features[(features.wafer_index <= wafer_index) &
                          (features.wafer_index > wafer_index - 5)]
        d = diagnoser.run(window).iloc[-1]
        cause, explained = d["diagnosis"], d["explained"]
        res.log.append(f"Step 4, diagnosis: {cause} (explains {explained:.0%} of the deviation).")

    # 5. Disposition
    lsl, usl = spec
    in_spec = lsl <= y <= usl
    if not in_spec:
        res.disposition = "HOLD LOT, chamber down (out of spec); escalate to process owner."
    elif cause == "none" and not hits:
        res.disposition = ("CONTINUE UNDER WATCH: single-chart signal with no confirming evidence "
                           "(likely a chart false alarm or a very early drift); sample every wafer "
                           "for the next lot and re-run this OCAP if the signal persists.")
    elif cause not in ("unknown", "none", "not run"):
        res.disposition = (f"CONTINUE with ticket: thickness {y:.1f} nm still in spec; "
                           f"maintenance ticket for {cause}; sample every wafer until fixed.")
    else:
        res.disposition = "HOLD LOT, chamber down (anomaly with no identified cause); escalate."
    res.log.append("Step 5, disposition: " + res.disposition)
    return res
