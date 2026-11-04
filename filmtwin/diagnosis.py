"""Module 3: detect chamber problems and diagnose their root cause.

Three steps, the same logic a good process engineer follows:

1. REFERENCE  Learn what "healthy" looks like from known-good runs, including
              normal context effects (rate drift between cleans, the low first
              wafer). A drop that happens every PM cycle is not a fault.
2. DETECT     Compare each wafer to the reference. Hotelling's T-squared
              combines all features (map shape, film properties, sensors) into
              one multivariate statistic that accounts for how they normally
              move together. A run rule (3 of the last 5 wafers out) cuts
              false alarms.
3. DIAGNOSE   Explain the deviation as a combination of known fault
              signatures. Each signature is computed from the physics in
              Module 1 (what one unit of that fault does to every feature),
              so the method is physics-informed rather than a black-box
              classifier. Because it is a least-squares decomposition, it can
              report two faults at once and estimate each one's severity.
              Deviations no signature explains are flagged "unknown".
"""
from __future__ import annotations

import dataclasses

import numpy as np
import pandas as pd

from .chamber import SensorModel
from .faults import Aging, FaultEvent, build_condition
from .features import FEATURES, features_from_arrays
from .physics import Condition, Hardware, Recipe, deposit, gas_state
from .wafer import Wafer, sites_49pt

# Signatures to fit. Signed ones can go either direction in real life.
SIGNATURES = {
    # name: (fault kind, unit severity, angle, signed)
    "edge_tc_drift": ("edge_tc_drift", 5.0, 0.0, True),
    "showerhead_clog": ("showerhead_clog", 0.04, 0.0, False),
    "mfc_drift": ("mfc_drift", 0.02, 0.0, True),
    "rf_match_degradation": ("rf_match_degradation", 0.03, 0.0, False),
    "edge_ring_erosion": ("edge_ring_erosion", 0.03, 0.0, False),
    "leveling_tilt_x": ("leveling_tilt", 0.005, 0.0, True),
    "leveling_tilt_y": ("leveling_tilt", 0.005, 90.0, True),
}
FAULT_GROUPS = {
    "edge_tc_drift": ["edge_tc_drift"],
    "showerhead_clog": ["showerhead_clog"],
    "mfc_drift": ["mfc_drift"],
    "rf_match_degradation": ["rf_match_degradation"],
    "edge_ring_erosion": ["edge_ring_erosion"],
    "leveling_tilt": ["leveling_tilt_x", "leveling_tilt_y"],
}


# ---------------------------------------------------------------------------
# 1. Reference model of healthy behaviour
# ---------------------------------------------------------------------------
def _context_matrix(wafers_since_clean):
    n = np.asarray(wafers_since_clean, dtype=float)
    x = n / 100.0
    return np.column_stack([np.ones_like(x), x, x ** 2, (n == 0).astype(float),
                            (n == 1).astype(float)])


@dataclasses.dataclass
class ReferenceModel:
    features: list = dataclasses.field(default_factory=lambda: list(FEATURES))
    coef: np.ndarray = None
    cov: np.ndarray = None
    cov_inv: np.ndarray = None
    chol: np.ndarray = None
    t2_limit: float = None

    def fit(self, healthy: pd.DataFrame, alpha: float = 0.001):
        X = _context_matrix(healthy["wafers_since_clean"])
        Y = healthy[self.features].to_numpy()
        self.coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
        R = Y - X @ self.coef
        self.cov = np.cov(R, rowvar=False)
        self.cov_inv = np.linalg.inv(self.cov)
        self.chol = np.linalg.cholesky(self.cov)
        self.t2_limit = chi2_quantile(1 - alpha, len(self.features))
        return self

    def expected(self, df):
        return _context_matrix(df["wafers_since_clean"]) @ self.coef

    def residuals(self, df):
        return df[self.features].to_numpy() - self.expected(df)

    def t2(self, df):
        R = self.residuals(df)
        return np.einsum("ij,jk,ik->i", R, self.cov_inv, R)


def chi2_quantile(p: float, k: int) -> float:
    """Wilson-Hilferty approximation of the chi-square quantile (no scipy needed)."""
    z = _norm_ppf(p)
    return k * (1 - 2 / (9 * k) + z * np.sqrt(2 / (9 * k))) ** 3


def _norm_ppf(p: float) -> float:
    # Acklam's rational approximation, accurate to ~1e-9
    a = [-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02,
         1.383577518672690e02, -3.066479806614716e01, 2.506628277459239e00]
    b = [-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02,
         6.680131188771972e01, -1.328068155288572e01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00,
         -2.549732539343734e00, 4.374664141464968e00, 2.938163982698783e00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00,
         3.754408661907416e00]
    if p < 0.02425:
        q = np.sqrt(-2 * np.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    if p > 1 - 0.02425:
        return -_norm_ppf(1 - p)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)


# ---------------------------------------------------------------------------
# 2. Physics-based fault signatures
# ---------------------------------------------------------------------------
class _NoNoise:
    """Stand-in random generator that returns the mean (no noise)."""

    def normal(self, loc=0.0, scale=1.0, size=None):
        return np.full(size, loc) if size is not None else loc


def _ideal_features(recipe, hw, cond, xy):
    res = deposit(xy, recipe, hw, cond, noise=False)
    gas = gas_state(recipe, hw, cond, float(np.mean(res.temperature_C) + 273.15))
    sens = SensorModel().read(recipe, cond, hw, gas, _NoNoise())
    return features_from_arrays(res.thickness_nm, res.refractive_index, res.werr, sens,
                                recipe.dep_time_s).iloc[0]


def build_signatures(recipe: Recipe = None, hw: Hardware = None, features=FEATURES):
    """Feature change per unit severity for each signature (central differences)."""
    recipe = recipe or Recipe()
    hw = hw or Hardware()
    xy = sites_49pt(Wafer())
    no_aging = Aging(wall_loss=0.0, first_wafer_dip=0.0)
    cols = {}
    for name, (kind, unit, angle, _signed) in SIGNATURES.items():
        def feats(sev):
            if sev == 0:
                cond = Condition()
            else:
                cond, _ = build_condition([FaultEvent(kind, 0, sev, angle_deg=angle)], 10, 10, 0, no_aging)
            return _ideal_features(recipe, hw, cond, xy)[features].to_numpy(dtype=float)

        if _signed:
            if kind == "leveling_tilt":
                # Tilt magnitude cannot be negative; flip the angle instead.
                plus = feats(unit)
                cond_m, _ = build_condition([FaultEvent(kind, 0, unit, angle_deg=angle + 180)], 10, 10, 0, no_aging)
                minus = _ideal_features(recipe, hw, cond_m, xy)[features].to_numpy(dtype=float)
            else:
                plus, minus = feats(unit), feats(-unit)
            cols[name] = (plus - minus) / (2 * unit)
        else:
            cols[name] = (feats(unit) - feats(0)) / unit
    return pd.DataFrame(cols, index=features)


# ---------------------------------------------------------------------------
# 3. Constrained least squares (Lawson-Hanson NNLS)
# ---------------------------------------------------------------------------
def nnls(A, b, max_iter: int = 200):
    """Solve min ||A x - b|| subject to x >= 0."""
    m, n = A.shape
    x = np.zeros(n)
    passive = np.zeros(n, dtype=bool)
    w = A.T @ (b - A @ x)
    for _ in range(max_iter):
        if passive.all() or np.all(w[~passive] <= 1e-10):
            break
        j = np.argmax(np.where(passive, -np.inf, w))
        passive[j] = True
        while True:
            z = np.zeros(n)
            z[passive], *_ = np.linalg.lstsq(A[:, passive], b, rcond=None)
            if np.all(z[passive] > 0):
                x = z
                break
            mask = passive & (z <= 0)
            alpha = np.min(x[mask] / (x[mask] - z[mask]))
            x = x + alpha * (z - x)
            passive &= x > 1e-12
            x[~passive] = 0.0
        w = A.T @ (b - A @ x)
    return x


# ---------------------------------------------------------------------------
# The diagnoser
# ---------------------------------------------------------------------------
@dataclasses.dataclass
class Diagnoser:
    reference: ReferenceModel
    signatures: pd.DataFrame
    min_significance: float = 6.0  # a fault must explain at least ~6 sigma of deviation
    relative_floor: float = 0.05  # ...and at least 5 % of the strongest fault's contribution
    run_rule: tuple = (3, 5)  # alarm when 3 of the last 5 wafers exceed the T2 limit
    min_explained: float = 0.7  # below this, the signature library can't explain the wafer

    @classmethod
    def train(cls, healthy_features: pd.DataFrame, recipe: Recipe = None, hw: Hardware = None, **kw):
        ref = ReferenceModel().fit(healthy_features)
        sig = build_signatures(recipe, hw, ref.features)
        return cls(ref, sig, **kw)

    def _whiten(self, M):
        return np.linalg.solve(self.reference.chol, M)

    def decompose(self, residual_row):
        """Severity estimate per signature for one wafer's residual vector."""
        S = self.signatures.to_numpy()
        names = list(self.signatures.columns)
        signed = np.array([SIGNATURES[n][3] for n in names])
        # Split signed signatures into +/- columns so everything is non-negative.
        A = np.column_stack([S, -S[:, signed]])
        Aw = self._whiten(A)
        rw = self._whiten(residual_row)
        x = nnls(Aw, rw)
        est = x[: len(names)].copy()
        est[signed] -= x[len(names):]
        fitted = self._whiten(S) @ est
        explained = 1.0 - np.sum((rw - fitted) ** 2) / max(np.sum(rw ** 2), 1e-12)
        Sw = self._whiten(S)
        contrib = {}
        for group, members in FAULT_GROUPS.items():
            idx = [names.index(m) for m in members]
            contrib[group] = float(np.linalg.norm(Sw[:, idx] @ est[idx]))
        return dict(zip(names, est)), contrib, float(explained)

    def run(self, feats: pd.DataFrame) -> pd.DataFrame:
        """Detect and diagnose every wafer. Returns one row per wafer."""
        ref = self.reference
        t2 = ref.t2(feats)
        raw = t2 > ref.t2_limit
        k, n = self.run_rule
        roll = pd.Series(raw.astype(int)).rolling(n, min_periods=1).sum().to_numpy()
        alarm = roll >= k
        R = ref.residuals(feats)

        rows = []
        for i in range(len(feats)):
            est, contrib, explained = self.decompose(R[i])
            top = max(contrib.values()) if contrib else 0.0
            present = [g for g, c in sorted(contrib.items(), key=lambda kv: -kv[1])
                       if c >= self.min_significance and c >= self.relative_floor * top]
            if not raw[i]:
                label = "none"
            elif not present or explained < self.min_explained:
                label = "unknown"
            else:
                label = "+".join(sorted(present))
            row = {"wafer_index": int(feats["wafer_index"].iloc[i]), "t2": float(t2[i]),
                   "t2_limit": ref.t2_limit, "out_of_limit": bool(raw[i]), "alarm": bool(alarm[i]),
                   "explained": explained, "diagnosis": label}
            row.update({f"est_{g}": v for g, v in est.items()})
            row.update({f"contrib_{g}": v for g, v in contrib.items()})
            rows.append(row)
        out = pd.DataFrame(rows)
        out["est_leveling_tilt"] = np.hypot(out["est_leveling_tilt_x"], out["est_leveling_tilt_y"])
        if "gt_label" in feats.columns:
            out["gt_label"] = feats["gt_label"].to_numpy() if "gt_label" in feats else None
        return out


# ---------------------------------------------------------------------------
# Scoring against ground truth
# ---------------------------------------------------------------------------
def _canonical(label: str) -> str:
    return "+".join(sorted(label.split("+"))) if label not in ("none", "unknown") else label


def score(results: pd.DataFrame, gt_labels, fault_start: int = None) -> dict:
    gt = pd.Series(list(gt_labels)).map(_canonical)
    pred = results["diagnosis"].map(_canonical)
    healthy = gt == "none"
    faulty = ~healthy
    out = {
        "false_alarm_rate": float(results.loc[healthy.values, "alarm"].mean()) if healthy.any() else np.nan,
        "detected": bool(results.loc[faulty.values, "alarm"].any()) if faulty.any() else False,
    }
    if faulty.any() and fault_start is not None:
        hits = results[(results.wafer_index >= fault_start) & results.alarm]
        out["detect_wafer"] = int(hits.wafer_index.iloc[0]) if len(hits) else None
        out["delay_wafers"] = (out["detect_wafer"] - fault_start) if len(hits) else None
    flagged = faulty.values & results["out_of_limit"].values
    out["diagnosis_accuracy"] = float((pred[flagged] == gt[flagged]).mean()) if flagged.any() else np.nan
    out["n_flagged_faulty"] = int(flagged.sum())
    return out


def shewhart_detection(series_train, series_test, wafer_index, start, run_rule=(3, 5)):
    """Classic single-variable 3-sigma chart, no context correction.

    Uses the same run rule as the multivariate method so the comparison is
    fair. Returns (first alarm wafer at or after `start`, alarm count before
    `start`).
    """
    mu, sd = float(np.mean(series_train)), float(np.std(series_train, ddof=1))
    x = np.asarray(series_test)
    out = ((x > mu + 3 * sd) | (x < mu - 3 * sd)).astype(int)
    k, n = run_rule
    alarm = pd.Series(out).rolling(n, min_periods=1).sum().to_numpy() >= k
    idx = np.asarray(wafer_index)
    hits = idx[(idx >= start) & alarm]
    false_alarms = int(alarm[idx < start].sum())
    return (int(hits[0]) if len(hits) else None), false_alarms
