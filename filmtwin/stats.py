"""Small statistics toolkit (numpy only): t and F distributions, regression with ANOVA.

Written out by hand so the project runs with numpy alone, and so every
number in a DOE table can be traced to a formula you can explain.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import lgamma

import numpy as np


def _betacf(a, b, x, max_iter=300, eps=3e-14):
    """Continued fraction for the incomplete beta function (Numerical Recipes)."""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c if abs(c) > 1e-300 else 1e300
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > 1e-300 else 1e-300)
        c = 1.0 + aa / c if abs(c) > 1e-300 else 1e300
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betainc(a, b, x):
    """Regularised incomplete beta I_x(a, b)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * np.log(x) + b * np.log(1 - x)
    bt = np.exp(lbt)
    if x < (a + 1) / (a + b + 2):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1 - x) / b


def t_sf2(t: float, df: float) -> float:
    """Two-sided p-value for a t statistic."""
    x = df / (df + t * t)
    return float(betainc(df / 2.0, 0.5, x))


def f_sf(f: float, d1: float, d2: float) -> float:
    """Upper-tail p-value for an F statistic."""
    if f <= 0:
        return 1.0
    x = d2 / (d2 + d1 * f)
    return float(betainc(d2 / 2.0, d1 / 2.0, x))


def t_ppf(p: float, df: float) -> float:
    """Inverse t CDF (two-sided critical value for p in (0.5, 1)) by bisection."""
    lo, hi = 0.0, 1000.0
    target = 2 * (1 - p)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if t_sf2(mid, df) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass
class OLSResult:
    terms: list
    coef: np.ndarray
    se: np.ndarray
    t: np.ndarray
    p: np.ndarray
    r2: float
    r2_adj: float
    rmse: float
    df_resid: int
    f_stat: float
    f_p: float
    XtX_inv: np.ndarray

    def table(self):
        import pandas as pd
        return pd.DataFrame({"term": self.terms, "estimate": self.coef, "std_error": self.se,
                             "t": self.t, "p_value": self.p})

    def predict(self, X, interval: float = None, n_new: int = 1):
        """Prediction, optionally with a prediction-interval half-width.

        n_new: number of new runs whose MEAN will be compared with the
        prediction (confirmation runs). n_new=1 is the interval for a single
        new wafer; larger n_new narrows the noise part, sqrt(1/n_new + leverage).
        """
        X = np.atleast_2d(X)
        yhat = X @ self.coef
        if interval is None:
            return yhat
        tcrit = t_ppf(interval, self.df_resid)
        lev = np.einsum("ij,jk,ik->i", X, self.XtX_inv, X)
        half = tcrit * self.rmse * np.sqrt(1.0 / n_new + lev)
        return yhat, half


def ols(X, y, terms):
    """Ordinary least squares with standard errors, t-tests, R^2 and overall F-test."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, k = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    coef = XtX_inv @ X.T @ y
    resid = y - X @ coef
    df = n - k
    sse = float(resid @ resid)
    sst = float(np.sum((y - y.mean()) ** 2))
    mse = sse / df if df > 0 else np.nan
    se = np.sqrt(np.clip(np.diag(XtX_inv) * mse, 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = coef / se
    p = np.array([t_sf2(abs(tt), df) if np.isfinite(tt) and df > 0 else np.nan for tt in t])
    r2 = 1 - sse / sst if sst > 0 else np.nan
    r2_adj = 1 - (1 - r2) * (n - 1) / df if df > 0 else np.nan
    ssr = sst - sse
    f_stat = (ssr / (k - 1)) / mse if (k > 1 and df > 0) else np.nan
    f_p = f_sf(f_stat, k - 1, df) if np.isfinite(f_stat) else np.nan
    return OLSResult(list(terms), coef, se, t, p, r2, r2_adj, float(np.sqrt(mse)), df, f_stat, f_p, XtX_inv)
