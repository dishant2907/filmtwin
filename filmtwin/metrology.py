"""Module 2: a simulated spectral reflectometer and the fitting engine behind it.

Two halves, the same as a real tool:

1. MEASURE  `Reflectometer.measure()` turns a true film (thickness, index)
   into a noisy reflectance spectrum, the raw signal the tool records.
2. FIT      `Reflectometer.fit()` recovers thickness and refractive index
   from that spectrum by matching it to the optical model:
     a. grid search over thickness (robust to the many local minima),
     b. Levenberg-Marquardt refinement of thickness, Cauchy A and a
        reflectance scale factor (absorbs reference calibration error),
     c. uncertainty and parameter correlation from the fit covariance.

Every real reflectometer or ellipsometer does some version of this. Knowing
where it breaks (wrong minimum, thin-film n-d correlation, wrong stack
model) is what separates someone who reads the number off the screen from
someone who can trust it.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .optics import cauchy_A_from_n633, cauchy_n, silicon_nk, stack_reflectance, CAUCHY_B_SIO2


@dataclass
class FitResult:
    thickness_nm: float
    n633: float
    scale: float
    gof: float  # goodness of fit, 1 = perfect
    sigma_thickness_nm: float
    sigma_n633: float
    corr_d_n: float  # correlation between thickness and index estimates
    sse: float
    iterations: int
    ok: bool
    message: str = ""


@dataclass
class Reflectometer:
    """Normal-incidence spectral reflectometer measuring a film on silicon."""

    wl_min_nm: float = 400.0
    wl_max_nm: float = 1000.0
    n_points: int = 301
    noise_sigma: float = 0.0015  # absolute reflectance noise per point
    calibration_sigma: float = 0.003  # per-measurement reference/scale error
    under_layers: list = field(default_factory=list)  # [(N, d_nm), ...] between film and Si
    seed: int = 0

    def __post_init__(self):
        self.wl = np.linspace(self.wl_min_nm, self.wl_max_nm, self.n_points)
        self.n_si = silicon_nk(self.wl)
        self.rng = np.random.default_rng(self.seed)

    # -- forward model -------------------------------------------------------
    def model(self, thickness_nm, A, scale=1.0, under_layers=None):
        """Model reflectance for a Cauchy film of the given thickness and A."""
        ul = self.under_layers if under_layers is None else under_layers
        N_film = cauchy_n(self.wl, A).astype(complex)
        d = np.asarray(thickness_nm, dtype=float)
        if d.ndim == 1:
            d = d[:, None]
        R = stack_reflectance(self.wl, [(N_film, d)] + list(ul), self.n_si)
        return scale * R

    # -- measurement ---------------------------------------------------------
    def measure(self, thickness_nm, n633, under_layers=None):
        """Simulated raw spectrum for a true film. Returns R (n_points,)."""
        A = cauchy_A_from_n633(n633)
        R = self.model(thickness_nm, A, under_layers=under_layers)
        scale = 1.0 + self.rng.normal(0.0, self.calibration_sigma)
        return R * scale + self.rng.normal(0.0, self.noise_sigma, size=R.shape)

    # -- fitting -------------------------------------------------------------
    def thickness_landscape(self, R_meas, d_grid, n633_guess=1.47, stride=1):
        """SSE vs thickness with index fixed and scale solved analytically."""
        A = cauchy_A_from_n633(n633_guess)
        idx = slice(None, None, stride)
        wl, nsi = self.wl[idx], self.n_si[idx]
        N_film = cauchy_n(wl, A).astype(complex)
        d = np.asarray(d_grid, dtype=float)[:, None]
        under = [(np.asarray(N)[idx] if np.ndim(N) else N, t) for N, t in self.under_layers]
        Rm = stack_reflectance(wl, [(N_film, d)] + under, nsi)
        Rx = np.asarray(R_meas)[idx]
        s = np.sum(Rm * Rx, axis=1) / np.sum(Rm * Rm, axis=1)
        return np.sum((Rx - s[:, None] * Rm) ** 2, axis=1)

    def fit(self, R_meas, d_min=10.0, d_max=2000.0, d_step=3.0, n633_guess=1.47,
            initial_thickness=None, fit_index="auto", max_iter=40):
        """Fit thickness (and optionally index) to a measured spectrum.

        fit_index: True fits thickness + index; False fixes the index at
        n633_guess; "auto" fits both, then falls back to fixed index if the
        index is not determinable (thin films, where n and d trade off).
        initial_thickness: if given, SKIP the grid search and start the local
        optimizer there. Useful to demonstrate how a bad starting point lands
        in the wrong minimum.
        """
        if fit_index == "auto":
            res = self.fit(R_meas, d_min, d_max, d_step, n633_guess, initial_thickness,
                           fit_index=True, max_iter=max_iter)
            if res.sigma_n633 > 0.005 or not np.isfinite(res.sigma_n633):
                fixed = self.fit(R_meas, d_min, d_max, d_step, n633_guess, initial_thickness,
                                 fit_index=False, max_iter=max_iter)
                fixed.message = ("index not determinable for this film (strong n-d correlation); "
                                 f"index fixed at {n633_guess:.4f}")
                return fixed
            return res

        if initial_thickness is None:
            # Coarse grid on every 2nd wavelength is plenty to find the right basin.
            grid = np.arange(d_min, d_max + d_step, d_step)
            sse = self.thickness_landscape(R_meas, grid, n633_guess, stride=2)
            d0 = float(grid[np.argmin(sse)])
        else:
            d0 = float(initial_thickness)

        A0 = cauchy_A_from_n633(n633_guess)
        p = np.array([d0, A0, 1.0]) if fit_index else np.array([d0, 1.0])

        def unpack(p):
            if fit_index:
                return p[0], p[1], p[2]
            return p[0], A0, p[1]

        def resid(p):
            d, A, s = unpack(p)
            return R_meas - self.model(d, A, s).ravel()

        steps = np.array([0.05, 1e-5, 1e-5]) if fit_index else np.array([0.05, 1e-5])
        lam = 1e-3
        r = resid(p)
        sse = float(r @ r)
        it = 0
        for it in range(1, max_iter + 1):
            J = np.empty((r.size, p.size))
            for j in range(p.size):
                dp = np.zeros_like(p)
                dp[j] = steps[j]
                J[:, j] = -(resid(p + dp) - r) / steps[j]
            JtJ = J.T @ J
            g = J.T @ r
            improved = False
            for _ in range(10):
                H = JtJ + lam * np.diag(np.diag(JtJ))
                try:
                    delta = np.linalg.solve(H, g)
                except np.linalg.LinAlgError:
                    lam *= 10
                    continue
                p_new = p + delta
                r_new = resid(p_new)
                sse_new = float(r_new @ r_new)
                if sse_new < sse:
                    p, r, sse = p_new, r_new, sse_new
                    lam = max(lam / 3, 1e-7)
                    improved = True
                    break
                lam *= 10
            if not improved or np.all(np.abs(delta) < np.array([1e-4, 1e-8, 1e-8])[: p.size]):
                break

        # Uncertainty from the local curvature: cov = s^2 (J^T J)^-1
        dof = max(r.size - p.size, 1)
        s2 = sse / dof
        J = np.empty((r.size, p.size))
        for j in range(p.size):
            dp = np.zeros_like(p)
            dp[j] = steps[j]
            J[:, j] = -(resid(p + dp) - r) / steps[j]
        try:
            cov = s2 * np.linalg.inv(J.T @ J)
        except np.linalg.LinAlgError:
            cov = np.full((p.size, p.size), np.nan)

        d, A, s = unpack(p)
        n633 = float(cauchy_n(633.0, A))
        sig_d = float(np.sqrt(max(cov[0, 0], 0)))
        if fit_index:
            sig_n = float(np.sqrt(max(cov[1, 1], 0)))  # dn633/dA = 1
            corr = float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])) if cov[0, 0] > 0 and cov[1, 1] > 0 else np.nan
        else:
            sig_n, corr = 0.0, 0.0
        sst = float(np.sum((R_meas - R_meas.mean()) ** 2))
        gof = 1.0 - sse / sst if sst > 0 else 0.0
        ok = gof > 0.99 and d_min <= d <= d_max and 1.3 < n633 < 2.2
        msg = "" if ok else "low goodness of fit or parameters out of range: check the stack model"
        return FitResult(float(d), n633, float(s), float(gof), sig_d, sig_n, corr, sse, it, ok, msg)


def gauge_study(meter: Reflectometer, thickness_nm: float, n633: float, repeats: int = 30,
                tolerance_nm: float = None):
    """Repeat the same measurement and report repeatability (a Gage R&R building block).

    tolerance_nm: full spec width (USL - LSL). If given, also returns the
    precision-to-tolerance ratio P/T = 6 * sigma / tolerance. Rule of thumb:
    P/T < 0.1 is good, 0.1-0.3 acceptable, > 0.3 the gauge is not capable.
    """
    fits = [meter.fit(meter.measure(thickness_nm, n633)) for _ in range(repeats)]
    d = np.array([f.thickness_nm for f in fits])
    n = np.array([f.n633 for f in fits])
    out = {
        "true_thickness_nm": thickness_nm,
        "true_n633": n633,
        "mean_thickness_nm": float(d.mean()),
        "bias_nm": float(d.mean() - thickness_nm),
        "repeatability_nm": float(d.std(ddof=1)),
        "mean_n633": float(n.mean()),
        "bias_n633": float(n.mean() - n633),
        "repeatability_n633": float(n.std(ddof=1)),
        "all_ok": all(f.ok for f in fits),
        "thickness_values": d,
        "n633_values": n,
    }
    if tolerance_nm:
        out["p_to_t"] = 6.0 * out["repeatability_nm"] / tolerance_nm
    return out


def measure_wafer(meter: Reflectometer, site_thickness_nm, site_n633, **fit_kwargs):
    """Measure every site of one wafer. Returns a list of FitResult."""
    return [meter.fit(meter.measure(d, n), **fit_kwargs)
            for d, n in zip(site_thickness_nm, site_n633)]
