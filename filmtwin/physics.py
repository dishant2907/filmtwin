"""Core deposition physics for a single-wafer showerhead PECVD chamber.

The model is the classic Grove model for CVD, applied locally at every point
on the wafer:

    flux to surface      F1 = h_g * (C_g - C_s)      (gas-phase mass transfer)
    flux consumed        F2 = k_s * C_s              (surface reaction)
    steady state F1 = F2 gives
    growth rate          G  = eta * V_m * [k_s * h_g / (k_s + h_g)] * C_g

    k_s   surface reaction rate constant (m/s), Arrhenius in temperature and
          scaled by plasma radical supply for PECVD
    h_g   mass transfer coefficient (m/s) = D / delta (diffusivity over
          boundary-layer thickness)
    C_g   precursor concentration in the gas (mol/m^3)
    V_m   molar volume of the deposited film (m^3/mol)
    eta   incorporation efficiency (fraction of reacted precursor that ends up
          in the film; a calibration knob)

The Damkohler number Da = k_s / h_g tells you the regime:
    Da << 1  reaction-limited: rate tracks temperature (Arrhenius)
    Da >> 1  transport-limited: rate tracks gas flow, pressure and spacing

Spatial non-uniformity comes from radial profiles of temperature, plasma
density and precursor depletion, plus tilt and a small random field.

IMPORTANT: the numbers are physically motivated but NOT calibrated to any
specific commercial tool. They are set so a default recipe lands in a
realistic window (~200-400 nm/min, ~1-2 % 1-sigma non-uniformity, mixed
regime). The film-property trends (refractive index, wet etch rate ratio)
are illustrative empirical forms, documented as such.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

R_GAS = 8.314  # J/(mol K)
P_STD = 101325.0  # Pa
T_STD = 273.15  # K
TORR = 133.322  # Pa per Torr


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------
@dataclass
class Recipe:
    """Process recipe knobs a process engineer would set on the tool."""

    temperature_C: float = 400.0  # heater setpoint
    pressure_torr: float = 2.5
    rf_power_W: float = 600.0
    precursor_sccm: float = 150.0  # e.g. SiH4
    oxidizer_sccm: float = 3000.0  # e.g. N2O
    spacing_mm: float = 10.0  # showerhead-to-wafer gap
    edge_temp_offset_C: float = 0.0  # edge heater zone trim
    dep_time_s: float = 103.0

    @property
    def total_sccm(self) -> float:
        return self.precursor_sccm + self.oxidizer_sccm

    @property
    def flow_ratio(self) -> float:
        """Precursor / oxidizer ratio (controls film stoichiometry)."""
        return self.precursor_sccm / self.oxidizer_sccm


@dataclass
class Hardware:
    """Fixed properties of one physical chamber.

    Two chambers of the same model differ slightly; make that variation with
    `Hardware.sample(rng)` for chamber matching studies.
    """

    name: str = "CH-A"
    # Surface kinetics
    Ea_kJ_mol: float = 8.0  # PECVD: weak temperature dependence (plasma supplies the activation)
    k_s_ref: float = 1.3  # m/s at T_ref and P_rf_ref
    T_ref_K: float = 673.15
    rf_sat_W: float = 900.0  # radical supply saturates with power
    rf_ref_W: float = 600.0
    # Transport
    D_ref_m2_s: float = 1.5e-5  # precursor diffusivity at 300 K, 1 atm
    delta_factor: float = 0.5  # boundary layer = factor * spacing
    showerhead_C: float = 150.0
    depletion_phi: float = 0.015  # fraction of flow that travels radially
    # Film
    molar_volume: float = 2.73e-5  # m^3/mol for SiO2
    eta: float = 0.06
    # Spatial signatures built into the hardware
    heater_split_rho: float = 0.70  # where centre zone hands over to edge zone
    heater_split_width: float = 0.06
    edge_heat_loss_C: float = 4.0  # wafer edge runs cooler
    plasma_center_peak: float = 0.008  # standing-wave style centre-high term
    plasma_edge_amp: float = 0.012  # edge-ring enhancement
    plasma_edge_width: float = 0.04
    tilt: float = 0.002  # residual leveling tilt (fraction per radius)
    tilt_angle_deg: float = 30.0
    # Noise
    wafer_to_wafer_sigma: float = 0.003  # fractional rate noise per wafer
    spatial_noise: float = 0.002  # fractional smooth random field
    gauge_sigma_nm: float = 0.4  # metrology repeatability
    ri_sigma: float = 0.0008  # refractive index repeatability
    werr_rel_sigma: float = 0.015  # wet etch test scatter (fractional)

    @staticmethod
    def sample(rng: np.random.Generator, name: str = "CH-X") -> "Hardware":
        """A chamber drawn from realistic tool-to-tool spread."""
        h = Hardware(name=name)
        h.k_s_ref *= rng.normal(1.0, 0.03)
        h.depletion_phi *= rng.normal(1.0, 0.10)
        h.edge_heat_loss_C *= rng.normal(1.0, 0.15)
        h.plasma_center_peak *= rng.normal(1.0, 0.20)
        h.plasma_edge_amp *= rng.normal(1.0, 0.20)
        h.tilt = abs(rng.normal(h.tilt, 0.001))
        h.tilt_angle_deg = rng.uniform(0, 360)
        return h


@dataclass
class Condition:
    """The chamber's *current* state, which drifts with use and faults.

    All values are multipliers or offsets relative to a healthy, freshly
    seasoned chamber. The simulator in `chamber.py` updates these wafer by
    wafer; you can also set them by hand to explore a single signature.
    """

    radical_factor: float = 1.0  # wall buildup / seasoning (k_s multiplier)
    edge_temp_error_C: float = 0.0  # actual minus intended edge temperature
    center_gas_loss: float = 0.0  # showerhead centre holes clogged (0-1)
    precursor_flow_error: float = 0.0  # MFC actual/setpoint - 1
    rf_delivered_factor: float = 1.0  # fraction of forward power delivered
    extra_edge_plasma: float = 0.0  # edge ring erosion
    extra_tilt: float = 0.0  # leveling fault
    extra_tilt_angle_deg: float = 0.0


# --------------------------------------------------------------------------
# Physics
# --------------------------------------------------------------------------
def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def wafer_temperature_K(rho, recipe: Recipe, hw: Hardware, cond: Condition):
    """Local wafer temperature from two heater zones plus edge heat loss."""
    edge_weight = _sigmoid((rho - hw.heater_split_rho) / hw.heater_split_width)
    T_C = (
        recipe.temperature_C
        + (recipe.edge_temp_offset_C + cond.edge_temp_error_C) * edge_weight
        - hw.edge_heat_loss_C * rho ** 4
    )
    return T_C + 273.15


def plasma_profile(rho, hw: Hardware, cond: Condition):
    """Relative radical density across the wafer (1.0 = nominal)."""
    centre = hw.plasma_center_peak * (1.0 - 2.0 * rho ** 2)
    edge_amp = hw.plasma_edge_amp + cond.extra_edge_plasma
    edge = edge_amp * np.exp(-(1.0 - np.clip(rho, 0, 1)) / hw.plasma_edge_width)
    return 1.0 + centre + edge


def rf_factor(power_W: float, hw: Hardware) -> float:
    """Radical generation vs delivered RF power, saturating, = 1 at ref power."""
    f = 1.0 - np.exp(-power_W / hw.rf_sat_W)
    f_ref = 1.0 - np.exp(-hw.rf_ref_W / hw.rf_sat_W)
    return f / f_ref


def gas_state(recipe: Recipe, hw: Hardware, cond: Condition, T_wafer_mean_K: float):
    """Bulk gas properties: precursor concentration, diffusivity, h_g, flow."""
    P = recipe.pressure_torr * TORR
    T_gas = 0.5 * (T_wafer_mean_K + hw.showerhead_C + 273.15)
    precursor_actual = recipe.precursor_sccm * (1.0 + cond.precursor_flow_error)
    total_actual = precursor_actual + recipe.oxidizer_sccm
    y = precursor_actual / total_actual
    C_g = y * P / (R_GAS * T_gas)  # mol/m^3
    D = hw.D_ref_m2_s * (T_gas / 300.0) ** 1.75 * (P_STD / P)
    delta = hw.delta_factor * recipe.spacing_mm * 1e-3
    h_g = D / delta
    # Actual volumetric flow at chamber conditions (m^3/s)
    Q_vol = total_actual * 1e-6 / 60.0 * (P_STD / P) * (T_gas / T_STD)
    return dict(P=P, T_gas=T_gas, y=y, C_g=C_g, D=D, h_g=h_g, Q_vol=Q_vol,
                precursor_actual=precursor_actual, total_actual=total_actual)


def surface_rate_constant(T_K, recipe: Recipe, hw: Hardware, cond: Condition, plasma):
    """k_s(T, plasma) in m/s."""
    arrhenius = np.exp(-hw.Ea_kJ_mol * 1e3 / R_GAS * (1.0 / T_K - 1.0 / hw.T_ref_K))
    power = recipe.rf_power_W * cond.rf_delivered_factor
    return hw.k_s_ref * arrhenius * rf_factor(power, hw) * plasma * cond.radical_factor


@dataclass
class DepositionResult:
    thickness_nm: np.ndarray
    rate_nm_min: np.ndarray
    refractive_index: np.ndarray
    werr: np.ndarray  # wet etch rate ratio vs thermal oxide
    temperature_C: np.ndarray
    damkohler: np.ndarray
    info: dict = field(default_factory=dict)


def deposit(
    xy_mm: np.ndarray,
    recipe: Recipe,
    hw: Hardware = None,
    cond: Condition = None,
    rng: np.random.Generator = None,
    usable_radius_mm: float = 147.0,
    noise: bool = True,
) -> DepositionResult:
    """Deposit one film and return local properties at the given sites.

    xy_mm: (N, 2) site positions. Use `noise=False` for the ideal map.
    """
    hw = hw or Hardware()
    cond = cond or Condition()
    rng = rng or np.random.default_rng()

    x, y = xy_mm[:, 0], xy_mm[:, 1]
    rho = np.sqrt(x ** 2 + y ** 2) / usable_radius_mm
    theta = np.arctan2(y, x)

    T = wafer_temperature_K(rho, recipe, hw, cond)
    plasma = plasma_profile(rho, hw, cond)
    gas = gas_state(recipe, hw, cond, float(np.mean(T)))
    k_s = surface_rate_constant(T, recipe, hw, cond, plasma)
    h_g = gas["h_g"]

    k_eff = k_s * h_g / (k_s + h_g)

    # Radial precursor depletion: consumption over the wafer vs supply.
    area = np.pi * (usable_radius_mm * 1e-3) ** 2
    beta = hw.depletion_phi * float(np.mean(k_eff)) * area / gas["Q_vol"]
    gas_profile = np.exp(-beta * rho ** 2)
    # Clogged centre showerhead holes starve the centre of the wafer.
    gas_profile = gas_profile * (1.0 - cond.center_gas_loss * np.exp(-(rho / 0.35) ** 2))

    rate_m_s = hw.eta * hw.molar_volume * k_eff * gas["C_g"] * gas_profile

    # Geometric tilt (chamber leveling) as a first-order planar term.
    tilt = hw.tilt * rho * np.cos(theta - np.deg2rad(hw.tilt_angle_deg))
    tilt += cond.extra_tilt * rho * np.cos(theta - np.deg2rad(cond.extra_tilt_angle_deg))
    rate_m_s = rate_m_s * (1.0 + tilt)

    if noise:
        rate_m_s = rate_m_s * rng.normal(1.0, hw.wafer_to_wafer_sigma)
        rate_m_s = rate_m_s * (1.0 + smooth_random_field(rho, theta, hw.spatial_noise, rng))

    rate_nm_min = rate_m_s * 1e9 * 60.0
    thickness = rate_m_s * 1e9 * recipe.dep_time_s
    if noise:
        thickness = thickness + rng.normal(0.0, hw.gauge_sigma_nm, size=thickness.shape)

    n, werr = film_properties(T, recipe, cond, plasma)
    if noise:
        # Ellipsometer RI repeatability and wet-etch test scatter.
        n = n + rng.normal(0.0, hw.ri_sigma, size=n.shape)
        werr = werr * rng.normal(1.0, hw.werr_rel_sigma, size=werr.shape)

    info = dict(gas, beta=beta, k_s_mean=float(np.mean(k_s)), k_eff_mean=float(np.mean(k_eff)))
    return DepositionResult(
        thickness_nm=thickness,
        rate_nm_min=rate_nm_min,
        refractive_index=n,
        werr=werr,
        temperature_C=T - 273.15,
        damkohler=k_s / h_g,
        info=info,
    )


def film_properties(T_K, recipe: Recipe, cond: Condition, plasma):
    """Illustrative film-property trends for PECVD SiO2 (SiH4 + N2O).

    - Refractive index rises as the film becomes Si-rich (higher SiH4/N2O).
      Stoichiometric thermal oxide is ~1.46.
    - Wet etch rate ratio (vs thermal oxide) falls with deposition
      temperature: hotter films are denser and contain less hydrogen.
    These are qualitative trends with plausible magnitudes, not a calibrated
    film model. Replace with fitted relationships once you have real data.
    """
    ratio = recipe.flow_ratio * (1.0 + cond.precursor_flow_error)
    T_C = T_K - 273.15
    n = 1.458 + 0.8 * np.maximum(ratio - 0.025, 0.0) - 0.00002 * (T_C - 400.0)
    n = np.clip(n, 1.45, 1.60) * np.ones_like(T_K)
    power_term = (600.0 / max(recipe.rf_power_W * cond.rf_delivered_factor, 50.0)) ** 0.2
    werr = 3.0 * np.exp(-(T_C - 300.0) / 150.0) * power_term / np.sqrt(plasma)
    return n, werr


def smooth_random_field(rho, theta, amplitude, rng, n_modes: int = 6):
    """Low-order random spatial variation (what's left after the systematic terms)."""
    field = np.zeros_like(rho, dtype=float)
    for _ in range(n_modes):
        m = rng.integers(0, 4)
        k = rng.uniform(0.5, 3.0)
        phase_r, phase_t = rng.uniform(0, 2 * np.pi, size=2)
        field += np.cos(k * np.pi * rho + phase_r) * np.cos(m * theta + phase_t)
    return amplitude * field / np.sqrt(n_modes / 2.0)
