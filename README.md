# FilmTwin

A physics-based digital twin of a single-wafer PECVD chamber, plus the analysis tools a thin-film process engineer uses every day: wafer maps, uniformity metrics, fault signatures, SPC, DOE and root-cause analysis.

FilmTwin answers one question: **when a wafer's thickness map goes out of spec, what in the chamber caused it, and what should be changed?**

Real fab data is proprietary, so almost no public datasets come with a known root cause. FilmTwin generates realistic wafer maps and tool sensor data for chambers that age and develop faults, with the true cause recorded for every wafer. That makes it possible to build and *score* diagnosis methods.

![Fault gallery](outputs/figures/02_fault_gallery.png)

## Status

| Module | What it does | Status |
|---|---|---|
| 1. Chamber twin | Grove-model deposition, radial profiles, aging, 6 fault types, sensors, datasets | **Done (v0.1)** |
| 2. Metrology engine | Transfer-matrix reflectance fitting for thickness and n, FTIR bond analysis | Next |
| 3. Root-cause diagnosis | Zernike decomposition of wafer maps, signature matching, scored against ground truth | Planned |
| 4. DOE, SPC, run-to-run | Factorial/RSM designs, control charts with OCAP, EWMA run-to-run control | Planned |
| 5. Virtual metrology | VM model on the public PHM 2016 CMP dataset | Planned |

## Quick start (Mac Terminal)

Requires Python 3.8 or newer.

```bash
cd filmtwin
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python3 tests/test_module1.py          # 11 physics and simulator checks
python3 examples/run_module1_demo.py   # figures + datasets (~10 s)
```

Figures land in `outputs/figures/`, datasets in `data/` (about 32 MB, so `data/` is git-ignored; anyone can regenerate it with the demo script).

## Project layout

```
filmtwin/
  physics.py    Grove model, temperature/plasma/gas profiles, film properties
  faults.py     Fault definitions and normal chamber aging
  chamber.py    Simulator: wafers, lots, cleans, sensors, metrology tables
  scenarios.py  8 ready-made fault stories with known root causes
  regime.py     Reaction- vs transport-limited analysis, knob sensitivities
  metrics.py    Uniformity metrics (1-sigma and range NU), radial profiles
  plots.py      Wafer maps, fault gallery, run charts, Grove curve
  wafer.py      300 mm wafer geometry, 49-site pattern, dense grids
examples/run_module1_demo.py
tests/test_module1.py
```

## The physics

At every point on the wafer the deposition rate follows the Grove model for CVD:

```
G = eta * V_m * [ k_s * h_g / (k_s + h_g) ] * C_g * g(r)
```

- `k_s` surface reaction rate constant: Arrhenius in local wafer temperature, scaled by RF-generated radical supply (saturating with power) and a radial plasma profile.
- `h_g = D / delta` gas-phase mass transfer coefficient. Diffusivity `D` scales as `T^1.75 / P`; boundary layer `delta` scales with showerhead-to-wafer spacing.
- `C_g` precursor concentration from partial pressure and gas temperature.
- `g(r)` radial precursor depletion: consumption over the wafer compared with gas supply, so more total flow means less depletion.
- `V_m`, `eta` film molar volume and incorporation efficiency.

The **Damkohler number** `Da = k_s / h_g` sets the regime. With `Da << 1` the process is reaction-limited and rate follows temperature. With `Da >> 1` it is transport-limited and rate follows flow, pressure and spacing. The default PECVD recipe sits at `Da ≈ 0.5` (mixed). `regime.py` contrasts it with a hypothetical thermal CVD chamber that shows the full textbook curve.

![Grove curve](outputs/figures/03_regime_grove_curve.png)

Spatial non-uniformity comes from:

- two heater zones plus edge heat loss,
- a centre-high plasma term and an edge-ring term,
- radial gas depletion,
- planar tilt,
- a small random field, wafer-to-wafer noise and gauge noise.

Film properties include refractive index (rises as the film turns Si-rich) and wet etch rate ratio (falls with deposition temperature). These are illustrative trends with plausible magnitudes.

Default recipe (SiH4 + N2O PECVD oxide): 400 C, 2.5 Torr, 600 W, 150/3000 sccm, 10 mm gap, 103 s. This gives about 500 nm at ~290 nm/min with 1.3 % 1-sigma non-uniformity.

## Scenarios

Each scenario is 600 wafers (24 lots of 25) with a chamber clean every 300 wafers.

| Scenario | Map signature | Sensor clue |
|---|---|---|
| `baseline` | Stable centre-thick bowl | Vpp creeps up between cleans |
| `edge_tc_drift` | Edge ~3 % thin, edge WERR up | Edge heater power falls |
| `showerhead_clog` | Centre thins; NU *improves* first, then worsens | None |
| `mfc_drift` | Mean shift, no shape change, RI drops | None (MFC reports setpoint) |
| `rf_match_degradation` | Mean shift, no shape change, RI flat | Reflected power up, Vpp down |
| `edge_ring_erosion` | Edge goes thin to thick; NU dips then rises | Weak |
| `leveling_tilt` | Step change after clean; one side thick | None (timing is the clue) |
| `mixed_clog_and_match` | Mean drop plus shape change | Reflected power explains only part |

Two lessons are built in on purpose. First, `mfc_drift` and `rf_match_degradation` produce the same thickness map, so you need refractive index or sensor data to tell them apart. Second, in `showerhead_clog` and `edge_ring_erosion` the 1-sigma NU looks *better* for a while, so an NU-only SPC chart misses a real fault.

## Data dictionary

`data/<scenario>_wafers.csv` has one row per wafer:

- **Context:** `wafer_id`, `lot_id`, `slot`, `chamber`, `timestamp`, `wafers_since_clean`, `rf_hours_since_clean`, `pm_event`, `sampled_for_metrology`
- `r_*` the recipe used
- `s_*` tool sensors as the tool reports them (pressure, throttle valve, RF forward/reflected, Vpp, heater powers, precursor flow reading, wall temperature)
- `m_*` metrology summaries (mean thickness, NU 1-sigma and range, RI, WERR, rate)
- `gt_*` **ground truth**: severity of each fault and a combined `gt_label`. A real fab never gives you these; use them only to score your methods.

`data/<scenario>_sites.csv` has one row per measured site: `site`, `x_mm`, `y_mm`, `thickness_nm`, `refractive_index`, `werr`, `temperature_C`.

## Limitations (stated plainly)

- The model is physically motivated but **not calibrated to any commercial tool**. Constants were set so the default recipe lands in a realistic window. Absolute numbers will differ from a real chamber; the directions and the relative signatures are the point.
- Film-property relationships (RI, WERR) are qualitative trends, not a validated film model.
- The gas model is 1D radial. There is no CFD, plasma chemistry or sheath physics.
- Faults are idealized single mechanisms. Real excursions are messier and often overlap, which is why `mixed_clog_and_match` exists.

Module 5 adds a real public dataset so the analysis methods are also proven on real fab data.

## Using it in your own code

```python
from filmtwin import ChamberSimulator, FaultEvent, Recipe, Hardware

sim = ChamberSimulator(
    recipe=Recipe(pressure_torr=3.0),
    faults=[FaultEvent("showerhead_clog", start_wafer=100, severity=0.08, ramp_wafers=300)],
    seed=1,
)
wafers, sites = sim.run(400)
print(wafers[["wafer_id", "m_thickness_mean_nm", "m_nu_1sigma_pct", "gt_label"]].tail())
```
