# Modules 2 and 3: what you should be able to explain

## Module 2: metrology

### 1. "How does a reflectometer measure film thickness?"

Light reflects off the top of the film and off the film/silicon interface. The two reflections interfere. Whether they add or cancel depends on the extra path the second one travels, 2 × n × d, compared with the wavelength. Sweep the wavelength and the reflectance oscillates. The tool fits a physical model of the stack to that spectrum, using the transfer-matrix method, and reads off the d and n that make the model match.

Show `20_spectrum_fit.png`.

### 2. "Why can the fit give a wrong answer?"

Three reasons. Each one has a figure.

- **Wrong valley** (`21_thickness_landscape.png`). Fit error vs thickness has many local minima, roughly one per interference order. A local optimizer started in the wrong one converges there. Real tools search globally first, or start from the recipe's target thickness.
- **Thin films** (`22_uncertainty_vs_thickness.png`). Below about 50 nm there isn't even one full fringe, so a thinner film with a higher index looks almost identical to a thicker one with a lower index. Tools fix n for thin films. If the fixed n is wrong, every thickness reading carries a bias of about the same percentage.
- **Wrong stack** (`25_wrong_stack_model.png`). The model assumes oxide directly on Si, but there's nitride underneath, so the answer is wrong. The goodness of fit drops, and that's the alarm.

The lesson: always look at GOF, not just the thickness number.

### 3. "What is a gauge study and why does it matter?"

You measure the same thing many times to see how much the gauge itself varies (repeatability). Then you compare that variation to the spec width using precision-to-tolerance: P/T = 6σ / (USL − LSL).

Here, P/T is 0.03, so the gauge uses only 3 % of the spec window. A common rule of thumb: below 0.1 is good, above 0.3 the gauge can't tell good wafers from bad. A full Gage R&R adds reproducibility: different operators, tools and days.

Tie-in: Module 1 assumed 0.4 nm gauge noise, and Module 2 measures 0.15 nm repeatability. The gap is reproducibility margin.

## Module 3: detection and diagnosis

### 4. "Why not just put an SPC chart on uniformity?"

- **Faults that make NU look better.** A clogged showerhead first makes uniformity *better*. If your NU chart only has an upper limit, which is common because "lower is better," you miss it.
- **Faults that don't change shape.** MFC drift and RF match faults barely change NU at all. The NU chart missed both.
- **No root cause.** A single chart tells you *that* something moved, never *what*.

Show `34_method_comparison.png`.

### 5. "What is Hotelling's T-squared?"

It's a multivariate distance from "healthy." It combines many features into one number, weighted by how they normally vary *together*. A small move in two features that normally move in opposite directions can be more suspicious than a big move in one noisy feature.

It's the backbone of many fab fault detection (FDC) systems.

**Be honest about the weakness:** pooling 15 features dilutes a shift in just one. That's why the plain mean-thickness chart caught MFC drift earlier. Real systems run both kinds of chart.

### 6. "How does your diagnosis work? Is it machine learning?"

It's physics-informed, not a black-box classifier. For each known fault, I compute from the chamber model what one small unit of that fault does to every feature. That's the fingerprint table, `31_signature_table.png`.

For a suspicious wafer, I solve for the non-negative mix of fingerprints that best explains its deviation. This is constrained least squares, weighted by the healthy covariance. That gives three things a classifier doesn't:

- two faults at once (`35_case_mixed.png`),
- a severity estimate for each fault (`33_severity_tracking.png`),
- "unknown" when nothing in the library explains the data (`36_case_unseen_fault.png`).

### 7. "Why do you correct for wafers-since-clean?"

The deposition rate drifts down between chamber cleans as the walls coat up, and the first wafer after a clean always runs low. Both are normal. Without the context model, either every PM cycle looks like a fault, or the limits get so wide that real faults slip through.

### 8. "How do you know it would work on a real tool?" (be ready for this)

Honest answer: I don't yet. It was tested on simulated data from the same physics that made the signatures, and that's the biggest caveat.

To reduce it, I tested on three chambers with different hardware than the signatures were built from (about 98 % accuracy) and on a fault type not in the library (correctly labelled "unknown"). On a real tool, I'd build signatures from a mix of physics, past incident data and engineer judgement, and validate them on historical excursions with known root causes.

That answer shows you understand validation, and interviewers value it more than a perfect score.
