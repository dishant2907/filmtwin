# Module 4: what you should be able to explain

## DOE

### 1. "Walk me through how you'd optimise a recipe."

Four stages, each with a figure.

1. **Screen.** A 2^(6-2) fractional factorial tests 6 knobs in 16 runs plus 4 centre points, instead of 64 runs for a full factorial. It finds the vital few (`40_doe_screening_pareto.png`).
2. **Model.** A face-centred central composite design on those 3 knobs fits curved response surfaces.
3. **Choose.** Desirability trades rate against uniformity against film quality, and the contour overlay shows the process window (`41_doe_process_window.png`).
4. **Confirm.** Run the chosen recipe and check that the results land inside the prediction intervals (`42_doe_confirmation.png`).

Result: uniformity went from 1.27 % to 0.20 % one-sigma, and rate from 292 to 336 nm/min, with refractive index held on target.

### 2. "What does resolution IV mean?"

Main effects are not confused with two-factor interactions, but two-factor interactions are confused with each other. That's fine for screening, because the goal is to find which knobs matter, not to untangle interactions. You can prove it from the design matrix: every main-effect column is orthogonal to every product of two other columns. `test_fractional_factorial_is_orthogonal_resolution_iv` checks exactly that.

### 3. "Why centre points?"

- **Pure error.** Repeated runs at the same settings measure run-to-run noise without assuming the model is right.
- **Curvature check.** If the centre-point average differs from the average of the corners, the response is curved and a linear model isn't enough. Here, wet etch rate showed curvature (p = 0.0001), which is one reason to move to a response surface.

### 4. "Why did you model the profile instead of uniformity?" (a strong answer to bring up yourself)

One-sigma non-uniformity is a standard deviation, so it can never go below zero. Near its minimum it folds into a V shape, and a quadratic model can't fit a V. The model over-predicted NU at the optimum (0.29 % predicted vs 0.20 % measured).

So I modelled the *signed* edge-minus-centre profile instead. It passes smoothly through zero, and the model fits it better (R² 0.996 vs 0.968). This is a general rule: model the physical quantity, not a summary statistic of it.

### 5. "One of your confirmation runs failed. Is the model wrong?"

Refractive index missed its prediction interval by 0.00017. That's statistically real, but it's 2 % of the ±0.008 window that matters for the film. Two lessons:

- **Statistical vs practical significance.** A very precise measurement makes tiny model errors detectable.
- **Multiple comparisons.** Five responses tested at 95 % each gives about a 23 % chance of at least one failure by chance alone.

I'd note it and move on. If a miss were a large share of the engineering window, I'd run more experiments in that region.

### 6. "Your optimum is at the edge of the design space. What now?"

Pressure ended at its upper limit (3.0 Torr). The true optimum may lie outside the tested region, so the next step is a follow-up experiment that extends pressure. I'd also watch for constraints the model doesn't include: particles from gas-phase reactions at high pressure, and film stress.

## SPC and OCAP

### 7. "I-MR, EWMA or CUSUM: which would you use?"

- **I-MR:** catches large, sudden shifts quickly.
- **EWMA and CUSUM:** built for small, sustained drifts.

On the MFC drift, the I chart with all four Western Electric rules signalled first (wafer 116), but it also has the most false alarms: 6 episodes per 600 healthy wafers, versus 0 for EWMA. A common practical setup is the I chart with rule 1 only, plus EWMA (`43_spc_charts.png`).

### 8. "What is an OCAP and why have one?"

An Out-of-Control Action Plan: a pre-agreed list of what to check, in what order, and who decides what. It makes the 3 a.m. response consistent and fast. See `docs/OCAP_PECVD_THICKNESS.md` for the flowchart.

Mine is executable. On the MFC excursion it made three different decisions at three points:

| Wafer | Situation | Decision |
|---|---|---|
| 116 | Weak chart signal, nothing confirms it | Continue under watch |
| 260 | In spec, MFC drift diagnosed | Continue with a maintenance ticket |
| 450 | Out of spec | Hold the lot, take the chamber down |

A good OCAP prevents both missed excursions and unnecessary tool downs.

### 9. "Why context-correct the chart?"

Rate drifts down between chamber cleans as the walls coat, and resets after each clean. On a raw chart, that sawtooth either triggers alarms every cycle or forces the limits so wide that real faults hide. Modelling the expected value against wafers-since-clean and charting the residual removes the known pattern and leaves the unknown one.

## Run-to-run control

### 10. "How does EWMA run-to-run control work?"

The controller assumes thickness = rate × time. After each measured wafer, it updates its rate estimate as a weighted blend: λ × the new observation + (1 − λ) × the old estimate. Then it sets the next deposition time to target ÷ rate.

On the MFC drift, out-of-spec wafers dropped from 43 % to 0 %, and Ppk rose from 0.06 to 2.74 (`44_r2r_control.png`, `46_capability.png`).

### 11. "How do you pick λ?"

- **Too low** and the controller lags the drift.
- **Too high** and it chases measurement noise.
- **Longer metrology delays** raise the error at every λ.

Here the best value is about 0.6–0.7 (`45_r2r_lambda_tuning.png`).

### 12. "What can go wrong with run-to-run control?" (the most important answer in this module)

Three things I found while building it.

1. **Resetting after a clean can hurt.** It throws away what the controller learned, and a clean doesn't fix a drifting MFC.
2. **The first wafer after a clean poisons feedback.** It always runs low, so the controller over-corrects the next wafers. The fix is feedforward, which adds known extra time to post-clean wafers, combined with converting that wafer's measurement to its settled equivalent before learning from it. This cut the worst-case error from 9.6 nm to 5.6 nm.
3. **The controller hides faults.** It held thickness in spec while the MFC drifted, so a thickness chart goes blind. The drift moved into the controller's output, so I put an engineering limit on deposition time: an alarm when it's more than 1 % from expected. It fired at wafer 192, with zero false alarms on five healthy runs.

Statistical 3-sigma limits fail on controller output, because the output is a smoothed average and highly autocorrelated. That's why fabs use engineering limits there.

I also changed the Module 3 diagnoser to use deposition *rate* instead of thickness. The controller changes time, not rate, so diagnosis still works under control: 95 % correct on the MFC drift.
