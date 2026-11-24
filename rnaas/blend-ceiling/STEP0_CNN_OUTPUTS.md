# Step 0: what the Salinas network outputs

Checked 2026-11-24. Result: **the single-transit CNN gives no transit or detection score.**
It is a regression network only. So, per the brief, Steps 1 to 3 run and the work stops
for Jon to pick a score (options at the end).

## Code inspected

- `helemysm/singletransiters` at commit `b7b9d49eded61643ceb715f9f99a0685fcbba7aa` ("print
  metrics", same as `origin/main` on 2026-11-24). Reused the existing clone at
  `salinas-second-transits/singletransiters` in the main checkout, with its own
  Python 3.10 env `.venv-cnn`. The clone has no local changes, nothing from it is copied
  here, and it has no license file.
- `helemysm/FII_TransformerNN` at commit `3f0ebcb6f860b33b76ba15a2bb2c05dbd59ab640`
  (2025-04-05), read through the GitHub API only and not cloned. It also has no license file.

## What each module outputs

| Module | Function / class | Output |
|---|---|---|
| `model/model.py` | `SingleTransitCNN.forward` | `(mu, log_sigma)`, each of shape (B, 3): depth, duration and `dtcorr` (time offset from `t0_hat`). There are no other heads. |
| `model/model.py` | `SoftArgmaxDtHead` | A softmax over time bins inside the window. It is used only to place the transit centre (`mu`, `sigma` in hours), and the weights are not returned. |
| `computeparams/param_inference_core.py` | `_run_inference_on_lc(lc, t0_hat)` | A dict with `depth_ppm`, `duration_h`, `dtcorr_h`, their calibrated sigmas, the predicted mid, start and end times, and the 99.5 h window (`t_rel_days`, `flux_window`). Our `step_cnn` wrapper calls this. |
| `computeparams/predict_csv.py`, `compute_param_single.py` | CLI wrappers | The same quantities, plus a radius from the TIC or Gaia R*. |
| `baseline/baseline_trap.py` | `fit_transit` | A trapezoid fit, used as the paper's comparison baseline. It is also regression. |
| `eval/*` | metrics | MAE/RMSE of the regression targets on synthetic sets. |

Checkpoint `models/baseline_conv1d.pt`:
- the state dict holds only `feat`, `mu_head_depth`, `sigma_head_depth`, `mu_head_dur`,
  `sigma_head_dur`, `dt_head` and `sigma_head_dt` (plus the `dt_mean` and `dt_std` buffers);
- the non-weight keys are `model_type=cnn`, `use_time_channel`, `center_flux` and the
  log-sigma clamps.

There is no classifier head anywhere in the code or the weights.

Training data (`data_gen/simgen_core.py`, `generate_data.py`):
- synthetic only, from batman transit models plus white noise, AR(1) red noise, a
  quasi-periodic stellar GP, a polynomial trend and outliers;
- **every training window contains a transit**: samples with integrated SNR < 3 are rejected,
  and b ≤ 0.90, so grazing geometries are excluded;
- no windows without a transit and no eclipsing binaries.

So the network has never seen a negative. Any output used as a score on FP/FA windows would be
off its training distribution.

Preprocessing, which the brief says to keep for any later scoring
(`param_inference_core.py`):
1. quality == 0, then remove NaNs;
2. bin to 30 min;
3. a fixed 99.5 h window (200 points) centred on `t0_hat`, linearly interpolated, with
   out-of-range points filled with the median;
4. divide by the window median and subtract 1;
5. the input is 2 channels: scaled time in [-1, 1] and flux.

PDCSAP is used if present, else SAP (`pick_lightcurve`).

## How the paper used it (arXiv 2610.24882 v1, HTML, read 2026-11-24)

- **The CNN is a characterisation step, not a detector.** The paper says detection "goes
  beyond the scope of this work" (Sec. 7). `t0_hat` is "provided by a prior detection
  pipeline" (Sec. 3.1).
- **Real-data validation:**
  - Sec. 5.1 uses 120 TOIs with TFOPWG disposition CP or KP and P > 10 d;
  - SPOC light curves, with QLP where SPOC is missing, and a random ±25 h jitter on the epoch;
  - metrics are depth, duration and Δt median/MAE (Table 3) and Δt failure rates (> 5 h:
    CNN 2% vs baseline 13%; > 10 h: 1% vs 10%);
  - no FP/FA set is used;
  - the 120 TICs are not listed in the paper, so a train/validation overlap check against our
    sample cannot be done from the paper. It would need the list from the authors.
- **Candidate search (Sec. 6.1):** the new candidates in Tables 4 to 6 came from a separate
  **Transformer classifier** (Salinas et al. 2025), retrained on ExoFOP labels with
  "negative labels to EBs and non-transit signals". Existing TOIs and CTOIs were excluded.
  BLS supplied periods for multi-transit systems, and the CNN then characterised each event.
- **Future work (Sec. 8):** the paper names false-positive vetting and joint detection plus
  characterisation as possible extensions only.

## The Transformer (the paper's actual detector)

`FII_TransformerNN`:
- `src/train_.py` trains `NNClassifier` with `nn.BCELoss()` and reports ROC AUC;
- `src/core/encoder_model.py` has a `ClassificationModule`;
- **so it does output a transit probability, but the repo has code only and no trained
  weights**;
- it predates the paper's retraining (last push 2025-04-05);
- it has no license.

Using it would mean one of two things:
- training it ourselves on labels we choose, which is close to option (b) below, with their
  architecture;
- asking the authors for the retrained weights. That is outreach and Jon's call. Those weights
  were trained on ExoFOP labels, which very likely include many of the TFOP FP/FA TOIs in our
  test set, so the test would leak unless the authors also supply their training list.

## Options for Jon (the brief's two, unchanged; nothing built yet)

- **(a) Proxy from the CNN's own outputs.** For example, the predicted depth over its
  calibrated sigma, `depth_ppm / sigma_depth`. It can be computed today with no training. The
  note would have to state plainly that it is a proxy: the network was trained only on windows
  that contain a planet transit, so the test becomes "does the regression's confidence differ
  between planets and FPs", not a classifier test.
- **(b) Our own baseline classifier** on the same 99.5 h, 30 min windows (scikit-learn on
  window features, or a small 1D CNN in torch), trained on labels kept apart from the test
  TOIs (for example, train on short-period TOIs and test on the long-period sample). This is a
  real classifier, but it is ours, not the Salinas network, so the note's framing changes to
  "a light-curve-only classifier".
