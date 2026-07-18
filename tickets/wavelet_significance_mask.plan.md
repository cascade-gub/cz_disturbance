# PLAN — wavelet red-noise significance mask (+ shared kernel)

Implementation plan for `wavelet_significance_mask.md`. Read that ticket first for the *why*;
this file is the *how*, grounded in the current state of the three notebooks.

## STATUS: IMPLEMENTED (both phases, output-preserving, verified)

- **Phase 1 (mask):** `R/wavelet_kernel.R` gained `ar1_spectrum`, `xpower_sig_level` (pointwise
  field), `xpower_sig_ref_tavg` (DOF-inflated time-average), `add_xpower_contour`. Wired into all
  three notebooks: a persistence-filtered significance **contour** on the coherence heatmaps **and**
  a red-noise 95% **reference** on the cross-power-vs-period figure (added that figure to 02/05).
  Additive: all tables byte-identical to baseline.
- **Phase 2 (kernel):** `make_wavelet_kernel(cfg)` now builds the WHOLE helper family; each notebook
  is `list2env(make_wavelet_kernel(cfg), environment())` + config. Output-preserving: re-knit 02/04/
  05, all tables byte-identical and numeric fingerprints unchanged except intended superset columns.
- **Deltas from this plan, all confirmed empirically:**
  - the analytic level needed a **`/wc$Scale`** factor the ticket's formula omitted (WaveletComp
    rectifies `Wave.xy` by `/Scale`) — calibrated against seeded AR(1) surrogates (`scratch calib*.R`);
  - the pointwise ν=2 field **speckles** on long daily records, so the contour keeps only coupling
    **sustained ≥3 cycles** (a duration/area reduction) — user-chosen;
  - both exhibits shipped (user chose "Both"); time-averaged reference validated (π/4 mean-const,
    T&C DOF), tight at event periods, mildly conservative at the annual band (documented in-kernel);
  - 02 prose corrected: Desync's collapse is time-resolved, not visible in a whole-record mean.
- **Verification layers used:** table byte-diff + numeric fingerprint of plot-input objects
  (`scratch fingerprint.R`) + per-plot-type visual spot-checks (the fingerprint/table diffs are
  blind to plot-only styling — season shading, axis units, lag label, wrap-break, subtitle/ylab
  codepoints — so those were checked by eye and by byte-comparing the exact strings).

## Headline: invert the ticket's order (deliberate deviation)

The ticket says **refactor the kernel first, add the mask second**. After reading all three
notebooks, I recommend the **opposite**, and the reasoning is the drift the ticket only partly
accounts for.

**What the ticket assumes drifted:** phase gating (04 only), plus `P`, `step`, `EF`.

**What actually drifted** (verified by reading `analysis_02/04/05`): the helpers differ
*structurally*, not just by the listed params —

| Helper | 02 (syn) | 04 (P→Q) | 05 (C→Q) |
|---|---|---|---|
| `wco` period bounds | `P/2^1.5 … P*2^1.5` | `4 … 730` | `4 … floor(nrow/3)` (data-dependent) |
| `block_stat` signature | `(t, v, circ)` | `(t, v, circ)` + NA-drop + returns `R` | `(t, v, blk, circ)` |
| `wave_summary` phases | uses `transition` span | `pre/disturbance/post` | `pre/disturbance/post` + `step` |
| `band_trace` | no gate, no `date`, no `driver_amp` | gate + `date` + `driver_amp` + `phase_ok` | `step` + `date`, no gate |
| `plot_wc` | plain image | `wc_marks` + annotated overlay | `spec.period.axis` step-axis |
| `plot_band_trace` | plain facets | `seg_id` wrap-breaking + antiphase rails | (date axis) |

A kernel configurable enough to reproduce **all three exactly** is a real amount of config, and
the divergent helpers are where the risk lives — so we sequence it, not skip it.

**The mask, by contrast, is new, pure, and identical across all three.** `ar1_spectrum` +
`xpower_sig_level` + one plot layer *add* a diagnostic; they touch no existing computation, so
they carry ~zero output-preservation risk. Doing the mask first delivers the **named feature** and
gets a real, low-risk `R/wavelet_kernel.R` on disk that Phase 2 then grows into the full kernel. So
the flip is about **sequencing**, not scope:

- **Phase 1 (do first): minimal kernel = the mask only.** Independently shippable; ~zero risk.
- **Phase 2 (committed, not optional): the full helper unification.** Every wavelet helper moves
  into the kernel, parameterized by config; each notebook becomes thin config + prose. Done
  helper-by-helper with a re-knit + diff after each.

**Why Phase 2 is required, not optional:** the drift *is* the problem this ticket exists to kill.
Leaving the structurally-divergent helpers notebook-local — even the awkward ones — is precisely
the mechanism that let 02/04/05 diverge in the first place. "Edit once" only holds if *all* the
helpers live in one place; a half-migrated kernel just relocates the drift. The divergent helpers
are unifiable via config (schema in Phase 2 below) — it costs more flags, and that cost is the
price of not drifting again.

---

## ⚠️ Correctness trap the ticket gets wrong — fix before coding

The ticket's Option A says: compute `xpower_sig_level` (a **point-wise, ν=2, single-estimate 95%**
level) and overlay it on the **`xpower_spectrum` figure**. But that figure plots
`rowMeans(wc$Power.xy)` — a **time-averaged** spectrum. Overlaying a single-estimate threshold on a
time-mean is a **degrees-of-freedom mismatch**: time-averaging shrinks the observable's variance
while the ν=2 threshold stays at single-estimate height, so the mean sits far below the line and
the test is drastically too conservative — genuine structure won't clear it. (Torrence & Compo
1998 handle time-averaged spectra with **inflated DOF, ν ≫ 2**, and a correspondingly *lower*
level.)

**Resolve by picking one — decide in Step 0:**

- **(A-field) Apply the ν=2 level to the *instantaneous* cross-power field** — i.e. contour
  `wc$Power.xy` on the coherence/`plot_wc` heatmap where ν=2 is exactly right. Simplest correct
  option; no DOF math. **Recommended.**
- **(A-avg) Derive a time-averaged threshold with proper DOF** for the `xpower_spectrum` figure —
  compute ν from the number of independent estimates in the time-mean (T&C 1998 §5) and lower the
  level accordingly. More math, more calibration surface.

The plan below is written for **A-field** and notes where A-avg diverges.

---

## Phase 1 — the minimal kernel + mask (recommended, ship this)

### Step 0 — decide the display quantity (blocks everything)
Choose **A-field** (instantaneous contour) vs **A-avg** (time-averaged threshold). Recommend
A-field. This choice fixes what the calibration in Step 4 must validate — they must be the **same
quantity**.

### Step 1 — settle two WaveletComp facts (cheap, do once, no knit)
1. **Is `wc$Power.xy` a modulus or a modulus²?** Decides whether observed compares to `level(f)`
   or `level(f)²`. Confirm from the installed WaveletComp source/docs — do **not** assume.
2. **What normalization is `Power.xy` in?** The `(Z₂/2)·√(P_x·P_y)` constant assumes unit-variance
   standardized channels (true here: `loess.span = 0` ⇒ `/sd`). Confirm `Power.xy` is the
   standardized cross-power, not re-multiplied by variances.

Capture both as one-line comments in the kernel next to `xpower_sig_level`.

### Step 2 — create `R/wavelet_kernel.R` (new `R/` dir at repo root)
Minimal contents — **only the new, pure, shared pieces**:

```r
# R/wavelet_kernel.R — shared wavelet significance helpers (see tickets/wavelet_significance_mask.md).
# Phase 1 scope: red-noise cross-power significance only. The per-notebook wco/band_trace/plot_*
# helpers are intentionally NOT unified here yet (they have drifted structurally — see the plan).

# AR(1) red-noise theoretical spectrum, unit-variance normalized (Torrence & Compo 1998, eq. 16).
# f = normalized frequency = 1 / Period, in the wavelet's own time unit (per step; dt = 1).
ar1_spectrum <- function(alpha, f) (1 - alpha^2) / (1 + alpha^2 - 2 * alpha * cos(2 * pi * f))

# 95% red-noise level for standardized cross-wavelet power |W_xy|/(sd_x sd_y), point-wise (nu = 2).
# alpha read straight off the RAW channels: lag-1 autocorrelation is standardization-invariant.
# NOTE (Step 1): confirm whether wc$Power.xy is |.| or |.|^2 → compare to level or level^2.
# NOTE (DOF): this nu=2 level is for the INSTANTANEOUS field. Do NOT overlay on a time-MEAN
#   spectrum without inflating nu (see plan, "correctness trap").
xpower_sig_level <- function(df, wc, conf = 0.95) {
  ax <- acf(df$climate,      1, plot = FALSE)$acf[2]
  ay <- acf(df$log_response, 1, plot = FALSE)$acf[2]
  f  <- 1 / wc$Period                              # per-step frequency (05: Period already in steps → correct)
  Z2 <- if (isTRUE(all.equal(conf, 0.95))) 3.999 else stop("tabulate Z_2 for other conf")
  (Z2 / 2) * sqrt(ar1_spectrum(ax, f) * ar1_spectrum(ay, f))   # one threshold per wc$Period row
}
```

Notes carried from the ticket that this respects:
- **05/`step`:** work in **step units** — `f = 1/wc$Period` with `Period` in steps, `α` from the
  binned channels. Do **not** multiply `Period` by `step_days` here. Correct as written (no `step`
  arg needed; the function reads the channels and periods of whatever `wc`/`df` it's handed).
- **Tracked, code not data:** `R/wavelet_kernel.R` is tracked R, not a data interface — no change
  to the git-ignore story.

### Step 3 — wire the mask into each notebook (3 small, additive diffs)
Each notebook's setup chunk gains one line:
```r
source("../R/wavelet_kernel.R")   # 02 in syn_data/, 04+05 in ms_data/ → both one dir below root
```
Path verified: all three live exactly one dir below root, and `build_docs.R` calls
`rmarkdown::render(rmd)` **without** overriding `knit_root_dir`, so knit wd = the `.Rmd`'s own
folder (same reason the existing `read_csv("…__series_daily.csv")` relative reads work). `../R/`
resolves from `syn_data/` and `ms_data/`. No new package (`here`/`rprojroot`) needed.

Then, **A-field** wiring — add a significance contour to `plot_wc`:
- In each `plot_wc`, after `wc.image(..., graphics.reset = FALSE)`, compute
  `lev <- xpower_sig_level(df, wc)` (pass `df` in — 02's `plot_wc` currently takes only `wc`/`main`,
  so add the `df` arg there), build the ratio field `wc$Power.xy / lev` (broadcast `lev` down each
  Period row), and `contour(..., levels = 1, add = TRUE)` on the live `wc.image` coordinates using
  the same `graphics.reset = FALSE` trick 04/05 already use for their guide lines. 04 and 05
  already keep the coordinate system live; 02 must add `graphics.reset = FALSE` to its `wc.image`
  call (small, additive).
- **05 axis:** the contour is drawn in native-step Period coordinates (`log2` y like the guides),
  so it lands correctly under the existing `spec.period.axis` day labels — no extra scaling.

If **A-avg** were chosen instead: add the threshold as a reference line/ribbon on the
`xpower_spectrum` figure (04 has it; add the same diagnostic to 02/05), using a DOF-inflated level.
Not recommended for the first pass.

### Step 4 — calibrate once against a seeded MC run (throwaway, not knitted)
On the **smallest** site only:
```r
set.seed(1); wc_mc <- analyze.coherency(<same df/args>, make.pval = TRUE, n.sim = 100)
```
Compare the MC significance contour to the analytic `level(f)` **on the same quantity chosen in
Step 0** — i.e. for A-field, compare MC's instantaneous-field significance to the analytic field
contour at the annual band and at one event band. Agreement ⇒ the modulus/constant/normalization
are right and MC can be dropped entirely. This run is **validation only**, never part of any knit
(keep it in a scratch script, not a chunk).

> The advisor's caveat: a seeded-MC contour on the instantaneous field validates the **field** use,
> not a time-averaged overlay. Match the calibration to the display quantity.

### Step 5 — re-knit all three, confirm additive-only
Because the mask only *adds* a contour/figure, every existing number and figure must be
**unchanged**. Re-knit 02, 04, 05; confirm the coherence heatmaps now carry the significance
contour and nothing else moved. Confirm the knit `../docs/` copy fired for all three.

### Step 6 — docs
One-line pointer in `CLAUDE.md` (root) + `ms_data/CLAUDE.md` + `syn_data/CLAUDE.md`: "shared
wavelet **significance** helpers live in `R/wavelet_kernel.R`." Add `R/wavelet_kernel.R` to the
toolchain note. (Phase 2 upgrades this to "**all** wavelet helpers are shared — edit once.")

### Phase 1 done when
- All three notebooks `source("../R/wavelet_kernel.R")` and show the red-noise significance
  contour on their coherence heatmaps; short-period shared power flags where it beats red noise,
  the annual band (correctly) does not at the degenerate sites.
- Analytic contour calibrated against one seeded MC run on the smallest site — they agree — and MC
  appears in **no** knit path (deterministic, no RNG, reproducible from tracked CSVs).
- Every pre-existing number/figure unchanged; HTML copied to `../docs/` for all three.
- Doc pointers added.

---

## Phase 2 — the full helper unification (committed)

**Goal:** every wavelet helper lives in `R/wavelet_kernel.R`, built once by a factory over a small
config; each notebook's `helpers` chunk shrinks to `K <- make_wavelet_kernel(cfg)` + the analysis
prose. Nothing wavelet-computational stays notebook-local — that is the only end-state that stops
the drift.

### The factory + config schema

```r
# R/wavelet_kernel.R (Phase 2 grows the Phase-1 file into this)
`%||%` <- function(a, b) if (is.null(a)) b else a

make_wavelet_kernel <- function(cfg) {
  # ---- config (every cross-notebook difference is ONE of these) ----------------------------
  step   <- cfg$step        %||% 1                 # native-cadence multiplier (05: step_days; else 1)
  Pd     <- cfg$P_days                             # band period in DAYS (02: 500, 04/05: 365)
  P      <- Pd / step                              # band period in the wavelet's own unit (steps)
  EF     <- round(sqrt(2) * Pd / step)             # COI / disturbance buffer, in steps
  blk    <- max(2, round(Pd / step))               # block size for block_stat, in steps
  lowerP <- cfg$lowerPeriod %||% 4                  # 02: P/2^1.5 ; 04/05: 4
  upperP <- cfg$upperPeriod                         # value OR function(df); 02: P*2^1.5, 04: 730, 05: \(df) floor(nrow(df)/3)
  gate   <- cfg$phase_gate  %||% FALSE              # 04: TRUE ; 02/05: FALSE (until re-baselined)
  A_MIN  <- cfg$PHASE_MIN_AMP %||% 0.5
  R_MIN  <- cfg$PHASE_MIN_R   %||% 0.5
  has_tr <- cfg$has_transition %||% FALSE           # 02: TRUE (post starts after `transition`); 04/05: FALSE
  marks  <- cfg$wc_marks                            # named period vector or NULL (02 currently none)
  xcol   <- cfg$trace_x %||% "time"                 # 02/04: "time" ; 05: "date"
  wrap   <- cfg$lag_wrap_break %||% FALSE           # 04: TRUE (break lag line at ±branch cut + antiphase rails)
  # ... build every helper below as a closure over these, then return them in a list.
}
```

Everything the three notebooks differ by is now exactly one config field. The **divergent
helpers are unified to their superset**, config selecting behavior:

| Helper | Unification strategy |
|---|---|
| `band_rows` | `abs(wc$Period * step - Pd) <= 0.05 * Pd` (05's formula; `step = 1` recovers 02/04). |
| `wco` | period bounds from `lowerP`/`upperP`; `upperP` accepts a number **or** `function(df)` so 05's `floor(nrow/3)` fits with no special case. |
| `block_stat` | adopt the **superset**: `(t, v, circ)` with internal `blk` (from config) + 04's NA-drop + returns `R`. 02/05 have no NAs (drop = no-op) and ignore `R` — behavior identical. |
| `band_trace` | always compute `driver_amp`, `phase_ok`, `date`; apply gating **iff `gate`** (`ok <- if (gate) ax >= A_MIN*median(ax[coi]) else TRUE`). `angle`/`lag` scale by `P`/`Pd` per unit. Extra columns are harmless where unplotted. |
| `wave_summary` | post-window start = `max(time[phase == (if (has_tr) "transition" else "disturbance")]) + EF`; report φ/lag gated by `R_MIN` **iff `gate`** (matches 04; off ⇒ 02/05 current behavior). |
| `plot_wc` | one renderer taking `df` (for the sig contour), `marks`, `step`. `marks = NULL` ⇒ 02's bare image; `step > 1` ⇒ 05's `spec.period.axis`. Significance contour (Phase 1) folds in here. |
| `plot_band_trace` | `xcol` picks time/date axis; `wrap = TRUE` adds 04's `seg_id` branch-cut breaking + antiphase rails, `FALSE` gives 02/05's plain facets. |
| easy wins | `smooth_ma`, `CMEAN`, `metrics_long`, `plot_wave_bars`, `emit_table`, `season_layer`, `phase_band_layers`, `fit_pre`, `resid_trace`, `xpower_spectrum`, `annual_r2` — near-identical; move verbatim. |

### Migration protocol (per helper, in this order)

Extract easy wins first, then the divergent ones, **one at a time**:

1. Move the helper into the factory as a closure; delete it from the notebook; call `K$helper(...)`.
2. Re-knit that notebook; **diff every number and figure against the pre-move baseline** (capture
   baselines by knitting all three *before* Phase 2 starts and stashing the HTML/PNGs).
3. Config is set to reproduce **today's** behavior — `gate = FALSE` for 02/05, `has_transition =
   TRUE` for 02, `marks = NULL` for 02, etc. The extraction must be **output-preserving**; any diff
   is a bug in the extraction, fixed before moving on.

### Deliberate behavior changes (separate, reviewed commits — only after extraction is clean)

These are *improvements the unified kernel now makes trivial*, but they change outputs, so they are
**not** silent:
- **Phase gating for 05** (`gate = TRUE`): the noisier `log Q` driver may trigger it → re-baseline
  05's figures/tables deliberately. Arguably correct, but a decision, not a side effect.
- **`xpower_spectrum` / `annual_r2` diagnostics** for 02/05 (04 has them): the kernel offers them
  to all three; adding them is a deliberate addition.
- **`plot_wc` guide marks / harmonized styling** for 02: if we decide 02 should gain the
  week/month/year guides, that is a visual re-baseline, flagged as such.

Keep these as follow-on commits so the extraction diff stays clean and reviewable.

### Phase 2 done when
- All three notebooks' `helpers` chunk is just `K <- make_wavelet_kernel(cfg)`; no wavelet helper
  is defined in any `.Rmd`.
- With config set to today's behavior, 02/04/05 re-knit **byte-for-behavior identical** to the
  pre-Phase-2 baselines (any intended change is a separate re-baselined commit).
- Doc pointers upgraded to "all wavelet helpers shared in `R/wavelet_kernel.R` — edit once" in root
  + both leg `CLAUDE.md`s.
- HTML copied to `../docs/` for all three.

---

## Out of scope

- **Option B** (cached MC coherence mask → tracked `<site>__wc_sigmask.csv` + offline
  `build_wc_sigmask.R`). Option A already meets the deterministic-knit constraint; Option B is a
  separate, optional exhibit for the coherence heatmap and can follow later if the stipple is
  wanted.
- Touching the tracked interface CSVs, `prep_ms_series.Rmd`, or `analysis_03` (per the ticket).

## Constraints honored throughout
Deterministic default path (no RNG, no `make.pval = TRUE` in any knit); reproducible from tracked
series CSVs alone; applied uniformly across sites and all three notebooks; minimal, additive diffs
in Phase 1; dark-page / white-matted theme and the knit → `../docs/` copy hook preserved.
