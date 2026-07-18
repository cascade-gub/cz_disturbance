# TICKET — shared wavelet kernel + red-noise significance mask

Two coupled goals, in the order they should be done:

1. **Refactor** the duplicated wavelet code into **one shared kernel** so a fix lands in one place.
2. **Add a red-noise significance mask** — cheaply and deterministically — as the first feature that
   kernel earns you (right now it would have to be written three times).

**Files in scope — all three wavelet notebooks:**

| Notebook | Leg | Key params | Has today |
|---|---|---|---|
| `syn_data/analysis_02_wavelet_synchrony.Rmd` | synthetic (idealised) | `P≈500`, `dt=1`, single-band, clean sine driver | no gating, no sig mask |
| `ms_data/analysis_04_wavelet_synchrony.Rmd` | P→Q (real) | `P=365`, `dt=1` day | **phase gating** (from `analysis_04_lag_fix.txt`), no sig mask |
| `ms_data/analysis_05_cq_wavelet_synchrony.Rmd` | C→Q (real) | `P_DAYS=365`, native **`step`** cadence | no gating, no sig mask, `step` scaling |

**Do not touch** the tracked interface CSVs (`*__series_daily.csv`, `*__cq_native.csv`),
`prep_ms_series.Rmd`, or the peak-method notebook (`analysis_03`).

**Status:** deferred / not started. Handoff explainer only.

---

## PART 1 — The shared wavelet kernel (do this first)

### The problem

All three notebooks define the same helpers, copied by hand. They have already **drifted**:
`analysis_04` gained phase gating (`PHASE_MIN_AMP`, `phase_ok`, lag→NA where the driver has no band
amplitude) via the lag-fix ticket; `analysis_02` and `analysis_05` never got it. Any new feature
(significance mask, an event band, a caption fix) currently means editing three files and keeping
them in sync by eyeball. That is the root cause this ticket removes.

### The design

Create **one tracked R file** — `R/wavelet_kernel.R` at repo root (new `R/` dir) — that defines the
whole helper family **once**, parameterized by the handful of things that actually differ between
notebooks. Each notebook `source()`s it in its setup chunk and supplies a small config.

**What differs across the three** (everything else is identical):

- **`P`** — the band period *in the wavelet's own time unit*: 500 (syn), 365 d (04), 365/`step`
  (05, so the band lands at 365 *days*).
- **`step`** — native-cadence multiplier; `1` for 02/04 (daily), `step_days` for 05. Generalize
  `band_rows` to `abs(wc$Period * step - P_days) <= 0.05 * P_days` with `step = 1` default; then
  02/04 are just the `step = 1` case of 05's formula.
- **`EF`** — COI buffer; `sqrt(2) * P_days / step` steps. Derives from the two above.
- **phase gating** — `PHASE_MIN_AMP` (04's `f = 0.5`), on/off flag. See reconciliation note below.
- **channel semantics** — all three use columns `climate` / `log_response`; only the *meaning*
  differs (05: `climate = log(Q)`, `log_response = log(C)`). **The code is identical**; no
  parameter needed, just a comment.

**Recommended shape — a factory (closures over config).** Cleanest R idiom; avoids the fragile
implicit-global coupling the notebooks use now:

```r
# R/wavelet_kernel.R
make_wavelet_kernel <- function(cfg) {
  # cfg = list(P_days=, step=1, phase_gate=TRUE, PHASE_MIN_AMP=0.5, lowerPeriod=4, upperPeriod=730)
  step <- cfg$step %||% 1
  P    <- cfg$P_days / step                 # band period in the wavelet's time unit
  EF   <- round(sqrt(2) * cfg$P_days / step)

  band_rows <- function(wc) which(abs(wc$Period * step - cfg$P_days) <= 0.05 * cfg$P_days)
  wco       <- function(df) analyze.coherency(as.data.frame(df[, c("climate","log_response")]),
                              my.pair = c("climate","log_response"), loess.span = 0, dt = 1,
                              dj = 1/20, lowerPeriod = cfg$lowerPeriod, upperPeriod = cfg$upperPeriod,
                              make.pval = FALSE, verbose = FALSE)
  band_trace <- function(wc, df) { ... }    # gating applied iff cfg$phase_gate
  block_stat <- function(...) { ... }
  wave_summary <- function(tr) { ... }
  plot_wc <- function(wc, main) { ... }
  xpower_sig_level <- function(df, wc, conf = 0.95) { ... }   # PART 2, lives here
  list(P=P, EF=EF, band_rows=band_rows, wco=wco, band_trace=band_trace,
       block_stat=block_stat, wave_summary=wave_summary, plot_wc=plot_wc,
       xpower_sig_level=xpower_sig_level)
}
```

Each notebook setup chunk then:

```r
source("../R/wavelet_kernel.R")            # 02 is in syn_data/, 04+05 in ms_data/ → both one dir below root
K <- make_wavelet_kernel(list(P_days = 365, step = 1, phase_gate = TRUE))
# then call K$wco(df), K$band_trace(wc, df), K$plot_wc(wc, main), K$P, K$EF, ...
```

A lighter first step, if the factory feels like too big a diff: move the function *bodies* verbatim
into the sourced file but keep them referencing the notebook's existing globals (`P`, `EF`,
`PHASE_MIN_AMP`, `step`). Smaller diff, but keeps the implicit-global coupling — treat as a
stepping stone, not the destination.

### Build / sourcing mechanics (verify these)

- Rmd knits with **working dir = the .Rmd's own folder** (both the RStudio Knit button and
  `rmarkdown::render` in `build_docs.R`). All three notebooks live exactly **one dir below repo
  root**, so `source("../R/wavelet_kernel.R")` resolves from `syn_data/` *and* `ms_data/`.
  Confirm `build_docs.R` doesn't override `knit_root_dir`; if it does, use an absolute/anchored path.
- `R/wavelet_kernel.R` is **code, tracked, published-independent** — it is *not* a data interface
  like the `*_daily.csv` files, so it doesn't change the git-ignore story. Add it to the repo and
  mention it in the toolchain notes.
- Avoid adding a new package dependency for path resolution (`here`/`rprojroot` are **not** in the
  current stack: `rmarkdown, knitr, dplyr, tidyr, readr, ggplot2, patchwork, WaveletComp`). The
  `../R/` relative path needs none.
- Note it in `CLAUDE.md` (repo root) and `ms_data/CLAUDE.md` / `syn_data/CLAUDE.md`: "wavelet
  helpers are shared in `R/wavelet_kernel.R`; edit once."

### ⚠️ Reconciliation — the refactor is NOT automatically output-preserving

Because the notebooks drifted, unifying them can *change outputs*. Handle deliberately:

- **Phase gating** is in `04` only. Extracting a shared `band_trace` means `02` and `05` would
  start running it too:
  - **`02` (synthetic):** the driver is a clean strong sine, so `ax >= 0.5·median(ax)` holds
    almost everywhere → gating is effectively a **no-op**; outputs should be unchanged. Verify.
  - **`05` (C→Q):** the driver (`log Q`) is noisier; gating **may trigger and change numbers**.
    That is arguably a *correct* improvement, but it is a behavior change — decide per-notebook via
    `cfg$phase_gate` and, if enabled on 05, re-baseline its published figures/tables intentionally,
    not silently.
- **Diagnostics** added by the lag-fix (`xpower_spectrum` figure, `annual_r2` table) live in `04`
  only. The kernel can offer them to all three; adding them to 02/05 is a deliberate addition.
- **`P` differs** (500 vs 365) and `02` has no `disturbance`-span buffering subtlety that 04/05 do
  (real phases). Keep those in config, don't hardcode.

**Protocol:** do the refactor as an **output-preserving extraction first** — same `cfg` per
notebook reproducing today's behavior (gating OFF for 02/05 initially), re-knit all three, confirm
**every number and figure is unchanged**. Only *then* flip on new behavior (gating for 05, the
significance mask) as separate, reviewed changes.

---

## PART 2 — Red-noise significance mask (the first feature the kernel earns)

### Why it's needed

The coherence heatmaps (`plot_wc`, `which.image = "wc"`) sit near **1 across the whole record** at
every site — the mathematical resting state of normalized wavelet coherence, not evidence of
coupling:

- Un-smoothed, coherence is **identically 1 everywhere** (`|W_xy|² = |W_x|²·|W_y|²` pointwise,
  Cauchy–Schwarz). Its entire dynamic range is manufactured by WaveletComp's time/scale smoothing.
- On **red (autocorrelated)** series — all of these, log-discharge/log-nitrate especially — phase
  barely moves within the smoothing window, so even weakly-related signals stay coherent.

So a bright coherence field is an *expected artifact*. The honest signal is in **un-normalized
cross-wavelet power** (`wc$Power.xy`) and band amplitude — already the notebooks' thesis. What's
missing is a significance reference: with no mask, every bright patch looks equally earned.
`make.pval = TRUE` would add one but was cut as **too expensive** (`n.sim` surrogates × full CWT
re-run × sites) and non-deterministic (RNG).

### The key realization: mask the right quantity

Coherence significance **genuinely needs Monte Carlo** — its distribution depends on the smoothing
and has no closed form (Torrence & Webster 1999; Grinsted et al. 2004). That is *why* it's
expensive. **Cross-wavelet power, by contrast, has an analytic red-noise significance level** —
no simulation (Torrence & Compo 1998) — and cross-power is the channel the notebooks *promote*. So
mask cross-power analytically instead of paying MC to mask the quantity you're demoting.

### OPTION A — analytic AR(1) cross-power significance *(recommended; deterministic, ~free)*

AR(1) red-noise theoretical spectrum, unit-variance normalized (T&C 1998, eq. 16):

```
P(f) = (1 − α²) / (1 + α² − 2α·cos(2πf))          f = normalized freq = 1 / wc$Period   (dt = 1)
```

95% red-noise level for cross-wavelet power `|W_xy|/(σ_x σ_y)`:

```
level(f) = (Z_ν(p)/ν) · √( P_x(f)·P_y(f) )          ν = 2,  Z_2(0.95) ≈ 3.999
```

Two simplifications here: `loess.span = 0` standardizes each series to **unit variance**
(`σ_x=σ_y=1`, and `wc$Power.xy` is already standardized); and `α` is **standardization-invariant**,
so read it straight off the raw channels. Kernel helper:

```r
ar1_spectrum <- function(alpha, f) (1 - alpha^2) / (1 + alpha^2 - 2 * alpha * cos(2 * pi * f))

xpower_sig_level <- function(df, wc, conf = 0.95) {          # one threshold per wc$Period row
  ax <- acf(df$climate,      1, plot = FALSE)$acf[2]
  ay <- acf(df$log_response, 1, plot = FALSE)$acf[2]
  f  <- 1 / wc$Period                                        # normalized freq (dt = 1, i.e. per step)
  Z2 <- if (isTRUE(all.equal(conf, 0.95))) 3.999 else stop("tabulate Z_2 for other conf")
  (Z2 / 2) * sqrt(ar1_spectrum(ax, f) * ar1_spectrum(ay, f))
}
```

**Where it plugs in (once, in the kernel; every notebook gets it):**
- `analysis_04`'s `xpower-spectrum` figure (mean `Power.xy` vs period): add the threshold as a
  reference line/ribbon and mark periods where observed cross-power beats it. This gives the
  "coupling lives at the event band, not the annual band" argument a red-noise reference.
- `analysis_02`/`05`: add the same cross-power-vs-period diagnostic (02 doesn't have it yet), and/or
  a per-band significance flag in `band_trace`/`wave_summary` alongside the existing gating columns.
- **`step` note (05):** work in **step units** — `f = 1/wc$Period` with `wc$Period` in steps, and
  `α` from the binned channels. Do **not** multiply Period by `step_days` here; `α` is the lag-1
  autocorrelation at the native cadence, so frequency must match. (This is automatic if the kernel
  keeps `step` in config and only converts to days for display.)

### ⚠️ IMPLEMENTATION CHECK (before trusting the contour)

1. **Is `wc$Power.xy` a modulus or a modulus²?** Decides whether you compare it to `level(f)` or
   `level(f)²`. Confirm from the installed WaveletComp's docs/source — don't assume.
2. **Calibrate once against a seeded MC run.** `set.seed(1); analyze.coherency(..., make.pval=TRUE,
   n.sim=100)` on the **smallest** site only, compare its significance contour to analytic
   `level(f)` at the annual and an event band. Agreement ⇒ the normalization/constant is right and
   you can drop MC entirely. This throwaway run is validation, not part of any knit.

### OPTION B — cached MC coherence mask *(optional; for the coherence heatmap exhibit only)*

If you want the classic stipple **on the coherence image**, don't pay at knit time — mirror the
repo's tracked-CSV interface pattern:

1. **Offline script** (`ms_data/build_wc_sigmask.R`) runs `wco(make.pval=TRUE)`, **seeded**, per
   site, writes a small **tracked** artifact (`<site>__wc_sigmask.csv`: p-field or contour coords).
2. Kernel's `plot_wc` reads it and overlays `contour(...)` on the live `wc.image` coordinates (same
   `graphics.reset = FALSE` trick the current annotations use).

Cost leaves the render loop; deterministic at knit (frozen mask). Document the offline rebuild in
`ms_data/CLAUDE.md`, same as `prep_ms_series`.

**Cost knobs if MC is ever run live** (stack them): **seed it** (`set.seed(k)` before the call —
base R RNG, restores determinism); **`n.sim = 50`** (30 for α=0.1); **restrict to the band scale
rows** you actually read; **right null** (`method = "AR"`/`"Fourier.rand"`, not `"white.noise"`
which is wrong for red data); **parallelize** across surrogates/sites (`parallel::parLapply` /
`future` `multisession` — Windows, no fork).

---

## Constraints (carry over from `analysis_04_lag_fix.txt`)

- **Deterministic default path:** no RNG, no `make.pval = TRUE` in what any knit runs. Option A is
  closed-form. Option B's MC lives in a separate offline script → tracked CSV.
- **Reproducible from the tracked series CSVs alone** for the default knit path.
- **Apply uniformly** across sites and, now, across all three wavelet notebooks via the kernel.
- **Minimal diffs; output-preserving refactor first.** Preserve the dark-page/white-matted theme
  and the knit → `../docs/` copy hook.

## Verify before done

- **Refactor:** all three notebooks source `R/wavelet_kernel.R`; with feature flags set to today's
  behavior, re-knit **02, 04, 05** and confirm every number/figure is **unchanged** (the
  extraction is output-preserving). Any intended change (05 gating) is re-baselined deliberately.
- **Significance:** cross-power figures show the red-noise threshold; short-period shared power
  flags significant where it beats red noise, annual band (correctly) does not; analytic contour
  was calibrated against one seeded MC run and they agree.
- If Option B done: cached sigmask overlays on the coherence heatmaps; offline rebuild documented.
- HTML copied to `../docs/` by the knit hook for all three.
- One-line pointers added to `CLAUDE.md` + both leg `CLAUDE.md`s that wavelet helpers now live in
  `R/wavelet_kernel.R` (edit once).

## References

- Torrence & Compo (1998), *A Practical Guide to Wavelet Analysis*, BAMS — AR(1) red-noise spectrum
  (eq. 16), cross-wavelet power significance (`Z_ν`).
- Torrence & Webster (1999) — coherence significance requires Monte Carlo (smoothing-dependent).
- Grinsted, Moore & Jevrejeva (2004), NPG — practical MC significance for wavelet coherence.
