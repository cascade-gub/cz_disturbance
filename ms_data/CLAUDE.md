# CLAUDE.md — Real-data (MacroSheds) disturbance–response analyses

Project context for `ms_data/`. This is the **real-data leg** of a two-part **methods
comparison**: how two different approaches to quantifying driver→response *synchrony* — the
central-tendency / peak-spacing approach and the Morlet-wavelet approach — behave on real
watersheds, versus on the idealised synthetic system in `../syn_data/`. The synthetic leg
(analyses I/II) frames the *theory and expectations* — what each method can and cannot see. This
leg runs the **same methods** on real, less-deterministic synchrony to see how they hold up. Each
method is a different, partial perspective on the system; some more complete than others. **The
real-data messiness is the stress-test bench, not a defect** — how each method responds (what it
recovers, what it loses, where it is fooled) is the finding.

Two real couplings are analyzed, ordered by how *causal* they are:

- **P→Q** (analyses III/IV) — precipitation → discharge. Near-mechanical (rain drives flow); the
  reference case. III = peak/moments, IV = wavelet.
- **C→Q** (analysis 05) — discharge → **nitrate concentration**. Deliberately *less causal /
  noisier*: nitrate is only partly set by flow (sources, hysteresis, dilution-vs-flushing, and at
  two sites a chemical manipulation acting directly on N). The wavelet's hardest real test.

> Companion doc: `README.md` (human orientation to the folder). Full generative spec of the
> *idealised* system: `../syn_data/CLAUDE.md`. Repo-wide context + toolchain: `../CLAUDE.md`.

## The datasets it analyzes

Two families of tracked CSVs, both built by `prep_ms_series.Rmd` from the local MacroSheds pull
(`../data/*_daily.csv`) + the screened disturbance windows (`../qualifying_disturbances.csv`).
**These CSVs are the only interface between the pull and the analyses** — as `synthetic_*.csv` are
in `syn_data/`. Analyses 03/04/05 need nothing but these (no `macrosheds`/feather dependency).

- **P→Q:** five `<site>__series_daily.csv` (one per disturbance) — daily, precip + discharge.
- **C→Q:** four `<site>__cq_native.csv` (MC06 dropped) — each on that site's **native NO3 cadence**
  (`step_days` bins), discharge + nitrate.

### Schema (all five files, identical)

Each row is one **day**. Columns:

- `time` — integer day index `1 … N` (the unit all peak/wavelet math runs in; `dt = 1` day).
- `date` — calendar date.
- `climate` — **daily precipitation (mm)**. The exogenous, *undisturbed* driver — the analog of
  the synthetic sine. Kept in **linear** units (non-negative, spiky, many genuine zero-rain days;
  no log taken on the driver).
- `response` — **daily discharge (L/s)**, floored strictly positive (see below). Analyzed in
  **log space** (`log(response)`) by 03/04 — strictly positive, right-skewed, multiplicative.
- `zero_flow` — logical; the day's raw discharge was ≤ 0 (no-flow) and was floored.
- `climate_interp` / `response_interp` — logical; the day is MacroSheds **gap-fill** (modelled,
  `ms_interp == 1`) rather than observed.
- `phase` — ordered factor `pre → disturbance → post` (there is **no** synthetic-style
  `transition` phase); `disturbance` is one day, or a short span at `fernow_WS-5`.
- `site`, `class` — the site tag and disturbance class (`pulse` / `chronic`).

## Channel mapping (the modelling choice)

`climate = precip`, `response = discharge`. Rationale, and how it maps onto the synthetic model:

- **Undisturbed driver.** The synthetic `climate` is never disturbed. **Precipitation is the only
  channel that preserves this** — it is exogenous. Discharge is itself altered by the disturbance
  (the *response*), so it can never be the driver here.
- **Regular, gap-free sampling.** Both peak-detection (III) and Morlet coherence (IV) need
  evenly-sampled, gap-free series. Daily precip & discharge qualify at all five sites; stream
  **chemistry (NO3) does not** (sparse, irregular grab samples) — so chemistry is deliberately
  **excluded** from this leg, even though it is the more canonical disturbance-biogeochem signal.
- **log space.** Discharge is strictly positive, right-skewed, seasonal → analyze `log(response)`,
  the direct analog of the synthetic `log(response)`.

Consequences (and the point): daily precip is **episodic, not a clean annual sine**, so the
seasonal band coupling is genuinely weaker/looser than the synthetic ideal — which is exactly the
"less-deterministic synchrony" the methods are being tested against.

## The five sites (source: `../qualifying_disturbances.csv`)

| Site | Disturbance | Class | Date | Series span (after trim) |
|------|-------------|-------|------|--------------------------|
| `hbef_w1` | chemical_addition (Ca) | pulse | 1999-11-01 | 1983-11 → 2022-11 |
| `hjandrews_GSWS10` | debris_flow | pulse | 1986-02-01 | 1976-02 → 2017-04 |
| `fernow_WS-3` | chemical_addition (acidification) | chronic | 1989-01-01 | 1971-01 → 2019-12 |
| `fernow_WS-5` | timber_harvest (clearcut) | pulse | 2007-01→04 | 1991-01 → 2019-04 |
| `santa_barbara_MC06` | burn (Jesusita fire) | pulse | 2009-05-10 | 2000-10 → 2021-05 |

Two are **chemistry** manipulations (`w1`, `WS-3`) that barely touch water yield → the
precip→discharge lens expects a near-null there (informative about what this mapping can/can't
see). `WS-5` (clearcut) is the classic water-yield disturbance — the site with the most for this
lens to find.

## How the series are built (`prep_ms_series.Rmd`)

1. **Trim** each site's screened window to the span where **both** precip and discharge exist —
   this trims `MC06` (precip starts 2000) and `hjandrews` (precip ends 2017).
2. **Daily grid.** Both channels are already gap-free inside the window (verified) — **no
   interpolation** is done here.
3. **Zero-flow floor.** `response = pmax(discharge, eps)`, `eps = ½ × 1st-percentile of positive
   flow`. A robust floor (the raw minimum can be a ~3e-5 L/s trace that would put no-flow days at
   an absurd log outlier). `response` is written to **6 dp** (4 dp would zero-out sub-1e-4 flows →
   `-Inf` in log). Floored days are flagged `zero_flow`.
4. **Phase.** `pre` before the disturbance start, `disturbance` within its span, `post` after.

## Provenance — `ms_interp` (read before trusting any coupling)

Discharge is **100% observed** at every site — the `response` channel is real throughout.
Precipitation is partly **MacroSheds modelled infill**:

| Site | precip % modelled |
|------|-------------------|
| hbef_w1 / fernow_WS-3 / fernow_WS-5 | ~0% |
| hjandrews_GSWS10 | ~13% |
| **santa_barbara_MC06** | **~88%** |

**`MC06`'s driver is almost entirely modelled** — any MC06 precip→discharge coupling is coupling
to a *modelled* driver, not observed rain. Analyses 03/04 surface this table and caveat MC06
throughout.

## C→Q leg (analysis 05) — the noisier coupling

Analysis 05 is Analysis IV's wavelet moved to **nitrate concentration → discharge**. Key
differences from P→Q, all forced by nitrate being **irregular grab samples** (not daily):

- **`<site>__cq_native.csv` schema:** `time` (bin index), `date`, `Q` (L/s), `C` (nitrate, mg N/L),
  `C_interp` / `Q_interp` (bin had no real sample → interpolated), `phase`, `step_days`, `site`.
  `Q`/`C` are raw; analysis 05 takes `climate = log(Q)`, `log_response = log(C)` (**log–log C–Q
  power law**; the pre-fit predictor slope is the C–Q exponent `b`).
- **Native-cadence grids, not a common grid.** Each site's bin = its median NO3 interval
  (`step_days`: hbef 7, hjandrews 21, fernow ×2 7). Each bin = **mean of the real samples in it**
  (aggregation, *not* interpolation — nothing invented where data exist); only empty bins are
  linearly filled and flagged (`C_interp`). This avoids fabricating the sparser records.
- **Shared physical-time axis.** The wavelet runs in native steps (`dt = 1`), then every period is
  rescaled to **days** (`Period × step_days`) and every lag to days (`angle/2π × 365`), so all four
  panels share one period-in-days axis and one annual band despite different grids. `EF` and the
  block size are `√2·365/step` and `365/step` steps. **A coarser grid can't resolve short periods**
  — hjandrews (21-day) goes blank below ~3 months; shown, not hidden. `lowerPeriod = 4` steps
  (WaveletComp's coherence-smoothing floor), so the low edge is ≈ `4·step_days`.
- **Four sites; MC06 dropped** (its nitrate is 44% empty months + a ~4-yr gap → a wavelet would
  read interpolation). Per-site **n, cadence, and % interpolated** are reported on every panel.
- Coherence saturates again (≈0.9–1) — the P→Q lesson carries over; C–Q signal, if any, is in
  cross-power / band amplitude, now against a genuinely weaker, noisier coupling.

## Analysis notes / gotchas

- **`P = 365` (annual band).** The single tuning constant in each analysis; everything derives
  from it (peak-smoothing/`gap` in III; the wavelet `EF`, band rows, `lower/upperPeriod`, and the
  `resid` lag scan in IV). Changed from the synthetic `P = 500`.
- **`log(response)` needs the floor.** Without the zero-flow floor + 6-dp rounding, `log(0)` →
  `-Inf` → `sd = NaN` → WaveletComp errors. Intermittent `MC06` (many floored days) shows a heavy
  left tail in log space — itself a real intermittent-stream signature.
- **No `transition` phase.** Real phases are `pre / disturbance / post`. The wavelet's settled
  windows are buffered `EF ≈ 516 d` from **both** the series ends and the disturbance, so ~1.4 yr
  of immediate post-disturbance recovery is excluded from the settled `post` average.
- **`loess.span = 0` is fine.** The analyzed wavelet range (4–730 d — wide enough to show the
  week/month/year guide lines) sits *below* any decadal trend, so no detrending is needed and the
  `×sd` amplitude rescale stays exact.
- **Peak method frays on the driver.** Smoothed daily precip gives an unstable "period" (fernow
  is bimodal ~280 d; `MC06` degenerates to a single peak → `period_climate = NA`). Discharge, by
  contrast, shows a clean ~annual peak (356–370 d). This split is a core finding of III.
- **Coherence saturates.** On real data WaveletComp's normalized, smoothed coherence sits ≈0.9–1
  across the whole record at every site — as the synthetic *Desync* case warned, coherence alone
  is not a usable discriminator; cross-power & band amplitude carry the real story (IV).
- **Deterministic.** Given the tracked CSVs, 03/04 reproduce every number on re-knit
  (`make.pval = FALSE`, no RNG). Rebuild the CSVs by knitting `prep_ms_series.Rmd` (needs the
  local `../data` pull).
