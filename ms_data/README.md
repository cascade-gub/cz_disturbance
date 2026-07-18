# `ms_data/` — Real-data (MacroSheds) disturbance–response analyses

The **real-data leg** of a two-part study. `../syn_data/` runs two different synchrony-quantifying
methods on an *idealised* system to set up the theory; this folder runs the **same two methods on
real MacroSheds watersheds** to see how they hold up under messy, less-deterministic coupling. It
is a **methods comparison** — peak/moments vs Morlet wavelet — not a hunt for a known answer. The
real-data messiness is the stress-test bench, and how each method copes is the result.

> Full spec + every gotcha: [`CLAUDE.md`](CLAUDE.md). Idealised leg: [`../syn_data/`](../syn_data/).
> Repo overview + R path: [`../CLAUDE.md`](../CLAUDE.md).

---

## At a glance

| File | What it is |
|------|------------|
| `prep_ms_series.Rmd` | **Prep** (analog of `syn_data/make_synthetic.Rmd`) — carves the real pull into the P→Q and C→Q datasets. |
| `*__series_daily.csv` (×5) | **P→Q** datasets — daily precip + discharge, one per screened disturbance (tracked). |
| `*__cq_native.csv` (×4) | **C→Q** datasets — discharge + nitrate on each site's native NO3 cadence (MC06 dropped; tracked). |
| `analysis_03_central_tendency.Rmd` | **Analysis III** — analysis I's method (log-space moments + peak-derived synchrony) on real precip→discharge. |
| `analysis_04_wavelet_synchrony.Rmd` | **Analysis IV** — analysis II's method (Morlet wavelet) on **P→Q** at the annual band. |
| `analysis_05_cq_wavelet_synchrony.Rmd` | **Analysis V** — the same wavelet on the noisier **C→Q** (nitrate→discharge) coupling. |
| `*.html` | Knitted output of each `.Rmd` (self-contained; also copied to `../docs/`). |
| `CLAUDE.md` | Authoritative spec for this folder. |
| `v2/` | The raw MacroSheds download (~1.8 GB feathers + shapefiles) — **git-ignored**. |

Each `.Rmd` reads its `*__series_daily.csv` (03/04) or `*__cq_native.csv` (05) from this directory
(`.`); those CSVs are the single interface between the local pull and the analyses (feathers are
never touched by the analyses).

---

## The datasets

Five daily CSVs, identical schema, one row per day (full schema in `CLAUDE.md`):

| Column | Meaning |
|--------|---------|
| `time` | integer day index `1 … N` (`dt = 1` day) |
| `date` | calendar date |
| `climate` | **daily precipitation (mm)** — the exogenous, undisturbed driver (linear units) |
| `response` | **daily discharge (L/s)**, floored positive — analyzed as `log(response)` |
| `zero_flow`, `climate_interp`, `response_interp` | no-flow / modelled-infill flags |
| `phase` | ordered `pre → disturbance → post` |
| `site`, `class` | site tag; `pulse` / `chronic` |

**Channel mapping:** `climate = precip` (the only *undisturbed* channel), `response = discharge`
(the disturbed variable). Stream chemistry (NO3) is deliberately excluded — too sparse/irregular
for peak-finding and wavelets. Daily precip is **episodic, not a clean sine**, so the seasonal
coupling is looser than the synthetic ideal — that is the point.

### The five disturbances

| Site | Disturbance | Class | Note |
|------|-------------|-------|------|
| `hbef_w1` | Ca addition | pulse | chemistry target → expect near-null in discharge |
| `hjandrews_GSWS10` | debris flow | pulse | physical channel disturbance; precip ~13% modelled |
| `fernow_WS-3` | chronic acidification | chronic | open-ended; chemistry target |
| `fernow_WS-5` | clearcut | pulse | classic water-yield disturbance — most for this lens to see |
| `santa_barbara_MC06` | burn | pulse | intermittent flow; **precip ~88% modelled** — hardest case |

---

## The analyses

The **same methods** as `syn_data`, on real watersheds, over two couplings of decreasing
causality. Each contrasts settled `pre` vs `post` windows with block standard errors (valid under
autocorrelation), now resting on **decades of annual cycles** rather than the synthetic ~2.

**Analysis III — P→Q central tendency (`analysis_03`).** Analysis I's peak/moments method on real
precip→discharge (five sites). Headline: discharge shows a clean ~annual peak, but the episodic
precip **driver** frays the method (unstable "period"; `MC06` degenerates), the amplitude ratio no
longer maps to a clean β, and pre/post differences tangle with decadal climate variability — a
**partial, noisy** read.

**Analysis IV — P→Q wavelet synchrony (`analysis_04`).** Analysis II's Morlet method at the
`P≈365 d` band (five sites). Headline: the synthetic lesson reappears **unprompted** — normalized
coherence **saturates near 1** everywhere (not a usable discriminator), and whatever real change
exists is read from **un-normalized cross-power and band amplitude**.

**Analysis V — C→Q wavelet synchrony (`analysis_05`).** The same wavelet moved to
**discharge→nitrate concentration** — a deliberately *less causal, noisier* coupling (four sites;
MC06 dropped). Nitrate is grab-sampled, so each site runs on its **native cadence** and periods are
rescaled to a **shared physical-time axis** (a coarse-sampled site simply can't resolve short
periods — shown, not hidden). Both channels in log space (C–Q power law); the pre-fit slope is the
C–Q exponent. Headline: coherence saturates as before, so the weaker C–Q signal — if seasonally
organised at all — lives in cross-power/amplitude, now most relevant at the two sites whose
disturbance hits nitrogen directly (hbef Ca, fernow WS-3 acidification).

Together: the methods on progressively harder couplings — idealised (I/II) → near-causal P→Q
(III/IV) → noisy C→Q (V) — each a different, partial view, incomplete in its own way.

---

## Provenance caveat (`ms_interp`)

Discharge is **fully observed** everywhere; precipitation is partly **MacroSheds modelled infill**
— negligible at hbef/fernow, ~13% at hjandrews, and **~88% at `MC06`** (a modelled driver). Read
MC06's coupling accordingly; 03/04 surface this table and repeat the caveat.

---

## Reproduce & build

Requirements: **R** (`../CLAUDE.md` has the path) with `rmarkdown, knitr, dplyr, tidyr, readr,
ggplot2, patchwork` (all notebooks) plus **`WaveletComp`** (analysis_04 / analysis_05).

- **Run an analysis** — knit `analysis_03/04_*.Rmd` (read the five `*__series_daily.csv`) or
  `analysis_05_*.Rmd` (reads the four `*__cq_native.csv`); they need nothing else (no local pull).
  IV takes ~4 min (five daily wavelet transforms); V is fast (four native-cadence transforms).
- **Regenerate the datasets** — knit `prep_ms_series.Rmd` (builds both the P→Q and C→Q CSVs). This
  one **needs the local MacroSheds pull** (`../data/*` + `../qualifying_disturbances.csv`), git-ignored.
- **Build everything to the published site** — from the **repo root**:
  ```
  Rscript build_docs.R
  ```
  renders both legs and copies the self-contained HTML into `../docs/` (it *skips* rather than
  fails on any notebook whose local inputs are absent). Knitting any `.Rmd` via the RStudio
  **Knit** button also copies its HTML to `../docs/` via the `knit:` YAML field.

The reports are published to GitHub Pages alongside the synthetic leg
(`../docs/index.html` links all of I–IV).
