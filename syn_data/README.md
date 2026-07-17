# `syn_data/` — Synthetic disturbance–response datasets & analyses

Idealized synthetic time series modeling an **exogenous climate driver and a lagged response
across a regime shift**, plus two self-contained analysis notebooks that characterize them. One
generator produces three CSVs; two notebooks analyze them; all three notebooks knit to
self-contained HTML that is also published for GitHub Pages.

> For the full generative spec (every parameter, exact phase bands, correlation values) see
> [`CLAUDE.md`](CLAUDE.md). This README is the human-facing orientation to the folder.

---

## At a glance

| File | What it is |
|------|------------|
| `make_synthetic.Rmd` | **Generator** — builds the three datasets (seed `42`, deterministic). |
| `synthetic_shift.csv` · `synthetic_desync.csv` · `synthetic_return.csv` | The three datasets (4250 / 4250 / 4500 rows). |
| `analysis_01_central_tendency.Rmd` | **Analysis I** — log-space moments + peak-derived synchrony, `pre` vs `post`. |
| `analysis_02_wavelet_synchrony.Rmd` | **Analysis II** — Morlet wavelet coherence / cross-power / phase-lag / amplitude at the `P≈500` band. |
| `*.html` | Knitted output of each `.Rmd` (self-contained; also copied to `../docs/`). |
| `CLAUDE.md` | Authoritative design spec for the datasets. |
| `pre_loglinear/` | Archived **old** datasets from the previous independent-sine model (superseded). |

Every `.Rmd` reads the CSVs from this directory (`.`); the datasets are the single interface
between the generator and the analyses.

---

## The datasets

Three CSVs, identical schema, one row per time step:

| Column | Meaning |
|--------|---------|
| `time` | integer step, `1 … N` |
| `climate` | exogenous driver — clean sine (period `P=500`) + constant **25%** Gaussian noise. Never disturbed. |
| `response` | **lagged log-linear transform of the realized noisy climate**: `log(response) = β·climate(t − lag) + noise`, i.e. `response = R0·exp(…)`. Strictly positive, right-skewed, multiplicative (log-space) noise. |
| `phase` | ordered factor `pre → disturbance → transition → post` (`disturbance` is the **single** step `td=2000`). |

`response` is a genuine product of `climate` — climate's 25% noise propagates (lagged) into it.
**Analyze in log space** (`log(response)`) to recover the clean lagged-sine structure.

### The three scenarios

All share one exogenous noisy `climate` and a single-step disturbance at `td=2000`; only
`response` is disturbed. They differ in what the disturbance does to the terminal regime:

| Scenario | N | Trajectory | Mechanism | Terminal state |
|----------|--:|-----------|-----------|----------------|
| **Shift**  | 4250 | regime n → **n+1**, permanent | `w` morphs 0→1 over `[2001,2250]`: sensitivity `β 1.0→0.5`, lag `50→10`; disturbance noise decays 30%→0 | settled regime **n+1** |
| **Desync** | 4250 | regime n → **noise**, permanent | `w=0`; signal gain fades 1→0 while noise grows to regime-n matched (log) variance | stationary **lognormal noise** (log-moments matched to regime n) |
| **Return** | 4500 | regime n → **noise** → regime n | `w=0`; gain dips 1→0→1 (trough at 2250, window `[2001,2500]`) — a **transient Desync** | recovered regime **n** |

Only **Shift** morphs parameters (`w>0`); **Desync** and **Return** hold `w=0` and decouple via
the fading signal gain (permanently vs transiently).

---

## The generative model

`response = R0·exp(β·climate(t − lag) + noise)` — a static, convex, strictly-positive transfer
of the lagged realized climate that is a **straight line of slope `β` in log space** (hence
*log-linear*). Regime n uses `β=1.0, lag=50`; regime n+1 uses `β=0.5, lag=10`. Transitions
crossfade two lagged copies of climate in log space via the weight `w`. `make_synthetic.Rmd`
opens with a **Transfer function** panel that plots this map in both real and log units.

Key parameters (`params` chunk): `P=500`, `Ac=1.0`, `beta_n/lag_n = 1.0/50`,
`beta_p/lag_p = 0.5/10`, `R0=1.0`, `clim_sd=0.25` (climate noise, linear), `base_sd=0.05`
(response baseline noise, log space), `dist_sd=0.30` (peak disturbance noise, log space),
`LEN_REG=2000`, `LEN_ADAPT=250`.

---

## The analyses

Two independent passes, deliberately using **different metrics**, each contrasting settled
`pre` vs `post` windows with **block standard errors** (valid under autocorrelation).

**Analysis I — central tendency (`analysis_01`).** The simplest thing that works: per-channel
moments and peak-spacing synchrony, all on `log(response)`. Headline reads:
- **Shift** — log-amplitude ratio `sd(logR)/sd(clim)` drops `1.0→0.5` (= β) and phase lag `50→10`
  while period holds ~P: a coherent **regime shift**, still coupled.
- **Desync** — `sd(log(response))` roughly holds but peak-spacing irregularity blows up and the
  response period detaches from climate's: **variance-preserving decoupling**.
- **Return** — every metric matches end-to-end: **no lasting signature**.

**Analysis II — wavelet synchrony (`analysis_02`).** The native signal-processing view: Morlet
wavelet (via **WaveletComp**) coherence, cross-power, phase→apparent-lag, and amplitude ratio at
the `P≈500` band, resolved through time (COI-buffered). Headline reads:
- **Shift** — coherence & cross-power hold; apparent lag pins `50→10` and amplitude ratio `1.0→0.5`
  far more sharply than peak offsets.
- **Desync** — the *normalized* coherence is **fooled** (stays ≈1 over noise); the decoupling
  registers only in **un-normalized cross-power and band amplitude collapsing to ≈0**. Lesson:
  coherence alone is not enough.
- **Return** — a flat null in the settled view; its transition is a brief, fully-recovering decouple.

---

## Reproduce & build

Requirements: **R** with `rmarkdown, knitr, dplyr, tidyr, readr, ggplot2, patchwork` (all three
notebooks) plus **`WaveletComp`** (analysis_02 only).

- **Regenerate the datasets** — knit `make_synthetic.Rmd` (writes the three CSVs to this
  directory). Output is deterministic given `set.seed(42)` and top-to-bottom execution.
- **Run an analysis** — knit `analysis_01_*.Rmd` or `analysis_02_*.Rmd` (they read the CSVs from
  `.`; regenerate the data first if you changed the generator).
- **Build everything to the published site** — from the **repo root**:
  ```
  Rscript build_docs.R
  ```
  renders all three reports and copies the self-contained HTML into `../docs/`.

Knitting any of these `.Rmd` (via the RStudio **Knit** button) also copies its HTML to `../docs/`
automatically, via the `knit:` field in the YAML header.

---

## Published site

The self-contained HTML reports are published to the repo-root `docs/` folder for **GitHub Pages**
(Settings → Pages → Deploy from a branch → `main` / `/docs`). `docs/index.html` is the landing
page linking the three reports. Once Pages is enabled they are served at
`https://cascade-gub.github.io/cz_disturbance/` (e.g. `…/analysis_02_wavelet_synchrony.html`).

---

## Notes & gotchas

- **Work in log space.** `log(response)` is a lag-shifted copy of `climate`; the matched-lag
  correlation is ≈ +0.99 and `sd(log(response))/sd(climate) ≈ β`. Raw `response` is right-skewed
  with multiplicative noise — don't assume Gaussian/additive residuals in linear units.
- **Desync's terminal block matches regime n in *log* moments only.** At equal log-variance the
  raw-unit variance still differs (terminal ≈ +40%), because regime n is `exp(sinusoid)` and the
  terminal block is `exp(Gaussian)`. Compare Desync in log space.
- **Only Shift changes its lag**, so only Shift's transition carries a real transient frequency
  excursion (a Doppler chirp, documented in the generator's Engine notes); Desync and Return hold
  the lag fixed. The analyses read settled regimes and do not measure the transient chirp.
- `pre_loglinear/` holds the **superseded** independent-sine datasets — kept for reference, not
  used by the current notebooks.
