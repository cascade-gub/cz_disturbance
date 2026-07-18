# CLAUDE.md — cz_disturbance (repo root)

Disturbance–response **synchrony** study. The repo is a **methods comparison**: two different
approaches to quantifying synchrony between a driver and a lagged response — central-tendency /
peak-spacing vs Morlet wavelet coherence — and their respective strengths and weaknesses. Each
approach is a different, partial perspective on the same coupled system.

- **`syn_data/`** — the *idealised, deterministic* leg. A synthetic generator
  (`make_synthetic.Rmd`) + two analyses (`analysis_01_central_tendency`,
  `analysis_02_wavelet_synchrony`) that establish the **theory and expectations** for what each
  method can and cannot see. See `syn_data/CLAUDE.md` for the full generative spec.
- **`ms_data/`** — the *real-data* leg. The same two approaches applied to real **MacroSheds**
  watershed data (`analysis_03_central_tendency`, `analysis_04_wavelet_synchrony`), testing how
  they hold up under **less-deterministic, real-world synchrony**. Channel mapping:
  `climate` = precipitation (undisturbed driver), `response` = discharge; daily resolution,
  seasonal period `P≈365`. See `ms_data/CLAUDE.md`.

The real-data messiness (episodic driver, trends, gaps, non-stationarity) is the **stress-test
bench** for the methods, not a defect — how each method responds is the finding.

## Toolchain

- **R:** `C:\Program Files\R\R-4.6.1\bin\Rscript.exe` (not on `PATH`; call by full path).
  - Bash tool: `"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" script.R`
  - PowerShell: `& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" script.R`
- Packages: `rmarkdown, knitr, dplyr, tidyr, readr, ggplot2, patchwork` (all notebooks) plus
  **`WaveletComp`** (the wavelet notebooks `analysis_02`/`04`/`05`). The `macrosheds` + `arrow`
  packages are needed only to (re)build the raw pull, not to run the analyses.
- **Shared wavelet kernel:** the three wavelet notebooks source **`R/wavelet_kernel.R`** (tracked
  R at repo root, `source("../R/wavelet_kernel.R")` from either leg dir) for the **whole** helper
  family — `wco`, `band_trace`, `block_stat`, `wave_summary`, `plot_*`, `fit_pre`, `resid_trace`,
  … plus the red-noise cross-power **significance** helpers. **Edit once, not three times.** Each
  notebook builds its copies via `make_wavelet_kernel(cfg)` and injects them with `list2env(...,`
  `environment())`, so call sites stay bare; `cfg` holds the handful of things that differ (band
  `P_days`, native `step`, phase gating, `transition` phase, axis, guide marks, …). The
  significance mask is analytic AR(1) red noise (Torrence & Compo 1998), so knits stay deterministic
  (no `make.pval`, no RNG). Add a feature or fix a bug in the kernel and all three inherit it.

## Data layout & git

- `ms_data/v2/` — raw MacroSheds download (~1.8 GB `.feather` + shapefiles). **Git-ignored.**
- `data/` — tidy per-site daily CSVs from `pull_data.R` (git-ignored working data).
- `qualifying_disturbances.csv` — the 5 screened disturbance events + pre/post windows (the
  real-data `phase` source of truth), produced by the `screen_events`/`verify`/`extend_windows`
  pipeline at repo root.
- The **tracked interface** for each analysis leg is small CSVs beside the notebooks:
  `syn_data/synthetic_*.csv` and `ms_data/*__series_daily.csv`. Raw feathers stay ignored.

## Build

From repo root: `"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" build_docs.R` renders every
notebook and copies the self-contained HTML into `docs/` for GitHub Pages
(`docs/index.html` is the landing page). Knitting an `.Rmd` via the RStudio **Knit** button
does the same copy via its `knit:` YAML field. `build_docs.R` renders **both legs** — `syn_data`
(I/II) and `ms_data` (III/IV) — plus `plots/eda_timeseries.Rmd`. The `ms_data` notebooks and their
`*__series_daily.csv` are **tracked** and published; only the raw `ms_data/v2/` pull and `plots/`
stay git-ignored. `prep_ms_series` and the EDA need the local pull (`data/`) to *re-render*, so
`build_docs.R` **skips (does not fail on)** any notebook whose local inputs are absent — the
tracked III/IV still build from the committed CSVs.

## Report styling — **standard: dark page + white-matted figures**

Every knitted report is a **dark-mode page** (`theme: darkly`, `highlight: breezedark`) with
**figures in their native (light) colour scheme, matted on white** so they stay legible against
the dark page. Do **not** re-theme plots dark to match the page — mat them instead. New notebooks
must follow this. The mechanism (see any `syn_data`/`ms_data` setup chunk):

- **PNG canvas:** `knitr::opts_chunk$set(dev.args = list(bg = "white"))` — native light figures.
- **ggplot theme:** plain `theme_set(theme_bw(base_size = 11) + theme(legend.position = "bottom"))`
  — native colours, no dark palette overrides.
- **White matting:** a `` ```{css figure-matting, echo=FALSE} `` chunk giving
  `.main-container img, div.figure img` a white background, padding, rounded corners and a soft
  shadow — the figures read as white cards on the dark page.
- **Base-R plots** (e.g. `WaveletComp` heatmaps): `par(bg = "white", fg = "black", col.axis =`
  `"black", col.lab = "black", col.main = "black")` so they mat the same way.
- **Error bars are black** (`colour = "black"`) for contrast against the white figure background;
  subtle in-figure reference lines (dashed disturbance markers) stay light grey.
