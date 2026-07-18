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
  **`WaveletComp`** (`analysis_02`/`analysis_04` only). The `macrosheds` + `arrow` packages are
  needed only to (re)build the raw pull, not to run the analyses.

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
does the same copy via its `knit:` YAML field.
