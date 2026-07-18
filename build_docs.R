# Render every published report and copy the self-contained HTML to ./docs for GitHub Pages.
# Run from the repo root:  Rscript build_docs.R
# (The per-file `knit:` YAML field does the same copy on the RStudio Knit button; this script is
#  the CLI / rebuild-all path, since rmarkdown::render() does not consult that field.)
# Note: the syn_data sources are tracked; plots/eda_timeseries.Rmd is a git-ignored notebook whose
# HTML is nonetheless published here — it needs the local data/ CSVs to render.
docs <- "docs"
dir.create(docs, showWarnings = FALSE, recursive = TRUE)

rmds <- c(# idealised leg (syn_data) — reproducible from git (deterministic generator + tracked CSVs)
          "syn_data/make_synthetic.Rmd",
          "syn_data/analysis_01_central_tendency.Rmd",
          "syn_data/analysis_02_wavelet_synchrony.Rmd",
          # real-data EDA — git-ignored notebook, needs the local data/ CSVs
          "plots/eda_timeseries.Rmd",
          # real-data leg (ms_data) — analyses 03/04 read the tracked *__series_daily.csv, so they
          # build from git alone; prep_ms_series rebuilds those CSVs and needs the LOCAL pull
          # (../data + ../qualifying_disturbances.csv), so it only renders where those exist.
          "ms_data/prep_ms_series.Rmd",
          "ms_data/analysis_03_central_tendency.Rmd",
          "ms_data/analysis_04_wavelet_synchrony.Rmd")

for (rmd in rmds) {
  # one failing render (e.g. a prep/EDA notebook without the local pull) must not block the rest.
  out <- tryCatch(rmarkdown::render(rmd, quiet = TRUE),   # renders next to the .Rmd (self_contained: true)
                  error = function(e) { message("SKIPPED ", rmd, " — ", conditionMessage(e)); NULL })
  if (is.null(out)) next
  file.copy(out, docs, overwrite = TRUE)
  message("published ", basename(out))
}
