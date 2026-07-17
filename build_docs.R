# Render every syn_data report and publish the self-contained HTML to ./docs for GitHub Pages.
# Run from the repo root:  Rscript build_docs.R
# (The per-file `knit:` YAML field does the same copy on the RStudio Knit button; this script is
#  the CLI / rebuild-all path, since rmarkdown::render() does not consult that field.)
docs <- "docs"
dir.create(docs, showWarnings = FALSE, recursive = TRUE)

rmds <- c("syn_data/make_synthetic.Rmd",
          "syn_data/analysis_01_central_tendency.Rmd",
          "syn_data/analysis_02_wavelet_synchrony.Rmd")

for (rmd in rmds) {
  out <- rmarkdown::render(rmd, quiet = TRUE)   # renders next to the .Rmd (self_contained: true)
  file.copy(out, docs, overwrite = TRUE)
  message("published ", basename(out))
}
