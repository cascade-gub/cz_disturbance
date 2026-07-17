pkgs <- c("rmarkdown","knitr","ggplot2","tidyr","dplyr","readr","lubridate",
          "scales","patchwork","stringr")
for (p in pkgs) cat(sprintf("%-12s %s\n", p, requireNamespace(p, quietly = TRUE)))
cat("pandoc:", rmarkdown::pandoc_available(), as.character(rmarkdown::pandoc_version()), "\n")
