suppressMessages({library(macrosheds)})
root <- "C:/Users/nic/cz_disturbance/ms_data"
dir.create(root, showWarnings = FALSE, recursive = TRUE)
doms <- c("fernow", "hbef", "hjandrews", "santa_barbara", "santee")
ms_download_core_data(macrosheds_root = root, domains = doms, quiet = FALSE)
cat("\nDOWNLOAD DONE\n")
print(list.files(root, recursive = FALSE))
