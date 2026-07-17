suppressMessages({library(macrosheds); library(dplyr)})
root <- "C:/Users/nic/cz_disturbance/ms_data"
x <- ms_load_product(root, prod = "stream_chemistry", domains = "hbef", warn = FALSE)
cat("cols:\n"); print(names(x)); print(head(as.data.frame(x), 3))
cat("\nms_status:"); print(table(x$ms_status, useNA = "ifany"))
cat("ms_interp:"); print(table(x$ms_interp, useNA = "ifany"))
q <- ms_load_product(root, prod = "discharge", domains = "hbef", warn = FALSE)
cat("\ndischarge cols:\n"); print(names(q)); print(head(as.data.frame(q), 2))
