suppressMessages({library(macrosheds); library(dplyr)})
root <- "C:/Users/nic/cz_disturbance/ms_data"
sites <- tribble(
  ~domain, ~site_code,
  "hbef", "w1", "hjandrews", "GSWS10", "fernow", "WS-3", "fernow", "WS-5",
  "santa_barbara", "MC06")

for (d in unique(sites$domain)) {
  sc <- ms_load_product(root, prod = "stream_chemistry", domains = d, warn = FALSE) %>%
    filter(ms_interp == 0, !is.na(val), site_code %in% sites$site_code[sites$domain == d])
  cat("\n====", d, "stream_chemistry vars by n obs ====\n")
  print(sc %>% group_by(site_code, var) %>%
    summarise(n = n(), first = min(date), last = max(date), .groups = "drop") %>%
    filter(grepl("NO3|spCond|SpCond|cond", var, ignore.case = TRUE) | n > 500) %>%
    arrange(site_code, desc(n)) %>% as.data.frame())
}
cat("\n==== products present per domain ====\n")
for (d in unique(sites$domain)) {
  cat(d, ":", paste(list.files(file.path(root, "v2", d)), collapse = ", "), "\n")
}
