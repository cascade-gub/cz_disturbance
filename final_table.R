suppressMessages({library(dplyr)})
YRS <- 10
res <- readRDS("verified.rds")

out <- res %>%
  filter(CHEM_OK | Q_OK) %>%
  # A chronic event with no end_date in the record (end coalesced to start) is an
  # ONGOING treatment: its "10 years post" would fall inside active disturbance,
  # so it cannot satisfy "10 yrs post-disturbance". Segregate rather than report.
  mutate(open_ended_chronic = grepl("chronic", classes) & grp_end == grp_start,
         qualifies_on = case_when(CHEM_OK & Q_OK ~ "stream chemistry + discharge",
                                  CHEM_OK ~ "stream chemistry only",
                                  TRUE ~ "discharge only"),
         window_start = grp_start - round(YRS * 365.25),
         window_end   = grp_end   + round(YRS * 365.25)) %>%
  arrange(desc(CHEM_OK), domain, site_code, grp_start) %>%
  select(network, domain, site_code,
         disturbance_start = grp_start, disturbance_end = grp_end,
         n_events_grouped = n_events, disturbance = types, class = classes,
         qualifies_on, window_start, window_end, open_ended_chronic,
         chem_days_pre = chem_pre_n, chem_days_post = chem_post_n,
         chem_pre_bins, chem_post_bins, chem_pre_minmo, chem_post_minmo,
         q_pre_bins, q_post_bins)

main <- out
excl <- out %>% filter(open_ended_chronic)

cat("=== QUALIFYING DISTURBANCES ===\n")
print(as.data.frame(main), width = 300)
cat("\n=== NOTE: chronic, no end date in record -> ongoing; post-window overlaps active treatment ===\n")
print(as.data.frame(excl %>% select(domain, site_code, disturbance_start, disturbance, class)))

safe_write <- function(d, f) {
  ok <- tryCatch({write.csv(d, f, row.names = FALSE); TRUE}, error = function(e) FALSE,
                 warning = function(w) FALSE)
  if (!ok) {
    f <- sub("\\.csv$", "_v2.csv", f)
    write.csv(d, f, row.names = FALSE)
  }
  cat("wrote", f, "\n")
}
safe_write(main, "qualifying_disturbances.csv")
safe_write(excl, "excluded_ongoing_chronic.csv")
cat(nrow(main), "disturbances at", n_distinct(paste(main$domain, main$site_code)), "sites\n")
