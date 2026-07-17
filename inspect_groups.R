suppressMessages({library(dplyr)})
x <- readRDS("screen.rds")
grp <- x$grp; cand <- x$cand

cat("=== groups that merged >1 event (within-1-month clustering) ===\n")
print(as.data.frame(grp %>% filter(n_events > 1) %>%
  select(domain, site_code, grp_start, grp_end, n_events, types)), width = 200)

cat("\n=== of the 48 disturbance-clear groups, why did 35 fail the data screen? ===\n")
f <- cand %>% filter(dist_clear, !(chem_ok | q_ok)) %>%
  mutate(across(c(chem_pre, chem_post, q_pre, q_post), ~round(.x, 1)))
print(as.data.frame(f %>% select(domain, site_code, grp_start, grp_end, types,
                                 chem_pre, chem_post, q_pre, q_post)), width = 200)

cat("\n=== near-miss: chem fails only on one side by <2 yrs ===\n")
print(as.data.frame(f %>% filter(!is.na(chem_pre), pmin(chem_pre, chem_post) > 8) %>%
  select(domain, site_code, grp_start, types, chem_pre, chem_post)))
