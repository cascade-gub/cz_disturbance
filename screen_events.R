# Screen disturbance events against MacroSheds data coverage (catalog extents).
suppressMessages({library(macrosheds); library(dplyr); library(readxl); library(lubridate)})

YRS <- 10          # required clear years pre and post
GAP_DAYS <- 31     # events within a month of each other are one disturbance

rec <- read_excel("disturbance_record.xlsx", sheet = "Sheet1") %>%
  select(network, domain, site_code, watershed_type, disturbance_source,
         disturbance_type, disturbance_def, start_date, end_date) %>%
  filter(!is.na(site_code))

cat("total rows:", nrow(rec), "\n")

# start_date is a character column: Excel serial numbers, except pre-1900 dates
# (santee/calhoun historic agriculture) which Excel stores as literal "YYYY-MM-DD" text.
xl_date <- function(x) {
  num <- suppressWarnings(as.numeric(x))
  out <- as.Date(num, origin = "1899-12-30")
  txt <- is.na(num) & !is.na(x)
  out[txt] <- as.Date(x[txt])
  out
}

ev <- rec %>%
  filter(!is.na(start_date)) %>%
  mutate(start = xl_date(start_date),
         end   = as.Date(end_date),
         end   = pmax(coalesce(end, start), start))
stopifnot(!any(is.na(ev$start)), !any(is.na(ev$end)))
cat("start range:", format(range(ev$start)), "\n")
cat("disturbance events (rows with start_date):", nrow(ev), "\n")
cat("sites with >=1 event:", n_distinct(paste(ev$domain, ev$site_code)), "\n")

# ---- group events within GAP_DAYS into single disturbances ----
grp <- ev %>%
  arrange(domain, site_code, start, end) %>%
  group_by(domain, site_code) %>%
  mutate(gap = as.numeric(start - lag(cummax(as.numeric(end)))),
         new_grp = is.na(gap) | gap > GAP_DAYS,
         grp_id = cumsum(new_grp)) %>%
  group_by(domain, site_code, grp_id) %>%
  summarise(network = first(network),
            watershed_type = first(watershed_type),
            grp_start = min(start),
            grp_end   = max(end),
            n_events  = n(),
            types     = paste(unique(disturbance_def), collapse = "; "),
            sources   = paste(unique(disturbance_source), collapse = "; "),
            classes   = paste(unique(disturbance_type), collapse = "; "),
            .groups = "drop")
cat("grouped disturbances:", nrow(grp), "\n")

# ---- neighbouring-disturbance clearance ----
grp <- grp %>%
  arrange(domain, site_code, grp_start) %>%
  group_by(domain, site_code) %>%
  mutate(prev_end   = lag(cummax(as.numeric(grp_end))),
         prev_end   = as.Date(prev_end, origin = "1970-01-01"),
         next_start = lead(grp_start),
         pre_clear_yrs  = ifelse(is.na(prev_end), Inf,
                                 as.numeric(grp_start - prev_end) / 365.25),
         post_clear_yrs = ifelse(is.na(next_start), Inf,
                                 as.numeric(next_start - grp_end) / 365.25)) %>%
  ungroup() %>%
  mutate(dist_clear = pre_clear_yrs >= YRS & post_clear_yrs >= YRS)

cat("disturbances with", YRS, "clear yrs of NO other disturbance both sides:",
    sum(grp$dist_clear), "\n")

# ---- data coverage from MacroSheds catalog (extent only; screening) ----
cov <- ms_var_catalog %>%
  filter(chem_category == "stream_conc" | variable_code == "discharge") %>%
  mutate(kind = ifelse(variable_code == "discharge", "q", "chem")) %>%
  group_by(domain, site_code, kind) %>%
  summarise(first = min(as.Date(first_record)), last = max(as.Date(last_record)),
            obs = sum(observations), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = kind, values_from = c(first, last, obs))

cand <- grp %>%
  left_join(cov, by = c("domain", "site_code")) %>%
  mutate(in_ms = !is.na(first_chem) | !is.na(first_q),
         chem_pre  = as.numeric(grp_start - first_chem) / 365.25,
         chem_post = as.numeric(last_chem - grp_end) / 365.25,
         q_pre     = as.numeric(grp_start - first_q) / 365.25,
         q_post    = as.numeric(last_q - grp_end) / 365.25,
         chem_ok = !is.na(chem_pre) & chem_pre >= YRS & chem_post >= YRS,
         q_ok    = !is.na(q_pre)    & q_pre    >= YRS & q_post    >= YRS)

# sites in disturbance record but not in MacroSheds at all
missing <- cand %>% filter(!in_ms) %>% distinct(network, domain, site_code)
cat("\nsites w/ events NOT in MacroSheds catalog:", nrow(missing), "\n")
print(as.data.frame(missing))

pass <- cand %>% filter(dist_clear, chem_ok | q_ok)
cat("\n=== CANDIDATES (catalog extent screen) ===\n")
cat("n =", nrow(pass), " sites =", n_distinct(paste(pass$domain, pass$site_code)), "\n")
print(as.data.frame(pass %>%
  select(network, domain, site_code, grp_start, grp_end, n_events, classes, types,
         chem_pre, chem_post, q_pre, q_post, chem_ok, q_ok) %>%
  mutate(across(c(chem_pre, chem_post, q_pre, q_post), ~round(.x, 1)))), width = 250)

saveRDS(list(grp = grp, cand = cand, pass = pass), "screen.rds")
write.csv(pass, "candidates_screen.csv", row.names = FALSE)
