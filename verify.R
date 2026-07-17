# Verify candidate disturbances against ACTUAL MacroSheds observations:
# catalog extent can hide multi-year gaps, so require every one of the 10 pre-
# and 10 post- year-bins to actually contain observations.
suppressMessages({library(macrosheds); library(dplyr); library(readxl); library(tidyr)})

YRS <- 10; GAP_DAYS <- 31
root <- "C:/Users/nic/cz_disturbance/ms_data"
doms <- c("fernow", "hbef", "hjandrews", "santa_barbara", "santee")

xl_date <- function(x) {
  num <- suppressWarnings(as.numeric(x)); out <- as.Date(num, origin = "1899-12-30")
  txt <- is.na(num) & !is.na(x); out[txt] <- as.Date(x[txt]); out
}

build_groups <- function(drop_pre1900 = FALSE) {
  ev <- read_excel("disturbance_record.xlsx", sheet = "Sheet1") %>%
    filter(!is.na(site_code), !is.na(start_date)) %>%
    mutate(start = xl_date(start_date), end = as.Date(end_date),
           end = pmax(coalesce(end, start), start))
  if (drop_pre1900) ev <- ev %>% filter(as.numeric(format(start, "%Y")) >= 1900)
  ev %>%
    arrange(domain, site_code, start, end) %>%
    group_by(domain, site_code) %>%
    mutate(gap = as.numeric(start - lag(cummax(as.numeric(end)))),
           grp_id = cumsum(is.na(gap) | gap > GAP_DAYS)) %>%
    group_by(network, domain, site_code, grp_id) %>%
    summarise(grp_start = min(start), grp_end = max(end), n_events = n(),
              types = paste(unique(disturbance_def), collapse = "; "),
              classes = paste(unique(disturbance_type), collapse = "; "),
              sources = paste(unique(disturbance_source), collapse = "; "),
              .groups = "drop") %>%
    arrange(domain, site_code, grp_start) %>%
    group_by(domain, site_code) %>%
    mutate(prev_end = as.Date(lag(cummax(as.numeric(grp_end))), origin = "1970-01-01"),
           next_start = lead(grp_start),
           pre_clear_yrs  = ifelse(is.na(prev_end), Inf, as.numeric(grp_start - prev_end)/365.25),
           post_clear_yrs = ifelse(is.na(next_start), Inf, as.numeric(next_start - grp_end)/365.25)) %>%
    ungroup() %>%
    mutate(dist_clear = pre_clear_yrs >= YRS & post_clear_yrs >= YRS)
}

grp <- build_groups(FALSE)

# ---- real observations ----
# ms_load_product returns no `domain` column, so load per-domain and tag.
# ms_interp == 1 rows are gap-filling placeholders (val is NA) - only real samples count.
load_obs <- function(prod) {
  bind_rows(lapply(doms, function(d) {
    ms_load_product(root, prod = prod, domains = d, warn = FALSE) %>%
      filter(ms_interp == 0, !is.na(val)) %>%
      mutate(domain = d, date = as.Date(date)) %>%
      distinct(domain, site_code, date)
  }))
}
chem <- load_obs("stream_chemistry")
q    <- load_obs("discharge")
cat("chem obs-days:", nrow(chem), " discharge obs-days:", nrow(q), "\n")

# site codes in the disturbance record must match MacroSheds site codes
cand_sites <- grp %>% filter(dist_clear, domain %in% doms) %>% distinct(domain, site_code)
unmatched <- cand_sites %>% anti_join(bind_rows(chem, q) %>% distinct(domain, site_code),
                                      by = c("domain", "site_code"))
cat("candidate sites with NO matching MacroSheds site_code:", nrow(unmatched), "\n")
if (nrow(unmatched)) print(as.data.frame(unmatched))

# Per group, per side: each of the 10 one-year bins must be "monthly-resolved" --
# at least one observation in >= MIN_MONTHS of the 12 calendar months in that bin.
MIN_MONTHS <- 11
bin_check <- function(dat, dom, site, anchor, direction) {
  d <- dat$date[dat$domain == dom & dat$site_code == site]
  empty <- tibble(bins_with_data = 0L, bins_monthly = 0L, n_days = 0L,
                  min_months = 0L, max_gap_d = NA_real_)
  if (!length(d)) return(empty)
  edges <- if (direction == "pre")
    seq(anchor, by = "-1 year", length.out = YRS + 1) else
    seq(anchor, by = "1 year", length.out = YRS + 1)
  lo <- min(edges); hi <- max(edges)
  w <- sort(d[d >= lo & d <= hi])
  if (!length(w)) return(empty)
  # 120 monthly slots across the window; bin k = slots (12k-11)..(12k). Slots are
  # anchored to the disturbance date, so "month" means an elapsed month of the
  # window, not a calendar month (a bin therefore has exactly 12, never 13).
  mo_edges <- seq(lo, by = "1 month", length.out = 12 * YRS + 1)
  slot <- cut(as.numeric(w), breaks = as.numeric(mo_edges),
              include.lowest = TRUE, labels = FALSE)
  filled <- unique(slot[!is.na(slot)])
  bin_of <- ceiling(filled / 12)
  months_per_bin <- tabulate(bin_of, nbins = YRS)   # distinct months filled in each bin
  tibble(bins_with_data = sum(months_per_bin > 0),
         bins_monthly   = sum(months_per_bin >= MIN_MONTHS),
         n_days         = length(w),
         min_months     = min(months_per_bin),      # weakest year in the window
         max_gap_d      = max(as.numeric(diff(c(lo, w, hi)))))
}

cands <- grp %>% filter(dist_clear, domain %in% doms)
res <- cands %>% rowwise() %>%
  mutate(
    c_pre  = list(bin_check(chem, domain, site_code, grp_start, "pre")),
    c_post = list(bin_check(chem, domain, site_code, grp_end,  "post")),
    q_pre  = list(bin_check(q,    domain, site_code, grp_start, "pre")),
    q_post = list(bin_check(q,    domain, site_code, grp_end,  "post"))) %>%
  mutate(chem_pre_bins  = c_pre$bins_monthly,  chem_pre_n  = c_pre$n_days,
         chem_post_bins = c_post$bins_monthly, chem_post_n = c_post$n_days,
         chem_pre_minmo = c_pre$min_months, chem_post_minmo = c_post$min_months,
         chem_pre_gap   = round(c_pre$max_gap_d/365.25, 2),
         chem_post_gap  = round(c_post$max_gap_d/365.25, 2),
         q_pre_bins  = q_pre$bins_monthly,  q_post_bins = q_post$bins_monthly,
         q_pre_minmo = q_pre$min_months, q_post_minmo = q_post$min_months,
         q_pre_n = q_pre$n_days, q_post_n = q_post$n_days) %>%
  select(-c_pre, -c_post, -q_pre, -q_post) %>% ungroup() %>%
  mutate(CHEM_OK = chem_pre_bins == YRS & chem_post_bins == YRS,
         Q_OK    = q_pre_bins == YRS & q_post_bins == YRS)

cat("\n=== VERIFIED: all 20 year-bins have data in >=", MIN_MONTHS, "of 12 months ===\n")
print(as.data.frame(res %>% filter(CHEM_OK | Q_OK) %>%
  select(network, domain, site_code, grp_start, grp_end, n_events, classes, types,
         chem_pre_bins, chem_post_bins, chem_pre_minmo, chem_post_minmo,
         chem_pre_n, chem_post_n, q_pre_bins, q_post_bins, CHEM_OK, Q_OK)), width = 300)

cat("\n=== candidates that FAILED the monthly-resolution test ===\n")
print(as.data.frame(res %>% filter(!(CHEM_OK | Q_OK)) %>%
  select(domain, site_code, grp_start, types, chem_pre_bins, chem_post_bins,
         chem_pre_minmo, chem_post_minmo, q_pre_bins, q_post_bins)), width = 220)

# ---- sensitivity: drop the 7 pre-1900 (historic agriculture) rows ----
g2 <- build_groups(TRUE)
key <- function(d) paste(d$domain, d$site_code, d$grp_start, d$grp_end)
a <- grp %>% filter(dist_clear); b <- g2 %>% filter(dist_clear)
cat("\n=== SENSITIVITY: dropping pre-1900 rows ===\n")
cat("disturbance-clear groups, with pre-1900:", nrow(a), " without:", nrow(b), "\n")
cat("groups that appear ONLY when pre-1900 rows dropped:\n")
print(as.data.frame(b %>% filter(!key(b) %in% key(a)) %>% select(domain, site_code, grp_start, grp_end, types)))
cat("groups that disappear when pre-1900 rows dropped:\n")
print(as.data.frame(a %>% filter(!key(a) %in% key(b)) %>% select(domain, site_code, grp_start, grp_end, types)))

saveRDS(res, "verified.rds")
