# The 10-yr rule is a MINIMUM. For each qualifying disturbance, extend the pre-
# and post- periods outward as far as they can go: keep adding one-year bins while
# (a) no other disturbance intrudes and (b) the bin still has data in >=11 months.
suppressMessages({library(macrosheds); library(dplyr); library(readxl); library(tidyr)})

YRS <- 10; GAP_DAYS <- 31; MIN_MONTHS <- 11
root <- "C:/Users/nic/cz_disturbance/ms_data"
doms <- c("fernow", "hbef", "hjandrews", "santa_barbara", "santee")

res <- readRDS("verified.rds")
qual <- res %>% filter(CHEM_OK | Q_OK)

load_obs <- function(prod) {
  bind_rows(lapply(doms, function(d) {
    ms_load_product(root, prod = prod, domains = d, warn = FALSE) %>%
      filter(ms_interp == 0, !is.na(val)) %>%
      mutate(domain = d, date = as.Date(date)) %>%
      distinct(domain, site_code, date)
  }))
}
chem <- load_obs("stream_chemistry"); q <- load_obs("discharge")

# months with >=1 observation in the one-year bin [a,b)
months_in <- function(d, a, b) {
  lo <- min(a, b); hi <- max(a, b)
  w <- d[d >= lo & d < hi]
  if (!length(w)) return(0L)
  edges <- seq(lo, by = "1 month", length.out = 13)
  length(unique(na.omit(cut(as.numeric(w), breaks = as.numeric(edges),
                            include.lowest = TRUE, labels = FALSE))))
}

# how many consecutive monthly-resolved years extend away from the anchor,
# stopping at `limit` (the neighbouring disturbance, or Inf if none)
shift_yrs <- function(d, n) if (n == 0) d else seq(d, by = paste(n, "years"), length.out = 2)[2]

extend <- function(dat, dom, site, anchor, direction, limit) {
  d <- dat$date[dat$domain == dom & dat$site_code == site]
  if (!length(d)) return(0L)
  sgn <- if (direction == "pre") -1L else 1L
  k <- 0L
  repeat {
    a <- shift_yrs(anchor, sgn * (k + 1L))   # outer edge of the next year-bin
    b <- shift_yrs(anchor, sgn * k)          # inner edge
    if (!is.na(limit) && ((direction == "pre" && a < limit) ||
                          (direction == "post" && a > limit))) break  # neighbouring disturbance
    if (months_in(d, a, b) < MIN_MONTHS) break
    k <- k + 1L
    if (k > 80L) break
  }
  k
}

out <- qual %>% rowwise() %>%
  mutate(
    lim_pre  = as.Date(ifelse(is.finite(pre_clear_yrs),  prev_end,   NA), origin = "1970-01-01"),
    lim_post = as.Date(ifelse(is.finite(post_clear_yrs), next_start, NA), origin = "1970-01-01"),
    chem_pre_yrs  = extend(chem, domain, site_code, grp_start, "pre",  lim_pre),
    chem_post_yrs = extend(chem, domain, site_code, grp_end,   "post", lim_post),
    q_pre_yrs     = extend(q,    domain, site_code, grp_start, "pre",  lim_pre),
    q_post_yrs    = extend(q,    domain, site_code, grp_end,   "post", lim_post)) %>%
  ungroup() %>%
  mutate(
    # usable window on the data type each site qualifies on
    use_pre  = ifelse(CHEM_OK, chem_pre_yrs,  q_pre_yrs),
    use_post = ifelse(CHEM_OK, chem_post_yrs, q_post_yrs),
    window_start = as.Date(mapply(function(s, k) as.character(shift_yrs(as.Date(s), -k)), grp_start, use_pre)),
    window_end   = as.Date(mapply(function(e, k) as.character(shift_yrs(as.Date(e),  k)), grp_end,   use_post)),
    qualifies_on = ifelse(CHEM_OK & Q_OK, "stream chemistry + discharge",
                   ifelse(CHEM_OK, "stream chemistry only", "discharge only")),
    open_ended_chronic = grepl("chronic", classes) & grp_end == grp_start) %>%
  arrange(desc(CHEM_OK), domain, site_code) %>%
  select(network, domain, site_code,
         disturbance_start = grp_start, disturbance_end = grp_end,
         disturbance = types, class = classes, qualifies_on,
         window_start, window_end, pre_yrs = use_pre, post_yrs = use_post,
         chem_pre_yrs, chem_post_yrs, q_pre_yrs, q_post_yrs, open_ended_chronic)

print(as.data.frame(out), width = 300)
stopifnot(all(out$pre_yrs >= YRS), all(out$post_yrs >= YRS))  # windows must still meet the minimum

f <- "qualifying_disturbances.csv"
ok <- tryCatch({write.csv(out, f, row.names = FALSE); TRUE}, error = function(e) FALSE)
if (!ok) { f <- "qualifying_disturbances_windows.csv"; write.csv(out, f, row.names = FALSE) }
cat("\nwrote", f, "\n")
