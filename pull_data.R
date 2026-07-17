# Pull daily discharge, daily precipitation, and stream chemistry (NO3 + spCond)
# for the 5 qualifying sites. Writes tidy per-site-per-product CSVs to ./data.
suppressMessages({library(macrosheds); library(dplyr); library(readr)})

root <- "C:/Users/nic/cz_disturbance/ms_data"
out  <- "C:/Users/nic/cz_disturbance/data"
dir.create(out, showWarnings = FALSE)

sites <- tibble::tribble(
  ~domain,          ~site_code, ~no3_var,     ~precip_prod,
  "hbef",           "w1",       "NO3_N",      "precipitation",
  "hjandrews",      "GSWS10",   "NO3_N",      "precipitation",
  "fernow",         "WS-3",     "NO3_N",      "CUSTOMprecipitation",
  "fernow",         "WS-5",     "NO3_N",      "CUSTOMprecipitation",
  "santa_barbara",  "MC06",     "NO3_NO2_N",  "precipitation")

# load a product for one domain, keep one site; return NULL if unavailable/empty
grab <- function(prod, dom, site) {
  x <- tryCatch(ms_load_product(root, prod = prod, domains = dom, warn = FALSE),
                error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  x <- x %>% filter(site_code == site)
  if (!nrow(x)) NULL else x
}

manifest <- list()
for (i in seq_len(nrow(sites))) {
  s <- sites[i, ]; tag <- paste0(s$domain, "_", s$site_code)
  cat("\n==", tag, "==\n")

  ## --- daily discharge (L/s) ---
  q <- grab("discharge", s$domain, s$site_code)
  if (!is.null(q)) {
    q <- q %>% transmute(domain = s$domain, site_code, date, var = "discharge",
                         val, unit = "L/s", ms_status, ms_interp) %>% arrange(date)
    f <- file.path(out, paste0(tag, "__discharge_daily.csv")); write_csv(q, f)
    cat("  discharge:", nrow(q), "days", format(min(q$date)), "->", format(max(q$date)),
        "| interp-filled:", sum(q$ms_interp == 1), "\n")
    manifest[[length(manifest)+1]] <- tibble(site = tag, product = "discharge_daily",
      n = nrow(q), first = min(q$date), last = max(q$date), file = basename(f))
  } else cat("  discharge: NONE\n")

  ## --- daily precipitation (mm) ---
  p <- grab(s$precip_prod, s$domain, s$site_code)
  if (!is.null(p)) {
    p <- p %>% transmute(domain = s$domain, site_code, date, var = "precipitation",
                         val, unit = "mm", ms_status, ms_interp) %>% arrange(date)
    f <- file.path(out, paste0(tag, "__precipitation_daily.csv")); write_csv(p, f)
    cat("  precip:", nrow(p), "days", format(min(p$date)), "->", format(max(p$date)),
        "(", s$precip_prod, ")\n")
    manifest[[length(manifest)+1]] <- tibble(site = tag, product = "precipitation_daily",
      n = nrow(p), first = min(p$date), last = max(p$date), file = basename(f))
  } else cat("  precip: NONE\n")

  ## --- stream chemistry: NO3 + specific conductance (grab samples) ---
  ch <- grab("stream_chemistry", s$domain, s$site_code)
  if (!is.null(ch)) {
    ch <- ch %>% filter(ms_interp == 0, !is.na(val),
                        var %in% c(s$no3_var, "spCond")) %>%
      transmute(domain = s$domain, site_code, date, var, val,
                ms_status) %>% arrange(date, var)
    f <- file.path(out, paste0(tag, "__chemistry_NO3_spCond.csv")); write_csv(ch, f)
    smry <- ch %>% count(var)
    cat("  chem:", paste(smry$var, smry$n, collapse = " | "), "\n")
    for (v in unique(ch$var)) {
      cv <- ch %>% filter(var == v)
      manifest[[length(manifest)+1]] <- tibble(site = tag, product = paste0("chem_", v),
        n = nrow(cv), first = min(cv$date), last = max(cv$date), file = basename(f))
    }
  } else cat("  chem: NONE\n")
}

man <- bind_rows(manifest)
write_csv(man, file.path(out, "MANIFEST.csv"))
cat("\n==== MANIFEST ====\n"); print(as.data.frame(man), width = 200)
cat("\nfiles in ./data:\n"); print(list.files(out))
