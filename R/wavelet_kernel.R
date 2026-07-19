# R/wavelet_kernel.R — shared wavelet helpers for the three synchrony notebooks
# (syn_data/analysis_02, ms_data/analysis_04, ms_data/analysis_05). EDIT ONCE.
#
# TWO parts:
#   (1) the red-noise cross-power SIGNIFICANCE helpers (ar1_spectrum / xpower_sig_level /
#       xpower_sig_ref_tavg / add_xpower_contour) — pure, identical across notebooks; and
#   (2) make_wavelet_kernel(cfg) — the factory that builds the WHOLE helper family (wco, band_trace,
#       block_stat, wave_summary, plot_*, fit_pre, resid_trace, ...) once, parameterized by the
#       handful of things that actually differ between the three notebooks. Each notebook does
#       `list2env(make_wavelet_kernel(cfg), environment())` in its helpers chunk, so its call sites
#       stay bare (wco(df), band_trace(wc, df), ...). This is the anti-drift mechanism: a fix or a
#       new feature lands in one place, not three. See tickets/wavelet_significance_mask.plan.md.
#
# WHY cross-power and not coherence: WaveletComp's normalized, smoothed coherence saturates near 1
# across the whole record (Cauchy–Schwarz resting state on red series), so the honest signal lives
# in un-normalized cross-wavelet power. Cross-power has an ANALYTIC AR(1) red-noise significance
# level (Torrence & Compo 1998) — no Monte Carlo, deterministic, reproducible at knit.

# ---- AR(1) red-noise theoretical spectrum (Torrence & Compo 1998, eq. 16) --------------------
# Unit-variance normalized discrete red-noise power at normalized frequency f = 1/Period (dt = 1,
# i.e. per wavelet step). alpha = lag-1 autocorrelation of the (standardization-invariant) channel.
ar1_spectrum <- function(alpha, f) (1 - alpha^2) / (1 + alpha^2 - 2 * alpha * cos(2 * pi * f))

# ---- lag-1 autocorrelations of the two analyzed channels ------------------------------------
# Read straight off the RAW channels: alpha is standardization-invariant, and loess.span = 0 makes
# WaveletComp standardize each series to unit variance, so sigma_x = sigma_y = 1 in the level below.
# (05: climate = log Q, log_response = log C — the code is identical, only the meaning differs.)
.chan_alpha <- function(df) c(
  x = acf(df$climate,      lag.max = 1, plot = FALSE)$acf[2],
  y = acf(df$log_response, lag.max = 1, plot = FALSE)$acf[2])

# ---- POINTWISE (instantaneous) cross-power red-noise 95% level -------------------------------
# One threshold per wc$Period row, for the INSTANTANEOUS cross-power field wc$Power.xy (nu = 2).
#
# Two WaveletComp facts, both verified against the installed source + a seeded AR(1) calibration
# (scratch: calib.R), NOT assumed:
#   (1) wc$Power.xy = Mod(Wave.xy) — a MODULUS, not a modulus^2 — so we compare it to `level`,
#       not `level^2`.
#   (2) Wave.xy = Wave.x * Conj(Wave.y) / Scale — WaveletComp rectifies the cross term by /Scale.
#       The T&C level therefore carries a matching / wc$Scale, or the contour's period-SHAPE is
#       wrong. Calibration confirmed the /Scale form (ratio emp_q95 / level ~ 1 across all periods;
#       the no-Scale form drifted by ~100x from short to long period).
# Level (T&C 1998, cross-wavelet power): |W_xy|/(sd_x sd_y) 95% = (Z_nu(p)/nu) sqrt(P_x P_y),
# nu = 2, Z_2(0.95) = 3.999; divided by wc$Scale for WaveletComp's rectified Power.xy.
xpower_sig_level <- function(df, wc, conf = 0.95) {
  a  <- .chan_alpha(df)
  f  <- 1 / wc$Period                                  # per-step frequency (05: Period already in steps)
  Z2 <- if (isTRUE(all.equal(conf, 0.95))) 3.999 else stop("tabulate Z_2 for conf != 0.95")
  (Z2 / 2) * sqrt(ar1_spectrum(a["x"], f) * ar1_spectrum(a["y"], f)) / wc$Scale
}

# ---- TIME-AVERAGED cross-power red-noise 95% reference (for the xpower-vs-period figure) ------
# The "where does shared power live" figure plots rowMeans(wc$Power.xy) — a TIME-AVERAGE over n_a
# columns. A pointwise (nu = 2) level is far too tall there: time-averaging shrinks the observable's
# variance, so the DOF must be inflated (Torrence & Compo 1998 eq. 23) and the reference lowered.
#   E[Power.xy | red-noise H0] = (pi/4) sqrt(P_x P_y) / Scale     (product of two Rayleigh moduli)
#   nu(s) = 2 sqrt( 1 + ( n_a dt / (gamma s) )^2 ),  gamma = 2.32 (Morlet),  s = wc$Scale
#   level_avg(f) = E[...] * qchisq(conf, nu) / nu
# Calibrated against the surrogate distribution of rowMeans(Power.xy), seeded (scratch: calib2.R),
# NOT assumed. Outcome: the (pi/4) mean-constant hits dead on (empirical Ebar/root = 0.793 vs
# pi/4 = 0.785), and the DOF reference matches emp. 95% to ~10% for periods up to ~3 months (the
# event bands, where the significance CALL carries the "coupling lives at short periods" argument).
# At the annual band and beyond the reference runs mildly CONSERVATIVE (sits ~20-50% above the
# surrogate 95%): a ~365 d wavelet fits only ~5 independent cycles in these records, so both the
# empirical target and the low-DOF (nu~5) formula are least certain there — and a conservative
# annual reference is the safe direction for judging whether an annual peak is genuinely present.
xpower_sig_ref_tavg <- function(df, wc, n_a, conf = 0.95, gamma = 2.32, dt = 1) {
  a    <- .chan_alpha(df)
  f    <- 1 / wc$Period
  Emean <- (pi / 4) * sqrt(ar1_spectrum(a["x"], f) * ar1_spectrum(a["y"], f)) / wc$Scale
  nu   <- 2 * sqrt(1 + (n_a * dt / (gamma * wc$Scale))^2)
  Emean * qchisq(conf, nu) / nu
}

# ---- persistence (area) filter for the pointwise significance field --------------------------
# The raw pointwise (nu = 2) test speckles badly on long high-resolution records: below ~a month
# it flags thousands of single-event crossings, an unreadable black texture. A genuine oscillation
# at scale s persists >~ one cycle in time (the Morlet footprint is ~s wide), whereas those isolated
# crossings are sub-footprint. So keep a significant point only where cross-power stays significant
# over a sustained run of >= `persist` cycles at that scale — a 1D-in-time duration reduction of the
# pointwise field (the cheap, defensible cousin of areawise significance; cf. Maraun & Kurths 2004,
# Schulte 2016). Removes the speckle, keeps sustained bands (the annual coupling, real event bands).
.persist_filter <- function(sig, period, persist = 1) {
  for (i in seq_len(nrow(sig))) {
    minlen <- max(1L, round(persist * period[i]))        # required run length, in time steps
    if (minlen <= 1L) next
    r <- rle(sig[i, ]); drop <- r$values & r$lengths < minlen
    if (any(drop)) { r$values[drop] <- FALSE; sig[i, ] <- inverse.rle(r) }
  }
  sig
}

# ---- overlay the significance contour on a live wc.image() -----------------------------------
# Call AFTER wc.image(..., graphics.reset = FALSE), same trick the notebooks use for guide lines.
# Outlines { instantaneous cross-power > red-noise 95%, sustained >= `persist` cycles } on the
# (time x log2 period) grid wc.image set up (wc$axis.1 = time index, wc$axis.2 = log2 Period). This
# marks where SHARED POWER is significant vs red noise — read against the saturated coherence colour
# field beneath it (coherence ~1 everywhere; significant cross-power is selective — the thesis).
add_xpower_contour <- function(df, wc, col = "black", lwd = 1, persist = 3) {
  op <- options(max.contour.segments = 1e5); on.exit(options(op))   # binary field → many short segments
  lev <- xpower_sig_level(df, wc)                        # length = nrow(Power.xy), one per Period row
  sig <- wc$Power.xy > lev                               # matrix > per-row vector recycles down columns → row i vs lev[i]
  sig <- .persist_filter(sig, wc$Period, persist)        # keep only coupling sustained >= persist cycles
  contour(wc$axis.1, wc$axis.2, t(sig) + 0, levels = 0.5, drawlabels = FALSE,
          add = TRUE, col = col, lwd = lwd)
}

# ==============================================================================================
# PART 2 — make_wavelet_kernel(cfg): the full shared helper family, config-parameterized.
# ==============================================================================================
# Config fields (defaults in parens) — every cross-notebook difference is exactly one of these:
#   P_days            band period in DAYS (02: 500, 04/05: 365)
#   phase_gate  (F)   blank phi/lag where the driver has no band amplitude (04 only)
#   PHASE_MIN_AMP(.5) trace gate: driver band-amp >= this x median(driver amp) (used iff phase_gate)
#   PHASE_MIN_R (.5)  summary gate: report phi/lag only where block-phase concentration R >= this
#   has_transition(F) settled `post` starts after a `transition` span (02) vs after `disturbance`
#   season      (F)   alternating calendar-quarter shading on the trace/resid plots (04/05)
#   trace_x  ("time") x column for trace/resid/bands — "time" (02/04) or "date" (05)
#   lag_wrap_break(F) break the lag line at +/-pi branch cuts and draw antiphase rails (04)
#   lag_lab ("app. lag") label of the lag facet ("app. lag (d)" for 04/05)
#   smooth_w   (51)   resid derivative smoothing window; number, or function(step) (05: 90/step)
#   lowerPeriod / upperPeriod   wco scale bounds; each a number OR function(P_wavelet, df)
#   wc_marks   (NULL) named period-guide vector for plot_wc (04/05); NULL -> bare image (02)
#   spec_axis   (F)   relabel the period axis to days via spec.period.axis (05, step != 1) vs manual
#   dist_vline ("grey70") colour of the dashed disturbance line
#   include_phase_R(F) add the phase_R column to emit_table (04)
#   table_note (NULL) extra text appended to every emit_table caption (04's phase_R explainer)
#   resp_lab / drv_lab  channel names in resid labels ("response"/"climate"; 05: "C"/"log(Q)")
#   resid_sub  a function(gain, L, step) -> the plot_resid subtitle string
#   gain_fmt  ("%.2f") sprintf format for the pre-fit gain in the resid subtitle
#   timelab / periodlab   base-graphics axis labels for plot_wc
# Channel semantics: all three analyze columns `climate` / `log_response`; only the MEANING differs
# (05: climate = log Q, log_response = log C). The code is identical — no config needed.

make_wavelet_kernel <- function(cfg) {
  `%||%`   <- function(a, b) if (is.null(a)) b else a
  P_days   <- cfg$P_days
  gate     <- isTRUE(cfg$phase_gate)
  A_MIN    <- cfg$PHASE_MIN_AMP %||% 0.5
  R_MIN    <- cfg$PHASE_MIN_R   %||% 0.5
  has_tr   <- isTRUE(cfg$has_transition)
  season_on<- isTRUE(cfg$season)
  trace_x  <- cfg$trace_x %||% "time"
  wrapbrk  <- isTRUE(cfg$lag_wrap_break)
  lag_lab  <- cfg$lag_lab %||% "app. lag"
  smooth_w <- cfg$smooth_w %||% 51
  lowerP   <- cfg$lowerPeriod
  upperP   <- cfg$upperPeriod
  wc_marks <- cfg$wc_marks
  spec_ax  <- isTRUE(cfg$spec_axis)
  dvcol    <- cfg$dist_vline %||% "grey70"
  incl_R   <- isTRUE(cfg$include_phase_R)
  tbl_note <- cfg$table_note
  resp_lab <- cfg$resp_lab %||% "response"
  drv_lab  <- cfg$drv_lab  %||% "climate"
  gain_fmt <- cfg$gain_fmt %||% "%.2f"
  timelab  <- cfg$timelab %||% "time"
  periodlab<- cfg$periodlab %||% "period"
  resid_sub<- cfg$resid_sub %||% function(g, L, step)
    sprintf(paste0("pre-fit predictor:  log(", resp_lab, ")(t) = ", gain_fmt,
                   " . ", drv_lab, "(t - %d)"), g, L)
  # resid axis labels live in the PLOT (not any table) — supply exact strings per notebook so the
  # extraction is output-preserving (02/04 use a unicode minus + two spaces; 05 uses one space).
  resid_ylab <- cfg$resid_ylab %||% paste0("residual = log(", resp_lab, ") − predicted  (green)")
  resid_sec  <- cfg$resid_sec  %||% "d(residual)/dt  (smoothed, orange)"
  # trace/resid cosmetics that drifted between notebooks (02 drew slightly thicker lines and a plain
  # "time" x-label; 04 used "time (days)"; 05 is on a date axis, no x-label). Kept as config so the
  # extraction is output-preserving down to line width and axis-label wording.
  x_lab      <- if (trace_x == "date") NULL else (cfg$x_lab %||% "time")
  trace_lw   <- cfg$trace_lw  %||% 0.3
  resid_lw   <- cfg$resid_lw  %||% 0.3
  dresid_lw  <- cfg$dresid_lw %||% 0.25
  metric_labs <- c(coherence = "coherence", xpower = "cross power",
                   amp_ratio = "amp ratio", lag = lag_lab)

  resolve <- function(x, P, df) if (is.function(x)) x(P, df) else x
  sm_win  <- function(step) { w <- if (is.function(smooth_w)) smooth_w(step) else smooth_w; max(3, round(w)) }

  # ---- small numerics -----------------------------------------------------------------------
  CMEAN     <- function(a) Arg(mean(exp(1i * a)))
  smooth_ma <- function(x, step = 1) { w <- sm_win(step); as.numeric(stats::filter(x, rep(1 / w, w))) }
  band_rows <- function(wc, step = 1) which(abs(wc$Period * step - P_days) <= 0.05 * P_days)

  # ---- WaveletComp wrapper (climate vs log_response; loess.span = 0) -------------------------
  wco <- function(df, step = 1) {
    P <- P_days / step
    analyze.coherency(as.data.frame(df[, c("climate", "log_response")]),
                      my.pair = c("climate", "log_response"), loess.span = 0, dt = 1, dj = 1/20,
                      lowerPeriod = resolve(lowerP, P, df), upperPeriod = resolve(upperP, P, df),
                      make.pval = FALSE, verbose = FALSE)
  }

  # ---- collapse the wc object to one time-resolved series at the P band ----------------------
  band_trace <- function(wc, df, step = 1) {
    bi   <- band_rows(wc, step); EF <- round(sqrt(2) * P_days / step)
    sd_c <- sd(df$climate); sd_r <- sd(df$log_response)
    ang  <- apply(wc$Angle[bi, , drop = FALSE], 2, CMEAN)
    coi  <- df$time >= EF & df$time <= nrow(df) - EF
    ax   <- colMeans(wc$Ampl.x[bi, , drop = FALSE])
    ok   <- if (gate) ax >= A_MIN * median(ax[coi]) else rep(TRUE, length(ax))
    tibble(time       = df$time,
           date       = if (is.null(df$date)) as.Date(NA) else df$date,
           phase      = df$phase,
           coherence  = colMeans(wc$Coherence[bi, , drop = FALSE]),
           xpower     = colMeans(wc$Power.xy[bi, , drop = FALSE]),
           amp_ratio  = colMeans(wc$Ampl.y[bi, , drop = FALSE]) /
                        colMeans(wc$Ampl.x[bi, , drop = FALSE]) * (sd_r / sd_c),
           driver_amp = ax,
           angle      = ifelse(ok, ang, NA_real_),
           lag        = ifelse(ok, ang / (2 * pi) * P_days, NA_real_),
           phase_ok   = ok,
           coi_ok     = coi)
  }

  # ---- settled pre/post summary with block SEs ----------------------------------------------
  # One estimate per P-long cycle in a COI-buffered window; NA (gated) steps dropped per block.
  # Returns R (phase concentration) too — used iff phase_gate for the summary phi/lag gate.
  block_stat <- function(t, v, step = 1, circ = FALSE) {
    blk <- max(2, round(P_days / step))
    b   <- (t - min(t)) %/% blk
    bv  <- vapply(split(v, b), function(vi) {
             vi <- vi[!is.na(vi)]
             if (!length(vi)) NA_real_ else if (circ) CMEAN(vi) else mean(vi)
           }, numeric(1))
    bv  <- bv[!is.na(bv)]; n <- length(bv)
    if (circ) {
      R    <- if (n) Mod(mean(exp(1i * bv))) else NA_real_
      disp <- if (n && R < 1) sqrt(-2 * log(R)) else 0
    } else { R <- NA_real_; disp <- if (n > 1) sd(bv) else NA_real_ }
    list(value = if (!n) NA_real_ else if (circ) CMEAN(bv) else mean(bv),
         se = if (n > 1) disp / sqrt(n) else NA_real_, n = n, R = R)
  }

  wave_summary <- function(tr, step = 1) {
    EF <- round(sqrt(2) * P_days / step)
    ds <- min(tr$time[tr$phase == "disturbance"])
    de <- if (has_tr) max(tr$time[tr$phase == "transition"]) else max(tr$time[tr$phase == "disturbance"])
    N  <- max(tr$time)
    wins <- list(pre  = tr$time >= EF      & tr$time <= ds - EF,
                 post = tr$time >= de + EF & tr$time <= N  - EF)
    one <- function(sel) {
      t  <- tr$time[sel]
      co <- block_stat(t, tr$coherence[sel], step); xp <- block_stat(t, tr$xpower[sel], step)
      ar <- block_stat(t, tr$amp_ratio[sel], step); an <- block_stat(t, tr$angle[sel], step, circ = TRUE)
      if (gate) {                                   # report phi/lag only where block phases concentrate
        ph_ok <- !is.na(an$R) && an$n >= 2 && an$R >= R_MIN
        phi <- if (ph_ok) an$value else NA_real_;  phise <- if (ph_ok) an$se else NA_real_
      } else { phi <- an$value; phise <- an$se }
      tibble(n_blocks = co$n,
             coherence = co$value, coherence_se = co$se,
             xpower    = xp$value, xpower_se    = xp$se,
             amp_ratio = ar$value, amp_ratio_se = ar$se,
             phi       = phi,                      phi_se = phise,
             lag       = phi / (2*pi) * P_days,    lag_se = phise / (2*pi) * P_days,
             phase_R   = an$R, phase_n = an$n)
    }
    bind_rows(pre = one(wins$pre), post = one(wins$post), .id = "phase") %>%
      mutate(phase = factor(phase, c("pre", "post")))
  }

  # ---- model-based residual (whole series, log space) ---------------------------------------
  fit_pre <- function(df, step = 1) {
    pre <- df[df$phase == "pre", ]; cl <- pre$climate; rs <- pre$log_response; n <- length(rs)
    Ls  <- 0:max(1, round(0.6 * P_days / step))
    cr  <- vapply(Ls, function(L) suppressWarnings(cor(rs[(L + 1):n], cl[seq_len(n - L)])), numeric(1))
    L   <- Ls[which.max(replace(cr, is.na(cr), -Inf))]
    x   <- cl[seq_len(n - L)]; y <- rs[(L + 1):n]
    list(lag = L, gain = sum(x * y) / sum(x * x))
  }
  resid_trace <- function(df, step = 1) {
    fp <- fit_pre(df, step); L <- fp$lag; N <- nrow(df)
    pred <- rep(NA_real_, N); if (L < N) pred[(L + 1):N] <- fp$gain * df$climate[seq_len(N - L)]
    resid <- df$log_response - pred
    out <- tibble(time = df$time,
                  date = if (is.null(df$date)) as.Date(NA) else df$date,
                  phase = df$phase, resid = resid,
                  dresid = c(NA, diff(smooth_ma(resid, step))))
    attr(out, "lag") <- L; attr(out, "gain") <- fp$gain; attr(out, "step") <- step; out
  }

  # ---- shared plotting furniture ------------------------------------------------------------
  phase_fill <- c(pre = "#9ecae1", disturbance = "grey45", transition = "#fcbba1", post = "#fdae6b")
  bar_fill   <- phase_fill[c("pre", "post")]
  phase_bands <- function(ph) ph %>% filter(phase != "disturbance") %>%
    group_by(phase) %>% summarise(xmin = min(.data[[trace_x]]), xmax = max(.data[[trace_x]]), .groups = "drop")
  phase_band_layers <- function(d) {
    bands <- phase_bands(d); td <- d[[trace_x]][d$phase == "disturbance"]
    list(geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase),
                   alpha = 0.30, inherit.aes = FALSE),
         geom_vline(xintercept = td, linetype = 2, colour = dvcol),
         scale_fill_manual(values = phase_fill, guide = "none"))
  }
  # alternating calendar-quarter shading (needs a `date` column); x placed on the trace_x axis.
  season_layer <- function(d) {
    if (!season_on) return(NULL)
    q  <- (as.integer(format(d$date, "%m")) - 1) %/% 3 + 1
    yr <- as.integer(format(d$date, "%Y"))
    r  <- tibble(x = d[[trace_x]], qi = yr * 4 + q, q = q) %>%
      group_by(qi) %>% summarise(xmin = min(x), xmax = max(x), odd = first(q) %% 2 == 1, .groups = "drop") %>%
      filter(odd)
    geom_rect(data = r, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey20", alpha = 0.14, inherit.aes = FALSE)
  }

  # ---- coherence heatmap (base graphics) + significance contour ------------------------------
  plot_wc <- function(wc, df, main, step = 1) {
    op <- par(bg = "white", fg = "black", col.axis = "black", col.lab = "black", col.main = "black")
    on.exit(par(op))
    at <- if (is.null(wc_marks)) NULL else wc_marks / step
    keep <- if (is.null(at)) NULL else which(at >= min(wc$Period) & at <= max(wc$Period))
    img_extra <- if (spec_ax && !is.null(at))
      list(spec.period.axis = list(at = at[keep], labels = names(wc_marks)[keep])) else list()
    do.call(wc.image, c(list(wc, which.image = "wc", plot.arrow = TRUE, plot.contour = FALSE,
      color.key = "interval", timelab = timelab, periodlab = periodlab,
      legend.params = list(lab = "wavelet coherence"), main = main, graphics.reset = FALSE), img_extra))
    if (!is.null(at)) {
      abline(h = log2(at[keep]), lty = 2, col = "black", lwd = 1.4)
      if (!spec_ax) {                              # 04: manual left-edge labels (default period axis kept)
        usr <- par("usr")
        text(x = usr[1] + 0.01 * (usr[2] - usr[1]), y = log2(at[keep]), labels = names(wc_marks)[keep],
             adj = c(0, -0.4), font = 2, cex = 0.9, col = "black")
      }
    }
    add_xpower_contour(df, wc)                     # red-noise 95% cross-power significance (Part 1)
    # Mark WHEN the disturbance occurred, on the time axis (wc$axis.1 = 1:nc = df rows). Use the
    # `disturbance`-phase column(s): one bold line for a single step, the onset/end pair for a span
    # (fernow_WS-5). The binned C-Q grids (05) carry NO disturbance bin for a single-instant event —
    # it is the pre→post boundary — so fall back to the midpoint of (last pre, first post). Drawn
    # last so it reads over the field, arrows and contour.
    td <- which(df$phase == "disturbance")
    if (!length(td)) {
      pre <- which(df$phase == "pre"); post <- which(df$phase == "post")
      if (length(pre) && length(post)) td <- (max(pre) + min(post)) / 2
    }
    if (length(td)) {
      abline(v = range(td), lty = 1, col = "black", lwd = 2)
      text(x = mean(range(td)), y = par("usr")[4], labels = "disturbance",
           adj = c(0.5, 1.3), font = 2, cex = 0.9, col = "black")
    }
  }

  # ---- time-resolved band metrics + settled bars + table -------------------------------------
  seg_id <- function(v, jump = Inf) { d <- abs(c(NA, diff(v))); cumsum(is.na(v) | (!is.na(d) & d > jump)) }
  plot_band_trace <- function(tr, title) {
    m <- names(metric_labs)
    wide <- tr %>% mutate(across(all_of(m), ~ ifelse(coi_ok, ., NA_real_)))
    long <- bind_rows(lapply(m, function(mm)
      tibble(x = wide[[trace_x]], metric = mm, value = wide[[mm]],
             seg = seg_id(wide[[mm]], if (wrapbrk && mm == "lag") P_days/2 else Inf)))) %>%
      mutate(metric = factor(metric, m))
    g <- ggplot() + phase_band_layers(tr) + season_layer(tr)
    if (wrapbrk) {                                 # antiphase rails + zero line + note on the lag facet
      lf <- factor("lag", m)
      g <- g +
        geom_hline(data = tibble(metric = lf, y = c(-P_days/2, P_days/2)), aes(yintercept = y),
                   linetype = 2, colour = "grey55", linewidth = 0.3) +
        geom_hline(data = tibble(metric = lf, y = 0), aes(yintercept = y),
                   linetype = 1, colour = "grey55", linewidth = 0.3) +
        geom_text(data = tibble(metric = lf, x = max(tr[[trace_x]]), y = P_days/2, lab = "antiphase (±½ yr)"),
                  aes(x, y, label = lab), hjust = 1, vjust = -0.4, size = 2.6, colour = "grey45")
    }
    g + geom_line(data = long, aes(x, value, group = interaction(metric, seg)), linewidth = trace_lw, na.rm = TRUE) +
      facet_wrap(~ metric, scales = "free_y", ncol = 1, strip.position = "left",
                 labeller = as_labeller(metric_labs)) +
      labs(title = title, x = x_lab, y = NULL)
  }

  plot_resid <- function(df, title, step = 1) {
    rt <- resid_trace(df, step); L <- attr(rt, "lag"); g <- attr(rt, "gain")
    bands <- phase_bands(df); ytop <- max(abs(rt$resid), na.rm = TRUE)
    sc <- ytop / max(abs(rt$dresid), na.rm = TRUE)
    ggplot() + phase_band_layers(df) + season_layer(df) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
      geom_line(data = rt, aes(.data[[trace_x]], dresid * sc), colour = "#d94801", linewidth = dresid_lw, na.rm = TRUE) +
      geom_line(data = rt, aes(.data[[trace_x]], resid),       colour = "#238b45", linewidth = resid_lw, na.rm = TRUE) +
      geom_text(data = mutate(bands, x = xmin + (xmax - xmin) / 2),   # Date-safe midpoint (== (xmin+xmax)/2 for numeric x)
                aes(x = x, y = ytop, label = phase), vjust = 1, size = 3, colour = "grey25") +
      scale_y_continuous(sec.axis = sec_axis(~ . / sc, name = resid_sec)) +
      labs(title = title, x = x_lab, y = resid_ylab, subtitle = resid_sub(g, L, step))
  }

  metrics_long <- function(sm, keys) {
    m <- names(metric_labs)
    val <- sm %>% select(all_of(c(keys, m))) %>% pivot_longer(all_of(m), names_to = "metric", values_to = "value")
    se  <- sm %>% select(all_of(c(keys, paste0(m, "_se")))) %>%
      rename_with(~ sub("_se$", "", .x), all_of(paste0(m, "_se"))) %>%
      pivot_longer(all_of(m), names_to = "metric", values_to = "se")
    left_join(val, se, by = c(keys, "metric")) %>% mutate(metric = factor(metric, m))
  }
  plot_wave_bars <- function(sm, title) {
    ggplot(metrics_long(sm, "phase"), aes(phase, value, fill = phase)) +
      geom_col(width = 0.7) +
      geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.25, colour = "black", na.rm = TRUE) +
      facet_wrap(~ metric, scales = "free_y", nrow = 1, labeller = as_labeller(metric_labs)) +
      scale_fill_manual(values = bar_fill, guide = "none") + labs(title = title, x = NULL, y = NULL)
  }
  emit_table <- function(sm, caption) {
    cols <- c("phase", "n_blocks", "coherence", "xpower", "amp_ratio", "phi", "lag")
    if (incl_R) cols <- c(cols, "phase_R")
    cap  <- if (is.null(tbl_note)) caption else paste(caption, tbl_note)
    cat(kable(sm %>% select(all_of(cols)), digits = 2, caption = cap), sep = "\n")
  }

  # ---- "where does the shared power live" diagnostic ----------------------------------------
  # Time-averaged cross-power vs period, plus the red-noise 95% reference (Part 1). `period` is in
  # the plotted unit: days (Period x step). n_a = the number of time columns rowMeans averages over.
  #   split = FALSE → one curve over the whole record (id, period, power, ref).
  #   split = TRUE  → one curve per settled phase (adds a `phase` pre/post column), each averaged
  #     over that phase's EF-BUFFERED window — the SAME windows wave_summary uses for the tables, so
  #     the pre/post spectra are directly comparable to the settled pre/post numbers and are free of
  #     COI/disturbance bleed (a band wavelet near the event would otherwise mix pre and post). The
  #     reference uses each window's own n_a (so a shorter window ⇒ higher, more conservative level);
  #     alpha stays whole-record (the null red-noise process is a series property).
  xpower_spectrum <- function(wc, df, id, step = 1, split = FALSE) {
    per <- wc$Period * step
    if (!split)
      return(tibble(id = id, period = per, power = rowMeans(wc$Power.xy),
                    ref = xpower_sig_ref_tavg(df, wc, n_a = ncol(wc$Power.xy))))
    EF <- round(sqrt(2) * P_days / step); N <- nrow(df)
    ds <- min(df$time[df$phase == "disturbance"])
    de <- if (has_tr) max(df$time[df$phase == "transition"]) else max(df$time[df$phase == "disturbance"])
    wins <- list(pre  = df$time >= EF      & df$time <= ds - EF,
                 post = df$time >= de + EF & df$time <= N  - EF)
    one <- function(sel, ph) tibble(id = id, phase = ph, period = per,
      power = rowMeans(wc$Power.xy[, sel, drop = FALSE]),
      ref   = xpower_sig_ref_tavg(df, wc, n_a = sum(sel)))
    bind_rows(one(wins$pre, "pre"), one(wins$post, "post")) %>%
      mutate(phase = factor(phase, c("pre", "post")))
  }

  list(P_days = P_days, metric_labs = metric_labs, phase_fill = phase_fill, bar_fill = bar_fill,
       CMEAN = CMEAN, smooth_ma = smooth_ma, band_rows = band_rows, wco = wco,
       band_trace = band_trace, block_stat = block_stat, wave_summary = wave_summary,
       fit_pre = fit_pre, resid_trace = resid_trace, phase_bands = phase_bands,
       phase_band_layers = phase_band_layers, season_layer = season_layer, plot_wc = plot_wc,
       seg_id = seg_id, plot_band_trace = plot_band_trace, plot_resid = plot_resid,
       metrics_long = metrics_long, plot_wave_bars = plot_wave_bars, emit_table = emit_table,
       xpower_spectrum = xpower_spectrum,
       xpower_sig_level = xpower_sig_level, xpower_sig_ref_tavg = xpower_sig_ref_tavg,
       add_xpower_contour = add_xpower_contour, ar1_spectrum = ar1_spectrum)
}
