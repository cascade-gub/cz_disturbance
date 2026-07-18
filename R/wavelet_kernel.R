# R/wavelet_kernel.R — shared wavelet helpers for the three synchrony notebooks
# (syn_data/analysis_02, ms_data/analysis_04, ms_data/analysis_05).
#
# Phase 1 (this file today): the red-noise cross-power SIGNIFICANCE helpers only — new, pure, and
# identical across all three notebooks, so they live here and are edited once. The per-notebook
# wco / band_trace / plot_* helpers are NOT unified here yet (they have drifted structurally); that
# is Phase 2 of tickets/wavelet_significance_mask.plan.md.
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

# ---- overlay the pointwise significance contour on a live wc.image() -------------------------
# Call AFTER wc.image(..., graphics.reset = FALSE), same trick the notebooks use for guide lines.
# Draws the boundary of { instantaneous cross-power > red-noise 95% } on the (time x log2 period)
# grid wc.image set up (wc$axis.1 = time index, wc$axis.2 = log2 Period). This marks where SHARED
# POWER is significant vs red noise — read against the saturated coherence colour field beneath it.
add_xpower_contour <- function(df, wc, col = "black", lwd = 1.6) {
  lev <- xpower_sig_level(df, wc)                       # length = nrow(Power.xy), one per Period row
  sig <- wc$Power.xy > lev                              # matrix > per-row vector recycles down columns → row i vs lev[i]
  contour(wc$axis.1, wc$axis.2, t(sig), levels = 1, drawlabels = FALSE,
          add = TRUE, col = col, lwd = lwd)
}
