# TICKET — period-distribution plots for the other wavelet parameters

Extend `analysis_04`'s **Figure 1** (time-averaged cross-wavelet power vs period, per site) into a
**family** of the same view for the *other* band parameters, so the reader can see how each wavelet
quantity is distributed across scale — not just read it at the single `P≈365` band.

**Notebooks:** `ms_data/analysis_04` first; then `ms_data/analysis_05` and `syn_data/analysis_02`
via the shared kernel. **Do not touch** the tracked series CSVs.

**Status:** deferred / not started.

---

## What "Figure 1" is (the template)

The `xpower-spectrum` chunk (`analysis_04`, ~line 398): `colMeans(wc$Power.xy)` — time-averaged
cross-wavelet power — plotted **vs period on a log-x axis**, faceted per site, with the annual band
(365 d) marked and week/month/year guides. It answers "**where across scale does the shared power
live?**" and shows it concentrating at short (event) periods, not the annual band. Helper:
`xpower_spectrum()`.

## The ask: the same "spectrum vs period" view for the other parameters

Each of these is already computed on the `wc` object; Figure 1 just happens to plot only one.
Add the analogous time-averaged-vs-period plot for:

- **Coherence** — `colMeans(wc$Coherence)` vs period. *This one is the honesty payoff:* it should
  show coherence sitting ≈1 across **all** periods, making the "coherence saturates everywhere"
  point visual and scale-resolved rather than asserted (ties directly to the significance-mask
  ticket).
- **Amplitude ratio** — band-amp ratio `Ampl.y / Ampl.x` (× `sd_r/sd_c` to physical units, as
  `band_trace` does) vs period. Shows at which scales the response amplifies/attenuates the driver.
- **Driver band amplitude** — `colMeans(wc$Ampl.x)` vs period — directly visualizes the aseasonal
  driver's near-absence of annual power (the root cause in the lag-oscillation ticket).
- **Phase / apparent lag** — mean phase (circular) vs period. Handle with care: only meaningful
  where the driver has amplitude (reuse the gating logic); likely show as phase, not day-lag, and
  annotate rather than draw a naïve line through undefined-phase scales.

Keep the exact styling of Figure 1: log-x period axis, week/month/year breaks, `P` marked,
per-site facets, native-light theme + white matting, `fig_n()` caption numbering.

## ⚠️ Confirm intent before building (one genuine ambiguity)

"Distribution plots" has two readings — pick with the user:

- **(A) Spectrum-over-period** — value vs period, exactly like Figure 1, for each parameter. This
  is the literal "like figure 1" reading and what this ticket is written around. **Assumed default.**
- **(B) Statistical distributions** — histograms / densities of each band parameter's *values*
  (e.g. distribution of the per-time coherence at `P≈365`, or of block estimates), to show spread
  and skew behind the settled point ± SE. Different figure entirely.

They are not mutually exclusive; (A) is the closer match to "figure 1". Confirm which (or both)
before implementing.

## Where the code goes

Add one generic helper to the **shared wavelet kernel** (`tickets/wavelet_significance_mask.md`,
Part 1) — e.g. `param_spectrum(wc, which = c("coherence","xpower","amp_ratio","driver_amp"))` —
so all three wavelet notebooks get the whole family from one definition, and Figure 1 becomes just
`param_spectrum(..., "xpower")`. Do this **after** (or as part of) the kernel extraction so it isn't
written three times. If the kernel refactor hasn't happened yet, still write it as a single
parameterized helper, not one chunk per parameter.

## Notes / gotchas

- **Circular quantities** (phase) must be averaged with `CMEAN`, never linearly — and gated where
  the driver lacks amplitude (see the lag-oscillation ticket). Don't plot a day-lag spectrum
  through scales where φ is undefined.
- **`analysis_05` `step` scaling:** the period axis is `wc$Period × step_days` (days); reuse the
  notebook's existing day-rescale so all panels share one period-in-days axis. Coarse-cadence sites
  (hjandrews, 21 d) legitimately go blank below ~3 months — show, don't hide.
- **`free_y` vs shared y:** Figure 1 uses `scales = "free_y"` (shape, not magnitude, per site).
  Decide per parameter — coherence should probably share a fixed `[0,1]` y so the saturation reads
  across sites; power/amplitude may keep free_y. State the choice in the caption (honesty about
  cross-site comparability, per the standardization discussion).

## Constraints

- Deterministic; reproducible from the tracked series CSVs alone.
- Native-light figures on white matting; `fig_n()` caption numbering; knit → `../docs/` copy hook.
- Minimal diffs; prefer one shared helper over per-notebook duplication.

## Verify before done

- Re-knit affected notebooks; new spectra render per site with correct period axis and `P` mark.
- The coherence spectrum visibly shows ≈1 across all periods (the saturation, made visual).
- If built in the kernel: Figure 1 still reproduces exactly (it's now one call into the helper).
- HTML copied to `../docs/`.
