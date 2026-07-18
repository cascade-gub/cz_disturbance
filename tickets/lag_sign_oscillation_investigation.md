# TICKET — investigate the ± apparent-lag oscillation at the 365 d band

**Type:** investigation / diagnosis (not a fix). Produce findings + a recommendation, not just code.
**Notebooks:** primarily `ms_data/analysis_04_wavelet_synchrony.Rmd` (P→Q); check whether the same
behavior appears in `ms_data/analysis_05` (C→Q) and `syn_data/analysis_02` (synthetic).
**Do not touch** the tracked series CSVs.

**Status:** deferred / not started.

---

## The phenomenon

In the band-metrics-through-time figures (the `*-trace` chunks, `plot_band_trace`), the
**apparent-lag panel at `P≈365`** oscillates between large positive and large negative values —
pinning near `±P/2` (`±182.5 d`) and flipping near-vertically, each jump ≈ one period.

## What is ALREADY known / done (read first — don't rediscover)

`tickets/analysis_04_lag_fix.txt` diagnosed and partially addressed this:

- The flips are **2π phase wraps, not real lags.** `band_trace` takes a circular mean of `wc$Angle`
  at the P band → result on `(−π, π]` → mapped linearly to days → trapped in `(−P/2, P/2]`. Where
  the true phase sits near `±π` (antiphase), noise flips it across the branch cut, so the plotted
  "lag" jumps `+182.5 ↔ −182.5` — the **same physical relationship** rendered as a sawtooth.
- **Root cause:** precip is nearly aseasonal (annual-harmonic R² ≈ 0.001), so at `P≈365` the
  **driver has essentially no amplitude**, its phase is undefined, and the phase *difference*
  wanders across `±π`.
- **Mitigation already in the notebook:** `band_trace` gates φ/lag to `NA` where driver band
  amplitude `ax < PHASE_MIN_AMP × median(ax[coi])` (`PHASE_MIN_AMP = 0.5`); `wave_summary`
  suppresses settled φ/lag where block-phase concentration `R < PHASE_MIN_R = 0.5`.

**So the mechanism is understood and the rendering is gated.** This ticket is the deeper pass the
lag-fix did not do: *characterize what remains, decide if the gating is right, and rule out real
structure.*

## Open questions to answer

1. **Does gating fully kill the oscillation, or do gated-IN regions still flip?** Plot the raw
   (un-gated) `ang → lag` alongside the gated version per site. Where `phase_ok == TRUE`, does the
   lag still change sign across the branch cut? If yes, the amplitude gate alone is insufficient.
2. **Is any surviving flip signal or artifact?** Distinguish three cases at each flip:
   (a) branch-cut wrap of a stable near-`±π` phase (artifact), (b) genuine drift of a well-defined
   phase through `±π` (real, but needs unwrapping to read), (c) noise where φ is undefined (should
   be gated). Use the per-time driver amplitude `ax` and a per-time phase-concentration measure to
   classify.
3. **Is `±P/2` here physically antiphase?** The lag-fix notes log-discharge's annual peak sits
   ~177 d (≈ ½ yr) from precip's. Confirm per site: is the surviving phase genuinely ~π (real
   antiphase between an aseasonal driver and a seasonal response), and should the panel therefore
   be *labeled* "antiphase (±½ yr)" rather than gated away?
4. **Threshold sensitivity.** Sweep `PHASE_MIN_AMP` (e.g. 0.25 / 0.5 / 0.75) and `PHASE_MIN_R`;
   show how much of the trace survives and whether the qualitative read is stable. Is `0.5`
   defensible or arbitrary?
5. **Unwrapping.** The lag-fix flagged optional unwrapping *within contiguous gated-in runs, never
   across gaps.* Test whether unwrapping makes any site's surviving phase interpretable, or whether
   it just manufactures apparent drift. Recommend for/against, per site.
6. **Cross-notebook:** does the same ± oscillation occur in `05` (C→Q — driver is `log Q`, which
   DOES have an annual cycle, so maybe not) and in `02` (clean synthetic sine — should be a stable,
   single-signed lag)? The contrast is itself the finding: the oscillation is a signature of an
   *aseasonal driver at the annual band*, and `02`/`05` are the controls that show it appearing
   only where the driver lacks band power.

## Suggested diagnostics (all deterministic, no RNG)

- Per site, overlay: raw lag, gated lag, driver band amplitude `ax` (with the `PHASE_MIN_AMP·median`
  threshold line), and a rolling phase-concentration `R` over ~1-cycle windows. One figure that
  shows "the flips happen exactly where `ax` collapses" makes the artifact self-evident.
- A small table: per site, fraction of COI-valid time that is phase-gated-out, number of sign
  flips before vs after gating, and the circular concentration `R` of the surviving phase.
- Tie back to `annual_r2` (already computed): flips ∝ how aseasonal the driver is.

## Deliverable

A short written finding (could live as a new section in `analysis_04`, or a memo in this ticket)
that states: the oscillation is/or-isn't fully explained by the aseasonal-driver branch-cut wrap;
whether the current gating is sufficient and the threshold defensible; whether surviving `±P/2`
should be **labeled as antiphase** rather than blanked; and a yes/no on unwrapping. If it prompts a
code change, that change should land in the **shared wavelet kernel**
(`tickets/wavelet_significance_mask.md`, Part 1) so all three notebooks get it once.

## Constraints

- Deterministic; reproducible from the tracked series CSVs alone.
- Don't invent coupling at the annual band — the goal is honest characterization, consistent with
  the notebook's thesis (the real P→Q coupling lives at the event band).
- Minimal, reversible diffs if the investigation turns into a notebook change.
