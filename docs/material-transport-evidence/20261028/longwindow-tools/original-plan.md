# Latest-source long-window isotope diagnosis

Frozen source: `c315c84973d0672d120b2b6926343a1b42b6ce21`.
Owner scope: this private directory only; no production edits or parameter tuning.

1. Reproduce all three original p9k windows with explicit NEC/AMZ 10x10 domains,
   actual IsoGSM forcing, fractionation enabled, and strict conservation gates.
2. Locate the current extrema by variable, date and cell. Do not assume old p9k
   values survive precipitation-fallback changes.
3. Hypothesis H1 (not yet a finding): strong evaporation tails are signed
   heavy-isotope exchange divided by positive outward water under high humidity,
   rather than extreme composition of a stored water pool. Test actual component
   gross-water/signed-mass ledgers and exact finite-pool inputs, never net host ET.
4. Privately instrument selected cells/days, replay the finite-pool formula at
   high precision, and compare complete history variable bits with the baseline
   to prove instrumentation did not alter the reproduced result.

Stop at evidence-backed numerical diagnosis. Do not claim observational accuracy,
global physical calibration, or all-region conclusions from selected samples.

## Primary-reference boundary

- Benettin et al. (2018), HESS 22, 2881, Eqs. 1 and 5–7:
  https://hess.copernicus.org/articles/22/2881/2018/
  The source provides the humidity-dependent Craig–Gordon flux and finite-volume
  residual-water solution. These equations support the model form, not this
  case's meteorological inputs or empirical accuracy.
- Dar et al. (2020), ACP 20, 11435, Eqs. 8–11:
  https://acp.copernicus.org/articles/20/11435/2020/
  The source separates escaping and captured vapor in the net flux budget.
  Inference from that budget: a negative heavy-isotope transfer can coexist with
  positive total-water loss if the heavy-isotope gradient points inward. This
  is a signed flux ratio, not a positive outgoing-parcel composition. Whether
  that condition occurs in each model tail still requires measured run inputs.
