# DSSAT-CSM rice Methane reference notes for CoLM

Date: 2026-05-15
Reference snapshot: `DSSAT/dssat-csm-os` `develop` commit `40e7738` cloned to `/tmp/dssat-csm-os` for inspection.
Primary files inspected: `Soil/GHG/Methane.for`, `Soil/GHG/Methmod.for`, `Soil/CERES_OrganicMatter/SoilOrg.for`, `Soil/GHG/GHG_mod.for`.

## 1. Source and license check

The user-provided table lists DSSAT-CSM as GPL-3.0. The current open-source repository snapshot inspected here ships `license.txt` with a BSD-3-Clause license notice. Treat this as **BSD-3-Clause at the inspected commit**, not GPL-3.0. Still, no DSSAT source code was copied into CoLM; these notes extract model structure and parameterization ideas only.

Official DSSAT page says the most recent `dssat-csm-os` source is available from GitHub and describes DSSAT-CSM as a modular crop-system model. The same page lists DSSAT v4.8.5 as the current released version.

## 2. DSSAT methane module shape

DSSAT's methane path is a rice/crop-system module, not a global wetland/lake methane module:

| Layer | DSSAT implementation | Role |
|---|---|---|
| Interface | `MethaneDynamics` in `Soil/GHG/Methane.for` | Couples CERES-Rice / soil organic matter / flood-water state to the methane solver. |
| Transport solver | Arah steady-state routines in `Soil/GHG/Methmod.for` | Solves O2 and CH4 concentrations and partitions fluxes among diffusion, root transport, ebullition, leaching, storage. |
| Carbon source | `newCO2` from soil organic matter | Daily decomposed carbon is partitioned into aerobic CO2 versus methanogenic substrate. |
| Water-management control | flood depth, drainage, soil water, air-filled porosity | Flooded or low-aeration layers permit CH4 production; drainage releases part of stored CH4. |
| Rice-plant pathway | root length density `RLV` -> root transmissivity | Plant-mediated transport is tied to crop root distribution. |
| Redox buffer | soil alternative electron acceptors `SAEA` plus fertilizer-added buffer | Suppresses methanogenesis until alternate electron acceptors are reduced; regenerates under aerated conditions. |

## 3. Key physical ideas worth carrying into CoLM

### 3.1 Methanogenesis is substrate- and redox-buffer-limited

DSSAT uses fresh decomposition carbon as the daily substrate supply and caps CH4-C formation by carbon availability. Under sufficiently aerated conditions, CH4 production is zero and substrate remains aerobic CO2. Under anaerobic conditions, alternate electron acceptors buffer methane formation before substrate can be routed to CH4.

CoLM mapping:
- Current CoLM `methane_prod` already links CH4 production to heterotrophic respiration and water-table/aeration scalars.
- The new optional microbial diagnostics also cap microbial production potential by available `hr_vr/catomw` substrate.
- Missing if rice-specific CH4 is desired: an explicit alternate-electron-acceptor buffer state analogous to DSSAT `SAEA` / reduced buffer.

### 3.2 Flooding and drainage are first-order rice controls

DSSAT strongly conditions CH4 production on flood/soil aeration and includes a drainage pulse: when floodwater is gone, part of previously stored CH4 is emitted and storage decreases.

CoLM mapping:
- Wetland/lake CH4 uses `finundated`, saturation split, surface water, ebullition, diffusion, and transport; it is not crop-management aware.
- For paddy/rice applications, CoLM needs explicit irrigation/flood-management coupling so drained paddy soil does not behave like a generic natural wetland.

### 3.3 Root-mediated transport is crop-root driven

DSSAT converts rice root length density into root transmissivity and sends CH4 through root transport in the Arah solver.

CoLM mapping:
- CoLM already has aerenchyma transport using root/LAI/NPP-related controls.
- For rice, a crop-specific porosity/transmissivity mode should use crop root distribution and rice phenology rather than only generic wetland vegetation controls.

### 3.4 Arah transport is more mechanistic for paddy columns, but overlaps CoLM transport

DSSAT's Arah solver solves a steady-state O2/CH4 diffusion-reaction problem with:
- gas and aqueous diffusion,
- Henry/phase partitioning,
- O2-limited CH4 oxidation,
- root-mediated gas movement,
- ebullition when dissolved CH4 exceeds solubility,
- leaching with drainage.

CoLM already has most analogous process families for natural wetland/lake columns. Therefore importing the entire Arah solver would duplicate architecture and create two competing CH4 transport systems. The safer path is to add rice-specific controls around CoLM's existing production/oxidation/transport, not to transplant DSSAT code.

## 4. Parameters observed in DSSAT snapshot

These are reference values, not CoLM defaults:

| Concept | DSSAT value / behavior | CoLM status |
|---|---:|---|
| Aeration cutoff for production | no CH4 when air-filled porosity exceeds about 30% of maximum pore air space | CoLM uses saturation/water-table/redox scalars, not this paddy threshold. |
| Buffer regeneration under aeration | about `0.070 d-1` | no explicit alternate-electron-acceptor buffer state. |
| Root transmissivity coefficient | `0.00015 m air / m root` | CoLM has separate aerenchyma geometry parameters. |
| CH4 oxidation potential in Arah column | about `1.5e-5 mol m-3 s-1` | CoLM default `vmax_methane_oxid=1.25e-5 mol m-3 s-1`, close in magnitude. |
| CH4 aqueous diffusion | about `1.49e-9 m2 s-1` | CoLM water diffusivity parameterization is temperature dependent and similar order. |
| O2 inhibition on CH4 production | high inhibition coefficient in Arah kinetics | CoLM has `oxinhib` and O2/redox controls. |

## 5. Recommended CoLM integration path

This is an optional **rice/paddy extension**, separate from the CTSM wetland/lake alignment work.

| Priority | Work package | Implementation direction | Default behavior |
|---|---|---|---|
| R0 | Documentation and comparison | Keep this DSSAT reference note and link it from the CH4 alignment plan. | No model change. |
| R1 | Paddy diagnostic mode | Add diagnostics that report DSSAT-like CH4 production, oxidation/consumption, plant flux, ebullition, diffusion, leaching/storage terms from existing CoLM fluxes where possible. | Off unless rice/paddy history requested. |
| R2 | Rice water-management scalar | Add optional paddy production scalar based on flooded/drained state, air-filled porosity, and crop/rice patch identity. | Off by default; natural wetland unchanged. |
| R3 | Alternate-electron-acceptor buffer | Add prognostic oxidized/reduced redox-buffer pools for paddy CH4 suppression/recovery. | Off by default; requires restart/LULCC support. |
| R4 | Crop-root transmissivity mode | Make rice aerenchyma transport use crop root distribution/phenology and optional DSSAT-like root transmissivity coefficient. | Off by default; generic aerenchyma unchanged. |
| R5 | Arah-style paddy column solver | Only if site validation shows CoLM transport cannot reproduce paddy CH4 dynamics. | Not recommended initially. |

## 6. Acceptance criteria before enabling any rice-specific mode

1. `igas_ch4 <= 0` remains a no-op.
2. Default namelist leaves existing wetland/lake CH4 unchanged.
3. Rice/paddy mode is gated by an explicit option and by rice/crop or paddy-water-management identity.
4. Carbon units are explicit: DSSAT reports kg C ha-1 d-1, while CoLM internal rates are mostly mol m-3 s-1 and mol m-2 s-1.
5. Drainage and storage pulses close the CH4-C budget.
6. Any new state has restart, LULCC/remap, history, and single-point smoke-test coverage.

## 7. Bottom line

DSSAT-CSM is useful as a **rice/paddy CH4 reference**, especially for flood/drainage management, alternate electron acceptor buffering, and crop-root transport. It should not replace the CoLM/CTSM-style natural wetland/lake CH4 core. The recommended next addition is a default-off paddy extension that reuses CoLM transport and adds DSSAT-inspired controls around substrate/redox/water-management and rice roots.
