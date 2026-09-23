# CaMa-Flood_v4
CaMa-Flood_v4 code on GitHub.

This GitHub repository is mainly used for those who want to make contributions to CaMa-Flood.
If you simply want to use CaMa-Flood, please visit the product webpage and register on Google Form. (the password for downloading is issued after registration). Then, you can download the CaMa-Flood package which contains data (map/input) to run the model.
http://hydro.iis.u-tokyo.ac.jp/~yamadai/cama-flood/

Now the codes of CaMa-Flood v4 is distributed under **Apache 2.0** license (so Open Source!). 
So you can feel free to use/modify the code. (note: data is distributed under difference license)

We are also happy to collaborate with external contributers for further development of CaMa-Flood.
We are now discussing how to merge contributions from multiple developpers.

Below is a tentative idea.

## Branches
- **master**        : main repository (latest stable code is here).
- **release_v4.XX** : the code archive released as v4.XX (Do not modify)
- **develop4XX**    : development branch based on v4.XX (please fork from here to make a contribution)
- **dev4XX_PROJ**   : working branch for merging a project (when you make pull request, we creat a specific branch for your project)


## When you make some update

0. Please contact the developper (Dai Yamazaki) in advance to make a discussion on the new scheme to be integrated to **CaMa-Food master**.

1. We recommend you to fork the **develop4XX** branch.

2. Make some modifications, following the below rules:

- The added scheme should be turned on/off by switch. Better if it does not affect the original codes when the new scheme is turned off.
- If the new scheme causes conflict, please discuss with the manager.
- We appreciate if the sample scripts is prepared to run CaMa-Flood with the new scheme. If possible, please design the scripts ready to be run by someone else. For example, avoid absolute pass, avoid environment dependent shebang (e.g. #!/usr/local/python), describe required library required.
- If the new scheme require additional data to run, please prepare a sample data required by the new sheme (better if the size of sample data is minimum). However, the data should be kept outside of GitHub. If your new scheme required data, please contact Yamazaki to discuss a better way of data management in CaMa-Flood package.

3. When your code is ready to be merged, please contact the repository manager (Dai Yamazaki). Then, we create a working branch **dev4XX_PROJ** for merging your commit. 
- Please fork **dev4XX_PROJ**, reflect your changes on this branch (only include nesessary changes).
- If you find a conflict which cannot be solved, please contact the code manager.

4. Please make a pull request to **dev4XX_PROJ** branch.
- If your commit cannot be merged by conflict, we will contact you.

5. Your new contribution is once merged to **dev4XX_PROJ**. Yamazaki will perform a test simulations to confirm your contribution works well. Then, the new scheme is merged to **develop4XX**  and **master** branch, by Yamazaki

## CoLM legacy two-way coupling

- Enable `CaMa_Flood` and `USEMPI`; `define.h` automatically disables
  `GridRiverLakeFlow` for this choice. Keep `CatchLateralFlow` disabled.
  These are alternative routing solvers, not additive ones.
- Use the **same main CoLM namelist input as GridRiverLakeFlow**:
  ```fortran
  DEF_UnitCatchment_file = './CoLMruntime/unitcatchment/grid_routing_data_15min.nc'
  DEF_USE_BIFURCATION = .false. ! enable only when bifurcation data are available
  ```
  Set the file to your external dataset. Relative paths use the process working
  directory, not the namelist directory. No dataset is bundled in the repository.
- When `DEF_UnitCatchment_file` is set, **no separate CaMa namelist is opened**:
  `DEF_CaMa_Namelist`, `CROUTINGNC`, `CDIMINFO`, all old static-map paths and
  format/endian selectors are bypassed. There is no AUTO/runtime directory scan.
  Invalid explicit inputs fail rather than falling back to another dataset.
- If the main NC path is empty / `null` / `NONE`, the original CaMa data and
  namelist path (`DEF_CaMa_Namelist`) is retained for backward compatibility,
  including old BIN and CaMa NC inputs. This legacy branch is not automatically
  converted into the new fixed-settings mode. Standalone CaMa calls are unchanged.
- In main-NC mode the fixed settings are local-inertial routing, floodplains and
  floodplain flow on, adaptive substeps, hourly coupling (`DT=3600`, `IFRQ_INP=1`),
  flood evaporation and infiltration on by default (`DEF_CaMa_FloodFeedback = .false.`
  makes the coupling one-way: the land then sees no flood water and CaMa takes
  nothing back), and the `inpmat_*` matrix from the same
  NC. Dates/leap years and history/restart scheduling come from CoLM.
  `DEF_USE_BIFURCATION` is reused. CoLM history includes basic channel/flood/sink
  diagnostics; CaMa does not open its own standalone output files.
- Additional CaMa groundwater-delay, tidal/mean-sea-level, slope-mixing,
  reservoir/irrigation, levee and legacy sediment/tracer processes are not enabled
  in the fixed main-NC mode. They require process-specific inputs/settings and
  remain available through legacy configuration. Main levee/reservoir operation
  requests in the NC mode fail explicitly rather than being silently ignored.
  CoLM MPI is supported; separate distributed-CaMa `UseMPI_CMF` is not supported
  by the bundled reader.
- Optional slope, groundwater-delay and mean-sea-level physics require their
  corresponding NC variables. Missing fields are errors, not implicit zeros.
  The bundled reader prefers the shared `topo_rivelv` (river-bed elevation) plus
  `topo_rivhgt`, converted to bank elevation as bed elevation + channel height;
  `topo_elevation` remains a legacy alternative. Downstream distance must be supplied
  as `topo_distance` or `topo_nxtdst`; it is not approximated from
  `topo_rivlen`.
- `seq_next` and bifurcation endpoints reference file rows, as in GridRiverLakeFlow.
  The reader derives upstream links and an internal river-first order in O(N),
  applying the same permutation to topology, topography, input matrix and
  bifurcations. Extra `seq`/upstream metadata is not required or authoritative.
  Matrix/profile dimension identities, not just lengths, distinguish transposes,
  including square arrays. CaMa retains its own channel/floodplain geometry:
  supplied `topo_rivstomax` must match length × width × height within rounding
  tolerance. `topo_fldstomax` is not used; CaMa derives its own total-storage curve.
  External bifurcation outlets (`bifurcation_down=0`) are not supported by this
  CaMa path and fail explicitly; ordinary river mouths remain supported.
- `DT` is the routing outer step; `IFRQ_INP` is the coupling interval in hours.
  The core fallback and supplied examples default to hourly coupling
  (`IFRQ_INP=1`, `DT=3600` s). The coupler no longer overwrites `DT`. Adaptive routing is enabled in the
  example. CoLM steps and routing outer steps must be whole minutes, matching
  the legacy CaMa calendar resolution.

### Water accounting and regional runs

Exchange is in **volumes**: covered patch runoff is integrated over intersection
areas, then assigned to routing cells by the input matrix. Missing regional area
is not filled by extrapolating the average runoff of present patches.
Sediment rainfall is an intensive mean: valid rainfall amounts are divided by
valid area-time, excluding missing patches and missing timesteps from both
numerator and denominator. This normalization is not applied to runoff or
floodwater withdrawals.

Flood storage is distributed using row-normalized input-matrix weights and true
regular-grid areas. A patch consumes a finite credit. Its withdrawals return to
the same supplying routing cells and are debited **before** routing; substeps
must not debit them again. Input-matrix weights do not replace `topo_area`.
Uncovered flood water remains in CaMa.

This does not identify the actually submerged vegetation patches: the allocation
assumes uniform inundation within each mapped grid intersection. Resolving which
patch floods first requires subgrid elevation/floodplain information.

The routing network is not cropped at the CoLM rectangle. By default
(`DEF_CaMa_StrictDomain = .false.`) a regional domain need not hold complete
drainage basins: CaMa warns how many river links or bifurcation ends cross the
domain edge, routes the runoff that has a routing cell, drops the runoff of grid
cells that no routing cell receives, and reports the dropped volume as a share of
the land runoff at the end of the run. Water entering across the domain edge is not
supplied; this legacy coupler does not invent external inflow. With
`DEF_CaMa_StrictDomain = .true.` such a run stops instead; supply complete drainage
basins through their outlets. Partial mapped cells contribute only their
covered volume.

### Land energy, restart and diagnostics

Flood evaporation is included before the ground-temperature solve. Only the dry
land contribution is removed from soil/snow; flood evaporation is taken from its
CaMa credit. This remains a single-ground-temperature approximation, not a
separate floodwater thermal/ice model. `LWEVAPFIX` / `LWINFILTFIX` are deprecated
compatibility options with no physical effect; enabling either prints a warning.
Enabled sinks always use the same credit/debit accounting regardless of these
flags. Negative open-water evaporation is set to zero (no CaMa condensation source). Flood infiltration bypasses the canopy.
Urban exchange has not been revised in this repair. Infiltration with active
water-transported tracers/isotopes is rejected: this legacy exchange does not
provide their floodwater composition. Compiling `TRACER` alone, or using
non-water-transported species, does not disable water-only coupling.
Re-infiltration uses the patch's single shared soil column; it does not track
separate saturation beneath the flooded fraction. Wet-area infiltration rates
therefore need site-scale validation rather than interpretation as a resolved
subgrid soil process.

Every CoLM restart and final step flushes the partial coupling window before
writing routing storage, so no pending land fluxes are lost at restart. Main-NC
mode writes NetCDF checkpoints. For a warm restart set the main namelist field
`DEF_CaMa_Restart_file` to the matching CaMa checkpoint; leaving it `null` means
an explicit cold start of routing and prints a warning. CoLM reads land initial
states for both cold and warm starts, so that read alone cannot identify a resume.
The legacy branch still uses `LRESTART=.TRUE.` and `CRESTSTO`. NetCDF restart
checks the checkpoint time and routing-grid coordinates before loading state,
and retains CaMa history accumulators across an intra-day restart. Older files
without those accumulators are accepted only for DAILY history at a day boundary
after spinup; otherwise they fail rather than silently truncate history. General land-history accumulators are saved separately before a forced partial
record is normalized, then restored after history initialization; the partial
record is still written. Legacy checkpoints without this sidecar warn that an
intra-period land record may be partial; new checkpoints require it. Checkpoint
files at the same timestamp must not be mixed across runs or overwritten during
saving. The sidecar is committed only after the physical restart and its
requirement marker are written; interrupted or markerless sidecars are rejected.
There is no automatic checkpoint search. Flood credits are published immediately after
restore. Saving a checkpoint can split an otherwise longer coupling interval;
compare restart continuity using the same checkpoint schedule. A history time
with no routing sample is left unwritten rather than dividing by zero.

Regression checks (Fortran compiler and NetCDF required for runtime probes):

```sh
python -m pytest -q tests/test_legacy_cama*.py
# Include actual 4-rank mapping checks with a fresh CaMa-enabled build:
COLM_CAMA_BUILD_DIR=/path/to/coupled-build/.bld python -m pytest -q tests/test_legacy_cama_mpi.py
```
