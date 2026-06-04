# Chloride wet-deposition forcing for TRACER

This folder builds an offline CoLM TRACER input for chloride (Cl-) wet deposition.
The model-side precip forcing expects a concentration-like value:

Model-side note: `DEF_TRACER_TYPES = "conservative"` is a generic
nonvolatile dissolved-solute path. It is not chloride-specific: precipitation
and runoff carry the tracer; evaporation, sublimation, transpiration, dew, and
frost carry zero tracer unless a separate deposition forcing is added.


```text
cl_precip_conc = chloride wet-deposition amount / precipitation water amount
units          = kg Cl kg-1 H2O
```

Recommended workflow:

1. Prepare chloride wet-deposition totals or fluxes on the forcing grid.
2. Prepare precipitation totals or fluxes for the same time records/grid.
3. Run `make_chloride_precip_conc.py` to write `CL_CONC_YYYY_MM.nc` files.
4. Point the Cl tracer parameter file at `forcing_role='precip'`,
   `forcing_fprefix='CL_CONC'`, `forcing_vname='cl_precip_conc'`,
   `forcing_input_mode='value'`.

Example for monthly totals:

```bash
python preprocess/Forcings/ChlorideWetDep/make_chloride_precip_conc.py \
  --wetdep nadp_cldep_monthly.nc --wetdep-var cl_wetdep --wetdep-unit mg_m-2 \
  --precip precip_monthly.nc --precip-var precip --precip-unit mm \
  --out-dir /path/to/CoLM_forcing/tracer \
  --fprefix CL_CONC --groupby month \
  --output-dtime-seconds 21600 --write-diagnostics
```

Output naming follows `main/TRACER/MOD_Tracer_Forcing.F90` for tracer-owned
forcing files, not the dataset-specific `MOD_UserSpecifiedForcing` names:

```text
month: CL_CONC_2006_04.nc
year : CL_CONC_2006.nc
day  : CL_CONC_2006_04_14.nc
```

Use monthly wetdep/monthly precipitation if the source product is monthly.
That distributes the monthly chloride input across the model precipitation
within that month through a fixed precipitation concentration. Do not divide a
monthly wetdep total by an instantaneous model rain rate offline.

For CoLM runtime reading, prefer `--output-dtime-seconds` equal to the tracer
forcing `forcing_dtime` (often 21600 for 6-hourly forcing). A monthly source
record is then repeated to every model-forcing time slot in that month, avoiding
a one-record-per-month file being indexed as if it contained submonthly records.
