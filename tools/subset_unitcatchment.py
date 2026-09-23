#!/usr/bin/env python3
"""Cut a regional unit-catchment routing network out of the global file.

mksrfdata does this itself when DEF_UnitCatchment_regional = .true. (module
mksrfdata/MOD_UnitCatchmentSubset.F90); this script is the reference
implementation the tests compare it against, and it also converts a routing
restart to a regional network (see below).

GridRiverLakeFlow routes the whole network in DEF_UnitCatchment_file.  For a
regional land domain almost all of it is dry and receives no runoff, yet every
worker still steps it.  This tool keeps only the river systems that receive
runoff from the land domain and writes a smaller file with the same layout, so
DEF_UnitCatchment_file can simply point at it.

A river system is everything that drains to one river mouth (or inland
depression).  Systems are kept whole, so no kept cell has a downstream or
upstream neighbour outside the file.  Bifurcation pathways may join different
systems; by default the selection is grown until every pathway with one end in
the set has both ends in it (--bif drop discards straddling pathways instead,
which is only valid for runs with DEF_USE_BIFURCATION = .false.).

Which cells receive runoff is decided with the file's own runoff-input matrix
(inpmat_x / inpmat_y): a unit catchment receives runoff when one of its input
cells contains part of the land domain.  The domain is given either as a
bounding box or as a 2-D mask file (positive cells belong to the domain).

Examples
    subset_unitcatchment.py network global.nc pearl.nc --bbox 102 115 21.5 27
    subset_unitcatchment.py network global.nc pearl.nc --mask mesh.nc --mask-var elmindex \\
        --lon-var longitude --lat-var latitude
    subset_unitcatchment.py restart pearl.nc global_restart_gridriver.nc pearl_restart_gridriver.nc

A run on the regional file needs its own routing restart: the restart is stored
per unit catchment and is refused on a network of a different size.  The
'restart' command carries a spun-up global state over to the regional file.
"""
import argparse
import sys

import numpy as np
import netCDF4 as nc

SEQ_DIM = "nseqmax"
PATH_DIM = "npthout"
DAM_DIM = "dam_ndams"


def domain_cells(src, args):
    """Boolean [nx+1, ny+1] array (1-based) of runoff-grid cells touched by the domain."""
    nx, ny = len(src.dimensions["nx"]), len(src.dimensions["ny"])
    west, north = float(src.west), float(src.north)
    dx = (float(src.east) - west) / nx
    dy = (north - float(src.south)) / ny

    cells = np.zeros((nx + 2, ny + 2), dtype=bool)

    if args.bbox:
        lon0, lon1, lat0, lat1 = args.bbox
        x0 = int(np.floor((lon0 - west) / dx)) + 1
        x1 = int(np.floor((lon1 - west) / dx)) + 1
        y0 = int(np.floor((north - lat1) / dy)) + 1
        y1 = int(np.floor((north - lat0) / dy)) + 1
        cells[max(x0, 1):min(x1, nx) + 1, max(y0, 1):min(y1, ny) + 1] = True
        return cells

    with nc.Dataset(args.mask) as mask:
        mask.set_auto_mask(False)
        field = mask[args.mask_var][:]
        lon = np.asarray(mask[args.lon_var][:], dtype=float)
        lat = np.asarray(mask[args.lat_var][:], dtype=float)
    if field.shape != (len(lat), len(lon)):
        sys.exit(f"mask variable {args.mask_var} has shape {field.shape}, "
                 f"expected ({len(lat)}, {len(lon)})")
    iy, ix = np.nonzero(field > 0)
    if iy.size == 0:
        sys.exit("the mask contains no positive cells")
    gx = np.floor((lon[ix] - west) / dx).astype(int) + 1
    gy = np.floor((north - lat[iy]) / dy).astype(int) + 1
    inside = (gx >= 1) & (gx <= nx) & (gy >= 1) & (gy <= ny)
    cells[gx[inside], gy[inside]] = True
    return cells


def river_mouth(seq_next):
    """0-based index of the mouth cell of every cell (seq_next > own index by construction)."""
    n = len(seq_next)
    mouth = np.arange(n)
    if np.any(seq_next[seq_next > 0] <= (np.nonzero(seq_next > 0)[0] + 1)):
        sys.exit("seq_next is not ordered upstream-to-downstream; unsupported file")
    for i in range(n - 1, -1, -1):
        if seq_next[i] > 0:
            mouth[i] = mouth[seq_next[i] - 1]
    return mouth


def select_systems(src, cells, bif_mode):
    seq_next = src["seq_next"][:].astype(np.int64)
    n = len(seq_next)
    mouth = river_mouth(seq_next)

    inp_x = src["inpmat_x"][:]
    inp_y = src["inpmat_y"][:]
    receives = np.zeros(n, dtype=bool)
    for k in range(inp_x.shape[1]):
        ok = (inp_x[:, k] > 0) & (inp_y[:, k] > 0)
        receives[ok] |= cells[inp_x[ok, k], inp_y[ok, k]]
    if not receives.any():
        sys.exit("no unit catchment receives runoff from the domain; check the domain and its coordinates")

    systems = set(np.unique(mouth[receives]).tolist())

    if bif_mode == "closure" and PATH_DIM in src.dimensions:
        up = mouth[src["bifurcation_upst"][:] - 1]
        dn = mouth[src["bifurcation_down"][:] - 1]
        while True:
            member = np.isin(up, list(systems)) | np.isin(dn, list(systems))
            grown = systems | set(up[member].tolist()) | set(dn[member].tolist())
            if grown == systems:
                break
            systems = grown

    keep = np.isin(mouth, list(systems))
    return keep, receives, systems


def copy_attrs(var_in, var_out):
    for name in var_in.ncattrs():
        if name != "_FillValue":
            var_out.setncattr(name, var_in.getncattr(name))


def write_subset(src, dst, keep, bif_mode):
    n = len(keep)
    new_index = np.zeros(n + 1, dtype=np.int64)          # 1-based old -> 1-based new, 0 = dropped
    new_index[1:][keep] = np.arange(1, keep.sum() + 1)

    def remap(values):
        """Renumber positive cell indices; non-positive codes (-9, -10, -9999) pass through."""
        values = np.asarray(values, dtype=np.int64)
        out = values.copy()
        pos = values > 0
        out[pos] = new_index[values[pos]]
        if np.any(out[pos] == 0):
            sys.exit("internal error: a kept cell references a dropped cell")
        return out

    keep_path = np.zeros(0, dtype=bool)
    if PATH_DIM in src.dimensions:
        up = src["bifurcation_upst"][:]
        dn = src["bifurcation_down"][:]
        both = keep[up - 1] & keep[dn - 1]
        straddle = keep[up - 1] ^ keep[dn - 1]
        if bif_mode == "closure" and straddle.any():
            sys.exit("internal error: bifurcation closure left straddling pathways")
        keep_path = both

    keep_dam = np.zeros(0, dtype=bool)
    if DAM_DIM in src.dimensions:
        keep_dam = keep[src["dam_seq"][:] - 1]

    sizes = {SEQ_DIM: int(keep.sum()), PATH_DIM: int(keep_path.sum()), DAM_DIM: int(keep_dam.sum())}
    selectors = {SEQ_DIM: keep, PATH_DIM: keep_path, DAM_DIM: keep_dam}

    for name, dim in src.dimensions.items():
        dst.createDimension(name, sizes.get(name, len(dim)))

    for name in src.ncattrs():
        dst.setncattr(name, src.getncattr(name))

    index_vars = {"seq": remap, "seq_next": remap, "seq_upst": remap,
                  "bifurcation_upst": remap, "bifurcation_down": remap, "dam_seq": remap}

    for name, var in src.variables.items():
        fill = getattr(var, "_FillValue", None)
        out = dst.createVariable(name, var.dtype, var.dimensions, fill_value=fill,
                                 zlib=var.ndim > 0 and var.size > 1000)
        copy_attrs(var, out)
        data = var[:]
        for axis, dim in enumerate(var.dimensions):
            if dim in selectors and selectors[dim].size:
                data = np.compress(selectors[dim], data, axis=axis)
        if name in index_vars:
            data = remap(data).astype(var.dtype)
        out[:] = data

    src_index = dst.createVariable("seq_src_index", "i4", (SEQ_DIM,))
    src_index.long_name = "index of this unit catchment in the source network"
    src_index.units = "1"
    src_index[:] = np.nonzero(keep)[0] + 1
    dst.setncattr("source_nseqmax", np.int64(n))

    n_new = int(keep.sum())
    n_river = int(np.count_nonzero(keep[:int(src.nseqriv)]))
    dst.setncattr("nseqall", np.int64(n_new))
    dst.setncattr("nseqmax", np.int64(n_new))
    dst.setncattr("nseqriv", np.int64(n_river))
    if PATH_DIM in src.dimensions:
        dst.setncattr("npthout", np.int64(sizes[PATH_DIM]))
    if DAM_DIM in src.dimensions:
        dst.setncattr("dam_ndams", np.int64(sizes[DAM_DIM]))
    dst.setncattr("title", str(src.title) + " (regional subset)")
    dst.setncattr("subset_bif_mode", bif_mode)
    return sizes


def verify(path):
    """Structural checks on the written file."""
    with nc.Dataset(path) as f:
        f.set_auto_mask(False)
        n = len(f.dimensions[SEQ_DIM])
        seq = f["seq"][:]
        nxt = f["seq_next"][:].astype(np.int64)
        upst = f["seq_upst"][:].astype(np.int64)          # (upnmax, n)
        upn = f["seq_upn"][:]
        assert np.array_equal(seq, np.arange(1, n + 1)), "seq is not 1..N"
        pos = nxt > 0
        assert np.all(nxt[pos] > np.nonzero(pos)[0] + 1), "seq_next not upstream-to-downstream"
        assert nxt[pos].max(initial=0) <= n, "seq_next points outside the file"
        assert upst.max(initial=0) <= n, "seq_upst points outside the file"
        listed = upst > 0
        assert np.array_equal(listed.sum(axis=0), upn), "seq_upn disagrees with seq_upst"
        cols = np.broadcast_to(np.arange(1, n + 1), upst.shape)
        assert np.all(nxt[upst[listed] - 1] == cols[listed]), "seq_upst is not the inverse of seq_next"
        if PATH_DIM in f.dimensions and len(f.dimensions[PATH_DIM]):
            assert f["bifurcation_upst"][:].max() <= n and f["bifurcation_down"][:].max() <= n
        if DAM_DIM in f.dimensions and len(f.dimensions[DAM_DIM]):
            assert f["dam_seq"][:].max() <= n


RESTART_DIM = "ucatch"
IDENTITY_VAR = "gridriver_ucatch_identity"


def convert_restart(subset_path, restart_in, restart_out):
    """Cut a GridRiverLake restart written on the source network down to the subset.

    Only restarts whose variables live on the unit-catchment dimension (plus the
    identity table) can be converted.  Levee, bifurcation, reservoir and tracer
    state sits on other dimensions and would have to be remapped as well, so those
    restarts are refused instead of being converted wrongly.
    """
    with nc.Dataset(subset_path) as sub:
        sub.set_auto_mask(False)
        if "seq_src_index" not in sub.variables:
            sys.exit(f"{subset_path} has no seq_src_index; regenerate it with this tool")
        src_index = sub["seq_src_index"][:].astype(np.int64) - 1
        source_n = int(sub.source_nseqmax)
        x, y, nxt = (sub[v][:].astype(float) for v in ("seq_x", "seq_y", "seq_next"))

    with nc.Dataset(restart_in) as rin:
        rin.set_auto_mask(False)
        if len(rin.dimensions[RESTART_DIM]) != source_n:
            sys.exit(f"restart has {len(rin.dimensions[RESTART_DIM])} unit catchments, "
                     f"the source network has {source_n}")
        for name, var in rin.variables.items():
            extra = [d for d in var.dimensions if d not in (RESTART_DIM, "gridriver_ucatch_identity_field")]
            if extra:
                sys.exit(f"variable {name} is on dimension {extra[0]}: restarts with levee, bifurcation, "
                         "reservoir or tracer state cannot be converted by this tool")

        with nc.Dataset(restart_out, "w", format="NETCDF4") as rout:
            for name, dim in rin.dimensions.items():
                rout.createDimension(name, len(src_index) if name == RESTART_DIM else len(dim))
            for name in rin.ncattrs():
                rout.setncattr(name, rin.getncattr(name))
            for name, var in rin.variables.items():
                fill = getattr(var, "_FillValue", None)
                out = rout.createVariable(name, var.dtype, var.dimensions, fill_value=fill)
                copy_attrs(var, out)
                data = var[:]
                if RESTART_DIM in var.dimensions:
                    data = np.take(data, src_index, axis=var.dimensions.index(RESTART_DIM))
                if name == IDENTITY_VAR:
                    data = np.array(data, dtype=float)
                    data[:, 1], data[:, 2], data[:, 3] = x, y, nxt
                out[:] = data
    print(f"wrote {restart_out} ({len(src_index)} of {source_n} unit catchments)")


def add_domain_arguments(ap):
    where = ap.add_mutually_exclusive_group(required=True)
    where.add_argument("--bbox", nargs=4, type=float, metavar=("LON0", "LON1", "LAT0", "LAT1"),
                       help="land domain as a bounding box (degrees)")
    where.add_argument("--mask", help="NetCDF file with a 2-D field whose positive cells are the land domain")
    ap.add_argument("--mask-var", default="elmindex")
    ap.add_argument("--lon-var", default="longitude")
    ap.add_argument("--lat-var", default="latitude")
    ap.add_argument("--bif", choices=("closure", "drop"), default="closure",
                    help="closure: keep systems joined by bifurcation pathways (default); "
                         "drop: discard straddling pathways (bifurcation must be off)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)

    net = sub.add_parser("network", help="write a regional unit-catchment network file")
    net.add_argument("source", help="global unit-catchment network file")
    net.add_argument("output", help="regional network file to write")
    add_domain_arguments(net)

    rst = sub.add_parser("restart", help="cut a GridRiverLake restart down to a regional network")
    rst.add_argument("subset", help="regional network file written by the 'network' command")
    rst.add_argument("restart_in", help="restart file written on the source network")
    rst.add_argument("restart_out")

    args = ap.parse_args(argv)

    if args.command == "restart":
        convert_restart(args.subset, args.restart_in, args.restart_out)
        return

    with nc.Dataset(args.source) as src:
        src.set_auto_mask(False)
        cells = domain_cells(src, args)
        keep, receives, systems = select_systems(src, cells, args.bif)
        with nc.Dataset(args.output, "w", format="NETCDF4") as dst:
            sizes = write_subset(src, dst, keep, args.bif)
        total = len(keep)

    verify(args.output)
    print(f"domain touches {int(cells.sum())} runoff-grid cells; {int(receives.sum())} unit catchments receive runoff")
    print(f"kept {len(systems)} river systems: {sizes[SEQ_DIM]} of {total} unit catchments "
          f"({100.0 * sizes[SEQ_DIM] / total:.2f}%), {sizes[PATH_DIM]} bifurcation pathways, {sizes[DAM_DIM]} dams")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
