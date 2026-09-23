"""Focused checks for the routing support, outgoing-face limiter, and volume."""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


FLOW = Path(__file__).resolve().parents[1] / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90"


def test_runoff_support_is_refreshed_before_both_remaps(tmp_path):
    flow = FLOW.read_text()
    step = flow.split("SUBROUTINE grid_riverlake_flow (", 1)[1]
    mask = step.index("filter_rnof = patchtype < 99 .and. patchmask")
    forcing = step.index("filter_rnof = filter_rnof .and. forcmask_pch", mask)
    water = step.index("remap_patch2inpm, rnof, rnof_gd", forcing)
    tracer = step.index("remap_patch2inpm, trc_rnof_step(itrc, :)", water)
    assert mask < forcing < water < tracer
    assert "filter = filter_rnof" in step[water:tracer]
    assert "filter = filter_rnof" in step[tracer:tracer + 150]
    fc = require_runnable_fortran_compiler(tmp_path)
    refresh = step[mask:water].split("\n", 2)[:2]
    code = """
program check
 implicit none
 integer, parameter :: r8=kind(1.d0)
 integer :: patchtype(2)=[1,1]
 logical :: patchmask(2)=.true., forcmask_pch(2), filter_rnof(2)
 real(r8) :: rnof(2)=[2._r8,3._r8], tracer(2)=[5._r8,7._r8]
 type forcing_type
  logical :: has_missing_value=.true.
 end type
 type(forcing_type) :: DEF_forcing
 forcmask_pch=[.true.,.true.]
@REFRESH@
 if (sum(rnof,mask=filter_rnof)/=5._r8) stop 1
 forcmask_pch(2)=.false.
@REFRESH@
 if (sum(rnof,mask=filter_rnof)/=2._r8) stop 2
 if (sum(tracer,mask=filter_rnof)/=5._r8) stop 3
 print *, 'HYDRO_MASK_OK'
end program
""".replace("@REFRESH@", "\n".join(refresh))
    source = tmp_path / "mask.f90"
    source.write_text(code)
    exe = tmp_path / "mask"
    compiled = subprocess.run([fc, str(source), "-o", str(exe)], capture_output=True, text=True)
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    ran = subprocess.run([str(exe)], capture_output=True, text=True)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert "HYDRO_MASK_OK" in ran.stdout


def test_particle_storage_uses_advanced_hydro_volume():
    flow = FLOW.read_text()
    block = flow.split("particle_floodarea(i) = total_floodarea(i)", 1)[1].split(
        "CALL tracer_lifecycle_route_diag_accumulate", 1
    )[0]
    assert "max(volresv(irsv), 0._r8)" in block
    assert "max(volwater_ucat(i), 0._r8)" in block
    assert "particle_water_storage(i) + max(levsto(i), 0._r8)" in block
    assert "topo_rivwth(i) * topo_rivlen(i)" not in block


def test_reactive_source_resets_outside_debug():
    flow = FLOW.read_text()
    step = flow.split("SUBROUTINE grid_riverlake_flow (", 1)[1]
    reset = step.index("IF (allocated(trc_reactive_source)) trc_reactive_source = 0._r8")
    debug = step.index("#ifdef CoLMDEBUG", reset)
    assert reset < debug


def test_inland_reservoir_overflow_updates_release_diagnostic():
    flow = FLOW.read_text()
    overflow = flow.split("! Remove excess water after transport has already run.", 1)[1]
    discharge = overflow.index("hflux_fc(i) = (volwater - topo_rivstomax(i)) / dt_all(irivsys(i))")
    diagnostic = overflow.index("IF (is_built_resv(i)) qresv_out(ucat2resv(i)) = hflux_fc(i)")
    assert discharge < diagnostic < overflow.index("volwater = topo_rivstomax(i)")


def test_empty_patch_worker_skips_unallocated_land_arrays_but_still_pushes():
    flow = FLOW.read_text().split("SUBROUTINE grid_riverlake_flow (", 1)[1]
    runoff = flow.split("! Forcing coverage can change", 1)[1].split("#ifdef TRACER", 1)[0]
    assert "IF (numpatch > 0) THEN" in runoff
    assert runoff.index("ENDIF") < runoff.index("CALL worker_push_data (push_inpm2ucat")
    rain = flow.split("allocate(prcp_pch(numpatch), filter_prcp(numpatch))", 1)[1]
    assert rain.index("IF (numpatch > 0) THEN") < rain.index("ieee_is_finite(forc_prc)")
    assert rain.index("ENDIF") < rain.index("CALL worker_push_data(push_inpm2ucat, prcp_push_fields")


def test_tracer_cold_start_uses_dam_construction_year():
    flow = FLOW.read_text()
    init = flow.split("SUBROUTINE grid_riverlake_flow_init (start_year)", 1)[1].split(
        "END SUBROUTINE grid_riverlake_flow_init", 1
    )[0]
    assert "IF (lake_type(i) /= 2) CYCLE" in init
    assert "IF (irsv < 1 .or. irsv > size(dam_build_year)) CYCLE" in init
    assert "is_built_resv_init(i) = start_year >= dam_build_year(irsv)" in init
    assert "trc_missing, is_built_resv_init)" in init


def test_small_donor_multi_face_limiter_conserves_reverse_transfer(tmp_path):
    fc = require_runnable_fortran_compiler(tmp_path)
    flow = FLOW.read_text()
    limiter = flow.split("            normal_outgoing_rate = 0._r8\n", 1)[1].split(
        "            IF (DEF_USE_BIFURCATION) THEN", 1
    )[0]
    momentum = limiter.split("            ! The face limiter also scales momentum fluxes.", 1)[1]
    momentum = momentum.split("            IF (DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT) THEN", 1)[1]
    momentum = "            IF (DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT) THEN" + momentum
    limiter = ("            normal_outgoing_rate = 0._r8\n" + limiter).split(
        "            ! The face limiter also scales momentum fluxes.", 1
    )[0]
    assert "ordinary_scale_next(i)" in limiter
    assert "sum_hflux_riv = hflux_fc - hflux_sumups" in limiter
    # Replace only communication calls with an in-process chain mapping; the
    # donor/face limiter and net-flux reconstruction are compiled from source.
    limiter = limiter.replace(
        "CALL worker_push_data(push_ups2ucat, hflux_sumups, mflux_sumups, fillvalue = 0._r8, mode = 'sum')",
        "CALL push_up(hflux_sumups, mflux_sumups)",
    ).replace(
        "CALL worker_push_data(push_next2ucat, ordinary_scale, ordinary_scale_next, fillvalue = 1._r8)",
        "CALL push_next(ordinary_scale, ordinary_scale_next)",
    ).replace(
        "CALL worker_push_data(push_ups2ucat, upstream_flux_fields, mode = 'sum')",
        "CALL push_up(hflux_fc, hflux_sumups)\n            CALL push_up(mflux_fc, mflux_sumups)",
    )
    code = """
module probe
 use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
 implicit none
 integer, parameter :: r8=kind(1.d0), numucat=3
 logical :: ucatfilter(3)=.true., is_built_resv(3)=.false.
 integer :: irivsys(3)=[1,1,1], ucat2resv(3)=[1,1,1]
 real(r8) :: hflux_fc(3), mflux_fc(3), hflux_sumups(3), mflux_sumups(3)
 real(r8) :: sum_hflux_riv(3), sum_mflux_riv(3), normal_outgoing_rate(3)
 real(r8) :: ordinary_scale(3), ordinary_scale_next(3), volwater_ucat(3)
 real(r8) :: volresv(1)=0, qresv_in(1)=0, qresv_out(1)=0
 real(r8) :: dt_all(1)=1, volwater, dt_this
 real(r8) :: veloc_riv(3)=0, momen_riv(3)=0, topo_rivare(3)=1, sum_zgrad_riv(3)=0
 logical :: DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT=.false.
 contains
 subroutine push_up(send,recv)
 real(r8), intent(in) :: send(:)
 real(r8), intent(out) :: recv(:)
 recv=[0._r8,send(1),send(2)]
 end subroutine
 subroutine push_next(send,recv)
 real(r8), intent(in) :: send(:)
 real(r8), intent(out) :: recv(:)
 recv=[send(2),send(3),1._r8]
 end subroutine
 subroutine CoLM_stop(message)
 character(*), intent(in) :: message
 print *, trim(message)
 error stop 13
 end subroutine
 subroutine limit_faces()
 integer :: i
@LIMITER@
@MOMENTUM@
 end subroutine
end module
program check
 use probe
 implicit none
 real(r8) :: after(3), before(3)
 character(8) :: mode
 call get_command_argument(1,mode)
 if (trim(mode)=='bad') then
  volwater_ucat=[1._r8,4._r8,1._r8]
  hflux_fc=[0._r8,3._r8,0._r8]
  mflux_fc=hflux_fc
  veloc_riv=[0._r8,1._r8,0._r8]
  momen_riv=0._r8
  DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT=.true.
  call limit_faces()
  stop 99
 endif
 ! Cell 2 loses through two faces: reverse on 1->2 and forward on 2->3.
 before=[1._r8,1.e-7_r8,1._r8]
 volwater_ucat=before
 hflux_fc=[-2._r8,3._r8,0._r8]
 mflux_fc=hflux_fc
 call limit_faces()
 after=before-sum_hflux_riv*dt_all(1)
 if (abs(ordinary_scale(2)-2.e-8_r8)>1.e-20_r8) stop 1
 if (after(2)<-1.e-20_r8) stop 2
 if (abs(sum(after)-sum(before))>1.e-12_r8) stop 3
 if (abs(hflux_fc(1)+4.e-8_r8)>1.e-20_r8) stop 4
 if (abs(hflux_fc(2)-6.e-8_r8)>1.e-20_r8) stop 5
 if (any(abs(mflux_fc-hflux_fc)>1.e-20_r8)) stop 6
 ! A dry downstream donor cannot reverse-release water despite upstream inflow.
 volwater_ucat=[1._r8,1._r8,0._r8]
 hflux_fc=[0._r8,-1._r8,0._r8]
 mflux_fc=hflux_fc
 call limit_faces()
 if (any(hflux_fc/=0._r8)) stop 7
 ! Reservoir inflow diagnostic must follow the finally capped upstream face.
 is_built_resv(2)=.true.
 volwater_ucat=[1.e-7_r8,0._r8,1._r8]
 hflux_fc=[2._r8,0._r8,0._r8]
 mflux_fc=hflux_fc
 call limit_faces()
 if (abs(qresv_in(1)-1.e-7_r8)>1.e-20_r8) stop 8
 ! Pre-cap momentum force is 3-2-.5=.5 (dt=1 acceptable), but
 ! capping the incoming face makes final force nearly 2.5.
 is_built_resv=.false.
 volwater_ucat=[1.e-7_r8,4._r8,1._r8]
 hflux_fc=[2._r8,3._r8,0._r8]
 mflux_fc=hflux_fc
 veloc_riv=[0._r8,1._r8,0._r8]
 momen_riv=veloc_riv
 sum_zgrad_riv=[0._r8,0.5_r8,0._r8]
 DEF_GRIDBASED_ROUTING_MOMENTUM_DT_LIMIT=.true.
 dt_all=1._r8
 call limit_faces()
 if (dt_all(1)>=0.5_r8) stop 9
 if (momen_riv(2)-(sum_mflux_riv(2)-sum_zgrad_riv(2))*dt_all(1)<-1.e-12_r8) stop 10
 print *, 'HYDRO_LIMITER_OK'
end program
""".replace("@LIMITER@", limiter).replace("@MOMENTUM@", momentum)
    source = tmp_path / "probe.f90"
    source.write_text(code)
    exe = tmp_path / "probe"
    compiled = subprocess.run(
        [fc, "-ffree-line-length-0", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow", str(source), "-o", str(exe)],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    ran = subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert "HYDRO_LIMITER_OK" in ran.stdout
    bad = subprocess.run([str(exe), "bad"], cwd=tmp_path, capture_output=True, text=True)
    assert bad.returncode != 0
    assert "non-positive final momentum dt" in bad.stdout
