from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
SUBPROCESS_TIMEOUT = 60


def test_malformed_sediment_namelist_reports_iomsg_field_name(tmp_path: Path) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)
    (tmp_path / "driver.f90").write_text(
        """
program sediment_namelist_driver
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
  integer, parameter :: MAX_SED_PARAM_CLASSES = 100
  type :: sediment_parameter_type
    integer  :: nsed = -1
    real(r8) :: grain_diameter(MAX_SED_PARAM_CLASSES) = -1._r8
    real(r8) :: grain_density = -1._r8
    real(r8) :: water_density = -1._r8
    real(r8) :: porosity = -1._r8
    integer  :: ndeposit_layers = -1
    real(r8) :: ignore_depth_m = -1._r8
    real(r8) :: active_layer_depth = -1._r8
    real(r8) :: viscosity = -1._r8
    real(r8) :: von_karman = -1._r8
    real(r8) :: settling_multiplier = -1._r8
    real(r8) :: yield_coefficient = -1._r8
    real(r8) :: slope_exponent = -1._r8
    real(r8) :: precipitation_exponent = -1._r8
    real(r8) :: unit_conversion = -1._r8
    real(r8) :: cfl_adv = -1._r8
    real(r8) :: max_timestep_s = -1._r8
    real(r8) :: bed_depth = -1._r8
  end type sediment_parameter_type
  type(sediment_parameter_type) :: DEF_SEDIMENT
  character(len=512) :: file_param, iomsg
  integer :: ierr, unit_nml
  namelist /nl_colm_sediment_parameter/ DEF_SEDIMENT

  call get_command_argument(1, file_param)
  open(newunit=unit_nml, status='OLD', file=trim(file_param), form='FORMATTED')
  iomsg = ''
  read(unit_nml, nml=nl_colm_sediment_parameter, iostat=ierr, iomsg=iomsg)
  close(unit_nml)
  if (ierr /= 0) then
    write(*,'(A,A)') 'ERROR: invalid &nl_colm_sediment_parameter in ', trim(file_param)
    write(*,'(A)') trim(iomsg)
    error stop 1
  endif
  write(*,'(A)') 'VALID'
end program sediment_namelist_driver
""",
        encoding="utf-8",
    )
    exe = tmp_path / "sediment_namelist_driver"
    compiled = subprocess.run(
        [compiler, "-ffree-line-length-0", str(tmp_path / "driver.f90"), "-o", str(exe)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr

    bad = tmp_path / "bad_sediment.nml"
    bad.write_text(
        """
&nl_colm_sediment_parameter
  DEF_SEDIMENT%unknown_field = 1
/
""",
        encoding="utf-8",
    )
    result = subprocess.run(
        [str(exe), str(bad)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert result.returncode != 0
    assert "invalid &nl_colm_sediment_parameter" in result.stdout
    assert "unknown_field" in result.stdout.lower()


def test_special_patch_evaporation_split_runtime_conserves_cap_and_branches(tmp_path: Path) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)
    (tmp_path / "define.h").write_text("#define TRACER\n", encoding="utf-8")
    (tmp_path / "precision.f90").write_text(
        "module MOD_Precision\n implicit none\n integer,parameter::r8=selected_real_kind(12)\nend module\n",
        encoding="utf-8",
    )
    (tmp_path / "defs.f90").write_text(
        """
module MOD_Tracer_Defs
 use MOD_Precision
 implicit none
 integer :: ntracers=3
 real(r8),parameter :: trc_tiny=1e-12_r8, trc_water_min_for_ratio=1e-12_r8
 ! Mirrors the real MOD_Tracer_Defs: parameter value and ref_ratio default
 ! must match production, or this driver would exercise different numbers.
 real(r8),parameter :: trc_delta_sanity_max=2.0e3_r8
 type :: tracer_info_type
   real(r8) :: ref_ratio = 1.0_r8
 end type tracer_info_type
 type(tracer_info_type) :: tracers(3)
 logical :: can_fixed(3)=(/.false.,.false.,.true./), uses_transport(3)=.true.
 logical :: nonvol(3)=(/.false.,.true.,.false./), has_limit(3)=.false.
 real(r8)::init_ratio(3)=(/0._r8,0._r8,0.25_r8/)
contains
 logical function tracer_can_use_fixed_signature(i); integer,intent(in)::i; tracer_can_use_fixed_signature=can_fixed(i); end
 logical function tracer_uses_land_water_transport(i); integer,intent(in)::i; tracer_uses_land_water_transport=uses_transport(i); end
 logical function tracer_is_nonvolatile_solute(i); integer,intent(in)::i; tracer_is_nonvolatile_solute=nonvol(i); end
 logical function tracer_has_dissolved_limit(i); integer,intent(in)::i; tracer_has_dissolved_limit=has_limit(i); end
 real(r8) function tracer_init_water_ratio(i); integer,intent(in)::i; tracer_init_water_ratio=init_ratio(i); end
 subroutine tracer_equilibrate_dissolved(i,w,d,s); integer,intent(in)::i; real(r8),intent(in)::w; real(r8),intent(inout)::d,s; end
end module
""",
        encoding="utf-8",
    )
    (tmp_path / "namelist.f90").write_text(
        """
module MOD_Namelist
 use MOD_Precision
 implicit none
 ! Real default from share/MOD_Namelist.F90; the sublimation skin mass
 ! changes the branch this test is here to check, so it is not arbitrary.
 real(r8) :: DEF_TRACER_SUBL_SKIN_MM = 5.0_r8
end module
""",
        encoding="utf-8",
    )
    (tmp_path / "forcing.f90").write_text(
        """
module MOD_Tracer_Forcing
 use MOD_Precision
 implicit none
contains
 real(r8) function tracer_forcing_precip_value(i,p); integer,intent(in)::i,p; tracer_forcing_precip_value=0._r8; end
 real(r8) function tracer_forcing_vapor_value(i,p); integer,intent(in)::i,p; tracer_forcing_vapor_value=0._r8; end
end module
""",
        encoding="utf-8",
    )
    (tmp_path / "frac.f90").write_text(
        """
module MOD_Tracer_Frac
 use MOD_Precision
 implicit none
contains
 logical function tracer_fractionation_active(i); integer,intent(in)::i; tracer_fractionation_active=.false.; end
 real(r8) function tracer_surface_relhum(a,b,c,d); real(r8),intent(in)::a,b,c; logical,intent(in)::d; tracer_surface_relhum=0._r8; end
 real(r8) function tracer_alpha_kinetic_craig_gordon(a,b); integer,intent(in)::a; logical,intent(in)::b; tracer_alpha_kinetic_craig_gordon=1._r8; end
 real(r8) function tracer_alpha_kinetic_open_water(a,b); integer,intent(in)::a; real(r8),intent(in)::b; tracer_alpha_kinetic_open_water=1._r8; end
 real(r8) function tracer_craig_gordon_evap_ratio(a,b,c,d,e,f,g); integer,intent(in)::a; real(r8),intent(in)::b,c,d,e,f; logical,intent(in)::g; tracer_craig_gordon_evap_ratio=b; end
 real(r8) function tracer_equilibrium_deposition_ratio(a,b,c,d); integer,intent(in)::a; real(r8),intent(in)::b,c; logical,intent(in)::d; tracer_equilibrium_deposition_ratio=b; end
end module
""",
        encoding="utf-8",
    )
    (tmp_path / "cons.f90").write_text(
        """
module MOD_Tracer_Conservation
 use MOD_Precision
 implicit none
contains
 subroutine tracer_save_storage(a,b,c); integer,intent(in)::a,b,c; end
 subroutine tracer_balance_check(ipatch,maxsnl,nl_soil,deltim,xerr_tracer,patchtype_in,water_err_in,water_dS_in,water_input_in,water_output_in,water_evap_in,water_rnof_in)
 integer,intent(in)::ipatch,maxsnl,nl_soil,patchtype_in; real(r8),intent(in)::deltim,water_err_in,water_dS_in,water_input_in,water_output_in,water_evap_in,water_rnof_in; real(r8),intent(out)::xerr_tracer; xerr_tracer=0._r8; end
 subroutine tracer_apply_reactive_processes(a,b,c,d); integer,intent(in)::a,b,c; real(r8),intent(in)::d; end
end module
""",
        encoding="utf-8",
    )
    (tmp_path / "hist.f90").write_text(
        """
module MOD_Tracer_Hist
 use MOD_Precision
 implicit none
contains
 subroutine tracer_hist_accumulate(ipatch,snl,maxsnl,nl_soil,a,b,c,d,e,f,g,h)
 integer,intent(in)::ipatch,snl,maxsnl,nl_soil; real(r8),intent(in)::a,b,c(:),d(:),e,f,g,h; end
end module
""",
        encoding="utf-8",
    )
    (tmp_path / "vars.f90").write_text(
        """
module MOD_Tracer_Vars
 use MOD_Precision
 implicit none
 integer,parameter::TRC_EVAP_KIND_SOILEVAP=1,TRC_EVAP_KIND_SUBL=2
 real(r8),allocatable::trc_wliq_soisno(:,:,:),trc_wice_soisno(:,:,:),trc_scv(:,:),trc_wdsrf(:,:),trc_ldew_rain(:,:),trc_ldew_snow(:,:),trc_rnof_step(:,:),a_trc_precip(:,:),a_trc_rsur(:,:),a_trc_rnof(:,:),trc_wetwat(:,:),trc_waterstorage(:,:),trc_storage_beg(:,:),trc_surface_residue(:,:),trc_subsurface_residue(:,:),trc_solid_soisno(:,:,:),trc_canopy_solid(:,:),trc_surface_solid(:,:),trc_subsurface_solid(:,:),trc_waterstorage_solid(:,:)
 logical,allocatable::trc_runtime_forced(:)
 real(r8)::booked(3,2)
contains
 subroutine tracer_book_evap_loss(i,p,m,w,k); integer,intent(in)::i,p,k; real(r8),intent(in)::m,w; booked(i,k)=booked(i,k)+m; end
 subroutine sync_tracer_patch_ratio(i,p,snl,maxsnl,nl_soil,wliq,wice,wa,wdsrf,scv,r); integer,intent(in)::i,p,snl,maxsnl,nl_soil; real(r8),intent(in)::wliq(:),wice(:),wa,wdsrf,scv,r; end
end module
""",
        encoding="utf-8",
    )
    (tmp_path / "driver.f90").write_text(
        """
program special_patch_driver
 use MOD_Precision
 use MOD_Tracer_Vars
 use MOD_Tracer_SpecialPatches
 implicit none
 real(r8) :: wliq(1:1), wice(1:1)
 allocate(trc_wliq_soisno(3,1,1),trc_wice_soisno(3,1,1),trc_scv(3,1),trc_wdsrf(3,1),trc_ldew_rain(3,1),trc_ldew_snow(3,1),trc_rnof_step(3,1),a_trc_precip(3,1),a_trc_rsur(3,1),a_trc_rnof(3,1),trc_wetwat(3,1),trc_storage_beg(3,1),trc_surface_residue(3,1),trc_subsurface_residue(3,1),trc_solid_soisno(3,1,1),trc_canopy_solid(3,1),trc_surface_solid(3,1),trc_subsurface_solid(3,1),trc_waterstorage_solid(3,1),trc_runtime_forced(3))
 trc_wliq_soisno=0; trc_wice_soisno=0; trc_scv=0; trc_wdsrf=0; trc_ldew_rain=0; trc_ldew_snow=0; trc_rnof_step=0; a_trc_precip=0; a_trc_rsur=0; a_trc_rnof=0; trc_wetwat=0; trc_storage_beg=0; trc_surface_residue=0; trc_subsurface_residue=0; trc_solid_soisno=0; trc_canopy_solid=0; trc_surface_solid=0; trc_subsurface_solid=0; trc_waterstorage_solid=0; trc_runtime_forced=.false.; booked=0
 trc_storage_beg(1,1)=4._r8
 wliq=0; wice=0
 call tracer_glacier_patch(1,0,1,1._r8, &
   0._r8,0._r8,0._r8,0._r8, &
   0._r8,6._r8,6._r8,0._r8,0._r8, &
   0._r8,0._r8,0._r8,0._r8, &
   10._r8,0._r8,273._r8,0._r8,100000._r8,wliq,wice)
 print '(6(ES20.12,1X))', booked(1,1), booked(1,2), booked(2,1), booked(2,2), booked(3,1), booked(3,2)
end program special_patch_driver
""",
        encoding="utf-8",
    )
    exe = tmp_path / "special_patch_driver"
    compiled = subprocess.run(
        [
            compiler,
            "-cpp",
            "-ffree-line-length-0",
            "-I",
            str(tmp_path),
            *(str(tmp_path / name) for name in ("precision.f90", "defs.f90", "namelist.f90", "forcing.f90", "frac.f90", "cons.f90", "hist.f90", "vars.f90")),
            # Real source, not a stub: the finite-pool limiter is the logic
            # under test here, so stubbing it would leave the assertions
            # checking the stub.
            str(ROOT / "main/TRACER/MOD_Tracer_EvapLimit.F90"),
            str(ROOT / "main/TRACER/MOD_Tracer_SpecialPatches.F90"),
            str(tmp_path / "driver.f90"),
            "-o",
            str(exe),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )
    assert compiled.returncode == 0, compiled.stdout + compiled.stderr
    result = subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True, timeout=SUBPROCESS_TIMEOUT)
    assert result.returncode == 0, result.stdout + result.stderr
    liq, ice, nonvol_liq, nonvol_ice, fixed_liq, fixed_ice = map(float, result.stdout.split())
    # The two branches are SEQUENTIAL, not an even split of the cap. Tracer 1
    # is well mixed over 10 mm carrying 4 units, so the 6 mm of liquid
    # evaporation takes 4 * 6/10 = 2.4; that leaves 1.6 units in 4 mm, and the
    # 6 mm sublimation demand exceeds what is left, so it takes all 1.6.
    assert liq == 2.4
    assert ice == 1.6
    # The cap: total booked loss can never exceed the storage it came from.
    assert liq + ice == 4.0
    # A non-volatile solute does not leave with the vapour at all.
    assert nonvol_liq == 0.0
    assert nonvol_ice == 0.0
    # A fixed-signature tracer loses ratio * water on each branch and does not
    # re-derive a ratio in between, so equal water losses give equal amounts:
    # 0.25 * 6 = 1.5 twice. This is what makes the 2.4/1.6 above a consequence
    # of the well-mixed re-derivation rather than of the branch order alone.
    assert fixed_liq == 1.5
    assert fixed_ice == 1.5
