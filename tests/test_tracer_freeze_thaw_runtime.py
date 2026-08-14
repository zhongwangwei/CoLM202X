"""Runtime check on the soil freeze/thaw tracer transfer in MOD_Tracer_Evapo.

The freeze step draws from the liquid pool that exists AFTER any thaw in the
same step, so its Rayleigh denominator must be wliq_soisno_bef + thaw_amt --
the same quantity wliq_post_phase is built from, and the same convention the
canopy path uses. Using the pre-thaw liquid overstates the freeze ratio by
(wliq_bef + thaw)/wliq_bef, and the numerator has already been incremented by
the thaw flux, so nothing downstream corrects it.

Reachable only when THERMAL reports thaw and freeze for the same layer in one
step, which the explicit soil_thaw_mass_th / soil_frzc_mass_th path allows.
The driver therefore uses that path rather than the sign-pattern heuristic,
whose two branches are mutually exclusive and cannot produce the case.

MOD_Tracer_Evapo and MOD_Tracer_EvapLimit are compiled from real source. The
Rayleigh helper is stubbed with the exact formula the real one uses when
fractionation is inactive, because what is under test is the denominator
Evapo passes in, not the Rayleigh math itself.
"""

from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
SUBPROCESS_TIMEOUT = 60


def _write_stubs(tmp_path: Path) -> None:
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
 integer :: ntracers=1
 real(r8),parameter :: trc_tiny=1e-12_r8, trc_delta_sanity_max=2.0e3_r8
 type :: tracer_info_type
   real(r8) :: ref_ratio = 1.0_r8
 end type tracer_info_type
 type(tracer_info_type) :: tracers(1)
contains
 logical function tracer_uses_land_water_transport(i); integer,intent(in)::i; tracer_uses_land_water_transport=.true.; end
 logical function tracer_is_nonvolatile_solute(i); integer,intent(in)::i; tracer_is_nonvolatile_solute=.false.; end
 real(r8) function tracer_init_water_ratio(i); integer,intent(in)::i; tracer_init_water_ratio=0._r8; end
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
 real(r8) :: DEF_TRACER_SUBL_SKIN_MM = 5.0_r8
 real(r8) :: DEF_TRACER_CANOPY_EQUILIBRATION = 0.0_r8
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
 real(r8) function tracer_forcing_vapor_value(i,p); integer,intent(in)::i,p; tracer_forcing_vapor_value=0._r8; end
end module
""",
        encoding="utf-8",
    )
    # tracer_rayleigh_freezing_loss reproduces the real function's
    # fractionation-inactive branch verbatim:
    #   min(freeze_water * pool_trc/pool_water, pool_trc)
    # including its pool_water/freeze_water guards, so the value it returns is
    # decided entirely by the denominator Evapo hands it.
    (tmp_path / "frac.f90").write_text(
        """
module MOD_Tracer_Frac
 use MOD_Precision
 implicit none
 real(r8),parameter :: eps=1e-12_r8
contains
 logical function tracer_fractionation_active(i); integer,intent(in)::i; tracer_fractionation_active=.false.; end
 real(r8) function tracer_rayleigh_freezing_loss(i,pool_trc,pool_water,freeze_water,temp_k)
  integer,intent(in)::i; real(r8),intent(in)::pool_trc,pool_water,freeze_water,temp_k
  tracer_rayleigh_freezing_loss=0._r8
  if (pool_water<=eps .or. freeze_water<=eps) return
  if (pool_trc<=eps) return
  if (freeze_water>=pool_water*(1._r8-1.e-12_r8)) then
    tracer_rayleigh_freezing_loss=max(pool_trc,0._r8); return
  endif
  tracer_rayleigh_freezing_loss=min(freeze_water*(max(pool_trc,0._r8)/pool_water),max(pool_trc,0._r8))
 end function
 real(r8) function tracer_surface_relhum(a,b,c,d); real(r8),intent(in)::a,b,c; logical,intent(in)::d; tracer_surface_relhum=0._r8; end
 real(r8) function tracer_alpha_kinetic_craig_gordon(a,b); integer,intent(in)::a; logical,intent(in)::b; tracer_alpha_kinetic_craig_gordon=1._r8; end
 real(r8) function tracer_alpha_liq_vap(a,b); integer,intent(in)::a; real(r8),intent(in)::b; tracer_alpha_liq_vap=1._r8; end
 real(r8) function tracer_craig_gordon_evap_ratio(a,b,c,d,e,f,g); integer,intent(in)::a; real(r8),intent(in)::b,c,d,e,f; logical,intent(in)::g; tracer_craig_gordon_evap_ratio=b; end
 real(r8) function tracer_equilibrium_deposition_ratio(a,b,c,d); integer,intent(in)::a; real(r8),intent(in)::b,c; logical,intent(in)::d; tracer_equilibrium_deposition_ratio=b; end
 real(r8) function tracer_equilibration_exchange(pool_trc,pool_water,vapor_ratio,alpha_eq,equilibration_fraction)
  real(r8),intent(in)::pool_trc,pool_water,vapor_ratio,alpha_eq,equilibration_fraction
  tracer_equilibration_exchange=0._r8
 end function
end module
""",
        encoding="utf-8",
    )
    (tmp_path / "vars.f90").write_text(
        """
module MOD_Tracer_Vars
 use MOD_Precision
 implicit none
 integer,parameter::TRC_EVAP_KIND_CANOPYEVAP=1,TRC_EVAP_KIND_SOILEVAP=2,TRC_EVAP_KIND_SUBL=3
 real(r8),allocatable::trc_wliq_soisno(:,:,:),trc_wice_soisno(:,:,:),trc_solid_soisno(:,:,:)
 real(r8),allocatable::trc_ldew_rain(:,:),trc_ldew_snow(:,:),trc_canopy_solid(:,:)
 real(r8),allocatable::a_trc_precip(:,:),a_trc_vapor_exchange(:,:),trc_numerical_residual_step(:,:)
contains
 subroutine tracer_book_evap_loss(i,p,m,w,k); integer,intent(in)::i,p,k; real(r8),intent(in)::m,w; end
end module
""",
        encoding="utf-8",
    )
    # One soil layer, no snow. Pre-THERMAL: 2 mm liquid holding 2 units of
    # tracer, 8 mm ice holding 8 -- both pools at ratio 1.0. THERMAL reports
    # the whole 8 mm thawing and 1 mm re-freezing in the same step, so the
    # liquid pool the freeze draws from is 2 + 8 = 10 mm holding 10 units,
    # still at ratio 1.0. Freezing 1 mm must therefore carry exactly 1 unit.
    # Against the pre-thaw 2 mm denominator it carries 10 * 1/2 = 5.
    (tmp_path / "driver.f90").write_text(
        """
program freeze_thaw_driver
 use MOD_Precision
 use MOD_Tracer_Vars
 use MOD_Tracer_Evapo
 implicit none
 real(r8) :: wliq(1:1), wice(1:1), wliq_bef(1:1), wice_bef(1:1)
 real(r8) :: thaw(1:1), frzc(1:1), tfrac(1:1)
 allocate(trc_wliq_soisno(1,1,1),trc_wice_soisno(1,1,1),trc_solid_soisno(1,1,1))
 allocate(trc_ldew_rain(1,1),trc_ldew_snow(1,1),trc_canopy_solid(1,1))
 allocate(a_trc_precip(1,1),a_trc_vapor_exchange(1,1),trc_numerical_residual_step(1,1))
 trc_ldew_rain=0; trc_ldew_snow=0; trc_canopy_solid=0
 a_trc_precip=0; a_trc_vapor_exchange=0; trc_numerical_residual_step=0
 trc_solid_soisno=0
 trc_wliq_soisno(1,1,1)=2._r8
 trc_wice_soisno(1,1,1)=8._r8
 wliq_bef=2._r8; wice_bef=8._r8
 thaw=8._r8; frzc=1._r8
 wliq=2._r8+8._r8-1._r8      ! 9
 wice=8._r8-8._r8+1._r8      ! 1
 tfrac=273._r8
 call tracer_evapo(1, 1._r8, 0, 1, &
   0._r8,0._r8,0._r8,0._r8, &
   wliq, wice, wliq_bef, wice_bef, &
   soil_thaw_mass_th=thaw, soil_frzc_mass_th=frzc, t_soisno_frac=tfrac)
 print '(2(ES20.12,1X))', trc_wliq_soisno(1,1,1), trc_wice_soisno(1,1,1)
end program freeze_thaw_driver
""",
        encoding="utf-8",
    )


def test_soil_freeze_uses_the_post_thaw_liquid_pool(tmp_path: Path) -> None:
    compiler = require_runnable_fortran_compiler(tmp_path)
    _write_stubs(tmp_path)
    exe = tmp_path / "freeze_thaw_driver"
    compiled = subprocess.run(
        [
            compiler,
            "-cpp",
            "-ffree-line-length-0",
            "-I",
            str(tmp_path),
            *(
                str(tmp_path / name)
                for name in (
                    "precision.f90",
                    "defs.f90",
                    "namelist.f90",
                    "forcing.f90",
                    "frac.f90",
                    "vars.f90",
                )
            ),
            str(ROOT / "main/TRACER/MOD_Tracer_EvapLimit.F90"),
            str(ROOT / "main/TRACER/MOD_Tracer_Evapo.F90"),
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

    result = subprocess.run(
        [str(exe)], cwd=tmp_path, capture_output=True, text=True, timeout=SUBPROCESS_TIMEOUT
    )
    assert result.returncode == 0, result.stdout + result.stderr
    liq, ice = map(float, result.stdout.split())

    # 10 units end up in the liquid pool after the thaw; freezing 1 mm out of
    # 10 mm at ratio 1.0 moves exactly 1 unit. The pre-thaw denominator gives
    # 5 instead, which would leave liq = 5.0 and ice = 5.0.
    assert abs(ice - 1.0) < 1e-12, f"freeze carried {ice}, expected 1.0"
    assert abs(liq - 9.0) < 1e-12, f"liquid left {liq}, expected 9.0"
    # Nothing is created or destroyed by an internal phase transfer.
    assert abs((liq + ice) - 10.0) < 1e-12
