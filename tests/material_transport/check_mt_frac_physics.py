#!/usr/bin/env python3
"""Focused regression checks for MOD_MT_Frac physics (F05 / F06 / F07).

The routines under test are extracted verbatim from
``main/MATERIAL/FAMILY_ISOTOPE/MOD_MT_Frac.F90`` (source anchors kept as
``routine()`` matches) and compiled with gfortran against a minimal stub host.
No model data, MPI or NetCDF is needed.

Covered invariants
------------------
F05  NSS leaf-water history term uses the stored OLD BULK delta:
     * the first-step increment must be independent of ``prev_delta_e``
       (old code used it, so this fails on the pre-fix source);
     * without NSS memory the result must be independent of ``prev_delta_b``;
     * with memory, changing ``prev_delta_b`` must change the result;
     * constant forcing iterates to a finite steady state whose
       transpiration ratio approaches the source ratio.
F06  water-vapour equivalent diffusion is isotope-only:
     * isotope -> positive diffusivity; nonvolatile solute -> exact 0;
     * the same holds for both the soil and the firn/snow helpers;
     * no air porosity (saturated, frozen) -> exact 0 for both families.
F07  frozen soil ice volume uses the ice density, not the liquid density:
     * a column whose liquid + ice exactly fill the pores has zero air
       porosity and therefore exactly zero vapour diffusivity (the old
       ``/1000`` conversion leaves ~2.9% air porosity here);
     * unfrozen and zero-liquid end members keep a positive pathway.
"""
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

root = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else \
    pathlib.Path(__file__).resolve().parents[2]
frac_path = root / 'main/MATERIAL/FAMILY_ISOTOPE/MOD_MT_Frac.F90'
source = frac_path.read_text()


def routine(name):
    pattern = (r'^ *[^\n]*\b(SUBROUTINE|FUNCTION) +' + name +
               r'\s*\([^\n]*(?:\n.*?)*?^ *END +(SUBROUTINE|FUNCTION) +' + name + r'\b[^\n]*')
    match = re.search(pattern, source, re.I | re.M)
    assert match, f'production routine not found: {name}'
    return match.group(0)


fixture = r'''
module fixture
use, intrinsic :: iso_fortran_env, only: real64
implicit none
integer, parameter :: r8 = real64
integer, parameter :: ntracers = 3
real(r8), parameter :: trc_tiny = 1.e-30_r8
real(r8), parameter :: tfrz = 273.15_r8
real(r8), parameter :: water_moles_per_mm = 1000._r8 / 18.01528_r8
real(r8), parameter :: liquid_water_molar_density = 55.5e3_r8
real(r8), parameter :: universal_gas_constant = 8.31446261815324_r8
real(r8), parameter :: DEF_ISO_NSS_LEAF_WATER_PER_LAI = 0.11_r8
real(r8), parameter :: DEF_ISO_NSS_LEAF_PATH_LENGTH = 0.02_r8
real(r8), parameter :: DEF_ISO_NSS_LEAF_RB = 120._r8
real(r8), parameter :: DEF_ISO_CG_RELHUM_MAX = 0.99_r8
real(r8), parameter :: ref_ratio_o18 = 0.0020052_r8
logical :: iso_family(ntracers) = [.true., .false., .false.]
contains
logical function tracer_fractionation_active(itrc)
  integer, intent(in) :: itrc
  tracer_fractionation_active = .true.
end function
logical function tracer_is_isotope(itrc)
  integer, intent(in) :: itrc
  tracer_is_isotope = iso_family(itrc)
end function
real(r8) function tracer_ratio_to_delta(itrc, ratio)
  integer, intent(in) :: itrc
  real(r8), intent(in) :: ratio
  tracer_ratio_to_delta = (ratio / ref_ratio_o18 - 1._r8) * 1000._r8
end function
real(r8) function tracer_delta_to_ratio(itrc, delta)
  integer, intent(in) :: itrc
  real(r8), intent(in) :: delta
  tracer_delta_to_ratio = ref_ratio_o18 * (1._r8 + delta / 1000._r8)
end function
real(r8) function tracer_alpha_liq_vap(itrc, temp_k)
  integer, intent(in) :: itrc
  real(r8), intent(in) :: temp_k
  tracer_alpha_liq_vap = 1.0098_r8
end function
real(r8) function tracer_alpha_ice_vap(itrc, temp_k)
  integer, intent(in) :: itrc
  real(r8), intent(in) :: temp_k
  tracer_alpha_ice_vap = 1.0117_r8
end function
real(r8) function tracer_diffusivity_ratio_air(itrc)
  integer, intent(in) :: itrc
  tracer_diffusivity_ratio_air = 1.0251_r8
end function
real(r8) function tracer_alpha_kinetic_leaf(itrc, ra, rb, rc)
  integer, intent(in) :: itrc
  real(r8), intent(in) :: ra, rb, rc
  tracer_alpha_kinetic_leaf = 1.0085_r8
end function
real(r8) function tracer_leaf_liquid_diffusivity(itrc, temp_k)
  integer, intent(in) :: itrc
  real(r8), intent(in) :: temp_k
  tracer_leaf_liquid_diffusivity = 2.0e-9_r8
end function
'''

extracted = [routine('tracer_saturation_vapor_pressure'),
             routine('tracer_soil_vapor_equivalent_diffusivity'),
             routine('tracer_isotope_soil_vapor_equivalent_diffusivity'),
             routine('tracer_snow_vapor_equivalent_diffusivity'),
             routine('tracer_isotope_snow_vapor_equivalent_diffusivity'),
             routine('tracer_transpiration_nss_ratio')]

program = r'''
end module fixture

program check_frac_physics
use fixture
implicit none
integer :: nfail, npass
real(r8), parameter :: psrf = 101325._r8
real(r8), parameter :: dz = 0.1_r8, porsl = 0.45_r8
real(r8), parameter :: rho_water = 1000._r8, rho_ice = 917._r8
real(r8) :: d_iso, d_solute, d_sat, d_unf, d_dry, d_ref
real(r8) :: ratio, src, vap
real(r8) :: de, db, pec, lm, tr
real(r8) :: de1, db1, pec1, lm1, tr1
real(r8) :: de2, db2, pec2, lm2, tr2
real(r8) :: de3, db3, pec3, lm3, tr3
real(r8) :: sat_ice_mm
integer :: k
logical :: ok

nfail = 0
npass = 0

! ---------------- F07: ice volume uses ice density ----------------
! liquid + ice exactly fill the pores: theta_l=0.1, theta_ice=0.35, air=0
sat_ice_mm = (porsl - 10._r8 / (rho_water * dz)) * rho_ice * dz
d_sat = tracer_soil_vapor_equivalent_diffusivity(10._r8, sat_ice_mm, dz, porsl, &
   268.15_r8, psrf, 1.0251_r8, 1.0117_r8)
ok = (d_sat == 0._r8)
call report('F07 saturated frozen soil has zero air porosity -> D_vap == 0', ok, nfail, npass)
if (.not. ok) write(*,*) '  sat_ice_mm=', sat_ice_mm, ' D=', d_sat

d_unf = tracer_soil_vapor_equivalent_diffusivity(10._r8, 0._r8, dz, porsl, &
   298.15_r8, psrf, 1.0251_r8, 1.0098_r8)
call report('F07 unfrozen partially wet soil keeps positive D_vap', d_unf > 0._r8, nfail, npass)

d_dry = tracer_soil_vapor_equivalent_diffusivity(0._r8, 0._r8, dz, porsl, &
   268.15_r8, psrf, 1.0251_r8, 1.0117_r8)
call report('F07 zero liquid/ice keeps positive D_vap', d_dry > 0._r8, nfail, npass)

! 10 mm liquid + 32.095 mm ice at dz=0.1 m, porsl=0.45 fills pores exactly
! (10/1000/0.1=0.10, 32.095/917/0.1=0.35).  The old /1000 conversion gives
! air porosity 0.02905 and a nonzero diffusivity.
d_ref = tracer_soil_vapor_equivalent_diffusivity(10._r8, 32.095_r8, dz, porsl, &
   268.15_r8, psrf, 1.0251_r8, 1.0117_r8)
call report('F07 ice density conversion: pore-filling column gives D_vap == 0', &
   d_ref == 0._r8, nfail, npass)

! ---------------- F06: vapour diffusion is isotope-only ----------------
d_iso = tracer_isotope_soil_vapor_equivalent_diffusivity(1, 10._r8, 0._r8, dz, porsl, &
   298.15_r8, psrf)
d_solute = tracer_isotope_soil_vapor_equivalent_diffusivity(2, 10._r8, 0._r8, dz, porsl, &
   298.15_r8, psrf)
call report('F06 soil vapour: isotope positive', d_iso > 0._r8, nfail, npass)
call report('F06 soil vapour: nonvolatile solute is exactly zero', d_solute == 0._r8, nfail, npass)

d_iso = tracer_isotope_snow_vapor_equivalent_diffusivity(1, 0._r8, 30._r8, 0.1_r8, &
   268.15_r8, psrf)
d_solute = tracer_isotope_snow_vapor_equivalent_diffusivity(2, 0._r8, 30._r8, 0.1_r8, &
   268.15_r8, psrf)
call report('F06 firn vapour: isotope positive', d_iso > 0._r8, nfail, npass)
call report('F06 firn vapour: nonvolatile solute is exactly zero', d_solute == 0._r8, nfail, npass)

d_iso = tracer_isotope_soil_vapor_equivalent_diffusivity(1, 10._r8, sat_ice_mm, dz, porsl, &
   268.15_r8, psrf)
d_solute = tracer_isotope_soil_vapor_equivalent_diffusivity(2, 10._r8, sat_ice_mm, dz, porsl, &
   268.15_r8, psrf)
call report('F06 saturated frozen soil: isotope zero', d_iso == 0._r8, nfail, npass)
call report('F06 saturated frozen soil: solute zero', d_solute == 0._r8, nfail, npass)

! ---------------- F05: NSS history term uses the old BULK delta ----------------
! Constant leaf geometry and forcing; old source -35 permil, old evaporation
! -5 permil, old bulk -20 permil; new source -50 permil.
ratio = tracer_delta_to_ratio(1, -20._r8)
src   = tracer_delta_to_ratio(1, -50._r8)
call tracer_transpiration_nss_ratio(1, src, ratio, 298.15_r8, 0.6_r8, psrf, &
   0.2_r8, 1800._r8, 2.0_r8, 30._r8, 200._r8, &
   -5._r8, -20._r8, 0.5_r8, 100._r8, tr, de, db, pec, lm)
call check1_invariants('F05 baseline', de, db, pec, lm, tr, nfail, npass)

! Same input but a wildly different OLD EVAPORATION delta: the first step must
! not depend on it.  This is the invariant the pre-fix source violated.
call tracer_transpiration_nss_ratio(1, src, ratio, 298.15_r8, 0.6_r8, psrf, &
   0.2_r8, 1800._r8, 2.0_r8, 30._r8, 200._r8, &
   -80._r8, -20._r8, 0.5_r8, 100._r8, tr2, de2, db2, pec2, lm2)
ok = same4(de, db, pec, tr, de2, db2, pec2, tr2)
call report('F05 first step is independent of prev_delta_e', ok, nfail, npass)

! Same input but a different OLD BULK delta must change the result.
call tracer_transpiration_nss_ratio(1, src, ratio, 298.15_r8, 0.6_r8, psrf, &
   0.2_r8, 1800._r8, 2.0_r8, 30._r8, 200._r8, &
   -5._r8, -10._r8, 0.5_r8, 100._r8, tr3, de3, db3, pec3, lm3)
ok = .not. same4(de, db, pec, tr, de3, db3, pec3, tr3)
call report('F05 with memory, prev_delta_b changes the first step', ok, nfail, npass)

! No NSS memory: prev_delta_b must be ignored (initialised to source).
call tracer_transpiration_nss_ratio(1, src, ratio, 298.15_r8, 0.6_r8, psrf, &
   0.2_r8, 1800._r8, 2.0_r8, 30._r8, 200._r8, &
   -5._r8, -20._r8, 0.5_r8, 0._r8, tr, de, db, pec, lm)
call tracer_transpiration_nss_ratio(1, src, ratio, 298.15_r8, 0.6_r8, psrf, &
   0.2_r8, 1800._r8, 2.0_r8, 30._r8, 200._r8, &
   -5._r8, -70._r8, 0.5_r8, 0._r8, tr2, de2, db2, pec2, lm2)
ok = same4(de, db, pec, tr, de2, db2, pec2, tr2)
call report('F05 without memory, prev_delta_b is ignored', ok, nfail, npass)

! Constant forcing must iterate to a finite steady state; at steady state the
! leaf storage tendency vanishes and the transpiration ratio approaches the
! source ratio.
vap = tracer_delta_to_ratio(1, -12._r8)
call tracer_transpiration_nss_ratio(1, src, vap, 298.15_r8, 0.6_r8, psrf, &
   0.2_r8, 1800._r8, 2.0_r8, 30._r8, 200._r8, &
   -20._r8, -20._r8, 0.5_r8, 100._r8, tr, de, db, pec, lm)
do k = 1, 400
   de1 = de; db1 = db; pec1 = pec; lm1 = lm
   call tracer_transpiration_nss_ratio(1, src, vap, 298.15_r8, 0.6_r8, psrf, &
      0.2_r8, 1800._r8, 2.0_r8, 30._r8, 200._r8, &
      de1, db1, pec1, lm1, tr, de, db, pec, lm)
end do
ok = (de == de) .and. (db == db) .and. (abs(de) < 1.e6_r8) .and. (abs(db) < 1.e6_r8) .and. &
     (abs(tr / src - 1._r8) < 0.05_r8)
call report('F05 constant forcing converges; steady transpiration ~= source', ok, nfail, npass)
if (.not. ok) write(*,*) '  converged: de=', de, ' db=', db, ' tr/src=', tr / src

write(*,'(A,I0,A,I0,A)') 'PASS: ', npass, ' checks; FAIL: ', nfail, ' checks'
if (nfail > 0) stop 1

contains

subroutine report(label, ok, nfail, npass)
character(len=*), intent(in) :: label
logical, intent(in) :: ok
integer, intent(inout) :: nfail, npass
if (ok) then
   write(*,'(A,A)') 'PASS: ', label
   npass = npass + 1
else
   write(*,'(A,A)') 'FAIL: ', label
   nfail = nfail + 1
end if
end subroutine

logical function same4(a1,a2,a3,a4,b1,b2,b3,b4)
real(r8), intent(in) :: a1,a2,a3,a4,b1,b2,b3,b4
same4 = (a1 == b1) .and. (a2 == b2) .and. (a3 == b3) .and. (a4 == b4)
end function

subroutine check1_invariants(label, de, db, pec, lm, tr, nfail, npass)
character(len=*), intent(in) :: label
real(r8), intent(in) :: de, db, pec, lm, tr
integer, intent(inout) :: nfail, npass
logical :: ok
ok = (de == de) .and. (db == db) .and. (tr == tr) .and. &
     (lm > 0._r8) .and. (pec >= 0._r8) .and. (pec <= 1._r8)
call report(trim(label)//' finite outputs and physical bounds', ok, nfail, npass)
end subroutine

end program check_frac_physics
'''

gfortran = shutil.which('gfortran')
if not gfortran:
    print('SKIP: gfortran not available; this check needs a Fortran compiler')
    sys.exit(0)

with tempfile.TemporaryDirectory() as tmp:
    tmp = pathlib.Path(tmp)
    src = tmp / 'check_frac_physics.F90'
    exe = tmp / 'check_frac_physics'
    src.write_text(fixture + '\n' +
                   '\n'.join(extracted) + '\n' +
                   program.lstrip('\n'))
    compile_cmd = [gfortran, '-std=f2008', '-O0', '-ffree-line-length-none',
                   '-fcheck=bounds', '-o', str(exe), str(src)]
    proc = subprocess.run(compile_cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print('FAIL: compile of extracted production routines failed')
        print(proc.stdout)
        print(proc.stderr)
        sys.exit(1)
    run = subprocess.run([str(exe)], capture_output=True, text=True)
    print(run.stdout, end='')
    if run.returncode != 0:
        print(run.stderr, end='')
        sys.exit(1)
