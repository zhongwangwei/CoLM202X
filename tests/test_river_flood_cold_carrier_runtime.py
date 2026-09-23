"""Compile the production pre-publication carrier repair, not a Python mirror."""
from pathlib import Path
import subprocess

from fortran_test_support import require_runnable_fortran_compiler

FLOW = Path(__file__).resolve().parents[1] / 'main/HYDRO/MOD_Grid_RiverLakeFlow.F90'


def test_first_flood_publication_materializes_wet_restart_placeholder(tmp_path):
    source = FLOW.read_text().split('SUBROUTINE grid_riverlake_flow_init (', 1)[1]
    # Compile the entire production boundary, including its guards. The old
    # location remains supported so the recorded pre-fix failure is runnable.
    guard = '      IF (DEF_GridRiverLake_FloodFeedback .and. p_is_worker) THEN'
    before_seed = source.split('         trc_restart_found = .false.', 1)[0]
    if guard in before_seed:
        block = guard + before_seed.split(guard, 1)[1].split('\n      ENDIF', 1)[0] + '\n      ENDIF'
    else:
        block = source.split('            flood_infil_period = 0._r8', 1)[1]
        block = block.split('         ENDIF\n      ENDIF', 1)[0]
    code = '''
module fixture
 implicit none
 integer, parameter :: r8=kind(1.d0), numucat=3
 real(r8), parameter :: RIVERMIN=1.e-5_r8, spval=-1.e36_r8
 type curve
 contains
  procedure :: volume
 end type
 type(curve) :: floodplain_curve(numucat)
 logical :: volwater_ucat_valid, DEF_USE_LEVEE, has_levee(numucat)
 logical :: DEF_GridRiverLake_FloodFeedback=.true., p_is_worker=.true.
 real(r8) :: volwater_ucat(numucat),wdsrf_ucat(numucat),levsto(numucat)
 real(r8), allocatable :: volresv(:)
 integer, allocatable :: ucat2resv(:), dam_build_year(:)
 integer :: lake_type(numucat)=0, start_year=2003
 contains
 real(r8) function volume(this,stage)
  class(curve) :: this
  real(r8), intent(in) :: stage
  volume=10._r8*stage
 end function
 real(r8) function levee_visible_volume_from_stage(i,stage,protected)
  integer, intent(in) :: i
  real(r8), intent(in) :: stage,protected
  levee_visible_volume_from_stage=max(0._r8,10._r8*stage-protected)
 end function
 subroutine materialize()
  integer :: i, irsv
''' + block + '''
 end subroutine
end module
program check
 use fixture
 implicit none
 integer :: levee_case,valid_case
 real(r8) :: expected(3),seeded_mass(3)
 do levee_case=0,1
  do valid_case=0,1
   DEF_USE_LEVEE=levee_case==1
   has_levee=[DEF_USE_LEVEE,.false.,.false.]
   volwater_ucat_valid=valid_case==1
   wdsrf_ucat=[1._r8,0._r8,1._r8]
   levsto=[2._r8,0._r8,0._r8]
   volwater_ucat=[0._r8,0._r8,7._r8]
   expected=[10._r8,0._r8,10._r8]
   if(DEF_USE_LEVEE) expected(1)=8._r8
   if(volwater_ucat_valid) expected(3)=7._r8
   ! Cold isotope seeding already recovers this wet-stage water. The first
   ! flood publication must expose the SAME carrier without deleting mass.
   seeded_mass=0.002_r8*expected
   call materialize()
   if(.not.volwater_ucat_valid) error stop 'carrier remains invalid'
   if(maxval(abs(volwater_ucat-expected))>1.e-13_r8) &
    error stop 'first publication lost wet carrier or replaced persisted volume'
   if(maxval(abs(seeded_mass-0.002_r8*volwater_ucat))>1.e-15_r8) &
    error stop 'uniform isotope amount no longer matches first published carrier'
  enddo
 enddo
 ! Reservoir publication reads its separate pool: recover only the missing
 ! sentinel for dams already built, never replace persisted or future pools.
 allocate(volresv(3), ucat2resv(3), dam_build_year(3))
 DEF_USE_LEVEE=.false.; has_levee=.false.
 volwater_ucat_valid=.true.; volwater_ucat=7._r8
 lake_type=2; ucat2resv=[1,2,3]; dam_build_year=[2000,2000,2100]
 wdsrf_ucat=[1._r8,2._r8,3._r8]; volresv=[spval,17._r8,spval]
 call materialize()
 if(volresv(1)/=10._r8) error stop 'built reservoir lost stage-derived carrier'
 if(volresv(2)/=17._r8) error stop 'persisted reservoir overwritten'
 if(volresv(3)/=spval) error stop 'future reservoir initialized prematurely'
 lake_type(1)=0; volresv(1)=spval
 call materialize()
 if(volresv(1)/=spval) error stop 'ordinary cell wrote reservoir pool'
 lake_type(1)=2; wdsrf_ucat(1)=0._r8
 call materialize()
 if(volresv(1)/=0._r8) error stop 'dry built reservoir gained water'
 print '(A)', 'COLD_FLOOD_CARRIER_PASS'
end program
'''
    path = tmp_path / 'carrier.f90'
    path.write_text(code)
    compiler = require_runnable_fortran_compiler(tmp_path)
    built = subprocess.run(
        [compiler, '-ffree-line-length-0', '-fcheck=all',
         '-ffpe-trap=invalid,zero,overflow', str(path), '-o', 'carrier'],
        cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run([str(tmp_path / 'carrier')], capture_output=True, text=True, timeout=30)
    assert ran.returncode == 0, ran.stdout + ran.stderr
    assert 'COLD_FLOOD_CARRIER_PASS' in ran.stdout
