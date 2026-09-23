"""Restart initialization must not add an extra wetland mixing operation."""
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_restart_preserves_pools_but_cold_start_enforces_solubility(tmp_path):
    """Execute the actual initialization tail for hot/cold and reservoir cases."""
    compiler = shutil.which('gfortran')
    if not compiler:
        pytest.skip('gfortran unavailable')
    source = (ROOT / 'main/TRACER/MOD_Tracer_LandPhase.F90').read_text()
    tail = source.split('IF (present(file_restart)) CALL tracer_lifecycle_land_read_restart(file_restart)', 1)[1]
    tail = tail.split('END SUBROUTINE tracer_init_from_arrays', 1)[0]
    program = '''program probe
implicit none
integer :: calls=0, mode
logical :: found_restart
real(8) :: pools(3), saved(3), water(1)
water=1
! Deliberately distinct surface/wetland/solid pools (not equilibrium).
saved=[0.2d0,0.7d0,0.1d0]
do mode=1,2
 found_restart=(mode==1)
 pools=saved
 call init()
 if (found_restart .and. any(pools/=saved)) stop 1
 if (.not.found_restart .and. calls/=1) stop 2
 pools=saved
 call init(water)
 if (found_restart .and. any(pools/=saved)) stop 3
 if (.not.found_restart .and. calls/=2) stop 4
enddo
contains
subroutine init(waterstorage)
real(8), optional :: waterstorage(1)
integer :: numpatch=1,maxsnl=0,nl_soil=1
real(8) :: ldew_rain(1)=0,wliq_soisno(1,1)=0,wa(1)=0,wdsrf(1)=1,wetwat(1)=2
''' + tail + '''end subroutine
subroutine tracer_enforce_solubility_from_water(n,m,l,rain,soil,aquifer,surface,wetland,reservoir)
integer :: n,m,l
real(8) :: rain(n),soil(m+1:l,n),aquifer(n),surface(n),wetland(n)
real(8),optional :: reservoir(n)
calls=calls+1
pools=sum(pools)/3d0
end subroutine
end program
'''
    path = tmp_path / 'probe.f90'
    path.write_text(program)
    exe = tmp_path / 'probe'
    subprocess.run([compiler, '-fcheck=all', '-ffree-line-length-none', str(path), '-o', str(exe)], check=True, capture_output=True)
    subprocess.run([str(exe)], check=True, capture_output=True)


def test_restart_compatibility_still_validates_all_phase_pools():
    source = (ROOT / 'main/TRACER/MOD_Tracer_Rest.F90').read_text()
    reader = source.split('SUBROUTINE read_land_tracer_restart', 1)[1].split('END SUBROUTINE read_land_tracer_restart', 1)[0]
    assert 'IF (.not. descriptor_matches) THEN' in reader
    assert 'whole-domain cold start' in reader
    assert "CALL CoLM_stop('incomplete or malformed committed generic land tracer restart')" in reader
    for pool in ('trc_wdsrf', 'trc_wetwat', 'trc_surface_solid', 'trc_solid_soisno'):
        assert f"tracer_dim_matches(file_restart, '{pool}'" in reader
    assert reader.index('CALL validate_land_tracer_restart_state(wa)') < reader.rindex('found_restart = .true.')
