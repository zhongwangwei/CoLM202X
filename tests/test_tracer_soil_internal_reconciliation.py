"""Execute production closed-column remapping and its legacy open-column fallback."""
from pathlib import Path
import os
import subprocess

import pytest
from fortran_test_support import require_runnable_fortran_compiler

ROOT = Path(__file__).resolve().parents[1]


def routine(source, kind, name):
    start = source.index(f"{kind} {name}")
    end = source.index(f"END {kind} {name}", start) + len(f"END {kind} {name}")
    return source[start:end]


CASES = [
    ("tiny_positive_donor_no_ratio_overflow", [1.e-310,0,0,0], [5.e-311,5.e-311,0,0], [1,0,0,0], [0]*4, [0.5,0.5,0,0], [0]*4, 0, False),
    ("forward_dry_passthrough", [4,0,0,0], [1,1,1,1], [8,0,0,0], [0]*4, [2]*4, [0]*4, 0, False),
    ("reverse_dry_passthrough", [0,0,0,4], [1,1,1,1], [0,0,0,8], [0]*4, [2]*4, [0]*4, 0, False),
    ("opposing_edges_current_donor", [0,4,0,4], [2,0,4,2], [0,4,0,12], [0]*4, [2,0,8,6], [0]*4, 0, False),
    ("gradient_updates_each_donor", [4,1,1,0], [1,1,1,3], [8,10,20,0], [0]*4, [2,4,8,24], [0]*4, 0, False),
    ("finite_solid_stays_at_origin", [4,0,0,0], [1,1,1,1], [4,0,0,0], [4,0,0,0], [1]*4, [4,0,0,0], 0, True),
    ("finite_dry_receiver_dissolves_locally", [4,0,0,0], [1,1,1,1], [0]*4, [0,5,0,0], [0,1,1,1], [0,2,0,0], 0, True),
    ("isotope_uniform_ratio", [4,0,0,0], [1,1,1,1], [0.008,0,0,0], [0]*4, [0.002]*4, [0]*4, 0, False),
    ("roundoff_closed_column", [4,0,0,0], [1,1,1,1+1.e-14], [8,0,0,0], [0]*4, [2]*4, [0]*4, 0, False),
    ("small_open_residual_not_internal_flow", [1,2,3,4], [1+1.e-8,2,3,4], [2,6,12,20], [0]*4, [2+2.e-8,6,12,20], [0]*4, 2.e-8, False),
    ("absent_mask_keeps_fallback", [4,0,0,0], [1,1,1,1], [8,0,0,0], [0]*4, [2,0,0,0], [0]*4, -6, False),
    ("disconnected_mask_no_cross_gap", [4,0,0,0], [1,1,1,1], [8,0,0,0], [0]*4, [2,0,0,0], [0]*4, -6, False),
    ("disconnected_closed_block_transports", [4,0,0,0], [2,2,0,0], [8,0,0,0], [0]*4, [4,4,0,0], [0]*4, 0, False),
    ("open_gain_keeps_numerical_ledger", [1,2,3,4], [2,2,3,4], [2,6,12,20], [0]*4, [4,6,12,20], [0]*4, 2, False),
    ("open_loss_keeps_numerical_ledger", [1,2,3,4], [0,2,3,4], [1,6,12,20], [0]*4, [0,6,12,20], [0]*4, -1, False),
]


def vector(values):
    return '[' + ','.join(f'{float(x)!r}_r8' for x in values) + ']'


@pytest.mark.parametrize('name,water,target,mass,solid,expected,expected_solid,ledger,finite', CASES, ids=[x[0] for x in CASES])
def test_closed_remap_and_open_fallback(tmp_path, name, water, target, mass, solid, expected, expected_solid, ledger, finite):
    source = Path(os.environ.get('COLM_SOIL_REMAP_SOURCE', ROOT / 'main/TRACER/MOD_Tracer_SoilWater.F90')).read_text()
    defs = (ROOT / 'main/TRACER/MOD_Tracer_Defs.F90').read_text()
    marker = 'SUBROUTINE reconcile_internal_soil_flow'
    helper = routine(source, 'SUBROUTINE', 'reconcile_internal_soil_flow') if marker in source else ''
    # Same behavioral test executes the original fallback for a baseline RED,
    # rather than failing because the new helper's name is absent.
    call = 'call reconcile_internal_soil_flow()' if helper else ''
    if helper:
        helper += '\n' + routine(source, 'SUBROUTINE', 'move_dissolved_face')
        assert source.index('CALL reconcile_internal_soil_flow()') < source.index('water_resid = wliq_soisno(j) - water_shadow(j)')
    start = source.rindex('            DO j = 1, nl_soil', 0, source.index('water_resid = wliq_soisno(j) - water_shadow(j)'))
    end = source.index('! The VSF layer residual', start)
    fallback = source[start:end]
    if name == 'tiny_positive_donor_no_ratio_overflow':
        # Exercise the actual face helper: an all-tiny column is intentionally
        # below the reconciliation eligibility threshold. M/W would overflow.
        assert helper, 'Tiny-donor test requires the production face helper'
        call = 'call move_dissolved_face(1, 2, 5.e-311_r8)'
        fallback = ''
    equilibrium = routine(defs, 'SUBROUTINE', 'tracer_equilibrate_dissolved')
    invocation = 'call run()' if name.startswith('absent_mask') else ('call run([.true.,.true.,.false.,.true.])' if name.startswith('disconnected') else 'call run([.true.,.true.,.true.,.true.])')
    program = f'''
module fixture
contains
subroutine run(permeable_soil)
 use, intrinsic :: ieee_arithmetic
 implicit none
 logical, optional, intent(in) :: permeable_soil(4)
 integer, parameter :: r8=kind(1.d0), nl_soil=4
 real(r8), parameter :: trc_tiny=1.e-30_r8, trc_water_min_for_ratio=1.e-12_r8
 integer :: itrc=1, ipatch=1, j
 real(r8) :: water_shadow(4),wliq_soisno(4),trc_wliq_soisno(1,4,1),trc_solid_soisno(1,4,1)
 real(r8) :: remap_face_water(nl_soil-1),remap_trial_water(nl_soil)
 real(r8) :: water_resid,water_shadow_ratio,trc_flux,soil_resid_trc,total_before
 type species
 real(r8) :: max_dissolved_conc=1._r8
 end type
 type(species) :: tracers(1)
 water_shadow={vector(water)}
 wliq_soisno={vector(target)}
 trc_wliq_soisno(1,:,1)={vector(mass)}
 trc_solid_soisno(1,:,1)={vector(solid)}
 total_before=sum(trc_wliq_soisno)+sum(trc_solid_soisno)
 soil_resid_trc=0._r8
 {call}
 {fallback}
 if(.not.all(ieee_is_finite(trc_wliq_soisno))) error stop 8
 if(.not.all(ieee_is_finite(trc_solid_soisno))) error stop 9
 if(.not.all(ieee_is_finite(water_shadow))) error stop 10
 if(.not.ieee_is_finite(soil_resid_trc)) error stop 11
 if(any(abs(trc_wliq_soisno(1,:,1)-{vector(expected)})>1.e-12_r8)) then
 print *, 'unexpected dissolved mass',trc_wliq_soisno
 error stop 1
 endif
 if(any(abs(trc_solid_soisno(1,:,1)-{vector(expected_solid)})>1.e-12_r8)) error stop 2
 if(abs(soil_resid_trc-({float(ledger)!r}_r8))>1.e-12_r8) error stop 3
 if(any(trc_wliq_soisno<0._r8).or.any(trc_solid_soisno<0._r8)) error stop 4
 if(abs(sum(trc_wliq_soisno)+sum(trc_solid_soisno)-total_before-soil_resid_trc)>1.e-12_r8) error stop 5
 if({'.true.' if ledger == 0 else '.false.'}) then
 if(any(abs(water_shadow-wliq_soisno)>1.e-12_r8)) error stop 6
 else
 if(any(water_shadow/={vector(water)})) error stop 7
 endif
contains
 {helper}
 {equilibrium}
 logical function tracer_has_dissolved_limit(index)
 integer,intent(in)::index
 tracer_has_dissolved_limit={'.true.' if finite else '.false.'}
 end function
 real(r8) function current_liq_ratio(index)
 integer,intent(in)::index
 current_liq_ratio=0._r8
 if(water_shadow(index)>trc_water_min_for_ratio) current_liq_ratio=trc_wliq_soisno(1,index,1)/water_shadow(index)
 end function
end subroutine
end module
program remap
 use fixture
 {invocation}
end program
'''
    compiler = require_runnable_fortran_compiler(tmp_path)
    src, exe = tmp_path / 'remap.f90', tmp_path / 'remap'
    src.write_text(program)
    built = subprocess.run([compiler,'-ffree-line-length-none','-fcheck=all','-ffpe-trap=invalid,zero,overflow',str(src),'-o',str(exe)], capture_output=True, text=True, timeout=60)
    assert built.returncode == 0, built.stdout + built.stderr
    ran = subprocess.run([str(exe)], capture_output=True, text=True, timeout=10)
    assert ran.returncode == 0, ran.stdout + ran.stderr
