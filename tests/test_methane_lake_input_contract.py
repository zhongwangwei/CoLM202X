from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STATE = (ROOT / "main/TRACER/MOD_Tracer_Reactive_Methane_State.F90").read_text()
INITIALIZE = (ROOT / "mkinidata/MOD_Initialize.F90").read_text()
AGGREGATION = (ROOT / "mksrfdata/Aggregation_LakeSoilC.F90").read_text()


def _subroutine(source: str, name: str) -> str:
    start = source.index(f"SUBROUTINE {name}")
    end = source.index(f"END SUBROUTINE {name}", start)
    return source[start:end]


def test_lulcc_remap_initializes_new_lakes_from_current_surface_data() -> None:
    remap_lulcc = _subroutine(STATE, "remap_methane_lulcc_state")
    no_snapshot = remap_lulcc.index("IF (.not. methane_lulcc_snapshot_valid) THEN")
    build_map = remap_lulcc.index("CALL build_lulcc_remap_map")
    assert "initialize_methane_lake_soilc_from_surface" in remap_lulcc
    assert "lake_soilc_srf" in remap_lulcc
    assert "initialize_lake_from_surface" in remap_lulcc
    assert no_snapshot < build_map


def test_enabled_lake_production_rejects_any_missing_lake_patch() -> None:
    initialize = _subroutine(STATE, "initialize_methane_lake_soilc_from_surface")
    assert "IF (lake_soilc_nlake == 0) RETURN" in initialize
    assert "lake_soilc_missing > 0" in initialize
    assert "CALL CoLM_stop" in initialize
    assert "WARNING: lake CH4 production is enabled" not in initialize


def test_single_point_lake_uses_existing_site_organic_matter_input() -> None:
    assert "lake_soilc_srf(:,ipatch)" in INITIALIZE
    assert "OM_density(:,ipatch)" in INITIALIZE


def test_gridded_lake_carbon_prefers_raw_area_weights_then_explicit_om_proxy() -> None:
    assert AGGREGATION.index("IF (lake_global > 0) THEN") < AGGREGATION.index(
        "IF (raw_exists) THEN"
    )
    assert AGGREGATION.index("IF (raw_exists) THEN") < AGGREGATION.index(
        "ncio_read_vector_complete"
    )
    assert "no lake patches in domain; writing zero lake sediment carbon" in AGGREGATION
    assert "area = area_one" in AGGREGATION
    assert "median(lake_soilc_one" not in AGGREGATION
    assert "raw lake_soilc.nc absent; using organic-matter proxy, not measured lake carbon" in AGGREGATION
    assert "n_source_layers = min(8, max(1, nl_soil-1))" in AGGREGATION
    assert "min(8, max(1, j-1)) == src_layer" in AGGREGATION
    assert "CALL ncio_read_vector_complete" in AGGREGATION
    assert "defval" not in AGGREGATION.split("CALL ncio_read_vector_complete", 1)[1].split(")", 1)[0]
    assert "ieee_is_finite(vf_om_value)" in AGGREGATION
    assert "vf_om_value < 0._r8 .or. vf_om_value > 1._r8" in AGGREGATION
    assert "CALL mpi_allreduce (invalid_local, invalid_global" in AGGREGATION
    assert "CALL CoLM_stop" in AGGREGATION


def test_lake_om_proxy_formula_maps_eight_source_layers_to_ten_model_layers() -> None:
    import re
    import pytest

    carbon = float(re.search(r"carbon_per_kg_om = ([0-9.]+)_r8", AGGREGATION).group(1))
    organic = float(re.search(r"organic_max_default = ([0-9.]+)_r8", AGGREGATION).group(1))
    assert carbon == 580.0 and organic == 130.0
    source_om = [0.01 * layer for layer in range(1, 9)]
    mapped = [carbon * organic * source_om[min(8, max(1, j - 1)) - 1]
              for j in range(1, 11)]
    assert mapped == pytest.approx([754.0, 754.0, 1508.0, 2262.0, 3016.0,
                                    3770.0, 4524.0, 5278.0, 6032.0, 6032.0])


def test_lake_carbon_initialization_runtime(tmp_path):
    """Execute production init selection/helper with absent surface data and LULCC masks."""
    import subprocess
    from fortran_test_support import require_runnable_fortran_compiler

    provider = (ROOT / 'main/TRACER/MOD_Tracer_Reactive_Methane.F90').read_text()
    selection = provider.split('      lake_restart_present = .false.', 1)[1].split(
        '      CALL allocate_methane_giems', 1)[0]
    selection = '      lake_restart_present = .false.' + selection
    helper = _subroutine(STATE, 'initialize_methane_lake_soilc_from_surface')
    helper += 'END SUBROUTINE initialize_methane_lake_soilc_from_surface\n'
    source = '''module MOD_SPMD_Task
contains
subroutine CoLM_stop(message)
 character(len=*),intent(in)::message
 print *,message
 error stop 1
end subroutine
end module
module probe
use, intrinsic :: ieee_arithmetic
implicit none
integer, parameter :: r8=kind(1.d0), nl_soil=2, PATCHTYPE_LAKE=4
real(r8), allocatable :: lake_soilc(:,:), lake_soilc_srf(:,:)
integer, allocatable :: patchtype(:)
logical :: p_is_worker=.true.
integer :: landpatch=0
type options
 logical :: allowlakeprod=.true.
end type
type(options) :: DEF_METHANE
contains
subroutine CoLM_stop(message)
 character(len=*),intent(in)::message
 print *,message
 error stop 1
end subroutine
 elemental logical function invalid_restart_value(x)
 real(r8),intent(in)::x
 invalid_restart_value = .not. ieee_is_finite(x) .or. abs(x)>1.d35
 end function
subroutine ncio_vector_group_presence(file, fields, patch, found)
 character(len=*),intent(in)::file,fields(:)
 integer,intent(in)::patch
 logical,intent(out)::found(:)
 found = file == 'valid' .or. file == 'zero'
end subroutine
subroutine initialize(file_restart)
 character(len=*),intent(in),optional::file_restart
 logical::lake_restart_present(2)
 character(len=32),parameter::lake_restart_fields(2) = &
 [character(len=32)::'ch4_conc_methane','ch4_lake_soilc']
''' + selection + '\nend subroutine\n' + helper + '''
end module
program main
use probe
implicit none
character(len=20)::mode
call get_command_argument(1,mode)
allocate(lake_soilc(2,2),patchtype(2))
lake_soilc=0
patchtype=PATCHTYPE_LAKE
select case(trim(mode))
case('valid','zero')
 if(trim(mode)=='zero') then
  allocate(lake_soilc_srf(2,2))
  lake_soilc_srf=5.d0
 endif
 call initialize(trim(mode)) ! Missing surface must not pre-empt restart loading.
 if(any(lake_soilc/=0)) error stop 2
case('cold')
 call initialize()
case('legacy')
 call initialize('legacy')
case('new_lake')
 allocate(lake_soilc_srf(2,2))
 lake_soilc_srf(:,1)=-1.d36
 lake_soilc_srf(:,2)=3.d0
 call initialize_methane_lake_soilc_from_surface(patchtype,lake_soilc_srf,.true.,[.false.,.true.])
 if(any(lake_soilc(:,1)/=0)) error stop 3
 if(any(lake_soilc(:,2)/=3)) error stop 4
case('old_lakes')
 call initialize_methane_lake_soilc_from_surface(patchtype,lake_soilc_srf,.true.,[.false.,.false.])
 if(any(lake_soilc/=0)) error stop 5
case('new_missing')
 call initialize_methane_lake_soilc_from_surface(patchtype,lake_soilc_srf,.true.,[.false.,.true.])
end select
end program
'''
    path = tmp_path / 'probe.f90'
    path.write_text(source)
    executable = tmp_path / 'probe'
    compiler = require_runnable_fortran_compiler(tmp_path)
    subprocess.run([compiler, '-ffree-line-length-none', '-fcheck=all', str(path),
                    '-o', str(executable)], cwd=tmp_path, check=True, capture_output=True)
    for mode in ('valid', 'zero', 'new_lake', 'old_lakes'):
        result = subprocess.run([str(executable), mode], capture_output=True, text=True)
        assert result.returncode == 0, (mode, result.stdout, result.stderr)
    for mode in ('cold', 'legacy', 'new_missing'):
        result = subprocess.run([str(executable), mode], capture_output=True, text=True)
        assert result.returncode != 0, mode
        assert 'requires lake_soilc surface data' in result.stdout


def test_lake_restart_probe_does_not_bypass_transaction_or_inventory_read():
    provider = (ROOT / 'main/TRACER/MOD_Tracer_Reactive_Methane.F90').read_text()
    reader = _subroutine(provider, 'ch4_reactive_read_restart')
    assert reader.index('CALL validate_methane_restart_transaction') < reader.index('CALL read_methane_restart')
    init = _subroutine(provider, 'ch4_reactive_init')
    assert "'ch4_conc_methane', 'ch4_lake_soilc'" in init
    assert init.index('CALL ncio_vector_group_presence') < init.index('IF (p_is_worker')
    assert 'file_restart)' in (ROOT / 'main/TRACER/MOD_Tracer_LandPhase.F90').read_text().split(
        'CALL tracer_lifecycle_land_init', 1)[1].split('\n', 1)[0]
