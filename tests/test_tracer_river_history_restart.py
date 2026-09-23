from pathlib import Path
import subprocess
import re

from fortran_test_support import require_runnable_fortran_compiler


SOURCE = (Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_RiverLake.F90").read_text()


def routine(name):
    start = re.search(rf"(?m)^\s*SUBROUTINE {name}\s*\(", SOURCE)
    assert start, name
    return SOURCE[start.end() - 1:].split(f"END SUBROUTINE {name}", 1)[0]


def test_mid_window_restart_keeps_numerator_and_own_elapsed_time():
    accumulate = routine("tracer_diag_accumulate_substep")
    write = routine("write_tracer_restart")
    read = routine("read_tracer_restart")
    history = routine("write_tracer_history")
    flush = routine("tracer_flush_acc")

    # The clock advances under the same ucatfilter, valid-system and positive-dt
    # guards as the tracer history numerators, including dry-cell flux records.
    assert accumulate.index("IF (.not. ucatfilter(i)) CYCLE") < accumulate.index(
        "a_trc_acctime(i) = a_trc_acctime(i) + dt_i"
    ) < accumulate.index("a_trc_out   (itrc, i) =")
    assert "CALL write_tracer_hist_acc (file_restart, 'trc_hist_acctime', acc1=a_trc_acctime)" in write
    assert "CALL read_tracer_hist_acc (file_restart, 'trc_hist_acctime', acc1=a_trc_acctime)" in read
    assert "IF (.not. hist_complete) THEN" in read
    assert "a_trc_acctime = 0._r8" in read
    assert "a_trc_acctime = 0._r8" in flush
    assert "a_trc_out(itrc, i) / a_trc_acctime(i)" in history
    assert "a_trc_storage_mass(itrc, i) / a_water_storage(i)" in history
    assert "a_water_storage(i) <= trc_delta_diag_vmin * a_trc_acctime(i)" in history



def test_reservoir_volume_uses_actual_built_mask_during_transport():
    volume = routine("get_cell_volume")
    substep = routine("tracer_substep")
    diagnostics = routine("tracer_diag_accumulate_substep")
    cold_start = routine("tracer_init_from_water")
    assert "IF (present(is_built_resv_cell)) use_reservoir = is_built_resv_cell" in volume
    assert "IF (use_reservoir .and. size(volresv_in) > 0" in volume
    assert "volwater_ucat(icell)" in volume  # unbuilt cells use routed river water
    assert "volwater, is_built_resv(i))" in substep
    assert "volwater, is_built_resv(i))" in diagnostics
    assert cold_start.count("volwater, is_built_resv(i))") == 2


def test_built_and_unbuilt_reservoir_volume_and_seeding_compiled(tmp_path):
    compiler = require_runnable_fortran_compiler(tmp_path)
    get_volume = "SUBROUTINE get_cell_volume" + routine("get_cell_volume") + "END SUBROUTINE get_cell_volume"
    cold_start = "SUBROUTINE tracer_init_from_water" + routine("tracer_init_from_water") + "END SUBROUTINE tracer_init_from_water"
    source = f"""
module MOD_Precision
 integer, parameter :: r8 = kind(1.d0)
end module
module MOD_Namelist
 logical :: DEF_USE_LEVEE = .false.
end module
module MOD_Vars_Global
 use MOD_Precision
 real(r8) :: spval = -9999._r8
end module
module MOD_Grid_RiverLakeNetwork
 use MOD_Precision
 integer :: numucat = 1
 type curve_type
 contains
   procedure :: volume
 end type
 type(curve_type), allocatable :: floodplain_curve(:)
 integer, allocatable :: lake_type(:)
contains
 real(r8) function volume(this, stage)
   class(curve_type), intent(in) :: this
   real(r8), intent(in) :: stage
   volume = 10._r8 * stage
 end function
end module
module MOD_Grid_RiverLakeLevee
 use MOD_Precision
 logical, allocatable :: has_levee(:)
 real(r8), allocatable :: levsto(:)
contains
 real(r8) function levee_visible_volume_from_stage(icell, stage, protected)
   integer, intent(in) :: icell
   real(r8), intent(in) :: stage, protected
   levee_visible_volume_from_stage = 10._r8 * stage - protected
 end function
end module
module MOD_Grid_RiverLakeTimeVars
 use MOD_Precision
 logical :: volwater_ucat_valid = .true.
 real(r8), allocatable :: volwater_ucat(:)
end module
module MOD_Tracer_Defs
 use MOD_Precision
 integer :: ntracers = 1
contains
 logical function tracer_uses_land_water_transport(itrc)
   integer, intent(in) :: itrc
   tracer_uses_land_water_transport = .true.
 end function
 logical function tracer_has_dissolved_limit(itrc)
   integer, intent(in) :: itrc
   tracer_has_dissolved_limit = .false.
 end function
 subroutine tracer_equilibrate_dissolved(itrc, water_mass, dissolved_mass, solid_mass)
   integer, intent(in) :: itrc
   real(r8), intent(in) :: water_mass
   real(r8), intent(inout) :: dissolved_mass, solid_mass
 end subroutine
 real(r8) function tracer_init_water_ratio(itrc)
   integer, intent(in) :: itrc
   tracer_init_water_ratio = 2._r8
 end function
end module
module actual_tracer_volume
 use MOD_Precision
 use MOD_Namelist
 use MOD_Tracer_Defs, only: ntracers, tracer_uses_land_water_transport, &
   tracer_has_dissolved_limit, tracer_equilibrate_dissolved
 implicit none
 logical :: p_is_worker = .true., p_is_master = .false.
 real(r8), allocatable :: trc_mass(:,:), trc_conc(:,:), trc_levsto(:,:)
 real(r8), allocatable :: trc_solid(:,:), trc_levsto_solid(:,:)
contains
{get_volume}
{cold_start}
subroutine update_tracer_concentration(itrc, icell, volwater)
 integer, intent(in) :: itrc, icell
 real(r8), intent(in) :: volwater
 trc_conc(itrc, icell) = trc_mass(itrc, icell) / volwater
end subroutine
end module
program check
 use MOD_Precision
 use MOD_Grid_RiverLakeNetwork
 use MOD_Grid_RiverLakeTimeVars
 use actual_tracer_volume
 implicit none
 real(r8) :: v
 allocate(lake_type(1), floodplain_curve(1), volwater_ucat(1))
 allocate(trc_mass(1,1), trc_conc(1,1), trc_levsto(1,1))
 lake_type = 2
 volwater_ucat = 7._r8
 call get_cell_volume(1, 2._r8, [50._r8], [1], v, .false.)
 if (v /= 7._r8) stop 1
 call get_cell_volume(1, 2._r8, [50._r8], [1], v, .true.)
 if (v /= 50._r8) stop 2
 call tracer_init_from_water([2._r8], [50._r8], [1], is_built_resv=[.false.])
 if (trc_mass(1,1) /= 14._r8 .or. trc_conc(1,1) /= 2._r8) stop 3
 call tracer_init_from_water([2._r8], [50._r8], [1], is_built_resv=[.true.])
 if (trc_mass(1,1) /= 100._r8 .or. trc_conc(1,1) /= 2._r8) stop 4
end program
"""
    path = tmp_path / "volume.f90"
    path.write_text(source)
    result = subprocess.run(
        [compiler, "-ffree-line-length-none", str(path), "-o", str(tmp_path / "volume")],
        cwd=tmp_path, text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    result = subprocess.run([str(tmp_path / "volume")], text=True, capture_output=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_production_restart_to_history_roundtrip_compiled(tmp_path):
    """Real MPI/NetCDF integration; rebuild the tested module, never fake I/O."""
    import os
    import shutil
    import re

    import pytest
    from fortran_test_support import netcdf_fortran_flags

    root = Path(__file__).resolve().parents[1]
    build = Path(os.environ.get("COLM_BLD_DIR", root / ".bld")).resolve()
    compiler, launcher = shutil.which("mpif90"), shutil.which("mpiexec")
    if not compiler or not launcher or not (build / "mod_tracer_riverlake.mod").exists():
        pytest.skip("requires MPI and a TRACER-enabled production build in COLM_BLD_DIR/.bld")
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1",
               OMPI_MCA_rmaps_base_oversubscribe="1")
    smoke = subprocess.run([launcher, "-n", "1", "/usr/bin/true"], env=env,
                           capture_output=True, text=True, timeout=30)
    if smoke.returncode:
        pytest.skip("MPI launcher cannot start processes")
    includes, libs = netcdf_fortran_flags()
    flags = ["-fopenmp", "-fdefault-real-8", "-ffree-line-length-none", "-cpp",
             "-fallow-argument-mismatch", "-fcheck=all", f"-I{tmp_path}",
             f"-I{build}", f"-I{root / 'include'}", f"-J{tmp_path}", *includes]

    def command(args):
        result = subprocess.run(args, cwd=tmp_path, text=True, capture_output=True, timeout=180)
        assert result.returncode == 0, result.stdout + result.stderr

    # Archive only in the temporary directory; shared build artifacts are read-only.
    archive = tmp_path / "model.a"
    objects = sorted(str(p) for p in build.glob("*.o")
                     if p.name not in {"MOD_Tracer_RiverLake.o", "MOD_Tracer_Lifecycle_Registrations_Stubs.o"})
    command(["ar", "rcs", str(archive), *objects])
    harness = root / "tests/tracer_river_history_restart_harness.F90"
    executable = tmp_path / "history_check"
    platform_libs = ["-framework", "Accelerate"] if os.uname().sysname == "Darwin" else ["-llapack", "-lblas"]

    def compile_source(source):
        path = tmp_path / "MOD_Tracer_RiverLake.F90"
        path.write_text(source)
        command([compiler, *flags, "-c", str(path), "-o", str(tmp_path / "tracer.o")])
        command([compiler, *flags, str(harness), str(tmp_path / "tracer.o"), str(archive),
                 *libs, *platform_libs, "-o", str(executable)])

    def run(mode, directory):
        directory.mkdir(exist_ok=True)
        return subprocess.run([launcher, "-n", "3", str(executable), mode], cwd=directory,
                              env=env, text=True, capture_output=True, timeout=60)

    compile_source(SOURCE)
    for mode in ("continuous", "split", "legacy", "partial", "clockless"):
        result = run(mode, tmp_path / mode)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "TRACER_HISTORY_RESTART_OK" in result.stdout

    result = run("finite_v2", tmp_path / "finite_v2")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "FINITE_RESTART_OK" in result.stdout

    # Generate committed old/new transactions with the production writer,
    # changing only the selected solid rows. Read both with the new reader.
    solid_write = """            write(varname, '(A,A)') 'trc_solid_', trim(tracer_names(itrc))
            CALL vector_gather_and_write ( &
               tmpvec, numucat, totalnumucat, ucat_data_address, file_restart, trim(varname), 'ucatch')
"""
    protected_write = solid_write.replace("'trc_solid_'", "'trc_levsto_solid_'")
    assert solid_write in SOURCE and protected_write in SOURCE
    old_writer = SOURCE.replace(
        "RIVER_TRACER_RESTART_SCHEMA_VERSION = 2",
        "RIVER_TRACER_RESTART_SCHEMA_VERSION = 1", 1,
    ).replace(solid_write, "", 1).replace(protected_write, "", 1)
    compile_source(old_writer)
    result = run("finite_write_only", tmp_path / "finite_v1")
    assert result.returncode == 0, result.stdout + result.stderr
    compile_source(SOURCE)
    result = run("finite_read_only", tmp_path / "finite_v1")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "FINITE_RESTART_OK" in result.stdout

    compile_source(SOURCE.replace(solid_write, "", 1))
    result = run("finite_write_only", tmp_path / "finite_missing")
    assert result.returncode == 0, result.stdout + result.stderr
    compile_source(SOURCE)
    result = run("finite_read_only", tmp_path / "finite_missing")
    assert result.returncode != 0
    assert "incomplete or malformed committed river tracer restart" in result.stdout + result.stderr

    # Negative control: reproduce the old missing-history writer, while keeping
    # valid prognostic state and descriptor metadata. It must not pass split.
    mutant, removed = re.subn(r"(?im)^.*CALL write_tracer_hist_acc \([^\n]+\n", "", SOURCE)
    assert removed >= 7
    compile_source(mutant)
    result = run("split", tmp_path / "missing_history_writer")
    assert result.returncode != 0
    assert "restart history mismatch" in result.stdout + result.stderr
