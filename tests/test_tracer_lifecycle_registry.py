from pathlib import Path
import subprocess
import tempfile

import pytest

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
DEFS = ROOT / "main" / "TRACER" / "MOD_Tracer_Defs.F90"
LIFECYCLE = ROOT / "main" / "TRACER" / "MOD_Tracer_Lifecycle.F90"
SUBPROCESS_TIMEOUT = 30

pytestmark = pytest.mark.skipif(
    not LIFECYCLE.exists(),
    reason="unified tracer lifecycle module has not been introduced yet",
)


@pytest.fixture(scope="module")
def lifecycle_driver():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        compiler = require_runnable_fortran_compiler(tmp)
        (tmp / "define.h").write_text(
            "#define TRACER\n#define BGC\n", encoding="utf-8"
        )
        (tmp / "precision.f90").write_text(
            """
module MOD_Precision
  implicit none
  integer, parameter :: r8 = selected_real_kind(12)
end module MOD_Precision
""",
            encoding="utf-8",
        )
        (tmp / "namelist.f90").write_text(
            """
module MOD_Namelist
  use MOD_Precision, only: r8
  implicit none
  integer :: DEF_TRACER_NUM = 0
  character(len=256) :: DEF_TRACER_NAMES = ''
  character(len=256) :: DEF_TRACER_TYPES = ''
  character(len=256) :: DEF_TRACER_MRAT = ''
  character(len=256) :: DEF_TRACER_REF_RATIO = ''
  character(len=256) :: DEF_TRACER_INIT_DELTA = ''
  character(len=256) :: DEF_TRACER_REACTIVE_DECAY_RATE = ''
  character(len=512) :: DEF_TRACER_PARAM_FILES = 'null'
  logical :: DEF_TRACER_USE_FRACTIONATION = .false.
end module MOD_Namelist
""",
            encoding="utf-8",
        )
        (tmp / "support.f90").write_text(
            """
module MOD_SPMD_Task
  implicit none
  logical :: p_is_master = .true.
contains
  subroutine CoLM_stop(message)
    character(len=*), intent(in), optional :: message
    if (present(message)) write(*,'(A)') trim(message)
    error stop 1
  end subroutine CoLM_stop
end module MOD_SPMD_Task

module MOD_Vars_Global
  use MOD_Precision, only: r8
  implicit none
  real(r8), parameter :: spval = -1.0e36_r8
end module MOD_Vars_Global

module MOD_DataType
  implicit none
  type :: block_data_real8_2d
    integer :: unused = 0
  end type block_data_real8_2d
end module MOD_DataType

module lifecycle_test_control
  use MOD_Precision, only: r8
  use MOD_DataType, only: block_data_real8_2d
  implicit none
  character(len=32) :: scenario = ''
  integer :: ch4_itrc = -1
  integer :: sediment_itrc = -1
  integer :: land_final_calls = 0
  integer :: route_final_calls = 0
contains
  subroutine count_land_final()
    land_final_calls = land_final_calls + 1
  end subroutine count_land_final

  subroutine count_route_final()
    route_final_calls = route_final_calls + 1
  end subroutine count_route_final

  subroutine noop_noarg()
  end subroutine noop_noarg

  subroutine noop_land_init(numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)
    integer, intent(in) :: numpatch, lc_year, jdate(3)
    character(len=*), intent(in) :: casename, dir_restart, dir_landdata
  end subroutine noop_land_init

  subroutine noop_restart(file_restart)
    character(len=*), intent(in) :: file_restart
  end subroutine noop_restart

  subroutine noop_land_write_restart(file_restart, compress)
    character(len=*), intent(in) :: file_restart
    integer, intent(in) :: compress
  end subroutine noop_land_write_restart

  subroutine noop_land_history(file_hist, itime_in_file, sumarea, filter, &
    nl_soil, forcing_has_missing_value, forcmask_pch)
    character(len=*), intent(in) :: file_hist
    integer, intent(in) :: itime_in_file, nl_soil
    type(block_data_real8_2d), intent(inout) :: sumarea
    logical, intent(inout) :: filter(:)
    logical, intent(in) :: forcing_has_missing_value, forcmask_pch(:)
  end subroutine noop_land_history

  subroutine noop_soil_step(istep_local, ipatch, idate, deltim)
    integer, intent(in) :: istep_local, ipatch, idate(3)
    real(r8), intent(in) :: deltim
  end subroutine noop_soil_step

  subroutine noop_route_calc(deltime)
    real(r8), intent(in) :: deltime
  end subroutine noop_route_calc

  subroutine noop_route_history(file_hist_ucat, itime_in_file_ucat)
    character(len=*), intent(in) :: file_hist_ucat
    integer, intent(in) :: itime_in_file_ucat
  end subroutine noop_route_history
end module lifecycle_test_control
""",
            encoding="utf-8",
        )
        (tmp / "driver.f90").write_text(
            """
#include <define.h>
subroutine register_all_tracer_providers()
  use MOD_Tracer_Defs, only: FAMILY_GAS, FAMILY_PARTICLE, &
    STATE_OWNER_PROVIDER, REACTION_NONE, REACTION_PROVIDER
  use MOD_Tracer_Lifecycle, only: tracer_lifecycle_hooks_type, &
    register_tracer_provider
  use lifecycle_test_control
  implicit none
  type(tracer_lifecycle_hooks_type) :: hooks
  integer :: duplicate_itrc

  ch4_itrc = 0
  sediment_itrc = 0

  if (trim(scenario) /= 'unknown' .and. &
      trim(scenario) /= 'ch4_without_provider') then
    hooks = tracer_lifecycle_hooks_type()
    hooks%land_final => count_land_final
    if (trim(scenario) /= 'gas_incomplete') then
      hooks%land_init => noop_land_init
      hooks%land_read_restart => noop_restart
      hooks%land_write_restart => noop_land_write_restart
      hooks%land_history => noop_land_history
      if (trim(scenario) /= 'reaction_incomplete') hooks%soil_step => noop_soil_step
    endif
    call register_tracer_provider('CH4', 'METHANE', 'methane', &
      FAMILY_GAS, STATE_OWNER_PROVIDER, REACTION_PROVIDER, hooks, ch4_itrc)
  endif

  hooks = tracer_lifecycle_hooks_type()
  hooks%route_final => count_route_final
  if (trim(scenario) /= 'particle_incomplete') then
    hooks%route_init => noop_noarg
    hooks%route_read_restart => noop_restart
    hooks%route_write_restart => noop_restart
    hooks%route_calc => noop_route_calc
    hooks%route_history => noop_route_history
    hooks%route_flush_history => noop_noarg
  endif
  if (trim(scenario) == 'mismatch') then
    call register_tracer_provider('SEDIMENT', 'SED', 'sediment', &
      FAMILY_GAS, STATE_OWNER_PROVIDER, REACTION_NONE, hooks, sediment_itrc)
  else
    call register_tracer_provider('SEDIMENT', 'SED', 'sediment', &
      FAMILY_PARTICLE, STATE_OWNER_PROVIDER, REACTION_NONE, hooks, sediment_itrc)
  endif

  if (trim(scenario) == 'duplicate') then
    call register_tracer_provider('SED', 'SEDIMENT', 'duplicate-sediment', &
      FAMILY_PARTICLE, STATE_OWNER_PROVIDER, REACTION_NONE, hooks, duplicate_itrc)
  endif
end subroutine register_all_tracer_providers

program lifecycle_driver
  use MOD_Namelist
  use MOD_Tracer_Defs
  use MOD_Tracer_Lifecycle
  use lifecycle_test_control
  implicit none
  character(len=512) :: ch4_parameter_file

  call get_command_argument(1, scenario)
  call get_command_argument(2, ch4_parameter_file)
  call configure_descriptors(trim(scenario), trim(ch4_parameter_file))
  call tracer_defs_init()
  call tracer_lifecycle_init()
  call tracer_lifecycle_validate()

  select case (trim(scenario))
  case ('happy')
    write(*,'(A,I0)') 'SIZE=', tracer_lifecycle_size()
    write(*,'(A,2(I0,1X))') 'INDEX=', ch4_itrc, sediment_itrc
    write(*,'(A,3(L1,1X))') 'HAS=', tracer_lifecycle_has_provider(1), &
      tracer_lifecycle_has_provider(2), tracer_lifecycle_has_provider(3)
    write(*,'(A)') 'CH4_PROVIDER=' // trim(tracer_lifecycle_provider(ch4_itrc))
    write(*,'(A)') 'SED_PROVIDER=' // trim(tracer_lifecycle_provider(sediment_itrc))
    write(*,'(A,3(A,1X))') 'FAMILY=', trim(tracers(1)%category), &
      trim(tracers(ch4_itrc)%category), trim(tracers(sediment_itrc)%category)
    write(*,'(A,3(I0,1X))') 'SOLUTE_DESC=', tracers(1)%family_id, &
      tracers(1)%state_owner, tracers(1)%reaction_mode
    write(*,'(A,3(I0,1X))') 'CH4_DESC=', tracers(ch4_itrc)%family_id, &
      tracers(ch4_itrc)%state_owner, tracers(ch4_itrc)%reaction_mode
    write(*,'(A,3(I0,1X))') 'SED_DESC=', tracers(sediment_itrc)%family_id, &
      tracers(sediment_itrc)%state_owner, tracers(sediment_itrc)%reaction_mode
    call tracer_lifecycle_land_final()
    write(*,'(A,2(I0,1X))') 'AFTER_LAND=', land_final_calls, route_final_calls
    call tracer_lifecycle_route_final()
    write(*,'(A,2(I0,1X))') 'AFTER_ROUTE=', land_final_calls, route_final_calls
  case ('absent')
    write(*,'(A,I0)') 'SIZE=', tracer_lifecycle_size()
    write(*,'(A,2(I0,1X))') 'INDEX=', ch4_itrc, sediment_itrc
    write(*,'(A,L1)') 'HAS=', tracer_lifecycle_has_provider(1)
  case ('zero')
    call tracer_lifecycle_land_final()
    call tracer_lifecycle_route_final()
    write(*,'(A,I0)') 'SIZE=', tracer_lifecycle_size()
    write(*,'(A,2(I0,1X))') 'CALLS=', land_final_calls, route_final_calls
  end select

  call tracer_lifecycle_reset()
  call tracer_defs_final()

contains
  subroutine configure_descriptors(action, parameter_file)
    character(len=*), intent(in) :: action, parameter_file

    select case (trim(action))
    case ('happy')
      DEF_TRACER_NUM = 3
      DEF_TRACER_NAMES = 'CL,METHANE,SED'
      DEF_TRACER_TYPES = 'solute,gas,particle'
      DEF_TRACER_MRAT = '35.453,16.04,1.0'
      DEF_TRACER_REF_RATIO = '1.0,1.0,1.0'
      DEF_TRACER_INIT_DELTA = '0.0,0.0,0.0'
      DEF_TRACER_REACTIVE_DECAY_RATE = '0.0,0.0,0.0'
      DEF_TRACER_PARAM_FILES = 'METHANE:' // trim(parameter_file)
    case ('absent')
      DEF_TRACER_NUM = 1
      DEF_TRACER_NAMES = 'CL'
      DEF_TRACER_TYPES = 'solute'
      DEF_TRACER_MRAT = '35.453'
      DEF_TRACER_REF_RATIO = '1.0'
      DEF_TRACER_INIT_DELTA = '0.0'
      DEF_TRACER_REACTIVE_DECAY_RATE = '0.0'
    case ('zero')
      DEF_TRACER_NUM = 0
    case ('unknown')
      DEF_TRACER_NUM = 1
      DEF_TRACER_NAMES = 'DUST'
      DEF_TRACER_TYPES = 'particle'
      DEF_TRACER_MRAT = '1.0'
      DEF_TRACER_REF_RATIO = '1.0'
      DEF_TRACER_INIT_DELTA = '0.0'
      DEF_TRACER_REACTIVE_DECAY_RATE = '0.0'
    case ('duplicate', 'mismatch', 'particle_incomplete')
      DEF_TRACER_NUM = 1
      DEF_TRACER_NAMES = 'SED'
      DEF_TRACER_TYPES = 'particle'
      DEF_TRACER_MRAT = '1.0'
      DEF_TRACER_REF_RATIO = '1.0'
      DEF_TRACER_INIT_DELTA = '0.0'
      DEF_TRACER_REACTIVE_DECAY_RATE = '0.0'
    case ('gas_incomplete', 'reaction_incomplete')
      DEF_TRACER_NUM = 1
      DEF_TRACER_NAMES = 'CH4'
      DEF_TRACER_TYPES = 'gas'
      DEF_TRACER_MRAT = '16.04'
      DEF_TRACER_REF_RATIO = '1.0'
      DEF_TRACER_INIT_DELTA = '0.0'
      DEF_TRACER_REACTIVE_DECAY_RATE = '0.0'
    case ('ch4_without_provider')
      DEF_TRACER_NUM = 1
      DEF_TRACER_NAMES = 'CH4'
      DEF_TRACER_TYPES = 'solute'
      DEF_TRACER_MRAT = '16.04'
      DEF_TRACER_REF_RATIO = '1.0'
      DEF_TRACER_INIT_DELTA = '0.0'
      DEF_TRACER_REACTIVE_DECAY_RATE = '0.0'
    end select
  end subroutine configure_descriptors
end program lifecycle_driver
""",
            encoding="utf-8",
        )
        ch4_parameter = tmp / "ch4_parameter.nml"
        ch4_parameter.write_text(
            """
&nl_colm_tracer_parameter
  DEF_TRACER%unit_kind = 'species_owned'
  DEF_TRACER%mol_weight = 16.04
  DEF_TRACER%ref_ratio = 1.0
  DEF_TRACER%init_delta = 0.0
  DEF_TRACER%reactive_decay_rate = 0.0
/
""",
            encoding="utf-8",
        )

        executable = tmp / "lifecycle_driver"
        command = [
            compiler,
            "-cpp",
            "-ffree-line-length-0",
            "-I",
            str(tmp),
            "-J",
            str(tmp),
            str(tmp / "precision.f90"),
            str(tmp / "namelist.f90"),
            str(tmp / "support.f90"),
            str(DEFS),
            str(LIFECYCLE),
            str(tmp / "driver.f90"),
            "-o",
            str(executable),
        ]
        compiled = subprocess.run(
            command,
            cwd=tmp,
            capture_output=True,
            text=True,
            timeout=SUBPROCESS_TIMEOUT,
        )
        assert compiled.returncode == 0, compiled.stdout + compiled.stderr
        yield executable, ch4_parameter, tmp


def run_driver(lifecycle_driver, scenario):
    executable, ch4_parameter, tmp = lifecycle_driver
    return subprocess.run(
        [str(executable), scenario, str(ch4_parameter)],
        cwd=tmp,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )


def test_route_activation_uses_registered_hooks_not_descriptor_family() -> None:
    source = LIFECYCLE.read_text(encoding="utf-8")
    body = source.split("logical FUNCTION tracer_lifecycle_route_has_active", 1)[1]
    body = body.split("END FUNCTION tracer_lifecycle_route_has_active", 1)[0]

    for hook in ("route_forcing_put", "route_diag_accumulate", "route_calc"):
        assert f"associated(lifecycle(i)%{hook})" in body
    assert "FAMILY_PARTICLE" not in body


def test_provider_aliases_resolve_once_to_descriptor_indices(lifecycle_driver):
    result = run_driver(lifecycle_driver, "happy")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "SIZE=3" in result.stdout
    assert "INDEX=2 3" in result.stdout
    assert "HAS=F T T" in result.stdout
    assert "ch4_provider=methane" in result.stdout.lower()
    assert "sed_provider=sediment" in result.stdout.lower()


def test_provider_registration_finalizes_ch4_and_sediment_descriptors(
    lifecycle_driver,
):
    result = run_driver(lifecycle_driver, "happy")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "FAMILY=solute gas particle" in result.stdout
    assert "SOLUTE_DESC=2 1 0" in result.stdout
    assert "CH4_DESC=4 2 2" in result.stdout
    assert "SED_DESC=3 2 0" in result.stdout


def test_startup_summary_reports_resolved_descriptor_provider_and_hook_groups(
    lifecycle_driver,
):
    result = run_driver(lifecycle_driver, "happy")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "TRACER lifecycle registry: count=3" in result.stdout
    assert (
        "[1] name=CL family=solute owner=generic_water reaction=none "
        "provider=- hooks=-"
    ) in result.stdout
    assert (
        "[2] name=METHANE family=gas owner=provider reaction=provider "
        "provider=methane hooks=land+reaction"
    ) in result.stdout
    assert (
        "[3] name=SED family=particle owner=provider reaction=none "
        "provider=sediment hooks=route"
    ) in result.stdout


def test_host_phase_dispatch_calls_only_the_matching_provider_hook(lifecycle_driver):
    result = run_driver(lifecycle_driver, "happy")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "AFTER_LAND=1 0" in result.stdout
    assert "AFTER_ROUTE=1 1" in result.stdout


def test_compiled_provider_absent_from_configuration_stays_inactive(lifecycle_driver):
    result = run_driver(lifecycle_driver, "absent")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "SIZE=1" in result.stdout
    assert "INDEX=0 0" in result.stdout
    assert "HAS=F" in result.stdout


def test_zero_tracers_allocates_a_safe_noop_lifecycle(lifecycle_driver):
    result = run_driver(lifecycle_driver, "zero")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "SIZE=0" in result.stdout
    assert "CALLS=0 0" in result.stdout


@pytest.mark.parametrize(
    ("scenario", "diagnostic"),
    [
        ("duplicate", "duplicate"),
        ("mismatch", "family"),
        ("unknown", "provider"),
        ("gas_incomplete", "required land"),
        ("particle_incomplete", "required route"),
        ("reaction_incomplete", "reaction-step"),
        ("ch4_without_provider", "compiled gas provider"),
    ],
)
def test_invalid_provider_contracts_fail_fast(
    lifecycle_driver, scenario, diagnostic
):
    result = run_driver(lifecycle_driver, scenario)
    assert result.returncode != 0
    assert diagnostic in (result.stdout + result.stderr).lower()


def test_lifecycle_registry_has_no_capacity_growth_or_refresh_registry():
    source = LIFECYCLE.read_text(encoding="utf-8").lower()
    assert "move_alloc" not in source
    assert "callback_capacity" not in source
    assert "refresh_dirty" not in source
    assert "class(*)" not in source
    assert "select type" not in source
    assert "tracer_reactive_" not in source
    assert "tracer_particle_" not in source
    assert "tracer_lifecycle_land_final" in source
    assert "tracer_lifecycle_route_final" in source
    assert "call tracer_lifecycle_report_registry" in source
    assert "tracer lifecycle registry: count=" in source
    for field in ("name=", "family=", "owner=", "reaction=", "provider=", "hooks="):
        assert field in source
