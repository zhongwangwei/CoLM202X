from pathlib import Path
import subprocess
import tempfile

import pytest

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
DEFS = ROOT / "main" / "TRACER" / "MOD_Tracer_Defs.F90"
EVAP_LIMIT = ROOT / "main" / "TRACER" / "MOD_Tracer_EvapLimit.F90"
STANDARD_CHLORIDE = ROOT / "run" / "standard_chloride_parameter.nml"
SUBPROCESS_TIMEOUT = 60


def test_descriptor_taxonomy_has_no_mutable_capability_booleans_or_setter():
    source = DEFS.read_text(encoding="utf-8")
    type_block = source.split("type :: tracer_info_type", 1)[1].split(
        "end type tracer_info_type", 1
    )[0]
    parameter_block = source.split("type :: tracer_parameter_type", 1)[1].split(
        "end type tracer_parameter_type", 1
    )[0]
    for field in (
        "uses_generic_land_water_transport",
        "is_nonvolatile",
        "uses_reaction",
        "has_fractionation",
    ):
        assert field not in type_block
        assert field not in parameter_block
    assert "tracer_set_land_water_transport" not in source


@pytest.fixture(scope="module")
def descriptor_driver():
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
  integer :: DEF_TRACER_NUM = 1
  character(len=256) :: DEF_TRACER_NAMES = 'CL'
  character(len=256) :: DEF_TRACER_TYPES = 'solute'
  character(len=256) :: DEF_TRACER_MRAT = '35.453'
  character(len=256) :: DEF_TRACER_REF_RATIO = '1.0'
  character(len=256) :: DEF_TRACER_INIT_DELTA = '0.0'
  character(len=256) :: DEF_TRACER_REACTIVE_DECAY_RATE = '0.0'
  character(len=512) :: DEF_TRACER_PARAM_FILES = 'null'
  logical :: DEF_TRACER_USE_FRACTIONATION = .false.
end module MOD_Namelist
""",
            encoding="utf-8",
        )
        (tmp / "spmd.f90").write_text(
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
""",
            encoding="utf-8",
        )
        (tmp / "vars_global.f90").write_text(
            """
module MOD_Vars_Global
  use MOD_Precision, only: r8
  implicit none
  real(r8), parameter :: spval = -1.0e36_r8
end module MOD_Vars_Global
""",
            encoding="utf-8",
        )
        (tmp / "driver.f90").write_text(
            """
#include <define.h>
program descriptor_driver
  use MOD_Namelist
  use MOD_Tracer_Defs
  use MOD_Tracer_EvapLimit, only: tracer_atmospheric_tracer_loss
  implicit none
  character(len=512) :: parameter_file
  character(len=256) :: tracer_name, tracer_category
  character(len=TRACER_DESCRIPTOR_IDENTITY_WIDTH) :: descriptor_text
  character(len=32) :: action, tracer_num_text
  integer :: i, k
  integer, allocatable :: descriptor_identity(:,:)
  real(r8) :: dissolved, solid, evap_loss

  call get_command_argument(1, parameter_file)
  call get_command_argument(2, tracer_name)
  call get_command_argument(3, tracer_category)
  call get_command_argument(4, tracer_num_text)
  call get_command_argument(5, action)
  if (len_trim(tracer_num_text) > 0) read(tracer_num_text, *) DEF_TRACER_NUM
  if (len_trim(tracer_name) > 0) DEF_TRACER_NAMES = trim(tracer_name)
  if (len_trim(tracer_category) > 0) DEF_TRACER_TYPES = trim(tracer_category)
  if (len_trim(parameter_file) > 0) then
    DEF_TRACER_PARAM_FILES = trim(DEF_TRACER_NAMES) // ':' // trim(parameter_file)
  endif
  call tracer_defs_init()
  if (ntracers == 0) then
    write(*,'(A,I0)') 'NTRACERS=', ntracers
    write(*,'(A,L1)') 'ALLOCATED=', allocated(tracers)
    write(*,'(A,10(L1,1X))') 'PRED=', tracer_is_isotope(1), tracer_is_solute(1), &
      tracer_is_gas(1), tracer_is_particle(1), tracer_is_conservative(1), &
      tracer_is_reactive(1), tracer_is_nonvolatile_solute(1), &
      tracer_uses_land_water_transport(1), tracer_uses_delta_diagnostics(1), &
      tracer_can_use_fixed_signature(1)
    write(*,'(A,3(I0,1X))') 'TAX=', tracer_get_family(1), &
      tracer_get_state_owner(1), tracer_get_reaction_mode(1)
    stop
  endif
  if (trim(action) == 'methane') i = tracer_index_for_name('CH4', 'METHANE')
  if (trim(action) == 'max') then
    write(*,'(A,ES24.16)') 'MAXC=', tracers(1)%max_dissolved_conc
    stop
  endif
  if (trim(action) == 'equilibrate') then
    dissolved = 2.0e-2_r8
    solid = 0.0_r8
    call tracer_equilibrate_dissolved(1, 1.0_r8, dissolved, solid)
    write(*,'(A,2(ES24.16,1X))') 'DRY=', dissolved, solid
    call tracer_equilibrate_dissolved(1, 2.0_r8, dissolved, solid)
    write(*,'(A,2(ES24.16,1X))') 'REWET=', dissolved, solid
    stop
  endif
  if (trim(action) == 'evaporation') then
    evap_loss = tracer_atmospheric_tracer_loss(0.5_r8, 10.0_r8, 2.0_r8, &
      293.15_r8, .false., unchanged_ratio, trc_tiny, -1.0_r8, &
      tracer_is_nonvolatile_solute(1))
    write(*,'(A,ES24.16)') 'EVAP_LOSS=', evap_loss
    stop
  endif
  if (trim(action) == 'identity' .or. trim(action) == 'transport_identity') then
    if (trim(action) == 'transport_identity') then
      call tracer_build_descriptor_identity(descriptor_identity, transport_only=.true.)
    else
      call tracer_build_descriptor_identity(descriptor_identity)
    endif
    descriptor_text = ''
    write(*,'(A,I0)') 'IDENTITY_COUNT=', size(descriptor_identity, 2)
    if (size(descriptor_identity, 2) > 0) then
      do k = 1, TRACER_DESCRIPTOR_IDENTITY_WIDTH
        if (descriptor_identity(k, 1) > 0) &
          descriptor_text(k:k) = achar(descriptor_identity(k, 1))
      enddo
      write(*,'(A)') 'IDENTITY=' // trim(descriptor_text)
    endif
    stop
  endif
  do i = 1, ntracers
    write(*,'(A)') 'NAME=' // trim(tracers(i)%name)
  enddo
  write(*,'(A,10(L1,1X))') 'PRED=', tracer_is_isotope(1), tracer_is_solute(1), &
    tracer_is_gas(1), tracer_is_particle(1), tracer_is_conservative(1), &
    tracer_is_reactive(1), tracer_is_nonvolatile_solute(1), &
      tracer_uses_land_water_transport(1), tracer_uses_delta_diagnostics(1), &
      tracer_can_use_fixed_signature(1)
  write(*,'(A)') 'CATEGORY=' // trim(tracers(1)%category)
  write(*,'(A,3(I0,1X))') 'TAX=', tracer_get_family(1), &
    tracer_get_state_owner(1), tracer_get_reaction_mode(1)
  write(*,'(A,3(ES24.16,1X))') 'CONC=', tracers(1)%init_conc, &
    tracers(1)%precip_default_conc, tracers(1)%vapor_default_conc
  write(*,'(A)') 'UNIT_KIND=' // trim(tracers(1)%unit_kind)
  write(*,'(A)') 'UNITS=' // trim(tracer_concentration_units(1))
  write(*,'(A,I0)') 'CHARGE=', tracers(1)%charge

contains
  real(r8) function unchanged_ratio(source_ratio, temp_k, from_ice)
    real(r8), intent(in) :: source_ratio, temp_k
    logical, intent(in) :: from_ice
    unchanged_ratio = source_ratio
  end function unchanged_ratio
end program descriptor_driver
""",
            encoding="utf-8",
        )

        executable = tmp / "descriptor_driver"
        command = [
            compiler,
            "-cpp",
            "-ffree-line-length-0",
            "-I",
            str(tmp),
            str(tmp / "precision.f90"),
            str(tmp / "namelist.f90"),
            str(tmp / "spmd.f90"),
            str(tmp / "vars_global.f90"),
            str(DEFS),
            str(EVAP_LIMIT),
            str(tmp / "driver.f90"),
            "-o",
            str(executable),
        ]
        subprocess.run(
            command,
            cwd=tmp,
            check=True,
            capture_output=True,
            text=True,
            timeout=SUBPROCESS_TIMEOUT,
        )

        (tmp / "define.h").write_text("#define TRACER\n", encoding="utf-8")
        executable_no_bgc = tmp / "descriptor_driver_no_bgc"
        command_no_bgc = list(command)
        command_no_bgc[command_no_bgc.index(str(executable))] = str(executable_no_bgc)
        subprocess.run(
            command_no_bgc,
            cwd=tmp,
            check=True,
            capture_output=True,
            text=True,
            timeout=SUBPROCESS_TIMEOUT,
        )
        yield executable, executable_no_bgc, tmp


def run_driver(
    descriptor_driver,
    parameter_text=None,
    tracer_name="CL",
    tracer_category="solute",
    tracer_num=1,
    action="none",
    bgc=True,
):
    executable, executable_no_bgc, tmp = descriptor_driver
    if not bgc:
        executable = executable_no_bgc
    parameter_arg = ""
    if parameter_text is not None:
        parameter_file = tmp / "parameter.nml"
        parameter_file.write_text(parameter_text, encoding="utf-8")
        parameter_arg = str(parameter_file)
    command = [
        str(executable),
        parameter_arg,
        tracer_name,
        tracer_category,
        str(tracer_num),
        action,
    ]
    return subprocess.run(
        command,
        cwd=tmp,
        capture_output=True,
        text=True,
        timeout=SUBPROCESS_TIMEOUT,
    )


def test_canonical_solute_defaults_are_derived(descriptor_driver):
    result = run_driver(descriptor_driver)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "PRED=F T F F T F T T F F" in result.stdout
    assert "CATEGORY=solute" in result.stdout
    assert "TAX=2 1 0" in result.stdout
    assert "UNIT_KIND=tracer_per_water" in result.stdout
    assert "UNITS=tracer/water" in result.stdout
    assert "CHARGE=0" in result.stdout


def test_legacy_conservative_alias_maps_to_solute(descriptor_driver):
    result = run_driver(descriptor_driver, tracer_category="conservative")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "PRED=F T F F T F T T F F" in result.stdout
    assert "CATEGORY=solute" in result.stdout
    assert "TAX=2 1 0" in result.stdout


def test_sanitized_tracer_names_are_unique_case_insensitively(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        tracer_name="CL,cl",
        tracer_category="solute,solute",
        tracer_num=2,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "NAME=CL\nNAME=cl_2\n" in result.stdout


def test_explicit_chloride_descriptor_is_loaded(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        STANDARD_CHLORIDE.read_text(encoding="utf-8"),
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "UNIT_KIND=mass_fraction" in result.stdout
    assert "UNITS=kg/kg water" in result.stdout
    assert "CHARGE=-1" in result.stdout
    assert "PRED=F T F F T F T T F F" in result.stdout
    assert "TAX=2 1 0" in result.stdout


def test_exact_restart_identity_changes_when_physical_descriptor_changes(
    descriptor_driver,
):
    baseline = run_driver(descriptor_driver, action="identity")
    chloride = run_driver(
        descriptor_driver,
        STANDARD_CHLORIDE.read_text(encoding="utf-8"),
        action="identity",
    )
    assert baseline.returncode == 0, baseline.stdout + baseline.stderr
    assert chloride.returncode == 0, chloride.stdout + chloride.stderr
    baseline_identity = next(
        line for line in baseline.stdout.splitlines() if line.startswith("IDENTITY=")
    )
    chloride_identity = next(
        line for line in chloride.stdout.splitlines() if line.startswith("IDENTITY=")
    )
    assert "IDENTITY_COUNT=1" in baseline.stdout
    assert "CL|solute|tracer_per_water|2|1|0|0|" in baseline_identity
    assert "CL|solute|mass_fraction|2|1|0|-1|" in chloride_identity
    assert chloride_identity != baseline_identity


def test_transport_identity_is_empty_for_provider_owned_particle(
    descriptor_driver,
):
    result = run_driver(
        descriptor_driver,
        tracer_name="SEDIMENT",
        tracer_category="particle",
        action="transport_identity",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines() == ["IDENTITY_COUNT=0"]


def test_zero_tracers_is_a_safe_noop(descriptor_driver):
    result = run_driver(descriptor_driver, tracer_num=0)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines() == [
        "NTRACERS=0",
        "ALLOCATED=F",
        "PRED=F F F F F F F F F F",
        "TAX=0 0 0",
    ]


def test_chloride_restores_the_historical_dissolved_concentration_limit(
    descriptor_driver,
):
    result = run_driver(
        descriptor_driver,
        STANDARD_CHLORIDE.read_text(encoding="utf-8"),
        action="max",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    line = next(line for line in result.stdout.splitlines() if line.startswith("MAXC="))
    assert float(line.removeprefix("MAXC=")) == pytest.approx(1.0e-2)


def test_dissolved_limit_precipitates_and_redissolves_without_losing_mass(
    descriptor_driver,
):
    result = run_driver(
        descriptor_driver,
        STANDARD_CHLORIDE.read_text(encoding="utf-8"),
        action="equilibrate",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    lines = result.stdout.splitlines()
    dry = [
        float(value)
        for value in next(line for line in lines if line.startswith("DRY="))
        .removeprefix("DRY=")
        .split()
    ]
    rewet = [
        float(value)
        for value in next(line for line in lines if line.startswith("REWET="))
        .removeprefix("REWET=")
        .split()
    ]
    assert dry == pytest.approx([1.0e-2, 1.0e-2])
    assert rewet == pytest.approx([2.0e-2, 0.0])
    assert sum(dry) == pytest.approx(2.0e-2)
    assert sum(rewet) == pytest.approx(2.0e-2)


@pytest.mark.parametrize(
    (
        "filename",
        "name",
        "category",
        "expected_category",
        "predicates",
        "taxonomy",
        "unit_kind",
        "units",
    ),
    [
        (
            "standard_O18_parameter.nml",
            "H2_18O",
            "isotope",
            "isotope",
            "PRED=T F F F F F F T T T",
            "TAX=1 1 0",
            "ratio",
            "R",
        ),
        (
            "standard_HDO_parameter.nml",
            "HDO",
            "isotope",
            "isotope",
            "PRED=T F F F F F F T T T",
            "TAX=1 1 0",
            "ratio",
            "R",
        ),
        (
            "standard_ch4_parameter.nml",
            "CH4",
            "gas",
            "gas",
            "PRED=F F T F T F F F F F",
            "TAX=4 2 0",
            "species_owned",
            "species-owned",
        ),
        (
            "standard_sediment_parameter.nml",
            "SEDIMENT",
            "particle",
            "particle",
            "PRED=F F F T F F F F F F",
            "TAX=3 2 0",
            "volume_fraction",
            "m3/m3 water",
        ),
    ],
)
def test_standard_species_descriptors_are_valid(
    descriptor_driver,
    filename,
    name,
    category,
    expected_category,
    predicates,
    taxonomy,
    unit_kind,
    units,
):
    result = run_driver(
        descriptor_driver,
        (ROOT / "run" / filename).read_text(encoding="utf-8"),
        tracer_name=name,
        tracer_category=category,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert f"CATEGORY={expected_category}" in result.stdout
    assert predicates in result.stdout
    assert taxonomy in result.stdout
    assert f"UNIT_KIND={unit_kind}" in result.stdout
    assert f"UNITS={units}" in result.stdout
    assert "CHARGE=0" in result.stdout


def test_malformed_parameter_file_fails_fast(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%unknown_field = 1
/
""",
    )
    assert result.returncode != 0
    assert "invalid &nl_colm_tracer_parameter" in result.stdout


@pytest.mark.parametrize(
    "field",
    [
        "uses_generic_land_water_transport",
        "is_nonvolatile",
        "uses_reaction",
        "has_fractionation",
    ],
)
def test_removed_free_boolean_fields_fail_fast(descriptor_driver, field):
    result = run_driver(
        descriptor_driver,
        f"""
&nl_colm_tracer_parameter
  DEF_TRACER%{field} = .true.
/
""",
    )
    assert result.returncode != 0
    assert field in result.stdout


def test_unsupported_unit_kind_is_rejected(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%unit_kind = 'mol_m3_water'
/
""",
    )
    assert result.returncode != 0
    assert "unit_kind" in result.stdout


def test_species_only_legacy_parameter_file_keeps_defaults(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_methane_parameter
/
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "CATEGORY=solute" in result.stdout
    assert "TAX=2 1 0" in result.stdout
    assert "UNIT_KIND=tracer_per_water" in result.stdout


def test_legacy_init_delta_remains_all_concentration_fallbacks(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%init_delta = 0.25
/
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    line = next(line for line in result.stdout.splitlines() if line.startswith("CONC="))
    values = [float(value) for value in line.removeprefix("CONC=").split()]
    assert values == pytest.approx([0.25, 0.25, 0.25])


def test_explicit_solute_concentrations_are_independent(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%init_delta = 0.25
  DEF_TRACER%init_conc = 1.0
  DEF_TRACER%precip_default_conc = 2.0
  DEF_TRACER%vapor_default_conc = 3.0
/
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    line = next(line for line in result.stdout.splitlines() if line.startswith("CONC="))
    values = [float(value) for value in line.removeprefix("CONC=").split()]
    assert values == pytest.approx([1.0, 2.0, 3.0])


def test_solute_positive_decay_derives_first_order_reaction(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%unit_kind = 'mass_fraction'
  DEF_TRACER%charge = -1
  DEF_TRACER%reactive_decay_rate = 1.0e-6
/
""",
        tracer_name="NO3",
        tracer_category="solute",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "CATEGORY=solute" in result.stdout
    assert "PRED=F T F F F T T T F F" in result.stdout
    assert "TAX=2 1 1" in result.stdout
    assert "UNIT_KIND=mass_fraction" in result.stdout
    assert "CHARGE=-1" in result.stdout


def test_generic_solute_is_retained_on_evaporation(descriptor_driver):
    result = run_driver(descriptor_driver, action="evaporation")
    assert result.returncode == 0, result.stdout + result.stderr
    line = next(
        line for line in result.stdout.splitlines() if line.startswith("EVAP_LOSS=")
    )
    assert float(line.removeprefix("EVAP_LOSS=")) == pytest.approx(0.0)


def test_bgc_off_rejects_species_owned_ch4_but_keeps_first_order_solute(
    descriptor_driver,
):
    methane = run_driver(
        descriptor_driver,
        (ROOT / "run" / "standard_ch4_parameter.nml").read_text(encoding="utf-8"),
        tracer_name="CH4",
        tracer_category="gas",
        bgc=False,
    )
    assert methane.returncode != 0
    assert "species-owned CH4 requires compiling with BGC" in methane.stdout

    generic = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%reactive_decay_rate = 1.0e-6
/
""",
        tracer_name="NO3",
        tracer_category="solute",
        bgc=False,
    )
    assert generic.returncode == 0, generic.stdout + generic.stderr
    assert "CATEGORY=solute" in generic.stdout
    assert "PRED=F T F F F T T T F F" in generic.stdout
    assert "TAX=2 1 1" in generic.stdout


def test_sanitized_names_remain_unique_after_suffix_collision(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        tracer_name="A,A_3,A",
        tracer_category="isotope,isotope,isotope",
        tracer_num=3,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    names = [line.removeprefix("NAME=") for line in result.stdout.splitlines() if line.startswith("NAME=")]
    assert names == ["A", "A_3", "A_4"]
    assert len(names) == len(set(names))


def test_ch4_aliases_cannot_both_be_registered(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        tracer_name="CH4,METHANE",
        tracer_category="gas,gas",
        tracer_num=2,
        action="methane",
    )
    assert result.returncode != 0
    assert "match multiple configured tracers" in result.stdout


@pytest.mark.parametrize(
    ("field", "value", "name", "category"),
    [
        ("ref_ratio", "0.0", "H2_18O", "isotope"),
        ("init_delta", "-1000.1", "H2_18O", "isotope"),
        ("init_conc", "-1.0", "CL", "solute"),
        ("precip_default_conc", "-1.0", "CL", "solute"),
        ("vapor_default_conc", "-1.0", "CL", "solute"),
        ("max_dissolved_conc", "0.0", "CL", "solute"),
        ("max_dissolved_conc", "NaN", "CL", "solute"),
        ("mol_weight", "NaN", "CL", "solute"),
        ("mol_weight", "Inf", "CL", "solute"),
        ("reactive_decay_rate", "NaN", "NO3", "solute"),
    ],
)
def test_nonphysical_or_nonfinite_descriptor_values_fail_fast(
    descriptor_driver, field, value, name, category
):
    result = run_driver(
        descriptor_driver,
        f"""
&nl_colm_tracer_parameter
  DEF_TRACER%{field} = {value}
/
""",
        tracer_name=name,
        tracer_category=category,
    )
    assert result.returncode != 0
    assert field in result.stdout
