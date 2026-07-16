from pathlib import Path
import subprocess
import tempfile

import pytest

from fortran_test_support import require_runnable_fortran_compiler


ROOT = Path(__file__).resolve().parents[1]
DEFS = ROOT / "main" / "TRACER" / "MOD_Tracer_Defs.F90"
METHANE_REGISTRY = (
    ROOT / "main" / "TRACER" / "MOD_Tracer_Reactive_Methane_Registry.F90"
)
STANDARD_CHLORIDE = ROOT / "run" / "standard_chloride_parameter.nml"
SUBPROCESS_TIMEOUT = 60


def test_species_owned_setter_keeps_descriptor_units_coherent():
    source = DEFS.read_text(encoding="utf-8")
    setter = source.split("SUBROUTINE tracer_set_land_water_transport", 1)[1].split(
        "END SUBROUTINE tracer_set_land_water_transport", 1
    )[0]
    assert "tracers(itrc)%unit_kind = 'species_owned'" in setter


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
  character(len=256) :: DEF_TRACER_TYPES = 'conservative'
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
#ifdef BGC
  use MOD_Tracer_Reactive_Methane_Registry, only: methane_registry_init
#endif
  implicit none
  character(len=512) :: parameter_file
  character(len=256) :: tracer_name, tracer_category
  character(len=32) :: action, tracer_num_text
  integer :: i

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
#ifdef BGC
  if (trim(action) == 'methane') call methane_registry_init()
#endif
  do i = 1, ntracers
    write(*,'(A)') 'NAME=' // trim(tracers(i)%name)
  enddo
  write(*,'(3(ES24.16,1X))') tracers(1)%init_conc, &
    tracers(1)%precip_default_conc, tracers(1)%vapor_default_conc
  write(*,'(A)') trim(tracers(1)%unit_kind)
  write(*,'(A)') trim(tracer_concentration_units(1))
  write(*,'(I0)') tracers(1)%charge
  write(*,'(3(L1,1X))') tracers(1)%uses_land_water_transport, &
    tracers(1)%is_nonvolatile, tracers(1)%uses_reaction
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
            str(METHANE_REGISTRY),
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
        command_no_bgc = [item for item in command if item != str(METHANE_REGISTRY)]
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
    tracer_category="conservative",
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


def test_legacy_conservative_defaults_are_preserved(descriptor_driver):
    result = run_driver(descriptor_driver)
    assert result.returncode == 0, result.stdout + result.stderr
    lines = result.stdout.splitlines()
    assert lines[-4:] == [
        "tracer_per_water",
        "tracer/water",
        "0",
        "T T F",
    ]


def test_explicit_chloride_descriptor_is_loaded(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        STANDARD_CHLORIDE.read_text(encoding="utf-8"),
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines()[-4:] == [
        "mass_fraction",
        "kg/kg water",
        "-1",
        "T T F",
    ]


@pytest.mark.parametrize(
    ("filename", "name", "category", "expected"),
    [
        ("standard_O18_parameter.nml", "H2_18O", "isotope", "ratio\nR\n0\nT F F"),
        ("standard_HDO_parameter.nml", "HDO", "isotope", "ratio\nR\n0\nT F F"),
        (
            "standard_ch4_parameter.nml",
            "CH4",
            "reactive",
            "species_owned\nspecies-owned\n0\nF F T",
        ),
        (
            "standard_sediment_parameter.nml",
            "SEDIMENT",
            "particle",
            "volume_fraction\nm3/m3 water\n0\nF F F",
        ),
    ],
)
def test_standard_species_descriptors_are_valid(
    descriptor_driver, filename, name, category, expected
):
    result = run_driver(
        descriptor_driver,
        (ROOT / "run" / filename).read_text(encoding="utf-8"),
        tracer_name=name,
        tracer_category=category,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert expected in result.stdout


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
    ("field", "value"),
    [
        ("uses_generic_land_water_transport", ".false."),
        ("is_nonvolatile", ".false."),
        ("uses_reaction", ".true."),
        ("unit_kind", "'mol_m3_water'"),
    ],
)
def test_incoherent_descriptor_is_rejected(
    descriptor_driver, field, value
):
    result = run_driver(
        descriptor_driver,
        f"""
&nl_colm_tracer_parameter
  DEF_TRACER%{field} = {value}
/
""",
    )
    assert result.returncode != 0
    assert field in result.stdout


def test_species_only_legacy_parameter_file_keeps_defaults(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_methane_parameter
/
""",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines()[-4:] == [
        "tracer_per_water",
        "tracer/water",
        "0",
        "T T F",
    ]


def test_reactive_category_cannot_disable_reaction_capability(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%unit_kind = 'species_owned'
  DEF_TRACER%uses_generic_land_water_transport = .false.
  DEF_TRACER%uses_reaction = .false.
/
""",
        tracer_name="CH4",
        tracer_category="reactive",
    )
    assert result.returncode != 0
    assert "uses_reaction" in result.stdout


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
    values = [float(value) for value in result.stdout.splitlines()[-5].split()]
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
    values = [float(value) for value in result.stdout.splitlines()[-5].split()]
    assert values == pytest.approx([1.0, 2.0, 3.0])


def test_generic_reactive_nonvolatile_solute_is_valid(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%unit_kind = 'mass_fraction'
  DEF_TRACER%charge = -1
  DEF_TRACER%is_nonvolatile = .true.
  DEF_TRACER%uses_reaction = .true.
  DEF_TRACER%reactive_decay_rate = 1.0e-6
/
""",
        tracer_name="NO3",
        tracer_category="reactive",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.splitlines()[-4:] == [
        "mass_fraction",
        "kg/kg water",
        "-1",
        "T T T",
    ]


def test_bgc_off_rejects_species_owned_ch4_but_keeps_generic_reactive_fallback(
    descriptor_driver,
):
    methane = run_driver(
        descriptor_driver,
        (ROOT / "run" / "standard_ch4_parameter.nml").read_text(encoding="utf-8"),
        tracer_name="CH4",
        tracer_category="reactive",
        bgc=False,
    )
    assert methane.returncode != 0
    assert "species-owned CH4 requires compiling with BGC" in methane.stdout

    generic = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%uses_reaction = .true.
  DEF_TRACER%reactive_decay_rate = 1.0e-6
/
""",
        tracer_name="NO3",
        tracer_category="reactive",
        bgc=False,
    )
    assert generic.returncode == 0, generic.stdout + generic.stderr
    assert generic.stdout.splitlines()[-1] == "T F T"


def test_species_owned_reactive_cannot_claim_nonvolatile_residue(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%unit_kind = 'species_owned'
  DEF_TRACER%uses_generic_land_water_transport = .false.
  DEF_TRACER%is_nonvolatile = .true.
  DEF_TRACER%uses_reaction = .true.
/
""",
        tracer_name="CH4",
        tracer_category="reactive",
    )
    assert result.returncode != 0
    assert "is_nonvolatile" in result.stdout


def test_isotope_cannot_claim_nonvolatile_solute_transport(descriptor_driver):
    result = run_driver(
        descriptor_driver,
        """
&nl_colm_tracer_parameter
  DEF_TRACER%is_nonvolatile = .true.
/
""",
        tracer_name="H2_18O",
        tracer_category="isotope",
    )
    assert result.returncode != 0
    assert "is_nonvolatile" in result.stdout


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
        tracer_category="reactive,reactive",
        tracer_num=2,
        action="methane",
    )
    assert result.returncode != 0
    assert "CH4 and METHANE aliases" in result.stdout


@pytest.mark.parametrize(
    ("field", "value", "name", "category"),
    [
        ("ref_ratio", "0.0", "H2_18O", "isotope"),
        ("init_delta", "-1000.1", "H2_18O", "isotope"),
        ("init_conc", "-1.0", "CL", "conservative"),
        ("precip_default_conc", "-1.0", "CL", "conservative"),
        ("vapor_default_conc", "-1.0", "CL", "conservative"),
        ("mol_weight", "NaN", "CL", "conservative"),
        ("mol_weight", "Inf", "CL", "conservative"),
        ("reactive_decay_rate", "NaN", "NO3", "reactive"),
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
