from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main" / "TRACER"


def source(name: str) -> str:
    return (TRACER / name).read_text()


def subroutine_body(text: str, name: str) -> str:
    match = re.search(
        rf"SUBROUTINE\s+{name}\b(?P<body>.*?)END\s+SUBROUTINE\s+{name}",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if match is None:
        raise AssertionError(f"subroutine {name} not found")
    return match.group("body")


def function_body(text: str, name: str) -> str:
    match = re.search(
        rf"(?:logical\s+)?FUNCTION\s+{name}\b(?P<body>.*?)END\s+FUNCTION\s+{name}",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if match is None:
        raise AssertionError(f"function {name} not found")
    return match.group("body")


class ReactiveArchitectureStaticTests(unittest.TestCase):
    def test_land_phase_facade_is_private_by_default(self):
        land_phase = source("MOD_Tracer_LandPhase.F90")
        declarations = land_phase[: land_phase.index("CONTAINS")]

        self.assertRegex(declarations, r"(?mi)^\s*PRIVATE\s*$")
        self.assertIn("PUBLIC :: land_tracer_init, land_tracer_final", declarations)
        self.assertIn("PUBLIC :: tracer_soil_step, tracer_report", declarations)

    def test_land_tracer_state_hides_lulcc_workspaces(self):
        variables = source("MOD_Tracer_Vars.F90")
        declarations = variables[: variables.index("CONTAINS")]

        self.assertRegex(declarations, r"(?mi)^\s*PRIVATE\s*$")
        self.assertIn("PUBLIC :: trc_wliq_soisno, trc_wice_soisno", declarations)
        self.assertIn(
            "PUBLIC :: save_land_tracer_lulcc_state, remap_land_tracer_lulcc_state",
            declarations,
        )
        self.assertNotRegex(declarations, r"(?mi)^\s*PUBLIC\s*::[^\n]*lulcc_trc_")

    def test_species_can_disable_generic_land_water_transport(self):
        defs = source("MOD_Tracer_Defs.F90")
        registry = source("MOD_Tracer_Reactive_Methane_Registry.F90")
        land_phase = source("MOD_Tracer_LandPhase.F90")

        self.assertIn("logical           :: uses_land_water_transport = .true.", defs)
        self.assertIn("PUBLIC :: tracer_set_land_water_transport", defs)
        uses_body = function_body(defs, "tracer_uses_land_water_transport")
        self.assertIn("tracers(itrc)%uses_land_water_transport", uses_body)
        self.assertIn("CALL tracer_set_land_water_transport (i, .false.)", registry)
        init_body = subroutine_body(land_phase, "tracer_init_from_arrays")
        self.assertLess(
            init_body.index("CALL tracer_reactive_init"),
            init_body.index("CALL allocate_Tracer_Vars"),
        )

    def test_particle_final_does_not_requery_activation(self):
        body = subroutine_body(source("MOD_Tracer_Particle.F90"), "tracer_particle_final")
        self.assertNotIn("particle_callback_enabled", body)
        self.assertIn("associated(particle_callbacks(i)%final)", body)

    def test_reactive_history_callbacks_receive_private_workspaces(self):
        body = subroutine_body(source("MOD_Tracer_Reactive.F90"), "tracer_reactive_history")
        self.assertIn("callback_sumarea = sumarea", body)
        self.assertIn("callback_filter = filter", body)
        self.assertIn("callback_sumarea, callback_filter", body)
        self.assertNotIn(
            "CALL reactive_callbacks(i)%history (file_hist, itime_in_file, sumarea, filter",
            body,
        )

    def test_generic_reactive_tracers_keep_generic_history(self):
        body = subroutine_body(source("MOD_Tracer_Hist.F90"), "tracer_hist_out")
        self.assertIn(
            "IF (.not. tracer_uses_land_water_transport(itrc_loc)) CYCLE",
            body,
        )
        self.assertNotIn("tracer_is_reactive(itrc_loc)", body)

    def test_methane_example_registers_only_implemented_species(self):
        example = (
            ROOT / "run" / "examples" / "Global_Grid_2x2_PFT_VG_BGC_CH4.nml"
        ).read_text()
        self.assertIn("DEF_TRACER_NUM         = 1", example)
        self.assertIn('DEF_TRACER_NAMES       = "CH4"', example)
        self.assertNotRegex(example, r"(?i)DEF_TRACER_NAMES\s*=.*\b(?:O2|CO2)\b")


if __name__ == "__main__":
    unittest.main()
