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
    def test_retired_family_registries_and_stubs_are_deleted(self):
        retired = (
            TRACER / "MOD_Tracer_Reactive.F90",
            TRACER / "MOD_Tracer_Particle.F90",
            ROOT / "include" / "tracer_reactive_species.inc",
            ROOT / "include" / "tracer_particle_species.inc",
            ROOT / "mkinidata" / "MOD_Tracer_Reactive_Registrations_Stubs.F90",
            ROOT / "mkinidata" / "MOD_Tracer_Particle_Registrations_Stubs.F90",
        )
        for path in retired:
            with self.subTest(path=path):
                self.assertFalse(path.exists())

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

    def test_state_owner_controls_generic_land_water_transport(self):
        defs = source("MOD_Tracer_Defs.F90")
        lifecycle = source("MOD_Tracer_Lifecycle.F90")
        land_phase = source("MOD_Tracer_LandPhase.F90")

        uses_body = function_body(defs, "tracer_uses_land_water_transport")
        self.assertIn(
            "tracers(itrc)%state_owner == STATE_OWNER_GENERIC_WATER", uses_body
        )
        registration = subroutine_body(lifecycle, "register_tracer_provider")
        self.assertIn("tracers(itrc)%state_owner = declared_owner", registration)
        init_body = subroutine_body(land_phase, "tracer_init_from_arrays")
        self.assertLess(
            init_body.index("CALL tracer_lifecycle_init()"),
            init_body.index("CALL allocate_Tracer_Vars"),
        )

    def test_provider_owned_rows_do_not_retain_generic_land_water_state(self):
        variables = source("MOD_Tracer_Vars.F90")
        restart = source("MOD_Tracer_Rest.F90")
        body = subroutine_body(
            variables, "zero_provider_owned_land_tracer_state"
        )

        self.assertIn("IF (tracer_uses_land_water_transport(itrc)) CYCLE", body)
        self.assertNotIn("tracer_is_particle", body)
        # Cold start, pre-read clearing, post-unpack clearing, and write all
        # enforce the same ownership boundary. Provider rows are never packed
        # into the compact generic restart transaction.
        self.assertEqual(
            restart.count("CALL zero_provider_owned_land_tracer_state()"), 4
        )
        self.assertIn("tracer_build_descriptor_identity(expected, transport_only=.true.)", restart)
        self.assertIn("IF (.not. tracer_uses_land_water_transport(itrc)) CYCLE", restart)
        self.assertNotIn("zero_particle_land_tracer_state", variables + restart)

    def test_route_final_uses_registered_rows_without_requerying_names(self):
        body = subroutine_body(
            source("MOD_Tracer_Lifecycle.F90"), "tracer_lifecycle_route_final"
        )
        self.assertIn("lifecycle_row_registered(i)", body)
        self.assertIn("associated(lifecycle(i)%route_final)", body)
        self.assertIn("CALL lifecycle(i)%route_final()", body)
        self.assertNotIn("tracer_index_for_name", body)

    def test_multiple_lifecycle_history_providers_receive_private_workspaces(self):
        body = subroutine_body(
            source("MOD_Tracer_Lifecycle.F90"), "tracer_lifecycle_land_history"
        )
        self.assertIn("callback_sumarea = sumarea", body)
        self.assertIn("callback_filter = filter", body)
        self.assertIn(
            "%land_history(file_hist, itime_in_file, callback_sumarea, callback_filter",
            body,
        )
        self.assertNotIn(
            "CALL lifecycle(i)%land_history(file_hist, itime_in_file, sumarea, filter",
            body,
        )

    def test_compiled_providers_share_one_manifest_and_registrar(self):
        manifest = (ROOT / "include" / "tracer_lifecycle_providers.inc").read_text()
        registrar = source("MOD_Tracer_Lifecycle_Registrations.F90")

        self.assertIn("ch4_register_tracer_provider", manifest)
        self.assertIn("register_sediment_tracer_provider", manifest)
        self.assertEqual(registrar.count('#include "tracer_lifecycle_providers.inc"'), 2)
        self.assertIn("CALL register_fn ()", registrar)

    def test_methane_parameter_lookup_accepts_both_name_aliases(self):
        body = subroutine_body(
            source("MOD_Tracer_Reactive_Methane.F90"),
            "ch4_reactive_resolve_param_file",
        )
        self.assertIn(
            "tracer_param_file_for_index (igas_ch4, 'CH4,METHANE'", body
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
