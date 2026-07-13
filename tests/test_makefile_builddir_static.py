from pathlib import Path
import re
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]


def make_variable(name: str) -> list[str]:
    probe = (
        "include Makefile\n"
        "print-variable:\n"
        f"\t@printf '%s\\n' '$({name})'\n"
    )
    result = subprocess.run(
        ["make", "-s", "-f", "-", "print-variable"],
        cwd=ROOT,
        input=probe,
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout.split()


class MakefileBuildDirectoryTests(unittest.TestCase):
    def test_every_object_recipe_waits_for_build_directory(self):
        makefile = (ROOT / "Makefile").read_text()
        rules = re.findall(
            r"(?m)^(?P<header>[^\n]*:[^\n]*)\n"
            r"\t(?P<recipe>[^\n]*-o \.bld/\$@[^\n]*)$",
            makefile,
        )
        self.assertGreater(len(rules), 0)
        for header, recipe in rules:
            with self.subTest(recipe=recipe):
                self.assertIn("| mkdir_build", header)

    def test_fortran_object_stages_have_predecessor_chains(self):
        makefile = (ROOT / "Makefile").read_text()
        self.assertIn("define chain_objects", makefile)
        for stage in (
            "OBJS_SHARED",
            "OBJS_MKSRFDATA",
            "OBJS_BASIC",
            "OBJS_MKINIDATA",
            "OBJECTS_CAMA",
            "OBJS_MAIN",
            "OBJS_POST1",
            "OBJS_POST2",
            "OBJS_POST3",
        ):
            with self.subTest(stage=stage):
                self.assertIn(f"$(call chain_objects,$({stage}))", makefile)

        # Do not trade the module race for a global parallel-build shutdown.
        self.assertNotRegex(makefile, r"(?m)^\.NOTPARALLEL\s*:")

    def test_cama_and_postprocess_wait_for_shared_modules(self):
        makefile = (ROOT / "Makefile").read_text()
        for target in (
            "$(OBJECTS_CAMA)",
            "$(OBJS_POST1)",
            "$(OBJS_POST2)",
            "$(OBJS_POST3)",
        ):
            match = re.search(
                rf"(?m)^{re.escape(target)}[^\n]*:[^\n]*$",
                makefile,
            )
            self.assertIsNotNone(match, target)
            with self.subTest(target=target):
                self.assertIn("${OBJS_SHARED}", match.group(0))

    def test_stage_order_respects_known_forward_module_dependencies(self):
        basic = make_variable("OBJS_BASIC")
        main = make_variable("OBJS_MAIN")

        self.assertLess(
            basic.index("MOD_Grid_RiverLakeTimeVars.o"),
            basic.index("MOD_Tracer_RiverLake.o"),
        )
        self.assertLess(
            basic.index("MOD_UserSpecifiedForcing.o"),
            basic.index("MOD_Tracer_Forcing.o"),
        )

        tracer_hist = main.index("MOD_Tracer_Hist.o")
        for provider in (
            "MOD_Vars_1DAccFluxes.o",
            "MOD_HistWriteBack.o",
            "MOD_HistGridded.o",
            "MOD_HistVector.o",
            "MOD_HistSingle.o",
        ):
            with self.subTest(provider=provider):
                self.assertLess(main.index(provider), tracer_hist)

    def test_parallel_dry_run_has_no_dependency_cycle(self):
        result = subprocess.run(
            ["make", "-n", "-j8", "colm.x"],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("Circular", result.stderr)


if __name__ == "__main__":
    unittest.main()
