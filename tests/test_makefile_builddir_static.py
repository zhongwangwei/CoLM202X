from pathlib import Path
import re
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRS = (
    "include",
    "share",
    "mksrfdata",
    "mkinidata",
    "main",
    "main/TRACER",
    "main/HYDRO",
    "main/BGC",
    "main/URBAN",
    "main/LULCC",
    "main/DA",
    "main/ParaOpt",
    "extends/CaMa/src",
    "postprocess",
)


def make_variable(name: str, overrides: tuple[str, ...] = ()) -> list[str]:
    probe = (
        "include Makefile\n"
        "print-variable:\n"
        f"\t@printf '%s\\n' '$({name})'\n"
    )
    result = subprocess.run(
        ["make", "-s", "-f", "-", *overrides, "print-variable"],
        cwd=ROOT,
        input=probe,
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout.split()


def source_module_graph() -> tuple[dict[str, str], dict[str, set[str]]]:
    providers: dict[str, str] = {}
    uses: dict[str, set[str]] = {}
    for directory in SOURCE_DIRS:
        for path in (ROOT / directory).glob("*.F90"):
            obj = f"{path.stem}.o"
            text = "\n".join(
                line.split("!", 1)[0]
                for line in path.read_text(errors="ignore").splitlines()
            )
            for match in re.finditer(
                r"(?mi)^\s*module\s+"
                r"(?!procedure\b|subroutine\b|function\b)([a-z]\w*)",
                text,
            ):
                providers[match.group(1).lower()] = obj
            uses[obj] = {
                match.group(1).lower()
                for match in re.finditer(
                    r"(?mi)^\s*use(?:\s*,[^:]*)?\s*"
                    r"(?:::\s*)?([a-z]\w*)",
                    text,
                )
            }
    return providers, uses


def explicit_object_graph() -> dict[str, set[str]]:
    text = (ROOT / "Makefile").read_text()
    text = re.sub(r"\\\n[ \t]*", " ", text)
    graph: dict[str, set[str]] = {}
    for match in re.finditer(
        r"(?m)^(?P<targets>(?:[A-Za-z0-9_]+\.o\s*)+):(?P<deps>[^\n]*)$",
        text,
    ):
        deps = set(re.findall(r"\b[A-Za-z0-9_]+\.o\b", match.group("deps")))
        for target in match.group("targets").split():
            graph.setdefault(target, set()).update(deps)
    return graph


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

    def test_fortran_modules_use_explicit_dependencies_not_stage_chains(self):
        makefile = (ROOT / "Makefile").read_text()
        self.assertNotIn("define chain_objects", makefile)
        self.assertNotIn("$(call chain_objects", makefile)

        # Representative forward edges must be encoded by the actual module
        # relationship, without serializing unrelated objects in a stage.
        for edge in (
            "MOD_Tracer_LandPhase.o: $(TRACER_BASIC_OBJS) MOD_Tracer_Reactive.o",
            "MOD_Grid_RiverLakeFlow.o: MOD_Grid_RiverLakeTimeVars.o",
            "MOD_Vars_TimeVariables.o: MOD_Tracer_Defs.o",
        ):
            with self.subTest(edge=edge):
                self.assertIn(edge, makefile)

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

    def test_all_same_stage_module_uses_have_a_make_ordering_path(self):
        providers, uses = source_module_graph()
        graph = explicit_object_graph()
        stages = {
            name: make_variable(name)
            for name in (
                "OBJS_SHARED",
                "OBJS_MKSRFDATA",
                "OBJS_BASIC",
                "OBJS_MKINIDATA",
                "OBJS_MAIN",
                "OBJS_POST1",
                "OBJS_POST2",
                "OBJS_POST3",
            )
        }
        stages["OBJECTS_CAMA"] = make_variable("OBJECTS_CAMA", ("CaMa=YES",))

        def reachable(target: str, provider: str) -> bool:
            pending = list(graph.get(target, ()))
            seen: set[str] = set()
            while pending:
                dependency = pending.pop()
                if dependency == provider:
                    return True
                if dependency not in seen:
                    seen.add(dependency)
                    pending.extend(graph.get(dependency, ()))
            return False

        missing: list[str] = []
        for stage, objects in stages.items():
            object_set = set(objects)
            for target in objects:
                for module in uses.get(target, ()):
                    provider = providers.get(module)
                    if (
                        provider is not None
                        and provider != target
                        and provider in object_set
                        and not reachable(target, provider)
                    ):
                        missing.append(f"{stage}: {target} USEs {module} from {provider}")

        self.assertEqual(missing, [], "\n".join(missing))

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

    def test_colm_alias_tracks_the_binary_for_true_noop_builds(self):
        makefile = (ROOT / "Makefile").read_text()

        self.assertIn(".PHONY: colm.x", makefile)
        self.assertIn("colm.x : run/colm.x", makefile)
        self.assertRegex(
            makefile,
            r"(?m)^run/colm\.x\s*:[^\n]*\| mkdir_build$",
        )


if __name__ == "__main__":
    unittest.main()
