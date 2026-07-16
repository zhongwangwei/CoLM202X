#!/usr/bin/env python3
"""Static regression checks for tracer accounting and runtime limits."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main/TRACER"


class TracerMiscStaticChecks(unittest.TestCase):
    def test_special_patch_evaporation_books_paired_water(self) -> None:
        source = (TRACER / "MOD_Tracer_SpecialPatches.F90").read_text(encoding="utf-8")
        calls = re.findall(
            r"CALL\s+tracer_book_evap_loss\s*\(\s*itrc\s*,\s*ipatch\s*,\s*"
            r"trc_evap\s*,\s*evap_mass\s*\)",
            source,
            re.I,
        )
        self.assertEqual(len(calls), 2)
        self.assertNotIn("a_trc_evap", source)

    def test_unused_phase1_sync_is_removed(self) -> None:
        sources = "\n".join(
            path.read_text(encoding="utf-8") for path in TRACER.glob("*.F90")
        )
        self.assertNotIn("sync_tracer_patch_phase1", sources)

    def test_abort_bad_point_limits_use_namelist_values(self) -> None:
        namelist = (ROOT / "share/MOD_Namelist.F90").read_text(encoding="utf-8")
        conservation = (TRACER / "MOD_Tracer_Conservation.F90").read_text(encoding="utf-8")
        variables = (TRACER / "MOD_Tracer_Vars.F90").read_text(encoding="utf-8")
        for name in (
            "DEF_TRACER_BALANCE_ABORT_NBAD",
            "DEF_TRACER_RESID_ABORT_NBAD",
            "DEF_TRACER_LULCC_ABORT_NBAD",
        ):
            self.assertGreaterEqual(namelist.count(name), 4)
        self.assertIn("nbad_total > DEF_TRACER_BALANCE_ABORT_NBAD", conservation)
        self.assertIn("resid_hard_nbad_total > DEF_TRACER_RESID_ABORT_NBAD", conservation)
        self.assertIn("nbad > DEF_TRACER_LULCC_ABORT_NBAD", variables)

    def test_cama_average_handles_zero_accumulated_duration(self) -> None:
        source = (
            ROOT / "extends/CaMa/src/cmf_ctrl_tracer_mod.F90"
        ).read_text(encoding="utf-8")
        routine = source.split("SUBROUTINE CMF_TRACER_DIAG_GETAVE", 1)[1].split(
            "END SUBROUTINE CMF_TRACER_DIAG_GETAVE", 1
        )[0]
        guard = routine.index("IF (NADD_out <= 0._JPRB) THEN")
        divide = routine.index("/ REAL(NADD_out,KIND=JPRB)")
        self.assertLess(guard, divide)
        self.assertIn("D2TRCOUT_oAVG(:,:)  = 0._JPRB", routine)
        self.assertIn("D2TRCDNS_oAVG(:,:)  = 0._JPRB", routine)
        self.assertIn("D2TRCPOUT_oAVG(:,:) = 0._JPRB", routine)
        self.assertIn("RETURN", routine[guard:divide])

    def test_tracer_history_writers_keep_only_trailing_field_barrier(self) -> None:
        source = (TRACER / "MOD_Tracer_Hist.F90").read_text(encoding="utf-8")
        for name in (
            "write_history_tracer_vector_2d",
            "write_history_tracer_ratio_vector_3d",
        ):
            routine = source.split(f"SUBROUTINE {name}", 1)[1].split(
                f"END SUBROUTINE {name}", 1
            )[0]
            barriers = list(re.finditer(r"(?i)\bmpi_barrier\b", routine))
            self.assertEqual(len(barriers), 1)
            self.assertGreater(barriers[0].start(), routine.index("mpi_recv"))

    def test_annavg_finrw_has_its_own_valid_sample_counter(self) -> None:
        accflux = (
            TRACER / "MOD_Tracer_Reactive_Methane_AccFlux.F90"
        ).read_text(encoding="utf-8")
        history = (
            TRACER / "MOD_Tracer_Reactive_Methane_Hist.F90"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "CALL acc_count1d (annavg_finrw, a_annavg_finrw_acc_num)",
            accflux,
        )
        self.assertIn("allocate (a_annavg_finrw_acc_num", accflux)
        self.assertIn("a_annavg_finrw_acc_num   (:) = 0._r8", accflux)
        self.assertIn("deallocate (a_annavg_finrw_acc_num", accflux)
        self.assertEqual(accflux.count("'ch4_a_annavg_finrw_acc_num'"), 2)
        finrw_write = history.split("mhist_on('f_annavg_finrw')", 1)[1].split(
            "CALL write_history_variable_2d", 1
        )[0]
        self.assertIn("acc_num=a_annavg_finrw_acc_num", finrw_write)
        self.assertNotIn("acc_num=a_methane_acc_num_extra", finrw_write)


if __name__ == "__main__":
    unittest.main()
