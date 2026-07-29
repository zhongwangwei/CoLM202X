#!/usr/bin/env python3
"""Static regression checks for tracer accounting and runtime limits."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
TRACER = ROOT / "main/TRACER"


class TracerMiscStaticChecks(unittest.TestCase):
    def test_special_patch_evaporation_books_liquid_and_ice_separately(self) -> None:
        source = (TRACER / "MOD_Tracer_SpecialPatches.F90").read_text(encoding="utf-8")
        for name in ("tracer_glacier_patch", "tracer_waterbody_patch"):
            routine = source.split(f"SUBROUTINE {name}", 1)[1].split(
                f"END SUBROUTINE {name}", 1
            )[0]
            self.assertRegex(
                routine,
                r"CALL\s+tracer_book_evap_loss\s*\(\s*itrc\s*,\s*ipatch\s*,\s*"
                r"trc_evap_liq\s*,\s*evap_liq_mass\s*,\s*TRC_EVAP_KIND_SOILEVAP\s*\)",
            )
            self.assertRegex(
                routine,
                r"CALL\s+tracer_book_evap_loss\s*\(\s*itrc\s*,\s*ipatch\s*,\s*"
                r"trc_evap_ice\s*,\s*evap_ice_mass\s*,\s*TRC_EVAP_KIND_SUBL\s*\)",
            )
            self.assertIn(
                "trc_evap_liq = tracer_atmospheric_tracer_loss",
                routine,
            )
            self.assertIn(
                "trc_evap_ice = tracer_atmospheric_tracer_loss(trc_after_liq, water_after_liq",
                " ".join(routine.replace("&", " ").split()),
            )
        self.assertNotIn("a_trc_evap", source)

    def test_glacier_consolidates_layered_snow_before_balance_snapshot(self) -> None:
        source = (TRACER / "MOD_Tracer_SpecialPatches.F90").read_text(encoding="utf-8")
        routine = source.split("SUBROUTINE tracer_glacier_patch", 1)[1].split(
            "END SUBROUTINE tracer_glacier_patch", 1
        )[0]
        snapshot = routine.index("CALL tracer_save_storage")
        consolidate = routine.index(
            "trc_scv(itrc, ipatch) = trc_scv(itrc, ipatch) +"
        )
        self.assertLess(consolidate, snapshot)
        self.assertIn(
            "sum(trc_wliq_soisno(itrc, maxsnl+1:0, ipatch))",
            routine[consolidate:snapshot],
        )
        self.assertIn(
            "sum(trc_wice_soisno(itrc, maxsnl+1:0, ipatch))",
            routine[consolidate:snapshot],
        )

    def test_soil_water_interface_omits_retired_raw_et_demand(self) -> None:
        source = (TRACER / "MOD_Tracer_SoilWater.F90").read_text(encoding="utf-8")
        self.assertIn("ieee_is_finite", source)
        self.assertIn(".not. ieee_is_finite(qcharge_eff)", source)
        signature = source.split("SUBROUTINE tracer_soil_water", 1)[1].split(
            "IMPLICIT NONE", 1
        )[0]
        self.assertIsNone(re.search(r"\betroot\b", signature, re.I))
        self.assertIn("etroot_actual", signature)
        self.assertIn("etroot_aquifer", signature)

    def test_layer_transport_uses_one_pretransport_ratio_snapshot(self) -> None:
        source = (TRACER / "MOD_Tracer_SoilWater.F90").read_text(encoding="utf-8")
        transport = source.split("! --- qlayer(j): layer j", 1)[1].split(
            "! 3. Groundwater", 1
        )[0]
        self.assertIn(
            "layer_transport_ratio(j) = current_liq_ratio(j)", transport
        )
        self.assertIn("ratio_src = layer_transport_ratio(j)", transport)
        self.assertIn("ratio_src = layer_transport_ratio(j+1)", transport)
        groundwater = source.split("! 3. Groundwater", 1)[1].split(
            "! 4. Subsurface runoff", 1
        )[0]
        self.assertIn("ratio_src = layer_transport_ratio(j)", groundwater)

        # A pre-step snapshot prevents fresh layer-1 inflow from crossing a
        # second interface in the same tracer update.
        ratios = [1.0, 0.0, 0.0]
        transfer_12 = 0.5 * ratios[0]
        transfer_23 = 0.5 * ratios[1]
        self.assertEqual(transfer_12, 0.5)
        self.assertEqual(transfer_23, 0.0)


    def test_colmmain_soil_water_calls_match_retired_etroot_interface(self) -> None:
        main = (ROOT / "main/CoLMMAIN.F90").read_text(encoding="utf-8")
        soil_water = (TRACER / "MOD_Tracer_SoilWater.F90").read_text(encoding="utf-8")
        signature = soil_water.split("SUBROUTINE tracer_soil_water", 1)[1].split(
            "IMPLICIT NONE", 1
        )[0]
        self.assertIn("etroot_actual", signature)
        self.assertIn("etroot_aquifer", signature)
        self.assertIsNone(re.search(r"\betroot\b", signature, re.I))

        # Delimit each call by balancing its parentheses rather than by keying
        # off whichever keyword argument happens to come last -- that anchor
        # breaks every time the argument list is extended.
        def call_text(text: str, start: int) -> str:
            depth = 0
            for i in range(text.index("(", start), len(text)):
                if text[i] == "(":
                    depth += 1
                elif text[i] == ")":
                    depth -= 1
                    if depth == 0:
                        return text[start : i + 1]
            raise AssertionError("unbalanced CALL tracer_soil_water argument list")

        calls = []
        offset = 0
        while True:
            start = main.find("CALL tracer_soil_water", offset)
            if start < 0:
                break
            call = call_text(main, start)
            calls.append(call)
            offset = start + len(call)
        self.assertEqual(len(calls), 2)
        for call in calls:
            self.assertIn("etroot_actual_trc", call)
            self.assertIn("etroot_aquifer_trc", call)
            self.assertNotIn("etroot_trc", call)

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

    def test_balance_uses_scale_aware_absolute_plus_relative_tolerance(self) -> None:
        source = (TRACER / "MOD_Tracer_Conservation.F90").read_text(encoding="utf-8")
        self.assertIn("trc_balance_abs_tol = 1.0e-12_r8", source)
        self.assertIn("trc_balance_rel_tol = 1.0e-10_r8", source)
        self.assertIn(
            "balance_tol = trc_balance_abs_tol + trc_balance_rel_tol * balance_scale",
            source,
        )
        self.assertNotIn("balance_scale = max(1._r8", source)

    def test_precipitation_and_runoff_tracer_signatures_are_written(self) -> None:
        variables = (TRACER / "MOD_Tracer_Vars.F90").read_text(encoding="utf-8")
        conservation = (TRACER / "MOD_Tracer_Conservation.F90").read_text(
            encoding="utf-8"
        )
        history = (TRACER / "MOD_Tracer_Hist.F90").read_text(encoding="utf-8")
        for name in ("a_water_precip", "a_water_rnof"):
            self.assertIn(f"allocate({name}", variables)
            self.assertIn(f"{name}(itrc, ipatch) = {name}(itrc, ipatch) +", conservation)
        for name in (
            "f_trc_delta_precip_",
            "f_trc_delta_runoff_",
            "f_trc_conc_precip_",
            "f_trc_conc_runoff_",
        ):
            self.assertIn(name, history)
        self.assertIn("a_trc_precip(itrc_loc, :), a_water_precip(itrc_loc, :)", history)
        self.assertIn("a_trc_rnof(itrc_loc, :), a_water_rnof(itrc_loc, :)", history)

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
