#!/usr/bin/env python3
"""Static regression checks for sediment MPI diagnostics and CFL safety."""

from pathlib import Path
import re
import unittest


SOURCE = (
    Path(__file__).resolve().parents[1]
    / "main/TRACER/MOD_Tracer_Particle_Sediment.F90"
)


class SedimentStaticChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = SOURCE.read_text(encoding="utf-8")
        cls.calc = cls.source.split("SUBROUTINE grid_sediment_calc", 1)[1].split(
            "END SUBROUTINE grid_sediment_calc", 1
        )[0]

    def test_diagnostic_reductions_are_batched(self) -> None:
        self.assertEqual(len(re.findall(r"CALL\s+mpi_allreduce", self.calc, re.I)), 5)
        self.assertEqual(self.calc.count("MPI_MIN"), 1)
        self.assertEqual(self.calc.count("MPI_MAX"), 2)
        self.assertEqual(self.calc.count("MPI_SUM"), 2)

    def test_all_diagnostic_values_are_reduced(self) -> None:
        expected_groups = {
            "precip_diag_global": [
                "max_sed_precip_local",
                "max_precip_rate_local",
                "max_slope_local",
            ],
            "diag_max_global": [
                "max_sedcon_local",
                "max_sedout_local",
                "max_bedout_local",
                "max_sedinp_local",
                "max_netflw_local",
                "max_shearvel_local",
                "max_flow_cancel_local",
                "max_es_raw_local",
                "max_d_raw_local",
                "max_es_eff_local",
                "max_d_eff_local",
            ],
            "diag_sum_global": [
                "sum_layer_local",
                "sum_seddep_local",
                "sum_sedsto_local",
                "sum_sedinp_local",
                "sum_sedout_down_local",
                "sum_sedout_up_local",
                "sum_sedout_abs_local",
                "sum_netflw_pos_local",
                "sum_netflw_neg_local",
                "sum_bed_bulk_local",
                "sum_bed_solid_local",
                "sum_rivout_signed_local",
                "sum_rivout_abs_local",
                "sum_es_raw_local",
                "sum_d_raw_local",
                "sum_es_eff_local",
                "sum_d_eff_local",
            ],
            "diag_count_global": [
                "n_wet_local",
                "n_shallow_local",
                "n_source_local",
                "n_susp_local",
                "n_bed_local",
                "n_exchange_pos_local",
                "n_exchange_neg_local",
                "n_es_raw_local",
                "n_d_raw_local",
                "n_es_eff_local",
                "n_d_eff_local",
                "n_flow_cancel_local",
            ],
        }
        for array, expected in expected_groups.items():
            match = re.search(rf"{array}\s*=\s*\(/(.*?)/\)", self.calc, re.S)
            self.assertIsNotNone(match, array)
            values = re.findall(r"\b\w+_(?:local)\b", match.group(1))
            self.assertEqual(values, expected, array)

    def test_cfl_is_never_raised_and_pathological_work_fails_fast(self) -> None:
        self.assertNotIn("SED_MIN_ADV_DT", self.source)
        self.assertNotRegex(self.calc, r"dt_cfl_global\s*=\s*min\(")
        self.assertIn("SED_MAX_ADV_SUBSTEPS = 100000", self.source)
        self.assertIn(
            "dt_cfl_global < dt_morph / real(SED_MAX_ADV_SUBSTEPS, r8)",
            self.calc,
        )
        self.assertIn(
            "CALL CoLM_stop('sediment advection CFL requires too many substeps')",
            self.calc,
        )
        self.assertIn(".not. (dt_cfl_global > 0._r8)", self.calc)
        self.assertIn(
            "CALL CoLM_stop('sediment advection CFL timestep must be valid and positive')",
            self.calc,
        )

    def test_rouse_profile_uses_near_bed_concentration_with_continuous_blend(self) -> None:
        exchange = self.source.split("SUBROUTINE calc_sediment_exchange", 1)[1].split(
            "END SUBROUTINE calc_sediment_exchange", 1
        )[0]
        exchange_no_continuations = exchange.replace("&", "")

        self.assertRegex(
            exchange,
            r"profile_factor\s*=\s*Zd\(ised\)\s*/\s*\(1\._r8\s*-\s*exp\(-Zd\(ised\)\)\)",
        )
        self.assertNotRegex(
            exchange,
            r"rouse_factor\s*=\s*\(1\._r8\s*-\s*exp\(-Zd\(ised\)\)\)\s*/\s*Zd\(ised\)",
        )
        self.assertIn(
            "EXCH_SHEARVEL_BLEND = 2._r8 * EXCH_SHEARVEL_MIN", exchange
        )
        self.assertRegex(
            exchange,
            r"IF\s*\(shearvel\(i\)\s*<=\s*EXCH_SHEARVEL_MIN\)\s*THEN\s*"
            r"rouse_factor\s*=\s*1\._r8",
        )
        self.assertRegex(
            exchange_no_continuations,
            r"transition_weight\s*=\s*min\(max\(\(shearvel\(i\)\s*-\s*"
            r"EXCH_SHEARVEL_MIN\)\s*/\s*\(EXCH_SHEARVEL_BLEND\s*-\s*"
            r"EXCH_SHEARVEL_MIN\),\s*0\._r8\),\s*1\._r8\)",
        )
        self.assertRegex(
            exchange_no_continuations,
            r"transition_weight\s*=\s*transition_weight\s*\*\s*transition_weight\s*\*\s*"
            r"\(3\._r8\s*-\s*2\._r8\s*\*\s*transition_weight\)",
        )
        self.assertRegex(
            exchange,
            r"rouse_factor\s*=\s*1\._r8\s*\+\s*transition_weight\s*\*\s*"
            r"\(profile_factor\s*-\s*1\._r8\)",
        )

    def test_bedload_solid_flux_is_converted_at_bulk_layer_boundary(self) -> None:
        advection = self.source.split("SUBROUTINE calc_sediment_advection", 1)[1].split(
            "END SUBROUTINE calc_sediment_advection", 1
        )[0]

        self.assertIn("Bedload solid-volume flux", advection)
        self.assertIn(
            "'bedload solid-volume flux, size class '", self.source
        )
        self.assertRegex(
            advection,
            r"bedout\(ised,i\)\s*=\s*min\(bedout\(ised,i\),\s*\(1\._r8\s*-\s*lambda\)\s*"
            r"\*\s*layer\(ised,i\)\s*/\s*dt\)",
        )
        self.assertRegex(
            advection,
            r"avail_bed_solid\(ised,i\)\s*=\s*max\(\(1\._r8\s*-\s*lambda\)\s*\*\s*"
            r"layer\(ised,i\)\s*-\s*max\(bedout\(ised,i\),\s*0\._r8\)\s*\*\s*dt,\s*0\._r8\)",
        )
        self.assertIn(
            "CALL limit_reverse_flux(bedout, avail_bed_solid, dt)", advection
        )
        self.assertRegex(
            advection,
            r"layer\(ised,i\)\s*=\s*layer\(ised,i\)\s*\+\s*"
            r"\(-bedout\(ised,i\)\s*\+\s*bed_ups\(ised,i\)\)\s*\*\s*dt\s*/\s*"
            r"\(1\._r8\s*-\s*lambda\)",
        )


if __name__ == "__main__":
    unittest.main()
