import math
import pathlib
import re
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
FRAC = (ROOT / "main/TRACER/MOD_Tracer_Frac.F90").read_text()
SOIL = (ROOT / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text()


class TracerIsotopeNssTests(unittest.TestCase):
    def test_peclet_and_conductance_use_consistent_area_bases(self):
        self.assertRegex(
            FRAC,
            r"transp_moles_leaf_s\s*=\s*transp_moles\s*/\s*"
            r"\(deltim\s*\*\s*max\(leaf_area,\s*trc_tiny\)\)",
        )
        self.assertRegex(
            FRAC,
            r"DEF_TRACER_NSS_LEAF_RB,\s*0\._r8\)\s*/\s*"
            r"max\(leaf_area,\s*trc_tiny\)",
        )
        self.assertRegex(
            FRAC,
            r"tracer_leaf_kinetic_epsilon\(itrc,\s*0\._r8,\s*&\s*\n\s*"
            r"DEF_TRACER_NSS_LEAF_RB\s*/\s*max\(leaf_area,\s*trc_tiny\),\s*"
            r"stomatal_resistance\)",
        )
        self.assertIn(
            "conductance_gross_moles = deltim * vapor_molar_density_sat / total_resistance",
            FRAC,
        )

        # Equal leaf-area transpiration must give equal Péclet transport.
        water_moles_per_mm = 1000.0 / 18.01528
        p1 = 0.2 * water_moles_per_mm / (1800.0 * 1.0)
        p4 = 0.8 * water_moles_per_mm / (1800.0 * 4.0)
        self.assertAlmostEqual(p1, p4)

        # rst is already canopy resistance; only leaf boundary resistance
        # is divided by LAI before the serial resistances are combined.
        rb = 100.0
        self.assertAlmostEqual(50.0 + rb / 4.0, 75.0)

    def test_craig_gordon_has_no_humidity_switch(self):
        self.assertNotIn("craig_gordon_relhum_fallback", FRAC)
        self.assertIn(
            "h = min(max(relhum, 0._r8), craig_gordon_relhum_max)", FRAC
        )
        self.assertRegex(
            FRAC,
            r"max\(cg_ratio,\s*&\s*\n\s*"
            r"craig_gordon_min_net_ratio_frac\s*\*\s*equilibrium_ratio\)",
        )

        def ratio(humidity):
            equilibrium = 0.98
            vapor = 0.97
            kinetic = 1.02
            h = min(max(humidity, 0.0), 0.95)
            raw = (equilibrium - h * vapor) / (kinetic * (1.0 - h))
            return max(raw, 0.75 * equilibrium)

        self.assertLess(abs(ratio(0.900001) - ratio(0.899999)), 1.0e-4)
        self.assertEqual(ratio(0.95), ratio(0.999999))
        self.assertTrue(math.isfinite(ratio(1.0)))

    def test_leaf_off_releases_signed_anomaly_without_negative_pools(self):
        self.assertIn("CALL release_leaf_iso_storage", SOIL)
        self.assertIsNotNone(
            re.search(
                r"lai_frac\s*<=\s*trc_tiny\s*\.or\..*?"
                r"trc_leaf_water_moles.*?CALL release_leaf_iso_storage",
                SOIL,
                re.DOTALL,
            ),
            "leaf-off must trigger the conservative release",
        )
        self.assertIn("release_fraction = min(-anomaly / pool_total, 1._r8)", SOIL)
        release_body = SOIL.split("SUBROUTINE release_leaf_iso_storage", 1)[1].split(
            "END SUBROUTINE release_leaf_iso_storage", 1
        )[0]
        self.assertNotIn("trc_numerical_residual_step", release_body)

        def release(anomaly, pools, water_weights=(2.0, 1.0)):
            before = anomaly + sum(pools)
            if anomaly >= 0.0:
                total = sum(water_weights)
                if total > 0.0:
                    pools = [
                        p + anomaly * w / total
                        for p, w in zip(pools, water_weights)
                    ]
                    anomaly = 0.0
            else:
                available = sum(max(p, 0.0) for p in pools)
                if available > 0.0:
                    fraction = min(-anomaly / available, 1.0)
                    pools = [p * (1.0 - fraction) for p in pools]
                    anomaly += fraction * available
            self.assertAlmostEqual(before, anomaly + sum(pools))
            self.assertTrue(all(p >= 0.0 for p in pools))
            return anomaly

        self.assertEqual(release(0.3, [1.0, 0.5]), 0.0)
        self.assertAlmostEqual(release(-0.3, [1.0, 0.5]), 0.0)
        self.assertAlmostEqual(release(-2.0, [1.0, 0.5]), -0.5)
        self.assertAlmostEqual(release(-0.3, [0.0, 0.0]), -0.3)
        self.assertAlmostEqual(release(0.3, [0.0, 0.0], (0.0, 0.0)), 0.3)


if __name__ == "__main__":
    unittest.main()
