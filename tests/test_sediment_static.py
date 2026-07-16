#!/usr/bin/env python3
"""Static regression checks for sediment MPI diagnostics and CFL safety."""

from pathlib import Path
import math
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

    def test_parameter_units_and_cfl_domain_fail_fast(self) -> None:
        reader = self.source.split("SUBROUTINE read_sediment_parameter_file", 1)[1].split(
            "END SUBROUTINE read_sediment_parameter_file", 1
        )[0]
        validator = self.source.split("SUBROUTINE validate_sediment_parameters", 1)[1].split(
            "END SUBROUTINE validate_sediment_parameters", 1
        )[0]

        self.assertIn("grain_density must be supplied in kg/m3", reader)
        self.assertIn("water_density must be supplied in kg/m3", reader)
        self.assertNotIn("IF (psedD > 100._r8)", reader)
        self.assertNotIn("IF (pwatD > 100._r8)", reader)
        self.assertIn("psedD = DEF_SEDIMENT%grain_density / 1000._r8", reader)
        self.assertIn("pwatD = DEF_SEDIMENT%water_density / 1000._r8", reader)
        self.assertIn("sed_cfl_adv <= 0._r8 .or. sed_cfl_adv > 1._r8", validator)

    def test_all_sediment_configuration_and_static_inputs_reject_nonfinite_values(self) -> None:
        reader = self.source.split("SUBROUTINE read_sediment_parameter_file", 1)[1].split(
            "END SUBROUTINE read_sediment_parameter_file", 1
        )[0]
        validator = self.source.split("SUBROUTINE validate_sediment_parameters", 1)[1].split(
            "END SUBROUTINE validate_sediment_parameters", 1
        )[0]
        initializer = self.source.split("SUBROUTINE grid_sediment_init", 1)[1].split(
            "END SUBROUTINE grid_sediment_init", 1
        )[0]
        static_reader = self.source.split("SUBROUTINE read_sediment_static_data", 1)[1].split(
            "END SUBROUTINE read_sediment_static_data", 1
        )[0]

        self.assertIn(
            "USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_is_finite",
            self.source,
        )
        scalar_fields = (
            "lambda",
            "lyrdph",
            "psedD",
            "pwatD",
            "visKin",
            "vonKar",
            "pset",
            "pyld",
            "pyldc",
            "pyldpc",
            "dsylunit",
            "sed_ignore_dph",
            "sed_cfl_adv",
            "sed_dt_max",
            "sed_bed_depth",
        )
        for field in scalar_fields:
            self.assertRegex(
                validator,
                rf"\.not\.\s*ieee_is_finite\({re.escape(field)}\)",
            )

        # Validate the namelist record before NaN can masquerade as an omitted
        # negative sentinel in the ordered override comparisons.
        for field in (
            "grain_density",
            "water_density",
            "porosity",
            "ignore_depth_m",
            "active_layer_depth",
            "viscosity",
            "von_karman",
            "settling_multiplier",
            "yield_coefficient",
            "slope_exponent",
            "precipitation_exponent",
            "unit_conversion",
            "cfl_adv",
            "max_timestep_s",
            "bed_depth",
            "grain_diameter",
        ):
            self.assertIn(f"ieee_is_finite(DEF_SEDIMENT%{field})", reader)

        self.assertIn("nsed <= 0", validator)
        self.assertIn("nlfp_sed <= 0", validator)
        for field, domain in (
            ("sDiam", "any(sDiam <= 0._r8)"),
            ("setvel", "any(setvel < 0._r8)"),
            ("sed_frc", "any(sed_frc < 0._r8)"),
            ("sed_slope", "any(sed_slope < 0._r8)"),
        ):
            self.assertIn(f"any(.not. ieee_is_finite({field}))", validator)
            self.assertIn(domain, validator)

        parse = initializer.index("CALL parse_grain_diameters()")
        settle = initializer.index("CALL calc_settling_velocities()")
        derived_validation = initializer.index(
            "CALL validate_sediment_parameters()", settle
        )
        static_read = initializer.index("CALL read_sediment_static_data", settle)
        self.assertLess(parse, settle)
        self.assertLess(settle, derived_validation)
        self.assertLess(derived_validation, static_read)
        self.assertLess(
            static_reader.index("CALL validate_sediment_parameters()"),
            static_reader.index("sed_slope = max(sed_slope, 0._r8)"),
        )
        self.assertLess(
            static_reader.index("CALL validate_sediment_parameters()"),
            static_reader.index("CALL normalize_sed_frc()"),
        )

    def test_finite_gate_closes_ordered_comparison_nan_and_infinity_holes(self) -> None:
        # Mirrors representative old guards.  NaN evades every ordered
        # comparison, while +Inf evades one-sided non-negative checks.
        old_porosity_invalid = lambda value: value < 0.0 or value >= 1.0
        old_yield_invalid = lambda value: value < 0.0

        self.assertFalse(old_porosity_invalid(math.nan))
        self.assertFalse(old_yield_invalid(math.nan))
        self.assertFalse(old_yield_invalid(math.inf))
        for value in (math.nan, math.inf, -math.inf):
            self.assertFalse(math.isfinite(value))

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
            "EXCH_SHEARVEL_BLEND = 2._r8 * EXCH_SHEARVEL_MIN", self.source
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

    def test_absolute_discharge_accumulator_has_complete_lifecycle(self) -> None:
        reader = self.source.split("SUBROUTINE read_sediment_restart", 1)[1].split(
            "END SUBROUTINE read_sediment_restart", 1
        )[0]
        writer = self.source.split("SUBROUTINE write_sediment_restart", 1)[1].split(
            "END SUBROUTINE write_sediment_restart", 1
        )[0]
        finalizer = self.source.split("SUBROUTINE grid_sediment_final", 1)[1].split(
            "END SUBROUTINE grid_sediment_final", 1
        )[0]

        self.assertIn("'sed_acc_abs_rivout'", reader)
        self.assertIn("'sed_acc_abs_rivout'", writer)
        self.assertIn("deallocate(sed_acc_abs_rivout)", finalizer.replace(" ", ""))

    def test_flow_cancellation_diagnostics_are_counted_once_per_routing_period(self) -> None:
        averaging = self.calc.split("DO WHILE (sed_time_remaining > 0._r8)", 1)[1]
        averaging = averaging.split("dt_cfl_local = dt_morph", 1)[0]
        self.assertIn("IF (iter_sed == 1) THEN", averaging)
        counted = averaging.split("IF (iter_sed == 1) THEN", 1)[1].split("ENDIF", 1)[0]
        self.assertIn("sum_rivout_signed_local", counted)
        self.assertIn("sum_rivout_abs_local", counted)
        self.assertIn("n_flow_cancel_local", counted)

    def test_cfl_cap_deposition_is_integrated_only_in_its_own_substep(self) -> None:
        advection = self.source.split("SUBROUTINE calc_sediment_advection", 1)[1].split(
            "END SUBROUTINE calc_sediment_advection", 1
        )[0]
        accumulator = self.source.split("SUBROUTINE accumulate_sediment_output", 1)[
            1
        ].split("END SUBROUTINE accumulate_sediment_output", 1)[0]

        self.assertIn("netflw_adv_step(:,:) = 0._r8", advection)
        self.assertIn("exch_d_adv_step(:,:) = 0._r8", advection)
        self.assertIn(
            "netflw_adv_step(:,i) = netflw_adv_step(:,i) - dTmp(:) / dt",
            advection,
        )
        self.assertIn(
            "netflw_adv_step(:,i) = netflw_adv_step(:,i) - sedsto(:,i) / dt",
            advection,
        )
        self.assertNotIn("netflw(:,i) = netflw(:,i) - dTmp(:) / dt", advection)
        self.assertNotIn("netflw(:,i) = netflw(:,i) - sedsto(:,i) / dt", advection)
        self.assertRegex(
            accumulator.replace("&", ""),
            r"a_netflw\(:,i\)\s*=\s*a_netflw\(:,i\)\s*\+\s*"
            r"\(netflw\(:,i\)\s*\+\s*netflw_adv_step\(:,i\)\)\s*\*\s*dt",
        )

        # A morphology-base rate is integrated for the whole interval, whereas
        # each cap/dry deposit is a one-substep mass.  The result must not depend
        # on the number of CFL partitions.
        total_dt = 8.0
        base_rate = 0.25
        total_deposit = 0.8
        expected = base_rate * total_dt - total_deposit
        for nsub in (1, 2, 8):
            sub_dt = total_dt / nsub
            deposits = [total_deposit / nsub] * nsub
            integrated = sum(
                (base_rate - deposit / sub_dt) * sub_dt for deposit in deposits
            )
            self.assertAlmostEqual(integrated, expected, places=14)

    def test_suspended_solid_volume_is_the_single_canonical_state(self) -> None:
        begin = self.source.split("SUBROUTINE begin_suspended_period", 1)[1].split(
            "END SUBROUTINE begin_suspended_period", 1
        )[0]
        advection = self.source.split("SUBROUTINE calc_sediment_advection", 1)[
            1
        ].split("END SUBROUTINE calc_sediment_advection", 1)[0]
        exchange = self.source.split("SUBROUTINE calc_sediment_exchange", 1)[
            1
        ].split("END SUBROUTINE calc_sediment_exchange", 1)[0]
        sediment_input = self.source.split("SUBROUTINE apply_sediment_input", 1)[
            1
        ].split("END SUBROUTINE apply_sediment_input", 1)[0]

        self.assertIn("Suspended solid volume [m3", self.source)
        self.assertIn("sedcon(:,i) = sedsto(:,i) / rivsto(i)", begin)
        self.assertNotIn("sedsto(:,i) = sedcon(:,i)", begin)
        self.assertNotIn("allocatable :: sedsto", advection)
        self.assertNotIn("allocate(sedsto", advection)
        self.assertNotRegex(advection, r"sedsto\(:,i\)\s*=\s*sedcon\(:,i\)")
        self.assertIn("sedsto(ised,i) = sedsto(ised,i) - sedout", advection)
        self.assertIn("sedsto(ised,i) = sedsto(ised,i) + netflw", exchange)
        self.assertIn("sedsto(:,i) = sedsto(:,i) + sedinp(:,i) * dt", sediment_input)
        self.assertIn("sum_sedsto_local = sum(sedsto)", self.calc)

    def test_clean_water_volume_changes_preserve_canonical_suspended_mass(self) -> None:
        mass_by_class = (0.2, 0.3, 0.5)
        total_mass = sum(mass_by_class)

        # Pure expansion and shrinkage change only concentration.  Since mass is
        # canonical, converting to concentration never feeds back into storage.
        for current_water in (0.5, 2.0, 8.0):
            concentration = tuple(mass / current_water for mass in mass_by_class)
            self.assertAlmostEqual(
                sum(c * current_water for c in concentration), total_mass, places=14
            )

        # A dry carrier deposits the old solid volume exactly once regardless of
        # how many CFL partitions the morphology interval uses.
        porosity = 0.4
        for nsub in (1, 2, 8):
            suspended = list(mass_by_class)
            bed_bulk = [0.0] * len(suspended)
            credited = 0.0
            for _ in range(nsub):
                deposit = suspended
                credited += sum(deposit)
                bed_bulk = [b + d / (1.0 - porosity) for b, d in zip(bed_bulk, deposit)]
                suspended = [0.0] * len(suspended)
            self.assertAlmostEqual(credited, total_mass, places=14)
            self.assertAlmostEqual(
                (1.0 - porosity) * sum(bed_bulk), total_mass, places=14
            )

    def test_restart_transaction_persists_mass_and_locks_physics_descriptor(self) -> None:
        reader = self.source.split("SUBROUTINE read_sediment_restart", 1)[1].split(
            "END SUBROUTINE read_sediment_restart", 1
        )[0]
        writer = self.source.split("SUBROUTINE write_sediment_restart(", 1)[1].split(
            "END SUBROUTINE write_sediment_restart", 1
        )[0]

        start = writer.index("'sed_restart_complete_meta', 0._r8")
        mass = writer.index("'sedsto_'")
        finish = writer.rindex("'sed_restart_complete_meta'")
        self.assertLess(start, mass)
        self.assertLess(mass, finish)
        self.assertLess(writer.index("validate_sediment_checkpoint_state('write')"), start)
        self.assertIn("validate_sediment_checkpoint_state('read')", reader)
        self.assertIn("'sed_restart_schema_meta'", reader)
        self.assertIn("'sed_restart_complete_meta', 1._r8", reader)
        self.assertIn("'sedsto_' // trim(cised)", reader)

        locked = (
            "sed_lambda_meta",
            "sed_lyrdph_meta",
            "sed_psedd_meta",
            "sed_pwatd_meta",
            "sed_viskin_meta",
            "sed_vonkar_meta",
            "sed_pset_meta",
            "sed_ignore_dph_meta",
            "sed_cfl_adv_meta",
            "sed_dt_max_meta",
            "sed_pyld_meta",
            "sed_pyldc_meta",
            "sed_pyldpc_meta",
            "sed_dsylunit_meta",
            "sed_max_conc_meta",
            "sed_bedload_coeff_meta",
            "sed_precip_threshold_meta",
            "sed_diam_meta_",
            "sed_setvel_meta_",
            "sed_frc_meta_",
            "sed_slope_meta_",
            "sed_rivwth_meta",
            "sed_rivlen_meta",
        )
        for field in locked:
            self.assertIn(field, reader)
            self.assertIn(field, writer)

        self.assertIn("legacy sediment restart has nonzero concentration", reader)
        self.assertIn("legacy_state_nonzero", reader)
        self.assertIn("migrated provably empty legacy state", reader)
        self.assertIn("partial sediment transaction has canonical fields", reader)

    def test_all_effective_deposition_paths_share_one_diagnostic_definition(self) -> None:
        sediment_input = self.source.split("SUBROUTINE apply_sediment_input", 1)[
            1
        ].split("END SUBROUTINE apply_sediment_input", 1)[0]

        self.assertIn(
            "exch_d_eff(:,i) = exch_d_eff(:,i) + dTmp(:) / dt", sediment_input
        )
        self.assertIn(
            "exch_d_eff(:,i) = exch_d_eff(:,i) + sedinp(:,i)", sediment_input
        )
        self.assertIn(
            "d_eff_seen = d_eff_seen .or.", self.calc
        )
        self.assertIn(
            "sum(exch_d_eff + exch_d_adv_step, dim=1) > 0._r8", self.calc
        )
        self.assertIn("n_d_eff_local = count(d_eff_seen)", self.calc)
        readme = (SOURCE.parent / "README.md").read_text(encoding="utf-8")
        readme_compact = " ".join(readme.split())
        self.assertIn("external hillslope source", readme_compact)
        self.assertIn("effective bed increment", readme_compact)

    def test_every_advected_wet_cell_participates_in_cfl(self) -> None:
        cfl = self.calc.split("dt_cfl_local = dt_morph", 1)[1].split(
            "#ifdef USEMPI", 1
        )[0]
        self.assertIn("IF (rivsto(i) <= 0._r8) CYCLE", cfl)
        self.assertNotIn("sed_ignore_dph", cfl)
        self.assertIn("dt_cell = sed_cfl_adv * rivsto(i) / abs(rivout(i))", cfl)

    def test_initial_named_bed_depth_cannot_silently_exceed_configuration(self) -> None:
        validator = self.source.split("SUBROUTINE validate_sediment_parameters", 1)[
            1
        ].split("END SUBROUTINE validate_sediment_parameters", 1)[0]
        self.assertIn(
            "sed_bed_depth < lyrdph * real(totlyrnum, r8)", validator
        )


if __name__ == "__main__":
    unittest.main()
