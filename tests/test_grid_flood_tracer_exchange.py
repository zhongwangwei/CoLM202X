"""Adjoint flood-tracer credit/debit check with evaporation before infiltration."""

from pathlib import Path


FLOW = (Path(__file__).resolve().parents[1] / "main/HYDRO/MOD_Grid_RiverLakeFlow.F90").read_text()
SOIL = (Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_SoilWater.F90").read_text()
MAIN = (Path(__file__).resolve().parents[1] / "main/CoLMMAIN.F90").read_text()
CONSERVATION = (Path(__file__).resolve().parents[1] / "main/TRACER/MOD_Tracer_Conservation.F90").read_text()


def test_mixed_visible_protected_signatures_and_evaporation_first_close():
    # Two catchments, two grids, two patches crossing both grids. Units: m2,m3,R*m3.
    overlap = ((10.0, 5.0), (10.0, 15.0))
    patches = ((12.0, 8.0), (8.0, 12.0))
    visible_water = (20.0, 15.0)
    protected_water = (10.0, 10.0)
    visible_tracer = (2.0, 7.5)  # deliberately unlike protected tracer
    protected_tracer = (8.0, 0.5)
    grid_area = tuple(sum(p[g] for p in patches) for g in range(2))
    donor_area = tuple(sum(row) for row in overlap)

    def publish(donor):
        grid = tuple(
            sum(donor[d] * overlap[d][g] / donor_area[d] for d in range(2)) / grid_area[g]
            for g in range(2)
        )
        return tuple(sum(grid[g] * patches[p][g] for g in range(2)) / sum(patches[p]) for p in range(2))

    water = tuple(a + b for a, b in zip(publish(visible_water), publish(protected_water)))
    tracer = tuple(a + b for a, b in zip(publish(visible_tracer), publish(protected_tracer)))
    evaporated = tuple(0.25 * w for w in water)
    infiltrated = tuple(0.75 * w for w in water)
    patch_fraction = tuple(i / (w - e) for w, e, i in zip(water, evaporated, infiltrated))
    patch_input = sum(tracer[p] * patch_fraction[p] * sum(patches[p]) for p in range(2))
    grid_fraction = tuple(
        sum(patch_fraction[p] * patches[p][g] for p in range(2)) / grid_area[g]
        for g in range(2)
    )
    donor_fraction = tuple(
        sum(grid_fraction[g] * overlap[d][g] for g in range(2)) / donor_area[d]
        for d in range(2)
    )
    donor_debit = sum(
        (visible_tracer[d] + protected_tracer[d]) * donor_fraction[d]
        for d in range(2)
    )
    assert abs(patch_input - donor_debit) < 1e-12
    assert patch_input == sum(visible_tracer) + sum(protected_tracer)


def test_exchange_is_wired_to_land_budget_and_donor_repartition():
    assert "infiltrated_fraction = flood_infil_acc(i)/(water_credit-flood_evap_acc(i))" in FLOW
    assert "ratio_patch(i) = evaporated_fraction + &" in FLOW
    assert "trc_mass(itrc,j) - coefficient_uc(j)*flood_visible_tracer_uc(itrc,j)" in FLOW
    assert "trc_levsto(itrc,j) - coefficient_uc(j)*flood_protected_tracer_uc(itrc,j)" in FLOW
    assert "flood_tracer_land_patch(itrc,i)*patch_area" in FLOW
    assert "land tracer input differs from donor exchange" in FLOW
    assert "CALL levee_tracer_repartition(j" in FLOW
    assert "late_tracer = late_tracer + &" in SOIL
    assert "flood_tracer_input(itrc) - flood_ground_evap_tracer" in SOIL
    assert "flood_destination_water = max(wdsrf,0._r8)+late_runoff_water-late_surface_water &" in SOIL
    assert "+ max(qinfl,0._r8)*deltim" in SOIL
    assert SOIL.index("early_runoff_water * ratio") < SOIL.index("late_tracer = late_tracer + &")
    assert SOIL.index("a_trc_rsur(itrc, ipatch) = a_trc_rsur") < SOIL.index(
        "late_tracer = late_tracer + &"
    )
    assert "a_trc_precip(itrc,ipatch) = a_trc_precip(itrc,ipatch) + flood_tracer_input(itrc)" in SOIL


def test_signed_isotope_vapor_uptake_adjoint_closes_with_distinct_sources():
    # C=0 for one donor and negative vapour loss: a nonnegative-only debit
    # or division by C would fail. Geometry matches the production adjoint.
    overlap = ((10.0, 5.0), (10.0, 15.0))
    patches = ((12.0, 8.0), (8.0, 12.0))
    donor_area = tuple(sum(row) for row in overlap)
    grid_area = tuple(sum(row[g] for row in patches) for g in range(2))
    water_src = (20.0, 30.0)
    tracer_src = (0.0, 0.06)

    def publish(src):
        grid = tuple(sum(src[d] * overlap[d][g] / donor_area[d]
                         for d in range(2)) / grid_area[g] for g in range(2))
        return tuple(sum(grid[g] * patches[p][g] for g in range(2))
                     / sum(patches[p]) for p in range(2))

    water, credit = publish(water_src), publish(tracer_src)
    evap = tuple(0.3 * w for w in water)
    infil = tuple(0.4 * w for w in water)
    vapor = tuple(-0.001 * w for w in water)
    f = tuple(i / (w - e) for w, e, i in zip(water, evap, infil))
    a = tuple(max(l, 0) / c if c else 0 for l, c in zip(vapor, credit))
    coeff = tuple(x + (1 - x) * y for x, y in zip(a, f))
    gain = tuple(max(-l, 0) * (1 - y) / w for l, y, w in zip(vapor, f, water))
    land = sum((c - l) * f[p] * sum(patches[p])
               for p, (c, l) in enumerate(zip(credit, vapor)))
    vapor_total = sum(l * sum(patches[p]) for p, l in enumerate(vapor))

    def adjoint(patch):
        grid = tuple(sum(patch[p] * patches[p][g] for p in range(2))
                     / grid_area[g] for g in range(2))
        return tuple(sum(grid[g] * overlap[d][g] for g in range(2))
                     / donor_area[d] for d in range(2))

    c_donor, g_donor = adjoint(coeff), adjoint(gain)
    donor_loss = sum(tracer_src[d] * c_donor[d] - water_src[d] * g_donor[d]
                     for d in range(2))
    assert abs(donor_loss - land - vapor_total) < 1e-12
    assert tracer_src[0] + water_src[0] * g_donor[0] >= 0
    assert vapor_total < 0


def test_late_flood_mixing_does_not_contaminate_prior_runoff_and_can_pond():
    # Ordinary water is already exposed to surface runoff when the flood
    # boundary input arrives. The VSF solver can divide the late mixture
    # between ponding and actual infiltration.
    ordinary_water, ordinary_tracer = 5.0, 0.5
    runoff, flood_water, flood_tracer = 2.0, 3.0, 1.2
    runoff_tracer = ordinary_tracer * runoff / ordinary_water
    late_water = ordinary_water - runoff + flood_water
    late_tracer = ordinary_tracer - runoff_tracer + flood_tracer
    pond, infiltrated = 2.0, 4.0
    assert pond + infiltrated == late_water
    assert abs(runoff_tracer - 0.2) < 1e-12  # independent of flood signature
    assert abs(pond * late_tracer / late_water + infiltrated * late_tracer / late_water
               + runoff_tracer - ordinary_tracer - flood_tracer) < 1e-12


def test_ground_evap_deficit_uses_flood_before_net_infiltration():
    # Host qgtop adds flood after ordinary qseva: no pond/rain, 1 mm qseva
    # and 2 mm flood gives net qinfl=1 mm. The first flood mm evaporates.
    flood_water, flood_tracer = 2.0, 0.004
    ordinary_qseva, net_qinfl = 1.0, 1.0
    surface_base_balance = net_qinfl - flood_water
    flood_ground_evap = min(flood_water, max(0, -surface_base_balance), ordinary_qseva)
    tracer_evap = flood_tracer * flood_ground_evap / flood_water
    late_ratio = (flood_tracer - tracer_evap) / net_qinfl
    assert flood_ground_evap == 1.0
    assert abs(late_ratio - 0.002) < 1e-15
    assert abs(tracer_evap + late_ratio * net_qinfl - flood_tracer) < 1e-15
    assert "min(flood_water, max(0._r8,-surface_base_balance), gwat_evap)" in SOIL
    assert "flood_tracer_input(itrc) - flood_ground_evap_tracer" in SOIL


def test_negative_qinfl_splits_evaporation_between_flood_and_soil():
    flood_water, ordinary_qseva, net_qinfl = 2.0, 3.0, -1.0
    top_soil_evap = max(-net_qinfl, 0)
    gwat_evap = ordinary_qseva - top_soil_evap
    balance = net_qinfl - flood_water
    flood_ground_evap = min(flood_water, max(0, -balance), gwat_evap)
    assert (top_soil_evap, flood_ground_evap) == (1.0, 2.0)
    assert abs(top_soil_evap + flood_ground_evap - ordinary_qseva) < 1e-15


def test_flood_only_rewet_dissolves_nonvolatile_residue_after_runoff():
    # The ordinary pool is dry, so its prior runoff must not receive this
    # residue. Flood water arrives at the late qgtop/ponding partition.
    ordinary_water, ordinary_runoff = 0.0, 0.0
    flood_water, flood_tracer, dry_residue = 2.0, 0.004, 0.006
    late_water = flood_water
    late_tracer = flood_tracer + dry_residue
    ratio = late_tracer / late_water
    assert ordinary_runoff == ordinary_water == 0.0
    assert abs(ratio - 0.005) < 1e-15
    assert abs(ratio * late_water - flood_tracer - dry_residue) < 1e-15
    runoff_at = SOIL.index("a_trc_rsur(itrc, ipatch) = a_trc_rsur")
    rewet_at = SOIL.index("late_tracer = late_tracer + trc_surface_residue(itrc,ipatch)")
    partition_at = SOIL.index("late_ratio = max(late_tracer,0._r8)/late_water")
    assert runoff_at < rewet_at < partition_at


def test_feedback_disables_fixed_signature_soil_assumption_after_drydown():
    # River isotopes stay in the soil after the instant qinfl_fld returns to
    # zero. Do not re-enable R_init diagnostics on the next dry/restart step.
    assert MAIN.count("flood_heterogeneous_in = LWINFILT .and. patchtype == 0") == 2
    assert "IF (flood_heterogeneous_in) THEN" in CONSERVATION
    assert "water_corrected_check = .false." in CONSERVATION
    assert "fixed_signature_step = .false." in CONSERVATION
