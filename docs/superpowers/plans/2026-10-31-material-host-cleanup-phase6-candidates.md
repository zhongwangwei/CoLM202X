# Further material stage consolidation candidates (read-only, not implemented)

Reuse existing MATERIAL modules; avoid new parallel facades or setter wrappers.

1. MAIN 1058–1079: bring irrigation capture, water diagnostic correction, flood
   input extraction, water transport and flood publication into the existing
   LandPhase water-stage entry. Keep the exact sequence; move the local material
   input vector inside MATERIAL and remove MAIN's direct nmts dependency if unused.
   Conservative union of all source USE edges found no reverse path from
   RiverLakeCoupling to LandPhase (56 reachable modules), hence no detected new
   cycle; this is static evidence, not a build result for that future change.
2. MAIN 1096–1105: use one SnowHost normalization entry, combine first and test
   updated snl before divide. Retain combine lb:1 vs divide lb:0 slices.
3. MAIN 972–978: absorb snowmelt carry and postthermal water snapshot into the
   existing LandPhase thermal-stage entry; preserve evap → carry → snapshot order.

These consolidate coherent material stages; they do not justify moving host
energy/water or SNICAR aerosol algorithms. Snow topology events, in-solver phase
sampling, and river update callbacks still need their host execution positions.
The three irrigation/glacier bookkeeping assignments could be wrapped individually,
but that adds indirection without removing meaningful complexity.

Freeze the verified phase-5 tree before any next implementation. Extend the
existing boundary/snow/lifecycle regressions and verify material on/off builds.
