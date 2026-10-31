# Phase 6: consolidate existing material stages

Frozen verified phase-5 tree: /tmp/colm-host-phase5-t4jdp4xy, durable copy with
SHA256 manifest in .omx/artifacts/material-host-cleanup-phase6-20261031/baseline.
Preserve all earlier uncommitted work. Scope is the three previously identified
stages, with no new dependencies or independent facade modules.

1. Extend the existing LandPhase water-stage entry to own irrigation capture,
   diagnostic flux correction, flood material extraction and publication. Preserve
   input -> water transport -> publication order. Move temporary species array
   inside MATERIAL. Caller supplies only necessary water state/forcing inputs.
2. Add one SnowHost normalization entry composing existing combine and divide;
   check updated snl after combine, preserving array bounds and slice differences.
3. Extend existing thermal-stage entry to own snowmelt carry and postthermal water
   snapshot, preserving evaporation -> carry -> snapshot order.

Implementation ownership: parent MAIN/build wiring; LandPhase executor handles
its two stages plus new differential coverage/affected tests; SnowHost executor
handles normalization plus snow adapter regression. Freeze before edits; existing
phase5 boundary/lifecycle/snow coverage locks behavior. All production changes
batched before final tests/builds. Separate read-only reviewer after integration.
Run actual extracted stage routines against frozen boundary calls with recording
services, compare all argument bits/order; actual snow topology and adapter tests;
existing lifecycle/soil/I/O checks when affected; generated registry, MT-off CPP
comparison, complete MT-on/off + PFT/CROP builds on T7920. Long simulations only if
focused proofs or integration reveal unresolved risk. No complete-model numerical
claim unless an actual model rerun is done.
