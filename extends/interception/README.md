# Extended canopy interception schemes

This directory keeps the non-default canopy interception parameterizations
(CLM4/CLM5/Noah-MP/MATSIRO/VIC/JULES/CoLM202x) outside the default `main/`
physics path.

Default `main/` keeps only the CoLM2014 interception path.  To re-enable or
port these extended schemes, do it under the `extend_interception` preprocessor
switch and keep the default CoLM2014 path unchanged.
