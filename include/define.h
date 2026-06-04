! 1. Spatial structure:
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land subgrid type classification:
!    Select one of the following options.
#undef LULC_USGS
#undef LULC_IGBP
#define LULC_IGBP_PFT
#undef LULC_IGBP_PC

! 2.1 3D Urban model (put it temporarily here):
#undef URBAN_MODEL
!    Dependence:  only LULC_IGBP subgrid type for
!    single point URBAN_MODEL right now.
#if (defined URBAN_MODEL && defined SinglePoint)
#define LULC_IGBP
#undef LULC_USGS
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
#endif

! 3. If defined, debug information is output.
#undef CoLMDEBUG
! 3.1 If defined, range of variables is checked.
#undef RangeCheck
! 3.1 If defined, surface data in vector is mapped to gridded data for checking.
#undef SrfdataDiag

! 4. If defined, MPI parallelization is enabled.
#define USEMPI
!    Conflict: not used when defined SingPoint.
#if (defined SinglePoint)
#undef USEMPI
#endif

! 5. Hydrological process options.
! 5.1 Two soil hydraulic models can be used.
#undef   Campbell_SOIL_MODEL
#define  vanGenuchten_Mualem_SOIL_MODEL
! 5.2 If defined, lateral flow is modeled.
#define CatchLateralFlow
!    Conflicts :
#ifndef CATCHMENT
#undef CatchLateralFlow
#endif

! 6. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood
#if (defined SinglePoint)
#undef CaMa_Flood
#endif
#ifndef USEMPI
#undef CaMa_Flood
#endif

#define GridRiverLakeFlow
!    Conflicts :
#if (defined CATCHMENT || defined SinglePoint)
#undef GridRiverLakeFlow
#endif

! NOTE: the former standalone river-lake sediment macro has been retired.
! Sediment is now a TRACER 'particle' species and is compiled/activated
! under #ifdef TRACER together with GridRiverLakeFlow.

! 7. If defined, BGC model is used.
#define BGC

!    Conflicts :  only used when LULC_IGBP_PFT is defined.
#ifndef LULC_IGBP_PFT
#ifndef LULC_IGBP_PC
#undef BGC
#endif
#endif
! 7.1 If defined, CROP model is used
#define CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif

! 8. If defined, open Land use and land cover change mode.
#undef LULCC

! 9. If defined, data assimilation is used.
#undef DataAssimilation
#if (defined DataAssimilation)
#define LULC_IGBP
#undef LULC_USGS
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
#endif

! 10. Interface to AI model.
#undef USESplitAI

! 11. External lake models.
#undef EXTERNAL_LAKE

! 12. Hyperspectral scheme.
#undef HYPERSPECTRAL

! 12b. If defined, extended canopy interception schemes are enabled.
#undef extend_interception

! 13. If defined, water tracer module is enabled (e.g. delta-18O, delta-D).
!     Default is OFF for production-safe builds; change to #define TRACER
!     only when water tracers / particle tracer species are explicitly needed.
#define TRACER
!    Conflicts: TRACER requires VariablySaturatedFlow soil hydrology
!    (vanGenuchten_Mualem_SOIL_MODEL). Disable when running with
!    Campbell_SOIL_MODEL.
#ifdef Campbell_SOIL_MODEL
#undef TRACER
#endif
!    Dependency: the TRACER subsystem (water isotopes + the sediment particle
!    species) routes through grid river/lake flow. Enabling TRACER REQUIRES
!    GridRiverLakeFlow.
#if (defined TRACER) && (!defined GridRiverLakeFlow)
#error "TRACER requires GridRiverLakeFlow to be defined in include/define.h"
#endif

! 13b. Methane reactive tracer.
!     Activation is runtime: register a tracer named "CH4" or "METHANE"
!     with category="reactive" in the &nl_colm DEF_TRACER_NAMES /
!     DEF_TRACER_TYPES namelist. The methane module is compiled whenever
!     both TRACER and BGC are defined; the registry resolves igas_ch4 at
!     run time and switches all methane logic on/off accordingly.
!     Additional dependency: requires LULC_IGBP_PFT or LULC_IGBP_PC for
!     pftfrac access (per-PFT NPP and root-respiration aggregation).
#if (defined TRACER) && (defined BGC)
#if (!defined LULC_IGBP_PFT && !defined LULC_IGBP_PC)
#error "Methane (TRACER+BGC) requires LULC_IGBP_PFT or LULC_IGBP_PC for pftfrac access."
#endif
#endif
