! 1. Spatial structure:
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land TYPE classification :
!    Select one of the following options.
#undef LULC_USGS
#define LULC_IGBP
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC
! 2.1 Urban model setting (put it temporarily here):
#undef URBAN_MODEL
#undef URBAN_LCZ

! 3. If defined, debug information is output.
#define CoLMDEBUG
! 3.1 If defined, range of variables is checked.
#define RangeCheck
! 3.2 If defined, surface data in vector is mapped to gridded data for checking.
#undef SrfdataDiag

! 4. If defined, MPI parallelization is enabled.
#define  USEMPI
!    Conflict: not used when defined SingPoint.
#if (defined SinglePoint)
#undef USEMPI
#endif

! 5. Hydrological process options.
! 5.1 Two soil hydraulic models can be used.
#undef Campbell_SOIL_MODEL
#define vanGenuchten_Mualem_SOIL_MODEL
! 5.2 If defined, lateral flow is modeled.
#define  LATERAL_FLOW
!    Conflicts :
#ifndef CATCHMENT
#undef LATERAL_FLOW
#endif

! 6. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood

#define GridRiverLakeFlow
!    Conflicts :
!    The river network is built per worker, so it needs MPI (as CaMa_Flood).
#ifndef USEMPI
#undef GridRiverLakeFlow
#endif
#if (defined CATCHMENT || defined SinglePoint || defined CaMa_Flood)
#undef GridRiverLakeFlow
#endif

! 7. If defined, BGC model is used.
#undef BGC

!    Conflicts : only used when LULC_IGBP_PFT or LULC_IGBP_PC is defined.
#ifndef LULC_IGBP_PFT
#ifndef LULC_IGBP_PC
#undef BGC
#endif
#endif

! 7.1 If defined, CROP model is used
#undef CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif

! 8. If defined, open Land use and land cover change mode.
#define LULCC

! 13. If defined, material transport is compiled in.
#define MaterialTransport
#if (defined MaterialTransport && defined Campbell_SOIL_MODEL)
#error "MaterialTransport requires vanGenuchten_Mualem_SOIL_MODEL"
#endif
