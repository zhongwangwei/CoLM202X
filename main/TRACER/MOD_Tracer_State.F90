#include <define.h>

MODULE MOD_Tracer_Soil_Erosion

USE MOD_Precision, only: r8
USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, tracer_info_type
USE MOD_Namelist, only: DEF_Tracer_Number
USE MOD_SPMD_Task, only: p_is_master
USE MOD_LandPatch, only:  numpatch 
USE MOD_Vars_Global, only: nl_soil 

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC :: Initialize_Soil_Erosion_Parameters
PUBLIC :: Calculate_Soil_Erosion_Tracers

! Parameters for a chosen erosion model could be stored here
! e.g., real(r8) :: critical_shear_stress_coeff, etc.

CONTAINS

SUBROUTINE Initialize_Soil_Erosion_Parameters()
  ! This subroutine would read any specific parameters needed for the
  ! chosen soil erosion model from a namelist or set defaults.
  IF (p_is_master) THEN
    WRITE(*,*) "MOD_Tracer_Soil_Erosion: Initializing parameters (conceptual)."
  ENDIF
  ! Example: Read erodibility base parameters, critical slope, etc.
END SUBROUTINE Initialize_Soil_Erosion_Parameters

SUBROUTINE Calculate_Soil_Erosion_Tracers( &
    patch_idx, &
    rainfall_intensity_mms, & ! Average intensity over dt, e.g., mm/s
    overland_flow_specific_discharge_ms, & ! q_overland (m^3/s / m_width = m^2/s) or depth*velocity (m/s)
    slope_decimal, &           ! e.g., m/m
    cover_management_factor, & ! C_USLE (unitless, 0-1)
    practice_factor, &         ! P_USLE (unitless, 0-1)
    dt_seconds, &              ! Timestep (s)
    eroded_sand_flux_kg_m2_s, &   ! Output: kg/m2/s
    eroded_silt_flux_kg_m2_s, &   ! Output: kg/m2/s
    eroded_clay_flux_kg_m2_s )    ! Output: kg/m2/s

  integer, intent(in) :: patch_idx
  real(r8), intent(in) :: rainfall_intensity_mms 
  real(r8), intent(in) :: overland_flow_specific_discharge_ms 
  real(r8), intent(in) :: slope_decimal 
  real(r8), intent(in) :: cover_management_factor 
  real(r8), intent(in) :: practice_factor 
  real(r8), intent(in) :: dt_seconds
  real(r8), intent(out) :: eroded_sand_flux_kg_m2_s
  real(r8), intent(out) :: eroded_silt_flux_kg_m2_s
  real(r8), intent(out) :: eroded_clay_flux_kg_m2_s

  ! --- Local Variables ---
  real(r8) :: R_factor_val       ! Rainfall erosivity factor for the event
  real(r8) :: K_factor_val       ! Soil erodibility factor
  real(r8) :: LS_factor_val      ! Slope length and steepness factor
  real(r8) :: C_factor_val       ! Cover and management factor (passed in)
  real(r8) :: P_factor_val       ! Support practice factor (passed in)
  
  real(r8) :: total_soil_loss_kg_m2_event ! Total soil mass eroded per unit area in this event (kg/m2)
  real(r8) :: top_soil_active_depth_m     = 0.01_r8 ! Active depth for erosion (e.g., 1 cm)
  real(r8) :: soil_bulk_density_kg_m3     = 1300.0_r8 ! Assumed/average bulk density

  integer :: sand_idx, silt_idx, clay_idx, i
  real(r8) :: sand_mass_frac, silt_mass_frac, clay_mass_frac
  real(r8) :: sand_frac, silt_frac, clay_frac ! Mass fractions of available sediment
  real(r8) :: available_sand_kg_m2, available_silt_kg_m2, available_clay_kg_m2
  real(r8) :: eroded_sand_kg_m2, eroded_silt_kg_m2, eroded_clay_kg_m2
  real(r8) :: total_available_soil_texture_kg_m2

  eroded_sand_flux_kg_m2_s = 0.0_r8
  eroded_silt_flux_kg_m2_s = 0.0_r8
  eroded_clay_flux_kg_m2_s = 0.0_r8

  IF (.NOT. allocated(DEF_Tracers) .OR. DEF_Tracer_Number == 0 .OR. nl_soil == 0) RETURN
  IF (.NOT. allocated(tracer_soil_concentration)) RETURN
  IF (patch_idx > numpatch .OR. patch_idx <= 0) RETURN ! Basic check

  ! 1. Identify indices for sand, silt, clay tracers
  sand_idx = -1; silt_idx = -1; clay_idx = -1
  DO i = 1, DEF_Tracer_Number
     IF (TRIM(ADJUSTL(DEF_Tracers(i)%name)) == 'sand') sand_idx = i
     IF (TRIM(ADJUSTL(DEF_Tracers(i)%name)) == 'silt') silt_idx = i
     IF (TRIM(ADJUSTL(DEF_Tracers(i)%name)) == 'clay') clay_idx = i
  ENDDO
  IF (sand_idx == -1 .OR. silt_idx == -1 .OR. clay_idx == -1) THEN
     IF (p_is_master) WRITE(*,*) "Warning: Sand, silt, or clay tracer not defined. Cannot calculate erosion for patch ", patch_idx
     RETURN
  ENDIF

  ! 2. Calculate Erosion Factors (Conceptual - placeholders for actual models)
  R_factor_val = rainfall_intensity_mms * dt_seconds * 0.01_r8 
                                                             
  sand_mass_frac = tracer_soil_concentration(sand_idx, patch_idx, 1) / soil_bulk_density_kg_m3
  silt_mass_frac = tracer_soil_concentration(silt_idx, patch_idx, 1) / soil_bulk_density_kg_m3
  clay_mass_frac = tracer_soil_concentration(clay_idx, patch_idx, 1) / soil_bulk_density_kg_m3
  K_factor_val = 0.03_r8 

  LS_factor_val = (slope_decimal / 0.09_r8)**1.4_r8 * (10.0_r8 / 22.13_r8)**0.4_r8 
  LS_factor_val = MAX(0.1_r8, MIN(LS_factor_val, 10.0_r8)) 

  C_factor_val = cover_management_factor
  P_factor_val = practice_factor

  ! 3. Calculate Total Potential Soil Loss (Conceptual)
  IF (overland_flow_specific_discharge_ms > 1.0e-9_r8) THEN
      total_soil_loss_kg_m2_event = overland_flow_specific_discharge_ms * slope_decimal * K_factor_val * C_factor_val * P_factor_val * dt_seconds * 500.0_r8 
  ELSE
      total_soil_loss_kg_m2_event = 0.0_r8
  ENDIF
  total_soil_loss_kg_m2_event = MAX(0.0_r8, total_soil_loss_kg_m2_event)

  ! 4. Apportion to Tracers and Update Concentrations
  IF (total_soil_loss_kg_m2_event > 1.0e-12_r8) THEN
      available_sand_kg_m2 = tracer_soil_concentration(sand_idx, patch_idx, 1) * top_soil_active_depth_m
      available_silt_kg_m2 = tracer_soil_concentration(silt_idx, patch_idx, 1) * top_soil_active_depth_m
      available_clay_kg_m2 = tracer_soil_concentration(clay_idx, patch_idx, 1) * top_soil_active_depth_m
      
      total_available_soil_texture_kg_m2 = available_sand_kg_m2 + available_silt_kg_m2 + available_clay_kg_m2

      IF (total_available_soil_texture_kg_m2 > 1.0e-9_r8) THEN
          sand_frac = available_sand_kg_m2 / total_available_soil_texture_kg_m2
          silt_frac = available_silt_kg_m2 / total_available_soil_texture_kg_m2
          clay_frac = available_clay_kg_m2 / total_available_soil_texture_kg_m2

          eroded_sand_kg_m2 = total_soil_loss_kg_m2_event * sand_frac
          eroded_silt_kg_m2 = total_soil_loss_kg_m2_event * silt_frac
          eroded_clay_kg_m2 = total_soil_loss_kg_m2_event * clay_frac
          
          eroded_sand_kg_m2 = MIN(eroded_sand_kg_m2, available_sand_kg_m2)
          eroded_silt_kg_m2 = MIN(eroded_silt_kg_m2, available_silt_kg_m2)
          eroded_clay_kg_m2 = MIN(eroded_clay_kg_m2, available_clay_kg_m2)

          tracer_soil_concentration(sand_idx, patch_idx, 1) = tracer_soil_concentration(sand_idx, patch_idx, 1) - eroded_sand_kg_m2 / top_soil_active_depth_m
          tracer_soil_concentration(silt_idx, patch_idx, 1) = tracer_soil_concentration(silt_idx, patch_idx, 1) - eroded_silt_kg_m2 / top_soil_active_depth_m
          tracer_soil_concentration(clay_idx, patch_idx, 1) = tracer_soil_concentration(clay_idx, patch_idx, 1) - eroded_clay_kg_m2 / top_soil_active_depth_m
          
          tracer_soil_concentration(sand_idx, patch_idx, 1) = MAX(0.0_r8, tracer_soil_concentration(sand_idx, patch_idx, 1))
          tracer_soil_concentration(silt_idx, patch_idx, 1) = MAX(0.0_r8, tracer_soil_concentration(silt_idx, patch_idx, 1))
          tracer_soil_concentration(clay_idx, patch_idx, 1) = MAX(0.0_r8, tracer_soil_concentration(clay_idx, patch_idx, 1))

          eroded_sand_flux_kg_m2_s = eroded_sand_kg_m2 / dt_seconds
          eroded_silt_flux_kg_m2_s = eroded_silt_kg_m2 / dt_seconds
          eroded_clay_flux_kg_m2_s = eroded_clay_kg_m2 / dt_seconds
      ENDIF
  ENDIF

END SUBROUTINE Calculate_Soil_Erosion_Tracers

END MODULE MOD_Tracer_Soil_Erosion
