#include <define.h>

MODULE MOD_Tracer_Overland_Flow_Seds

USE MOD_Precision, only: r8
USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, DEF_Tracer_Number, tracer_info_type
USE MOD_Tracer_Sediment_Util, only: DIAM_SAND_COLM, DIAM_SILT_COLM, DIAM_CLAY_COLM 
! For particle properties. Assumes these are in meters.
USE MOD_SPMD_Task, only: p_is_master ! For conditional writes

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC :: Initialize_Overland_Sed_Params
PUBLIC :: Transport_Sediment_Overland_Flow

! Parameters for transport equations (e.g., coefficients for a chosen TC formula)
! real(r8) :: tc_coeff_a, tc_coeff_b, etc.

CONTAINS

SUBROUTINE Initialize_Overland_Sed_Params()
  IF (p_is_master) THEN
    WRITE(*,*) "MOD_Tracer_Overland_Flow_Seds: Initializing parameters (conceptual)."
  ENDIF
  ! This subroutine would read any specific parameters needed for the
  ! chosen overland flow sediment transport model from a namelist or set defaults.
END SUBROUTINE Initialize_Overland_Sed_Params

SUBROUTINE Transport_Sediment_Overland_Flow( &
    patch_idx, &
    eroded_sand_flux_kg_m2_s,  & ! Input from MOD_Tracer_Soil_Erosion
    eroded_silt_flux_kg_m2_s,  & ! Input from MOD_Tracer_Soil_Erosion
    eroded_clay_flux_kg_m2_s,  & ! Input from MOD_Tracer_Soil_Erosion
    overland_flow_depth_m,     & ! From CoLM Hydrology
    overland_flow_velocity_ms, & ! From CoLM Hydrology
    patch_area_m2,             & ! From CoLM Patch data
    flow_path_length_m,        & ! From CoLM Patch data (or derived)
    slope_decimal,             & ! From CoLM Patch data (e.g., m/m)
    dt_seconds,                & ! Model timestep
    yield_sand_kgs,            & ! Output: kg/s delivered to channel
    yield_silt_kgs,            & ! Output: kg/s delivered to channel
    yield_clay_kgs )             ! Output: kg/s delivered to channel
    
  integer, intent(in) :: patch_idx
  real(r8), intent(in) :: eroded_sand_flux_kg_m2_s 
  real(r8), intent(in) :: eroded_silt_flux_kg_m2_s 
  real(r8), intent(in) :: eroded_clay_flux_kg_m2_s 
  real(r8), intent(in) :: overland_flow_depth_m     
  real(r8), intent(in) :: overland_flow_velocity_ms 
  real(r8), intent(in) :: patch_area_m2             
  real(r8), intent(in) :: flow_path_length_m        
  real(r8), intent(in) :: slope_decimal             
  real(r8), intent(in) :: dt_seconds                
  real(r8), intent(out) :: yield_sand_kgs            
  real(r8), intent(out) :: yield_silt_kgs            
  real(r8), intent(out) :: yield_clay_kgs            

  ! --- Local Variables ---
  real(r8) :: input_sand_kg, input_silt_kg, input_clay_kg ! Total mass eroded onto surface for this dt
  real(r8) :: tc_sand_kg, tc_silt_kg, tc_clay_kg          ! Transport capacity for the patch over dt (kg)
  ! real(r8) :: deposited_sand_kg, deposited_silt_kg, deposited_clay_kg ! (Optional to track deposition amount)

  real(r8), parameter :: WATER_DENSITY_KG_M3 = 1000.0_r8
  real(r8), parameter :: GRAVITY_MS2 = 9.81_r8
  
  ! Placeholder efficiency factors for conceptual TC calculation.
  ! These would be replaced by physically-based model coefficients.
  real(r8) :: efficiency_factor_sand = 0.002_r8 
  real(r8) :: efficiency_factor_silt = 0.01_r8  
  real(r8) :: efficiency_factor_clay = 0.02_r8  ! Clay might stay in suspension longer

  yield_sand_kgs = 0.0_r8
  yield_silt_kgs = 0.0_r8
  yield_clay_kgs = 0.0_r8

  IF (overland_flow_velocity_ms < 1.0e-6_r8 .OR. overland_flow_depth_m < 1.0e-6_r8) THEN
      ! No overland flow, so no transport from this patch (eroded material stays or handled by infiltration)
      RETURN
  ENDIF
  IF (patch_area_m2 < 1.0e-6_r8) RETURN
  IF (dt_seconds < 1.0e-6_r8) RETURN


  ! 1. Calculate total eroded mass entering overland flow during this timestep (kg)
  input_sand_kg = eroded_sand_flux_kg_m2_s * patch_area_m2 * dt_seconds
  input_silt_kg = eroded_silt_flux_kg_m2_s * patch_area_m2 * dt_seconds
  input_clay_kg = eroded_clay_flux_kg_m2_s * patch_area_m2 * dt_seconds

  ! 2. Calculate Transport Capacity (TC) for each sediment type (kg for the patch over dt)
  !    This is a highly conceptual placeholder using stream power concept.
  !    Stream Power per unit area (omega) = rho_w * g * depth * velocity * slope  (W/m^2 or J/s/m^2)
  !    TC_mass = efficiency * (omega / (particle_specific_weight * settling_velocity_or_diameter_term)) * area * dt
  !    For simplicity, let's use a very basic form: TC ~ efficiency * stream_power * area * dt / (particle_diameter_metric)
  
  real(r8) :: stream_power_total_W ! Total stream power for the patch (Watts or J/s)
  ! (overland_flow_depth_m * overland_flow_velocity_ms * patch_area_m2) is approx Q (m3/s) * (m_width / m_width) if specific discharge is per unit width
  ! If overland_flow_specific_discharge_ms is q (m^2/s = depth*velocity), then Q = q * width.
  ! Assuming overland_flow_velocity_ms is average velocity and overland_flow_depth_m is average depth.
  ! Q_patch = overland_flow_depth_m * overland_flow_velocity_ms * flow_width_m (conceptual width)
  ! For a patch, total stream power = rho_w * g * Q_patch * slope_decimal
  ! Q_patch * slope = (depth * velocity * width) * slope
  ! Stream power per unit area = rho_w * g * depth * velocity * slope
  ! Total stream power = (Stream power per unit area) * patch_area_m2
  
  stream_power_total_W = WATER_DENSITY_KG_M3 * GRAVITY_MS2 * overland_flow_depth_m * overland_flow_velocity_ms * slope_decimal * patch_area_m2

  ! Conceptual transport capacity (kg that can be transported over dt)
  ! The denominator involving diameter is to make TC smaller for larger particles.
  ! This is NOT a physically validated formula, just a placeholder structure.
  IF (DIAM_SAND_COLM > 1.0e-9_r8) THEN
    tc_sand_kg = efficiency_factor_sand * stream_power_total_W * dt_seconds / (DIAM_SAND_COLM * GRAVITY_MS2)
  ELSE
    tc_sand_kg = 1.0e9_r8 ! Effectively infinite TC if diameter is zero/tiny (for non-settling part)
  ENDIF
  
  IF (DIAM_SILT_COLM > 1.0e-9_r8) THEN
    tc_silt_kg = efficiency_factor_silt * stream_power_total_W * dt_seconds / (DIAM_SILT_COLM * GRAVITY_MS2)
  ELSE
    tc_silt_kg = 1.0e9_r8
  ENDIF

  IF (DIAM_CLAY_COLM > 1.0e-9_r8) THEN
    tc_clay_kg = efficiency_factor_clay * stream_power_total_W * dt_seconds / (DIAM_CLAY_COLM * GRAVITY_MS2)
  ELSE
    tc_clay_kg = 1.0e9_r8
  ENDIF
  
  tc_sand_kg = MAX(0.0_r8, tc_sand_kg)
  tc_silt_kg = MAX(0.0_r8, tc_silt_kg)
  tc_clay_kg = MAX(0.0_r8, tc_clay_kg)

  ! 3. Calculate Sediment Yield to Channel (kg/s)
  !    Simplified: Yield is the minimum of input material and transport capacity.
  !    (Assumes no significant pre-existing load in overland flow from upslope patches for this simple model)
  
  yield_sand_kgs = MIN(input_sand_kg, tc_sand_kg) / dt_seconds
  yield_silt_kgs = MIN(input_silt_kg, tc_silt_kg) / dt_seconds
  yield_clay_kgs = MIN(input_clay_kg, tc_clay_kg) / dt_seconds
  
  ! Optional: Calculate deposition within the patch
  ! deposited_sand_kg = input_sand_kg - (yield_sand_kgs * dt_seconds)
  ! deposited_silt_kg = input_silt_kg - (yield_silt_kgs * dt_seconds)
  ! deposited_clay_kg = input_clay_kg - (yield_clay_kgs * dt_seconds)
  ! This deposited material could be added to a conceptual 'surface loose sediment' store
  ! or simply considered "lost" from the cascade for this timestep if not transported.

END SUBROUTINE Transport_Sediment_Overland_Flow

END MODULE MOD_Tracer_Overland_Flow_Seds
