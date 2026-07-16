#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_EvapLimit
!=======================================================================
! Shared finite-pool evaporation/sublimation limiter for tracer losses.
! Callers provide the current isotope/reactive evaporation ratio callback;
! this helper owns the stiff-substep constants so canopy, soil and wetland
! paths cannot silently diverge.
!=======================================================================

   USE MOD_Precision

   IMPLICIT NONE
   PRIVATE

   real(r8), parameter :: evaplimit_default_max_loss_fraction = 0.10_r8
   integer,  parameter :: evaplimit_default_max_substeps = 80

   abstract interface
      real(r8) FUNCTION tracer_evap_ratio_callback (source_ratio, temp_k, from_ice)
         IMPORT :: r8
         real(r8), intent(in) :: source_ratio
         real(r8), intent(in) :: temp_k
         logical,  intent(in) :: from_ice
      END FUNCTION tracer_evap_ratio_callback
   end interface

   PUBLIC :: tracer_evaporative_tracer_loss
   PUBLIC :: tracer_atmospheric_tracer_loss

CONTAINS

   real(r8) FUNCTION tracer_atmospheric_tracer_loss (pool_trc, pool_water, water_loss, &
      temp_k, from_ice, evap_ratio_fn, trc_tiny, r_max, is_nonvolatile)

      real(r8), intent(in) :: pool_trc
      real(r8), intent(in) :: pool_water
      real(r8), intent(in) :: water_loss
      real(r8), intent(in) :: temp_k
      logical,  intent(in) :: from_ice
      procedure(tracer_evap_ratio_callback) :: evap_ratio_fn
      real(r8), intent(in) :: trc_tiny
      real(r8), intent(in) :: r_max
      logical,  intent(in) :: is_nonvolatile

      IF (is_nonvolatile) THEN
         tracer_atmospheric_tracer_loss = 0._r8
      ELSE
         tracer_atmospheric_tracer_loss = tracer_evaporative_tracer_loss(pool_trc, pool_water, &
            water_loss, temp_k, from_ice, evap_ratio_fn, trc_tiny, r_max)
      ENDIF
   END FUNCTION tracer_atmospheric_tracer_loss

   real(r8) FUNCTION tracer_evaporative_tracer_loss (pool_trc, pool_water, water_loss, &
      temp_k, from_ice, evap_ratio_fn, trc_tiny, r_max)

      real(r8), intent(in) :: pool_trc
      real(r8), intent(in) :: pool_water
      real(r8), intent(in) :: water_loss
      real(r8), intent(in) :: temp_k
      logical,  intent(in) :: from_ice
      procedure(tracer_evap_ratio_callback) :: evap_ratio_fn
      real(r8), intent(in) :: trc_tiny
      ! NUMERICAL/PHYSICAL safety ceiling on the residual pool concentration R
      ! -- this is NOT part of the Craig-Gordon fractionation itself. Evaporative
      ! fractionation enriches the residual without bound as a pool nears
      ! drydown; over repeated wet/dry cycles the in-pool tracer can run away to
      ! non-physical values (delta up to ~+1e13 permil). Once source_ratio
      ! reaches r_max, evaporation switches to NON-FRACTIONATING
      ! (flux_ratio = source_ratio): the vapour carries water_loss*R, which is
      ! water-matched and exactly conservative, so the residual R freezes near
      ! r_max instead of running away -- without shedding any water-less excess.
      ! Pass r_max<=0 to DISABLE; nonvolatile solutes should be handled by
      ! callers as zero evaporative tracer loss.
      real(r8), intent(in) :: r_max

      integer  :: isub
      real(r8) :: source_ratio, flux_ratio
      real(r8) :: remaining_trc, remaining_water, loss_left, step_loss, trc_loss_step

      tracer_evaporative_tracer_loss = 0._r8
      IF (pool_water <= trc_tiny .or. water_loss <= trc_tiny) RETURN
      IF (pool_trc <= trc_tiny) RETURN

      IF (water_loss >= pool_water * (1._r8 - 1.e-12_r8)) THEN
         ! ALL_POOL_EVAPORATION: when the water-side flux removes the whole
         ! finite pool, the integrated outgoing tracer is the whole tracer
         ! inventory.  Applying a last instantaneous fractionation ratio would
         ! leave non-zero tracer in a zero-water pool and break conservation.
         tracer_evaporative_tracer_loss = max(pool_trc, 0._r8)
         RETURN
      ENDIF

      source_ratio = max(pool_trc, 0._r8) / pool_water
      flux_ratio = evap_ratio_fn(source_ratio, temp_k, from_ice)
      ! Stop evaporative enrichment at the physical ceiling. Above r_max the
      ! vapour leaves with the pool's own ratio, so tracer loss remains matched
      ! to water loss and the residual pool R stops increasing.
      IF (r_max > 0._r8 .and. source_ratio > r_max) flux_ratio = source_ratio
      IF (water_loss <= evaplimit_default_max_loss_fraction * pool_water .or. &
          abs(flux_ratio - source_ratio) <= 1.e-12_r8 * max(source_ratio, trc_tiny)) THEN
         tracer_evaporative_tracer_loss = min(water_loss * flux_ratio, max(pool_trc, 0._r8))
         RETURN
      ENDIF

      remaining_trc = max(pool_trc, 0._r8)
      remaining_water = pool_water
      loss_left = water_loss
      DO isub = 1, evaplimit_default_max_substeps
         IF (loss_left <= trc_tiny) EXIT
         IF (remaining_water <= trc_tiny .or. remaining_trc <= trc_tiny) EXIT

         step_loss = min(loss_left, evaplimit_default_max_loss_fraction * remaining_water)
         IF (isub == evaplimit_default_max_substeps) step_loss = loss_left
         step_loss = min(step_loss, remaining_water)

         source_ratio = remaining_trc / remaining_water
         flux_ratio = evap_ratio_fn(source_ratio, temp_k, from_ice)
         ! Stop enrichment at the ceiling: non-fractionating evap above r_max
         ! freezes the residual R near r_max while preserving water-matched loss.
         IF (r_max > 0._r8 .and. source_ratio > r_max) flux_ratio = source_ratio
         trc_loss_step = min(step_loss * flux_ratio, remaining_trc)
         remaining_trc = remaining_trc - trc_loss_step
         remaining_water = remaining_water - step_loss
         loss_left = loss_left - step_loss
      ENDDO

      tracer_evaporative_tracer_loss = min(max(pool_trc, 0._r8), &
         max(pool_trc, 0._r8) - max(remaining_trc, 0._r8))

   END FUNCTION tracer_evaporative_tracer_loss

END MODULE MOD_Tracer_EvapLimit
#endif
