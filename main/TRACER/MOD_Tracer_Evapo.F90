#include <define.h>

MODULE MOD_Tracer_Evapo

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars, only: trc_ldew_rain, trc_ldew_snow, &
      trc_wliq_soisno, trc_wice_soisno, &
      a_trc_evap, a_trc_trans, a_trc_precip

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_evapo (ipatch, deltim, snl, nl_soil, &
      fevpl, etr, qseva, qsdew, qsubl, qfros, &
      qseva_soil, qsdew_soil, qsubl_soil, qfros_soil, &
      qseva_snow, qsdew_snow, qsubl_snow, qfros_snow, &
      rootflux, wliq_soisno, wice_soisno, &
      ldew_rain, ldew_snow)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: snl, nl_soil
      real(r8), intent(in) :: fevpl, etr
      real(r8), intent(in) :: qseva, qsdew, qsubl, qfros
      real(r8), intent(in) :: qseva_soil, qsdew_soil, qsubl_soil, qfros_soil
      real(r8), intent(in) :: qseva_snow, qsdew_snow, qsubl_snow, qfros_snow
      real(r8), intent(in) :: rootflux(1:nl_soil)
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: ldew_rain, ldew_snow

      integer  :: itrc, j, lb
      real(r8) :: ratio, trc_flux, R_atm
      real(r8) :: canopy_evap_liq

      IF (ntracers <= 0) RETURN
      lb = snl + 1
      canopy_evap_liq = max(fevpl - etr, 0._r8)

      DO itrc = 1, ntracers
         R_atm = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! Canopy liquid evaporation
         IF (canopy_evap_liq > trc_tiny .and. ldew_rain > trc_tiny) THEN
            ratio = trc_ldew_rain(itrc, ipatch) / ldew_rain
            trc_flux = canopy_evap_liq * ratio * deltim
            trc_flux = min(trc_flux, trc_ldew_rain(itrc, ipatch))
            trc_ldew_rain(itrc, ipatch) = trc_ldew_rain(itrc, ipatch) - trc_flux
            a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
         ENDIF

         ! Soil evaporation
         IF (qseva_soil > trc_tiny) THEN
            j = max(lb, 1)
            IF (wliq_soisno(j) > trc_tiny) THEN
               ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno(j)
               trc_flux = qseva_soil * ratio * deltim
               trc_flux = min(trc_flux, trc_wliq_soisno(itrc, j, ipatch))
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

         ! Soil dew
         IF (qsdew_soil > trc_tiny) THEN
            j = max(lb, 1)
            trc_flux = qsdew_soil * R_atm * deltim
            trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF

         ! Snow evaporation
         IF (qseva_snow > trc_tiny .and. snl < 0) THEN
            j = lb
            IF (wliq_soisno(j) > trc_tiny) THEN
               ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno(j)
               trc_flux = qseva_snow * ratio * deltim
               trc_flux = min(trc_flux, trc_wliq_soisno(itrc, j, ipatch))
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

         ! Soil sublimation
         IF (qsubl_soil > trc_tiny) THEN
            j = max(lb, 1)
            IF (wice_soisno(j) > trc_tiny) THEN
               ratio = trc_wice_soisno(itrc, j, ipatch) / wice_soisno(j)
               trc_flux = qsubl_soil * ratio * deltim
               trc_flux = min(trc_flux, trc_wice_soisno(itrc, j, ipatch))
               trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

         ! Snow sublimation
         IF (qsubl_snow > trc_tiny .and. snl < 0) THEN
            j = lb
            IF (wice_soisno(j) > trc_tiny) THEN
               ratio = trc_wice_soisno(itrc, j, ipatch) / wice_soisno(j)
               trc_flux = qsubl_snow * ratio * deltim
               trc_flux = min(trc_flux, trc_wice_soisno(itrc, j, ipatch))
               trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ENDIF
         ENDIF

         ! Frost on soil
         IF (qfros_soil > trc_tiny) THEN
            j = max(lb, 1)
            trc_flux = qfros_soil * R_atm * deltim
            trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF

         ! Frost on snow
         IF (qfros_snow > trc_tiny .and. snl < 0) THEN
            j = lb
            trc_flux = qfros_snow * R_atm * deltim
            trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF

         ! Snow dew
         IF (qsdew_snow > trc_tiny .and. snl < 0) THEN
            j = lb
            trc_flux = qsdew_snow * R_atm * deltim
            trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
            a_trc_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch) + trc_flux
         ENDIF

         ! Transpiration
         DO j = 1, nl_soil
            IF (rootflux(j) > trc_tiny .and. wliq_soisno(j) > trc_tiny) THEN
               ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno(j)
               trc_flux = rootflux(j) * ratio * deltim
               trc_flux = min(trc_flux, trc_wliq_soisno(itrc, j, ipatch))
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
               a_trc_trans(itrc, ipatch) = a_trc_trans(itrc, ipatch) + trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) + trc_flux
            ELSEIF (rootflux(j) < -trc_tiny .and. wliq_soisno(j) > trc_tiny) THEN
               ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno(j)
               trc_flux = abs(rootflux(j)) * ratio * deltim
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
               a_trc_trans(itrc, ipatch) = a_trc_trans(itrc, ipatch) - trc_flux
               a_trc_evap(itrc, ipatch) = a_trc_evap(itrc, ipatch) - trc_flux
            ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE tracer_evapo

END MODULE MOD_Tracer_Evapo
