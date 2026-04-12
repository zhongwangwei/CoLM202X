#include <define.h>

MODULE MOD_Tracer_SoilWater

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars, only: trc_wliq_soisno, trc_wice_soisno, &
      trc_wa, trc_wdsrf, trc_wetwat, &
      a_trc_qinfl, a_trc_qcharge

   IMPLICIT NONE

CONTAINS

   SUBROUTINE tracer_soil_water (ipatch, deltim, snl, nl_soil, &
      qinfl, qlayer, qcharge, &
      wliq_soisno, wice_soisno, &
      wliq_soisno_old, wice_soisno_old, &
      wa, wa_old, wdsrf, wdsrf_old, &
      wetwat, wetwat_old, pg_rain, pg_snow, &
      use_vsf)

      IMPLICIT NONE
      integer,  intent(in) :: ipatch
      real(r8), intent(in) :: deltim
      integer,  intent(in) :: snl, nl_soil
      real(r8), intent(in) :: qinfl
      real(r8), intent(in) :: qlayer(0:nl_soil)
      real(r8), intent(in) :: qcharge
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wliq_soisno_old(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno_old(snl+1:nl_soil)
      real(r8), intent(in) :: wa, wa_old, wdsrf, wdsrf_old
      real(r8), intent(in) :: wetwat, wetwat_old
      real(r8), intent(in) :: pg_rain, pg_snow
      logical,  intent(in) :: use_vsf

      integer  :: itrc, j, lb, j_src
      real(r8) :: ratio, trc_flux, R_precip
      real(r8) :: dwliq, dwice

      IF (ntracers <= 0) RETURN
      lb = snl + 1

      DO itrc = 1, ntracers
         R_precip = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)

         ! Throughfall to surface water
         trc_flux = (pg_rain + pg_snow) * R_precip * deltim
         trc_wdsrf(itrc, ipatch) = trc_wdsrf(itrc, ipatch) + trc_flux

         ! Infiltration
         IF (qinfl > trc_tiny) THEN
            IF (wdsrf > trc_tiny) THEN
               ratio = trc_wdsrf(itrc, ipatch) / wdsrf
            ELSE
               ratio = R_precip
            ENDIF
            trc_flux = qinfl * ratio * deltim
            trc_flux = min(trc_flux, trc_wdsrf(itrc, ipatch))
            trc_wdsrf(itrc, ipatch) = trc_wdsrf(itrc, ipatch) - trc_flux
            j = max(lb, 1)
            trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
            a_trc_qinfl(itrc, ipatch) = a_trc_qinfl(itrc, ipatch) + trc_flux
         ENDIF

         ! Inter-layer flux
         DO j = 1, nl_soil - 1
            IF (qlayer(j) > trc_tiny) THEN
               j_src = j
               IF (wliq_soisno_old(j_src) > trc_tiny) THEN
                  ratio = trc_wliq_soisno(itrc, j_src, ipatch) / max(wliq_soisno(j_src), trc_tiny)
                  trc_flux = qlayer(j) * ratio * deltim
                  trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j_src, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j,   ipatch) = trc_wliq_soisno(itrc, j,   ipatch) - trc_flux
                  trc_wliq_soisno(itrc, j+1, ipatch) = trc_wliq_soisno(itrc, j+1, ipatch) + trc_flux
               ENDIF
            ELSEIF (qlayer(j) < -trc_tiny) THEN
               j_src = j + 1
               IF (wliq_soisno_old(j_src) > trc_tiny) THEN
                  ratio = trc_wliq_soisno(itrc, j_src, ipatch) / max(wliq_soisno(j_src), trc_tiny)
                  trc_flux = abs(qlayer(j)) * ratio * deltim
                  trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j_src, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j+1, ipatch) = trc_wliq_soisno(itrc, j+1, ipatch) - trc_flux
                  trc_wliq_soisno(itrc, j,   ipatch) = trc_wliq_soisno(itrc, j,   ipatch) + trc_flux
               ENDIF
            ENDIF
         ENDDO

         ! Groundwater recharge
         IF (qcharge > trc_tiny) THEN
            j = nl_soil
            IF (wliq_soisno(j) > trc_tiny) THEN
               ratio = trc_wliq_soisno(itrc, j, ipatch) / wliq_soisno(j)
               trc_flux = qcharge * ratio * deltim
               trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
               trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) + trc_flux
               a_trc_qcharge(itrc, ipatch) = a_trc_qcharge(itrc, ipatch) + trc_flux
            ENDIF
         ELSEIF (qcharge < -trc_tiny) THEN
            j = nl_soil
            IF (wa > trc_tiny) THEN
               ratio = trc_wa(itrc, ipatch) / max(wa, trc_tiny)
               trc_flux = abs(qcharge) * ratio * deltim
               trc_flux = min(trc_flux, max(trc_wa(itrc, ipatch), 0._r8))
               trc_wa(itrc, ipatch) = trc_wa(itrc, ipatch) - trc_flux
               trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
               a_trc_qcharge(itrc, ipatch) = a_trc_qcharge(itrc, ipatch) - trc_flux
            ENDIF
         ENDIF

         ! Freeze/thaw
         DO j = lb, nl_soil
            dwliq = wliq_soisno(j) - wliq_soisno_old(j)
            dwice = wice_soisno(j) - wice_soisno_old(j)
            IF (dwliq < -trc_tiny .and. dwice > trc_tiny) THEN
               IF (wliq_soisno_old(j) > trc_tiny) THEN
                  ratio = trc_wliq_soisno(itrc, j, ipatch) / max(wliq_soisno(j) - dwliq, trc_tiny)
                  trc_flux = abs(dwliq) * ratio
                  trc_flux = min(trc_flux, max(trc_wliq_soisno(itrc, j, ipatch), 0._r8))
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) - trc_flux
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ELSEIF (dwice < -trc_tiny .and. dwliq > trc_tiny) THEN
               IF (wice_soisno_old(j) > trc_tiny) THEN
                  ratio = trc_wice_soisno(itrc, j, ipatch) / max(wice_soisno(j) - dwice, trc_tiny)
                  trc_flux = abs(dwice) * ratio
                  trc_flux = min(trc_flux, max(trc_wice_soisno(itrc, j, ipatch), 0._r8))
                  trc_wice_soisno(itrc, j, ipatch) = trc_wice_soisno(itrc, j, ipatch) - trc_flux
                  trc_wliq_soisno(itrc, j, ipatch) = trc_wliq_soisno(itrc, j, ipatch) + trc_flux
               ENDIF
            ENDIF
         ENDDO

         ! Wetland (VSF mode)
         IF (use_vsf) THEN
            IF (wetwat > trc_tiny .and. wetwat_old > trc_tiny) THEN
               IF (wetwat > wetwat_old) THEN
                  trc_flux = (wetwat - wetwat_old) * R_precip
                  trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) + trc_flux
               ELSE
                  ratio = trc_wetwat(itrc, ipatch) / wetwat_old
                  trc_flux = (wetwat_old - wetwat) * ratio
                  trc_flux = min(trc_flux, max(trc_wetwat(itrc, ipatch), 0._r8))
                  trc_wetwat(itrc, ipatch) = trc_wetwat(itrc, ipatch) - trc_flux
               ENDIF
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE tracer_soil_water

END MODULE MOD_Tracer_SoilWater
