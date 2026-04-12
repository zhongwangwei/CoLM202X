#include <define.h>

MODULE MOD_Tracer_Conservation

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny
   USE MOD_Tracer_Vars

   IMPLICIT NONE

   real(r8), parameter :: trc_balance_tol = 1.0e-10_r8

CONTAINS

   SUBROUTINE tracer_save_storage (ipatch, snl, nl_soil)
      IMPLICIT NONE
      integer, intent(in) :: ipatch, snl, nl_soil
      integer :: itrc, j

      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         trc_storage_beg(itrc, ipatch) = 0._r8
         trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
            + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
         DO j = snl + 1, nl_soil
            trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
         trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
            + trc_wa(itrc, ipatch) + trc_wdsrf(itrc, ipatch) + trc_wetwat(itrc, ipatch)
      ENDDO
   END SUBROUTINE tracer_save_storage

   SUBROUTINE tracer_balance_check (ipatch, snl, nl_soil, deltim, xerr_tracer)
      IMPLICIT NONE
      integer,  intent(in)  :: ipatch, snl, nl_soil
      real(r8), intent(in)  :: deltim
      real(r8), intent(out) :: xerr_tracer

      integer  :: itrc, j
      real(r8) :: storage_end, total_input, total_output, err

      xerr_tracer = 0._r8
      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         storage_end = 0._r8
         storage_end = storage_end + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
         DO j = snl + 1, nl_soil
            storage_end = storage_end &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
         storage_end = storage_end &
            + trc_wa(itrc, ipatch) + trc_wdsrf(itrc, ipatch) + trc_wetwat(itrc, ipatch)

         total_input = a_trc_precip(itrc, ipatch)
         total_output = a_trc_evap(itrc, ipatch) + a_trc_rnof(itrc, ipatch)

         err = storage_end - trc_storage_beg(itrc, ipatch) - total_input + total_output
         trc_balance_err(itrc, ipatch) = err

         IF (abs(err) > trc_balance_tol) THEN
            WRITE(*,'(A,A,A,I8,A,E15.7)') &
               'WARNING: Tracer balance error [', TRIM(tracers(itrc)%name), &
               '] patch=', ipatch, ' err=', err
         ENDIF

         xerr_tracer = max(xerr_tracer, abs(err) / deltim)
      ENDDO
   END SUBROUTINE tracer_balance_check

END MODULE MOD_Tracer_Conservation
