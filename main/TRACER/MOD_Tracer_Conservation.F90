#include <define.h>

MODULE MOD_Tracer_Conservation

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny
   USE MOD_Tracer_Vars

   IMPLICIT NONE

   real(r8), parameter :: trc_balance_tol = 1.0e-10_r8

   ! Per-step accumulator snapshots (saved at start of each timestep)
   real(r8), allocatable, save :: snap_precip(:,:)  ! (ntracers, numpatch)
   real(r8), allocatable, save :: snap_evap  (:,:)
   real(r8), allocatable, save :: snap_rnof  (:,:)

CONTAINS

   SUBROUTINE tracer_save_storage (ipatch, snl, nl_soil)
      IMPLICIT NONE
      integer, intent(in) :: ipatch, snl, nl_soil
      integer :: itrc, j

      IF (ntracers <= 0) RETURN

      ! Allocate snapshots on first call
      IF (.not. allocated(snap_precip)) THEN
         allocate(snap_precip(ntracers, size(trc_storage_beg,2))); snap_precip = 0._r8
         allocate(snap_evap  (ntracers, size(trc_storage_beg,2))); snap_evap   = 0._r8
         allocate(snap_rnof  (ntracers, size(trc_storage_beg,2))); snap_rnof   = 0._r8
      ENDIF

      DO itrc = 1, ntracers
         ! Save current storage
         trc_storage_beg(itrc, ipatch) = 0._r8
         trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
            + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
         DO j = snl + 1, nl_soil
            trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
         trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
            + trc_wa(itrc, ipatch) + trc_wdsrf(itrc, ipatch) + trc_wetwat(itrc, ipatch)

         ! Snapshot accumulators at start of this step
         snap_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch)
         snap_evap  (itrc, ipatch) = a_trc_evap  (itrc, ipatch)
         snap_rnof  (itrc, ipatch) = a_trc_rnof  (itrc, ipatch)
      ENDDO
   END SUBROUTINE tracer_save_storage

   SUBROUTINE tracer_balance_check (ipatch, snl, nl_soil, deltim, xerr_tracer)
      IMPLICIT NONE
      integer,  intent(in)  :: ipatch, snl, nl_soil
      real(r8), intent(in)  :: deltim
      real(r8), intent(out) :: xerr_tracer

      integer  :: itrc, j
      real(r8) :: storage_end, step_input, step_output, err

      xerr_tracer = 0._r8
      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         ! Current total storage
         storage_end = 0._r8
         storage_end = storage_end + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
         DO j = snl + 1, nl_soil
            storage_end = storage_end &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
         storage_end = storage_end &
            + trc_wa(itrc, ipatch) + trc_wdsrf(itrc, ipatch) + trc_wetwat(itrc, ipatch)

         ! Per-step fluxes = current accumulator - snapshot at step start
         step_input  = a_trc_precip(itrc, ipatch) - snap_precip(itrc, ipatch)
         step_output = (a_trc_evap(itrc, ipatch) - snap_evap(itrc, ipatch)) &
                     + (a_trc_rnof(itrc, ipatch) - snap_rnof(itrc, ipatch))

         ! Conservation: Δstorage = input - output
         err = storage_end - trc_storage_beg(itrc, ipatch) - step_input + step_output
         trc_balance_err(itrc, ipatch) = err

         IF (abs(err) > trc_balance_tol .and. ipatch <= 3) THEN
            WRITE(*,'(A,A,A,I6,A,E12.5,A,E12.5,A,E12.5,A,E12.5,A,E12.5)') &
               'TRC_ERR [', TRIM(tracers(itrc)%name), &
               '] p=', ipatch, &
               ' dS=', storage_end - trc_storage_beg(itrc, ipatch), &
               ' in=', step_input, &
               ' out=', step_output, &
               ' err=', err, &
               ' S_end=', storage_end
         ENDIF

         xerr_tracer = max(xerr_tracer, abs(err) / deltim)
      ENDDO
   END SUBROUTINE tracer_balance_check

END MODULE MOD_Tracer_Conservation
