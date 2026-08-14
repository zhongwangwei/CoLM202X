PROGRAM methane_host_water_harness

! ---------------------------------------------------------------------------
! Numerical regression for the methane saturated/unsaturated water partition.
!
! The tolerance change that widened the host-water guard was covered only by
! source-string assertions, and those cannot see a conservation error. They
! did not: inside the tolerated band the saturated column was over-allocated
! and the area-weighted total came out at finundated*pore_volume instead of
! the host's vtot, creating (finundated*pore_volume - vtot) of water per layer.
!
! This exercises the arithmetic itself. The invariant is the one the routine's
! own comment promises:
!
!    finundated*sat + (1-finundated)*unsat == host value, exactly,
!
! for liquid and ice separately, plus sat <= pore_volume and unsat >= 0.
!
! The partition is reproduced here rather than called, because it sits in the
! middle of methane() and has no separable entry point. Any change to those
! lines must be mirrored here -- which is the point: the mirror is what makes
! the invariant checkable at all.
! ---------------------------------------------------------------------------

   USE MOD_Precision

   IMPLICIT NONE

   integer  :: ncase, nfail
   real(r8) :: pore, f, vliq, vice

      nfail = 0
      ncase = 0

      ! pore_volume, finundated, vliq, vice
      CALL check ('deficit 1.6%% of pore, f=0.8', 0.04_r8, 0.8_r8, 0.03136_r8, 0.0_r8)
      CALL check ('deficit 2%% of pore,   f=0.5', 0.04_r8, 0.5_r8, 0.01920_r8, 0.0_r8)
      CALL check ('deficit with ice',             0.04_r8, 0.8_r8, 0.02000_r8, 0.01136_r8)
      CALL check ('deficit 4.9%% (just inside)',  0.04_r8, 0.8_r8, 0.03004_r8, 0.0_r8)
      CALL check ('exact vtot = f*pore',          0.04_r8, 0.8_r8, 0.03200_r8, 0.0_r8)
      CALL check ('normal',                       0.04_r8, 0.8_r8, 0.03600_r8, 0.0_r8)
      CALL check ('excess 1.6%% over pore',       0.04_r8, 0.8_r8, 0.04064_r8, 0.0_r8)
      CALL check ('fully inundated f=1',          0.04_r8, 1.0_r8, 0.03000_r8, 0.0_r8)
      CALL check ('dry column',                   0.04_r8, 0.8_r8, 0.0_r8,     0.0_r8)
      CALL check ('all ice',                      0.04_r8, 0.7_r8, 0.0_r8,     0.02500_r8)
      CALL check ('tiny finundated',              0.04_r8, 1.e-6_r8, 0.03_r8,  0.0_r8)
      CALL check ('thin layer',                   1.e-4_r8, 0.9_r8, 5.e-5_r8,  0.0_r8)

      write(*,'(A,I0,A,I0,A)') 'MHW ', ncase-nfail, '/', ncase, ' cases conserved'
      IF (nfail == 0) THEN
         write(*,'(A)') 'MHW PASS'
      ELSE
         write(*,'(A,I0)') 'MHW FAIL nfail=', nfail
         STOP 1
      ENDIF

CONTAINS

   SUBROUTINE check (name, pore_volume, finundated, vliq_in, vice_in)

   character(len=*), intent(in) :: name
   real(r8), intent(in) :: pore_volume, finundated, vliq_in, vice_in

   real(r8) :: vliq, vice, vtot, scale
   real(r8) :: vliq_sat_alloc, vice_sat_alloc, uliq, uice
   real(r8) :: wliq, wice
   real(r8), parameter :: tol = 1.e-14_r8
   logical :: ok

      ncase = ncase + 1
      vliq = vliq_in
      vice = vice_in
      vtot = vliq + vice

      ! ---- mirrors MOD_Tracer_Reactive_Methane_Physics.F90 ----
      scale = min(pore_volume / max(vtot, 1.e-12_r8), &
                  1._r8 / max(finundated, 1.e-12_r8))
      vliq_sat_alloc = vliq * scale
      vice_sat_alloc = vice * scale
      IF (finundated < 1._r8) THEN
         uliq = max(0._r8, (vliq - finundated*vliq_sat_alloc) / (1._r8 - finundated))
         uice = max(0._r8, (vice - finundated*vice_sat_alloc) / (1._r8 - finundated))
      ELSE
         uliq = 0._r8
         uice = 0._r8
      ENDIF
      ! ---------------------------------------------------------

      wliq = finundated*vliq_sat_alloc + (1._r8 - finundated)*uliq
      wice = finundated*vice_sat_alloc + (1._r8 - finundated)*uice

      ok = .true.
      IF (abs(wliq - vliq) > tol * max(vliq, 1._r8)) ok = .false.
      IF (abs(wice - vice) > tol * max(vice, 1._r8)) ok = .false.
      IF (vliq_sat_alloc + vice_sat_alloc > pore_volume * (1._r8 + 1.e-12_r8)) ok = .false.
      IF (uliq < -tol .or. uice < -tol) ok = .false.

      IF (.not. ok) THEN
         nfail = nfail + 1
         write(*,'(A,A)')        'MHW FAIL case: ', trim(name)
         write(*,'(A,3E22.14)')  '    vliq host/weighted/created : ', vliq, wliq, wliq-vliq
         write(*,'(A,3E22.14)')  '    vice host/weighted/created : ', vice, wice, wice-vice
         write(*,'(A,2E22.14)')  '    sat total / pore_volume    : ', &
            vliq_sat_alloc+vice_sat_alloc, pore_volume
      ENDIF

   END SUBROUTINE check

END PROGRAM methane_host_water_harness
