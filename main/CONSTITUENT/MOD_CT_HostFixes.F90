#include <define.h>

#ifdef ConstituentTransport
MODULE MOD_CT_HostFixes

!-----------------------------------------------------------------------
! Fixes for two upstream problems that break restart reproducibility.
! Constituent runs need a split run (restart) to reproduce an unbroken
! run exactly, so both fixes are applied here, and ONLY when at least one
! constituent is switched on (ct_any_enabled). With constituents off, or
! with ConstituentTransport not defined, the host keeps its upstream
! behaviour bit for bit: the host calls below are inside
! #ifdef ConstituentTransport and both functions return the upstream
! value when ct_any_enabled() is false.
!
! 1. Start year of the BGC nitrogen-deposition / fire data
!    (ct_bgc_data_init_year). CoLM.F90 passes sdate(1) after
!    adj2end(sdate). A run starting at 00:00 on 1 January has its start
!    date moved to the end of the previous year (day 365/366, 86400 s),
!    so the data of the previous year are read and kept for the whole
!    first year (the annual update in the time loop only fires when the
!    year changes). A run restarted later in the same year reads the
!    right year, so the split run and the unbroken run differ (ndep,
!    sminn_vr, ...). The year that the simulation actually starts in is
!    DEF_simulation_time%start_year (s_year), as CoLM202X uses.
!
! 2. History file continued after a restart (ct_hist_continue_file).
!    The host creates a new history file when the file name differs from
!    the last file written in this run; that name is 'null' at start-up,
!    so a run restarted inside a history group (e.g. on day 2 with
!    monthly files) creates the file again and loses the records written
!    before the restart. When the run starts inside the group of the file
!    it is about to write and the file exists, it is marked as the last
!    file, so the host appends to it (ncio_write_time finds the time
!    record or appends a new one).
!-----------------------------------------------------------------------

   USE MOD_CT_Namelist, only: ct_any_enabled
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: ct_bgc_data_init_year
   PUBLIC :: ct_hist_continue_file

CONTAINS

   integer FUNCTION ct_bgc_data_init_year (sdate_year, start_year)
   ! sdate_year: sdate(1) after adj2end (upstream argument)
   ! start_year: DEF_simulation_time%start_year
   IMPLICIT NONE
   integer, intent(in) :: sdate_year, start_year

      IF (ct_any_enabled()) THEN
         ct_bgc_data_init_year = start_year
      ELSE
         ct_bgc_data_init_year = sdate_year
      ENDIF

   END FUNCTION ct_bgc_data_init_year


   SUBROUTINE ct_hist_continue_file (file_hist, idate, file_last)
   ! file_hist: the history file about to be written
   ! idate:     current model time (end of the time step)
   ! file_last: the host's last-written file name; set to file_hist when
   !            the existing file has to be continued
   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_simulation_time, DEF_HIST_groupby, DEF_HIST_mode
   USE MOD_TimeManager, only: julian2monthday, monthday2julian, minutes_since_1900
   IMPLICIT NONE
   character(len=*), intent(in)    :: file_hist
   integer,          intent(in)    :: idate(3)
   character(len=*), intent(inout) :: file_last

   integer :: month, mday, jday_group, jday_start, min_group, min_start
   logical :: lcontinue

      IF (.not. ct_any_enabled()) RETURN
      ! only the first write of this run
      IF (trim(file_last) /= 'null') RETURN
      ! per-block files ('block' mode) carry other names; not handled
      IF (trim(DEF_HIST_mode) /= 'one') RETURN

      ! beginning of the history group that idate belongs to
      ! (idate is end-of-step, so 86400 s still belongs to that day)
      SELECT CASE (trim(DEF_HIST_groupby))
      CASE ('YEAR')
         jday_group = 1
      CASE ('MONTH')
         CALL julian2monthday (idate(1), idate(2), month, mday)
         CALL monthday2julian (idate(1), month, 1, jday_group)
      CASE ('DAY')
         jday_group = idate(2)
      CASE DEFAULT
         RETURN
      END SELECT
      min_group = minutes_since_1900 (idate(1), jday_group, 0)

      CALL monthday2julian (DEF_simulation_time%start_year, DEF_simulation_time%start_month, &
         DEF_simulation_time%start_day, jday_start)
      min_start = minutes_since_1900 (DEF_simulation_time%start_year, jday_start, &
         DEF_simulation_time%start_sec)

      ! A run starting at the beginning of the group writes a new file, as upstream.
      IF (min_start <= min_group) RETURN

      ! The file is opened by the master only; decide there so that all
      ! processes agree.
      lcontinue = .false.
      IF (p_is_master) inquire (file=trim(file_hist), exist=lcontinue)
#ifdef USEMPI
      CALL mpi_bcast (lcontinue, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif

      IF (lcontinue) THEN
         file_last = file_hist
         IF (p_is_master) write(*,'(A)') 'Constituent restart: continuing history file ' // trim(file_hist)
      ENDIF

   END SUBROUTINE ct_hist_continue_file

END MODULE MOD_CT_HostFixes
#endif
