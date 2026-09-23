PROGRAM river_hist_spinup_harness
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_TimeManager, only: timestamp
   USE MOD_Hist, only: hist_out
   USE MOD_Grid_RiverLakeHistState, only: acctime_ucat, a_wdsrf_ucat, a_discharge
   USE MOD_Grid_RiverLakeTimeVars, only: wdsrf_ucat, acctime_rnof, acc_rnof_uc
   IMPLICIT NONE
   type(timestamp) :: stamp

   CALL spmd_init ()
   stamp = timestamp(2000, 1, 0)
   allocate (acctime_ucat(1), a_wdsrf_ucat(1), a_discharge(1))
   allocate (wdsrf_ucat(1), acc_rnof_uc(1))
   acctime_ucat = 10._r8
   a_wdsrf_ucat = 20._r8
   a_discharge = 30._r8
   wdsrf_ucat = 40._r8
   acc_rnof_uc = 50._r8
   acctime_rnof = 60._r8

   CALL hist_out ((/2000,1,0/), 1800._r8, stamp, stamp, stamp, '/tmp', 'spinup')

   IF (any(acctime_ucat /= 0._r8) .or. any(a_wdsrf_ucat /= 0._r8) .or. &
       any(a_discharge /= 0._r8)) STOP 1
   IF (any(wdsrf_ucat /= 40._r8) .or. any(acc_rnof_uc /= 50._r8) .or. &
       acctime_rnof /= 60._r8) STOP 2
   write(*,'(A)') 'SPINUP_ROUTE_HISTORY_OK'
   CALL spmd_exit ()
END PROGRAM river_hist_spinup_harness
