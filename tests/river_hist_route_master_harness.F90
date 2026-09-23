PROGRAM river_hist_route_master_harness

! ---------------------------------------------------------------------------
! Drives the REAL route-history dispatcher (route_hist_begin, route_hist_write_*,
! route_hist_end/final in MOD_Grid_RiverLakeHistRoute) in DEF_HIST_mode='block'
! under a real CoLM SPMD layout that includes the MASTER.
!
! Unlike river_hist_shard_harness.F90, which calls the shard writers directly
! and wraps them in its own rank guard, this program calls the dispatcher on
! EVERY rank exactly as MOD_Grid_RiverLakeHist / the tracer writers do.  So it
! fails if the master reaches a group collective, and it writes two history
! files in one process so that per-file state (the bifurcation dimensions) is
! exercised as well.
!
! Rank layout (real CoLM): the last rank is master; each group starts with an
! IO rank. Worker 0 is left empty on purpose (numucat = 0).
!
! Environment: RH_OUTDIR (required)  RH_NGROUP  RH_TOTALNUMUCAT  RH_TOTALNPTHOUT
!              RH_TOTALNUMRESV  RH_EMPTY_LAST_GROUP
! Output: shard files RH_OUTDIR/hrtest_hist_unitcat_<day>_seg..._shardNNNNN.nc
! Prints HISTROUTE_OK on the master when every rank got through.
! ---------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_ucid, griducat
   USE MOD_Grid_Reservoir, only: numresv, totalnumresv, resv_global_id
   USE MOD_Grid_RiverLakeHistRoute

   IMPLICIT NONE

   integer, parameter :: NLON = 12, NLAT = 8, NPTHLEV = 3, NDAYS = 2
   integer :: ngroup, totalnpthout, iday, iw, nactive, k, empty_last_group
   character(len=256) :: outdir, fbase
   real(r8), allocatable :: lon(:), lat(:)
   real(r8), allocatable :: ucat_val(:), resv_val(:), zero(:)
   real(r8), allocatable :: mat(:,:), mat_zero(:,:)
   integer,  allocatable :: pth_gid(:), pth_none(:)
   integer,  allocatable :: seq_xy(:)
   integer :: npth_local, itime
   logical :: first

      CALL spmd_init ()

      CALL read_env_str ('RH_OUTDIR', '', outdir)
      IF (len_trim(outdir) == 0) THEN
         IF (p_is_master) write(*,*) 'RH_OUTDIR is required'
         CALL CoLM_stop ()
      ENDIF
      CALL read_env_int ('RH_NGROUP',        1, ngroup)
      CALL read_env_int ('RH_TOTALNUMUCAT', 37, totalnumucat)
      CALL read_env_int ('RH_TOTALNPTHOUT', 11, totalnpthout)
      CALL read_env_int ('RH_TOTALNUMRESV',  5, totalnumresv)
      CALL read_env_int ('RH_EMPTY_LAST_GROUP', 0, empty_last_group)

      CALL divide_processes_into_groups (ngroup)

      ! ---- what the model would have read from the parameter file ----
      DEF_HIST_mode = 'block'
      DEF_HIST_FREQ = 'DAILY'
      DEF_CASE_NAME = 'hrtest'
      DEF_HIST_CompressLevel = 1
      DEF_UnitCatchment_file = trim(outdir) // '/ucat_seq.nc'
      griducat%nlon = NLON
      griducat%nlat = NLAT

      IF (p_is_master) THEN
         allocate (seq_xy(totalnumucat))
         CALL ncio_create_file (trim(DEF_UnitCatchment_file))
         CALL ncio_define_dimension (trim(DEF_UnitCatchment_file), 'seq', totalnumucat)
         seq_xy = (/ (mod(k-1, NLON) + 1, k = 1, totalnumucat) /)
         CALL ncio_write_serial (trim(DEF_UnitCatchment_file), 'seq_x', seq_xy, 'seq')
         seq_xy = (/ (mod((k-1)/NLON, NLAT) + 1, k = 1, totalnumucat) /)
         CALL ncio_write_serial (trim(DEF_UnitCatchment_file), 'seq_y', seq_xy, 'seq')
      ENDIF
      CALL mpi_barrier (p_comm_glb, p_err)

      allocate (lon(NLON), lat(NLAT))
      lon = (/ (real(k,r8), k = 1, NLON) /)
      lat = (/ (real(k,r8), k = 1, NLAT) /)

      CALL build_local_data ()

      DO iday = 1, NDAYS
         write(fbase,'(A,A,I3.3,A)') trim(outdir), '/hrtest_hist_unitcat_', iday, '.nc'
         first = .true.

         ! Every rank, master included, enters the dispatcher.
         CALL route_hist_begin (trim(fbase), (/2000, iday, 0/), first, lon, lat, itime)

         IF (first) CALL route_hist_write_ucat (ucat_val, 'mask_test', &
            longname = 'static mask', units = '1', no_time = .true.)

         CALL route_hist_write_ucat (ucat_val + 1000._r8*iday, 'f_test', &
            longname = 'synthetic unitcat field', units = 'm')
         CALL route_hist_write_resv (resv_val + 1000._r8*iday, 'volresv', &
            longname = 'synthetic reservoir field', units = 'm^3')
         CALL route_hist_write_bif_matrix (mat + 1000._r8*iday, NPTHLEV, npth_local, pth_gid, &
            totalnpthout, 'f_bifflw_lev', longname = 'synthetic pathway flow', units = 'm^3/s')

         CALL route_hist_end ()
      ENDDO

      CALL route_hist_final ()

      CALL mpi_barrier (p_comm_glb, p_err)
      IF (p_is_master) write(*,'(A)') 'HISTROUTE_OK'

      CALL spmd_exit ()

CONTAINS

   SUBROUTINE build_local_data ()

   integer :: iwk, ir, ip

      numucat = 0
      numresv = 0
      npth_local = 0
      allocate (zero(0), mat_zero(NPTHLEV,0), pth_none(0))
      allocate (ucat_val(0), resv_val(0), mat(NPTHLEV,0), pth_gid(0))

      IF (p_is_worker) THEN
         ! worker 0 stays empty; the others deal ids round-robin so consecutive
         ! global ids live on different shards
         nactive = p_np_worker - 1
         IF (empty_last_group /= 0) nactive = 1
         iw = p_iam_worker - 1
         IF (iw >= 0 .and. nactive > 0) THEN
            numucat = count((/ (mod(k-1, nactive) == iw, k = 1, totalnumucat) /))
            deallocate (ucat_val); allocate (ucat_val(numucat))
            IF (allocated(ucat_ucid)) deallocate (ucat_ucid)
            allocate (ucat_ucid(numucat))
            iwk = 0
            DO k = 1, totalnumucat
               IF (mod(k-1, nactive) == iw) THEN
                  iwk = iwk + 1
                  ucat_ucid(iwk) = k
                  ucat_val(iwk) = real(k, r8) + 0.5_r8
               ENDIF
            ENDDO

            numresv = count((/ (mod(k-1, nactive) == iw, k = 1, totalnumresv) /))
            deallocate (resv_val); allocate (resv_val(numresv))
            IF (allocated(resv_global_id)) deallocate (resv_global_id)
            allocate (resv_global_id(numresv))
            ir = 0
            DO k = 1, totalnumresv
               IF (mod(k-1, nactive) == iw) THEN
                  ir = ir + 1
                  resv_global_id(ir) = k
                  resv_val(ir) = real(k, r8) + 0.25_r8
               ENDIF
            ENDDO

            npth_local = count((/ (mod(k-1, nactive) == iw, k = 1, totalnpthout) /))
            deallocate (mat, pth_gid); allocate (mat(NPTHLEV,npth_local), pth_gid(npth_local))
            ip = 0
            DO k = 1, totalnpthout
               IF (mod(k-1, nactive) == iw) THEN
                  ip = ip + 1
                  pth_gid(ip) = k
                  mat(:,ip) = real(k, r8) + (/ 0.1_r8, 0.2_r8, 0.3_r8 /)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF (.not. allocated(ucat_ucid)) allocate (ucat_ucid(0))
      IF (.not. allocated(resv_global_id)) allocate (resv_global_id(0))

   END SUBROUTINE build_local_data

   SUBROUTINE read_env_int (name, default, value)
   character(len=*), intent(in) :: name
   integer, intent(in)  :: default
   integer, intent(out) :: value
   character(len=64) :: buf
   integer :: ierr, st
      value = default
      CALL get_environment_variable (name, buf, status = st)
      IF (st == 0) THEN
         read (buf, *, iostat = ierr) value
         IF (ierr /= 0) value = default
      ENDIF
   END SUBROUTINE read_env_int

   SUBROUTINE read_env_str (name, default, value)
   character(len=*), intent(in)  :: name, default
   character(len=*), intent(out) :: value
   integer :: st
      value = default
      CALL get_environment_variable (name, value, status = st)
      IF (st /= 0) value = default
   END SUBROUTINE read_env_str

END PROGRAM river_hist_route_master_harness
