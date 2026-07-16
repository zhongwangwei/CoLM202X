PROGRAM river_restart_mpi_harness
   USE MOD_Precision, only: r8
   USE MOD_SPMD_Task
   USE MOD_DataType, only: pointer_int32_1d
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir
   USE MOD_Grid_RiverLakeLevee, only: levsto
   USE MOD_Grid_RiverLakeTimeVars
   IMPLICIT NONE

   character(len=64) :: mode
   character(len=1024) :: restart_file
   integer :: expected_ranks

   CALL get_command_argument(1, mode)
   CALL get_command_argument(2, restart_file)
   IF (len_trim(mode) == 0 .or. len_trim(restart_file) == 0) THEN
      ERROR STOP 'usage: river_restart_mpi_harness MODE RESTART_FILE'
   ENDIF

   CALL spmd_init()
   CALL configure_roles()

   SELECT CASE (trim(mode))
   CASE ('write', 'write-incomplete', 'write-empty')
      expected_ranks = 3
   CASE ('read', 'read-legacy-prev', 'read-empty')
      expected_ranks = 5
   CASE DEFAULT
      CALL CoLM_stop('unknown river restart harness mode')
   END SELECT
   IF (p_np_glb /= expected_ranks) CALL CoLM_stop('unexpected MPI rank count for restart harness mode')

   SELECT CASE (trim(mode))
   CASE ('write')
      CALL run_write(trim(restart_file), 6, .true.)
   CASE ('write-incomplete')
      CALL run_write(trim(restart_file), 6, .false.)
   CASE ('write-empty')
      CALL run_write(trim(restart_file), 0, .true.)
   CASE ('read')
      CALL run_read(trim(restart_file), 6, .false.)
   CASE ('read-legacy-prev')
      CALL run_read(trim(restart_file), 6, .false.)
   CASE ('read-empty')
      CALL run_read(trim(restart_file), 0, .false.)
   END SELECT

   IF (p_is_master) write(*,'(A,1X,A,1X,A,I0,A)') &
      'river restart MPI harness: PASS', trim(mode), '(', p_np_glb, ' ranks)'
   CALL cleanup_test_state()
   CALL spmd_exit()

CONTAINS

   SUBROUTINE configure_roles()
      integer :: global_rank, role_count, worker_rank

      ! Keep one rank independent of both the master and worker roles.  This
      ! rank must still enter every global restart collective and barrier.
      p_is_io = p_iam_glb == p_address_master - 1
      p_is_worker = .not. p_is_master .and. .not. p_is_io

      IF (p_is_worker) THEN
         CALL mpi_comm_split(p_comm_glb, 1, p_iam_glb, p_comm_worker, p_err)
         CALL mpi_comm_rank(p_comm_worker, p_iam_worker, p_err)
      ELSE
         CALL mpi_comm_split(p_comm_glb, MPI_UNDEFINED, p_iam_glb, p_comm_worker, p_err)
         p_iam_worker = -1
      ENDIF
      allocate(p_itis_worker(0:p_np_glb-1))
      CALL mpi_allgather(p_iam_worker, 1, MPI_INTEGER, p_itis_worker, 1, &
         MPI_INTEGER, p_comm_glb, p_err)
      p_np_worker = count(p_itis_worker >= 0)
      allocate(p_address_worker(0:p_np_worker-1))
      p_address_worker = -1
      DO global_rank = 0, p_np_glb-1
         worker_rank = p_itis_worker(global_rank)
         IF (worker_rank >= 0) p_address_worker(worker_rank) = global_rank
      ENDDO

      IF (p_is_io) THEN
         CALL mpi_comm_split(p_comm_glb, 1, p_iam_glb, p_comm_io, p_err)
         CALL mpi_comm_rank(p_comm_io, p_iam_io, p_err)
      ELSE
         CALL mpi_comm_split(p_comm_glb, MPI_UNDEFINED, p_iam_glb, p_comm_io, p_err)
         p_iam_io = -1
      ENDIF
      allocate(p_itis_io(0:p_np_glb-1))
      CALL mpi_allgather(p_iam_io, 1, MPI_INTEGER, p_itis_io, 1, &
         MPI_INTEGER, p_comm_glb, p_err)
      p_np_io = count(p_itis_io >= 0)
      allocate(p_address_io(0:p_np_io-1))
      p_address_io = -1
      DO global_rank = 0, p_np_glb-1
         IF (p_itis_io(global_rank) >= 0) &
            p_address_io(p_itis_io(global_rank)) = global_rank
      ENDDO

      role_count = merge(1, 0, p_is_io .and. .not. p_is_master .and. .not. p_is_worker)
      CALL mpi_allreduce(MPI_IN_PLACE, role_count, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      IF (role_count /= 1 .or. p_np_io /= 1) &
         CALL CoLM_stop('restart harness requires exactly one IO-only rank')
   END SUBROUTINE configure_roles

   SUBROUTINE setup_test_network(ntotal)
      integer, intent(in) :: ntotal
      integer :: i, iw, gid, nlocal

      totalnumucat = ntotal
      numucat = 0
      IF (p_is_worker) THEN
         DO gid = p_iam_worker + 1, ntotal, p_np_worker
            numucat = numucat + 1
         ENDDO
      ENDIF

      allocate(ucat_ucid(numucat), x_ucat(numucat), y_ucat(numucat), ucat_next(numucat))
      i = 0
      IF (p_is_worker) THEN
         DO gid = p_iam_worker + 1, ntotal, p_np_worker
            i = i + 1
            ucat_ucid(i) = gid
            x_ucat(i) = 100 + gid
            y_ucat(i) = 200 + 2 * gid
            ucat_next(i) = merge(gid + 1, 0, gid < ntotal)
         ENDDO
      ENDIF

      allocate(ucat_data_address(0:p_np_worker-1))
      allocate(resv_data_address(0:p_np_worker-1))
      DO iw = 0, p_np_worker-1
         nlocal = 0
         DO gid = iw + 1, ntotal, p_np_worker
            nlocal = nlocal + 1
         ENDDO
         allocate(ucat_data_address(iw)%val(nlocal))
         i = 0
         DO gid = iw + 1, ntotal, p_np_worker
            i = i + 1
            ucat_data_address(iw)%val(i) = gid
         ENDDO
      ENDDO

      numresv = 0
      totalnumresv = merge(4, 0, ntotal > 0)
      IF (p_is_worker) THEN
         DO gid = p_iam_worker + 1, totalnumresv, p_np_worker
            numresv = numresv + 1
         ENDDO
      ENDIF
      DO iw = 0, p_np_worker-1
         nlocal = 0
         DO gid = iw + 1, totalnumresv, p_np_worker
            nlocal = nlocal + 1
         ENDDO
         allocate(resv_data_address(iw)%val(nlocal))
         i = 0
         DO gid = iw + 1, totalnumresv, p_np_worker
            i = i + 1
            resv_data_address(iw)%val(i) = gid
         ENDDO
      ENDDO
      allocate(resv_global_id(numresv), resv_ucid(numresv))
      IF (p_is_worker) THEN
         DO i = 1, numresv
            gid = p_iam_worker + 1 + (i - 1) * p_np_worker
            resv_global_id(i) = gid
            resv_ucid(i) = 10000 + 7 * gid
         ENDDO
      ENDIF
      allocate(levsto(numucat))
      levsto = 0._r8
   END SUBROUTINE setup_test_network

   SUBROUTINE initialize_expected_state()
      integer :: i, gid

      acctime_rnof = 45.5_r8
      DO i = 1, numucat
         gid = ucat_ucid(i)
         wdsrf_ucat(i) = expected_current(gid)
         wdsrf_ucat_prev(i) = expected_previous(gid)
         veloc_riv(i) = expected_velocity(gid)
         acc_rnof_uc(i) = expected_runoff(gid)
         volwater_ucat(i) = expected_volume(gid)
      ENDDO
      DO i = 1, numresv
         gid = p_iam_worker + 1 + (i - 1) * p_np_worker
         volresv(i) = expected_reservoir_volume(gid)
      ENDDO
      wdsrf_ucat_prev_valid = .true.
      volwater_ucat_valid = .true.
   END SUBROUTINE initialize_expected_state

   SUBROUTINE run_write(file_restart, ntotal, commit_restart)
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: ntotal
      logical, intent(in) :: commit_restart

      CALL setup_test_network(ntotal)
      CALL allocate_GridRiverLakeTimeVars()
      CALL initialize_expected_state()
      CALL write_GridRiverLakeTimeVars(file_restart)
      IF (commit_restart) CALL commit_GridRiverLakeRestart(file_restart)
   END SUBROUTINE run_write

   SUBROUTINE run_read(file_restart, ntotal, expect_distinct_previous)
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: ntotal
      logical, intent(in) :: expect_distinct_previous
      integer :: failures, i, gid

      CALL setup_test_network(ntotal)
      CALL allocate_GridRiverLakeTimeVars()
      CALL read_GridRiverLakeTimeVars(file_restart)

      failures = 0
      IF (abs(acctime_rnof - 45.5_r8) > 1.e-12_r8) failures = failures + 1
      IF (.not. wdsrf_ucat_prev_valid) failures = failures + 1
      IF (ntotal > 0 .and. .not. volwater_ucat_valid) failures = failures + 1
      IF (ntotal == 0 .and. volwater_ucat_valid) failures = failures + 1
      IF (size(wdsrf_ucat) /= numucat .or. size(wdsrf_ucat_prev) /= numucat) failures = failures + 1
      IF (size(volresv) /= numresv) failures = failures + 1
      IF (ntotal > 0 .and. totalnumresv /= 4) failures = failures + 1
      IF (ntotal == 0 .and. (totalnumresv /= 0 .or. numresv /= 0)) failures = failures + 1
      IF (p_is_io .and. size(volresv) /= 0) failures = failures + 1

      DO i = 1, numucat
         gid = ucat_ucid(i)
         IF (abs(wdsrf_ucat(i) - expected_current(gid)) > 1.e-12_r8) failures = failures + 1
         IF (expect_distinct_previous) THEN
            IF (abs(wdsrf_ucat_prev(i) - expected_previous(gid)) > 1.e-12_r8) failures = failures + 1
         ELSE
            IF (abs(wdsrf_ucat_prev(i) - expected_current(gid)) > 1.e-12_r8) failures = failures + 1
         ENDIF
         IF (abs(veloc_riv(i) - expected_velocity(gid)) > 1.e-12_r8) failures = failures + 1
         IF (abs(acc_rnof_uc(i) - expected_runoff(gid)) > 1.e-12_r8) failures = failures + 1
         IF (abs(volwater_ucat(i) - expected_volume(gid)) > 1.e-12_r8) failures = failures + 1
      ENDDO
      DO i = 1, numresv
         gid = p_iam_worker + 1 + (i - 1) * p_np_worker
         IF (abs(volresv(i) - expected_reservoir_volume(gid)) > 1.e-12_r8) failures = failures + 1
      ENDDO

      CALL mpi_allreduce(MPI_IN_PLACE, failures, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      IF (failures /= 0) THEN
         IF (p_is_master) write(*,'(A,I0)') 'river restart verification failures=', failures
         CALL CoLM_stop('river restart round-trip mismatch')
      ENDIF
   END SUBROUTINE run_read

   real(r8) FUNCTION expected_current(gid)
      integer, intent(in) :: gid
      expected_current = 1000._r8 + real(gid, r8)
   END FUNCTION expected_current

   real(r8) FUNCTION expected_previous(gid)
      integer, intent(in) :: gid
      expected_previous = 2000._r8 + 3._r8 * real(gid, r8)
   END FUNCTION expected_previous

   real(r8) FUNCTION expected_velocity(gid)
      integer, intent(in) :: gid
      expected_velocity = -real(gid, r8) / 8._r8
   END FUNCTION expected_velocity

   real(r8) FUNCTION expected_runoff(gid)
      integer, intent(in) :: gid
      expected_runoff = real(gid, r8) / 1000._r8
   END FUNCTION expected_runoff

   real(r8) FUNCTION expected_volume(gid)
      integer, intent(in) :: gid
      expected_volume = 50000._r8 + 11._r8 * real(gid, r8)
   END FUNCTION expected_volume

   real(r8) FUNCTION expected_reservoir_volume(gid)
      integer, intent(in) :: gid
      expected_reservoir_volume = 700000._r8 + 101._r8 * real(gid, r8)
   END FUNCTION expected_reservoir_volume

   SUBROUTINE cleanup_test_state()
      integer :: i

      CALL deallocate_GridRiverLakeTimeVars()
      IF (allocated(ucat_ucid)) deallocate(ucat_ucid)
      IF (allocated(x_ucat)) deallocate(x_ucat)
      IF (allocated(y_ucat)) deallocate(y_ucat)
      IF (allocated(ucat_next)) deallocate(ucat_next)
      IF (allocated(levsto)) deallocate(levsto)
      IF (allocated(ucat_data_address)) THEN
         DO i = lbound(ucat_data_address, 1), ubound(ucat_data_address, 1)
            IF (allocated(ucat_data_address(i)%val)) deallocate(ucat_data_address(i)%val)
         ENDDO
         deallocate(ucat_data_address)
      ENDIF
      IF (allocated(resv_data_address)) THEN
         DO i = lbound(resv_data_address, 1), ubound(resv_data_address, 1)
            IF (allocated(resv_data_address(i)%val)) deallocate(resv_data_address(i)%val)
         ENDDO
         deallocate(resv_data_address)
      ENDIF
      IF (allocated(resv_global_id)) deallocate(resv_global_id)
      IF (allocated(resv_ucid)) deallocate(resv_ucid)
   END SUBROUTINE cleanup_test_state

END PROGRAM river_restart_mpi_harness
