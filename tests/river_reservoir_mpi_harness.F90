PROGRAM river_reservoir_mpi_harness
   USE MOD_Precision, only: r8
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Namelist, only: DEF_ReservoirPara_file
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir
   IMPLICIT NONE

   integer :: failures, failures_global

   CALL get_command_argument(1, DEF_ReservoirPara_file)
   IF (len_trim(DEF_ReservoirPara_file) == 0) &
      ERROR STOP 'usage: river_reservoir_mpi_harness PARAMETER_FILE'

   CALL spmd_init ()
   CALL configure_roles ()
   IF (p_np_glb /= 5 .or. p_np_worker /= 3) &
      CALL CoLM_stop ('reservoir harness requires 5 ranks and 3 workers')

   IF (p_is_master) CALL write_parameter_file (trim(DEF_ReservoirPara_file))
   CALL mpi_barrier (p_comm_glb, p_err)

   CALL setup_local_network ()
   CALL reservoir_init ()

   failures = 0
   IF (totalnumresv /= 2) failures = failures + 1
   IF (.not. allocated(dam_GRAND_ID)) THEN
      failures = failures + 1
   ELSEIF (size(dam_GRAND_ID) /= 2) THEN
      failures = failures + 1
   ELSEIF (any(dam_GRAND_ID /= [1001, 1003])) THEN
      failures = failures + 1
   ENDIF

   IF (p_is_worker) CALL check_worker (failures)
   IF (p_is_master) CALL check_master_addresses (failures)
   CALL mpi_allreduce (failures, failures_global, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)

   IF (p_is_master) THEN
      IF (failures_global == 0) THEN
         write(*,'(A)') 'river reservoir sparse-catalogue MPI harness: PASS'
      ELSE
         write(*,'(A,I0)') 'river reservoir sparse-catalogue MPI harness: FAIL count=', failures_global
      ENDIF
   ENDIF

   CALL reservoir_final ()
   IF (allocated(ucat_ucid)) deallocate(ucat_ucid)
   IF (allocated(lake_type)) deallocate(lake_type)
   IF (failures_global /= 0) CALL CoLM_stop ('reservoir sparse-catalogue MPI regression')
   CALL spmd_exit ()

CONTAINS

   SUBROUTINE configure_roles ()
      integer :: global_rank, worker_rank

      p_is_io = p_iam_glb == p_address_master - 1
      p_is_worker = .not. p_is_master .and. .not. p_is_io
      IF (p_is_worker) THEN
         CALL mpi_comm_split (p_comm_glb, 1, p_iam_glb, p_comm_worker, p_err)
         CALL mpi_comm_rank (p_comm_worker, p_iam_worker, p_err)
      ELSE
         CALL mpi_comm_split (p_comm_glb, MPI_UNDEFINED, p_iam_glb, p_comm_worker, p_err)
         p_iam_worker = -1
      ENDIF

      allocate (p_itis_worker(0:p_np_glb-1))
      CALL mpi_allgather (p_iam_worker, 1, MPI_INTEGER, p_itis_worker, 1, &
         MPI_INTEGER, p_comm_glb, p_err)
      p_np_worker = count(p_itis_worker >= 0)
      allocate (p_address_worker(0:p_np_worker-1))
      p_address_worker = -1
      DO global_rank = 0, p_np_glb-1
         worker_rank = p_itis_worker(global_rank)
         IF (worker_rank >= 0) p_address_worker(worker_rank) = global_rank
      ENDDO
   END SUBROUTINE configure_roles

   SUBROUTINE write_parameter_file (filename)
      character(len=*), intent(in) :: filename
      logical :: exists

      inquire(file=filename, exist=exists)
      IF (exists) THEN
         OPEN(unit=97, file=filename, status='old')
         CLOSE(unit=97, status='delete')
      ENDIF
      CALL ncio_create_file (filename)
      CALL ncio_define_dimension (filename, 'reservoir', 5)
      CALL ncio_write_serial (filename, 'dam_GRAND_ID', [1001,1002,1003,1004,1005], 'reservoir')
      CALL ncio_write_serial (filename, 'dam_seq',      [100,200,300,400,500], 'reservoir')
      CALL ncio_write_serial (filename, 'dam_year',     [1901,1902,1903,1904,1905], 'reservoir')
      CALL ncio_write_serial (filename, 'dam_TotalVol_mcm', [10,20,30,40,50]*1._r8, 'reservoir')
      CALL ncio_write_serial (filename, 'dam_ConVol_mcm',   [5,10,15,20,25]*1._r8, 'reservoir')
      CALL ncio_write_serial (filename, 'dam_Qn',           [1,2,3,4,5]*1._r8, 'reservoir')
      CALL ncio_write_serial (filename, 'dam_Qf',           [10,20,30,40,50]*1._r8, 'reservoir')
   END SUBROUTINE write_parameter_file

   SUBROUTINE setup_local_network ()
      numucat = merge(1, 0, p_is_worker)
      allocate (ucat_ucid(numucat), lake_type(numucat))
      lake_type = 0
      IF (.not. p_is_worker) RETURN
      SELECT CASE (p_iam_worker)
      CASE (0)
         ucat_ucid = 300
      CASE (1)
         ucat_ucid = 100
      CASE DEFAULT
         ucat_ucid = 999
      END SELECT
   END SUBROUTINE setup_local_network

   SUBROUTINE check_worker (nfail)
      integer, intent(inout) :: nfail
      integer :: expected_gid, expected_row

      SELECT CASE (p_iam_worker)
      CASE (0)
         expected_gid = 2; expected_row = 3
      CASE (1)
         expected_gid = 1; expected_row = 1
      CASE DEFAULT
         IF (numresv /= 0 .or. lake_type(1) /= 0) nfail = nfail + 1
         IF (.not. allocated(resv_global_id)) nfail = nfail + 1
         IF (allocated(resv_global_id)) THEN
            IF (size(resv_global_id) /= 0) nfail = nfail + 1
         ENDIF
         RETURN
      END SELECT

      IF (numresv /= 1 .or. lake_type(1) /= 2) nfail = nfail + 1
      IF (size(resv_global_id) /= 1 .or. resv_global_id(1) /= expected_gid) nfail = nfail + 1
      IF (resv_ucid(1) /= ucat_ucid(1)) nfail = nfail + 1
      IF (dam_build_year(1) /= 1900 + expected_row) nfail = nfail + 1
      IF (abs(volresv_total(1) - real(10*expected_row, r8)*1.e6_r8) > 1._r8) nfail = nfail + 1
      IF (abs(volresv_normal(1) - real(5*expected_row, r8)*1.e6_r8) > 1._r8) nfail = nfail + 1
      IF (qresv_normal(1) /= real(expected_row, r8)) nfail = nfail + 1
      IF (qresv_flood(1) /= real(10*expected_row, r8)) nfail = nfail + 1
   END SUBROUTINE check_worker

   SUBROUTINE check_master_addresses (nfail)
      integer, intent(inout) :: nfail

      IF (.not. allocated(resv_data_address)) THEN
         nfail = nfail + 1
         RETURN
      ENDIF
      IF (.not. allocated(resv_data_address(0)%val)) nfail = nfail + 1
      IF (.not. allocated(resv_data_address(1)%val)) nfail = nfail + 1
      IF (allocated(resv_data_address(0)%val)) THEN
         IF (size(resv_data_address(0)%val) /= 1 .or. resv_data_address(0)%val(1) /= 2) nfail = nfail + 1
      ENDIF
      IF (allocated(resv_data_address(1)%val)) THEN
         IF (size(resv_data_address(1)%val) /= 1 .or. resv_data_address(1)%val(1) /= 1) nfail = nfail + 1
      ENDIF
      IF (allocated(resv_data_address(2)%val)) nfail = nfail + 1
   END SUBROUTINE check_master_addresses

END PROGRAM river_reservoir_mpi_harness
