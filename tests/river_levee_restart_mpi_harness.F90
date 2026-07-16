PROGRAM river_levee_restart_mpi_harness
   USE MOD_Precision, only: r8
   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_USE_LEVEE, DEF_USE_BIFURCATION
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir, only: numresv, totalnumresv, ucat2resv, &
      resv_data_address, reservoir_final
   USE MOD_Grid_RiverLakeLevee, only: has_levee, levsto, levee_init, &
      read_levee_restart, levee_final
   USE MOD_Grid_RiverLakeTimeVars
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE

   integer, parameter :: NCELL_GLOBAL = 6
   character(len=64) :: mode
   character(len=1024) :: input_file, output_file

   CALL get_command_argument(1, mode)
   CALL get_command_argument(2, input_file)
   CALL get_command_argument(3, output_file)
   IF (len_trim(mode) == 0 .or. len_trim(input_file) == 0) THEN
      ERROR STOP 'usage: river_levee_restart_mpi_harness MODE INPUT [OUTPUT]'
   ENDIF

   CALL spmd_init()
   CALL configure_roles()
   DEF_USE_BIFURCATION = .false.

   SELECT CASE (trim(mode))
   CASE ('write-enabled')
      IF (p_np_glb /= 3) CALL CoLM_stop('levee restart write requires 3 ranks')
      DEF_USE_LEVEE = .true.
      CALL configure_network()
      CALL allocate_GridRiverLakeTimeVars()
      CALL initialize_time_state()
      CALL levee_init()
      CALL initialize_levee_state()
      CALL write_GridRiverLakeTimeVars(trim(input_file))
      CALL commit_GridRiverLakeRestart(trim(input_file))

   CASE ('read-enabled', 'read-manifest-disabled', 'read-fold-disabled')
      IF (p_np_glb /= 5) CALL CoLM_stop('levee restart read requires 5 ranks')
      DEF_USE_LEVEE = trim(mode) /= 'read-fold-disabled'
      CALL configure_network()
      CALL allocate_GridRiverLakeTimeVars()
      CALL read_GridRiverLakeTimeVars(trim(input_file))
      CALL levee_init()
      CALL read_levee_restart(trim(input_file), &
         restart_transaction_validated_in = restart_transaction_validated, &
         restart_feature_manifest_present_in = restart_feature_manifest_present, &
         restart_levee_enabled_in = restart_levee_enabled, &
         fold_protected_to_visible = .not. DEF_USE_LEVEE, &
         volwater_ucat_io = volwater_ucat, &
         volwater_ucat_valid_io = volwater_ucat_valid, &
         wdsrf_ucat_in = wdsrf_ucat)
      CALL verify_loaded_state(trim(mode))

      IF (trim(mode) == 'read-enabled') THEN
         IF (len_trim(output_file) == 0) &
            CALL CoLM_stop('levee restart read-enabled requires output file')
         CALL write_GridRiverLakeTimeVars(trim(output_file))
         CALL commit_GridRiverLakeRestart(trim(output_file))
      ENDIF

   CASE DEFAULT
      CALL CoLM_stop('unknown levee restart harness mode')
   END SELECT

   IF (p_is_master) write(*,'(A,1X,A,1X,A,I0,A)') &
      'river levee restart MPI harness: PASS', trim(mode), '(', p_np_glb, ' ranks)'

   CALL deallocate_GridRiverLakeTimeVars()
   CALL levee_final()
   CALL reservoir_final()
   CALL riverlake_network_final()
   CALL spmd_exit()

CONTAINS

   SUBROUTINE configure_roles()
      integer :: global_rank, role_count, worker_rank

      p_is_io = p_iam_glb == p_address_master - 1
      p_is_worker = .not. p_is_master .and. .not. p_is_io

      IF (p_is_worker) THEN
         CALL mpi_comm_split(p_comm_glb, 1, p_iam_glb, p_comm_worker, p_err)
         CALL mpi_comm_rank(p_comm_worker, p_iam_worker, p_err)
      ELSE
         CALL mpi_comm_split(p_comm_glb, MPI_UNDEFINED, p_iam_glb, p_comm_worker, p_err)
         p_iam_worker = -1
      ENDIF
      allocate (p_itis_worker(0:p_np_glb-1))
      CALL mpi_allgather(p_iam_worker, 1, MPI_INTEGER, p_itis_worker, 1, &
         MPI_INTEGER, p_comm_glb, p_err)
      p_np_worker = count(p_itis_worker >= 0)
      allocate (p_address_worker(0:p_np_worker-1))
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
      allocate (p_itis_io(0:p_np_glb-1))
      CALL mpi_allgather(p_iam_io, 1, MPI_INTEGER, p_itis_io, 1, &
         MPI_INTEGER, p_comm_glb, p_err)
      p_np_io = count(p_itis_io >= 0)
      allocate (p_address_io(0:p_np_io-1))
      p_address_io = -1
      DO global_rank = 0, p_np_glb-1
         IF (p_itis_io(global_rank) >= 0) &
            p_address_io(p_itis_io(global_rank)) = global_rank
      ENDDO

      role_count = merge(1, 0, p_is_io .and. .not. p_is_master .and. .not. p_is_worker)
      CALL mpi_allreduce(MPI_IN_PLACE, role_count, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      IF (role_count /= 1 .or. p_np_io /= 1) &
         CALL CoLM_stop('levee restart harness requires exactly one IO-only rank')
   END SUBROUTINE configure_roles


   SUBROUTINE configure_network()
      integer :: i, iw, gid, nlocal

      totalnumucat = NCELL_GLOBAL
      totalnpthout = 0
      npthout_local = 0
      npthlev_bif = 0
      numucat = 0
      IF (p_is_worker) THEN
         DO gid = p_iam_worker + 1, NCELL_GLOBAL, p_np_worker
            numucat = numucat + 1
         ENDDO
      ENDIF

      allocate (ucat_ucid(numucat), x_ucat(numucat), y_ucat(numucat), ucat_next(numucat))
      allocate (topo_rivelv(numucat), topo_rivhgt(numucat), topo_rivlen(numucat))
      allocate (topo_rivman(numucat), topo_rivwth(numucat), topo_rivare(numucat))
      allocate (topo_rivstomax(numucat), topo_area(numucat), topo_fldhgt(4,numucat))
      allocate (lake_type(numucat), levee_frc_data(numucat), levee_hgt_data(numucat))
      allocate (floodplain_curve(numucat), ucat2resv(numucat))

      i = 0
      IF (p_is_worker) THEN
         DO gid = p_iam_worker + 1, NCELL_GLOBAL, p_np_worker
            i = i + 1
            ucat_ucid(i) = gid
            x_ucat(i) = 300 + gid
            y_ucat(i) = 400 + 2*gid
            ucat_next(i) = merge(gid + 1, 0, gid < NCELL_GLOBAL)
            topo_rivelv(i) = 0._r8
            topo_rivhgt(i) = 2._r8
            topo_rivlen(i) = 100._r8
            topo_rivman(i) = 0.03_r8
            topo_rivwth(i) = 10._r8
            topo_rivare(i) = 1000._r8
            topo_rivstomax(i) = 2000._r8
            topo_area(i) = 4000._r8
            topo_fldhgt(:,i) = [0.5_r8, 1._r8, 1.5_r8, 2._r8]
            lake_type(i) = 0
            levee_frc_data(i) = 0.5_r8
            levee_hgt_data(i) = 1.25_r8
            ucat2resv(i) = 0
            CALL configure_curve(i)
         ENDDO
      ENDIF

      allocate (ucat_data_address(0:p_np_worker-1))
      allocate (resv_data_address(0:p_np_worker-1))
      DO iw = 0, p_np_worker-1
         nlocal = 0
         DO gid = iw + 1, NCELL_GLOBAL, p_np_worker
            nlocal = nlocal + 1
         ENDDO
         allocate (ucat_data_address(iw)%val(nlocal))
         allocate (resv_data_address(iw)%val(0))
         i = 0
         DO gid = iw + 1, NCELL_GLOBAL, p_np_worker
            i = i + 1
            ucat_data_address(iw)%val(i) = gid
         ENDDO
      ENDDO
      numresv = 0
      totalnumresv = 0
   END SUBROUTINE configure_network


   SUBROUTINE configure_curve(i)
      integer, intent(in) :: i
      integer :: level

      floodplain_curve(i)%nlfp = 4
      floodplain_curve(i)%rivhgt = 2._r8
      floodplain_curve(i)%rivare = 1000._r8
      floodplain_curve(i)%rivstomax = 2000._r8
      allocate (floodplain_curve(i)%flphgt(0:4))
      allocate (floodplain_curve(i)%flparea(0:4))
      allocate (floodplain_curve(i)%flpaccare(0:4))
      allocate (floodplain_curve(i)%flpstomax(0:4))
      floodplain_curve(i)%flphgt = [0._r8, 0.5_r8, 1._r8, 1.5_r8, 2._r8]
      floodplain_curve(i)%flparea = [0._r8, 1000._r8, 1000._r8, 1000._r8, 1000._r8]
      floodplain_curve(i)%flpaccare(0) = 0._r8
      floodplain_curve(i)%flpstomax(0) = 0._r8
      DO level = 1, 4
         floodplain_curve(i)%flpaccare(level) = floodplain_curve(i)%flpaccare(level-1) &
            + floodplain_curve(i)%flparea(level)
         floodplain_curve(i)%flpstomax(level) = floodplain_curve(i)%flpstomax(level-1) &
            + 0.5_r8 * (floodplain_curve(i)%flpaccare(level) &
            + floodplain_curve(i)%flpaccare(level-1)) &
            * (floodplain_curve(i)%flphgt(level)-floodplain_curve(i)%flphgt(level-1))
      ENDDO
   END SUBROUTINE configure_curve


   SUBROUTINE initialize_time_state()
      integer :: i, gid

      acctime_rnof = 23.5_r8
      DO i = 1, numucat
         gid = ucat_ucid(i)
         wdsrf_ucat(i) = floodplain_curve(i)%depth(expected_total(gid))
         veloc_riv(i) = real(gid, r8) / 20._r8
         acc_rnof_uc(i) = real(gid, r8) / 1000._r8
         volwater_ucat(i) = expected_visible(gid)
      ENDDO
      volwater_ucat_valid = .true.
   END SUBROUTINE initialize_time_state


   SUBROUTINE initialize_levee_state()
      integer :: i, gid

      DO i = 1, numucat
         gid = ucat_ucid(i)
         IF (.not. has_levee(i)) CALL CoLM_stop('test levee geometry did not create a levee')
         levsto(i) = expected_protected(gid)
      ENDDO
   END SUBROUTINE initialize_levee_state


   SUBROUTINE verify_loaded_state(read_mode)
      character(len=*), intent(in) :: read_mode
      integer :: failures, i, gid
      real(r8) :: actual_total, expected_global_total, expected_local_total, expected_volume

      failures = 0
      actual_total = 0._r8
      expected_local_total = 0._r8
      IF (.not. volwater_ucat_valid) failures = failures + 1

      DO i = 1, numucat
         gid = ucat_ucid(i)
         SELECT CASE (trim(read_mode))
         CASE ('read-enabled')
            expected_volume = expected_visible(gid)
            IF (abs(levsto(i) - expected_protected(gid)) > 1.e-10_r8) failures = failures + 1
            expected_local_total = expected_local_total + expected_total(gid)
         CASE ('read-manifest-disabled')
            expected_volume = expected_visible(gid)
            IF (levsto(i) /= 0._r8) failures = failures + 1
            expected_local_total = expected_local_total + expected_visible(gid)
         CASE ('read-fold-disabled')
            expected_volume = expected_total(gid)
            IF (levsto(i) /= 0._r8) failures = failures + 1
            expected_local_total = expected_local_total + expected_total(gid)
         CASE DEFAULT
            CALL CoLM_stop('unknown levee verification mode')
         END SELECT

         IF (.not. ieee_is_finite(volwater_ucat(i))) failures = failures + 1
         IF (volwater_ucat(i) < 0._r8) failures = failures + 1
         IF (abs(volwater_ucat(i) - expected_volume) > 1.e-10_r8) failures = failures + 1
         actual_total = actual_total + volwater_ucat(i) + levsto(i)
      ENDDO

      expected_global_total = expected_local_total
      CALL mpi_allreduce(MPI_IN_PLACE, actual_total, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, expected_global_total, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
      CALL mpi_allreduce(MPI_IN_PLACE, failures, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      IF (abs(actual_total - expected_global_total) > 1.e-9_r8) failures = failures + 1
      IF (failures /= 0) THEN
         IF (p_is_master) write(*,'(A,I0,A,2ES18.8)') &
            'levee restart verification failures=', failures, ', totals=', &
            actual_total, expected_global_total
         CALL CoLM_stop('levee restart state mismatch')
      ENDIF
   END SUBROUTINE verify_loaded_state


   real(r8) FUNCTION expected_visible(gid)
      integer, intent(in) :: gid
      expected_visible = 5000._r8 + 17._r8 * real(gid, r8)
   END FUNCTION expected_visible


   real(r8) FUNCTION expected_protected(gid)
      integer, intent(in) :: gid
      expected_protected = 250._r8 + 3._r8 * real(gid, r8)
   END FUNCTION expected_protected


   real(r8) FUNCTION expected_total(gid)
      integer, intent(in) :: gid
      expected_total = expected_visible(gid) + expected_protected(gid)
   END FUNCTION expected_total

END PROGRAM river_levee_restart_mpi_harness
