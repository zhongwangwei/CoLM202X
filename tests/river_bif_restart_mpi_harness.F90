PROGRAM river_bif_restart_mpi_harness
   USE MOD_Precision, only: r8
   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_USE_LEVEE, DEF_USE_BIFURCATION
   USE MOD_WorkerPushData, only: build_worker_pushdata
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir, only: numresv, totalnumresv, ucat2resv, &
      resv_data_address, reservoir_final
   USE MOD_Grid_RiverLakeLevee, only: levee_init, levee_final
   USE MOD_Grid_RiverLakeBifurcation, only: bifurcation_init, &
      bifurcation_calc, bifurcation_invalidate_static_dn, &
      read_bifurcation_restart, bifurcation_final
   USE MOD_Grid_RiverLakeTimeVars
   IMPLICIT NONE

   integer, parameter :: NCELL_GLOBAL = 6
   character(len=64) :: mode
   character(len=1024) :: input_file, output_file
   logical :: expect_loaded

   CALL get_command_argument(1, mode)
   CALL get_command_argument(2, input_file)
   CALL get_command_argument(3, output_file)
   IF (len_trim(mode) == 0 .or. len_trim(input_file) == 0) THEN
      ERROR STOP 'usage: river_bif_restart_mpi_harness MODE INPUT [OUTPUT]'
   ENDIF

   CALL spmd_init()
   CALL configure_roles()
   DEF_USE_LEVEE = .false.
   DEF_USE_BIFURCATION = .true.

   SELECT CASE (trim(mode))
   CASE ('write')
      IF (p_np_glb /= 3) CALL CoLM_stop('BIF restart write requires 3 ranks')
      CALL configure_network()
      CALL allocate_GridRiverLakeTimeVars()
      CALL initialize_time_state()
      CALL levee_init()
      CALL bifurcation_init()
      CALL generate_nonzero_bif_state()
      CALL write_GridRiverLakeTimeVars(trim(input_file))
      CALL commit_GridRiverLakeRestart(trim(input_file))

   CASE ('read-valid', 'read-cold')
      IF (p_np_glb /= 5) CALL CoLM_stop('BIF restart read requires 5 ranks')
      IF (len_trim(output_file) == 0) CALL CoLM_stop('BIF restart read requires output file')
      expect_loaded = trim(mode) == 'read-valid'
      CALL configure_network()
      CALL allocate_GridRiverLakeTimeVars()
      CALL read_GridRiverLakeTimeVars(trim(input_file))
      CALL levee_init()
      CALL bifurcation_init()
      CALL load_bif_state(expect_loaded)
      CALL write_GridRiverLakeTimeVars(trim(output_file))
      CALL commit_GridRiverLakeRestart(trim(output_file))

   CASE DEFAULT
      CALL CoLM_stop('unknown BIF restart harness mode')
   END SELECT

   IF (p_is_master) write(*,'(A,1X,A,1X,A,I0,A)') &
      'river BIF restart MPI harness: PASS', trim(mode), '(', p_np_glb, ' ranks)'

   CALL deallocate_GridRiverLakeTimeVars()
   CALL bifurcation_final()
   CALL levee_final()
   CALL reservoir_final()
   CALL riverlake_network_final()
   CALL spmd_exit()

CONTAINS

   SUBROUTINE configure_roles()
      integer :: global_rank, role_count, worker_rank

      ! Reserve the rank below master as a true IO-only participant.  It owns
      ! no river/BIF state but must enter all global restart collectives.
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
         CALL CoLM_stop('BIF restart harness requires exactly one IO-only rank')
   END SUBROUTINE configure_roles


   SUBROUTINE configure_network()
      integer :: i, iw, gid, down_gid, nlocal

      totalnumucat = NCELL_GLOBAL
      totalnpthout = NCELL_GLOBAL
      npthlev_bif = 2
      numucat = 0
      IF (p_is_worker) THEN
         DO gid = p_iam_worker + 1, NCELL_GLOBAL, p_np_worker
            numucat = numucat + 1
         ENDDO
      ENDIF
      npthout_local = numucat

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
            x_ucat(i) = 100 + gid
            y_ucat(i) = 200 + 2*gid
            ucat_next(i) = modulo(gid, NCELL_GLOBAL) + 1
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
            levee_frc_data(i) = 1._r8
            levee_hgt_data(i) = 0._r8
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

      allocate (pth_upst_local(npthout_local), pth_down_local(npthout_local))
      allocate (pth_down_ucid(npthout_local), pth_global_id(npthout_local))
      allocate (pth_dst(npthout_local), pth_elv(2,npthout_local), pth_wth(2,npthout_local))
      allocate (pth_man(2), bif_incoming_pths(1,numucat), bif_incoming_wts(1,numucat))
      pth_man = [0.03_r8, 0.05_r8]
      DO i = 1, npthout_local
         gid = ucat_ucid(i)
         down_gid = modulo(gid, NCELL_GLOBAL) + 1
         pth_upst_local(i) = i
         pth_down_ucid(i) = down_gid
         pth_down_local(i) = local_index(down_gid)
         pth_global_id(i) = gid
         pth_dst(i) = 1000._r8 + real(gid, r8)
         pth_elv(:,i) = [0._r8, 2.5_r8]
         pth_wth(:,i) = [10._r8, 5._r8]
         bif_incoming_pths(1,i) = modulo(gid-2+NCELL_GLOBAL, NCELL_GLOBAL) + 1
         bif_incoming_wts(1,i) = 1._r8
      ENDDO
      max_bif_incoming = 1

      IF (p_is_worker) THEN
         CALL build_worker_pushdata(numucat, ucat_ucid, npthout_local, &
            pth_down_ucid, push_bif_dn2pth)
         CALL build_worker_pushdata(npthout_local, pth_global_id, numucat, &
            bif_incoming_pths, bif_incoming_wts, push_bif_influx)
      ENDIF
   END SUBROUTINE configure_network


   integer FUNCTION local_index(gid)
      integer, intent(in) :: gid
      integer :: i

      local_index = -1
      DO i = 1, numucat
         IF (ucat_ucid(i) == gid) THEN
            local_index = i
            RETURN
         ENDIF
      ENDDO
   END FUNCTION local_index


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

      acctime_rnof = 37.5_r8
      DO i = 1, numucat
         gid = ucat_ucid(i)
         wdsrf_ucat(i) = expected_current(gid)
         wdsrf_ucat_prev(i) = expected_previous(gid)
         veloc_riv(i) = -real(gid, r8)/10._r8
         acc_rnof_uc(i) = real(gid, r8)/1000._r8
         volwater_ucat(i) = floodplain_curve(i)%volume(wdsrf_ucat(i))
      ENDDO
      wdsrf_ucat_prev_valid = .true.
      volwater_ucat_valid = .true.
   END SUBROUTINE initialize_time_state


   SUBROUTINE generate_nonzero_bif_state()
      logical, allocatable :: active(:), is_reservoir(:)
      integer, allocatable :: river_system(:)
      real(r8), allocatable :: dt(:), normal_outgoing(:), empty_reservoir(:)

      IF (.not. p_is_worker) RETURN
      allocate (active(numucat), is_reservoir(numucat), river_system(numucat))
      allocate (dt(1), normal_outgoing(numucat), empty_reservoir(0))
      active = .true.
      is_reservoir = .false.
      river_system = 1
      dt = 30._r8
      normal_outgoing = 0._r8
      CALL bifurcation_invalidate_static_dn()
      CALL bifurcation_calc(wdsrf_ucat, wdsrf_ucat_prev, volwater_ucat, &
         volwater_ucat_valid, empty_reservoir, is_reservoir, dt, &
         river_system, active, normal_outgoing)
      deallocate (active, is_reservoir, river_system, dt, normal_outgoing, empty_reservoir)
   END SUBROUTINE generate_nonzero_bif_state


   SUBROUTINE load_bif_state(expect_restart_loaded)
      logical, intent(in) :: expect_restart_loaded
      logical :: restart_loaded
      integer :: failures

      CALL read_bifurcation_restart(trim(input_file), &
         wdsrf_ucat_prev_restart_found, restart_loaded, &
         restart_transaction_validated, restart_feature_manifest_present, &
         restart_bifurcation_enabled)
      failures = merge(0, 1, restart_loaded .eqv. expect_restart_loaded)
      IF (.not. restart_loaded) THEN
         wdsrf_ucat_prev = wdsrf_ucat
         wdsrf_ucat_prev_valid = .true.
      ENDIF
      CALL mpi_allreduce(MPI_IN_PLACE, failures, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      IF (failures /= 0) CALL CoLM_stop('unexpected BIF restart_loaded result')
   END SUBROUTINE load_bif_state


   real(r8) FUNCTION expected_current(gid)
      integer, intent(in) :: gid
      expected_current = merge(1.8_r8, 0.7_r8, modulo(gid, 2) == 1)
   END FUNCTION expected_current


   real(r8) FUNCTION expected_previous(gid)
      integer, intent(in) :: gid
      expected_previous = 0.5_r8 * expected_current(gid) + 0.05_r8*real(gid, r8)
   END FUNCTION expected_previous

END PROGRAM river_bif_restart_mpi_harness
