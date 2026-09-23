! Production tracer_substep against real MPI push mappings. Run with 2-4 ranks.
PROGRAM river_bif_tracer_harness
   USE mpi, only: MPI_Init, MPI_Finalize, MPI_Comm_rank, MPI_Comm_size, MPI_Allreduce, &
      MPI_COMM_WORLD, MPI_REAL8, MPI_INTEGER, MPI_SUM, MPI_MAX
   USE MOD_Precision, only: r8
   USE MOD_SPMD_Task, only: p_err, p_comm_glb, p_comm_worker, p_iam_glb, &
      p_iam_worker, p_np_glb, p_np_worker, p_np_io, p_address_master, &
      p_is_master, p_is_worker, p_is_io, p_is_writeback, p_itis_worker, p_address_worker
   USE MOD_Namelist, only: DEF_USE_LEVEE, DEF_TRACER_NUM, DEF_TRACER_NAMES, &
      DEF_TRACER_TYPES, DEF_TRACER_MRAT, DEF_TRACER_REF_RATIO, &
      DEF_TRACER_INIT_DELTA, DEF_TRACER_REACTIVE_DECAY_RATE
   USE MOD_WorkerPushData, only: build_worker_pushdata
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, npthout_local, npthlev_bif, &
      ucat_ucid, ucat_next, ucat_ups, wts_ups, lake_type, floodplain_curve, pth_upst_local, &
      pth_down_local, pth_down_ucid, pth_global_id, bif_incoming_pths, &
      bif_incoming_wts, push_next2ucat, push_ups2ucat, push_bif_dn2pth, push_bif_influx
   USE MOD_Grid_RiverLakeLevee, only: has_levee, levsto
   USE MOD_Grid_RiverLakeTimeVars, only: volwater_ucat, volwater_ucat_valid
   USE MOD_Tracer_RiverLake, only: river_lake_tracer_init, tracer_substep, &
      tracer_limiter_stats, trc_mass, trc_levsto, trc_solid, trc_levsto_solid, &
      trc_inp_buf, trc_flux_out, tracer_init_from_water
   USE MOD_Tracer_Defs, only: tracers
   IMPLICIT NONE
   integer :: rank, nranks, i, scenario, calls, iterations, peak, soft, peak_global
   integer :: failures, total_failures
   real(r8) :: initial(2), final(2), global_initial(2), global_final(2)
   real(r8) :: runoff(2), outlet(2), global_runoff(2), global_outlet(2)
   real(r8) :: main_flux(1), bif_flux(2,1), stage(1), residual(2)
   real(r8) :: mobile_before(2), solid_before(2), moved(2), global_moved(2)
   real(r8) :: uniform_ratio(2), water_after, incoming_water, uniform_error

   CALL MPI_Init(p_err)
   CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, p_err)
   CALL MPI_Comm_size(MPI_COMM_WORLD, nranks, p_err)
   IF (nranks < 2 .or. nranks > 4) ERROR STOP 'run with 2-4 MPI ranks'

   p_comm_glb = MPI_COMM_WORLD
   p_comm_worker = MPI_COMM_WORLD
   p_iam_glb = rank
   p_iam_worker = rank
   p_np_glb = nranks
   p_np_worker = nranks
   p_np_io = 0
   p_address_master = 0
   p_is_master = rank == 0
   p_is_worker = .true.
   p_is_io = .false.
   p_is_writeback = .false.
   allocate(p_itis_worker(0:nranks-1), p_address_worker(0:nranks-1))
   DO i = 0, nranks-1
      p_itis_worker(i) = i
      p_address_worker(i) = i
   ENDDO

   numucat = 1
   totalnumucat = nranks
   npthout_local = 1
   npthlev_bif = 2
   allocate(ucat_ucid(1), ucat_next(1), ucat_ups(1,1), wts_ups(1,1), lake_type(1))
   ucat_ucid = rank + 1
   ucat_next = rank + 2
   IF (rank == nranks-1) ucat_next = 0  ! sea outlet
   ucat_ups = rank                     ! first rank has no ordinary upstream
   wts_ups = 1._r8
   lake_type = 0
   CALL build_worker_pushdata(1, ucat_ucid, 1, ucat_next, push_next2ucat)
   CALL build_worker_pushdata(1, ucat_ucid, 1, ucat_ups, wts_ups, push_ups2ucat)

   allocate(pth_upst_local(1), pth_down_local(1), pth_down_ucid(1), &
      pth_global_id(1), bif_incoming_pths(1,1), bif_incoming_wts(1,1))
   pth_upst_local = 1
   pth_down_local = -1
   pth_down_ucid = modulo(rank+1,nranks)+1  ! BIF forms a cross-rank ring
   pth_global_id = rank+1
   bif_incoming_pths = modulo(rank-1+nranks,nranks)+1
   bif_incoming_wts = 1._r8
   CALL build_worker_pushdata(1, ucat_ucid, 1, pth_down_ucid, push_bif_dn2pth)
   CALL build_worker_pushdata(1, pth_global_id, 1, bif_incoming_pths, &
      bif_incoming_wts, push_bif_influx)

   DEF_USE_LEVEE = .true.
   allocate(has_levee(1), levsto(1), volwater_ucat(1))
   has_levee = .true.
   volwater_ucat_valid = .true.
   DEF_TRACER_NUM = 2
   DEF_TRACER_NAMES = 'A,B'
   DEF_TRACER_TYPES = 'solute,solute'
   DEF_TRACER_MRAT = '18,18'
   DEF_TRACER_REF_RATIO = '1,1'
   DEF_TRACER_INIT_DELTA = '0,0'
   DEF_TRACER_REACTIVE_DECAY_RATE = '0,0'
   CALL river_lake_tracer_init()

   failures = 0
   DO scenario = 1, 4
      main_flux = 0._r8
      bif_flux = 0._r8
      volwater_ucat = 1._r8
      levsto = 0.5_r8
      stage = 1._r8
      trc_inp_buf = 0._r8
      DO i = 1, 2
         trc_mass(i,:) = (0.5_r8 + 0.1_r8*rank)*i
         trc_levsto(i,:) = (0.1_r8 + 0.03_r8*rank)*(3-i)
      ENDDO
      SELECT CASE (scenario)
      CASE (1)  ! Reverse ordinary flow; positive visible and negative protected BIF.
         IF (rank < nranks-1) main_flux = -0.1_r8
         bif_flux(1,:) = 0.05_r8
         bif_flux(2,:) = -0.03_r8
      CASE (2)  ! Dry queued donor with two exits must enter MPI Jacobi fallback.
         IF (rank == 0) THEN
            volwater_ucat = 0._r8
            stage = 0._r8
            trc_mass(:,1) = [0.2_r8, 0.1_r8]
            main_flux = 0.1_r8
            bif_flux(1,:) = 0.1_r8
         ENDIF
      CASE (3)  ! Sub-dry-threshold protected pool exports its own signature.
         IF (rank == 0) THEN
            levsto = 1.e-9_r8
            trc_levsto(:,1) = [0.4e-9_r8, 0.7e-9_r8]
            bif_flux(2,:) = 1.e-9_r8
         ENDIF
      CASE (4)  ! Open budget: runoff at head, ordinary outlet at sea.
         main_flux = 0.1_r8
         IF (rank == 0) trc_inp_buf(:,1) = [0.07_r8, 0.11_r8]
      END SELECT

      initial = [sum(trc_mass(1,:))+sum(trc_levsto(1,:)), &
                 sum(trc_mass(2,:))+sum(trc_levsto(2,:))]
      runoff = 0._r8
      IF (rank == 0) runoff = trc_inp_buf(:,1)
      CALL MPI_Allreduce(initial, global_initial, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, p_err)
      CALL MPI_Allreduce(runoff, global_runoff, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, p_err)
      CALL tracer_limiter_stats(calls, iterations, peak, soft, reset=.true.)
      CALL tracer_substep(1._r8, [1._r8], [1], main_flux, [0._r8], &
         stage, [.true.], [real(r8)::], [integer::], [.false.], &
         do_bif=.true., bif_hflux_lev_in=bif_flux, npthout_local_in=1)
      final = [sum(trc_mass(1,:))+sum(trc_levsto(1,:)), &
               sum(trc_mass(2,:))+sum(trc_levsto(2,:))]
      outlet = 0._r8
      IF (rank == nranks-1) outlet = trc_flux_out(:,1)
      CALL MPI_Allreduce(final, global_final, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, p_err)
      CALL MPI_Allreduce(outlet, global_outlet, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, p_err)
      CALL tracer_limiter_stats(calls, iterations, peak, soft)
      CALL MPI_Allreduce(peak, peak_global, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, p_err)
      residual = global_final + global_outlet - global_initial - global_runoff
      IF (maxval(abs(residual)) > 2.e-12_r8) failures = failures+1
      IF (minval(trc_mass) < -1.e-12_r8 .or. minval(trc_levsto) < -1.e-12_r8) failures = failures+1
      IF (scenario == 1 .and. peak_global /= 0) failures = failures+1
      IF (scenario == 2 .and. peak_global < 1) failures = failures+1
      IF (rank == 0) WRITE(*,'(A,I0,A,2ES13.5,A,I0)') 'TRACER_DYNAMIC scenario=', &
         scenario, ' residual=', residual, ' peak_iter=', peak_global
      IF (rank == 0 .and. scenario == 4) WRITE(*,'(A,2ES13.5,A,2ES13.5,A,2ES13.5,A,2ES13.5)') &
         'TRACER_BUDGET initial=', global_initial, ' runoff=', global_runoff, &
         ' outlet=', global_outlet, ' final=', global_final
   ENDDO

   ! Compile-time coverage of the actual cold-start helper and transport,
   ! not a Python mirror: uniform composition must travel with its carrier.
   uniform_ratio = [0.002_r8, 0.00015_r8]
   tracers(1)%init_conc = uniform_ratio(1)
   tracers(2)%init_conc = uniform_ratio(2)
   allocate(floodplain_curve(1))
   floodplain_curve(1)%rivhgt = 2._r8
   floodplain_curve(1)%rivare = 10._r8
   floodplain_curve(1)%rivstomax = 20._r8
   DEF_USE_LEVEE = .false.
   has_levee = .false.
   stage = 1._r8
   volwater_ucat = 0._r8 ! Legacy cold checkpoint: wet stage, placeholder volume.
   volwater_ucat_valid = .true.
   CALL tracer_init_from_water(stage, [real(r8)::], [integer::], is_built_resv=[.false.])
   uniform_error = maxval(abs(trc_mass(:,1)-10._r8*uniform_ratio))
   IF (uniform_error > 1.e-14_r8) failures = failures+1
   IF (rank == 0) WRITE(*,'(A,ES13.5)') 'TRACER_UNIFORM cold_placeholder error=', uniform_error
   DEF_USE_LEVEE = .true.
   has_levee = .true.

   DO scenario = 1, 3
      stage = 1._r8
      volwater_ucat = 1._r8
      levsto = 0.5_r8
      trc_mass(:,1) = uniform_ratio
      trc_levsto(:,1) = 0.5_r8*uniform_ratio
      trc_inp_buf = 0._r8
      main_flux = 0._r8
      bif_flux = 0._r8
      water_after = 1._r8
      SELECT CASE (scenario)
      CASE (1)
         ! Closed cross-rank ring: gross export exceeds local storage, but
         ! equal incoming water keeps both visible and protected pools wet.
         bif_flux(1,:) = 1.4_r8
         bif_flux(2,:) = -0.6_r8
      CASE (2)
         main_flux = 0.2_r8
         incoming_water = 0.2_r8
         IF (rank == 0) incoming_water = 0._r8
         water_after = 1._r8 + incoming_water - main_flux(1)
      CASE (3)
         main_flux = -0.2_r8
         IF (rank == nranks-1) main_flux = 0._r8
         incoming_water = -0.2_r8
         IF (rank == 0) incoming_water = 0._r8
         water_after = 1._r8 + incoming_water - main_flux(1)
      END SELECT
      CALL tracer_substep(1._r8, [1._r8], [1], main_flux, [1._r8-water_after], &
         stage, [.true.], [real(r8)::], [integer::], [.false.], &
         do_bif=.true., bif_hflux_lev_in=bif_flux, npthout_local_in=1)
      uniform_error = maxval(abs(trc_mass(:,1)-uniform_ratio*water_after))
      uniform_error = max(uniform_error, maxval(abs(trc_levsto(:,1)-0.5_r8*uniform_ratio)))
      IF (uniform_error > 2.e-12_r8) failures = failures+1
      IF (rank == 0) WRITE(*,'(A,I0,A,ES13.5)') 'TRACER_UNIFORM transport=', scenario, ' error=', uniform_error
   ENDDO

   ! The same production BIF kernel with a finite-solubility row. Positive
   ! visible and negative protected fluxes both cross ranks, while solids
   ! stay attached to their source catchments.
   tracers(1)%max_dissolved_conc = 0.2_r8
   CALL river_lake_tracer_init()
   IF (.not. allocated(trc_solid)) failures = failures+1
   main_flux = 0._r8
   bif_flux(1,:) = 0.05_r8
   bif_flux(2,:) = -0.03_r8
   stage = 1._r8
   volwater_ucat = 1._r8
   levsto = 0.5_r8
   trc_inp_buf = 0._r8
   trc_mass = 0._r8
   trc_levsto = 0._r8
   trc_solid = 0._r8
   trc_levsto_solid = 0._r8
   trc_mass(1,1) = 0.06_r8 + 0.02_r8*rank
   trc_levsto(1,1) = 0.04_r8 + 0.01_r8*rank
   IF (rank == 0) THEN
      trc_mass(1,1) = 0.2_r8
      trc_levsto(1,1) = 0.1_r8
      trc_solid(1,1) = 1._r8
      trc_levsto_solid(1,1) = 0.5_r8
   ENDIF
   mobile_before = [trc_mass(1,1), trc_levsto(1,1)]
   solid_before = [trc_solid(1,1), trc_levsto_solid(1,1)]
   initial(1) = sum(mobile_before) + sum(solid_before)
   initial(2) = 0._r8
   CALL MPI_Allreduce(initial, global_initial, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, p_err)
   CALL tracer_substep(1._r8, [1._r8], [1], main_flux, [0._r8], &
      stage, [.true.], [real(r8)::], [integer::], [.false.], &
      do_bif=.true., bif_hflux_lev_in=bif_flux, npthout_local_in=1)
   final(1) = trc_mass(1,1) + trc_levsto(1,1) + trc_solid(1,1) + trc_levsto_solid(1,1)
   final(2) = trc_mass(2,1) + trc_levsto(2,1)
   CALL MPI_Allreduce(final, global_final, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, p_err)
   moved = [abs(trc_mass(1,1)-mobile_before(1)), abs(trc_levsto(1,1)-mobile_before(2))]
   CALL MPI_Allreduce(moved, global_moved, 2, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, p_err)
   IF (minval(global_moved) <= 1.e-8_r8) failures = failures+1
   IF (maxval(abs([trc_solid(1,1),trc_levsto_solid(1,1)]-solid_before)) > 1.e-12_r8) failures = failures+1
   IF (trc_mass(1,1) > 0.2_r8*volwater_ucat(1)+1.e-12_r8 .or. &
       trc_levsto(1,1) > 0.2_r8*levsto(1)+1.e-12_r8) failures = failures+1
   IF (maxval(abs(global_final-global_initial)) > 2.e-12_r8) failures = failures+1
   IF (rank == 0) WRITE(*,'(A,2ES13.5,A,ES13.5)') 'TRACER_FINITE_BIF moved=', &
      global_moved, ' residual=', global_final(1)-global_initial(1)

   CALL MPI_Allreduce(failures, total_failures, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, p_err)
   IF (rank == 0) WRITE(*,'(A,I0,A,I0)') 'TRACER_DYNAMIC ranks=', nranks, ' failures=', total_failures
   CALL MPI_Finalize(p_err)
   IF (total_failures /= 0) ERROR STOP 1
END PROGRAM river_bif_tracer_harness
