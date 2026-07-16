PROGRAM river_bif_physics_harness
   USE mpi, only: MPI_Init, MPI_Finalize, MPI_Comm_rank, MPI_Comm_size, &
      MPI_Allreduce, MPI_COMM_WORLD_F => MPI_COMM_WORLD, &
      MPI_INTEGER_F => MPI_INTEGER, MPI_REAL8_F => MPI_REAL8, &
      MPI_SUM_F => MPI_SUM, MPI_MAX_F => MPI_MAX
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   USE MOD_Precision, only: r8
   USE MOD_SPMD_Task, only: p_err, p_comm_glb, p_comm_worker, p_iam_glb, &
      p_iam_worker, p_np_glb, p_np_worker, p_np_io, p_address_master, &
      p_is_master, p_is_worker, p_is_io, p_is_writeback, p_itis_worker, &
      p_address_worker
   USE MOD_Namelist, only: DEF_USE_LEVEE, DEF_USE_BIFURCATION
   USE MOD_WorkerPushData, only: build_worker_pushdata
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir, only: ucat2resv, reservoir_final
   USE MOD_Grid_RiverLakeLevee, only: levee_init, levee_fldstg, &
      levee_visible_volume_from_stage, levee_apply_protected_flux, &
      levee_repartition_storage, levee_final, levsto, levdph
   USE MOD_Grid_RiverLakeBifurcation, only: bifurcation_init, &
      bifurcation_calc, bifurcation_final, bifurcation_invalidate_static_dn, &
      bif_hflux_sum, bif_hflux_lev, bif_lev_hflux_sum
   IMPLICIT NONE

   integer :: rank, nranks, local_failures, global_failures
   integer :: ncell, npath
   real(r8), allocatable :: visible_volume(:), total_initial(:)
   real(r8), allocatable, target :: stage(:), stage_previous(:)
   logical, allocatable :: is_reservoir(:), active(:)
   integer, allocatable :: river_system(:)
   real(r8), allocatable :: dt_all(:), normal_outgoing(:), empty_reservoir(:)

   CALL MPI_Init(p_err)
   CALL MPI_Comm_rank(MPI_COMM_WORLD_F, rank, p_err)
   CALL MPI_Comm_size(MPI_COMM_WORLD_F, nranks, p_err)
   CALL configure_spmd(rank, nranks)

   local_failures = 0
   DEF_USE_LEVEE = .true.
   DEF_USE_BIFURCATION = .true.

   CALL configure_network(rank, nranks, ncell, npath)
   CALL levee_init()
   CALL bifurcation_init()

   allocate (visible_volume(ncell), total_initial(ncell), stage(ncell), stage_previous(ncell))
   allocate (is_reservoir(ncell), active(ncell), river_system(ncell))
   allocate (dt_all(1), normal_outgoing(ncell), empty_reservoir(0))
   is_reservoir = .false.
   active = .true.
   river_system = 1
   dt_all = 60._r8
   normal_outgoing = 0._r8

   CALL test_levee_curve_inverse_and_boundaries(local_failures)
   CALL test_channel_previous_depth_geometric_mean(local_failures)

   CALL bifurcation_final()
   DEF_USE_LEVEE = .false.
   CALL bifurcation_init()
   CALL test_nonlevee_multilayer_semi_implicit(local_failures)

   CALL bifurcation_final()
   DEF_USE_LEVEE = .true.
   CALL bifurcation_init()
   CALL test_tiny_residual_zero_donor(local_failures)

   CALL bifurcation_final()
   CALL bifurcation_init()
   CALL test_multilayer_conservation_and_wet_dry(local_failures)

   CALL MPI_Allreduce(local_failures, global_failures, 1, MPI_INTEGER_F, MPI_SUM_F, &
      MPI_COMM_WORLD_F, p_err)
   IF (rank == 0) THEN
      IF (global_failures == 0) THEN
         write(*,'(A,I0,A)') 'river BIF physics harness: PASS (', nranks, ' ranks)'
      ELSE
         write(*,'(A,I0)') 'river BIF physics harness: FAIL count=', global_failures
      ENDIF
   ENDIF

   CALL bifurcation_final()
   CALL levee_final()
   CALL reservoir_final()
   CALL riverlake_network_final()
   CALL release_spmd()
   CALL MPI_Finalize(p_err)
   IF (global_failures /= 0) ERROR STOP 1

CONTAINS

   SUBROUTINE configure_spmd(rank_in, nranks_in)
      integer, intent(in) :: rank_in, nranks_in
      integer :: i

      p_comm_glb = MPI_COMM_WORLD_F
      p_comm_worker = MPI_COMM_WORLD_F
      p_iam_glb = rank_in
      p_iam_worker = rank_in
      p_np_glb = nranks_in
      p_np_worker = nranks_in
      p_np_io = 0
      p_address_master = 0
      p_is_master = rank_in == 0
      p_is_worker = .true.
      p_is_io = .false.
      p_is_writeback = .false.

      allocate (p_itis_worker(0:nranks_in-1))
      allocate (p_address_worker(0:nranks_in-1))
      DO i = 0, nranks_in-1
         p_itis_worker(i) = i
         p_address_worker(i) = i
      ENDDO
   END SUBROUTINE configure_spmd


   SUBROUTINE release_spmd()
      IF (allocated(p_itis_worker)) deallocate (p_itis_worker)
      IF (allocated(p_address_worker)) deallocate (p_address_worker)
   END SUBROUTINE release_spmd


   SUBROUTINE configure_network(rank_in, nranks_in, ncell_out, npath_out)
      integer, intent(in) :: rank_in, nranks_in
      integer, intent(out) :: ncell_out, npath_out
      integer :: i, gid, down_gid, incoming_gid

      IF (nranks_in == 1) THEN
         ncell_out = 2
         npath_out = 2
         totalnumucat = 2
         totalnpthout = 2
      ELSE
         ncell_out = 1
         npath_out = 1
         totalnumucat = nranks_in
         totalnpthout = nranks_in
      ENDIF
      numucat = ncell_out
      npthout_local = npath_out
      npthlev_bif = 2

      allocate (ucat_ucid(ncell_out), x_ucat(ncell_out), y_ucat(ncell_out))
      allocate (topo_rivelv(ncell_out), topo_rivhgt(ncell_out), topo_rivlen(ncell_out))
      allocate (topo_rivwth(ncell_out), topo_rivare(ncell_out), topo_rivstomax(ncell_out))
      allocate (topo_area(ncell_out), topo_fldhgt(4,ncell_out))
      allocate (lake_type(ncell_out), levee_frc_data(ncell_out), levee_hgt_data(ncell_out))
      allocate (floodplain_curve(ncell_out), ucat2resv(ncell_out))

      DO i = 1, ncell_out
         IF (nranks_in == 1) THEN
            gid = i
         ELSE
            gid = rank_in + 1
         ENDIF
         ucat_ucid(i) = gid
         x_ucat(i) = gid
         y_ucat(i) = 1
         topo_rivelv(i) = 0._r8
         topo_rivhgt(i) = 2._r8
         topo_rivlen(i) = 100._r8
         topo_rivwth(i) = 10._r8
         topo_rivare(i) = 1000._r8
         topo_rivstomax(i) = 2000._r8
         topo_area(i) = 4000._r8
         topo_fldhgt(:,i) = [0.5_r8, 1._r8, 1.5_r8, 2._r8]
         lake_type(i) = 0
         levee_frc_data(i) = 0.5_r8
         levee_hgt_data(i) = 1.5_r8
         ucat2resv(i) = 0
         CALL configure_curve(i)
      ENDDO

      allocate (pth_upst_local(npath_out), pth_down_local(npath_out))
      allocate (pth_down_ucid(npath_out), pth_global_id(npath_out))
      allocate (pth_dst(npath_out), pth_elv(2,npath_out), pth_wth(2,npath_out))
      allocate (pth_man(2), bif_incoming_pths(1,ncell_out), bif_incoming_wts(1,ncell_out))

      IF (nranks_in == 1) THEN
         pth_upst_local = [1, 2]
         pth_down_local = [2, 1]
         pth_down_ucid = [2, 1]
         pth_global_id = [1, 2]
         bif_incoming_pths(1,:) = [2, 1]
      ELSE
         gid = rank_in + 1
         down_gid = modulo(rank_in + 1, nranks_in) + 1
         incoming_gid = modulo(rank_in - 1 + nranks_in, nranks_in) + 1
         pth_upst_local(1) = 1
         pth_down_local(1) = -1
         pth_down_ucid(1) = down_gid
         pth_global_id(1) = gid
         bif_incoming_pths(1,1) = incoming_gid
      ENDIF
      pth_dst = 1000._r8
      pth_elv(1,:) = 0._r8
      pth_elv(2,:) = 2._r8
      pth_wth(1,:) = 10._r8
      pth_wth(2,:) = 5._r8
      pth_man = [0.03_r8, 0.05_r8]
      bif_incoming_wts = 1._r8

      CALL build_worker_pushdata(numucat, ucat_ucid, npthout_local, &
         pth_down_ucid, push_bif_dn2pth)
      CALL build_worker_pushdata(npthout_local, pth_global_id, numucat, &
         bif_incoming_pths, bif_incoming_wts, push_bif_influx)
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
      floodplain_curve(i)%flparea(0) = 0._r8
      floodplain_curve(i)%flparea(1:) = 1000._r8
      floodplain_curve(i)%flpaccare(0) = 0._r8
      floodplain_curve(i)%flpstomax(0) = 0._r8
      DO level = 1, 4
         floodplain_curve(i)%flpaccare(level) = &
            floodplain_curve(i)%flpaccare(level-1) + floodplain_curve(i)%flparea(level)
         floodplain_curve(i)%flpstomax(level) = floodplain_curve(i)%flpstomax(level-1) &
            + 0.5_r8 * (floodplain_curve(i)%flpaccare(level) &
            + floodplain_curve(i)%flpaccare(level-1)) &
            * (floodplain_curve(i)%flphgt(level) - floodplain_curve(i)%flphgt(level-1))
      ENDDO
   END SUBROUTINE configure_curve


   SUBROUTINE test_levee_curve_inverse_and_boundaries(failures)
      integer, intent(inout) :: failures
      real(r8), parameter :: samples(10) = [0._r8, 1000._r8, 2000._r8, &
         3000._r8, 4000._r8, 5250._r8, 5500._r8, 5625._r8, 5750._r8, 6500._r8]
      real(r8), parameter :: boundaries(4) = [2000._r8, 4000._r8, 5500._r8, 5750._r8]
      real(r8) :: volume, depth, protected, protected_depth, flood_fraction
      real(r8) :: visible, depth_left, depth_right, protected_left, protected_right
      real(r8) :: dummy_depth, eps, tol
      integer :: i

      DO i = 1, size(samples)
         volume = samples(i)
         CALL levee_fldstg(1, volume, depth, protected, protected_depth, flood_fraction)
         visible = levee_visible_volume_from_stage(1, depth, protected)
         tol = 2.e-10_r8 * max(1._r8, volume)
         IF (abs(visible + protected - volume) > tol) failures = failures + 1
         IF (.not. ieee_is_finite(depth) .or. .not. ieee_is_finite(protected_depth)) failures = failures + 1
         IF (depth < 0._r8 .or. protected < 0._r8 .or. protected_depth < 0._r8) failures = failures + 1
         IF (flood_fraction < 0._r8 .or. flood_fraction > 1._r8) failures = failures + 1
      ENDDO

      DO i = 1, size(boundaries)
         eps = 1.e-7_r8 * max(1._r8, boundaries(i))
         CALL levee_fldstg(1, boundaries(i)-eps, depth_left, protected_left, &
            dummy_depth, flood_fraction)
         CALL levee_fldstg(1, boundaries(i)+eps, depth_right, protected_right, &
            dummy_depth, flood_fraction)
         IF (abs(depth_right-depth_left) > 2.e-5_r8) failures = failures + 1
         IF (abs(protected_right-protected_left) > 4._r8*eps) failures = failures + 1
      ENDDO
   END SUBROUTINE test_levee_curve_inverse_and_boundaries


   SUBROUTINE test_channel_previous_depth_geometric_mean(failures)
      integer, intent(inout) :: failures
      real(r8), parameter :: previous_scale = 0.25_r8
      real(r8) :: explicit_peak, geometric_peak, fallback_peak
      real(r8) :: local_peak, expected_ratio, measured_ratio
      integer :: i

      ! With zero initial momentum and no active storage limiter, the first
      ! channel-layer flux is proportional to the face depth used by the
      ! local-inertial pressure-gradient update.  This makes the actual BIF
      ! flux ratio a direct executable check of CaMa's
      ! sqrt(current_depth * previous_depth) channel-depth rule.
      DO i = 1, ncell
         IF (modulo(ucat_ucid(i), 2) == 1) THEN
            stage(i) = 1._r8
         ELSE
            stage(i) = 0.9_r8
         ENDIF
         stage_previous(i) = stage(i)
         visible_volume(i) = floodplain_curve(i)%volume(stage(i))
         levsto(i) = 0._r8
         levdph(i) = 0._r8
      ENDDO
      CALL bifurcation_invalidate_static_dn()
      CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., &
         empty_reservoir, is_reservoir, dt_all, river_system, active, normal_outgoing)
      local_peak = maxval(abs(bif_hflux_lev(1,:)))
      CALL MPI_Allreduce(local_peak, explicit_peak, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)

      CALL bifurcation_final()
      CALL bifurcation_init()
      stage_previous = previous_scale * stage
      CALL bifurcation_invalidate_static_dn()
      CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., &
         empty_reservoir, is_reservoir, dt_all, river_system, active, normal_outgoing)
      local_peak = maxval(abs(bif_hflux_lev(1,:)))
      CALL MPI_Allreduce(local_peak, geometric_peak, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)

      expected_ratio = sqrt(previous_scale)
      IF (explicit_peak <= 0._r8) THEN
         failures = failures + 1
      ELSE
         measured_ratio = geometric_peak / explicit_peak
         IF (abs(measured_ratio-expected_ratio) > 2.e-11_r8) failures = failures + 1
      ENDIF

      ! CaMa's dry-previous-depth fallback is the current face depth.
      CALL bifurcation_final()
      CALL bifurcation_init()
      stage_previous = 0._r8
      CALL bifurcation_invalidate_static_dn()
      CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., &
         empty_reservoir, is_reservoir, dt_all, river_system, active, normal_outgoing)
      local_peak = maxval(abs(bif_hflux_lev(1,:)))
      CALL MPI_Allreduce(local_peak, fallback_peak, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)
      IF (abs(fallback_peak-explicit_peak) > 2.e-11_r8 * max(1._r8, explicit_peak)) &
         failures = failures + 1
   END SUBROUTINE test_channel_previous_depth_geometric_mean


   SUBROUTINE test_nonlevee_multilayer_semi_implicit(failures)
      integer, intent(inout) :: failures
      real(r8) :: baseline_channel, baseline_overbank
      real(r8) :: geometric_channel, geometric_overbank
      real(r8) :: dry_channel, dry_overbank
      real(r8) :: local_channel, local_overbank
      real(r8) :: measured_ratio, expected_ratio
      integer :: i

      ! CaMa selects a different pathway-depth rule when levees are disabled:
      ! every layer is semi-implicit, and a 0.01 m previous-depth floor avoids
      ! a zero-conductance cold start.  Use water surfaces that make the
      ! current face depths exactly 3 m (channel) and 1 m (overbank).
      DO i = 1, ncell
         IF (modulo(ucat_ucid(i), 2) == 1) THEN
            stage(i) = 3._r8
            stage_previous(i) = 3._r8
         ELSE
            stage(i) = 2.9_r8
            stage_previous(i) = 2.9_r8
         ENDIF
         visible_volume(i) = 1.e8_r8
         levsto(i) = 0._r8
         levdph(i) = 0._r8
      ENDDO
      CALL bifurcation_invalidate_static_dn()
      CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., &
         empty_reservoir, is_reservoir, dt_all, river_system, active, normal_outgoing)
      local_channel = maxval(abs(bif_hflux_lev(1,:)))
      local_overbank = maxval(abs(bif_hflux_lev(2,:)))
      CALL MPI_Allreduce(local_channel, baseline_channel, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)
      CALL MPI_Allreduce(local_overbank, baseline_overbank, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)

      CALL bifurcation_final()
      CALL bifurcation_init()
      DO i = 1, ncell
         IF (modulo(ucat_ucid(i), 2) == 1) THEN
            stage_previous(i) = 2.25_r8
         ELSE
            stage_previous(i) = 2.15_r8
         ENDIF
      ENDDO
      CALL bifurcation_invalidate_static_dn()
      CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., &
         empty_reservoir, is_reservoir, dt_all, river_system, active, normal_outgoing)
      local_channel = maxval(abs(bif_hflux_lev(1,:)))
      local_overbank = maxval(abs(bif_hflux_lev(2,:)))
      CALL MPI_Allreduce(local_channel, geometric_channel, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)
      CALL MPI_Allreduce(local_overbank, geometric_overbank, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)

      IF (baseline_channel <= 0._r8 .or. baseline_overbank <= 0._r8) THEN
         failures = failures + 1
      ELSE
         measured_ratio = geometric_channel / baseline_channel
         expected_ratio = sqrt(2.25_r8 / 3._r8)
         IF (abs(measured_ratio-expected_ratio) > 2.e-11_r8) failures = failures + 1
         measured_ratio = geometric_overbank / baseline_overbank
         expected_ratio = sqrt(0.25_r8 / 1._r8)
         IF (abs(measured_ratio-expected_ratio) > 2.e-11_r8) failures = failures + 1
      ENDIF

      CALL bifurcation_final()
      CALL bifurcation_init()
      stage_previous = 0._r8
      CALL bifurcation_invalidate_static_dn()
      CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., &
         empty_reservoir, is_reservoir, dt_all, river_system, active, normal_outgoing)
      local_channel = maxval(abs(bif_hflux_lev(1,:)))
      local_overbank = maxval(abs(bif_hflux_lev(2,:)))
      CALL MPI_Allreduce(local_channel, dry_channel, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)
      CALL MPI_Allreduce(local_overbank, dry_overbank, 1, MPI_REAL8_F, MPI_MAX_F, &
         MPI_COMM_WORLD_F, p_err)

      IF (baseline_channel > 0._r8) THEN
         measured_ratio = dry_channel / baseline_channel
         expected_ratio = sqrt(0.01_r8 / 3._r8)
         IF (abs(measured_ratio-expected_ratio) > 2.e-11_r8) failures = failures + 1
      ENDIF
      IF (baseline_overbank > 0._r8) THEN
         measured_ratio = dry_overbank / baseline_overbank
         expected_ratio = sqrt(0.01_r8 / 1._r8)
         IF (abs(measured_ratio-expected_ratio) > 2.e-11_r8) failures = failures + 1
      ENDIF
   END SUBROUTINE test_nonlevee_multilayer_semi_implicit


   SUBROUTINE test_tiny_residual_zero_donor(failures)
      integer, intent(inout) :: failures
      real(r8), parameter :: dry_stage = 1.e-5_r8
      real(r8), parameter :: stage_delta = 8.5e-7_r8
      real(r8) :: local_peak, global_peak
      integer :: i

      DO i = 1, ncell
         IF (modulo(ucat_ucid(i), 2) == 1) THEN
            stage(i) = dry_stage + stage_delta
         ELSE
            stage(i) = dry_stage
         ENDIF
         stage_previous(i) = stage(i)
         visible_volume(i) = floodplain_curve(i)%volume(stage(i))
         levsto(i) = 0._r8
         levdph(i) = 0._r8
      ENDDO
      CALL bifurcation_invalidate_static_dn()
      CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., empty_reservoir, &
         is_reservoir, dt_all, river_system, active, normal_outgoing)

      local_peak = 0._r8
      IF (npath > 0) local_peak = maxval(abs(bif_hflux_lev(1,:)))
      CALL MPI_Allreduce(local_peak, global_peak, 1, MPI_REAL8_F, MPI_MAX_F, MPI_COMM_WORLD_F, p_err)
      IF (global_peak <= 0._r8 .or. global_peak >= 1.e-9_r8) failures = failures + 1

      stage = dry_stage
      stage_previous = stage
      visible_volume = 0._r8
      levsto = 0._r8
      levdph = 0._r8
      CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., empty_reservoir, &
         is_reservoir, dt_all, river_system, active, normal_outgoing)

      local_peak = 0._r8
      IF (npath > 0) local_peak = maxval(abs(bif_hflux_lev))
      CALL MPI_Allreduce(local_peak, global_peak, 1, MPI_REAL8_F, MPI_MAX_F, MPI_COMM_WORLD_F, p_err)
      IF (global_peak > 1.e-25_r8) failures = failures + 1
      IF (any(visible_volume < 0._r8) .or. any(levsto < 0._r8)) failures = failures + 1
   END SUBROUTINE test_tiny_residual_zero_donor


   SUBROUTINE test_multilayer_conservation_and_wet_dry(failures)
      integer, intent(inout) :: failures
      integer :: i, step
      integer :: local_channel_positive, local_channel_negative
      integer :: global_channel_positive, global_channel_negative
      integer :: local_overbank_positive, local_overbank_negative
      integer :: global_overbank_positive, global_overbank_negative
      integer :: local_wet, local_dry, global_wet, global_dry
      real(r8) :: flood_fraction, clipped, visible_flux
      real(r8) :: local_flux_sum, global_flux_sum, local_layer_sum, global_layer_sum
      real(r8) :: local_total, global_total, initial_global_total, tolerance

      DO i = 1, ncell
         IF (modulo(ucat_ucid(i), 2) == 1) THEN
            total_initial(i) = 6500._r8
         ELSE
            total_initial(i) = 10._r8
         ENDIF
         CALL levee_fldstg(i, total_initial(i), stage(i), levsto(i), levdph(i), flood_fraction)
         stage_previous(i) = stage(i)
         visible_volume(i) = total_initial(i) - levsto(i)
      ENDDO

      local_total = sum(visible_volume + levsto)
      CALL MPI_Allreduce(local_total, initial_global_total, 1, MPI_REAL8_F, MPI_SUM_F, &
         MPI_COMM_WORLD_F, p_err)
      CALL bifurcation_invalidate_static_dn()

      DO step = 1, 8
         CALL bifurcation_calc(stage, stage_previous, visible_volume, .true., empty_reservoir, &
            is_reservoir, dt_all, river_system, active, normal_outgoing)

         local_flux_sum = sum(bif_hflux_sum)
         local_layer_sum = sum(bif_lev_hflux_sum)
         CALL MPI_Allreduce(local_flux_sum, global_flux_sum, 1, MPI_REAL8_F, MPI_SUM_F, &
            MPI_COMM_WORLD_F, p_err)
         CALL MPI_Allreduce(local_layer_sum, global_layer_sum, 1, MPI_REAL8_F, MPI_SUM_F, &
            MPI_COMM_WORLD_F, p_err)
         tolerance = 2.e-12_r8 * max(1._r8, sum(abs(bif_hflux_sum)))
         IF (abs(global_flux_sum) > tolerance) failures = failures + 1
         IF (abs(global_layer_sum) > tolerance) failures = failures + 1

         IF (step == 1) THEN
            local_channel_positive = count(bif_hflux_lev(1,:) > 0._r8)
            local_channel_negative = count(bif_hflux_lev(1,:) < 0._r8)
            local_overbank_positive = count(bif_hflux_lev(2,:) > 0._r8)
            local_overbank_negative = count(bif_hflux_lev(2,:) < 0._r8)
            CALL MPI_Allreduce(local_channel_positive, global_channel_positive, 1, &
               MPI_INTEGER_F, MPI_SUM_F, MPI_COMM_WORLD_F, p_err)
            CALL MPI_Allreduce(local_channel_negative, global_channel_negative, 1, &
               MPI_INTEGER_F, MPI_SUM_F, MPI_COMM_WORLD_F, p_err)
            CALL MPI_Allreduce(local_overbank_positive, global_overbank_positive, 1, &
               MPI_INTEGER_F, MPI_SUM_F, MPI_COMM_WORLD_F, p_err)
            CALL MPI_Allreduce(local_overbank_negative, global_overbank_negative, 1, &
               MPI_INTEGER_F, MPI_SUM_F, MPI_COMM_WORLD_F, p_err)
            IF (global_channel_positive == 0 .or. global_channel_negative == 0) &
               failures = failures + 1
            IF (global_overbank_positive == 0 .or. global_overbank_negative == 0) &
               failures = failures + 1
         ENDIF

         ! Match the production substep ordering: preserve the current depth
         ! after flux evaluation but before the stage update below.
         stage_previous = stage
         DO i = 1, ncell
            visible_flux = bif_hflux_sum(i) - bif_lev_hflux_sum(i)
            visible_volume(i) = visible_volume(i) - visible_flux * dt_all(1)
            CALL levee_apply_protected_flux(i, bif_lev_hflux_sum(i), dt_all(1), clipped)
            IF (clipped > 1.e-12_r8) failures = failures + 1
            CALL levee_repartition_storage(i, visible_volume(i), stage(i), flood_fraction)
         ENDDO

         IF (any(.not. ieee_is_finite(visible_volume)) .or. &
             any(.not. ieee_is_finite(levsto)) .or. any(.not. ieee_is_finite(stage))) &
            failures = failures + 1
         IF (any(visible_volume < -1.e-12_r8) .or. any(levsto < -1.e-12_r8) .or. &
             any(stage < -1.e-12_r8)) failures = failures + 1

         local_total = sum(visible_volume + levsto)
         CALL MPI_Allreduce(local_total, global_total, 1, MPI_REAL8_F, MPI_SUM_F, &
            MPI_COMM_WORLD_F, p_err)
         tolerance = 2.e-10_r8 * max(1._r8, initial_global_total)
         IF (abs(global_total-initial_global_total) > tolerance) failures = failures + 1
      ENDDO

      local_wet = 0
      local_dry = 0
      DO i = 1, ncell
         IF (modulo(ucat_ucid(i), 2) == 0 .and. &
             visible_volume(i)+levsto(i) > total_initial(i)) local_wet = local_wet + 1
         IF (modulo(ucat_ucid(i), 2) == 1 .and. &
             visible_volume(i)+levsto(i) < total_initial(i)) local_dry = local_dry + 1
      ENDDO
      CALL MPI_Allreduce(local_wet, global_wet, 1, MPI_INTEGER_F, MPI_SUM_F, MPI_COMM_WORLD_F, p_err)
      CALL MPI_Allreduce(local_dry, global_dry, 1, MPI_INTEGER_F, MPI_SUM_F, MPI_COMM_WORLD_F, p_err)
      IF (global_wet == 0 .or. global_dry == 0) failures = failures + 1
   END SUBROUTINE test_multilayer_conservation_and_wet_dry

END PROGRAM river_bif_physics_harness
