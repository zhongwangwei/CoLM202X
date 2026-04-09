#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeBifurcation
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Bifurcation (multi-channel flow) module for grid-based river-lake routing.
!   Computes water exchange through bifurcation pathways using a Riemann-solver-
!   based shallow water formulation consistent with the main channel solver in
!   MOD_Grid_RiverLakeFlow.
!
! Created by CoLM team, April 2026
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist, only: DEF_USE_LEVEE
   USE MOD_SPMD_Task
   USE MOD_Grid_Reservoir, only: ucat2resv
   USE MOD_Grid_RiverLakeLevee, only: has_levee, levdph, levsto, levee_hgt, levee_topsto, &
      volwater_ucat, volwater_ucat_valid
   USE MOD_Grid_RiverLakeNetwork, only: floodplain_curve
   IMPLICIT NONE

   real(r8), parameter :: BIFMIN = 1.e-5_r8

   ! ----- State variables -----
   real(r8), allocatable :: pth_veloc     (:,:)  ! velocity (npthlev, npthout_local) [m/s]
   real(r8), allocatable :: pth_momen     (:,:)  ! momentum (npthlev, npthout_local) [m^2/s]
   real(r8), allocatable :: bif_hflux_lev (:,:)  ! effective volume flux per pathway layer [m^3/s]
   real(r8), allocatable :: bif_hflux_sum (:)    ! net volume flux per ucat [m^3/s]

   PUBLIC :: bifurcation_init
   PUBLIC :: bifurcation_calc
   PUBLIC :: read_bifurcation_restart
   PUBLIC :: write_bifurcation_restart
   PUBLIC :: bifurcation_final

CONTAINS

   ! =========================================================================
   SUBROUTINE bifurcation_init ()
   ! =========================================================================
   !
   ! Allocate and initialize bifurcation state arrays.
   !
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: numucat, npthout_local, npthlev_bif

   IMPLICIT NONE

      IF (.not. p_is_worker) RETURN

      allocate (pth_veloc     (npthlev_bif, npthout_local))
      allocate (pth_momen     (npthlev_bif, npthout_local))
      allocate (bif_hflux_lev (npthlev_bif, npthout_local))
      allocate (bif_hflux_sum (numucat))

      pth_veloc     (:,:) = 0._r8
      pth_momen     (:,:) = 0._r8
      bif_hflux_lev (:,:) = 0._r8
      bif_hflux_sum (:)   = 0._r8

   END SUBROUTINE bifurcation_init


   ! =========================================================================
   SUBROUTINE bifurcation_calc (wdsrf_ucat, veloc_riv, volresv, is_built_resv, dt_all, irivsys, ucatfilter)
   ! =========================================================================
   !
   ! Compute bifurcation fluxes for one sub-timestep using a two-rarefaction
   ! Riemann solver. The net volume flux per unit catchment is accumulated
   ! into bif_hflux_sum(:).
   !
   ! =========================================================================

   USE MOD_Const_Physical, only: grav
   USE MOD_WorkerPushData
   USE MOD_Grid_RiverLakeNetwork, only: &
      numucat, npthout_local, npthlev_bif, &
      pth_upst_local, pth_down_local, &
      pth_dst, pth_elv, pth_wth, pth_man, &
      pth_global_id, &
      push_bif_dn2pth, push_bif_influx, &
      topo_rivelv

   IMPLICIT NONE

   real(r8), intent(in) :: wdsrf_ucat (:)   ! water depth above riverbed [m]
   real(r8), intent(in) :: veloc_riv  (:)   ! river velocity [m/s]
   real(r8), intent(in) :: volresv    (:)   ! reservoir water volume [m^3]
   logical,  intent(in) :: is_built_resv (:) ! reservoir-built mask
   real(r8), intent(in) :: dt_all     (:)   ! timestep per ucat [s]
   integer,  intent(in) :: irivsys    (:)   ! river system id per ucat
   logical,  intent(in) :: ucatfilter (:)   ! active ucat filter

   ! Local variables
   real(r8), allocatable :: wdsrf_dn_pth     (:)  ! downstream water depth at pathway
   real(r8), allocatable :: rivelv_dn_pth    (:)  ! downstream riverbed elevation at pathway
   real(r8), allocatable :: levdph_dn_pth    (:)  ! downstream protected-side depth at pathway
   real(r8), allocatable :: has_levee_dn_pth (:)  ! downstream levee flag at pathway [0/1]
   real(r8), allocatable :: storage_dn_pth   (:)  ! downstream available storage at pathway [m^3]
   real(r8), allocatable :: storage_ucat     (:)  ! local available storage per ucat [m^3]
   real(r8), allocatable :: pth_hflux_total (:) ! total volume flux per pathway [m^3/s]
   real(r8), allocatable :: bif_influx (:)      ! incoming flux from remote pathways [m^3/s]
   real(r8), allocatable :: has_levee_r8   (:)  ! local levee flag in real form for push

   real(r8) :: rivelv_up, rivelv_dn, wdsrf_up, wdsrf_dn, wdsrf_up_eff, wdsrf_dn_eff
   real(r8) :: bedelv_pth, height_up, height_dn
   real(r8) :: v_up, v_dn
   real(r8) :: veloct_fc, height_fc
   real(r8) :: vwave_up, vwave_dn
   real(r8) :: hflux_up, hflux_dn, mflux_up, mflux_dn
   real(r8) :: hflux_lev, mflux_lev  ! flux for this layer
   real(r8) :: width_pth, pth_are
   real(r8) :: friction
   real(r8) :: dt, storage_up, storage_dn, storage_ref, rate
   integer  :: ipth, ilev, i_up, i_dn

      IF (.not. p_is_worker) RETURN

      ! Reset net flux accumulator
      bif_hflux_sum(:) = 0._r8
      bif_hflux_lev(:,:) = 0._r8

      IF (npthout_local <= 0 .and. numucat <= 0) THEN
         ! Still must participate in MPI push calls below
         allocate (wdsrf_dn_pth  (0))
         allocate (rivelv_dn_pth (0))
         allocate (levdph_dn_pth (0))
         allocate (has_levee_dn_pth (0))
         allocate (pth_hflux_total (0))
         allocate (bif_influx (0))
         allocate (has_levee_r8 (0))
         allocate (storage_ucat (0))
         allocate (storage_dn_pth (0))
      ELSE
         allocate (wdsrf_dn_pth  (npthout_local))
         allocate (rivelv_dn_pth (npthout_local))
         allocate (levdph_dn_pth (npthout_local))
         allocate (has_levee_dn_pth (npthout_local))
         allocate (storage_dn_pth  (npthout_local))
         allocate (pth_hflux_total (npthout_local))
         allocate (bif_influx (numucat))
         allocate (has_levee_r8 (numucat))
         allocate (storage_ucat (numucat))
      ENDIF

      wdsrf_dn_pth     (:) = 0._r8
      rivelv_dn_pth    (:) = 0._r8
      levdph_dn_pth    (:) = 0._r8
      has_levee_dn_pth (:) = 0._r8
      storage_dn_pth   (:) = 0._r8
      pth_hflux_total(:) = 0._r8
      storage_ucat    (:) = 0._r8

      IF (DEF_USE_LEVEE .and. allocated(has_levee_r8)) THEN
         has_levee_r8(:) = 0._r8
         IF (allocated(has_levee)) THEN
            WHERE (has_levee)
               has_levee_r8 = 1._r8
            END WHERE
         ENDIF
      ENDIF

      IF (numucat > 0) THEN
         DO i_up = 1, numucat
            storage_ucat(i_up) = available_storage_ucat(i_up, wdsrf_ucat(i_up), volresv, is_built_resv)
         ENDDO
      ENDIF

      ! ----- Step 2: Get downstream cell state via push objects -----
      ! Use -9999 as fillvalue to mark pathways with no valid downstream cell
      ! (domain boundary or unresolved remote). These are skipped in Step 3.
      CALL worker_push_data (push_bif_dn2pth, wdsrf_ucat,  wdsrf_dn_pth,  fillvalue = -9999._r8)
      CALL worker_push_data (push_bif_dn2pth, topo_rivelv, rivelv_dn_pth, fillvalue = -9999._r8)
      IF (DEF_USE_LEVEE) THEN
         CALL worker_push_data (push_bif_dn2pth, levdph,        levdph_dn_pth,    fillvalue = 0._r8)
         CALL worker_push_data (push_bif_dn2pth, has_levee_r8,  has_levee_dn_pth, fillvalue = 0._r8)
      ENDIF
      CALL worker_push_data (push_bif_dn2pth, storage_ucat, storage_dn_pth, fillvalue = -9999._r8)

      ! ----- Step 3: Riemann solver for each pathway / layer -----
      DO ipth = 1, npthout_local

         i_up = pth_upst_local(ipth)

         ! Skip if upstream ucat is not active
         IF (i_up < 1 .or. i_up > numucat) CYCLE
         IF (.not. ucatfilter(i_up)) CYCLE

         ! Skip if downstream cell is outside domain (CaMa: CYCLE when JSEQP<=0)
         IF (wdsrf_dn_pth(ipth) < -9000._r8 .or. rivelv_dn_pth(ipth) < -9000._r8) CYCLE

         rivelv_up = topo_rivelv(i_up)
         wdsrf_up  = wdsrf_ucat(i_up)
         rivelv_dn = rivelv_dn_pth(ipth)
         wdsrf_dn  = wdsrf_dn_pth(ipth)
         dt        = dt_all(irivsys(i_up))

         pth_hflux_total(ipth) = 0._r8

         i_dn = pth_down_local(ipth)

         DO ilev = 1, npthlev_bif

            width_pth = pth_wth(ilev, ipth)
            IF (width_pth <= 0._r8) CYCLE

            IF (ilev == 1) THEN
               wdsrf_up_eff = wdsrf_up
            ELSEIF (ilev > 1) THEN
               IF (DEF_USE_LEVEE .and. allocated(has_levee) .and. has_levee(i_up)) THEN
                  wdsrf_up_eff = levdph(i_up)
               ELSE
                  wdsrf_up_eff = wdsrf_up
               ENDIF
            ELSE
               wdsrf_up_eff = wdsrf_up
            ENDIF

            IF (ilev == 1) THEN
               wdsrf_dn_eff = wdsrf_dn
            ELSEIF (ilev > 1) THEN
               IF (DEF_USE_LEVEE .and. allocated(has_levee) .and. i_dn > 0 .and. i_dn <= numucat &
                  .and. has_levee(i_dn)) THEN
                  wdsrf_dn_eff = levdph(i_dn)
               ELSEIF (DEF_USE_LEVEE .and. has_levee_dn_pth(ipth) > 0.5_r8) THEN
                  wdsrf_dn_eff = levdph_dn_pth(ipth)
               ELSE
                  wdsrf_dn_eff = wdsrf_dn
               ENDIF
            ELSE
               wdsrf_dn_eff = wdsrf_dn
            ENDIF

            ! Effective bed elevation for this pathway layer
            ! Use max of upstream bed, downstream bed, and pathway layer elevation
            ! (same logic as main channel: bedelv_fc = max(rivelv_i, bedelv_next_i))
            bedelv_pth = max(rivelv_up, rivelv_dn, pth_elv(ilev, ipth))

            ! Flow depths at both ends relative to pathway bed
            height_up = max(0._r8, wdsrf_up_eff + rivelv_up - bedelv_pth)
            height_dn = max(0._r8, wdsrf_dn_eff + rivelv_dn - bedelv_pth)

            ! Skip if both ends are dry
            IF (height_up < BIFMIN .and. height_dn < BIFMIN) THEN
               hflux_lev = 0._r8
               pth_momen(ilev, ipth) = 0._r8
               pth_veloc(ilev, ipth) = 0._r8
               CYCLE
            ENDIF

            ! Upstream velocity from pathway state; downstream velocity = 0
            v_up = pth_veloc(ilev, ipth)
            v_dn = 0._r8

            ! --- Two-rarefaction Riemann solver ---

            ! Middle-state velocity and height
            veloct_fc = 0.5_r8 * (v_up + v_dn) &
               + sqrt(grav * height_up) - sqrt(grav * height_dn)

            height_fc = 1._r8/grav * (0.5_r8*(sqrt(grav*height_up) + sqrt(grav*height_dn)) &
               + 0.25_r8 * (v_up - v_dn)) ** 2

            ! Wave speeds
            IF (height_up > 0._r8) THEN
               vwave_up = min(v_up - sqrt(grav*height_up), veloct_fc - sqrt(grav*height_fc))
            ELSE
               vwave_up = v_dn - 2.0_r8 * sqrt(grav*height_dn)
            ENDIF

            IF (height_dn > 0._r8) THEN
               vwave_dn = max(v_dn + sqrt(grav*height_dn), veloct_fc + sqrt(grav*height_fc))
            ELSE
               vwave_dn = v_up + 2.0_r8 * sqrt(grav*height_up)
            ENDIF

            ! Fluxes at left and right states
            hflux_up = v_up  * height_up
            hflux_dn = v_dn  * height_dn
            mflux_up = v_up**2  * height_up + 0.5_r8*grav * height_up**2
            mflux_dn = v_dn**2  * height_dn + 0.5_r8*grav * height_dn**2

            ! Select flux based on wave structure, scaled by pathway width
            IF (vwave_up >= 0._r8) THEN
               hflux_lev = width_pth * hflux_up
               mflux_lev = width_pth * mflux_up
            ELSEIF (vwave_dn <= 0._r8) THEN
               hflux_lev = width_pth * hflux_dn
               mflux_lev = width_pth * mflux_dn
            ELSE
               hflux_lev = width_pth * (vwave_dn*hflux_up - vwave_up*hflux_dn &
                  + vwave_up*vwave_dn*(height_dn - height_up)) / (vwave_dn - vwave_up)
               mflux_lev = width_pth * (vwave_dn*mflux_up - vwave_up*mflux_dn &
                  + vwave_up*vwave_dn*(hflux_dn - hflux_up)) / (vwave_dn - vwave_up)
            ENDIF

            ! --- Update pathway momentum (semi-implicit friction, same as main channel) ---
            ! pth_are = pathway cross-section area (width * representative depth)
            ! Momentum equation: d(h*v)/dt = -d(h*v^2 + gh^2/2)/dx - friction
            ! Discrete: momen_new = (momen_old - mflux/are*dt) / (1 + friction*dt)
            pth_are = width_pth * max(height_up, height_dn, BIFMIN)

            IF (max(height_up, height_dn) >= BIFMIN) THEN
               friction = grav * pth_man(ilev)**2 &
                  / max(height_up, height_dn)**(7._r8/3._r8) * abs(pth_momen(ilev, ipth))
               pth_momen(ilev, ipth) = (pth_momen(ilev, ipth) &
                  - mflux_lev / pth_are * dt) &
                  / (1._r8 + friction * dt)
               pth_veloc(ilev, ipth) = pth_momen(ilev, ipth) / max(height_up, height_dn)
            ELSE
               pth_momen(ilev, ipth) = 0._r8
               pth_veloc(ilev, ipth) = 0._r8
            ENDIF

            ! Clamp velocity
            pth_veloc(ilev, ipth) = max(-20._r8, min(20._r8, pth_veloc(ilev, ipth)))

            ! Accumulate total flux for this pathway [m^3/s]
            pth_hflux_total(ipth) = pth_hflux_total(ipth) + hflux_lev
            bif_hflux_lev(ilev, ipth) = hflux_lev

         ENDDO  ! ilev

         ! ----- Step 4: 5% storage limiter -----
         IF (abs(pth_hflux_total(ipth)) > 1.e-10_r8 .and. dt > 0._r8) THEN
            storage_up = storage_ucat(i_up)
            storage_dn = storage_dn_pth(ipth)
            storage_ref = min(storage_up, storage_dn)
            storage_ref = max(storage_ref, 1._r8)
            rate = 0.05_r8 * storage_ref / (abs(pth_hflux_total(ipth)) * dt)
            IF (rate < 1._r8) THEN
               pth_hflux_total(ipth) = pth_hflux_total(ipth) * rate
               bif_hflux_lev(:, ipth) = bif_hflux_lev(:, ipth) * rate
               pth_momen(:, ipth) = pth_momen(:, ipth) * rate
               pth_veloc(:, ipth) = pth_veloc(:, ipth) * rate
               pth_hflux_total(ipth) = SUM(bif_hflux_lev(:, ipth))
            ENDIF
         ENDIF

         ! ----- Step 5: Accumulate to upstream ucat (local) -----
         bif_hflux_sum(i_up) = bif_hflux_sum(i_up) + pth_hflux_total(ipth)

         ! Also accumulate directly to downstream ucat if local
         i_dn = pth_down_local(ipth)
         IF (i_dn > 0 .and. i_dn <= numucat) THEN
            bif_hflux_sum(i_dn) = bif_hflux_sum(i_dn) - pth_hflux_total(ipth)
         ENDIF

      ENDDO  ! ipth

      ! ----- Step 6: Scatter flux to remote downstream ucats -----
      IF (numucat > 0) THEN
         bif_influx(:) = 0._r8
      ENDIF
      CALL worker_push_data (push_bif_influx, pth_hflux_total, bif_influx, &
         fillvalue = 0._r8, mode = 'sum')

      ! ----- Step 7: Subtract remote inflow, avoiding double-counting -----
      ! Remove contributions from local pathways (already counted in step 5)
      DO ipth = 1, npthout_local
         i_dn = pth_down_local(ipth)
         IF (i_dn > 0 .and. i_dn <= numucat) THEN
            bif_influx(i_dn) = bif_influx(i_dn) - pth_hflux_total(ipth)
         ENDIF
      ENDDO

      ! Apply remote inflow to bif_hflux_sum
      DO i_up = 1, numucat
         bif_hflux_sum(i_up) = bif_hflux_sum(i_up) - bif_influx(i_up)
      ENDDO

      ! Clean up
      deallocate (wdsrf_dn_pth)
      deallocate (rivelv_dn_pth)
      deallocate (levdph_dn_pth)
      deallocate (has_levee_dn_pth)
      deallocate (storage_dn_pth)
      deallocate (pth_hflux_total)
      deallocate (bif_influx)
      deallocate (has_levee_r8)
      deallocate (storage_ucat)

   END SUBROUTINE bifurcation_calc


   ! =========================================================================
   FUNCTION available_storage_ucat (i, wdsrf, volresv_in, is_built_resv_in) RESULT(storage)
   ! =========================================================================
   !
   ! Current host-model available storage semantics for limiter reference.
   !
   ! =========================================================================

   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: wdsrf
   real(r8), intent(in) :: volresv_in(:)
   logical,  intent(in) :: is_built_resv_in(:)
   real(r8)             :: storage

      storage = 0._r8

      IF (i < 1) RETURN

      IF (allocated(ucat2resv) .and. is_built_resv_in(i)) THEN
         IF (i <= size(ucat2resv)) THEN
            IF (ucat2resv(i) > 0 .and. ucat2resv(i) <= size(volresv_in)) THEN
               storage = volresv_in(ucat2resv(i))
               RETURN
            ENDIF
         ENDIF
      ENDIF

      IF (DEF_USE_LEVEE .and. allocated(has_levee) .and. has_levee(i)) THEN
         IF (volwater_ucat_valid .and. allocated(volwater_ucat)) THEN
            storage = volwater_ucat(i)
         ELSEIF (allocated(levsto) .and. levsto(i) > 0._r8) THEN
            IF (wdsrf <= floodplain_curve(i)%rivhgt + levee_hgt(i) + 1.e-6_r8) THEN
               storage = levee_topsto(i)
            ELSE
               storage = floodplain_curve(i)%volume (wdsrf) - levsto(i)
            ENDIF
            storage = max(storage, 0._r8)
         ELSE
            storage = floodplain_curve(i)%volume (wdsrf)
         ENDIF
      ELSE
         storage = floodplain_curve(i)%volume (wdsrf)
      ENDIF

   END FUNCTION available_storage_ucat


   ! =========================================================================
   SUBROUTINE read_bifurcation_restart (file_restart)
   ! =========================================================================
   !
   ! Read bifurcation pathway state from restart in global pathway order.
   !
   ! =========================================================================

   USE MOD_NetCDFSerial, only: ncio_var_exist
   USE MOD_Vector_ReadWrite, only: vector_read_matrix_and_scatter
   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, totalnpthout, pth_global_id
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical :: has_pth_veloc, has_pth_momen

      IF (.not. p_is_worker) RETURN
      IF (npthlev_bif <= 0 .or. npthout_local <= 0 .or. totalnpthout <= 0) RETURN

      has_pth_veloc = ncio_var_exist(file_restart, 'pth_veloc', readflag = .false.)
      has_pth_momen = ncio_var_exist(file_restart, 'pth_momen', readflag = .false.)

      if (has_pth_veloc .and. has_pth_momen) then
         CALL vector_read_matrix_and_scatter ( &
            file_restart, pth_veloc, npthlev_bif, npthout_local, 'pth_veloc', pth_global_id, totalnpthout)
         CALL vector_read_matrix_and_scatter ( &
            file_restart, pth_momen, npthlev_bif, npthout_local, 'pth_momen', pth_global_id, totalnpthout)
      else
         ! Missing pathway variables mean the restart file predates bifurcation
         ! support, so keep the zero-initialized cold-start state for compatibility.
         pth_veloc(:,:) = 0._r8
         pth_momen(:,:) = 0._r8
      endif

   END SUBROUTINE read_bifurcation_restart


   ! =========================================================================
   SUBROUTINE write_bifurcation_restart (file_restart)
   ! =========================================================================
   !
   ! Write bifurcation pathway state in global pathway order.
   !
   ! =========================================================================

   USE MOD_Namelist, only: DEF_REST_CompressLevel
   USE MOD_NetCDFSerial, only: ncio_define_dimension, ncio_write_serial
   USE MOD_Vector_ReadWrite, only: vector_gather_matrix_to_master
   USE MOD_Grid_RiverLakeNetwork, only: npthlev_bif, npthout_local, totalnpthout, pth_global_id
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   integer :: ncol_local_write
   integer, allocatable :: pth_global_id_write(:)
   real(r8), allocatable :: wdata(:,:)
   real(r8), allocatable :: pth_veloc_write(:,:), pth_momen_write(:,:)

      IF (totalnpthout <= 0 .or. npthlev_bif <= 0) RETURN

      ! Guard: bifurcation state may not be initialised (e.g. mkinidata).
      ! Write zero state so the restart file has valid variables.
      IF (p_is_worker) THEN
         IF (allocated(pth_veloc) .and. allocated(pth_momen)) THEN
            ncol_local_write = npthout_local
            pth_veloc_write = pth_veloc
            pth_momen_write = pth_momen
            allocate (pth_global_id_write (ncol_local_write))
            pth_global_id_write(:) = pth_global_id(:)
         ELSE
            ! Not initialised: contribute zero-state columns
            ncol_local_write = npthout_local
            allocate (pth_veloc_write     (npthlev_bif, ncol_local_write))
            allocate (pth_momen_write     (npthlev_bif, ncol_local_write))
            allocate (pth_global_id_write (ncol_local_write))
            pth_veloc_write = 0._r8
            pth_momen_write = 0._r8
            pth_global_id_write(:) = pth_global_id(:)
         ENDIF
      ELSE
         ncol_local_write = 0
         allocate (pth_veloc_write     (npthlev_bif, 0))
         allocate (pth_momen_write     (npthlev_bif, 0))
         allocate (pth_global_id_write (0))
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_define_dimension(file_restart, 'bifurcation_level',   npthlev_bif)
         CALL ncio_define_dimension(file_restart, 'bifurcation_pathway', totalnpthout)
      ENDIF

      ! Gather local (npthlev_bif, npthout_local) state into global
      ! (npthlev_bif, totalnpthout) pathway order using pth_global_id.
      CALL vector_gather_matrix_to_master ( &
         pth_veloc_write, npthlev_bif, ncol_local_write, totalnpthout, pth_global_id_write, wdata)

      IF (p_is_master) THEN
         CALL ncio_write_serial (file_restart, 'pth_veloc', wdata, &
            'bifurcation_level', 'bifurcation_pathway', DEF_REST_CompressLevel)
         deallocate (wdata)
      ENDIF

      CALL vector_gather_matrix_to_master ( &
         pth_momen_write, npthlev_bif, ncol_local_write, totalnpthout, pth_global_id_write, wdata)

      IF (p_is_master) THEN
         CALL ncio_write_serial (file_restart, 'pth_momen', wdata, &
            'bifurcation_level', 'bifurcation_pathway', DEF_REST_CompressLevel)
         deallocate (wdata)
      ENDIF

      deallocate (pth_veloc_write)
      deallocate (pth_momen_write)
      deallocate (pth_global_id_write)

   END SUBROUTINE write_bifurcation_restart


   ! =========================================================================
   SUBROUTINE bifurcation_final ()
   ! =========================================================================
   !
   ! Deallocate all bifurcation arrays.
   !
   ! =========================================================================

   IMPLICIT NONE

      IF (allocated(pth_veloc))     deallocate (pth_veloc)
      IF (allocated(pth_momen))     deallocate (pth_momen)
      IF (allocated(bif_hflux_lev)) deallocate (bif_hflux_lev)
      IF (allocated(bif_hflux_sum)) deallocate (bif_hflux_sum)

   END SUBROUTINE bifurcation_final

END MODULE MOD_Grid_RiverLakeBifurcation
#endif
