#include <define.h>

MODULE MOD_Tracer_Conservation

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, trc_tiny, delta_to_R
   USE MOD_Tracer_Vars

   IMPLICIT NONE

   real(r8), parameter :: trc_balance_tol = 1.0e-10_r8

   ! Per-step accumulator snapshots (saved at start of each timestep)
   real(r8), allocatable, save :: snap_precip(:,:)  ! (ntracers, numpatch)
   real(r8), allocatable, save :: snap_evap  (:,:)
   real(r8), allocatable, save :: snap_rnof  (:,:)

   ! Worker-local worst-patch tracker for tracer_balance_check. Previously
   ! the routine only printed when `ipatch <= 1 .and. itrc == 1`, so every
   ! other failing patch was invisible. These counters accumulate across
   ! patches within a step; caller drains them via tracer_balance_report
   ! at step end.
   real(r8), save :: balance_worst_err   = 0._r8
   integer,  save :: balance_worst_ipatch = 0
   integer,  save :: balance_worst_itrc   = 0
   integer,  save :: balance_worst_ptype  = -1
   integer,  save :: balance_nbad         = 0

   PUBLIC :: deallocate_tracer_conservation
   PUBLIC :: tracer_balance_report

CONTAINS

   !---------------------------------------------------------------
   ! Release the save-level snap_* snapshots. Must be called from
   ! tracer_final so that a subsequent tracer_init on the same
   ! process (regression harness, embedded driver) does not re-use
   ! stale arrays sized to the previous numpatch — the lazy
   ! `IF (.not. allocated) allocate` guard in tracer_save_storage
   ! only fires on the first call ever, so without this deallocate
   ! a changed numpatch would either truncate or overrun snap_*.
   !---------------------------------------------------------------
   SUBROUTINE deallocate_tracer_conservation ()
      IMPLICIT NONE
      IF (allocated(snap_precip)) deallocate(snap_precip)
      IF (allocated(snap_evap  )) deallocate(snap_evap  )
      IF (allocated(snap_rnof  )) deallocate(snap_rnof  )
   END SUBROUTINE deallocate_tracer_conservation

   SUBROUTINE tracer_save_storage (ipatch, snl, nl_soil, waterstorage)
      IMPLICIT NONE
      integer, intent(in) :: ipatch, snl, nl_soil
      ! Optional: patch-level irrigation reservoir (mm of water equivalent).
      ! Present only when DEF_USE_IRRIGATION is on under CROP. When passed,
      ! trc_waterstorage is re-synced to waterstorage*R_init (Phase-1
      ! invariant) and included in the storage balance; this mirrors the
      ! water-side convention at CoLMMAIN:806 that treats waterstorage as
      ! part of the column inventory (so irrigation becomes an internal
      ! transfer instead of an external atmospheric input).
      real(r8), intent(in), optional :: waterstorage
      integer :: itrc, j
      real(r8) :: R_init

      IF (ntracers <= 0) RETURN

      ! Allocate snapshots on first call
      IF (.not. allocated(snap_precip)) THEN
         allocate(snap_precip(ntracers, size(trc_storage_beg,2))); snap_precip = 0._r8
         allocate(snap_evap  (ntracers, size(trc_storage_beg,2))); snap_evap   = 0._r8
         allocate(snap_rnof  (ntracers, size(trc_storage_beg,2))); snap_rnof   = 0._r8
      ENDIF

      DO itrc = 1, ntracers
         trc_rnof_step(itrc, ipatch) = 0._r8

         ! Phase-1 re-sync of the irrigation reservoir tracer to current
         ! waterstorage. Under Phase 1 all refills arrive at R_init so
         ! the ratio is constant; Phase 2 (time-varying refill R) will
         ! have to track refills explicitly and drop this re-sync.
         IF (present(waterstorage) .and. allocated(trc_waterstorage)) THEN
            R_init = delta_to_R(tracers(itrc)%init_delta, tracers(itrc)%ref_ratio)
            trc_waterstorage(itrc, ipatch) = max(waterstorage, 0._r8) * R_init
         ENDIF

         ! Save current storage
         trc_storage_beg(itrc, ipatch) = 0._r8
         trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
            + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
         DO j = snl + 1, nl_soil
            trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
         trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
            + trc_wa(itrc, ipatch) + trc_wdsrf(itrc, ipatch) + trc_wetwat(itrc, ipatch) &
            + trc_scv(itrc, ipatch)
         IF (allocated(trc_waterstorage)) THEN
            trc_storage_beg(itrc, ipatch) = trc_storage_beg(itrc, ipatch) &
               + trc_waterstorage(itrc, ipatch)
         ENDIF

         ! Snapshot accumulators at start of this step
         snap_precip(itrc, ipatch) = a_trc_precip(itrc, ipatch)
         snap_evap  (itrc, ipatch) = a_trc_evap  (itrc, ipatch)
         snap_rnof  (itrc, ipatch) = a_trc_rnof  (itrc, ipatch)
      ENDDO
   END SUBROUTINE tracer_save_storage

   SUBROUTINE tracer_balance_check (ipatch, snl, nl_soil, deltim, xerr_tracer, patchtype_in)
      IMPLICIT NONE
      integer,  intent(in)  :: ipatch, snl, nl_soil
      real(r8), intent(in)  :: deltim
      real(r8), intent(out) :: xerr_tracer
      ! Optional patchtype passthrough so TRC_BAL report can identify the
      ! offending branch (0=soil, 1=urban, 2=wetland, 3=glacier, 4=waterbody).
      ! Kept optional so old call sites (if any) still compile.
      integer,  intent(in), optional :: patchtype_in

      integer  :: itrc, j
      real(r8) :: storage_end, step_input, step_output, err
      real(r8) :: R_input, trc_expected, trc_actual, max_drift
      real(r8) :: w_pool, t_pool, drift

      xerr_tracer = 0._r8
      IF (ntracers <= 0) RETURN

      DO itrc = 1, ntracers
         ! Current total storage
         storage_end = 0._r8
         storage_end = storage_end + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
         DO j = snl + 1, nl_soil
            storage_end = storage_end &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
         storage_end = storage_end &
            + trc_wa(itrc, ipatch) + trc_wdsrf(itrc, ipatch) + trc_wetwat(itrc, ipatch) &
            + trc_scv(itrc, ipatch)
         IF (allocated(trc_waterstorage)) THEN
            storage_end = storage_end + trc_waterstorage(itrc, ipatch)
         ENDIF

         ! Per-step fluxes = current accumulator - snapshot at step start
         step_input  = a_trc_precip(itrc, ipatch) - snap_precip(itrc, ipatch)
         step_output = (a_trc_evap(itrc, ipatch) - snap_evap(itrc, ipatch)) &
                     + (a_trc_rnof(itrc, ipatch) - snap_rnof(itrc, ipatch))

         ! Conservation: Δstorage = input - output
         err = storage_end - trc_storage_beg(itrc, ipatch) - step_input + step_output
         trc_balance_err(itrc, ipatch) = err

         IF (abs(err) > trc_balance_tol) THEN
            ! Track worst-ever failure within this step (worker-local).
            ! The existing per-patch WRITE below stays narrow (patch 1,
            ! tracer 1) to avoid log flooding; tracer_balance_report
            ! surfaces the hidden bulk at end of step.
            IF (abs(err) > abs(balance_worst_err)) THEN
               balance_worst_err    = err
               balance_worst_ipatch = ipatch
               balance_worst_itrc   = itrc
               IF (present(patchtype_in)) THEN
                  balance_worst_ptype = patchtype_in
               ELSE
                  balance_worst_ptype = -1
               ENDIF
            ENDIF
            balance_nbad = balance_nbad + 1
            IF (ipatch <= 1 .and. itrc == 1) THEN
               WRITE(*,'(A,I6,A,E12.5,A,E12.5,A,E12.5,A,E12.5)') &
                  'TRC_ERR p=', ipatch, &
                  ' dS=', storage_end - trc_storage_beg(itrc, ipatch), &
                  ' in=', step_input, ' out=', step_output, ' err=', err
            ENDIF
         ENDIF

         xerr_tracer = max(xerr_tracer, abs(err) / deltim)
      ENDDO
   END SUBROUTINE tracer_balance_check

   !---------------------------------------------------------------
   ! Drain worker-local balance counters. Call once per CoLM step
   ! after all patches have run tracer_balance_check. Reduces the
   ! worst-err across workers so a single-patch failure elsewhere
   ! still lands in the per-rank log. Zeros the state so next step
   ! starts fresh.
   !
   ! This runs from inside CoLMDRIVER (worker-only gated at
   ! CoLM.F90:528). Reducing over p_comm_glb here would hang — master
   ! and IO never enter. Use p_comm_worker so the collective stays on
   ! the same ranks that produced the counts. Worker rank 0 prints the
   ! summary (not p_is_master, which lives on a different comm).
   !---------------------------------------------------------------
   SUBROUTINE tracer_balance_report ()
#ifdef USEMPI
      USE MOD_SPMD_Task, only: p_comm_worker, p_iam_worker, p_err, p_is_worker, &
                               p_np_worker
#endif
      IMPLICIT NONE
      real(r8) :: worst_abs, reduced_abs
      integer  :: nbad_total
      integer  :: winner_rank, local_winner
      logical  :: print_me
      ! MAXLOC reduction: pair |err| with rank so we can ship ipatch/itrc
      ! from the rank that owns the worst patch. Standard MPI_MAXLOC on
      ! real+int works with MPI_2DOUBLE_PRECISION wrapping (rank as real).
      real(r8) :: in_pair(2), out_pair(2)
#ifdef USEMPI
      include 'mpif.h'
#endif

      worst_abs = abs(balance_worst_err)
      nbad_total  = balance_nbad

#ifdef USEMPI
      reduced_abs = worst_abs
      IF (p_is_worker) THEN
         CALL mpi_reduce(worst_abs,    reduced_abs, 1, MPI_REAL8,   MPI_MAX, 0, p_comm_worker, p_err)
         CALL mpi_reduce(balance_nbad, nbad_total,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_worker, p_err)
         ! Broadcast reduced_abs back so every worker can decide whether
         ! to contribute its ipatch/itrc (only the owner of worst_abs
         ! ships it via another reduce).
         CALL mpi_bcast(reduced_abs, 1, MPI_REAL8, 0, p_comm_worker, p_err)

         ! Two-stage: pack (|err|, rank_as_real) and MAXLOC to find owner.
         in_pair(1) = worst_abs
         in_pair(2) = real(p_iam_worker, r8)
         out_pair   = in_pair
         CALL mpi_reduce(in_pair, out_pair, 1, MPI_2DOUBLE_PRECISION, &
                         MPI_MAXLOC, 0, p_comm_worker, p_err)
         winner_rank = nint(out_pair(2))
         ! Broadcast winner rank to everyone so the owner can ship
         ! ipatch/itrc to worker 0 via a point-to-point (cheap vs another reduce).
         CALL mpi_bcast(winner_rank, 1, MPI_INTEGER, 0, p_comm_worker, p_err)
         IF (winner_rank /= 0) THEN
            IF (p_iam_worker == winner_rank) THEN
               CALL mpi_send(balance_worst_ipatch, 1, MPI_INTEGER, 0, 91, p_comm_worker, p_err)
               CALL mpi_send(balance_worst_itrc,   1, MPI_INTEGER, 0, 92, p_comm_worker, p_err)
               CALL mpi_send(balance_worst_ptype,  1, MPI_INTEGER, 0, 93, p_comm_worker, p_err)
            ELSEIF (p_iam_worker == 0) THEN
               CALL mpi_recv(balance_worst_ipatch, 1, MPI_INTEGER, winner_rank, 91, &
                             p_comm_worker, MPI_STATUS_IGNORE, p_err)
               CALL mpi_recv(balance_worst_itrc,   1, MPI_INTEGER, winner_rank, 92, &
                             p_comm_worker, MPI_STATUS_IGNORE, p_err)
               CALL mpi_recv(balance_worst_ptype,  1, MPI_INTEGER, winner_rank, 93, &
                             p_comm_worker, MPI_STATUS_IGNORE, p_err)
            ENDIF
         ENDIF
      ENDIF
      print_me = (p_is_worker .and. p_iam_worker == 0)
#else
      reduced_abs = worst_abs
      print_me    = .true.
#endif

      IF (print_me .and. nbad_total > 0) THEN
         ! "nbad_entries" = (patch, tracer) pairs exceeding tol.
         ! Upper bound on unique patches is nbad_entries / ntracers
         ! (each patch contributes at most ntracers entries).
         ! ptype: 0=soil 1=urban 2=wetland 3=glacier 4=waterbody, -1=unknown.
         WRITE(*,'(A,I8,A,E12.5,A,I8,A,I4,A,I3)') &
            'TRC_BAL step report: nbad_entries=', nbad_total, &
            ' worst_abs_err=', reduced_abs, &
            ' @ipatch=', balance_worst_ipatch, &
            ' itrc=', balance_worst_itrc, &
            ' ptype=', balance_worst_ptype
      ENDIF

      balance_worst_err    = 0._r8
      balance_worst_ipatch = 0
      balance_worst_itrc   = 0
      balance_worst_ptype  = -1
      balance_nbad         = 0
   END SUBROUTINE tracer_balance_report

END MODULE MOD_Tracer_Conservation
