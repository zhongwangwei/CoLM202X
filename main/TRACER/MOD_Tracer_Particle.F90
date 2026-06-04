#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Particle
!=======================================================================
! Particle-tracer dispatch layer.
!
! HYDRO owns water routing and calls this generic layer. Particle species
! such as suspended sediment register TRACER-owned callbacks below it, so
! HYDRO does not depend on species-private sediment modules.
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Tracer_Defs, only: tracer_lower
   IMPLICIT NONE

   interface
      SUBROUTINE register_all_particle_callbacks ()
      END SUBROUTINE register_all_particle_callbacks
   end interface

   abstract interface
      logical FUNCTION particle_has_if ()
      END FUNCTION particle_has_if

      SUBROUTINE particle_noarg_if ()
      END SUBROUTINE particle_noarg_if

      SUBROUTINE particle_refresh_if ()
      END SUBROUTINE particle_refresh_if

      SUBROUTINE particle_read_restart_if (file_restart)
         character(len=*), intent(in) :: file_restart
      END SUBROUTINE particle_read_restart_if

      SUBROUTINE particle_write_restart_if (file_restart)
         character(len=*), intent(in) :: file_restart
      END SUBROUTINE particle_write_restart_if

      SUBROUTINE particle_forcing_put_if (precip, dt)
         USE MOD_Precision
         real(r8), intent(in) :: precip(:)
         real(r8), intent(in) :: dt
      END SUBROUTINE particle_forcing_put_if

      SUBROUTINE particle_diag_accumulate_if (dt_all, irivsys, ucatfilter, &
         veloc, wdsrf, rivsto_input, rivout_fc, floodarea)
         USE MOD_Precision
         real(r8), intent(in) :: dt_all(:)
         integer,  intent(in) :: irivsys(:)
         logical,  intent(in) :: ucatfilter(:)
         real(r8), intent(in) :: veloc(:)
         real(r8), intent(in) :: wdsrf(:)
         real(r8), intent(in) :: rivsto_input(:)
         real(r8), intent(in) :: rivout_fc(:)
         real(r8), intent(in) :: floodarea(:)
      END SUBROUTINE particle_diag_accumulate_if

      SUBROUTINE particle_calc_if (deltime)
         USE MOD_Precision
         real(r8), intent(in) :: deltime
      END SUBROUTINE particle_calc_if

      SUBROUTINE particle_history_if (file_hist_ucat, itime_in_file_ucat)
         character(len=*), intent(in) :: file_hist_ucat
         integer, intent(in) :: itime_in_file_ucat
      END SUBROUTINE particle_history_if

      SUBROUTINE particle_remap_lulcc_if (patchclass_new, eindex_new, patchclass_old, eindex_old, &
         lccpct_patches, old_patch_area)
         USE MOD_Precision
         integer, intent(in) :: patchclass_new(:), patchclass_old(:)
         integer*8, intent(in) :: eindex_new(:), eindex_old(:)
         real(r8), intent(in), optional :: lccpct_patches(:,:)
         real(r8), intent(in), optional :: old_patch_area(:)
      END SUBROUTINE particle_remap_lulcc_if
   end interface

   PRIVATE

   integer, parameter :: PARTICLE_CALLBACKS_INITIAL_CAPACITY = 2

   type :: particle_callbacks_type
      character(len=32) :: name = ''
      procedure(particle_has_if),             pointer, nopass :: has => null()
      procedure(particle_refresh_if),         pointer, nopass :: refresh => null()
      procedure(particle_noarg_if),           pointer, nopass :: init => null()
      procedure(particle_read_restart_if),    pointer, nopass :: read_restart => null()
      procedure(particle_forcing_put_if),     pointer, nopass :: forcing_put => null()
      procedure(particle_diag_accumulate_if), pointer, nopass :: diag_accumulate => null()
      procedure(particle_calc_if),            pointer, nopass :: calc => null()
      procedure(particle_history_if),         pointer, nopass :: history => null()
      procedure(particle_noarg_if),           pointer, nopass :: flush_history => null()
      procedure(particle_write_restart_if),   pointer, nopass :: write_restart => null()
      procedure(particle_noarg_if),           pointer, nopass :: save_lulcc => null()
      procedure(particle_remap_lulcc_if),     pointer, nopass :: remap_lulcc => null()
      procedure(particle_noarg_if),           pointer, nopass :: final => null()
   end type particle_callbacks_type

   type(particle_callbacks_type), allocatable, save :: particle_callbacks(:)
   integer, save :: n_particle_callbacks = 0
   logical, save :: particle_callbacks_ready = .false.
   logical, save :: particle_refresh_dirty = .true.

   PUBLIC :: register_particle_callbacks
   PUBLIC :: mark_particle_callbacks_dirty
   PUBLIC :: tracer_particle_has_active
   PUBLIC :: tracer_particle_init, tracer_particle_final
   PUBLIC :: tracer_particle_read_restart, tracer_particle_write_restart
   PUBLIC :: tracer_particle_forcing_put, tracer_particle_diag_accumulate
   PUBLIC :: tracer_particle_calc
   PUBLIC :: tracer_particle_write_history, tracer_particle_flush_history
   PUBLIC :: tracer_particle_save_lulcc_state, tracer_particle_remap_lulcc_state

CONTAINS

   SUBROUTINE ensure_particle_callbacks_registered ()

      IMPLICIT NONE

      IF (particle_callbacks_ready) RETURN
      particle_callbacks_ready = .true.

      CALL register_all_particle_callbacks ()

   END SUBROUTINE ensure_particle_callbacks_registered

   SUBROUTINE prepare_particle_dispatch ()

      IMPLICIT NONE

      CALL ensure_particle_callbacks_registered ()
      IF (particle_refresh_dirty) CALL refresh_particle_callback_states_if_dirty ()

   END SUBROUTINE prepare_particle_dispatch

   SUBROUTINE mark_particle_callbacks_dirty ()

      IMPLICIT NONE

      particle_refresh_dirty = .true.

   END SUBROUTINE mark_particle_callbacks_dirty

   SUBROUTINE register_particle_callbacks (name, has_fn, refresh_fn, init_fn, read_restart_fn, &
      forcing_put_fn, diag_accumulate_fn, calc_fn, history_fn, flush_history_fn, &
      write_restart_fn, save_lulcc_fn, remap_lulcc_fn, final_fn)

      IMPLICIT NONE
      character(len=*), intent(in) :: name
      procedure(particle_has_if),             optional :: has_fn
      procedure(particle_refresh_if),         optional :: refresh_fn
      procedure(particle_noarg_if),           optional :: init_fn
      procedure(particle_read_restart_if),    optional :: read_restart_fn
      procedure(particle_forcing_put_if),     optional :: forcing_put_fn
      procedure(particle_diag_accumulate_if), optional :: diag_accumulate_fn
      procedure(particle_calc_if),            optional :: calc_fn
      procedure(particle_history_if),         optional :: history_fn
      procedure(particle_noarg_if),           optional :: flush_history_fn
      procedure(particle_write_restart_if),   optional :: write_restart_fn
      procedure(particle_noarg_if),           optional :: save_lulcc_fn
      procedure(particle_remap_lulcc_if),     optional :: remap_lulcc_fn
      procedure(particle_noarg_if),           optional :: final_fn

      integer :: idx

      IF (len_trim(name) <= 0) THEN
         CALL CoLM_stop ('MOD_Tracer_Particle: cannot register particle tracer with empty name')
      ENDIF
      IF (.not. present(has_fn)) THEN
         CALL CoLM_stop ('MOD_Tracer_Particle: particle tracer registration requires has_fn')
      ENDIF
      IF (particle_registration_conflicts(name)) THEN
         CALL CoLM_stop ('MOD_Tracer_Particle: duplicate particle tracer name registration')
      ENDIF

      CALL ensure_particle_callback_capacity (n_particle_callbacks + 1)
      n_particle_callbacks = n_particle_callbacks + 1
      idx = n_particle_callbacks

      particle_callbacks(idx)%name = trim(name)
      particle_callbacks(idx)%has => has_fn
      IF (present(refresh_fn))         particle_callbacks(idx)%refresh => refresh_fn
      IF (present(init_fn))            particle_callbacks(idx)%init => init_fn
      IF (present(read_restart_fn))    particle_callbacks(idx)%read_restart => read_restart_fn
      IF (present(forcing_put_fn))     particle_callbacks(idx)%forcing_put => forcing_put_fn
      IF (present(diag_accumulate_fn)) particle_callbacks(idx)%diag_accumulate => diag_accumulate_fn
      IF (present(calc_fn))            particle_callbacks(idx)%calc => calc_fn
      IF (present(history_fn))         particle_callbacks(idx)%history => history_fn
      IF (present(flush_history_fn))   particle_callbacks(idx)%flush_history => flush_history_fn
      IF (present(write_restart_fn))   particle_callbacks(idx)%write_restart => write_restart_fn
      IF (present(save_lulcc_fn))      particle_callbacks(idx)%save_lulcc => save_lulcc_fn
      IF (present(remap_lulcc_fn))     particle_callbacks(idx)%remap_lulcc => remap_lulcc_fn
      IF (present(final_fn))           particle_callbacks(idx)%final => final_fn
      CALL mark_particle_callbacks_dirty ()

   END SUBROUTINE register_particle_callbacks

   SUBROUTINE refresh_particle_callback_states_if_dirty ()

      IMPLICIT NONE
      integer :: i

      IF (.not. particle_refresh_dirty) RETURN
      particle_refresh_dirty = .false.
      IF (.not. allocated(particle_callbacks)) RETURN

      DO i = 1, n_particle_callbacks
         IF (associated(particle_callbacks(i)%refresh)) CALL particle_callbacks(i)%refresh ()
      ENDDO

   END SUBROUTINE refresh_particle_callback_states_if_dirty

   SUBROUTINE ensure_particle_callback_capacity (needed)

      IMPLICIT NONE
      integer, intent(in) :: needed

      type(particle_callbacks_type), allocatable :: grown(:)
      integer :: old_size, new_size

      IF (.not. allocated(particle_callbacks)) THEN
         allocate(particle_callbacks(max(PARTICLE_CALLBACKS_INITIAL_CAPACITY, needed)))
         RETURN
      ENDIF

      old_size = size(particle_callbacks)
      IF (needed <= old_size) RETURN

      new_size = max(needed, old_size * 2)
      allocate(grown(new_size))
      IF (n_particle_callbacks > 0) grown(1:n_particle_callbacks) = particle_callbacks(1:n_particle_callbacks)
      CALL move_alloc(grown, particle_callbacks)

   END SUBROUTINE ensure_particle_callback_capacity

   SUBROUTINE clear_particle_callbacks ()

      IMPLICIT NONE

      IF (allocated(particle_callbacks)) deallocate(particle_callbacks)
      n_particle_callbacks = 0
      particle_callbacks_ready = .false.
      particle_refresh_dirty = .true.

   END SUBROUTINE clear_particle_callbacks

   logical FUNCTION particle_registration_conflicts (name)

      IMPLICIT NONE
      character(len=*), intent(in) :: name
      integer :: i

      particle_registration_conflicts = .false.
      IF (.not. allocated(particle_callbacks)) RETURN

      DO i = 1, n_particle_callbacks
         IF (trim(tracer_lower(particle_callbacks(i)%name)) == trim(tracer_lower(name))) THEN
            particle_registration_conflicts = .true.
            RETURN
         ENDIF
      ENDDO

   END FUNCTION particle_registration_conflicts

   logical FUNCTION particle_callback_enabled (idx)

      IMPLICIT NONE
      integer, intent(in) :: idx

      particle_callback_enabled = .false.
      IF (idx < 1 .or. idx > n_particle_callbacks) RETURN
      IF (.not. associated(particle_callbacks(idx)%has)) RETURN
      particle_callback_enabled = particle_callbacks(idx)%has()

   END FUNCTION particle_callback_enabled

   logical FUNCTION tracer_particle_has_active ()

      IMPLICIT NONE
      integer :: i

      CALL prepare_particle_dispatch ()
      tracer_particle_has_active = .false.

      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i)) THEN
            tracer_particle_has_active = .true.
            RETURN
         ENDIF
      ENDDO

   END FUNCTION tracer_particle_has_active

   SUBROUTINE tracer_particle_init ()

      IMPLICIT NONE
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%init)) &
            CALL particle_callbacks(i)%init ()
      ENDDO

   END SUBROUTINE tracer_particle_init

   SUBROUTINE tracer_particle_read_restart (file_restart)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%read_restart)) &
            CALL particle_callbacks(i)%read_restart (file_restart)
      ENDDO

   END SUBROUTINE tracer_particle_read_restart

   SUBROUTINE tracer_particle_write_restart (file_restart)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%write_restart)) &
            CALL particle_callbacks(i)%write_restart (file_restart)
      ENDDO

   END SUBROUTINE tracer_particle_write_restart

   SUBROUTINE tracer_particle_forcing_put (precip, dt)

      IMPLICIT NONE
      real(r8), intent(in) :: precip(:)
      real(r8), intent(in) :: dt
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%forcing_put)) &
            CALL particle_callbacks(i)%forcing_put (precip, dt)
      ENDDO

   END SUBROUTINE tracer_particle_forcing_put

   SUBROUTINE tracer_particle_diag_accumulate (dt_all, irivsys, ucatfilter, &
      veloc, wdsrf, rivsto_input, rivout_fc, floodarea)

      IMPLICIT NONE
      real(r8), intent(in) :: dt_all(:)
      integer,  intent(in) :: irivsys(:)
      logical,  intent(in) :: ucatfilter(:)
      real(r8), intent(in) :: veloc(:)
      real(r8), intent(in) :: wdsrf(:)
      real(r8), intent(in) :: rivsto_input(:)
      real(r8), intent(in) :: rivout_fc(:)
      real(r8), intent(in) :: floodarea(:)
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%diag_accumulate)) &
            CALL particle_callbacks(i)%diag_accumulate (dt_all, irivsys, ucatfilter, &
               veloc, wdsrf, rivsto_input, rivout_fc, floodarea)
      ENDDO

   END SUBROUTINE tracer_particle_diag_accumulate

   SUBROUTINE tracer_particle_calc (deltime)

      IMPLICIT NONE
      real(r8), intent(in) :: deltime
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%calc)) &
            CALL particle_callbacks(i)%calc (deltime)
      ENDDO

   END SUBROUTINE tracer_particle_calc

   SUBROUTINE tracer_particle_write_history (file_hist_ucat, itime_in_file_ucat)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_hist_ucat
      integer, intent(in) :: itime_in_file_ucat
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%history)) &
            CALL particle_callbacks(i)%history (file_hist_ucat, itime_in_file_ucat)
      ENDDO

   END SUBROUTINE tracer_particle_write_history

   SUBROUTINE tracer_particle_flush_history ()

      IMPLICIT NONE
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%flush_history)) &
            CALL particle_callbacks(i)%flush_history ()
      ENDDO

   END SUBROUTINE tracer_particle_flush_history

   SUBROUTINE tracer_particle_save_lulcc_state ()

      IMPLICIT NONE
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%save_lulcc)) &
            CALL particle_callbacks(i)%save_lulcc ()
      ENDDO

   END SUBROUTINE tracer_particle_save_lulcc_state

   SUBROUTINE tracer_particle_remap_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
      lccpct_patches, old_patch_area)

      IMPLICIT NONE
      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8), intent(in), optional :: lccpct_patches(:,:)
      real(r8), intent(in), optional :: old_patch_area(:)
      integer :: i

      CALL prepare_particle_dispatch ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%remap_lulcc)) THEN
            CALL particle_callbacks(i)%remap_lulcc (patchclass_new, eindex_new, patchclass_old, eindex_old, &
               lccpct_patches, old_patch_area)
         ENDIF
      ENDDO
      CALL mark_particle_callbacks_dirty ()

   END SUBROUTINE tracer_particle_remap_lulcc_state

   SUBROUTINE tracer_particle_final ()

      IMPLICIT NONE
      integer :: i

      CALL ensure_particle_callbacks_registered ()
      DO i = 1, n_particle_callbacks
         IF (particle_callback_enabled(i) .and. associated(particle_callbacks(i)%final)) &
            CALL particle_callbacks(i)%final ()
      ENDDO
      CALL clear_particle_callbacks ()

   END SUBROUTINE tracer_particle_final

END MODULE MOD_Tracer_Particle
#endif
