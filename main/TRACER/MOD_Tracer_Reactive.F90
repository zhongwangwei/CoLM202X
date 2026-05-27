#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Reactive
!=======================================================================
! Reactive-tracer dispatch layer.
!
! The generic TRACER lifecycle calls this module. Species-specific
! reactive implementations register callbacks below this layer; the
! generic entry points only iterate the registry.
!=======================================================================

   USE MOD_Precision
   USE MOD_DataType, only: block_data_real8_2d
   USE MOD_SPMD_Task, only: CoLM_stop
   USE MOD_Tracer_Defs, only: ntracers, tracers, tracer_is_reactive
   IMPLICIT NONE

   interface
      SUBROUTINE register_all_reactive_callbacks ()
      END SUBROUTINE register_all_reactive_callbacks
   end interface

   abstract interface
      logical FUNCTION reactive_has_if (tracer_name)
         character(len=*), intent(in) :: tracer_name
      END FUNCTION reactive_has_if

      SUBROUTINE reactive_init_if (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)
         integer, intent(in) :: numpatch
         integer, intent(in) :: lc_year
         integer, intent(in) :: jdate(3)
         character(len=*), intent(in) :: casename
         character(len=*), intent(in) :: dir_restart
         character(len=*), intent(in) :: dir_landdata
      END SUBROUTINE reactive_init_if

      SUBROUTINE reactive_read_restart_if (file_restart)
         character(len=*), intent(in) :: file_restart
      END SUBROUTINE reactive_read_restart_if

      SUBROUTINE reactive_write_restart_if (file_restart, compress)
         character(len=*), intent(in) :: file_restart
         integer, intent(in) :: compress
      END SUBROUTINE reactive_write_restart_if

      SUBROUTINE reactive_lake_step_if (istep_local, ipatch, idate, deltim_phy, isub, nsub)
         USE MOD_Precision
         integer,  intent(in) :: istep_local
         integer,  intent(in) :: ipatch
         integer,  intent(in) :: idate(3)
         real(r8), intent(in) :: deltim_phy
         integer,  intent(in) :: isub
         integer,  intent(in) :: nsub
      END SUBROUTINE reactive_lake_step_if

      SUBROUTINE reactive_wetland_decomp_if (ipatch)
         integer, intent(in) :: ipatch
      END SUBROUTINE reactive_wetland_decomp_if

      SUBROUTINE reactive_soil_step_if (istep_local, ipatch, idate, deltim)
         USE MOD_Precision
         integer,  intent(in) :: istep_local
         integer,  intent(in) :: ipatch
         integer,  intent(in) :: idate(3)
         real(r8), intent(in) :: deltim
      END SUBROUTINE reactive_soil_step_if

      SUBROUTINE reactive_noarg_if ()
      END SUBROUTINE reactive_noarg_if

      SUBROUTINE reactive_history_if (file_hist, itime_in_file, sumarea, filter, &
         nl_soil, forcing_has_missing_value, forcmask_pch)
         USE MOD_DataType, only: block_data_real8_2d
         character(len=*), intent(in) :: file_hist
         integer, intent(in) :: itime_in_file
         type(block_data_real8_2d), intent(inout) :: sumarea
         logical, intent(inout) :: filter(:)
         integer, intent(in) :: nl_soil
         logical, intent(in) :: forcing_has_missing_value
         logical, intent(in) :: forcmask_pch(:)
      END SUBROUTINE reactive_history_if

      SUBROUTINE reactive_remap_lulcc_if (patchclass_new, eindex_new, patchclass_old, eindex_old, &
         lccpct_patches, old_patch_area)
         USE MOD_Precision
         integer, intent(in) :: patchclass_new(:), patchclass_old(:)
         integer*8, intent(in) :: eindex_new(:), eindex_old(:)
         real(r8), intent(in), optional :: lccpct_patches(:,:)
         real(r8), intent(in), optional :: old_patch_area(:)
      END SUBROUTINE reactive_remap_lulcc_if

      SUBROUTINE reactive_publish_levee_flood_if (fldfrc_patch)
         USE MOD_Precision
         real(r8), intent(in) :: fldfrc_patch(:)
      END SUBROUTINE reactive_publish_levee_flood_if

      SUBROUTINE reactive_publish_flood_if (fldfrc_patch, flddph_patch)
         USE MOD_Precision
         real(r8), intent(in) :: fldfrc_patch(:)
         real(r8), intent(in) :: flddph_patch(:)
      END SUBROUTINE reactive_publish_flood_if
   end interface

   PRIVATE

   integer, parameter :: REACTIVE_CALLBACKS_INITIAL_CAPACITY = 4

   type :: reactive_callbacks_type
      character(len=32)  :: name = ''
      character(len=128) :: aliases = ''
      procedure(reactive_has_if),             pointer, nopass :: has => null()
      procedure(reactive_noarg_if),           pointer, nopass :: refresh => null()
      procedure(reactive_init_if),            pointer, nopass :: init => null()
      procedure(reactive_read_restart_if),    pointer, nopass :: read_restart => null()
      procedure(reactive_write_restart_if),   pointer, nopass :: write_restart => null()
      procedure(reactive_lake_step_if),       pointer, nopass :: lake_step => null()
      procedure(reactive_wetland_decomp_if),  pointer, nopass :: wetland_decomp => null()
      procedure(reactive_soil_step_if),       pointer, nopass :: soil_step => null()
      procedure(reactive_noarg_if),           pointer, nopass :: report => null()
      procedure(reactive_noarg_if),           pointer, nopass :: flush_acc_fluxes => null()
      procedure(reactive_noarg_if),           pointer, nopass :: accumulate_fluxes => null()
      procedure(reactive_history_if),         pointer, nopass :: history => null()
      procedure(reactive_noarg_if),           pointer, nopass :: save_lulcc => null()
      procedure(reactive_remap_lulcc_if),     pointer, nopass :: remap_lulcc => null()
      procedure(reactive_publish_levee_flood_if), pointer, nopass :: publish_levee_flood => null()
      procedure(reactive_publish_flood_if),    pointer, nopass :: publish_flood => null()
      procedure(reactive_noarg_if),           pointer, nopass :: final => null()
   end type reactive_callbacks_type

   type(reactive_callbacks_type), allocatable, save :: reactive_callbacks(:)
   integer, save :: n_reactive_callbacks = 0
   logical, save :: reactive_callbacks_ready = .false.
   logical, save :: reactive_refresh_dirty = .true.

   PUBLIC :: tracer_reactive_has
   PUBLIC :: tracer_reactive_init, tracer_reactive_final
   PUBLIC :: tracer_reactive_resolve_step, tracer_reactive_lake_step
   PUBLIC :: tracer_reactive_wetland_decomp, tracer_reactive_soil_step
   PUBLIC :: tracer_reactive_report
   PUBLIC :: tracer_reactive_write_restart, tracer_reactive_read_restart
   PUBLIC :: tracer_reactive_history
   PUBLIC :: tracer_reactive_flush_acc_fluxes, tracer_reactive_accumulate_fluxes
   PUBLIC :: tracer_reactive_save_lulcc_state, tracer_reactive_remap_lulcc_state
   PUBLIC :: tracer_reactive_publish_levee_flood_patch, tracer_reactive_publish_flood_patch
   PUBLIC :: register_reactive_callbacks

CONTAINS

   SUBROUTINE ensure_reactive_callbacks_registered ()

      IMPLICIT NONE

      IF (reactive_callbacks_ready) RETURN
      reactive_callbacks_ready = .true.

      CALL register_all_reactive_callbacks ()

   END SUBROUTINE ensure_reactive_callbacks_registered

   SUBROUTINE prepare_reactive_dispatch ()

      IMPLICIT NONE

      CALL ensure_reactive_callbacks_registered ()
      IF (reactive_refresh_dirty) CALL refresh_reactive_callback_states_if_dirty ()

   END SUBROUTINE prepare_reactive_dispatch

   SUBROUTINE mark_reactive_callbacks_dirty ()

      IMPLICIT NONE

      reactive_refresh_dirty = .true.

   END SUBROUTINE mark_reactive_callbacks_dirty

   SUBROUTINE register_reactive_callbacks (name, aliases, has_fn, refresh_fn, init_fn, &
      read_restart_fn, write_restart_fn, lake_step_fn, wetland_decomp_fn, &
      soil_step_fn, report_fn, flush_acc_fluxes_fn, accumulate_fluxes_fn, &
      history_fn, save_lulcc_fn, remap_lulcc_fn, publish_levee_flood_fn, &
      publish_flood_fn, final_fn)

      IMPLICIT NONE
      character(len=*), intent(in) :: name
      character(len=*), intent(in), optional :: aliases
      procedure(reactive_has_if),            optional :: has_fn
      procedure(reactive_noarg_if),          optional :: refresh_fn
      procedure(reactive_init_if),           optional :: init_fn
      procedure(reactive_read_restart_if),   optional :: read_restart_fn
      procedure(reactive_write_restart_if),  optional :: write_restart_fn
      procedure(reactive_lake_step_if),      optional :: lake_step_fn
      procedure(reactive_wetland_decomp_if), optional :: wetland_decomp_fn
      procedure(reactive_soil_step_if),      optional :: soil_step_fn
      procedure(reactive_noarg_if),          optional :: report_fn
      procedure(reactive_noarg_if),          optional :: flush_acc_fluxes_fn
      procedure(reactive_noarg_if),          optional :: accumulate_fluxes_fn
      procedure(reactive_history_if),        optional :: history_fn
      procedure(reactive_noarg_if),          optional :: save_lulcc_fn
      procedure(reactive_remap_lulcc_if),    optional :: remap_lulcc_fn
      procedure(reactive_publish_levee_flood_if), optional :: publish_levee_flood_fn
      procedure(reactive_publish_flood_if),  optional :: publish_flood_fn
      procedure(reactive_noarg_if),          optional :: final_fn

      integer :: idx

      IF (len_trim(name) <= 0) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive: cannot register reactive tracer with empty name')
      ENDIF
      IF (.not. present(has_fn)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive: reactive tracer registration requires has_fn')
      ENDIF
      IF (present(aliases)) THEN
         IF (alias_list_has_duplicates(name, aliases)) THEN
            CALL CoLM_stop ('MOD_Tracer_Reactive: duplicate reactive tracer alias in registration')
         ENDIF
      ENDIF
      IF (reactive_registration_conflicts(name, aliases)) THEN
         CALL CoLM_stop ('MOD_Tracer_Reactive: duplicate reactive tracer name or alias registration')
      ENDIF

      idx = find_reactive_callback (name)
      IF (idx <= 0) THEN
         n_reactive_callbacks = n_reactive_callbacks + 1
         CALL ensure_reactive_callback_capacity (n_reactive_callbacks)
         idx = n_reactive_callbacks
      ENDIF

      reactive_callbacks(idx)%name = upper_name(name)
      reactive_callbacks(idx)%aliases = ''
      IF (present(aliases)) reactive_callbacks(idx)%aliases = aliases
      reactive_callbacks(idx)%has => has_fn
      IF (present(refresh_fn))          reactive_callbacks(idx)%refresh => refresh_fn
      IF (present(init_fn))              reactive_callbacks(idx)%init => init_fn
      IF (present(read_restart_fn))      reactive_callbacks(idx)%read_restart => read_restart_fn
      IF (present(write_restart_fn))     reactive_callbacks(idx)%write_restart => write_restart_fn
      IF (present(lake_step_fn))         reactive_callbacks(idx)%lake_step => lake_step_fn
      IF (present(wetland_decomp_fn))    reactive_callbacks(idx)%wetland_decomp => wetland_decomp_fn
      IF (present(soil_step_fn))         reactive_callbacks(idx)%soil_step => soil_step_fn
      IF (present(report_fn))            reactive_callbacks(idx)%report => report_fn
      IF (present(flush_acc_fluxes_fn))  reactive_callbacks(idx)%flush_acc_fluxes => flush_acc_fluxes_fn
      IF (present(accumulate_fluxes_fn)) reactive_callbacks(idx)%accumulate_fluxes => accumulate_fluxes_fn
      IF (present(history_fn))           reactive_callbacks(idx)%history => history_fn
      IF (present(save_lulcc_fn))        reactive_callbacks(idx)%save_lulcc => save_lulcc_fn
      IF (present(remap_lulcc_fn))       reactive_callbacks(idx)%remap_lulcc => remap_lulcc_fn
      IF (present(publish_levee_flood_fn)) reactive_callbacks(idx)%publish_levee_flood => publish_levee_flood_fn
      IF (present(publish_flood_fn))     reactive_callbacks(idx)%publish_flood => publish_flood_fn
      IF (present(final_fn))             reactive_callbacks(idx)%final => final_fn
      CALL mark_reactive_callbacks_dirty ()

   END SUBROUTINE register_reactive_callbacks

   SUBROUTINE ensure_reactive_callback_capacity (required)

      IMPLICIT NONE
      integer, intent(in) :: required

      type(reactive_callbacks_type), allocatable :: grown(:)
      integer :: old_size, new_size

      IF (required <= 0) RETURN
      IF (.not. allocated(reactive_callbacks)) THEN
         allocate (reactive_callbacks(max(REACTIVE_CALLBACKS_INITIAL_CAPACITY, required)))
         RETURN
      ENDIF
      IF (required <= size(reactive_callbacks)) RETURN

      old_size = size(reactive_callbacks)
      new_size = max(required, max(REACTIVE_CALLBACKS_INITIAL_CAPACITY, 2 * old_size))
      allocate (grown(new_size))
      IF (old_size > 0) grown(1:old_size) = reactive_callbacks(1:old_size)
      CALL move_alloc (grown, reactive_callbacks)

   END SUBROUTINE ensure_reactive_callback_capacity

   SUBROUTINE refresh_reactive_callback_states_if_dirty ()

      IMPLICIT NONE

      IF (.not. reactive_refresh_dirty) RETURN
      IF (.not. allocated(tracers)) RETURN
      CALL refresh_reactive_callback_states ()
      reactive_refresh_dirty = .false.

   END SUBROUTINE refresh_reactive_callback_states_if_dirty

   SUBROUTINE refresh_reactive_callback_states ()

      IMPLICIT NONE
      integer :: i

      DO i = 1, n_reactive_callbacks
         IF (associated(reactive_callbacks(i)%refresh)) CALL reactive_callbacks(i)%refresh ()
      ENDDO

   END SUBROUTINE refresh_reactive_callback_states

   logical FUNCTION reactive_registration_conflicts (name, aliases)

      IMPLICIT NONE
      character(len=*), intent(in) :: name
      character(len=*), intent(in), optional :: aliases

      integer :: i, start_pos, end_pos, list_len
      character(len=32) :: item, target

      reactive_registration_conflicts = .false.
      target = upper_name(name)

      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_matches(i, target)) THEN
            reactive_registration_conflicts = .true.
            RETURN
         ENDIF
      ENDDO

      IF (.not. present(aliases)) RETURN

      list_len = len_trim(aliases)
      start_pos = 1
      DO WHILE (start_pos <= list_len)
         end_pos = index(aliases(start_pos:list_len), ',')
         IF (end_pos == 0) THEN
            end_pos = list_len
         ELSE
            end_pos = start_pos + end_pos - 2
         ENDIF
         item = upper_name(adjustl(trim(aliases(start_pos:end_pos))))
         IF (len_trim(item) > 0) THEN
            DO i = 1, n_reactive_callbacks
               IF (reactive_callback_matches(i, item)) THEN
                  reactive_registration_conflicts = .true.
                  RETURN
               ENDIF
            ENDDO
         ENDIF
         start_pos = end_pos + 2
      ENDDO

   END FUNCTION reactive_registration_conflicts

   logical FUNCTION alias_list_has_duplicates (name, aliases)

      IMPLICIT NONE
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: aliases

      integer, parameter :: MAX_ALIASES_PER_REACTIVE = 32
      character(len=32) :: seen(MAX_ALIASES_PER_REACTIVE)
      character(len=32) :: item, target
      integer :: nseen, start_pos, end_pos, list_len, i

      alias_list_has_duplicates = .false.
      seen = ''
      nseen = 0
      target = upper_name(name)
      list_len = len_trim(aliases)
      start_pos = 1

      DO WHILE (start_pos <= list_len)
         end_pos = index(aliases(start_pos:list_len), ',')
         IF (end_pos == 0) THEN
            end_pos = list_len
         ELSE
            end_pos = start_pos + end_pos - 2
         ENDIF

         item = upper_name(adjustl(trim(aliases(start_pos:end_pos))))
         IF (len_trim(item) > 0) THEN
            IF (trim(item) == trim(target)) THEN
               alias_list_has_duplicates = .true.
               RETURN
            ENDIF
            DO i = 1, nseen
               IF (trim(item) == trim(seen(i))) THEN
                  alias_list_has_duplicates = .true.
                  RETURN
               ENDIF
            ENDDO
            IF (nseen >= MAX_ALIASES_PER_REACTIVE) THEN
               CALL CoLM_stop ('MOD_Tracer_Reactive: too many aliases for one reactive tracer')
            ENDIF
            nseen = nseen + 1
            seen(nseen) = item
         ENDIF

         start_pos = end_pos + 2
      ENDDO

   END FUNCTION alias_list_has_duplicates

   integer FUNCTION find_reactive_callback (name)

      IMPLICIT NONE
      character(len=*), intent(in) :: name
      integer :: i
      character(len=32) :: target

      find_reactive_callback = 0
      target = upper_name(name)
      DO i = 1, n_reactive_callbacks
         IF (trim(reactive_callbacks(i)%name) == trim(target)) THEN
            find_reactive_callback = i
            RETURN
         ENDIF
      ENDDO

   END FUNCTION find_reactive_callback

   logical FUNCTION tracer_reactive_has (tracer_name)

      IMPLICIT NONE
      character(len=*), intent(in) :: tracer_name

      integer :: i
      character(len=32) :: target

      CALL prepare_reactive_dispatch ()

      tracer_reactive_has = .false.
      target = upper_name(tracer_name)

      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_matches(i, target)) THEN
            tracer_reactive_has = reactive_callback_active(i, tracer_name)
            RETURN
         ENDIF
      ENDDO

      IF (.not. allocated(tracers)) RETURN
      DO i = 1, ntracers
         IF (upper_name(tracers(i)%name) == target .and. tracer_is_reactive(i)) THEN
            tracer_reactive_has = .true.
            RETURN
         ENDIF
      ENDDO

   END FUNCTION tracer_reactive_has

   SUBROUTINE tracer_reactive_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)

      IMPLICIT NONE
      integer, intent(in) :: numpatch
      integer, intent(in) :: lc_year
      integer, intent(in) :: jdate(3)
      character(len=*), intent(in) :: casename
      character(len=*), intent(in) :: dir_restart
      character(len=*), intent(in) :: dir_landdata

      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%init)) THEN
            CALL reactive_callbacks(i)%init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_init

   SUBROUTINE tracer_reactive_read_restart (file_restart)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%read_restart)) THEN
            CALL reactive_callbacks(i)%read_restart (file_restart)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_read_restart

   SUBROUTINE tracer_reactive_write_restart (file_restart, compress)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: compress
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%write_restart)) THEN
            CALL reactive_callbacks(i)%write_restart (file_restart, compress)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_write_restart

   SUBROUTINE tracer_reactive_resolve_step (istep_in, istep_local)

      IMPLICIT NONE
      integer, intent(in),  optional :: istep_in
      integer, intent(out)           :: istep_local

      IF (present(istep_in)) THEN
         istep_local = istep_in
      ELSE
         istep_local = 1
      ENDIF

   END SUBROUTINE tracer_reactive_resolve_step

   SUBROUTINE tracer_reactive_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim_phy
      integer,  intent(in) :: isub
      integer,  intent(in) :: nsub
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%lake_step)) THEN
            CALL reactive_callbacks(i)%lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_lake_step

   SUBROUTINE tracer_reactive_wetland_decomp (ipatch)

      IMPLICIT NONE
      integer, intent(in) :: ipatch
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%wetland_decomp)) THEN
            CALL reactive_callbacks(i)%wetland_decomp (ipatch)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_wetland_decomp

   SUBROUTINE tracer_reactive_soil_step (istep_local, ipatch, idate, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%soil_step)) THEN
            CALL reactive_callbacks(i)%soil_step (istep_local, ipatch, idate, deltim)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_soil_step

   SUBROUTINE tracer_reactive_report ()

      IMPLICIT NONE
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%report)) THEN
            CALL reactive_callbacks(i)%report ()
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_report

   SUBROUTINE tracer_reactive_flush_acc_fluxes ()

      IMPLICIT NONE
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%flush_acc_fluxes)) THEN
            CALL reactive_callbacks(i)%flush_acc_fluxes ()
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_flush_acc_fluxes

   SUBROUTINE tracer_reactive_accumulate_fluxes ()

      IMPLICIT NONE
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%accumulate_fluxes)) THEN
            CALL reactive_callbacks(i)%accumulate_fluxes ()
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_accumulate_fluxes

   SUBROUTINE tracer_reactive_history (file_hist, itime_in_file, sumarea, filter, &
      nl_soil, forcing_has_missing_value, forcmask_pch)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_hist
      integer, intent(in) :: itime_in_file
      type(block_data_real8_2d), intent(inout) :: sumarea
      logical, intent(inout) :: filter(:)
      integer, intent(in) :: nl_soil
      logical, intent(in) :: forcing_has_missing_value
      logical, intent(in) :: forcmask_pch(:)
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%history)) THEN
            CALL reactive_callbacks(i)%history (file_hist, itime_in_file, sumarea, filter, &
               nl_soil, forcing_has_missing_value, forcmask_pch)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_history

   SUBROUTINE tracer_reactive_save_lulcc_state ()

      IMPLICIT NONE
      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%save_lulcc)) THEN
            CALL reactive_callbacks(i)%save_lulcc ()
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_save_lulcc_state

   SUBROUTINE tracer_reactive_remap_lulcc_state (patchclass_new, eindex_new, patchclass_old, eindex_old, &
      lccpct_patches, old_patch_area)

      IMPLICIT NONE
      integer, intent(in) :: patchclass_new(:), patchclass_old(:)
      integer*8, intent(in) :: eindex_new(:), eindex_old(:)
      real(r8), intent(in), optional :: lccpct_patches(:,:)
      real(r8), intent(in), optional :: old_patch_area(:)

      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%remap_lulcc)) THEN
            CALL reactive_callbacks(i)%remap_lulcc (patchclass_new, eindex_new, patchclass_old, eindex_old, &
               lccpct_patches, old_patch_area)
         ENDIF
      ENDDO
      CALL mark_reactive_callbacks_dirty ()

   END SUBROUTINE tracer_reactive_remap_lulcc_state

   SUBROUTINE tracer_reactive_publish_levee_flood_patch (fldfrc_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)

      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%publish_levee_flood)) THEN
            CALL reactive_callbacks(i)%publish_levee_flood (fldfrc_patch)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_publish_levee_flood_patch

   SUBROUTINE tracer_reactive_publish_flood_patch (fldfrc_patch, flddph_patch)

      IMPLICIT NONE
      real(r8), intent(in) :: fldfrc_patch(:)
      real(r8), intent(in) :: flddph_patch(:)

      integer :: i

      CALL prepare_reactive_dispatch ()
      DO i = 1, n_reactive_callbacks
         IF (reactive_callback_enabled(i) .and. &
             associated(reactive_callbacks(i)%publish_flood)) THEN
            CALL reactive_callbacks(i)%publish_flood (fldfrc_patch, flddph_patch)
         ENDIF
      ENDDO

   END SUBROUTINE tracer_reactive_publish_flood_patch

   SUBROUTINE tracer_reactive_final ()

      IMPLICIT NONE
      integer :: i

      CALL ensure_reactive_callbacks_registered ()
      DO i = 1, n_reactive_callbacks
         IF (associated(reactive_callbacks(i)%final)) CALL reactive_callbacks(i)%final ()
      ENDDO
      CALL clear_reactive_callbacks ()

   END SUBROUTINE tracer_reactive_final

   logical FUNCTION reactive_callback_enabled (idx)

      IMPLICIT NONE
      integer, intent(in) :: idx

      reactive_callback_enabled = .false.
      IF (idx < 1 .or. idx > n_reactive_callbacks) RETURN
      IF (.not. allocated(reactive_callbacks)) RETURN

      reactive_callback_enabled = reactive_callback_active(idx, reactive_callbacks(idx)%name)

   END FUNCTION reactive_callback_enabled

   logical FUNCTION reactive_callback_active (idx, tracer_name)

      IMPLICIT NONE
      integer, intent(in) :: idx
      character(len=*), intent(in) :: tracer_name

      reactive_callback_active = .false.
      IF (idx < 1 .or. idx > n_reactive_callbacks) RETURN
      IF (.not. allocated(reactive_callbacks)) RETURN
      IF (associated(reactive_callbacks(idx)%has)) THEN
         reactive_callback_active = reactive_callbacks(idx)%has (tracer_name)
      ELSE
         CALL CoLM_stop ('MOD_Tracer_Reactive: registered reactive tracer is missing has callback')
      ENDIF

   END FUNCTION reactive_callback_active

   logical FUNCTION reactive_callback_matches (idx, target)

      IMPLICIT NONE
      integer, intent(in) :: idx
      character(len=*), intent(in) :: target

      reactive_callback_matches = .false.
      IF (idx < 1 .or. idx > n_reactive_callbacks) RETURN
      IF (.not. allocated(reactive_callbacks)) RETURN
      IF (trim(upper_name(reactive_callbacks(idx)%name)) == trim(upper_name(target))) THEN
         reactive_callback_matches = .true.
         RETURN
      ENDIF
      reactive_callback_matches = alias_list_contains (reactive_callbacks(idx)%aliases, target)

   END FUNCTION reactive_callback_matches

   logical FUNCTION alias_list_contains (aliases, target)

      IMPLICIT NONE
      character(len=*), intent(in) :: aliases
      character(len=*), intent(in) :: target

      integer :: start_pos, end_pos, list_len
      character(len=32) :: item, target_upper

      alias_list_contains = .false.
      list_len = len_trim(aliases)
      IF (list_len <= 0) RETURN

      target_upper = upper_name(target)
      start_pos = 1
      DO WHILE (start_pos <= list_len)
         end_pos = index(aliases(start_pos:list_len), ',')
         IF (end_pos == 0) THEN
            end_pos = list_len
         ELSE
            end_pos = start_pos + end_pos - 2
         ENDIF
         item = upper_name(adjustl(trim(aliases(start_pos:end_pos))))
         IF (trim(item) == trim(target_upper)) THEN
            alias_list_contains = .true.
            RETURN
         ENDIF
         start_pos = end_pos + 2
      ENDDO

   END FUNCTION alias_list_contains

   SUBROUTINE clear_reactive_callbacks ()

      IMPLICIT NONE
      integer :: i

      IF (allocated(reactive_callbacks)) THEN
         DO i = 1, n_reactive_callbacks
            reactive_callbacks(i)%name = ''
            reactive_callbacks(i)%aliases = ''
            nullify(reactive_callbacks(i)%has)
            nullify(reactive_callbacks(i)%refresh)
            nullify(reactive_callbacks(i)%init)
            nullify(reactive_callbacks(i)%read_restart)
            nullify(reactive_callbacks(i)%write_restart)
            nullify(reactive_callbacks(i)%lake_step)
            nullify(reactive_callbacks(i)%wetland_decomp)
            nullify(reactive_callbacks(i)%soil_step)
            nullify(reactive_callbacks(i)%report)
            nullify(reactive_callbacks(i)%flush_acc_fluxes)
            nullify(reactive_callbacks(i)%accumulate_fluxes)
            nullify(reactive_callbacks(i)%history)
            nullify(reactive_callbacks(i)%save_lulcc)
            nullify(reactive_callbacks(i)%remap_lulcc)
            nullify(reactive_callbacks(i)%publish_levee_flood)
            nullify(reactive_callbacks(i)%publish_flood)
            nullify(reactive_callbacks(i)%final)
         ENDDO
         deallocate (reactive_callbacks)
      ENDIF
      n_reactive_callbacks = 0
      reactive_callbacks_ready = .false.
      reactive_refresh_dirty = .true.

   END SUBROUTINE clear_reactive_callbacks

   character(len=32) FUNCTION upper_name (value)

      IMPLICIT NONE
      character(len=*), intent(in) :: value

      integer :: i, code

      upper_name = adjustl(value)
      DO i = 1, len_trim(upper_name)
         code = iachar(upper_name(i:i))
         IF (code >= iachar('a') .and. code <= iachar('z')) THEN
            upper_name(i:i) = achar(code - iachar('a') + iachar('A'))
         ENDIF
      ENDDO

   END FUNCTION upper_name

END MODULE MOD_Tracer_Reactive
#endif
