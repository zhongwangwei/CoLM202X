#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_ForcingInput

   ! ------------------------------------------------------------------
   ! Per-species tracer forcing input.
   !
   ! Each tracer's parameter file (the one mapped by DEF_TRACER_PARAM_FILES
   ! and already read for &nl_colm_tracer_parameter) may carry a second,
   ! optional namelist group describing that tracer's OWN forcing inputs:
   !
   !   &nl_colm_tracer_forcing
   !      forcing_num        = 2
   !      forcing_role       = 'precip', 'vapor'
   !      forcing_fprefix    = 'IsoGSM_prate', 'IsoGSM_Q'
   !      forcing_vname      = 'prate1sfc',    'spfh12m'
   !      forcing_tintalgo   = 'nearest',      'linear'
   !      forcing_dtime      = 21600, 21600
   !      forcing_offset     = 10800, 10800
   !      forcing_input_mode = 'normalized_over_total', 'normalized_over_total'
   !   /
   !
   ! `role` is the extension point. Water-isotope species use 'precip' and
   ! 'vapor' (consumed by MOD_Tracer_Forcing). Reactive species such as CH4
   ! reserve their own roles (e.g. 'inundation', 'atm'); those entries are
   ! loaded and stored here, ready for the CH4 module to query once it is
   ! wired to consume them. Unknown roles are kept but ignored by the
   ! water-isotope forcing path.
   !
   ! This module replaces the former global DEF_forcing%tracer_* / legacy
   ! per-species DEF_forcing%precipitation_O18_* namelist path: tracer
   ! forcing now lives entirely in the per-species parameter files and is
   ! read by the TRACER subsystem itself.
   ! ------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Tracer_Defs, only: ntracers, tracers, tracer_param_file_for_index, tracer_lower

   IMPLICIT NONE
   SAVE
   PRIVATE

   integer, parameter, public :: TRACER_FORCING_MAX = 8

   type, public :: tracer_forcing_spec_type
      character(len=32)  :: role       = 'none'
      character(len=256) :: fprefix    = 'null'
      character(len=64)  :: vname      = 'null'
      character(len=32)  :: tintalgo   = 'linear'
      integer            :: dtime      = 21600
      integer            :: offset     = 0
      character(len=32)  :: input_mode = 'normalized_over_total'
   end type tracer_forcing_spec_type

   integer, allocatable, public :: tracer_forcing_n(:)                 ! (ntracers)
   type(tracer_forcing_spec_type), allocatable, public :: &
        tracer_forcing_specs(:,:)                                      ! (TRACER_FORCING_MAX, ntracers)

   PUBLIC :: tracer_forcing_input_load
   PUBLIC :: tracer_forcing_input_final
   PUBLIC :: tracer_forcing_input_count
   PUBLIC :: tracer_forcing_input_get
   PUBLIC :: tracer_forcing_input_find

CONTAINS

   SUBROUTINE tracer_forcing_input_load ()
      USE MOD_SPMD_Task, only: p_is_master, CoLM_stop
      IMPLICIT NONE
      integer :: itrc, k, nf, ierr, unit_nml
      logical :: found, fexists
      character(len=256) :: nlfile
      character(len=512) :: iomsg

      ! namelist scratch (parallel arrays for cross-compiler robustness)
      integer            :: forcing_num
      character(len=32)  :: forcing_role(TRACER_FORCING_MAX)
      character(len=256) :: forcing_fprefix(TRACER_FORCING_MAX)
      character(len=64)  :: forcing_vname(TRACER_FORCING_MAX)
      character(len=32)  :: forcing_tintalgo(TRACER_FORCING_MAX)
      integer            :: forcing_dtime(TRACER_FORCING_MAX)
      integer            :: forcing_offset(TRACER_FORCING_MAX)
      character(len=32)  :: forcing_input_mode(TRACER_FORCING_MAX)
      namelist /nl_colm_tracer_forcing/ forcing_num, forcing_role, forcing_fprefix, &
         forcing_vname, forcing_tintalgo, forcing_dtime, forcing_offset, forcing_input_mode

      IF (ntracers <= 0) RETURN
      IF (.not. allocated(tracer_forcing_n))     allocate(tracer_forcing_n(ntracers))
      IF (.not. allocated(tracer_forcing_specs)) allocate(tracer_forcing_specs(TRACER_FORCING_MAX, ntracers))
      tracer_forcing_n(:) = 0

      DO itrc = 1, ntracers
         CALL tracer_param_file_for_index (itrc, '', nlfile, found)
         IF (.not. found) CYCLE
         INQUIRE (file=trim(nlfile), exist=fexists)
         IF (.not. fexists) CYCLE
         IF (.not. tracer_forcing_group_present(nlfile)) CYCLE

         forcing_num            = 0
         forcing_role(:)        = 'none'
         forcing_fprefix(:)     = 'null'
         forcing_vname(:)       = 'null'
         forcing_tintalgo(:)    = 'linear'
         forcing_dtime(:)       = 21600
         forcing_offset(:)      = 0
         forcing_input_mode(:)  = 'normalized_over_total'

         open (newunit=unit_nml, status='OLD', file=trim(nlfile), form='FORMATTED')
         iomsg = ''
         read (unit_nml, nml=nl_colm_tracer_forcing, iostat=ierr, iomsg=iomsg)
         close (unit_nml)
         IF (ierr /= 0) THEN
            IF (p_is_master) THEN
               write(*,'(3A)') 'ERROR tracer_forcing_input_load: invalid &nl_colm_tracer_forcing in ', &
                  trim(nlfile), ': '
               write(*,'(A)') trim(iomsg)
            ENDIF
            CALL CoLM_stop()
         ENDIF

         IF (forcing_num < 0 .or. forcing_num > TRACER_FORCING_MAX) THEN
            IF (p_is_master) THEN
               write(*,'(A,I0,A,I0,A,A,A)') 'ERROR tracer_forcing_input_load: forcing_num=', forcing_num, &
                  ' must be between 0 and ', TRACER_FORCING_MAX, ' in ', trim(nlfile), '.'
            ENDIF
            CALL CoLM_stop()
         ENDIF
         nf = forcing_num
         tracer_forcing_n(itrc) = nf
         DO k = 1, nf
            IF (forcing_dtime(k) <= 0) THEN
               IF (p_is_master) THEN
                  write(*,'(A,I0,A,A,A,I0,A)') 'ERROR tracer_forcing_input_load: forcing_dtime(', k, &
                     ') for tracer "', trim(tracers(itrc)%name), '" must be > 0, got ', forcing_dtime(k), '.'
               ENDIF
               CALL CoLM_stop()
            ENDIF
            tracer_forcing_specs(k,itrc)%role       = adjustl(forcing_role(k))
            tracer_forcing_specs(k,itrc)%fprefix    = adjustl(forcing_fprefix(k))
            tracer_forcing_specs(k,itrc)%vname      = adjustl(forcing_vname(k))
            tracer_forcing_specs(k,itrc)%tintalgo   = adjustl(forcing_tintalgo(k))
            tracer_forcing_specs(k,itrc)%dtime      = forcing_dtime(k)
            tracer_forcing_specs(k,itrc)%offset     = forcing_offset(k)
            tracer_forcing_specs(k,itrc)%input_mode = adjustl(forcing_input_mode(k))
         ENDDO

         IF (p_is_master) write(*,'(A,I0,A,A,A,I0,A)') &
            ' [tracer-forcing] tracer ', itrc, ' (', trim(tracers(itrc)%name), '): ', nf, ' forcing input(s)'
      ENDDO
   END SUBROUTINE tracer_forcing_input_load

   integer FUNCTION tracer_forcing_input_count (itrc)
      IMPLICIT NONE
      integer, intent(in) :: itrc
      tracer_forcing_input_count = 0
      IF (.not. allocated(tracer_forcing_n)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      tracer_forcing_input_count = tracer_forcing_n(itrc)
   END FUNCTION tracer_forcing_input_count

   FUNCTION tracer_forcing_input_get (itrc, k) RESULT(spec)
      IMPLICIT NONE
      integer, intent(in) :: itrc, k
      type(tracer_forcing_spec_type) :: spec
      ! defaults from the type definition apply if out of range
      IF (.not. allocated(tracer_forcing_specs)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      IF (k < 1 .or. k > TRACER_FORCING_MAX) RETURN
      spec = tracer_forcing_specs(k,itrc)
   END FUNCTION tracer_forcing_input_get

   integer FUNCTION tracer_forcing_input_find (itrc, role)
      IMPLICIT NONE
      integer, intent(in) :: itrc
      character(len=*), intent(in) :: role
      integer :: k
      tracer_forcing_input_find = 0
      IF (.not. allocated(tracer_forcing_n)) RETURN
      IF (itrc < 1 .or. itrc > ntracers) RETURN
      DO k = 1, tracer_forcing_n(itrc)
         IF (tracer_forcing_role_equal(tracer_forcing_specs(k,itrc)%role, role)) THEN
            tracer_forcing_input_find = k
            RETURN
         ENDIF
      ENDDO
   END FUNCTION tracer_forcing_input_find

   logical FUNCTION tracer_forcing_role_equal (a, b)
      IMPLICIT NONE
      character(len=*), intent(in) :: a, b
      tracer_forcing_role_equal = (tracer_lower(a) == tracer_lower(b))
   END FUNCTION tracer_forcing_role_equal

   logical FUNCTION tracer_forcing_group_present (nlfile)
      IMPLICIT NONE
      character(len=*), intent(in) :: nlfile
      character(len=1024) :: line
      integer :: ierr, unit_nml

      tracer_forcing_group_present = .false.
      open(newunit=unit_nml, status='OLD', file=trim(nlfile), form='FORMATTED')
      DO
         read(unit_nml, '(A)', iostat=ierr) line
         IF (ierr /= 0) EXIT
         line = adjustl(line)
         IF (line(1:1) == '!') CYCLE
         IF (index(tracer_lower(line), '&nl_colm_tracer_forcing') == 1 .or. &
             index(tracer_lower(line), '$nl_colm_tracer_forcing') == 1) THEN
            tracer_forcing_group_present = .true.
            EXIT
         ENDIF
      ENDDO
      close(unit_nml)
   END FUNCTION tracer_forcing_group_present

   SUBROUTINE tracer_forcing_input_final ()
      IMPLICIT NONE
      IF (allocated(tracer_forcing_n))     deallocate(tracer_forcing_n)
      IF (allocated(tracer_forcing_specs)) deallocate(tracer_forcing_specs)
   END SUBROUTINE tracer_forcing_input_final

END MODULE MOD_Tracer_ForcingInput
#endif
