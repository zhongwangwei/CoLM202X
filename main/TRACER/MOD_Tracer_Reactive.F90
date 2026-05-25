#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Reactive
!=======================================================================
! Reactive-tracer dispatch layer.
!
! The generic TRACER lifecycle calls this module. Species-specific
! reactive implementations stay below this layer and are activated only
! when their tracer is registered as reactive.
!=======================================================================

   USE MOD_Precision
   USE MOD_DataType, only: block_data_real8_2d
   USE MOD_Tracer_Defs, only: ntracers, tracers, tracer_is_reactive
   IMPLICIT NONE
#ifdef BGC
   logical, external :: ch4_reactive_has
   external :: ch4_reactive_init, ch4_reactive_final
   external :: ch4_reactive_lake_step, ch4_reactive_wetland_decomp
   external :: ch4_reactive_soil_step, ch4_reactive_report
   external :: ch4_reactive_write_restart, ch4_reactive_read_restart
   external :: ch4_reactive_flush_acc_fluxes, ch4_reactive_accumulate_fluxes
   interface
      SUBROUTINE methane_reactive_history (file_hist, itime_in_file, sumarea, filter, &
         nl_soil, forcing_has_missing_value, forcmask_pch)
         USE MOD_DataType, only: block_data_real8_2d
         character(len=*), intent(in) :: file_hist
         integer, intent(in) :: itime_in_file
         type(block_data_real8_2d), intent(inout) :: sumarea
         logical, intent(inout) :: filter(:)
         integer, intent(in) :: nl_soil
         logical, intent(in) :: forcing_has_missing_value
         logical, intent(in) :: forcmask_pch(:)
      END SUBROUTINE methane_reactive_history
   end interface
#endif
   PRIVATE

   PUBLIC :: tracer_reactive_has
   PUBLIC :: tracer_reactive_init, tracer_reactive_final
   PUBLIC :: tracer_reactive_resolve_step, tracer_reactive_lake_step
   PUBLIC :: tracer_reactive_wetland_decomp, tracer_reactive_soil_step
   PUBLIC :: tracer_reactive_report
   PUBLIC :: tracer_reactive_write_restart, tracer_reactive_read_restart
   PUBLIC :: tracer_reactive_history
   PUBLIC :: tracer_reactive_flush_acc_fluxes, tracer_reactive_accumulate_fluxes

CONTAINS

   logical FUNCTION tracer_reactive_has (tracer_name)

      IMPLICIT NONE
      character(len=*), intent(in) :: tracer_name

      integer :: i
      character(len=32) :: target

      tracer_reactive_has = .false.
      target = upper_name(tracer_name)

#ifdef BGC
      IF (target == 'CH4' .or. target == 'METHANE') THEN
         tracer_reactive_has = ch4_reactive_has()
         RETURN
      ENDIF
#endif

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

#ifdef BGC
      CALL ch4_reactive_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)
#endif

   END SUBROUTINE tracer_reactive_init

   SUBROUTINE tracer_reactive_read_restart (file_restart)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart

#ifdef BGC
      CALL ch4_reactive_read_restart (file_restart)
#endif

   END SUBROUTINE tracer_reactive_read_restart

   SUBROUTINE tracer_reactive_write_restart (file_restart, compress)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: compress

#ifdef BGC
      CALL ch4_reactive_write_restart (file_restart, compress)
#endif

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

#ifdef BGC
      CALL ch4_reactive_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)
#endif

   END SUBROUTINE tracer_reactive_lake_step

   SUBROUTINE tracer_reactive_wetland_decomp (ipatch)

      IMPLICIT NONE
      integer, intent(in) :: ipatch

#ifdef BGC
      CALL ch4_reactive_wetland_decomp (ipatch)
#endif

   END SUBROUTINE tracer_reactive_wetland_decomp

   SUBROUTINE tracer_reactive_soil_step (istep_local, ipatch, idate, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

#ifdef BGC
      CALL ch4_reactive_soil_step (istep_local, ipatch, idate, deltim)
#endif

   END SUBROUTINE tracer_reactive_soil_step

   SUBROUTINE tracer_reactive_report ()

      IMPLICIT NONE

#ifdef BGC
      CALL ch4_reactive_report ()
#endif

   END SUBROUTINE tracer_reactive_report

   SUBROUTINE tracer_reactive_flush_acc_fluxes ()

      IMPLICIT NONE

#ifdef BGC
      CALL ch4_reactive_flush_acc_fluxes ()
#endif

   END SUBROUTINE tracer_reactive_flush_acc_fluxes

   SUBROUTINE tracer_reactive_accumulate_fluxes ()

      IMPLICIT NONE

#ifdef BGC
      CALL ch4_reactive_accumulate_fluxes ()
#endif

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

#ifdef BGC
      CALL methane_reactive_history (file_hist, itime_in_file, sumarea, filter, &
         nl_soil, forcing_has_missing_value, forcmask_pch)
#endif

   END SUBROUTINE tracer_reactive_history

   SUBROUTINE tracer_reactive_final ()

      IMPLICIT NONE

#ifdef BGC
      CALL ch4_reactive_final ()
#endif

   END SUBROUTINE tracer_reactive_final

   character(len=32) FUNCTION upper_name (value)

      IMPLICIT NONE
      character(len=*), intent(in) :: value

      integer :: i, code

      upper_name = adjustl(value)
      DO i = 1, len(upper_name)
         code = iachar(upper_name(i:i))
         IF (code >= iachar('a') .and. code <= iachar('z')) THEN
            upper_name(i:i) = achar(code - 32)
         ENDIF
      ENDDO
      upper_name = trim(upper_name)

   END FUNCTION upper_name

END MODULE MOD_Tracer_Reactive

SUBROUTINE tracer_reactive_write_restart (file_restart, compress)

   USE MOD_Tracer_Reactive, only: tracer_reactive_write_restart_mod => tracer_reactive_write_restart

   IMPLICIT NONE
   character(len=*), intent(in) :: file_restart
   integer, intent(in) :: compress

   CALL tracer_reactive_write_restart_mod(file_restart, compress)

END SUBROUTINE tracer_reactive_write_restart
#endif
