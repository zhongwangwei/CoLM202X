#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane
!=======================================================================
! CH4 reactive-tracer implementation behind the generic reactive dispatch.
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task, only: p_is_worker, CoLM_stop
   USE MOD_Namelist, only: DEF_USE_METHANE_para, DEF_file_METHANE_para, &
      DEF_file_GIEMS, DEF_wetland_finundation_scheme, DEF_USE_SpatialpH, &
      DEF_METHANE_inundation_mode, DEF_TRACER_REACTIVE_PARAM_FILES, &
      DEF_TRACER_REACTIVE_OPTIONS
   USE MOD_Vars_TimeInvariants, only: patchtype, lake_soilc_srf, patchlatr, patchlonr
   USE MOD_Tracer_Methane_Registry, only: methane_registry_init, igas_ch4
   USE MOD_Tracer_Methane_State,    only: allocate_methane_state, &
      init_methane_wetland_fraction_cache, deallocate_methane_state, &
      read_methane_restart, write_methane_restart, initialize_methane_lake_soilc_from_surface
   USE MOD_Tracer_Methane_AccFlux,  only: allocate_methane_acc_fluxes, &
      deallocate_methane_acc_fluxes, flush_methane_acc_fluxes, accumulate_methane_fluxes
   USE MOD_Tracer_Methane_Microbes, only: allocate_methane_microbes_state, &
      deallocate_methane_microbes_state, read_methane_microbes_restart, &
      write_methane_microbes_restart
   USE MOD_Tracer_Methane_Const,    only: read_methane_namelist, &
      configure_methane_inundation_mode, DEF_METHANE
   USE MOD_Tracer_Methane_GIEMS,    only: allocate_methane_giems, &
      deallocate_methane_giems, read_methane_giems, giems_active
   USE MOD_Tracer_Methane_pH,       only: allocate_methane_ph, &
      deallocate_methane_ph, read_methane_ph_patch
   USE MOD_Tracer_Methane_VegOverride, only: allocate_wetland_aere_overrides, &
      deallocate_wetland_aere_overrides
   USE MOD_Tracer_Reactive_Methane_Impl, only: ch4_impl_lake_step => ch4_reactive_lake_step, &
      ch4_impl_wetland_decomp => ch4_reactive_wetland_decomp, &
      ch4_impl_soil_step => ch4_reactive_soil_step, &
      ch4_impl_report => ch4_reactive_report

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: ch4_reactive_has
   PUBLIC :: ch4_reactive_init, ch4_reactive_final
   PUBLIC :: ch4_reactive_lake_step, ch4_reactive_wetland_decomp
   PUBLIC :: ch4_reactive_soil_step, ch4_reactive_report
   PUBLIC :: ch4_reactive_write_restart, ch4_reactive_read_restart
   PUBLIC :: ch4_reactive_flush_acc_fluxes, ch4_reactive_accumulate_fluxes

CONTAINS

   logical FUNCTION ch4_reactive_has ()

      IMPLICIT NONE

      ch4_reactive_has = (igas_ch4 > 0)

   END FUNCTION ch4_reactive_has

   SUBROUTINE ch4_reactive_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)

      IMPLICIT NONE
      integer, intent(in) :: numpatch
      integer, intent(in) :: lc_year
      integer, intent(in) :: jdate(3)
      character(len=*), intent(in) :: casename
      character(len=*), intent(in) :: dir_restart
      character(len=*), intent(in) :: dir_landdata

      character(len=256) :: file_restart_trc
      character(len=14)  :: cdate_restart
      character(len=256) :: cyear_restart
      character(len=256) :: file_param
      logical :: use_param
      real(r8), allocatable :: giems_dummy_patch(:)

      CALL methane_registry_init()
      IF (.not. ch4_reactive_has()) RETURN

      CALL ch4_reactive_resolve_param_file (use_param, file_param)
      IF (use_param) THEN
         CALL read_methane_namelist (file_param)
      END IF
      CALL ch4_reactive_apply_options ()
      CALL configure_methane_inundation_mode ()
      CALL allocate_methane_state (numpatch)
      CALL init_methane_wetland_fraction_cache (numpatch)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL allocate_methane_microbes_state (numpatch)
      ENDIF
      CALL allocate_methane_acc_fluxes (numpatch)

      write(cyear_restart,'(i4.4)') lc_year
      write(cdate_restart,'(i4.4,"-",i3.3,"-",i5.5)') jdate(1), jdate(2), jdate(3)
      file_restart_trc = trim(dir_restart)//'/'//trim(cdate_restart)//'/'//trim(casename)// &
                         '_restart_'//trim(cdate_restart)//'_lc'//trim(cyear_restart)//'.nc'
      CALL ch4_reactive_read_restart (file_restart_trc)
      CALL initialize_methane_lake_soilc_from_surface (patchtype, lake_soilc_srf, DEF_METHANE%allowlakeprod)

      CALL allocate_methane_giems (numpatch)
      IF (DEF_wetland_finundation_scheme == 5) THEN
         IF (trim(DEF_file_GIEMS) == 'null') THEN
            CALL CoLM_stop (' ***** ERROR: DEF_wetland_finundation_scheme=5 requires DEF_file_GIEMS.')
         ENDIF
         IF (p_is_worker .and. allocated(patchlatr) .and. allocated(patchlonr) .and. &
             size(patchlatr) >= numpatch .and. size(patchlonr) >= numpatch) THEN
            CALL read_methane_giems (DEF_file_GIEMS, patchlatr, patchlonr, numpatch)
         ELSE
            allocate(giems_dummy_patch(0))
            CALL read_methane_giems (DEF_file_GIEMS, giems_dummy_patch, giems_dummy_patch, 0)
            deallocate(giems_dummy_patch)
         ENDIF
         IF (.not. giems_active) THEN
            CALL CoLM_stop (' ***** ERROR: GIEMS file could not be loaded for methane scheme 5.')
         ENDIF
      ENDIF

      CALL allocate_methane_ph (numpatch)
      IF (DEF_USE_SpatialpH) THEN
         CALL read_methane_ph_patch (trim(dir_landdata)//'/soil/'//trim(cyear_restart)// &
            '/methane_ph_patches.nc', numpatch)
      ENDIF

      CALL allocate_wetland_aere_overrides (numpatch)

   END SUBROUTINE ch4_reactive_init


   SUBROUTINE ch4_reactive_resolve_param_file (use_param, file_param)

      IMPLICIT NONE
      logical, intent(out) :: use_param
      character(len=*), intent(out) :: file_param

      character(len=512) :: mapped
      logical :: found

      use_param = DEF_USE_METHANE_para
      file_param = DEF_file_METHANE_para

      CALL ch4_reactive_find_species_entry (DEF_TRACER_REACTIVE_PARAM_FILES, mapped, found)
      IF (found) THEN
         file_param = trim(mapped)
         use_param = len_trim(file_param) > 0 .and. trim(ch4_lower(file_param)) /= 'null'
      ENDIF

   END SUBROUTINE ch4_reactive_resolve_param_file

   SUBROUTINE ch4_reactive_apply_options ()

      IMPLICIT NONE
      character(len=512) :: options
      logical :: found

      CALL ch4_reactive_find_species_entry (DEF_TRACER_REACTIVE_OPTIONS, options, found)
      IF (found) CALL ch4_reactive_apply_option_pairs (options)

   END SUBROUTINE ch4_reactive_apply_options

   SUBROUTINE ch4_reactive_find_species_entry (raw_list, value, found)

      IMPLICIT NONE
      character(len=*), intent(in) :: raw_list
      character(len=*), intent(out) :: value
      logical, intent(out) :: found

      integer :: start_pos, end_pos, list_len, colon_pos
      character(len=512) :: entry, key

      value = ''
      found = .false.
      list_len = len_trim(raw_list)
      IF (list_len <= 0) RETURN
      IF (trim(ch4_lower(raw_list)) == 'null') RETURN

      start_pos = 1
      DO WHILE (start_pos <= list_len)
         end_pos = ch4_next_entry_end(raw_list, start_pos, list_len)
         entry = adjustl(trim(raw_list(start_pos:end_pos)))
         colon_pos = index(entry, ':')
         IF (colon_pos > 0) THEN
            key = adjustl(trim(entry(:colon_pos-1)))
            IF (ch4_reactive_is_alias(key)) THEN
               value = adjustl(trim(entry(colon_pos+1:)))
               found = .true.
               RETURN
            ENDIF
         ENDIF
         start_pos = end_pos + 2
      ENDDO

   END SUBROUTINE ch4_reactive_find_species_entry

   integer FUNCTION ch4_next_entry_end (raw_list, start_pos, list_len)

      IMPLICIT NONE
      character(len=*), intent(in) :: raw_list
      integer, intent(in) :: start_pos
      integer, intent(in) :: list_len

      integer :: i

      ch4_next_entry_end = list_len
      DO i = start_pos, list_len
         IF (raw_list(i:i) == ';') THEN
            ch4_next_entry_end = i - 1
            RETURN
         ENDIF
      ENDDO

   END FUNCTION ch4_next_entry_end

   SUBROUTINE ch4_reactive_apply_option_pairs (options)

      IMPLICIT NONE
      character(len=*), intent(in) :: options

      integer :: start_pos, end_pos, opt_len, eq_pos
      character(len=256) :: pair, key, value

      opt_len = len_trim(options)
      start_pos = 1
      DO WHILE (start_pos <= opt_len)
         end_pos = ch4_next_option_end(options, start_pos, opt_len)
         pair = adjustl(trim(options(start_pos:end_pos)))
         eq_pos = index(pair, '=')
         IF (eq_pos > 0) THEN
            key = ch4_lower(adjustl(trim(pair(:eq_pos-1))))
            value = adjustl(trim(pair(eq_pos+1:)))
            SELECT CASE (trim(key))
            CASE ('inundation_mode', 'methane_inundation_mode', 'ch4_inundation_mode')
               IF (len_trim(value) > 0) THEN
                  DEF_METHANE_inundation_mode = value(:min(len(DEF_METHANE_inundation_mode), len_trim(value)))
               ENDIF
            END SELECT
         ENDIF
         start_pos = end_pos + 2
      ENDDO

   END SUBROUTINE ch4_reactive_apply_option_pairs

   integer FUNCTION ch4_next_option_end (options, start_pos, opt_len)

      IMPLICIT NONE
      character(len=*), intent(in) :: options
      integer, intent(in) :: start_pos
      integer, intent(in) :: opt_len

      integer :: i

      ch4_next_option_end = opt_len
      DO i = start_pos, opt_len
         IF (options(i:i) == ',') THEN
            ch4_next_option_end = i - 1
            RETURN
         ENDIF
      ENDDO

   END FUNCTION ch4_next_option_end

   logical FUNCTION ch4_reactive_is_alias (raw_name)

      IMPLICIT NONE
      character(len=*), intent(in) :: raw_name

      character(len=32) :: name

      name = ch4_upper(raw_name)
      ch4_reactive_is_alias = trim(name) == 'CH4' .or. trim(name) == 'METHANE'

   END FUNCTION ch4_reactive_is_alias

   FUNCTION ch4_upper (raw) RESULT(out)

      IMPLICIT NONE
      character(len=*), intent(in) :: raw
      character(len=max(1,len_trim(raw))) :: out
      integer :: i, ia

      out = adjustl(trim(raw))
      DO i = 1, len_trim(out)
         ia = iachar(out(i:i))
         IF (ia >= iachar('a') .and. ia <= iachar('z')) THEN
            out(i:i) = achar(ia - iachar('a') + iachar('A'))
         ENDIF
      ENDDO

   END FUNCTION ch4_upper

   FUNCTION ch4_lower (raw) RESULT(out)

      IMPLICIT NONE
      character(len=*), intent(in) :: raw
      character(len=max(1,len_trim(raw))) :: out
      integer :: i, ia

      out = adjustl(trim(raw))
      DO i = 1, len_trim(out)
         ia = iachar(out(i:i))
         IF (ia >= iachar('A') .and. ia <= iachar('Z')) THEN
            out(i:i) = achar(ia - iachar('A') + iachar('a'))
         ENDIF
      ENDDO

   END FUNCTION ch4_lower

   SUBROUTINE ch4_reactive_read_restart (file_restart)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart

      IF (.not. ch4_reactive_has()) RETURN
      CALL read_methane_restart (file_restart)
      IF (DEF_METHANE%use_microbial_pools) THEN
         CALL read_methane_microbes_restart (file_restart)
      ENDIF

   END SUBROUTINE ch4_reactive_read_restart

   SUBROUTINE ch4_reactive_write_restart (file_restart, compress)

      IMPLICIT NONE
      character(len=*), intent(in) :: file_restart
      integer, intent(in) :: compress

      IF (.not. ch4_reactive_has()) RETURN
      CALL write_methane_restart(file_restart, compress)
      CALL write_methane_microbes_restart(file_restart, compress)

   END SUBROUTINE ch4_reactive_write_restart

   SUBROUTINE ch4_reactive_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim_phy
      integer,  intent(in) :: isub
      integer,  intent(in) :: nsub

      IF (ch4_reactive_has()) THEN
         CALL ch4_impl_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)
      ENDIF

   END SUBROUTINE ch4_reactive_lake_step

   SUBROUTINE ch4_reactive_wetland_decomp (ipatch)

      IMPLICIT NONE
      integer, intent(in) :: ipatch

      IF (ch4_reactive_has()) CALL ch4_impl_wetland_decomp (ipatch)

   END SUBROUTINE ch4_reactive_wetland_decomp

   SUBROUTINE ch4_reactive_soil_step (istep_local, ipatch, idate, deltim)

      IMPLICIT NONE
      integer,  intent(in) :: istep_local
      integer,  intent(in) :: ipatch
      integer,  intent(in) :: idate(3)
      real(r8), intent(in) :: deltim

      IF (ch4_reactive_has()) CALL ch4_impl_soil_step (istep_local, ipatch, idate, deltim)

   END SUBROUTINE ch4_reactive_soil_step

   SUBROUTINE ch4_reactive_report ()

      IMPLICIT NONE

      IF (ch4_reactive_has()) CALL ch4_impl_report ()

   END SUBROUTINE ch4_reactive_report

   SUBROUTINE ch4_reactive_flush_acc_fluxes ()

      IMPLICIT NONE

      IF (ch4_reactive_has()) CALL flush_methane_acc_fluxes ()

   END SUBROUTINE ch4_reactive_flush_acc_fluxes

   SUBROUTINE ch4_reactive_accumulate_fluxes ()

      IMPLICIT NONE

      IF (ch4_reactive_has()) CALL accumulate_methane_fluxes ()

   END SUBROUTINE ch4_reactive_accumulate_fluxes

   SUBROUTINE ch4_reactive_final ()

      IMPLICIT NONE

      IF (.not. ch4_reactive_has()) RETURN

      CALL deallocate_methane_acc_fluxes ()
      CALL deallocate_wetland_aere_overrides ()
      CALL deallocate_methane_ph ()
      CALL deallocate_methane_giems ()
      CALL deallocate_methane_microbes_state ()
      CALL deallocate_methane_state ()

   END SUBROUTINE ch4_reactive_final


END MODULE MOD_Tracer_Reactive_Methane

logical FUNCTION ch4_reactive_has ()

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_has_mod => ch4_reactive_has

   IMPLICIT NONE

   ch4_reactive_has = ch4_reactive_has_mod()

END FUNCTION ch4_reactive_has

SUBROUTINE ch4_reactive_init (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_init_mod => ch4_reactive_init

   IMPLICIT NONE
   integer, intent(in) :: numpatch
   integer, intent(in) :: lc_year
   integer, intent(in) :: jdate(3)
   character(len=*), intent(in) :: casename
   character(len=*), intent(in) :: dir_restart
   character(len=*), intent(in) :: dir_landdata

   CALL ch4_reactive_init_mod (numpatch, lc_year, jdate, casename, dir_restart, dir_landdata)

END SUBROUTINE ch4_reactive_init

SUBROUTINE ch4_reactive_read_restart (file_restart)

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_read_restart_mod => ch4_reactive_read_restart

   IMPLICIT NONE
   character(len=*), intent(in) :: file_restart

   CALL ch4_reactive_read_restart_mod (file_restart)

END SUBROUTINE ch4_reactive_read_restart

SUBROUTINE ch4_reactive_write_restart (file_restart, compress)

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_write_restart_mod => ch4_reactive_write_restart

   IMPLICIT NONE
   character(len=*), intent(in) :: file_restart
   integer, intent(in) :: compress

   CALL ch4_reactive_write_restart_mod (file_restart, compress)

END SUBROUTINE ch4_reactive_write_restart

SUBROUTINE ch4_reactive_lake_step (istep_local, ipatch, idate, deltim_phy, isub, nsub)

   USE MOD_Precision
   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_lake_step_mod => ch4_reactive_lake_step

   IMPLICIT NONE
   integer,  intent(in) :: istep_local
   integer,  intent(in) :: ipatch
   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim_phy
   integer,  intent(in) :: isub
   integer,  intent(in) :: nsub

   CALL ch4_reactive_lake_step_mod (istep_local, ipatch, idate, deltim_phy, isub, nsub)

END SUBROUTINE ch4_reactive_lake_step

SUBROUTINE ch4_reactive_wetland_decomp (ipatch)

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_wetland_decomp_mod => ch4_reactive_wetland_decomp

   IMPLICIT NONE
   integer, intent(in) :: ipatch

   CALL ch4_reactive_wetland_decomp_mod (ipatch)

END SUBROUTINE ch4_reactive_wetland_decomp

SUBROUTINE ch4_reactive_soil_step (istep_local, ipatch, idate, deltim)

   USE MOD_Precision
   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_soil_step_mod => ch4_reactive_soil_step

   IMPLICIT NONE
   integer,  intent(in) :: istep_local
   integer,  intent(in) :: ipatch
   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim

   CALL ch4_reactive_soil_step_mod (istep_local, ipatch, idate, deltim)

END SUBROUTINE ch4_reactive_soil_step

SUBROUTINE ch4_reactive_report ()

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_report_mod => ch4_reactive_report

   IMPLICIT NONE

   CALL ch4_reactive_report_mod ()

END SUBROUTINE ch4_reactive_report

SUBROUTINE ch4_reactive_flush_acc_fluxes ()

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_flush_acc_fluxes_mod => ch4_reactive_flush_acc_fluxes

   IMPLICIT NONE

   CALL ch4_reactive_flush_acc_fluxes_mod ()

END SUBROUTINE ch4_reactive_flush_acc_fluxes

SUBROUTINE ch4_reactive_accumulate_fluxes ()

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_accumulate_fluxes_mod => ch4_reactive_accumulate_fluxes

   IMPLICIT NONE

   CALL ch4_reactive_accumulate_fluxes_mod ()

END SUBROUTINE ch4_reactive_accumulate_fluxes

SUBROUTINE ch4_reactive_final ()

   USE MOD_Tracer_Reactive_Methane, only: ch4_reactive_final_mod => ch4_reactive_final

   IMPLICIT NONE

   CALL ch4_reactive_final_mod ()

END SUBROUTINE ch4_reactive_final
#endif

