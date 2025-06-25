#include <define.h>

MODULE MOD_Tracer_Hist

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!     Write out tracer forcing data to history files.
!     This module mirrors the structure and functionality of MOD_Hist.F90
!
!----------------------------------------------------------------------------

   USE MOD_Tracer_Vars_1DAccFluxes
   USE MOD_Vars_Global, only: spval
   USE MOD_NetCDFSerial

   USE MOD_HistGridded
#if (defined UNSTRUCTURED || defined CATCHMENT)
   USE MOD_HistVector
#endif
#ifdef SinglePoint
   USE MOD_HistSingle
#endif

   PUBLIC :: tracer_hist_init
   PUBLIC :: tracer_hist_out
   PUBLIC :: tracer_hist_final

   character(len=10) :: TracerHistForm ! 'Gridded', 'Vector', 'Single'
   character(len=256) :: tracer_file_last = 'null'

CONTAINS

   SUBROUTINE tracer_hist_init (dir_hist, lulcc_call)

   USE MOD_SPMD_Task
   USE MOD_Namelist

   IMPLICIT NONE

   character(len=*), intent(in) :: dir_hist
   logical, optional, intent(in) :: lulcc_call

   CALL allocate_tracer_acc_fluxes ()
   CALL FLUSH_tracer_acc_fluxes ()

#if (defined UNSTRUCTURED || defined CATCHMENT)
      TracerHistForm = 'Vector'
#elif (defined SinglePoint)
      TracerHistForm = 'Single'
#else
      TracerHistForm = 'Gridded'
#endif
      
      CALL tracer_hist_gridded_init (dir_hist, lulcc_call)

   END SUBROUTINE tracer_hist_init


   SUBROUTINE tracer_hist_final ()

   IMPLICIT NONE

      ! Simplified finalization - no special cleanup needed
      ! tracer_file_last = 'null'

   END SUBROUTINE tracer_hist_final


   SUBROUTINE tracer_hist_out (idate, deltim, itstamp, etstamp, ptstamp, dir_hist, casename)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_TimeManager
   USE MOD_LandPatch
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Tracer_Namelist_Defs, only: DEF_Tracer_Forcings_NL
   USE MOD_Hist, only: hist_write_time, ghist, mp2g_hist
   USE MOD_HistGridded, only: landfraction
   
   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   type(timestamp), intent(in) :: itstamp
   type(timestamp), intent(in) :: etstamp  
   type(timestamp), intent(in) :: ptstamp
   character(len=*), intent(in) :: dir_hist
   character(len=*), intent(in) :: casename

   ! Local variables - mirror main hist module
   logical :: lwrite
   character(len=256) :: file_hist_tracer
   integer :: itime_in_file
   integer :: month, day
   character(len=10) :: cdate

   type(block_data_real8_2d) :: sumarea
   real(r8), allocatable :: vecacc(:)
   logical,  allocatable :: filter(:)
   
   integer :: i, j

      ! Check if tracers are enabled and properly initialized
      IF (.NOT. DEF_USE_Tracer .OR. num_tracers_acc == 0) THEN
         IF (p_is_master) THEN
            WRITE(*,*) 'Tracer history skipped: DEF_USE_Tracer=', DEF_USE_Tracer, ', num_tracers_acc=', num_tracers_acc
         ENDIF
         RETURN
      ENDIF

      ! Mirror the accumulation logic from main hist module
      IF (itstamp <= ptstamp) THEN
         CALL FLUSH_tracer_acc_fluxes ()
         RETURN
      ELSE
         CALL accumulate_tracer_fluxes ()
      ENDIF
      ! Mirror the writing frequency logic
      select CASE (trim(adjustl(DEF_HIST_FREQ)))
      CASE ('TIMESTEP')
         lwrite = .true.
      CASE ('HOURLY')
         lwrite = isendofhour (idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE ('DAILY')
         lwrite = isendofday  (idate, deltim) .or. (.not. (itstamp < etstamp))

      CASE ('MONTHLY')
         lwrite = isendofmonth(idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE ('YEARLY')
         lwrite = isendofyear (idate, deltim) .or. (.not. (itstamp < etstamp))
      CASE default
         lwrite = .false.
         write(*,*) &
         'Warning : Please USE one of TIMESTEP/HOURLY/DAILY/MONTHLY/YEARLY for tracer history frequency.'
         write(*,*) &
         '          Set to FALSE by default.                                                     '
      END select

      IF (lwrite) THEN

         ! Mirror the date formatting logic
         CALL julian2monthday(idate(1), idate(2), month, day)

         IF ( trim(DEF_HIST_groupby) == 'YEAR' ) THEN
            write(cdate,'(i4.4)') idate(1)
         ELSEIF ( trim(DEF_HIST_groupby) == 'MONTH' ) THEN
            write(cdate,'(i4.4,"-",i2.2)') idate(1), month
         ELSEIF ( trim(DEF_HIST_groupby) == 'DAY' ) THEN
            write(cdate,'(i4.4,"-",i2.2,"-",i2.2)') idate(1), month, day
         ELSE
            write(*,*) 'Warning : Please USE one of DAY/MONTH/YEAR for tracer history group.'
         ENDIF

         ! Create tracer history filename - mirror main hist pattern
         file_hist_tracer = trim(dir_hist) // '/' // trim(casename) //'_hist_tracer_'//trim(cdate)//'.nc'

         ! Write time coordinate - mirror main hist module
         CALL hist_write_time (file_hist_tracer, tracer_file_last, 'time', idate, itime_in_file)

         ! Allocate arrays - mirror main hist module
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               allocate (filter (numpatch))
               allocate (vecacc (numpatch))
            ENDIF
         ENDIF

         IF (TracerHistForm == 'Gridded') THEN
            IF (p_is_io) THEN
               CALL allocate_block_data (ghist, sumarea)
            ENDIF
         ENDIF

         ! Set up filter - mirror main hist module
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               filter(:) = patchtype < 99
               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (TracerHistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! Write land area and fraction for first time in file - mirror main hist
         IF (TracerHistForm == 'Gridded') THEN
            IF (trim(file_hist_tracer) /= trim(tracer_file_last)) THEN
               CALL hist_write_var_real8_2d (file_hist_tracer, 'landarea', ghist, -1, sumarea, &
                  compress = 1, longname = 'land area', units = 'km2')
               CALL hist_write_var_real8_2d (file_hist_tracer, 'landfraction', ghist, -1, landfraction, &
                  compress = 1, longname = 'land fraction', units = '-')
            ENDIF
         ENDIF

         ! Write tracer forcing variables - mirror main hist pattern
         CALL tracer_hist_gridded_out (file_hist_tracer, itime_in_file, sumarea, filter)

         ! Update file_last - mirror main hist module
         tracer_file_last = file_hist_tracer

         ! Clean up - mirror main hist module
         CALL FLUSH_tracer_acc_fluxes ()

         IF (p_is_worker) THEN
            IF (allocated(vecacc)) deallocate (vecacc)
            IF (allocated(filter)) deallocate (filter)
         ENDIF

      ENDIF

   END SUBROUTINE tracer_hist_out


   SUBROUTINE tracer_hist_gridded_init (dir_hist, lulcc_call)

   USE MOD_Precision
   USE MOD_SPMD_Task

   IMPLICIT NONE

   character(len=*) , intent(in) :: dir_hist
   logical, optional, intent(in) :: lulcc_call

      IF (p_is_master) THEN
         write(*,'(A)') 'Initializing tracer history (gridded mode)'
      ENDIF

   END SUBROUTINE tracer_hist_gridded_init


   SUBROUTINE tracer_hist_gridded_out (file_hist_tracer, itime_in_file, sumarea, filter)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Tracer_Namelist_Defs, only: DEF_Tracer_Forcings_NL

   IMPLICIT NONE

   character(len=256), intent(in) :: file_hist_tracer
   integer, intent(in) :: itime_in_file
   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

   ! Local variables
   integer :: i, j
   character(len=64) :: varname
   character(len=256) :: longname
   character(len=64) :: units

      ! Write tracer forcing variables
      DO i = 1, num_tracers_acc
         DO j = 1, DEF_Tracer_Forcings_NL(i)%NVAR
            IF (j <= max_vars_acc .AND. trim(DEF_Tracer_Forcings_NL(i)%vname(j)) /= 'NULL') THEN
               
               ! Create variable name
               varname = 'tracer_forc_' // &
                  trim(DEF_Tracer_Forcings_NL(i)%tracer_name) // '_' // &
                  trim(DEF_Tracer_Forcings_NL(i)%vname(j))
               
               ! Set metadata
               longname = trim(DEF_Tracer_Forcings_NL(i)%tracer_name) // ' tracer forcing: ' // &
                  trim(DEF_Tracer_Forcings_NL(i)%vname(j))
               units = 'permil_or_concentration'
               
               CALL write_tracer_history_variable_2d (.true., &
                  a_tracer_forc(i, j, :), file_hist_tracer, varname, itime_in_file, &
                  sumarea, filter, longname, units)
            ENDIF
         ENDDO
      ENDDO

      ! Write specific common variables
      !DO i = 1, num_tracers_acc
         ! Precipitation tracer
      !   varname = 'tracer_prate_' // trim(DEF_Tracer_Forcings_NL(i)%tracer_name)
      !   longname = trim(DEF_Tracer_Forcings_NL(i)%tracer_name) // ' in precipitation'
      !   units = 'permil'
         
      !   CALL write_tracer_history_variable_2d (.true., &
      !      a_tracer_prate(i, :), file_hist_tracer, varname, itime_in_file, &
      !      sumarea, filter, longname, units)
         
         ! Humidity tracer
      !   varname = 'tracer_spfh_' // trim(DEF_Tracer_Forcings_NL(i)%tracer_name)
      !   longname = trim(DEF_Tracer_Forcings_NL(i)%tracer_name) // ' in specific humidity'
      !   units = 'permil'
         
      !   CALL write_tracer_history_variable_2d (.true., &
      !      a_tracer_spfh(i, :), file_hist_tracer, varname, itime_in_file, &
      !      sumarea, filter, longname, units)
      !ENDDO

   END SUBROUTINE tracer_hist_gridded_out


   ! Mirror write_history_variable_2d from main hist module
   SUBROUTINE write_tracer_history_variable_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

   IMPLICIT NONE

   logical, intent(in) :: is_hist
   real(r8), intent(inout) :: acc_vec(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units
   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

      IF (.not. is_hist) RETURN

      select CASE (TracerHistForm)
      CASE ('Gridded')
         CALL tracer_flux_map_and_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_tracer_history_variable_2d


   ! Tracer-specific version of flux_map_and_write_2d that uses tracer_nac instead of nac
   SUBROUTINE tracer_flux_map_and_write_2d ( &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

   USE MOD_Block
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(inout) :: acc_vec(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

   ! Local variables
   type(block_data_real8_2d) :: flux_xy_2d
   integer :: iblkme, xblk, yblk, xloc, yloc
   integer :: compress

      ! Use tracer_nac instead of nac, and add safety check
      IF (p_is_worker) THEN
         IF (tracer_nac > 0.0) THEN
            WHERE (acc_vec /= spval)  acc_vec = acc_vec / tracer_nac
         ELSE
            WHERE (acc_vec /= spval)  acc_vec = 0.0  ! Set to zero if no accumulations
         ENDIF
      ENDIF
      
      IF (p_is_io)      CALL allocate_block_data (ghist, flux_xy_2d)

      CALL mp2g_hist%pset2grid (acc_vec, flux_xy_2d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)

                  IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                     IF (flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           = flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                     ENDIF
                  ELSE
                     flux_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, &
         flux_xy_2d, compress, longname, units)

   END SUBROUTINE tracer_flux_map_and_write_2d


END MODULE MOD_Tracer_Hist 