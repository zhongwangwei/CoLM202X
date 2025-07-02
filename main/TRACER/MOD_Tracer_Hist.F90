#include <define.h>

MODULE MOD_Tracer_Hist

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!     Write out 3D tracer forcing data to history files.
!     This module uses only the 3D tracer system.
!
!----------------------------------------------------------------------------

   USE MOD_Tracer_Vars_3DAccFluxes
   USE MOD_Vars_Global, only: spval
   USE MOD_NetCDFSerial
   USE MOD_Tracer_Vars_3DForcing, only: tracer_forc, num_tracers_3d, max_vars_3d
   USE MOD_Hist
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

   CALL allocate_tracer_3d_acc_fluxes ()
   CALL FLUSH_tracer_3d_acc_fluxes ()

   TracerHistForm = 'Gridded'
#if (defined UNSTRUCTURED || defined CATCHMENT)
   IF (DEF_HISTORY_IN_VECTOR) THEN
      TracerHistForm = 'Vector'
   ENDIF
#endif
#ifdef SinglePoint
         TracerHistForm = 'Single'
#endif
   
   ! Note: ghist and mp2g_hist are already initialized by main hist_init
   ! We reuse the existing grid configuration for tracer history
   
#ifdef SinglePoint
   IF (TracerHistForm == 'Single') THEN
      CALL hist_single_init  ()
   ENDIF
#endif
   
#ifdef CatchLateralFlow
         CALL hist_basin_init ()
#endif
   
   END SUBROUTINE tracer_hist_init


   SUBROUTINE tracer_hist_final ()

   IMPLICIT NONE
   CALL deallocate_tracer_3d_acc_fluxes ()

#ifdef SinglePoint
      CALL tracer_hist_single_final ()
#endif

#ifdef CatchLateralFlow
      CALL tracer_hist_basin_final ()
#endif


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
   USE MOD_Tracer_Forcing, only: forcmask_pch
   
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
   integer :: days_month(1:12)

   integer :: i, j, itrace, ivar

   character(len=64) :: varname
   character(len=256) :: longname
   character(len=64) :: units


   IF (itstamp <= ptstamp) THEN
       CALL FLUSH_tracer_3d_acc_fluxes ()
       RETURN
    ELSE
       CALL accumulate_tracer_3d_fluxes ()
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

         CALL julian2monthday(idate(1), idate(2), month, day)

         days_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         IF (isleapyear(idate(1))) days_month(2) = 29

         IF ( trim(DEF_HIST_groupby) == 'YEAR' ) THEN
            write(cdate,'(i4.4)') idate(1)
#ifdef SinglePoint
            IF (USE_SITE_HistWriteBack) THEN
               memory_to_disk = isendofyear(idate,deltim) .or. (.not. (itstamp < etstamp))
            ENDIF
#endif
         ELSEIF ( trim(DEF_HIST_groupby) == 'MONTH' ) THEN
            write(cdate,'(i4.4,"-",i2.2)') idate(1), month
#ifdef SinglePoint
            IF (USE_SITE_HistWriteBack) THEN
               memory_to_disk = isendofmonth(idate,deltim) .or. (.not. (itstamp < etstamp))
            ENDIF
#endif
         ELSEIF ( trim(DEF_HIST_groupby) == 'DAY' ) THEN
            write(cdate,'(i4.4,"-",i2.2,"-",i2.2)') idate(1), month, day
#ifdef SinglePoint
            IF (USE_SITE_HistWriteBack) THEN
               memory_to_disk = isendofday(idate,deltim) .or. (.not. (itstamp < etstamp))
            ENDIF
#endif
         ELSE
            write(*,*) 'Warning : Please USE one of DAY/MONTH/YEAR for history group.'
         ENDIF

         file_hist_tracer = trim(dir_hist) // '/' // trim(casename) //'_hist_tracer_'//trim(cdate)//'.nc'

         CALL hist_write_time (file_hist_tracer, tracer_file_last, 'time', idate, itime_in_file)

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

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN
               filter(:) = patchtype < 99
               do i = 1, num_tracers_3d
                  if (DEF_Tracer_Forcings_NL(i)%has_missing_value) then
                     filter = filter .and. forcmask_pch
                  endif
               enddo
               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (TracerHistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         IF (TracerHistForm == 'Gridded') THEN
            IF (trim(file_hist_tracer) /= trim(tracer_file_last)) THEN
               CALL hist_write_var_real8_2d (file_hist_tracer, 'landarea', ghist, -1, sumarea, &
                  compress = 1, longname = 'land area', units = 'km2')
               CALL hist_write_var_real8_2d (file_hist_tracer, 'landfraction', ghist, -1, landfraction, &
                  compress = 1, longname = 'land fraction', units = '-')
            ENDIF
         ENDIF
         
         ! ------------------------------------------------------------------------------------------
         ! Mapping the fluxes and state variables at patch [numpatch] to grid
         ! ------------------------------------------------------------------------------------------
         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               filter(:) = patchtype < 99

               do i = 1, num_tracers_3d
                  if (DEF_Tracer_Forcings_NL(i)%has_missing_value) then
                     filter = filter .and. forcmask_pch
                  endif
               enddo

               filter = filter .and. patchmask
            ENDIF
         ENDIF

         IF (HistForm == 'Gridded') THEN
            CALL mp2g_hist%get_sumarea (sumarea, filter)
         ENDIF

         ! Write 3D tracer forcing data from accumulated arrays
         DO itrace = 1, num_tracers_3d
            DO ivar = 1, max_vars_3d
               CALL tracer_flux_map_and_write_2d( &
               a_tracer_forc_3d(itrace, ivar, :), file_hist_tracer, &
               DEF_Tracer_Forcings_NL(itrace)%vname(ivar), itime_in_file, sumarea, filter, &
               '-', '-')
            ENDDO
         ENDDO

         IF (allocated(filter)) deallocate (filter)

         CALL FLUSH_tracer_3d_acc_fluxes ()
         tracer_file_last = file_hist_tracer

         IF (p_is_worker) THEN
            IF (allocated(vecacc)) deallocate (vecacc)
            IF (allocated(filter)) deallocate (filter)
         ENDIF

      ENDIF

   END SUBROUTINE tracer_hist_out


   SUBROUTINE tracer_flux_map_and_write_2d ( &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

   USE MOD_Block
   USE MOD_Tracer_Vars_3DAccFluxes, only: tracer_3d_nac
   USE MOD_Vars_Global, only: spval
   USE MOD_Hist, only: ghist
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

      ! Use tracer's own accumulation counter
      IF (p_is_worker) THEN
         IF (tracer_3d_nac > 0) THEN
            WHERE (acc_vec /= spval)  acc_vec = acc_vec / tracer_3d_nac
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