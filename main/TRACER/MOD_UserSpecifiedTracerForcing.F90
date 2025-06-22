#include <define.h>

MODULE MOD_UserSpecifiedTracerForcing
   USE MOD_Precision
   USE MOD_Tracer_Namelist_Defs, only: DEF_Tracers, nl_tracer_forcing_type, DEF_Tracer_Forcings_NL, MAX_TRACER_FORCING_VARS
 

   IMPLICIT NONE

   character(len=256) :: dataset

   integer  :: NVAR_tracer_forcing                ! variable number of forcing data
   integer  :: startyr_tracer_forcing               ! start year of forcing data
   integer  :: startmo_tracer_forcing               ! start month of forcing data
   integer  :: endyr_tracer_forcing                 ! END year of forcing data
   integer  :: endmo_tracer_forcing                 ! END month of forcing data
   integer, allocatable :: dtime_tracer_forcing(:)  ! time interval of forcing data
   integer, allocatable :: offset_tracer_forcing(:) ! offset of forcing data

   logical :: leapyear_tracer_forcing               ! leapyear calendar
   logical :: data2d_tracer_forcing                 ! data in 2 dimension (lon, lat)
   logical :: hightdim_tracer_forcing                ! have "z" dimension (height)
   logical :: dim2d_tracer_forcing                  ! lat/lon value in 2 dimension (lon, lat)

   character(len=256) :: latname_tracer_forcing                  ! dimension name of latitude
   character(len=256) :: lonname_tracer_forcing                  ! dimension name of longitude

   character(len=256) :: groupby_tracer_forcing                  ! file grouped by year/month

   character(len=256), allocatable :: fprefix_tracer_forcing(:)  ! file prefix
   character(len=256), allocatable :: vname_tracer_forcing(:)    ! variable name
   character(len=256), allocatable :: timelog_tracer_forcing(:)  ! variable time log info
   character(len=256), allocatable :: tintalgo_tracer_forcing(:) ! interpolation algorithm

   ! ----- public subroutines -----
   PUBLIC :: init_user_specified_tracer_forcing ! initialization of the selected forcing dataset
   PUBLIC :: tracerfilename                 ! identify the forcing file name
   PUBLIC :: tracerpreprocess               ! preprocess the forcing data

CONTAINS

   ! ----------------
   SUBROUTINE init_user_specified_tracer_forcing(i)
      USE MOD_Namelist
      USE MOD_SPMD_Task
      IMPLICIT NONE

      integer, intent(in) :: i

   ! Local variables
      integer :: ivar
      integer :: debug_alloc_status

      IF (p_is_master) THEN
         WRITE(*,*) "    === DEBUG: init_user_specified_tracer_forcing(", i, ") START ==="
         WRITE(*,*) "      Input tracer index: ", i
      ENDIF

      NVAR_tracer_forcing = DEF_Tracer_Forcings_NL(i)%NVAR
      IF (p_is_master) WRITE(*,*) "      NVAR_tracer_forcing = ", NVAR_tracer_forcing

      IF (p_is_master) WRITE(*,*) "      Allocating time-related arrays..."
      IF (allocated(dtime_tracer_forcing )) THEN
         IF (p_is_master) WRITE(*,*) "        Deallocating existing dtime_tracer_forcing"
         deallocate(dtime_tracer_forcing)
      ENDIF
      IF (allocated(offset_tracer_forcing)) THEN
         IF (p_is_master) WRITE(*,*) "        Deallocating existing offset_tracer_forcing"
         deallocate(offset_tracer_forcing)
      ENDIF
      
      allocate (dtime_tracer_forcing  (NVAR_tracer_forcing), stat=debug_alloc_status)
      IF (debug_alloc_status /= 0) THEN
         IF (p_is_master) WRITE(*,*) "        ERROR: Failed to allocate dtime_tracer_forcing, status = ", debug_alloc_status
      ENDIF
      
      allocate (offset_tracer_forcing (NVAR_tracer_forcing), stat=debug_alloc_status)
      IF (debug_alloc_status /= 0) THEN
         IF (p_is_master) WRITE(*,*) "        ERROR: Failed to allocate offset_tracer_forcing, status = ", debug_alloc_status
      ENDIF

      IF (p_is_master) WRITE(*,*) "      Allocating variable-related arrays..."
      IF (allocated(fprefix_tracer_forcing )) THEN
         IF (p_is_master) WRITE(*,*) "        Deallocating existing fprefix_tracer_forcing"
         deallocate(fprefix_tracer_forcing )
      ENDIF
      IF (allocated(vname_tracer_forcing   )) THEN
         IF (p_is_master) WRITE(*,*) "        Deallocating existing vname_tracer_forcing"
         deallocate(vname_tracer_forcing   )
      ENDIF
      IF (allocated(timelog_tracer_forcing )) THEN
         IF (p_is_master) WRITE(*,*) "        Deallocating existing timelog_tracer_forcing"
         deallocate(timelog_tracer_forcing )
      ENDIF
      IF (allocated(tintalgo_tracer_forcing)) THEN
         IF (p_is_master) WRITE(*,*) "        Deallocating existing tintalgo_tracer_forcing"
         deallocate(tintalgo_tracer_forcing)
      ENDIF
      
      allocate (fprefix_tracer_forcing  (NVAR_tracer_forcing), stat=debug_alloc_status)
      IF (debug_alloc_status /= 0) THEN
         IF (p_is_master) WRITE(*,*) "        ERROR: Failed to allocate fprefix_tracer_forcing, status = ", debug_alloc_status
      ENDIF
      
      allocate (vname_tracer_forcing    (NVAR_tracer_forcing), stat=debug_alloc_status)
      IF (debug_alloc_status /= 0) THEN
         IF (p_is_master) WRITE(*,*) "        ERROR: Failed to allocate vname_tracer_forcing, status = ", debug_alloc_status
      ENDIF
      
      allocate (timelog_tracer_forcing  (NVAR_tracer_forcing), stat=debug_alloc_status)
      IF (debug_alloc_status /= 0) THEN
         IF (p_is_master) WRITE(*,*) "        ERROR: Failed to allocate timelog_tracer_forcing, status = ", debug_alloc_status
      ENDIF
      
      allocate (tintalgo_tracer_forcing (NVAR_tracer_forcing), stat=debug_alloc_status)
      IF (debug_alloc_status /= 0) THEN
         IF (p_is_master) WRITE(*,*) "        ERROR: Failed to allocate tintalgo_tracer_forcing, status = ", debug_alloc_status
      ENDIF

      IF (p_is_master) WRITE(*,*) "      Copying configuration from DEF_Tracer_Forcings_NL(", i, ")..."

      startyr_tracer_forcing          = DEF_Tracer_Forcings_NL(i)%startyr          ! start year of forcing data
      startmo_tracer_forcing          = DEF_Tracer_Forcings_NL(i)%startmo          ! start month of forcing data
      endyr_tracer_forcing            = DEF_Tracer_Forcings_NL(i)%endyr            ! end year of forcing data
      endmo_tracer_forcing            = DEF_Tracer_Forcings_NL(i)%endmo            ! end month of forcing data
      dtime_tracer_forcing(1:NVAR_tracer_forcing)  = DEF_Tracer_Forcings_NL(i)%dtime(1:NVAR_tracer_forcing)   ! time interval of forcing data
      offset_tracer_forcing(1:NVAR_tracer_forcing) = DEF_Tracer_Forcings_NL(i)%offset(1:NVAR_tracer_forcing)  ! offset of forcing data

      leapyear_tracer_forcing         = DEF_Tracer_Forcings_NL(i)%leapyear         ! whether leapyear calendar
      data2d_tracer_forcing           = DEF_Tracer_Forcings_NL(i)%data2d           ! whether data in 2 dimension (lon, lat)
      hightdim_tracer_forcing         = DEF_Tracer_Forcings_NL(i)%hightdim         ! whether have "z" dimension (height)
      dim2d_tracer_forcing            = DEF_Tracer_Forcings_NL(i)%dim2d            ! whether lat/lon in 2 dimension (lon,lat)

      latname_tracer_forcing          = DEF_Tracer_Forcings_NL(i)%latname          ! dimension name of latitude
      lonname_tracer_forcing          = DEF_Tracer_Forcings_NL(i)%lonname          ! dimension name of longitude

      groupby_tracer_forcing          = DEF_Tracer_Forcings_NL(i)%groupby          ! file grouped by year/month

      IF (p_is_master) THEN
         WRITE(*,*) "        Time period: ", startyr_tracer_forcing, "/", startmo_tracer_forcing, " to ", endyr_tracer_forcing, "/", endmo_tracer_forcing
         WRITE(*,*) "        Leap year: ", leapyear_tracer_forcing
         WRITE(*,*) "        Data 2D: ", data2d_tracer_forcing
         WRITE(*,*) "        High dim: ", hightdim_tracer_forcing
         WRITE(*,*) "        Dim 2D: ", dim2d_tracer_forcing
         WRITE(*,*) "        Lat name: ", TRIM(latname_tracer_forcing)
         WRITE(*,*) "        Lon name: ", TRIM(lonname_tracer_forcing)
         WRITE(*,*) "        Group by: ", TRIM(groupby_tracer_forcing)
         WRITE(*,*) "        Time intervals: ", (dtime_tracer_forcing(ivar), ivar=1, min(NVAR_tracer_forcing, 5))
         IF (NVAR_tracer_forcing > 5) WRITE(*,*) "          ... and ", NVAR_tracer_forcing-5, " more"
         WRITE(*,*) "        Offsets: ", (offset_tracer_forcing(ivar), ivar=1, min(NVAR_tracer_forcing, 5))
         IF (NVAR_tracer_forcing > 5) WRITE(*,*) "          ... and ", NVAR_tracer_forcing-5, " more"
      ENDIF

      IF (p_is_master) WRITE(*,*) "      Copying variable-specific information..."
      DO ivar = 1, NVAR_tracer_forcing
         fprefix_tracer_forcing (ivar) = DEF_Tracer_Forcings_NL(i)%fprefix(ivar)   ! file prefix
         vname_tracer_forcing   (ivar) = DEF_Tracer_Forcings_NL(i)%vname(ivar)     ! variable name
         timelog_tracer_forcing (ivar) = DEF_Tracer_Forcings_NL(i)%timelog(ivar)   ! variable name
         tintalgo_tracer_forcing(ivar) = DEF_Tracer_Forcings_NL(i)%tintalgo(ivar)  ! interpolation algorithm
         
         IF (p_is_master .AND. ivar <= 3) THEN
            WRITE(*,*) "        Variable ", ivar, ":"
            WRITE(*,*) "          Prefix: ", TRIM(fprefix_tracer_forcing(ivar))
            WRITE(*,*) "          Name: ", TRIM(vname_tracer_forcing(ivar))
            WRITE(*,*) "          Time log: ", TRIM(timelog_tracer_forcing(ivar))
            WRITE(*,*) "          Interpolation: ", TRIM(tintalgo_tracer_forcing(ivar))
         ENDIF
      ENDDO
      IF (p_is_master .AND. NVAR_tracer_forcing > 3) THEN
         WRITE(*,*) "        ... and ", NVAR_tracer_forcing-3, " more variables"
      ENDIF

      IF (p_is_master) THEN
         WRITE(*,*) "    === DEBUG: init_user_specified_tracer_forcing(", i, ") COMPLETED ==="
      ENDIF

   END SUBROUTINE init_user_specified_tracer_forcing

   ! ----------------
   FUNCTION tracerfilename(year, month, day, var_i) RESULT(metfilename)

   USE MOD_Namelist
   USE MOD_SPMD_Task
   IMPLICIT NONE

   integer, intent(in) :: year
   integer, intent(in) :: month
   integer, intent(in) :: day
   integer, intent(in) :: var_i
   character(len=256)  :: metfilename
   character(len=256)  :: yearstr
   character(len=256)  :: monthstr

      IF (p_is_master) THEN
         WRITE(*,*) "      === DEBUG: tracerfilename called ==="
         WRITE(*,*) "        Input: year=", year, " month=", month, " day=", day, " var_i=", var_i
      ENDIF

      write(yearstr, '(I4.4)') year
      write(monthstr, '(I2.2)') month

      IF (p_is_master) THEN
         WRITE(*,*) "        Formatted strings: yearstr='", TRIM(yearstr), "' monthstr='", TRIM(monthstr), "'"
         WRITE(*,*) "        Dataset name: '", TRIM(DEF_Tracer_Forcings_NL(1)%dataset_name), "'"
      ENDIF

      select CASE (trim(DEF_Tracer_Forcings_NL(1)%dataset_name))
      
      CASE ('isogsm')
         !DESCRIPTION
         !===========
            !--- Isotopes-incorporated Global Spectral Model (IsoGSM)
   
         !data source:
         !-------------------
            !---https://isotope.iis.u-tokyo.ac.jp/about-our-lab?lang=en
   
         !References:
         !-------------------
            !---Bong, H., Cauquoin, A., Okazaki, A., Chang, E.-C., Werner, M., Wei, Z., et al. (2024). 
            !   Process-based intercomparison of water isotope-enabled models and reanalysis nudging effects. 
            !   Journal of Geophysical Research: Atmospheres, 129, e2023JD038719. 
            !   https://doi.org/10.1029/2023JD038719
   
         !REVISION HISTORY
         !----------------
            !---2025.03.23   Zhongwang Wei @ SYSU: add the isotope forcing data
   
            metfilename = '/'//trim(fprefix_tracer_forcing(var_i))//'_'//trim(yearstr)//'.nc'
            IF (p_is_master) THEN
               WRITE(*,*) "        CASE IsoGSM selected"
               WRITE(*,*) "        fprefix_tracer_forcing(", var_i, ") = '", TRIM(fprefix_tracer_forcing(var_i)), "'"
               WRITE(*,*) "        Generated filename: '", TRIM(metfilename), "'"
            ENDIF

      
      CASE ('POINT')
         metfilename = '/'//trim(fprefix_tracer_forcing(1))
         IF (p_is_master) THEN
            WRITE(*,*) "        CASE POINT selected"
            WRITE(*,*) "        fprefix_tracer_forcing(1) = '", TRIM(fprefix_tracer_forcing(1)), "'"
            WRITE(*,*) "        Generated filename: '", TRIM(metfilename), "'"
         ENDIF
         
      CASE DEFAULT
         IF (p_is_master) THEN
            WRITE(*,*) "        WARNING: Unknown dataset name '", TRIM(DEF_Tracer_Forcings_NL(1)%dataset_name), "'"
            WRITE(*,*) "        Using default POINT format"
         ENDIF
         metfilename = '/'//trim(fprefix_tracer_forcing(1))
      END select
      
      IF (p_is_master) THEN
         WRITE(*,*) "      === DEBUG: tracerfilename returning '", TRIM(metfilename), "' ==="
      ENDIF
      
      ! IF (DEF_USE_CBL_HEIGHT) THEN
      !    select CASE (var_i)
      !    CASE (9)
      !       metfilename = '/'//trim(fprefix_tracer_forcing(9))//'_'//trim(yearstr)//'_'//trim(monthstr)//&
      !          '_boundary_layer_height.nc4'
      !    END select
      ! ENDIF
   END FUNCTION tracerfilename

 ! preprocess for forcing data [not applicable yet for PRINCETON]
 ! ------------------------------------------------------------
   SUBROUTINE tracerpreprocess(grid, forcn)

   USE MOD_Const_Physical
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Qsadv
   IMPLICIT NONE
   type(grid_type), intent(in) :: grid
   type(block_data_real8_2d), intent(inout) :: forcn(:)

   integer  :: iblkme, ib, jb, i, j
   real(r8) :: es, esdT, qsat_tmp, dqsat_tmpdT, e, ea

      !----------------------------------------------------------------------------
      ! use polynomials to calculate saturation vapor pressure and derivative with
      ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
      ! required to convert relative humidity to specific humidity
      !----------------------------------------------------------------------------
      IF (trim(DEF_Tracer_Forcings_NL(1)%dataset_name) == 'POINT') THEN
#ifdef SinglePoint
         CALL qsadv(forcn(1)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1), &
                    forcn(3)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1), &
                    es,esdT,qsat_tmp,dqsat_tmpdT)
         IF (qsat_tmp < forcn(2)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1)) THEN
            forcn(2)%blk(gblock%xblkme(1),gblock%yblkme(1))%val(1,1) = qsat_tmp
         ENDIF
#endif
      ELSE
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)

            DO j = 1, grid%ycnt(jb)
               DO i = 1, grid%xcnt(ib)

                  select CASE (trim(DEF_Tracer_Forcings_NL(1)%dataset_name))


                  CASE ('IsoGSM') ! IsoGSM forcing
                     CALL qsadv (forcn(1)%blk(ib,jb)%val(i,j), forcn(3)%blk(ib,jb)%val(i,j), &
                        es,esdT,qsat_tmp,dqsat_tmpdT)
                     IF (qsat_tmp < forcn(2)%blk(ib,jb)%val(i,j)) THEN
                        forcn(2)%blk(ib,jb)%val(i,j) = qsat_tmp
                     ENDIF

                  END select

               ENDDO
            ENDDO
         ENDDO
      ENDIF

   END SUBROUTINE tracerpreprocess

END MODULE MOD_UserSpecifiedTracerForcing
! ---------- EOP ------------
