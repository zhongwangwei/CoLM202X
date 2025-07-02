#include <define.h>

!-----------------------------------------------------------------------
! Created by Yongjiu Dai, 03/2014
!-----------------------------------------------------------------------

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
MODULE MOD_Tracer_Vars_PFTimeVariables
!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Define PFT time variables
!
!  Added by Hua Yuan, 08/2019
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_TimeManager
#ifdef BGC
   USE MOD_BGC_Vars_PFTimeVariables
#endif

   IMPLICIT NONE
   SAVE
!-----------------------------------------------------------------------
! Time-varying state variables which required by restart run

   ! for LULC_IGBP_PFT or LULC_IGBP_PC

   real(r8), allocatable :: ldew_p_O18       (:) !depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_p_O18  (:) !depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_p_O18  (:) !depth of snow on foliage [mm]
   real(r8), allocatable :: qref_p_O18       (:) !2 m height air specific humidity
   real(r8), allocatable :: vegwp_p_O18    (:,:) !vegetation water potential [mm]


   real(r8), allocatable :: ldew_p_H2       (:) !depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_p_H2  (:) !depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_p_H2  (:) !depth of snow on foliage [mm]
   real(r8), allocatable :: qref_p_H2       (:) !2 m height air specific humidity
   real(r8), allocatable :: vegwp_p_H2    (:,:) !vegetation water potential [mm]

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_PFTimeVariables
   PUBLIC :: Tracer_deallocate_PFTimeVariables
   PUBLIC :: Tracer_READ_PFTimeVariables
   PUBLIC :: Tracer_WRITE_PFTimeVariables
#ifdef RangeCheck
   PUBLIC :: Tracer_check_PFTimeVariables
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Tracer_allocate_PFTimeVariables ()
   !--------------------------------------------------------------------
   ! Allocates memory for CoLM 1d [numpft] variables
   !--------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_LandPFT
   USE MOD_Vars_Global
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN

            allocate (ldew_p_O18       (numpft)) ; ldew_p_O18       (:) = spval !depth of water on foliage [mm]
            allocate (ldew_rain_p_O18  (numpft)) ; ldew_rain_p_O18  (:) = spval !depth of rain on foliage [mm]
            allocate (ldew_snow_p_O18  (numpft)) ; ldew_snow_p_O18  (:) = spval !depth of snow on foliage [mm]
            allocate (ldew_p_H2        (numpft)) ; ldew_p_H2        (:) = spval !depth of water on foliage [mm]
            allocate (ldew_rain_p_H2   (numpft)) ; ldew_rain_p_H2   (:) = spval !depth of rain on foliage [mm]
            allocate (ldew_snow_p_H2   (numpft)) ; ldew_snow_p_H2   (:) = spval !depth of snow on foliage [mm]

            allocate (qref_p_O18       (numpft)) ; qref_p_O18       (:) = spval !2 m height air specific humidity
            allocate (qref_p_H2       (numpft)) ; qref_p_H2       (:) = spval !2 m height air specific humidity

         ENDIF
      ENDIF



   END SUBROUTINE Tracer_allocate_PFTimeVariables

   SUBROUTINE Tracer_READ_PFTimeVariables (file_restart)

   USE MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, DEF_USE_IRRIGATION
   USE MOD_NetCDFVector
   USE MOD_LandPFT
   USE MOD_Vars_Global

   IMPLICIT NONE

      character(len=*), intent(in) :: file_restart


      CALL ncio_read_vector (file_restart, 'ldew_p_O18',landpft, ldew_p_O18 )
      CALL ncio_read_vector (file_restart, 'ldew_rain_p_O18',landpft, ldew_rain_p_O18 )
      CALL ncio_read_vector (file_restart, 'ldew_snow_p_O18',landpft, ldew_snow_p_O18 )
      CALL ncio_read_vector (file_restart, 'ldew_p_H2',landpft, ldew_p_H2 )
      CALL ncio_read_vector (file_restart, 'ldew_rain_p_H2',landpft, ldew_rain_p_H2 )
      CALL ncio_read_vector (file_restart, 'ldew_snow_p_H2',landpft, ldew_snow_p_H2 )

      CALL ncio_read_vector (file_restart, 'qref_p_O18',  landpft, qref_p_O18      )
      CALL ncio_read_vector (file_restart, 'qref_p_H2',  landpft, qref_p_H2      )

      CALL ncio_read_vector (file_restart, 'vegwp_p_O18  ',  nvegwcs, landpft, vegwp_p_O18 )
      CALL ncio_read_vector (file_restart, 'vegwp_p_H2  ',  nvegwcs, landpft, vegwp_p_H2 )


   END SUBROUTINE Tracer_READ_PFTimeVariables

   SUBROUTINE Tracer_WRITE_PFTimeVariables (file_restart)

   USE MOD_Namelist, only: DEF_REST_CompressLevel, DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, &
                           DEF_USE_IRRIGATION
   USE MOD_LandPFT
   USE MOD_NetCDFVector
   USE MOD_Vars_Global
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   ! Local variables
   integer :: compress

      compress = DEF_REST_CompressLevel

      CALL ncio_create_file_vector (file_restart, landpft)
      CALL ncio_define_dimension_vector (file_restart, landpft, 'pft')
      CALL ncio_define_dimension_vector (file_restart, landpft, 'band', 2)
      CALL ncio_define_dimension_vector (file_restart, landpft, 'rtyp', 2)
IF(DEF_USE_PLANTHYDRAULICS)THEN
      CALL ncio_define_dimension_vector (file_restart, landpft, 'vegnodes', nvegwcs)
ENDIF

      CALL ncio_write_vector (file_restart, 'ldew_p_O18       ', 'pft', landpft, ldew_p_O18       , compress)
      CALL ncio_write_vector (file_restart, 'ldew_rain_p_O18  ', 'pft', landpft, ldew_rain_p_O18  , compress)
      CALL ncio_write_vector (file_restart, 'ldew_snow_p_O18  ', 'pft', landpft, ldew_snow_p_O18  , compress)
      CALL ncio_write_vector (file_restart, 'ldew_p_H2        ', 'pft', landpft, ldew_p_H2        , compress)
      CALL ncio_write_vector (file_restart, 'ldew_rain_p_H2   ', 'pft', landpft, ldew_rain_p_H2   , compress)
      CALL ncio_write_vector (file_restart, 'ldew_snow_p_H2   ', 'pft', landpft, ldew_snow_p_H2   , compress)

      CALL ncio_write_vector (file_restart, 'qref_p_O18', 'pft', landpft, qref_p_O18, compress)
      CALL ncio_write_vector (file_restart, 'qref_p_H2', 'pft', landpft, qref_p_H2, compress)

      CALL ncio_write_vector (file_restart, 'vegwp_p_O18  ', 'vegnodes', nvegwcs,  'pft', landpft, vegwp_p_O18, compress)     
      CALL ncio_write_vector (file_restart, 'vegwp_p_H2  ', 'vegnodes', nvegwcs,  'pft', landpft, vegwp_p_H2, compress)

   END SUBROUTINE Tracer_WRITE_PFTimeVariables


   SUBROUTINE Tracer_deallocate_PFTimeVariables
   !--------------------------------------------------------------------
   ! Deallocates memory for CoLM 1d [numpft/numpc] variables
   !--------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_LandPFT

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN

            deallocate (ldew_p_O18       )  ! depth of water on foliage [mm]
            deallocate (ldew_rain_p_O18  )  ! depth of rain on foliage [mm]
            deallocate (ldew_snow_p_O18  )  ! depth of snow on foliage [mm]
            deallocate (ldew_p_H2        )  ! depth of water on foliage [mm]
            deallocate (ldew_rain_p_H2   )  ! depth of rain on foliage [mm]
            deallocate (ldew_snow_p_H2   )  ! depth of snow on foliage [mm]

            deallocate (qref_p_O18       )  ! 2 m height air specific humidity
            deallocate (qref_p_H2       )  ! 2 m height air specific humidity

            deallocate (vegwp_p_O18    )  ! vegetation water potential [mm]
            deallocate (vegwp_p_H2    )  ! vegetation water potential [mm]
! Ozone Stress variables
         ENDIF
      ENDIF


   END SUBROUTINE Tracer_deallocate_PFTimeVariables

#ifdef RangeCheck
   SUBROUTINE Tracer_check_PFTimeVariables

   USE MOD_RangeCheck
   USE MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, DEF_USE_IRRIGATION

   IMPLICIT NONE


      CALL check_vector_data ('        ldew_p_O18', ldew_p_O18       )
      CALL check_vector_data ('   ldew_rain_p_O18', ldew_rain_p_O18  )
      CALL check_vector_data ('   ldew_snow_p_O18', ldew_snow_p_O18  )
      CALL check_vector_data ('        ldew_p_H2', ldew_p_H2        )
      CALL check_vector_data ('   ldew_rain_p_H2', ldew_rain_p_H2   )
      CALL check_vector_data ('   ldew_snow_p_H2', ldew_snow_p_H2   )

      CALL check_vector_data ('        qref_p_O18', qref_p_O18       )
      CALL check_vector_data ('        qref_p_H2', qref_p_H2       )

   END SUBROUTINE Tracer_check_PFTimeVariables
#endif

END MODULE MOD_Tracer_Vars_PFTimeVariables
#endif


MODULE MOD_Tracer_Vars_TimeVariables
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

   USE MOD_Precision
   USE MOD_TimeManager

   IMPLICIT NONE
   SAVE
! -----------------------------------------------------------------
! Time-varying state variables which required by restart run

   real(r8), allocatable :: wliq_soisno_O18 (:,:) ! liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_soisno_O18 (:,:) ! ice lens in layers [kg/m2]
   real(r8), allocatable :: wliq_soisno_H2 (:,:) ! liquid water in layers [kg/m2]
   real(r8), allocatable :: wice_soisno_H2 (:,:) ! ice lens in layers [kg/m2]

   real(r8), allocatable :: rootr_O18   (:,:) ! transpiration contribution fraction from different layers
   real(r8), allocatable :: rootflux_O18 (:,:) ! water exchange between soil and root. Positive: soil->root [?]
   real(r8), allocatable :: rootr_H2   (:,:) ! transpiration contribution fraction from different layers
   real(r8), allocatable :: rootflux_H2 (:,:) ! water exchange between soil and root. Positive: soil->root [?]

   real(r8), allocatable :: ldew_O18       (:) ! depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_O18  (:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_O18  (:) ! depth of snow on foliage [mm]
   real(r8), allocatable :: ldew_H2        (:) ! depth of water on foliage [mm]
   real(r8), allocatable :: ldew_rain_H2   (:) ! depth of rain on foliage [mm]
   real(r8), allocatable :: ldew_snow_H2   (:) ! depth of snow on foliage [mm]

   real(r8), allocatable :: scv_O18       (:) ! snow cover, water equivalent [mm]
   real(r8), allocatable :: scv_H2       (:) ! snow cover, water equivalent [mm]

   real(r8), allocatable :: wa_O18        (:) ! water storage in aquifer [mm]
   real(r8), allocatable :: wa_H2         (:) ! water storage in aquifer [mm]
   real(r8), allocatable :: wetwat_O18    (:) ! water storage in wetland [mm]
   real(r8), allocatable :: wetwat_H2    (:) ! water storage in wetland [mm]
   real(r8), allocatable :: wat_O18      (:) ! total water storage [mm]
   real(r8), allocatable :: wat_H2       (:) ! total water storage [mm]
   real(r8), allocatable :: wdsrf_O18    (:) ! depth of surface water [mm]
   real(r8), allocatable :: wdsrf_H2    (:) ! depth of surface water [mm]

   real(r8), allocatable :: qref_O18      (:) ! 2 m height air specific humidity of O18
   real(r8), allocatable :: qref_H2       (:) ! 2 m height air specific humidity of H2

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Tracer_allocate_TimeVariables
   PUBLIC :: Tracer_deallocate_TimeVariables
   PUBLIC :: Tracer_READ_TimeVariables
   PUBLIC :: Tracer_WRITE_TimeVariables
#ifdef RangeCheck
   PUBLIC :: Tracer_check_TimeVariables
#endif


!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE Tracer_allocate_TimeVariables
   !--------------------------------------------------------------------
   ! Allocates memory for CoLM 1d [numpatch] variables
   !--------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE


      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN


            allocate (wliq_soisno_O18(maxsnl+1:nl_soil,numpatch)); wliq_soisno_O18 (:,:) = spval
            allocate (wice_soisno_O18(maxsnl+1:nl_soil,numpatch)); wice_soisno_O18 (:,:) = spval
            allocate (wliq_soisno_H2(maxsnl+1:nl_soil,numpatch)); wliq_soisno_H2 (:,:)   = spval
            allocate (wice_soisno_H2(maxsnl+1:nl_soil,numpatch)); wice_soisno_H2 (:,:)   = spval

            allocate (rootr_O18         (1:nl_soil,numpatch)); rootr_O18   (:,:)         = spval
            allocate (rootr_H2         (1:nl_soil,numpatch)); rootr_H2   (:,:)           = spval
            allocate (rootflux_O18     (1:nl_soil,numpatch)); rootflux_O18(:,:)          = spval
            allocate (rootflux_H2     (1:nl_soil,numpatch)); rootflux_H2(:,:)            = spval

            allocate (ldew_O18                   (numpatch)); ldew_O18     (:)           = spval
            allocate (ldew_rain_O18              (numpatch)); ldew_rain_O18     (:)      = spval
            allocate (ldew_snow_O18              (numpatch)); ldew_snow_O18     (:)      = spval
            allocate (ldew_H2                    (numpatch)); ldew_H2     (:)            = spval
            allocate (ldew_rain_H2               (numpatch)); ldew_rain_H2     (:)       = spval
            allocate (ldew_snow_H2               (numpatch)); ldew_snow_H2     (:)       = spval

            allocate (scv_O18                     (numpatch)); scv_O18       (:)         = spval
            allocate (scv_H2                     (numpatch)); scv_H2       (:) = spval

            allocate (wa_O18                    (numpatch)); wa_O18        (:) = spval
            allocate (wa_H2                    (numpatch)); wa_H2        (:) = spval
            allocate (wetwat_O18                 (numpatch)); wetwat_O18    (:) = spval
            allocate (wetwat_H2                 (numpatch)); wetwat_H2    (:) = spval
            allocate (wat_O18                   (numpatch)); wat_O18      (:) = spval
            allocate (wat_H2                    (numpatch)); wat_H2       (:) = spval
            allocate (wdsrf_O18                 (numpatch)); wdsrf_O18    (:) = spval
            allocate (wdsrf_H2                 (numpatch)); wdsrf_H2    (:) = spval

            allocate (qref_O18                   (numpatch)); qref_O18      (:) = spval
            allocate (qref_H2                    (numpatch)); qref_H2       (:) = spval

         ENDIF
      ENDIF




   END SUBROUTINE Tracer_allocate_TimeVariables



   SUBROUTINE Tracer_deallocate_TimeVariables ()

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch
   IMPLICIT NONE

   !--------------------------------------------------------------------
   ! Deallocates memory for CoLM 1d [numpatch] variables
   !--------------------------------------------------------------------

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN

            deallocate (wliq_soisno_O18       )
            deallocate (wice_soisno_O18       )
            deallocate (wliq_soisno_H2         )
            deallocate (wice_soisno_H2         )

            deallocate (rootr_O18             )
            deallocate (rootr_H2             )
            deallocate (rootflux_O18         )
            deallocate (rootflux_H2         )

            deallocate (ldew_O18              )
            deallocate (ldew_rain_O18         )
            deallocate (ldew_snow_O18         )
            deallocate (ldew_H2               )
            deallocate (ldew_rain_H2          )
            deallocate (ldew_snow_H2          )

            deallocate (scv_O18               )
            deallocate (scv_H2               )

            deallocate (wa_O18                 )
            deallocate (wa_H2                  )
            deallocate (wetwat_O18             )
            deallocate (wetwat_H2              )
            deallocate (wat_O18                )
            deallocate (wat_H2                 )
            deallocate (wdsrf_O18              )
            deallocate (wdsrf_H2               )

            deallocate (qref_O18              )
            deallocate (qref_H2               )
       endIF
      ENDIF

   END SUBROUTINE Tracer_deallocate_TimeVariables


   !---------------------------------------
   FUNCTION Tracer_save_to_restart (idate, deltim, itstamp, ptstamp) result(rwrite)

   USE MOD_Namelist
   IMPLICIT NONE

   logical :: rwrite

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   type(timestamp), intent(in) :: itstamp, ptstamp


      ! added by yuan, 08/31/2014
      SELECTCASE (trim(adjustl(DEF_WRST_FREQ)))
      CASE ('TIMESTEP')
         rwrite = .true.
      CASE ('HOURLY')
         rwrite = isendofhour (idate, deltim)
      CASE ('DAILY')
         rwrite = isendofday(idate, deltim)
      CASE ('MONTHLY')
         rwrite = isendofmonth(idate, deltim)
      CASE ('YEARLY')
         rwrite = isendofyear(idate, deltim)
      CASE default
         rwrite = .false.
         write(*,*) 'Warning: Please USE one of TIMESTEP/HOURLY/DAILY/MONTHLY/YEARLY for restart frequency.'
         write(*,*) '         Set to FALSE by default.                                                     '
      ENDSELECT

      IF (rwrite) THEN
         rwrite = ((ptstamp <= itstamp) .or. isendofyear(idate,deltim))
      ENDIF

   END FUNCTION Tracer_save_to_restart


   SUBROUTINE Tracer_WRITE_TimeVariables (idate, lc_year, site, dir_restart)

   !====================================================================
   ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
   !====================================================================

   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_REST_CompressLevel, DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, &
                           DEF_USE_IRRIGATION, DEF_USE_Dynamic_Lake, SITE_landtype
   USE MOD_LandPatch
   USE MOD_NetCDFVector
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants, only: dz_lake
   USE MOD_Const_LC, only: patchtypes
   IMPLICIT NONE

   integer, intent(in) :: idate(3)
   integer, intent(in) :: lc_year      !year of land cover type data
   character(len=*), intent(in) :: site
   character(len=*), intent(in) :: dir_restart

   ! Local variables
   character(len=256) :: file_restart
   character(len=14)  :: cdate
   character(len=256) :: cyear         !character for lc_year
   integer :: compress

      compress = DEF_REST_CompressLevel

      ! land cover type year
      write(cyear,'(i4.4)') lc_year
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)

      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(dir_restart)//'/'//trim(cdate))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_tracer_restart_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'


      CALL ncio_create_file_vector (file_restart, landpatch)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')

      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake)

IF(DEF_USE_PLANTHYDRAULICS)THEN
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'vegnodes', nvegwcs)
ENDIF

      CALL ncio_define_dimension_vector (file_restart, landpatch, 'band', 2)
      CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp', 2)

      ! Time-varying state variables which required by restart run
 
      CALL ncio_write_vector (file_restart, 'wliq_soisno_O18', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wliq_soisno_O18, compress) ! liquid water in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'wice_soisno_O18', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wice_soisno_O18, compress) ! ice lens in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'wliq_soisno_H2', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wliq_soisno_H2, compress) ! liquid water in layers [kg/m2]
      CALL ncio_write_vector (file_restart, 'wice_soisno_H2', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wice_soisno_H2, compress) ! ice lens in layers [kg/m2]

      CALL ncio_write_vector (file_restart, 'ldew_O18'     , 'patch', landpatch, ldew_O18 ,      compress )
      CALL ncio_write_vector (file_restart, 'ldew_rain_O18', 'patch', landpatch, ldew_rain_O18 , compress )
      CALL ncio_write_vector (file_restart, 'ldew_snow_O18', 'patch', landpatch, ldew_snow_O18 , compress )
      CALL ncio_write_vector (file_restart, 'ldew_H2'      , 'patch', landpatch, ldew_H2 ,       compress )
      CALL ncio_write_vector (file_restart, 'ldew_rain_H2' , 'patch', landpatch, ldew_rain_H2 ,  compress )
      CALL ncio_write_vector (file_restart, 'ldew_snow_H2' , 'patch', landpatch, ldew_snow_H2 ,  compress )

      CALL ncio_write_vector (file_restart, 'scv_O18'   , 'patch', landpatch, scv_O18       , compress)                    ! snow cover, water equivalent [mm]
      CALL ncio_write_vector (file_restart, 'scv_H2'   , 'patch', landpatch, scv_H2       , compress)                    ! snow cover, water equivalent [mm]

      CALL ncio_write_vector (file_restart, 'wa_O18  '   , 'patch', landpatch, wa_O18    , compress)                    ! water storage in aquifer [mm]
      CALL ncio_write_vector (file_restart, 'wa_H2  '   , 'patch', landpatch, wa_H2    , compress)                    ! water storage in aquifer [mm]
      CALL ncio_write_vector (file_restart, 'wetwat_O18  '   , 'patch', landpatch, wetwat_O18    , compress)                    ! water storage in wetland [mm]
      CALL ncio_write_vector (file_restart, 'wetwat_H2  '   , 'patch', landpatch, wetwat_H2    , compress)                    ! water storage in wetland [mm]

      CALL ncio_write_vector (file_restart, 'qref_O18 ', 'patch', landpatch, qref_O18 , compress) ! 2 m height air specific humidity of O18
      CALL ncio_write_vector (file_restart, 'qref_H2 ', 'patch', landpatch, qref_H2 , compress) ! 2 m height air specific humidity of H2


   END SUBROUTINE Tracer_WRITE_TimeVariables


   SUBROUTINE Tracer_READ_TimeVariables (idate, lc_year, site, dir_restart)

   !====================================================================
   ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
   !====================================================================

   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_NetCDFVector
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_LandPatch
   USE MOD_Vars_Global
   USE MOD_Vars_TimeInvariants, only: dz_lake
   USE MOD_Const_LC, only: patchtypes

   IMPLICIT NONE

   integer, intent(in) :: idate(3)
   integer, intent(in) :: lc_year      !year of land cover type data
   character(len=*), intent(in) :: site
   character(len=*), intent(in) :: dir_restart

   ! Local variables
   character(len=256) :: file_restart
   character(len=14)  :: cdate, cyear

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write(*,'(/,A26)') 'Loading Time Variables ...'
      ENDIF

      ! land cover type year
      write(cyear,'(i4.4)') lc_year

      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
      file_restart = trim(dir_restart)// '/'//trim(cdate)//'/' // trim(site) //'_tracer_restart_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

      ! Time-varying state variables which required by restart run
 
      CALL ncio_read_vector (file_restart, 'wliq_soisno_O18', nl_soil-maxsnl, landpatch, wliq_soisno_O18) ! liquid water in layers [kg/m2]
      CALL ncio_read_vector (file_restart, 'wice_soisno_O18', nl_soil-maxsnl, landpatch, wice_soisno_O18) ! ice lens in layers [kg/m2]
      CALL ncio_read_vector (file_restart, 'wliq_soisno_H2', nl_soil-maxsnl, landpatch, wliq_soisno_H2)   ! liquid water in layers [kg/m2]
      CALL ncio_read_vector (file_restart, 'wice_soisno_H2', nl_soil-maxsnl, landpatch, wice_soisno_H2)   ! ice lens in layers [kg/m2]

      CALL ncio_read_vector (file_restart, 'ldew_O18',      landpatch, ldew_O18 )
      CALL ncio_read_vector (file_restart, 'ldew_rain_O18', landpatch, ldew_rain_O18 )
      CALL ncio_read_vector (file_restart, 'ldew_snow_O18', landpatch, ldew_snow_O18 )
      CALL ncio_read_vector (file_restart, 'ldew_H2',       landpatch, ldew_H2 )
      CALL ncio_read_vector (file_restart, 'ldew_rain_H2',  landpatch, ldew_rain_H2 )
      CALL ncio_read_vector (file_restart, 'ldew_snow_H2',  landpatch, ldew_snow_H2 )

      CALL ncio_read_vector (file_restart, 'scv_O18'   , landpatch, scv_O18       ) ! snow cover, water equivalent [mm]
      CALL ncio_read_vector (file_restart, 'scv_H2'   , landpatch, scv_H2       ) ! snow cover, water equivalent [mm]

      CALL ncio_read_vector (file_restart, 'wa_O18  '   , landpatch, wa_O18     ) ! water storage in aquifer [mm]
      CALL ncio_read_vector (file_restart, 'wa_H2  '   , landpatch, wa_H2     ) ! water storage in aquifer [mm]
      CALL ncio_read_vector (file_restart, 'wetwat_O18  '   , landpatch, wetwat_O18     ) ! water storage in wetland [mm]
      CALL ncio_read_vector (file_restart, 'wetwat_H2  '   , landpatch, wetwat_H2     ) ! water storage in wetland [mm]

      CALL ncio_read_vector (file_restart, 'qref_O18 ', landpatch, qref_O18 ) ! 2 m height air specific humidity of O18
      CALL ncio_read_vector (file_restart, 'qref_H2 ', landpatch, qref_H2 ) ! 2 m height air specific humidity of H2


#ifdef RangeCheck
      CALL Tracer_check_TimeVariables
#endif

      IF (p_is_master) THEN
         write(*,*) 'Loading Tracer Time Variables done.'
      ENDIF

   END SUBROUTINE Tracer_READ_TimeVariables


#ifdef RangeCheck
   SUBROUTINE Tracer_check_TimeVariables ()

   USE MOD_SPMD_Task
   USE MOD_RangeCheck
   USE MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS, DEF_USE_IRRIGATION, &
                           DEF_USE_SNICAR, DEF_USE_Dynamic_Lake
   USE MOD_Vars_TimeInvariants, only: dz_lake

   IMPLICIT NONE

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/,A27)') 'Checking Time Variables ...'
      ENDIF

      CALL check_vector_data ('ldew_O18       [mm]   ', ldew_O18       ) ! depth of water on foliage [mm]
      CALL check_vector_data ('ldew_rain_O18  [mm]   ', ldew_rain_O18  ) ! depth of rain on foliage [mm]
      CALL check_vector_data ('ldew_snow_O18  [mm]   ', ldew_snow_O18  ) ! depth of snow on foliage [mm]
      CALL check_vector_data ('ldew_H2       [mm]   ', ldew_H2       ) ! depth of water on foliage [mm]
      CALL check_vector_data ('ldew_rain_H2  [mm]   ', ldew_rain_H2  ) ! depth of rain on foliage [mm]
      CALL check_vector_data ('ldew_snow_H2  [mm]   ', ldew_snow_H2  ) ! depth of snow on foliage [mm]

      CALL check_vector_data ('scv_O18         [mm]   ', scv_O18        ) ! snow cover, water equivalent [mm]
      CALL check_vector_data ('scv_H2         [mm]   ', scv_H2        ) ! snow cover, water equivalent [mm]

      CALL check_vector_data ('wa_O18      [mm]   ', wa_O18     ) ! water storage in aquifer [mm]
      CALL check_vector_data ('wa_H2      [mm]   ', wa_H2       ) ! water storage in aquifer [mm]
      CALL check_vector_data ('wetwat_O18  [mm]   ', wetwat_O18 ) ! water storage in wetland [mm]
      CALL check_vector_data ('wetwat_H2  [mm]   ', wetwat_H2   ) ! water storage in wetland [mm]

      CALL check_vector_data ('wliq_soisno_O18 [kg/m2]', wliq_soisno_O18) ! liquid water in layers [kg/m2]
      CALL check_vector_data ('wice_soisno_O18 [kg/m2]', wice_soisno_O18) ! ice lens in layers [kg/m2]
      CALL check_vector_data ('wliq_soisno_H2 [kg/m2]', wliq_soisno_H2) ! liquid water in layers [kg/m2]
      CALL check_vector_data ('wice_soisno_H2 [kg/m2]', wice_soisno_H2) ! ice lens in layers [kg/m2]

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

  END SUBROUTINE Tracer_check_TimeVariables
#endif


END MODULE MOD_Tracer_Vars_TimeVariables
! ---------- EOP ------------
