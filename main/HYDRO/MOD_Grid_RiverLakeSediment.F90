#include <define.h>

#ifdef GridRiverLakeSediment
MODULE MOD_Grid_RiverLakeSediment
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Sediment transport module for GridRiverLakeFlow.
!   Ported from CaMa-Flood sediment module (cmf_ctrl_sed_mod.F90)
!
! Created by: Claude Code, Dec 2025
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   IMPLICIT NONE

   !-------------------------------------------------------------------------------------
   ! Module Parameters
   !-------------------------------------------------------------------------------------
   integer,  save :: nsed           ! Number of sediment size classes
   integer,  save :: totlyrnum      ! Number of deposition layers
   integer,  save :: nlfp_sed       ! Number of floodplain layers for slope

   real(r8), save :: lambda         ! Porosity [-]
   real(r8), save :: lyrdph         ! Active layer depth [m]
   real(r8), save :: psedD          ! Sediment density [g/cm3]
   real(r8), save :: pwatD          ! Water density [g/cm3]
   real(r8), save :: visKin         ! Kinematic viscosity [m2/s]
   real(r8), save :: vonKar         ! Von Karman coefficient [-]

   ! Sediment yield parameters
   real(r8), save :: pyld           ! Yield coefficient
   real(r8), save :: pyldc          ! Slope exponent
   real(r8), save :: pyldpc         ! Precipitation exponent
   real(r8), save :: dsylunit       ! Unit conversion factor

   !-------------------------------------------------------------------------------------
   ! Static Data (read from DEF_UnitCatchment_file)
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sed_frc   (:,:)    ! Sediment fraction [nsed, numucat]
   real(r8), allocatable :: sed_slope (:,:)    ! Floodplain slope [nlfp_sed, numucat]
   real(r8), allocatable :: sDiam     (:)      ! Grain diameter [nsed]
   real(r8), allocatable :: setvel    (:)      ! Settling velocity [nsed]

   !-------------------------------------------------------------------------------------
   ! State Variables
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sedcon  (:,:)      ! Suspended sediment concentration [numucat, nsed]
   real(r8), allocatable :: layer   (:,:)      ! Active layer storage [numucat, nsed]
   real(r8), allocatable :: seddep  (:,:,:)    ! Deposition layer storage [numucat, totlyrnum, nsed]

   !-------------------------------------------------------------------------------------
   ! Diagnostic Variables
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: sedout  (:,:)      ! Suspended sediment outflow [numucat, nsed]
   real(r8), allocatable :: bedout  (:,:)      ! Bedload outflow [numucat, nsed]
   real(r8), allocatable :: sedinp  (:,:)      ! Erosion input [numucat, nsed]
   real(r8), allocatable :: netflw  (:,:)      ! Net exchange flux [numucat, nsed]
   real(r8), allocatable :: shearvel(:)        ! Shear velocity [numucat]
   real(r8), allocatable :: critshearvel(:,:)  ! Critical shear velocity [numucat, nsed]
   real(r8), allocatable :: susvel  (:,:)      ! Suspension velocity [numucat, nsed]

   !-------------------------------------------------------------------------------------
   ! Accumulated Variables for Sediment Time-stepping
   !-------------------------------------------------------------------------------------
   real(r8), save :: sed_acc_time              ! Accumulated time for averaging
   real(r8), allocatable :: sed_acc_veloc(:)   ! Accumulated velocity [numucat]
   real(r8), allocatable :: sed_acc_wdsrf(:)   ! Accumulated water depth [numucat]
   real(r8), allocatable :: sed_precip(:)      ! Precipitation for sediment yield [numucat]

   !-------------------------------------------------------------------------------------
   ! Accumulated Variables for History Output
   !-------------------------------------------------------------------------------------
   real(r8), allocatable :: a_sedcon  (:,:)    ! Accumulated sedcon
   real(r8), allocatable :: a_sedout  (:,:)    ! Accumulated sedout
   real(r8), allocatable :: a_bedout  (:,:)    ! Accumulated bedout
   real(r8), allocatable :: a_sedinp  (:,:)    ! Accumulated sedinp
   real(r8), allocatable :: a_netflw  (:,:)    ! Accumulated netflw
   real(r8), allocatable :: a_layer   (:,:)    ! Accumulated layer
   real(r8), allocatable :: a_shearvel(:)      ! Accumulated shearvel

   !-------------------------------------------------------------------------------------
   ! Public Subroutines
   !-------------------------------------------------------------------------------------
   PUBLIC :: grid_sediment_init
   PUBLIC :: grid_sediment_calc
   PUBLIC :: grid_sediment_final
   PUBLIC :: sediment_diag_accumulate
   PUBLIC :: sediment_forcing_put

CONTAINS

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_init()
   ! Initialize sediment module
   !-------------------------------------------------------------------------------------
   USE MOD_NetCDFSerial
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, &
      ucat_data_address, topo_rivwth, topo_rivlen
   IMPLICIT NONE

   character(len=256) :: parafile
   integer :: i

      IF (.not. DEF_USE_SEDIMENT) RETURN

      IF (p_is_io) THEN
         WRITE(*,*) 'Initializing sediment module...'
      ENDIF

      ! Set parameters from namelist
      lambda    = DEF_SED_LAMBDA
      lyrdph    = DEF_SED_LYRDPH
      psedD     = DEF_SED_DENSITY
      pwatD     = DEF_SED_WATER_DENSITY
      visKin    = DEF_SED_VISKIN
      vonKar    = DEF_SED_VONKAR
      totlyrnum = DEF_SED_TOTLYRNUM
      pyld      = DEF_SED_PYLD
      pyldc     = DEF_SED_PYLDC
      pyldpc    = DEF_SED_PYLDPC
      dsylunit  = DEF_SED_DSYLUNIT

      parafile = DEF_UnitCatchment_file

      ! Read dimensions from NetCDF file
      IF (p_is_master) THEN
         CALL ncio_inquire_length(parafile, 'sed_n', nsed)
         CALL ncio_inquire_length(parafile, 'slope_layers', nlfp_sed)
      ENDIF

#ifdef USEMPI
      CALL mpi_bcast(nsed, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast(nlfp_sed, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      IF (p_is_io) THEN
         WRITE(*,*) 'Sediment module: nsed=', nsed, ' nlfp_sed=', nlfp_sed
      ENDIF

      ! Parse grain diameters from string
      CALL parse_grain_diameters()

      ! Calculate settling velocities
      CALL calc_settling_velocities()

      ! Read static data
      CALL read_sediment_static_data(parafile)

      ! Allocate and initialize state variables
      CALL allocate_sediment_vars()

      ! Initialize state from sed_frc
      CALL initialize_sediment_state()

      IF (p_is_io) THEN
         WRITE(*,*) 'Sediment module initialized successfully.'
      ENDIF

   END SUBROUTINE grid_sediment_init

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_calc(deltime)
   ! Main sediment calculation routine
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: deltime
      ! Placeholder - to be implemented
      WRITE(*,*) 'MOD_Grid_RiverLakeSediment: grid_sediment_calc called, dt=', deltime
   END SUBROUTINE grid_sediment_calc

   !-------------------------------------------------------------------------------------
   SUBROUTINE sediment_diag_accumulate(dt, veloc, wdsrf)
   ! Accumulate water flow variables for sediment calculation
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: dt
   real(r8), intent(in) :: veloc(:)
   real(r8), intent(in) :: wdsrf(:)
      ! Placeholder - to be implemented
   END SUBROUTINE sediment_diag_accumulate

   !-------------------------------------------------------------------------------------
   SUBROUTINE sediment_forcing_put(precip)
   ! Store precipitation for sediment yield calculation
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: precip(:)
      ! Placeholder - to be implemented
   END SUBROUTINE sediment_forcing_put

   !-------------------------------------------------------------------------------------
   SUBROUTINE parse_grain_diameters()
   ! Parse grain diameters from comma-separated string
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   character(len=256) :: str
   integer :: i, j, k, n
   integer :: iostat

      allocate(sDiam(nsed))
      str = trim(adjustl(DEF_SED_DIAMETER))

      n = 0
      j = 1
      DO i = 1, len_trim(str)
         IF (str(i:i) == ',' .or. i == len_trim(str)) THEN
            n = n + 1
            IF (i == len_trim(str)) THEN
               k = i
            ELSE
               k = i - 1
            ENDIF
            IF (n <= nsed) THEN
               read(str(j:k), *, iostat=iostat) sDiam(n)
               IF (iostat /= 0) THEN
                  IF (p_is_io) THEN
                     WRITE(*,*) 'ERROR: Failed to parse grain diameter at position', n
                     WRITE(*,*) '  String fragment: "', str(j:k), '"'
                  ENDIF
                  STOP
               ENDIF
            ENDIF
            j = i + 1
         ENDIF
      ENDDO

      IF (n /= nsed) THEN
         IF (p_is_io) THEN
            WRITE(*,*) 'ERROR: Number of diameters does not match nsed'
            WRITE(*,*) '  Parsed:', n, ' Expected:', nsed
         ENDIF
         STOP
      ENDIF

      ! Validate that all parsed diameter values are positive
      DO i = 1, nsed
         IF (sDiam(i) <= 0._r8) THEN
            IF (p_is_io) THEN
               WRITE(*,*) 'ERROR: Grain diameter must be positive'
               WRITE(*,*) '  Class:', i, ' Value:', sDiam(i)
            ENDIF
            STOP
         ENDIF
      ENDDO

      IF (p_is_io) THEN
         WRITE(*,*) 'Grain diameters (m):', sDiam
      ENDIF

   END SUBROUTINE parse_grain_diameters

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_settling_velocities()
   ! Calculate settling velocity using Stokes-Rubey formula
   !-------------------------------------------------------------------------------------
   USE MOD_Const_Physical, only: grav
   IMPLICIT NONE
   real(r8) :: sTmp
   integer :: i

      allocate(setvel(nsed))

      DO i = 1, nsed
         sTmp = 6.0_r8 * visKin / sDiam(i)
         setvel(i) = sqrt(2.0_r8/3.0_r8 * (psedD-pwatD)/pwatD * grav * sDiam(i) &
                    + sTmp*sTmp) - sTmp
      ENDDO

      IF (p_is_io) THEN
         WRITE(*,*) 'Settling velocities (m/s):', setvel
      ENDIF

   END SUBROUTINE calc_settling_velocities

   !-------------------------------------------------------------------------------------
   SUBROUTINE read_sediment_static_data(parafile)
   ! Read sediment static data from NetCDF file
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, readin_riverlake_parameter
   IMPLICIT NONE
   character(len=*), intent(in) :: parafile

      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            allocate(sed_frc  (nsed, numucat))
            allocate(sed_slope(nlfp_sed, numucat))
         ENDIF
      ENDIF

      ! Read sediment fraction
      CALL readin_riverlake_parameter(parafile, 'sed_frc', rdata2d=sed_frc)

      ! Read sediment slope
      CALL readin_riverlake_parameter(parafile, 'sed_slope', rdata2d=sed_slope)

      ! Normalize sed_frc (ensure sum = 1)
      CALL normalize_sed_frc()

      IF (p_is_io) THEN
         WRITE(*,*) 'Sediment static data read successfully.'
      ENDIF

   END SUBROUTINE read_sediment_static_data

   !-------------------------------------------------------------------------------------
   SUBROUTINE normalize_sed_frc()
   ! Normalize sediment fractions to ensure they sum to 1.0 for each unit catchment
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE
   integer :: i
   real(r8) :: frc_sum

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         frc_sum = sum(sed_frc(:,i))
         IF (frc_sum > 0._r8) THEN
            sed_frc(:,i) = sed_frc(:,i) / frc_sum
         ELSE
            ! If all fractions are zero or negative, set equal fractions
            sed_frc(:,i) = 1._r8 / real(nsed, r8)
         ENDIF
      ENDDO

      IF (p_is_io) THEN
         WRITE(*,*) 'Sediment fractions normalized.'
      ENDIF

   END SUBROUTINE normalize_sed_frc

   !-------------------------------------------------------------------------------------
   SUBROUTINE allocate_sediment_vars()
   ! Allocate sediment state and diagnostic variables
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      ! State variables
      allocate(sedcon(nsed, numucat))
      allocate(layer (nsed, numucat))
      allocate(seddep(nsed, totlyrnum, numucat))

      ! Diagnostic variables
      allocate(sedout      (nsed, numucat))
      allocate(bedout      (nsed, numucat))
      allocate(sedinp      (nsed, numucat))
      allocate(netflw      (nsed, numucat))
      allocate(shearvel    (numucat))
      allocate(critshearvel(nsed, numucat))
      allocate(susvel      (nsed, numucat))

      ! Accumulation variables
      allocate(sed_acc_veloc(numucat))
      allocate(sed_acc_wdsrf(numucat))
      allocate(sed_precip   (numucat))

      ! History output variables
      allocate(a_sedcon  (nsed, numucat))
      allocate(a_sedout  (nsed, numucat))
      allocate(a_bedout  (nsed, numucat))
      allocate(a_sedinp  (nsed, numucat))
      allocate(a_netflw  (nsed, numucat))
      allocate(a_layer   (nsed, numucat))
      allocate(a_shearvel(numucat))

      ! Initialize to zero
      sedcon       = 0._r8
      layer        = 0._r8
      seddep       = 0._r8
      sedout       = 0._r8
      bedout       = 0._r8
      sedinp       = 0._r8
      netflw       = 0._r8
      shearvel     = 0._r8
      critshearvel = 0._r8
      susvel       = 0._r8
      sed_acc_veloc = 0._r8
      sed_acc_wdsrf = 0._r8
      sed_acc_time  = 0._r8
      sed_precip    = 0._r8
      a_sedcon     = 0._r8
      a_sedout     = 0._r8
      a_bedout     = 0._r8
      a_sedinp     = 0._r8
      a_netflw     = 0._r8
      a_layer      = 0._r8
      a_shearvel   = 0._r8

   END SUBROUTINE allocate_sediment_vars

   !-------------------------------------------------------------------------------------
   SUBROUTINE initialize_sediment_state()
   ! Initialize sediment state from sed_frc
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat, topo_rivwth, topo_rivlen
   IMPLICIT NONE
   integer :: i, ilyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      ! Initialize active layer based on sed_frc
      DO i = 1, numucat
         layer(:,i) = lyrdph * topo_rivwth(i) * topo_rivlen(i) * sed_frc(:,i)

         ! Initialize deposition layers
         DO ilyr = 1, totlyrnum - 1
            seddep(:,ilyr,i) = layer(:,i)
         ENDDO
         ! Bottom layer gets extra depth
         seddep(:,totlyrnum,i) = max(10._r8 - lyrdph*totlyrnum, 0._r8) &
            * topo_rivwth(i) * topo_rivlen(i) * sed_frc(:,i)
      ENDDO

   END SUBROUTINE initialize_sediment_state

   !-------------------------------------------------------------------------------------
   SUBROUTINE grid_sediment_final()
   ! Cleanup sediment module
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
      IF (allocated(sed_frc     )) deallocate(sed_frc     )
      IF (allocated(sed_slope   )) deallocate(sed_slope   )
      IF (allocated(sDiam       )) deallocate(sDiam       )
      IF (allocated(setvel      )) deallocate(setvel      )
      IF (allocated(sedcon      )) deallocate(sedcon      )
      IF (allocated(layer       )) deallocate(layer       )
      IF (allocated(seddep      )) deallocate(seddep      )
      IF (allocated(sedout      )) deallocate(sedout      )
      IF (allocated(bedout      )) deallocate(bedout      )
      IF (allocated(sedinp      )) deallocate(sedinp      )
      IF (allocated(netflw      )) deallocate(netflw      )
      IF (allocated(shearvel    )) deallocate(shearvel    )
      IF (allocated(critshearvel)) deallocate(critshearvel)
      IF (allocated(susvel      )) deallocate(susvel      )
      IF (allocated(sed_acc_veloc)) deallocate(sed_acc_veloc)
      IF (allocated(sed_acc_wdsrf)) deallocate(sed_acc_wdsrf)
      IF (allocated(sed_precip  )) deallocate(sed_precip  )
      IF (allocated(a_sedcon    )) deallocate(a_sedcon    )
      IF (allocated(a_sedout    )) deallocate(a_sedout    )
      IF (allocated(a_bedout    )) deallocate(a_bedout    )
      IF (allocated(a_sedinp    )) deallocate(a_sedinp    )
      IF (allocated(a_netflw    )) deallocate(a_netflw    )
      IF (allocated(a_layer     )) deallocate(a_layer     )
      IF (allocated(a_shearvel  )) deallocate(a_shearvel  )
   END SUBROUTINE grid_sediment_final

END MODULE MOD_Grid_RiverLakeSediment
#endif
