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
   FUNCTION calc_shear_velocity(rivvel, rivdph, rivman) RESULT(svel)
   ! Calculate shear velocity using Manning's equation
   !-------------------------------------------------------------------------------------
   USE MOD_Const_Physical, only: grav
   IMPLICIT NONE
   real(r8), intent(in) :: rivvel   ! River velocity [m/s]
   real(r8), intent(in) :: rivdph   ! River depth [m]
   real(r8), intent(in) :: rivman   ! Manning coefficient
   real(r8) :: svel

      IF (rivdph > 0._r8) THEN
         svel = sqrt(grav * rivman**2 * rivvel**2 * rivdph**(-1._r8/3._r8))
      ELSE
         svel = 0._r8
      ENDIF

   END FUNCTION calc_shear_velocity

   !-------------------------------------------------------------------------------------
   FUNCTION calc_critical_shear_velocity(diam) RESULT(csvel)
   ! Calculate critical shear velocity using Shields curve
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in) :: diam    ! Grain diameter [m]
   real(r8) :: csvel               ! [(cm/s)^2]
   real(r8) :: cA, cB

      cB = 1._r8
      IF (diam >= 0.00303_r8) THEN
         cA = 80.9_r8
      ELSEIF (diam >= 0.00118_r8) THEN
         cA = 134.6_r8
         cB = 31._r8 / 32._r8
      ELSEIF (diam >= 0.000565_r8) THEN
         cA = 55._r8
      ELSEIF (diam >= 0.000065_r8) THEN
         cA = 8.41_r8
         cB = 11._r8 / 32._r8
      ELSE
         cA = 226._r8
      ENDIF

      csvel = cA * (diam * 100._r8) ** cB

   END FUNCTION calc_critical_shear_velocity

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_critical_shear_egiazoroff(i, svel, csvel_out)
   ! Calculate critical shear velocity using Egiazoroff equation for mixed-size sediment
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   integer,  intent(in)  :: i           ! Unit catchment index
   real(r8), intent(in)  :: svel        ! Shear velocity [m/s]
   real(r8), intent(out) :: csvel_out(nsed)

   real(r8) :: dMean, csVel0
   integer  :: ised
   real(r8) :: layer_sum

      layer_sum = sum(layer(:,i))

      IF (layer_sum <= 0._r8) THEN
         csvel_out(:) = 1.e20_r8
         RETURN
      ENDIF

      ! Calculate mean diameter
      dMean = 0._r8
      DO ised = 1, nsed
         dMean = dMean + sDiam(ised) * layer(ised,i) / layer_sum
      ENDDO

      csVel0 = calc_critical_shear_velocity(dMean)

      DO ised = 1, nsed
         IF (sDiam(ised) / dMean >= 0.4_r8) THEN
            csvel_out(ised) = sqrt(csVel0 * sDiam(ised) / dMean) * &
               (log10(19._r8) / log10(19._r8 * sDiam(ised) / dMean)) * 0.01_r8
         ELSE
            csvel_out(ised) = sqrt(0.85_r8 * csVel0) * 0.01_r8
         ENDIF
      ENDDO

   END SUBROUTINE calc_critical_shear_egiazoroff

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_suspend_velocity(csvel, svel, susvel_out)
   ! Calculate suspension velocity using Uchida & Fukuoka (2019) Eq.44
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   real(r8), intent(in)  :: csvel(nsed)  ! Critical shear velocity [m/s]
   real(r8), intent(in)  :: svel         ! Shear velocity [m/s]
   real(r8), intent(out) :: susvel_out(nsed)

   real(r8) :: alpha, a, cB, sTmp
   integer  :: ised

      alpha = vonKar / 6._r8
      a = 0.08_r8
      cB = 1._r8 - lambda

      susvel_out(:) = 0._r8

      DO ised = 1, nsed
         IF (csvel(ised) > svel) CYCLE
         IF (svel <= 0._r8) CYCLE

         sTmp = setvel(ised) / alpha / svel
         susvel_out(ised) = max(setvel(ised) * cB / (1._r8 + sTmp) * &
            (1._r8 - a*sTmp) / (1._r8 + (1._r8-a)*sTmp), 0._r8)
      ENDDO

   END SUBROUTINE calc_suspend_velocity

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_sediment_advection(dt, rivout, rivsto)
   ! Calculate suspended sediment and bedload advection
   !-------------------------------------------------------------------------------------
   USE MOD_Const_Physical, only: grav
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_next, topo_rivwth, push_next2ucat
   USE MOD_WorkerPushData
   IMPLICIT NONE

   real(r8), intent(in) :: dt           ! Time step [s]
   real(r8), intent(in) :: rivout(:)    ! River outflow [m3/s]
   real(r8), intent(in) :: rivsto(:)    ! River storage [m3]

   real(r8), allocatable :: sedsto(:,:)       ! Sediment storage [nsed, numucat]
   real(r8), allocatable :: bOut(:,:), sOut(:,:)
   real(r8), allocatable :: brate(:,:), srate(:,:)
   real(r8), allocatable :: sedcon_next(:,:)

   integer  :: i, ised, i0, i1
   real(r8) :: plusVel, minusVel
   real(r8) :: layer_sum

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      allocate(sedsto(nsed, numucat))
      allocate(bOut  (nsed, numucat))
      allocate(sOut  (nsed, numucat))
      allocate(brate (nsed, numucat))
      allocate(srate (nsed, numucat))
      allocate(sedcon_next(nsed, numucat))

      ! Calculate sediment storage from concentration
      DO i = 1, numucat
         sedsto(:,i) = sedcon(:,i) * max(rivsto(i), 0._r8)
      ENDDO

      ! Get downstream sediment concentration
      DO ised = 1, nsed
         CALL worker_push_data(push_next2ucat, sedcon(ised,:), sedcon_next(ised,:), &
            fillvalue = 0._r8)
      ENDDO

      bOut(:,:) = 0._r8
      sOut(:,:) = 0._r8

      ! Calculate outflows
      DO i = 1, numucat
         IF (rivout(i) >= 0._r8) THEN
            i0 = i
            i1 = ucat_next(i)
         ELSE
            i0 = ucat_next(i)
            i1 = i
         ENDIF

         IF (rivout(i) == 0._r8) THEN
            sedout(:,i) = 0._r8
            bedout(:,i) = 0._r8
            CYCLE
         ENDIF

         ! Suspended sediment outflow
         IF (i0 < 0) THEN
            sedout(:,i) = sedcon(:,i1) * rivout(i)
         ELSE
            sedout(:,i) = sedcon(:,i0) * rivout(i)
            sOut(:,i0) = sOut(:,i0) + abs(sedout(:,i)) * dt
         ENDIF

         ! Bedload outflow
         layer_sum = sum(layer(:,i))
         IF (all(critshearvel(:,i) >= shearvel(i)) .or. layer_sum == 0._r8 .or. i0 < 0) THEN
            bedout(:,i) = 0._r8
         ELSE
            DO ised = 1, nsed
               IF (critshearvel(ised,i) >= shearvel(i) .or. layer(ised,i) == 0._r8) THEN
                  bedout(ised,i) = 0._r8
                  CYCLE
               ENDIF
               plusVel = shearvel(i) + critshearvel(ised,i)
               minusVel = shearvel(i) - critshearvel(ised,i)
               bedout(ised,i) = 17._r8 * topo_rivwth(i) * plusVel * minusVel * minusVel &
                  / ((psedD-pwatD)/pwatD) / grav * layer(ised,i) / layer_sum
               IF (i0 > 0) bOut(ised,i0) = bOut(ised,i0) + bedout(ised,i) * dt
            ENDDO
         ENDIF
      ENDDO

      ! Adjust outflow if larger than available storage
      brate(:,:) = 1._r8
      srate(:,:) = 1._r8

      DO i = 1, numucat
         DO ised = 1, nsed
            IF (sOut(ised,i) > 1.e-8_r8) THEN
               srate(ised,i) = min(sedsto(ised,i) / sOut(ised,i), 1._r8)
            ENDIF
            IF (bOut(ised,i) > 1.e-8_r8) THEN
               brate(ised,i) = min(layer(ised,i) / bOut(ised,i), 1._r8)
            ENDIF
         ENDDO
      ENDDO

      ! Apply adjusted outflows and update storage
      DO i = 1, numucat
         IF (rivout(i) >= 0._r8) THEN
            i0 = i
            i1 = ucat_next(i)
         ELSE
            i0 = ucat_next(i)
            i1 = i
         ENDIF

         IF (i0 > 0) THEN
            sedout(:,i) = sedout(:,i) * srate(:,i0)
            sedsto(:,i0) = max(sedsto(:,i0) - abs(sedout(:,i)) * dt, 0._r8)
            bedout(:,i) = bedout(:,i) * brate(:,i0)
            layer(:,i0) = max(layer(:,i0) - abs(bedout(:,i)) * dt, 0._r8)
         ENDIF

         IF (i1 > 0) THEN
            sedsto(:,i1) = max(sedsto(:,i1) + abs(sedout(:,i)) * dt, 0._r8)
            layer(:,i1) = max(layer(:,i1) + abs(bedout(:,i)) * dt, 0._r8)
         ENDIF
      ENDDO

      ! Update concentration from storage
      DO i = 1, numucat
         IF (rivsto(i) > 0._r8) THEN
            sedcon(:,i) = sedsto(:,i) / rivsto(i)
         ENDIF
      ENDDO

      deallocate(sedsto, bOut, sOut, brate, srate, sedcon_next)

   END SUBROUTINE calc_sediment_advection

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_sediment_exchange(dt, rivsto, rivwth, rivlen)
   ! Calculate suspension-deposition exchange
   !-------------------------------------------------------------------------------------
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: dt
   real(r8), intent(in) :: rivsto(:)
   real(r8), intent(in) :: rivwth(:)
   real(r8), intent(in) :: rivlen(:)

   real(r8) :: Es(nsed), D(nsed), Zd(nsed)
   real(r8) :: sedsto(nsed), dTmp(nsed)
   real(r8) :: layer_sum, area, dTmp1
   integer  :: i, ised
   real(r8), parameter :: IGNORE_DPH = 0.05_r8

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         IF (rivsto(i) < rivwth(i) * rivlen(i) * IGNORE_DPH) THEN
            netflw(:,i) = 0._r8
            CYCLE
         ENDIF

         layer_sum = sum(layer(:,i))
         area = rivwth(i) * rivlen(i)

         ! Calculate entrainment (suspension)
         IF (layer_sum == 0._r8 .or. all(susvel(:,i) == 0._r8)) THEN
            Es(:) = 0._r8
         ELSE
            Es(:) = susvel(:,i) * (1._r8 - lambda) * area * layer(:,i) / layer_sum
            Es(:) = max(Es(:), 0._r8)
         ENDIF

         ! Calculate deposition
         IF (shearvel(i) == 0._r8 .or. all(setvel(:) == 0._r8)) THEN
            D(:) = 0._r8
         ELSE
            DO ised = 1, nsed
               Zd(ised) = 6._r8 * setvel(ised) / vonKar / shearvel(i)
               D(ised) = setvel(ised) * area * sedcon(ised,i) * &
                  Zd(ised) / (1._r8 - exp(-Zd(ised)))
            ENDDO
            D(:) = max(D(:), 0._r8)
         ENDIF

         netflw(:,i) = Es(:) - D(:)

         ! Apply exchange with mass conservation
         sedsto(:) = sedcon(:,i) * rivsto(i)

         DO ised = 1, nsed
            IF (netflw(ised,i) == 0._r8) THEN
               CYCLE
            ELSEIF (netflw(ised,i) > 0._r8) THEN
               ! Suspension: transfer from layer to water
               dTmp1 = netflw(ised,i) * dt / (1._r8 - lambda)
               IF (dTmp1 < layer(ised,i)) THEN
                  layer(ised,i) = layer(ised,i) - dTmp1
               ELSE
                  netflw(ised,i) = layer(ised,i) * (1._r8 - lambda) / dt
                  layer(ised,i) = 0._r8
               ENDIF
               sedsto(ised) = sedsto(ised) + netflw(ised,i) * dt
            ELSE
               ! Deposition: transfer from water to layer
               IF (abs(netflw(ised,i)) * dt < sedsto(ised)) THEN
                  sedsto(ised) = max(sedsto(ised) - abs(netflw(ised,i)) * dt, 0._r8)
               ELSE
                  netflw(ised,i) = -sedsto(ised) / dt
                  sedsto(ised) = 0._r8
               ENDIF
               layer(ised,i) = layer(ised,i) + abs(netflw(ised,i)) * dt / (1._r8 - lambda)
            ENDIF
         ENDDO

         ! Add erosion input
         sedsto(:) = sedsto(:) + sedinp(:,i) * dt

         ! Limit concentration to 1%
         IF (sum(sedsto(:)) > rivsto(i) * 0.01_r8) THEN
            dTmp(:) = (sum(sedsto(:)) - rivsto(i) * 0.01_r8) * sedsto(:) / sum(sedsto(:))
            netflw(:,i) = netflw(:,i) - dTmp(:) / dt
            sedsto(:) = sedsto(:) - dTmp(:)
            layer(:,i) = layer(:,i) + dTmp(:) / (1._r8 - lambda)
         ENDIF

         ! Update concentration
         IF (rivsto(i) > 0._r8) THEN
            sedcon(:,i) = sedsto(:) / rivsto(i)
         ENDIF
      ENDDO

   END SUBROUTINE calc_sediment_exchange

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_layer_redistribution(rivwth, rivlen)
   ! Redistribute sediment into vertical bed layers
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: rivwth(:)
   real(r8), intent(in) :: rivlen(:)

   real(r8) :: lyrvol, diff
   real(r8) :: layerP(nsed), seddepP(totlyrnum+1, nsed), tmp(nsed)
   integer  :: i, ilyr, jlyr, slyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      DO i = 1, numucat
         lyrvol = lyrdph * rivwth(i) * rivlen(i)

         ! Ensure non-negative values
         layer(:,i) = max(layer(:,i), 0._r8)
         seddep(:,:,i) = max(seddep(:,:,i), 0._r8)

         ! If total bed storage less than layer volume
         IF (sum(layer(:,i)) + sum(seddep(:,:,i)) <= lyrvol) THEN
            layer(:,i) = layer(:,i) + sum(seddep(:,:,i), dim=2)
            seddep(:,:,i) = 0._r8
            CYCLE
         ENDIF

         ! Distribute into top exchange layer
         layerP(:) = layer(:,i)
         IF (sum(layerP(:)) >= lyrvol) THEN
            layer(:,i) = layerP(:) * min(lyrvol / sum(layerP(:)), 1._r8)
            layerP(:) = max(layerP(:) - layer(:,i), 0._r8)
            slyr = 0
         ELSEIF (sum(seddep(:,:,i)) > 0._r8) THEN
            layerP(:) = 0._r8
            DO ilyr = 1, totlyrnum
               diff = lyrvol - sum(layer(:,i))
               IF (diff <= 0._r8) EXIT
               IF (sum(seddep(:,ilyr,i)) <= diff) THEN
                  layer(:,i) = layer(:,i) + seddep(:,ilyr,i)
                  seddep(:,ilyr,i) = 0._r8
                  slyr = ilyr + 1
               ELSE
                  tmp(:) = diff * seddep(:,ilyr,i) / sum(seddep(:,ilyr,i))
                  layer(:,i) = layer(:,i) + tmp(:)
                  seddep(:,ilyr,i) = max(seddep(:,ilyr,i) - tmp(:), 0._r8)
                  slyr = ilyr
                  EXIT
               ENDIF
            ENDDO
         ELSE
            seddep(:,:,i) = 0._r8
            CYCLE
         ENDIF

         IF (sum(seddep(:,:,i)) == 0._r8) CYCLE

         ! Distribute remaining bedload into vertical deposition layers
         seddepP(1,:) = layerP(:)
         seddepP(2:,:) = seddep(:,:,i)
         seddep(:,:,i) = 0._r8

         DO ilyr = 1, totlyrnum - 1
            IF (sum(seddep(:,ilyr,i)) == lyrvol) CYCLE
            DO jlyr = slyr + 1, totlyrnum + 1
               diff = lyrvol - sum(seddep(:,ilyr,i))
               IF (diff <= 0._r8) EXIT
               IF (sum(seddepP(jlyr,:)) <= diff) THEN
                  seddep(:,ilyr,i) = seddep(:,ilyr,i) + seddepP(jlyr,:)
                  seddepP(jlyr,:) = 0._r8
               ELSE
                  tmp(:) = diff * seddepP(jlyr,:) / sum(seddepP(jlyr,:))
                  seddep(:,ilyr,i) = seddep(:,ilyr,i) + tmp(:)
                  seddepP(jlyr,:) = max(seddepP(jlyr,:) - tmp(:), 0._r8)
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         IF (sum(seddepP) > 0._r8) THEN
            seddep(:,totlyrnum,i) = sum(seddepP, dim=1)
         ENDIF
      ENDDO

   END SUBROUTINE calc_layer_redistribution

   !-------------------------------------------------------------------------------------
   SUBROUTINE calc_sediment_yield(fldfrc, grarea)
   ! Calculate sediment yield from precipitation
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   IMPLICIT NONE

   real(r8), intent(in) :: fldfrc(:)   ! Flooded fraction
   real(r8), intent(in) :: grarea(:)   ! Grid area [m2]

   real(r8) :: precip_mm
   integer  :: i, ilyr

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      sedinp(:,:) = 0._r8

      DO i = 1, numucat
         ! Convert precip from mm/s or kg/m2/s to mm/day for threshold
         precip_mm = sed_precip(i) * 86400._r8

         IF (precip_mm <= 10._r8) CYCLE

         ! Calculate erosion for each floodplain layer
         DO ilyr = 1, nlfp_sed
            IF (fldfrc(i) * nlfp_sed > real(ilyr, r8)) CYCLE  ! No erosion if submerged

            sedinp(:,i) = sedinp(:,i) + &
               pyld * (sed_precip(i) * 3600._r8)**pyldpc * sed_slope(ilyr,i)**pyldc / 3600._r8 &
               * grarea(i) * min(real(ilyr, r8)/real(nlfp_sed, r8) - fldfrc(i), 1._r8/real(nlfp_sed, r8)) &
               * dsylunit * sed_frc(:,i)
         ENDDO
      ENDDO

   END SUBROUTINE calc_sediment_yield

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
