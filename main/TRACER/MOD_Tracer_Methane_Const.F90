#include <define.h>
#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Methane_Const
!=======================================================================
! methane constants
!=======================================================================
	USE MOD_Precision

	IMPLICIT NONE

	PUBLIC

	!------------------------------------------------------------------
	! Constants for the Methane reactive tracer
	!------------------------------------------------------------------
	! Note some of these constants are also used in CNNitrifDenitrifMod

	integer :: iloop  ! loop index

	integer, public, parameter :: ngases      =   3     ! CH4, O2, & CO2

	!------------------------------------------------------------------

	real(r8), public, parameter :: catomw = 12.011_r8 ! molar mass of C atoms (g/mol)
	real(r8), public, parameter :: methane_atomw = 16.04_r8 ! molar mass of CH4 atoms (g/mol)

	real(r8), public :: s_con(ngases,4)    ! Schmidt # calculation constants (spp, #)
	data (s_con(1,iloop),iloop=1,4) /1898.0_r8, -110.1_r8, 2.834_r8, -0.02791_r8/ ! CH4
	data (s_con(2,iloop),iloop=1,4) /1801.0_r8, -120.1_r8, 3.7818_r8, -0.047608_r8/ ! O2
	data (s_con(3,iloop),iloop=1,4) /1911.0_r8, -113.7_r8, 2.967_r8, -0.02943_r8/ ! CO2

	real(r8), public :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (*10^-9 m2/s)
	data (d_con_w(1,iloop),iloop=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
	data (d_con_w(2,iloop),iloop=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
	data (d_con_w(3,iloop),iloop=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

	real(r8), public :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (*10^-4 m2/s)
	data (d_con_g(1,iloop),iloop=1,2) /0.1875_r8, 0.0013_r8/ ! CH4
	data (d_con_g(2,iloop),iloop=1,2) /0.1759_r8, 0.00117_r8/ ! O2
	data (d_con_g(3,iloop),iloop=1,2) /0.1325_r8, 0.0009_r8/ ! CO2

	real(r8), public :: c_h(ngases)    ! constant (K) for Henry's law (4.12, Wania)
	data c_h(1:3) /1600._r8, 1500._r8, 2400._r8/ ! CH4, O2, CO2

	real(r8), public :: kh_theta(ngases)    ! Henry's constant (mol/L/atm) at standard temperature (298K)
   data kh_theta(1:3) /1.4e-3_r8, 1.3e-3_r8, 3.4e-2_r8/ ! CH4, O2, CO2

	real(r8), public :: kh_tbase = 298.15_r8 ! base temperature for calculation of Henry's constant (K)

!------------------------------------------------------------------

	real(r8), public, parameter :: rgasm = 8.31446261815324_r8 ! Universal gas constant [J mol-1 K-1]
                                 ![J/K/mol]=[N*m/K/mol]=[Pa*m3/K/mol]
	!!! rgas Different from CoLM
	!!! rgas in CoLM is gas constant for dry air [J/kg/K]
	!!! not Universal gas constant
	real(r8), public, parameter :: rgasLatm = 0.08206_r8 ! L*atm/mol/K
   ! rgasLatm      = rgasm         / 101325 / 0.001
   ! [L*atm/mol/K] = [Pa*m3/K/mol] / [Pa/atm] / [m3/L]

	real(r8), public, parameter :: secspday = 86400._r8 ! Seconds per day

   type Methane_type
      ! ---------------------------------------------------- Variable parameter ----------------------------------------------------
      ! methane production constants
      real(r8) :: q10methane =2._r8             ! additional Q10 for methane production ABOVE the soil decomposition temperature relationship (doc:Q10 Baseline:2 Range:1.5~4) (params:2.0)
      real(r8) :: f_methane = 0.2            ! ratio of CH4 production to total C mineralization (Baseline:0.2 Range:NA params:0.2)

      ! methane oxidation constants
	      ! Oxidation vmax is per aqueous/active-water volume; methane_oxid
	      ! multiplies by vol_aqu to produce bulk-soil mol m-3 s-1 rates.
	      real(r8) :: vmax_methane_oxid = 1.25e-5_r8       ! Unit: mol m-3-aqueous s-1
	      real(r8) :: vmax_oxid_unsat = 1.25e-6_r8         ! Unit: mol m-3-aqueous s-1
      real(r8) :: k_m = 5.e-3_r8                ! Michaelis-Menten oxidation rate constant for CH4 concentration (params:5e-3 code:5.e-6_r8 * 1000._r8) (doc:KCH4 Baseline:5e-3 Range:5e-4~5e-2)
      real(r8) :: k_m_unsat = 5.e-4_r8           ! Michaelis-Menten oxidation rate constant for CH4 concentration (params:5e-4 code:5.e-6_r8 * 1000._r8 / 10._r8) (doc:KCH4 Baseline:5e-3 Range:5e-4~5e-2)
      real(r8) :: k_m_o2 =2.e-2_r8             ! Michaelis-Menten oxidation rate constant for O2 concentration (params:2e-2 code:20.e-6_r8 * 1000._r8) (doc:KO2 Baseline:2e-2 Range:2e-3~2e-1)
      real(r8) :: q10_methane_oxid = 1.9_r8         ! Q10 oxidation constant (? params:1.9)

      ! optional microbial-pool dynamics (default off; flux impact only when override is explicitly enabled)
      logical  :: use_microbial_pools = .false.
      logical  :: use_microbial_flux_override = .false.
      logical  :: use_microbial_dormancy = .true.
      real(r8) :: B_init_methanogen = 1.0_r8
      real(r8) :: B_init_methanotroph = 1.0_r8
      real(r8) :: B_min_methanogen = 1.0e-3_r8
      real(r8) :: B_min_methanotroph = 1.0e-3_r8
      real(r8) :: mu_max_methanogen = 0.2_r8
      real(r8) :: mu_max_methanotroph = 0.5_r8
      real(r8) :: gamma_methanogen = 0.05_r8
      real(r8) :: gamma_methanotroph = 0.10_r8
      real(r8) :: gamma_microbial_dormant = 0.005_r8
      real(r8) :: gamma_microbial_freeze = 0.001_r8
      ! DEPRECATED: legacy rate-based substrate half-saturation [mol C m-3 s-1].
      ! Current microbial-pool code uses K_substrate_methanogen_pool instead.
      real(r8) :: K_substrate_methanogen_rate = 1.0e-7_r8
      real(r8) :: K_substrate_methanogen_pool = 0.04_r8   ! substrate-pool half-saturation [mol C m-3]
      real(r8) :: K_inh_O2_methanogen = 1.0_r8
      real(r8) :: kappa_m_methanogen = 1.0e-2_r8
      real(r8) :: kappa_m_methanotroph = 1.0e-2_r8
      real(r8) :: q10_microbe_growth = 3.0_r8
      real(r8) :: T_ref_microbe = 298.15_r8
      real(r8) :: dormancy_rate_active = 0.1_r8
      real(r8) :: dormancy_rate_revive = 1.0_r8
      ! DEPRECATED: legacy shared dormancy threshold; split fields below are used.
      real(r8) :: dormancy_threshold_f_S = 0.1_r8
      real(r8) :: dormancy_threshold_methanogen_fS = 0.1_r8
      real(r8) :: dormancy_threshold_methanogen_fO2 = 0.1_r8
      real(r8) :: dormancy_threshold_methanotroph_fS = 0.1_r8
      real(r8) :: dormancy_threshold_methanotroph_fO2 = 0.1_r8

      ! methane ebullition constants
      real(r8) :: vgc_max  =0.15_r8            ! gas-volume fraction scale for ebullition threshold/target (unitless)

      ! methane aerenchyma constants
      real(r8) :: nongrassporosratio = 1._r8/3._r8   ! Ratio of root porosity in non-grass to grass, used for aerenchyma transport (params:0.33)
      real(r8) :: poros_tiller = 0.3_r8    ! Porosity for grass tiller
      real(r8) :: unsat_aere_ratio = 0.05_r8/0.3_r8    ! Ratio to multiply upland vegetation aerenchyma porosity by compared to inundated systems (params:0.1666666667 code:0.05_r8 / 0.3_r8)
      real(r8) :: porosmin = 0.05_r8            ! minimum aerenchyma porosity (unitless)(params:0.05 code:0.05_r8)
      real(r8) :: aere_radius = 2.9e-3_r8 ! Aerenchyma radius
      real(r8) :: rob = 3._r8                 ! ratio of root length to vertical depth ("root obliquity") (params:3. code:3._r8)
      real(r8) :: scale_factor_aere = 1._r8   ! scale factor on the aerenchyma area for sensitivity tests (1 params:1.) (doc:Fa Baseline:1 Range:0.5~1.5)

      ! methane transport constants
      real(r8) :: scale_factor_gasdiff = 1._r8   ! For sensitivity tests; convection would allow this to be > 1(? params:1.) (doc:fD0 Basline:1 Range:1,10 Unit:m2 s-1)
      real(r8) :: scale_factor_liqdiff = 1._r8   ! For sensitivity tests; convection would allow this to be > 1(? params:1.) (doc:fD0 Basline:1 Range:1,10 Unit:m2 s-1)

      ! -------------------------------------------------- Invariant parameter ----------------------------------------------------
      ! methane production constants
      real(r8) :: mino2lim = 0.2_r8         ! minimum anaerobic decomposition rate as a fraction of potential aerobic rate (0.2+ params:0.2)
      real(r8) :: q10methane_base = 295._r8 ! temperature at which the effective f_methane actually equals the constant f_methane (295+ params:295)
      real(r8) :: q10lake =3._r8       ! lake sediment decomposition Q10 used directly by methane physics.
      real(r8) :: q10lakebase = 298._r8       ! (K) base temperature for lake CH4 production (params:298. code:298._r8)
      real(r8) :: cnscalefactor=1.        ! scale factor on CN decomposition for assigning methane flux (?- params:1.)

      real(r8) :: redoxlag =30.           ! Number of days to lag in the calculation of finundated_lag (30+ params:30.)
      real(r8) :: lake_decomp_fact =9e-11    ! Base decomposition rate (1/s) at 25C (1 params:9e-11)
      real(r8) :: redoxlag_vertical=30._r8   ! time lag (days) to inhibit production for newly unsaturated layers (30+ params:0.)
      real(r8) :: pHmax = 9._r8          ! maximum pH for methane production(params:9. code:9.)
      real(r8) :: pHmin = 2.2_r8         ! minimum pH for methane production(params:2.2 code:2.2)
      real(r8) :: oxinhib = 400._r8          ! inhibition of methane production by oxygen (m^3/mol) (400+? params:400.)

      ! methane oxidation constants
      real(r8) :: smp_crit =-2.4e5_r8            ! Critical soil moisture potential (mm) (params:-2.4e5)

      ! methane ebullition constants
      real(r8) :: bubble_f =0.57_r8            ! CH4 content in gas bubbles (Kellner et al. 2006)

      ! methane aerenchyma constants
      real(r8) :: aereoxid =0._r8            ! fraction of methane flux entering aerenchyma rhizosphere that will be(? params:0.)
      real(r8) :: wet_lai = 4._r8 ! Legacy namelist field retained for compatibility; not used by current physics.
      real(r8) :: tiller_C = 0.22_r8 ! Per tiller 0.22 g C [g C/tiller]

      ! methane transport constants
      real(r8) :: satpow  =2._r8             ! exponent on watsat for saturated soil solute diffusion (2? params:2.)
      real(r8) :: capthick = 100._r8         ! min thickness before assuming h2osfc is impermeable (mm) (params:100.code:100._r8)

      ! additional constants
      real(r8) :: atm_methane  = 1.7e-6_r8         ! Atmospheric CH4 mixing ratio to prescribe if not provided by the atmospheric model (params:1.7e-6 code:1.7e-6_r8 search:1.9e-6) (mol/mol) could change with year
      real(r8) :: om_frac_sf = 1._r8          ! Scale factor for organic matter fraction (unitless)(? params:NA)

      ! -------------------------------------------------- Switch ----------------------------------------------------
      logical :: use_aereoxid_prog = .true. ! if false then aereoxid is read off of
      ! the parameter file and may be modifed by the user (default aereoxid on the
      ! file is 0.0).

      logical :: transpirationloss = .true. ! switch for activating CH4 loss from transpiration
                                    ! Transpiration loss assumes that the methane concentration in dissolved soil
                                    ! water remains constant through the plant and is released when the water evaporates
                                    ! from the stomata.
                                    ! Currently hard-wired to true; impact is < 1 Tg CH4/yr

      logical :: allowlakeprod = .false. ! Switch to allow production under lakes based on soil carbon dataset
                                 ! (Methane can be produced, and CO2 produced from methane oxidation,
                                 ! which will slowly reduce the available carbon stock, if ! replenishlakec, but no other biogeochem is done.)
                                 ! Note: switching this off turns off ALL lake methane biogeochem. However, 0 values
                                 ! will still be averaged into the concentration _sat history fields.

	      logical :: usephfact = .false. ! Switch to use pH factor in methane production; defaults to neutral pH if no spatial pH field is wired.

      logical :: replenishlakec = .true. ! Switch for keeping carbon storage under lakes constant
                                    ! so that lakes do not affect the carbon balance
                                    ! Good for long term rather than transient warming experiments
            ! NOTE SWITCHING THIS OFF ASSUMES TRANSIENT CARBON SUPPLY FROM LAKES; COUPLED MODEL WILL NOT CONSERVE CARBON
            ! IN THIS MODE.

      logical :: methane_offline = .true.    ! true --> Methane is not passed between the land & atmosphere.
                                 ! NEM is not added to NEE flux to atm. to correct for methane production,
                                 ! and ambient CH4 is set to constant 2009 value.

      logical :: methane_rmcnlim = .false.   ! Remove the N and low moisture limitations on SOM HR when calculating
                                 ! methanogenesis.
                                 ! Note: this option has not been extensively tested.
                                 ! Currently hardwired off.

      logical :: anoxicmicrosites = .false. ! Use Arah & Stephen 1998 expression to allow production above the water table
                                    ! Currently hardwired off; expression is crude.

      logical :: methane_frzout = .false.    ! Exclude CH4 from frozen fraction of soil pore H2O, to simulate "freeze-out" pulse
                                 ! as in Mastepanov 2008.
                                 ! Causes slight increase in emissions in the fall and decrease in the spring.
                                 ! Currently hardwired off; small impact.

      ! public :: methane_conrd ! Read and initialize CH4 constants

      logical :: use_nitrif_denitrif = .true.

      logical :: anoxia  = .true. ! true => anoxia is applied to heterotrophic respiration also considered in CH4 model
                                    ! default value reset in controlMod
                                    ! Whether to enable the anoxia for the seasonally induatded zones
      logical :: use_vertical_redoxlag = .true. ! Whether to enable the vertical redox lag effect

      ! SIF is enabled only when coupled BGC does not already apply anoxia limits.
      logical :: bgc_anoxia_limits_decomp = .false.  ! true: BGC o_scalar already limits decomposition
      logical :: use_ch4_sif              = .true.   ! true: apply CH4 seasonal inundation factor
      ! If BGC later applies real o_scalar limits, set (true, false).
   END type Methane_type

   type Methane_hydrology_type
            real(r8) :: vdcf = 2._r8
      real(r8) :: slopebeta = -3._r8
      real(r8) :: slopemax = 0.4_r8
      real(r8) :: pc = 0.4_r8
   END type Methane_hydrology_type

   type (Methane_type) :: DEF_METHANE
   type (Methane_hydrology_type) :: DEF_METHANE_hydrology

CONTAINS

   SUBROUTINE read_methane_namelist (nlfile)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   IMPLICIT NONE

   character(len=*), intent(in) :: nlfile
   ! Local variables
   logical :: fexists
   integer :: ierr
   integer :: unit_nml

   namelist /nl_colm_methane_parameter/ DEF_METHANE,DEF_METHANE_hydrology

      ! Read on every rank so all workers use the same methane constants.
      ! The previous master-only read left non-master ranks at default values
      ! whenever users overrode DEF_METHANE in the parameter namelist.
      INQUIRE (file=trim(nlfile), exist=fexists)
      IF (.not. fexists) THEN
         CALL CoLM_Stop (' ***** ERROR: methane parameter file does not exist: '// trim(nlfile))
      ENDIF
      open(newunit=unit_nml, status='OLD', file=trim(nlfile), form="FORMATTED")
      read(unit_nml, nml=nl_colm_methane_parameter, iostat=ierr)
      IF (ierr /= 0) THEN
         close(unit_nml)
         CALL CoLM_Stop (' ***** ERROR: Problem reading namelist: '// trim(nlfile))
      ENDIF
      close(unit_nml)

   END SUBROUTINE read_methane_namelist
END MODULE MOD_Tracer_Methane_Const
#endif
