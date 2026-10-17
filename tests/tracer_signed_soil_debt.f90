#include <define.h>
program tracer_signed_soil_debt
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
#ifdef USEMPI
   use MOD_SPMD_Task, only: spmd_init, spmd_exit
#endif
   use MOD_Precision
   use MOD_Tracer_Defs
   use MOD_Tracer_Vars
   use MOD_Tracer_Conservation, only: tracer_save_storage
   use MOD_Tracer_SoilWater
   use MOD_Namelist, only: DEF_TRACER_SOIL_DIFFUSION, DEF_TRACER_SOIL_VAPOR_DIFFUSION, &
      DEF_TRACER_USE_FRACTIONATION
   implicit none
   real(r8), allocatable :: water(:), zero(:), ice(:)
   real(r8) :: initial_mass(2), flux
   integer :: i, cycle, hold, layers
   character(len=32) :: mode

#ifdef USEMPI
   call spmd_init()
#endif
   call get_command_argument(1, mode)
   if (len_trim(mode) == 0) mode = 'uniform'
   ntracers = 2
   allocate(tracers(ntracers))
   do i = 1, ntracers
      tracers(i)%name = 'signed-debt-test'
      tracers(i)%category = 'conservative'
      tracers(i)%family_id = FAMILY_SOLUTE
      tracers(i)%state_owner = STATE_OWNER_GENERIC_WATER
      tracers(i)%ref_ratio = 1._r8
      tracers(i)%init_delta = 0._r8
      tracers(i)%init_conc = 1._r8
      tracers(i)%precip_default_conc = 1._r8
      tracers(i)%vapor_default_conc = 1._r8
      tracers(i)%mol_weight = 18._r8
      tracers(i)%reactive_decay_rate = 0._r8
   enddo
   tracers(1)%family_id = FAMILY_ISOTOPE
   DEF_TRACER_SOIL_DIFFUSION = .false.
   DEF_TRACER_SOIL_VAPOR_DIFFUSION = .false.
   DEF_TRACER_USE_FRACTIONATION = .false.
   layers = 2
   if (mode == 'throughdown' .or. mode == 'throughup' .or. &
       mode == 'gradientdown' .or. mode == 'gradientup') layers = 3
   allocate(water(layers), zero(layers), ice(layers))
   call allocate_Tracer_Vars(1, 0, layers)
   zero = 0._r8
   ice = 0._r8
   if (layers == 2) then
      water = [.01_r8, 10._r8]
   else
      water = [1._r8, .01_r8, 10._r8]
      if (mode == 'throughup' .or. mode == 'gradientup') water = [10._r8, .01_r8, 1._r8]
   endif
   if (mode == 'reverse') water = [10._r8, .01_r8]
   if (mode == 'atlimit') water(1) = -.00025_r8
   if (mode == 'enriched' .or. mode == 'overflow') water(1) = -.000125_r8
   do i = 1, ntracers
      trc_wliq_soisno(i,:,1) = water
      if (mode == 'gradientdown') trc_wliq_soisno(i,1,1) = 2._r8*water(1)
      if (mode == 'gradientup') trc_wliq_soisno(i,3,1) = 2._r8*water(3)
      if (mode == 'enriched') trc_wliq_soisno(i,2,1) = 2._r8*water(2)
      if (mode == 'different') trc_wliq_soisno(i,2,1) = .01_r8
      initial_mass(i) = sum(trc_wliq_soisno(i,:,1))
   enddo

   select case (trim(mode))
   case ('uniform', 'nondebt')
      do cycle = 1, 4
         flux = water(1) + .000125_r8
         if (mode == 'nondebt') flux = .005_r8
         call advance(flux, .true.)
         do hold = 1, 3
            call advance(0._r8, .true.)
         enddo
         call advance(-.01_r8, .true.)
      enddo
   case ('gradientdown', 'gradientup')
      if (mode == 'gradientdown') then
         call advance(.02_r8, .false., .030125_r8)
      else
         call advance(-.030125_r8, .false., -.02_r8)
      endif
      call check_gradient()
      do hold = 1, 3
         call advance(0._r8, .false.)
         call check_gradient()
      enddo
   case ('throughdown', 'throughup')
      ! The middle layer passes more than its snapshot stock and needs upstream mixing.
      if (mode == 'throughdown') then
         call advance(.02_r8, .true., .030125_r8)
      else
         call advance(-.030125_r8, .true., -.02_r8)
      endif
      do hold = 1, 3
         call advance(0._r8, .true.)
      enddo
      call advance(0._r8, .true., -.01_r8)
      if (water(2) <= 0._r8) error stop 'middle layer did not rewet'
   case ('lifecycle')
      call advance(.010125_r8, .false.)
      call tracer_save_storage(1, 0, layers)
      error stop 'full timestep accepted experimental debt'
   case ('nan-coupling')
      call advance(.010125_r8, .false.)
      error stop 'nonfinite coupling was not rejected'
   case ('atlimit')
      call advance(0._r8, .true.)
      call advance(-.01_r8, .true.)
   case ('partial')
      call advance(.010125_r8, .true.)
      call advance(-.0000625_r8, .true.)
      call advance(0._r8, .true.)
      call advance(-.01_r8, .true.)
   case ('reverse')
      call advance(-.010125_r8, .true.)
      call advance(0._r8, .true.)
      call advance(.01_r8, .true.)
   case ('enriched')
      ! Positive repayment must not exceed the initial concentration range [1,2].
      call advance(-.0002_r8, .false.)
      error stop 'enriched repayment overshoot was not rejected'
   case ('nan')
      trc_wliq_soisno(1,1,1) = ieee_value(0._r8, ieee_quiet_nan)
      call advance(.010125_r8, .false.)
      error stop 'nonfinite inventory was not rejected'
   case ('overflow')
      trc_wliq_soisno(1,1,1) = -huge(1._r8)/2._r8
      call advance(0._r8, .false.)
      error stop 'ratio overflow was not rejected'
   case ('phase', 'diffusion', 'reaction', 'provider', 'solubility')
      if (mode == 'phase') ice(1) = 1._r8
      if (mode == 'diffusion') DEF_TRACER_SOIL_DIFFUSION = .true.
      if (mode == 'reaction') tracers(2)%reaction_mode = REACTION_FIRST_ORDER
      if (mode == 'provider') tracers(2)%state_owner = STATE_OWNER_PROVIDER
      if (mode == 'solubility') tracers(2)%max_dissolved_conc = 2._r8
      call advance(.010125_r8, .false.)
      error stop 'unsupported process was not rejected'
   case ('different')
      call advance(.010125_r8, .false.)
      call advance(0._r8, .false.)
      ! A low-concentration refill cannot pay the existing signed tracer debt.
      ! A supported implementation must reject, not return positive water / negative mass.
      call advance(-.0002_r8, .false.)
      error stop 'unsafe different-concentration refill was not rejected'
   case ('overlimit')
      ! With two soil layers, solver tol_z is .001/(2*2) = .00025 mm.
      call advance(.0105_r8, .false.)
      error stop 'overlimit soil water debt was not rejected'
   case default
      error stop 'unknown case'
   end select
   print *, 'PASS ', trim(mode)
   call deallocate_Tracer_Vars()
#ifdef USEMPI
   call spmd_exit()
#endif
contains
   subroutine check_gradient()
      real(r8) :: expected(3)
      integer :: species
      expected = [1.96_r8, -.000125_r8*(5._r8/3._r8), 10._r8+.030125_r8*(5._r8/3._r8)]
      if (mode == 'gradientup') expected = expected(3:1:-1)
      if (abs(water(2)+.000125_r8) > 1.e-12_r8) error stop 'gradient middle water'
      do species = 1, ntracers
         if (any(abs(trc_wliq_soisno(species,:,1)-expected) > 1.e-11_r8)) &
            error stop 'gradient throughflow mixing mass'
         if (abs(trc_wliq_soisno(species,2,1)/water(2)-5._r8/3._r8) > 1.e-9_r8) &
            error stop 'gradient throughflow mixing ratio'
      enddo
   end subroutine check_gradient

   subroutine advance(face_flux, uniform, second_face)
      real(r8), intent(in) :: face_flux
      logical, intent(in) :: uniform
      real(r8), intent(in), optional :: second_face
      real(r8) :: after(layers), q(0:layers), mass(layers), infiltration
      integer :: species
      q = 0._r8
      q(1) = face_flux
      if (present(second_face)) q(2) = second_face
      after = water + q(0:layers-1) - q(1:layers)
      infiltration = 0._r8
      if (mode == 'nan-coupling') infiltration = ieee_value(0._r8, ieee_quiet_nan)
      call tracer_soil_water(ipatch=1, deltim=1._r8, snl=0, nl_soil=layers, &
         qlayer=q, qinfl=infiltration, qcharge=0._r8, rsur=0._r8, rsub=0._r8, &
         qseva_in=0._r8, qsdew_in=0._r8, qsubl_in=0._r8, qfros_in=0._r8, &
         qseva_soil=0._r8, qsdew_soil=0._r8, qsubl_soil=0._r8, qfros_soil=0._r8, &
         qseva_snow=0._r8, qsdew_snow=0._r8, qsubl_snow=0._r8, qfros_snow=0._r8, &
         sm=0._r8, fsno=0._r8, split_soilsnow=.false., &
         wliq_soisno=after, wice_soisno=ice, wliq_soisno_bef=water, wice_soisno_bef=ice, &
         wa=0._r8, wa_bef=0._r8, wdsrf=0._r8, wdsrf_bef=0._r8, &
         wetwat=0._r8, wetwat_bef=0._r8, pg_rain=0._r8, pg_snow=0._r8, &
         wblc_ice_sink=zero, etroot_actual=zero, etroot_aquifer=0._r8, qflx_irrig_ground=0._r8, &
         qgtop_solver=0._r8, permeable_soil=spread(.true., 1, layers))
      do species = 1, ntracers
         mass = trc_wliq_soisno(species,:,1)
         write(*,*) &
            'case=', trim(mode), ' species=', species, ' water=', after, ' signed_mass=', mass, &
            ' balance/ledger/waterledger=', sum(mass)-initial_mass(species), &
            trc_numerical_residual_step(species,1), trc_numerical_water_step(species,1)
         if (.not.all(ieee_is_finite(mass))) error stop 'nonfinite tracer mass'
         if (.not.all(ieee_is_finite(after))) error stop 'nonfinite water'
         if (.not.ieee_is_finite(trc_numerical_residual_step(species,1))) error stop 'nonfinite residual'
         if (.not.ieee_is_finite(trc_numerical_water_step(species,1))) error stop 'nonfinite water ledger'
         ! Signed inventory must close without phantom compensation or clipping.
         if (abs(sum(mass)-initial_mass(species)) > 1.e-11_r8) error stop 'signed inventory not conserved'
         if (abs(trc_numerical_residual_step(species,1)) > 1.e-11_r8) error stop 'phantom numerical mass'
         if (abs(trc_numerical_water_step(species,1)) > 1.e-11_r8) error stop 'phantom numerical water'
         if (uniform) then
            if (any(abs(mass-after) > 1.e-11_r8)) error stop 'uniform signature changed'
         endif
         if (any(after > 1.e-12_r8 .and. mass < -1.e-12_r8)) &
            error stop 'unphysical positive carrier with negative tracer'
      enddo
      water = after
   end subroutine advance
end program tracer_signed_soil_debt
