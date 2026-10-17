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
   real(r8), allocatable :: water(:), zero(:), ice(:), profile(:)
   real(r8) :: initial_mass(2), base(2), lowest(2), highest(2), expected(2), flux, debt
   integer :: i, cycle, hold, layers
   integer(8) :: draw_state
   logical :: check_expected
   character(len=32) :: mode

#ifdef USEMPI
   call spmd_init()
#endif
   call get_command_argument(1, mode)
   if (len_trim(mode) == 0) mode = 'uniform'
   ! Species 1 has an exact unit ratio; species 2 has an O18-like ratio, so rounding is exercised.
   base = [1._r8, 2.0052e-3_r8]
   if (mode == 'drainlarge') base = [.37_r8, 50._r8]
   ! Finite, but ratio * water overflows unless the kernel rejects it up front.
   if (mode == 'hugeratio') base(1) = huge(1._r8)/16._r8
   ntracers = 2
   allocate(tracers(ntracers))
   do i = 1, ntracers
      tracers(i)%name = 'signed-debt-test'
      tracers(i)%category = 'conservative'
      tracers(i)%family_id = FAMILY_ISOTOPE
      tracers(i)%state_owner = STATE_OWNER_GENERIC_WATER
      tracers(i)%ref_ratio = base(i)
      tracers(i)%init_delta = 0._r8
      tracers(i)%init_conc = base(i)
      tracers(i)%precip_default_conc = base(i)
      tracers(i)%vapor_default_conc = base(i)
      tracers(i)%mol_weight = 18._r8
      tracers(i)%reactive_decay_rate = 0._r8
   enddo
   if (mode == 'solute' .or. mode == 'solubility') tracers(2)%family_id = FAMILY_SOLUTE
   DEF_TRACER_SOIL_DIFFUSION = .false.
   DEF_TRACER_SOIL_VAPOR_DIFFUSION = .false.
   DEF_TRACER_USE_FRACTIONATION = .false.
   layers = 2
   select case (trim(mode))
   case ('throughdown', 'throughup', 'gradientdown', 'gradientup', 'debtthrough', &
         'emptybackflow', 'convergent', 'randomwalk')
      layers = 3
   end select
   allocate(water(layers), zero(layers), ice(layers), profile(layers))
   call allocate_Tracer_Vars(1, 0, layers)
   zero = 0._r8
   ice = 0._r8
   profile = 1._r8
   if (layers == 2) then
      water = [.01_r8, 10._r8]
   else
      water = [1._r8, .01_r8, 10._r8]
      if (mode == 'throughup' .or. mode == 'gradientup') water = [10._r8, .01_r8, 1._r8]
   endif
   select case (trim(mode))
   case ('reverse')
      water = [10._r8, .01_r8]
   case ('atlimit')
      water(1) = -.00025_r8
   case ('enriched', 'overflow')
      water(1) = -.000125_r8
   case ('hugeratio')
      water = [-.0001_r8, 10._r8]
   case ('partialsettle')
      water(1) = -.0002_r8
   case ('drainlarge')
      water = [10._r8, 10._r8]
   case ('debtthrough', 'convergent')
      water(2) = -.0001_r8
   case ('emptybackflow')
      water(1) = 0._r8
      water(2) = 1._r8
   case ('randomwalk')
      water = [.0004_r8, .0004_r8, 5._r8]
   end select
   select case (trim(mode))
   case ('gradientdown')
      profile(1) = 2._r8
   case ('gradientup')
      profile(3) = 2._r8
   case ('enriched', 'partialsettle')
      profile(2) = 2._r8
   case ('different')
      profile(2) = .001_r8
   case ('debtthrough')
      profile = [2._r8, 1.5_r8, 1.5_r8]
   case ('convergent', 'randomwalk')
      profile = [2._r8, 1.5_r8, 1._r8]
   case ('emptybackflow')
      profile = [1._r8, 2._r8, 1._r8]
   end select
   do i = 1, ntracers
      trc_wliq_soisno(i,:,1) = water*base(i)*profile
      initial_mass(i) = sum(trc_wliq_soisno(i,:,1))
      lowest(i) = base(i)*minval(profile, mask=abs(water) > trc_water_min_for_ratio)
      highest(i) = base(i)*maxval(profile, mask=abs(water) > trc_water_min_for_ratio)
   enddo
   check_expected = .false.
   expected = 0._r8

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
   case ('repayexact', 'repaynear')
      ! Repaying to (almost) zero water divides rounding by a tiny carrier.
      call advance(.010125_r8, .true.)
      debt = water(1)
      if (mode == 'repayexact') then
         call advance(debt, .true.)
         if (water(1) /= 0._r8) error stop 'exact repayment left water'
         if (any(trc_wliq_soisno(:,1,1) /= 0._r8)) error stop 'exact repayment left tracer'
      else
         call advance(debt*(1._r8+1.e-6_r8), .true.)
      endif
      call advance(-.01_r8, .true.)
   case ('drainlarge')
      ! Draining a large stock into debt cancels most of its mass.
      call advance(water(1) + .000125_r8, .true.)
      call advance(0._r8, .true.)
      call advance(-.01_r8, .true.)
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
      ! Richer water repays a debt: the remainder keeps the supplier's ratio, the excess is booked.
      call advance(-.0002_r8, .false.)
      check_expected = .true.
      expected = -.000125_r8*(2._r8-1._r8)*base
      call check_ratio(1, 2._r8)
   case ('different')
      call advance(.010125_r8, .false.)
      call advance(0._r8, .false.)
      ! Leaner water repays a richer debt; the shortfall is booked as a numerical source.
      check_expected = .true.
      do i = 1, ntracers
         expected(i) = water(1)*(trc_wliq_soisno(i,2,1)/water(2) - trc_wliq_soisno(i,1,1)/water(1))
      enddo
      call advance(-.0002_r8, .false.)
      call check_ratio(1, trc_wliq_soisno(1,2,1)/water(2)/base(1))
   case ('partialsettle')
      ! Inflow smaller than the debt: the remaining debt keeps its own ratio.
      call advance(-.0001_r8, .false.)
      if (abs(water(1)+.0001_r8) > 1.e-15_r8) error stop 'partial settlement water'
      check_expected = .true.
      expected = -.0001_r8*(2._r8-1._r8)*base
      call check_ratio(1, 1._r8)
   case ('debtthrough')
      ! A debt layer in throughflow repays first, then passes pure supplier water on.
      call advance(.0003_r8, .false., .00015_r8)
      check_expected = .true.
      expected = -.0001_r8*(2._r8-1.5_r8)*base
      call check_ratio(2, 2._r8)
   case ('convergent')
      ! Both neighbours repay one debt; it settles at their volume-weighted ratio.
      call advance(.0001_r8, .false., -.0002_r8)
      check_expected = .true.
      expected = -.0001_r8*(4._r8/3._r8-1.5_r8)*base
      call check_ratio(2, 4._r8/3._r8)
   case ('emptybackflow')
      ! An empty layer drains, then the same receiver refills it.
      call advance(.0001_r8, .false.)
      call check_ratio(1, 2._r8)
      call check_ratio(2, 2._r8)
      call advance(-.0002_r8, .false.)
      call check_ratio(1, 2._r8)
   case ('randomwalk')
      draw_state = 20261017_8
      do cycle = 1, 400
         call random_faces()
      enddo
   case ('nan')
      trc_wliq_soisno(1,1,1) = ieee_value(0._r8, ieee_quiet_nan)
      call advance(.010125_r8, .false.)
      error stop 'nonfinite inventory was not rejected'
   case ('overflow')
      trc_wliq_soisno(1,1,1) = -huge(1._r8)/2._r8
      call advance(0._r8, .false.)
      error stop 'ratio overflow was not rejected'
   case ('hugeratio')
      call advance(-5._r8, .false.)
      error stop 'ratio product overflow was not rejected'
   case ('phase', 'diffusion', 'reaction', 'provider', 'solubility', 'solute')
      if (mode == 'phase') ice(1) = 1._r8
      if (mode == 'diffusion') DEF_TRACER_SOIL_DIFFUSION = .true.
      if (mode == 'reaction') tracers(2)%reaction_mode = REACTION_FIRST_ORDER
      if (mode == 'provider') tracers(2)%state_owner = STATE_OWNER_PROVIDER
      if (mode == 'solubility') tracers(2)%max_dissolved_conc = 2._r8
      call advance(.010125_r8, .false.)
      error stop 'unsupported process was not rejected'
   case ('overlimit')
      ! With two soil layers, solver tol_z is .001/(2*2) = .00025 mm.
      call advance(.0105_r8, .false.)
      error stop 'overlimit soil water debt was not rejected'
   case default
      error stop 'unknown case'
   end select
   if (check_expected) then
      do i = 1, ntracers
         if (abs(trc_numerical_residual_step(i,1)-expected(i)) > mass_tol(i)+1.e-9_r8*abs(expected(i))) &
            error stop 'settlement residual differs from expected'
      enddo
   endif
   print *, 'PASS ', trim(mode)
   call deallocate_Tracer_Vars()
#ifdef USEMPI
   call spmd_exit()
#endif
contains
   real(r8) function mass_tol(species)
      integer, intent(in) :: species
      mass_tol = 1.e-12_r8*base(species)*max(1._r8, sum(abs(water)))
   end function mass_tol

   subroutine check_ratio(level, relative)
      integer, intent(in) :: level
      real(r8), intent(in) :: relative
      integer :: species
      do species = 1, ntracers
         if (abs(trc_wliq_soisno(species,level,1)/water(level) - relative*base(species)) > &
             1.e-9_r8*relative*base(species)) error stop 'settled layer ratio'
      enddo
   end subroutine check_ratio

   subroutine check_gradient()
      real(r8) :: expected_mass(3)
      integer :: species
      expected_mass = [1.96_r8, -.000125_r8*(5._r8/3._r8), 10._r8+.030125_r8*(5._r8/3._r8)]
      if (mode == 'gradientup') expected_mass = expected_mass(3:1:-1)
      if (abs(water(2)+.000125_r8) > 1.e-12_r8) error stop 'gradient middle water'
      do species = 1, ntracers
         if (any(abs(trc_wliq_soisno(species,:,1)-expected_mass*base(species)) > mass_tol(species))) &
            error stop 'gradient throughflow mixing mass'
         if (abs(trc_wliq_soisno(species,2,1)/water(2)-5._r8/3._r8*base(species)) > 1.e-9_r8*base(species)) &
            error stop 'gradient throughflow mixing ratio'
      enddo
   end subroutine check_gradient

   subroutine random_faces()
      ! Deterministic Park-Miller draws; keep only steps whose debt stays within the solver limit.
      real(r8) :: q1, q2, limit
      integer :: attempt
      limit = 1.e-3_r8/(2._r8*real(layers,r8))
      do attempt = 1, 50
         draw_state = mod(16807_8*draw_state, 2147483647_8)
         q1 = (real(draw_state,r8)/2147483647._r8 - .5_r8)*.0012_r8
         draw_state = mod(16807_8*draw_state, 2147483647_8)
         q2 = (real(draw_state,r8)/2147483647._r8 - .5_r8)*.0012_r8
         if (water(1)-q1 >= -limit .and. water(2)+q1-q2 >= -limit .and. water(3)+q2 >= -limit) then
            call advance(q1, .false., q2)
            return
         endif
      enddo
   end subroutine random_faces

   subroutine advance(face_flux, uniform, second_face)
      real(r8), intent(in) :: face_flux
      logical, intent(in) :: uniform
      real(r8), intent(in), optional :: second_face
      real(r8) :: after(layers), q(0:layers), mass(layers), infiltration, inflow(layers)
      real(r8) :: residual_before(ntracers), step_bound(ntracers)
      integer :: species, level
      q = 0._r8
      q(1) = face_flux
      if (present(second_face)) q(2) = second_face
      after = water + q(0:layers-1) - q(1:layers)
      infiltration = 0._r8
      if (mode == 'nan-coupling') infiltration = ieee_value(0._r8, ieee_quiet_nan)
      residual_before = trc_numerical_residual_step(:,1)
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
      ! Repayment books at most the repaid water times the concentration range.
      inflow = max(q(0:layers-1), 0._r8) + max(-q(1:layers), 0._r8)
      step_bound = sum(min(inflow, max(-water, 0._r8)))*(highest-lowest)
      water = after
      do species = 1, ntracers
         mass = trc_wliq_soisno(species,:,1)
         if (mode /= 'randomwalk') write(*,*) &
            'case=', trim(mode), ' species=', species, ' water=', after, ' signed_mass=', mass, &
            ' balance/ledger/waterledger=', sum(mass)-initial_mass(species), &
            trc_numerical_residual_step(species,1), trc_numerical_water_step(species,1)
         if (.not.all(ieee_is_finite(mass))) error stop 'nonfinite tracer mass'
         if (.not.all(ieee_is_finite(after))) error stop 'nonfinite water'
         if (.not.ieee_is_finite(trc_numerical_residual_step(species,1))) error stop 'nonfinite residual'
         if (.not.ieee_is_finite(trc_numerical_water_step(species,1))) error stop 'nonfinite water ledger'
         ! Signed inventory closes once the booked settlement residual is included.
         if (abs(sum(mass)-initial_mass(species)-trc_numerical_residual_step(species,1)) > mass_tol(species)) &
            error stop 'signed inventory not conserved'
         ! Checked per call, so bookings of opposite sign cannot cancel; a call without
         ! repayment may book rounding only.
         if (abs(trc_numerical_residual_step(species,1)-residual_before(species)) > &
             step_bound(species)*(1._r8+1.e-9_r8)+mass_tol(species)) &
            error stop 'residual exceeds repaid debt times concentration range'
         if (abs(trc_numerical_water_step(species,1)) > 1.e-11_r8) error stop 'phantom numerical water'
         if (uniform) then
            if (any(abs(mass-after*base(species)) > mass_tol(species))) error stop 'uniform signature changed'
         endif
         if (any(after > 0._r8 .and. mass < 0._r8) .or. any(after < 0._r8 .and. mass > 0._r8)) &
            error stop 'carrier and tracer signs differ'
         do level = 1, layers
            if (abs(after(level)) <= trc_water_min_for_ratio) cycle
            if (mass(level)/after(level) < lowest(species)*(1._r8-1.e-12_r8) .or. &
                mass(level)/after(level) > highest(species)*(1._r8+1.e-12_r8)) &
               error stop 'ratio left the initial range'
         enddo
      enddo
   end subroutine advance
end program tracer_signed_soil_debt
