#include <define.h>

#ifdef TRACER
MODULE MOD_Tracer_Hist

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval
   USE MOD_NetCDFSerial
   USE MOD_HistGridded
#if (defined UNSTRUCTURED || defined CATCHMENT)
   USE MOD_HistVector
#endif
#ifdef SinglePoint
   USE MOD_HistSingle
#endif
   USE MOD_Tracer_Defs, only: ntracers, tracers, mass_to_delta, trc_tiny, &
      trc_delta_sanity_max, trc_flux_water_min_for_delta, trc_water_min_for_delta, &
      tracer_uses_delta_diagnostics
   USE MOD_Tracer_Vars

   IMPLICIT NONE

   character(len=10) :: HistForm = 'Gridded'

CONTAINS

   SUBROUTINE tracer_hist_accumulate (ipatch, snl, maxsnl, nl_soil, &
      ldew_rain, ldew_snow, wliq_soisno, wice_soisno, &
      wa, wdsrf, wetwat, scv)
      IMPLICIT NONE
      integer,  intent(in) :: ipatch, snl, maxsnl, nl_soil
      real(r8), intent(in) :: ldew_rain, ldew_snow
      real(r8), intent(in) :: wliq_soisno(snl+1:nl_soil)
      real(r8), intent(in) :: wice_soisno(snl+1:nl_soil)
      ! Post-step pool states paired with the existing canopy/soil/snow args.
      ! These feed four additional per-pool tracer-ratio diagnostics in
      ! MOD_Hist (f_trc_conc_wa / wdsrf / wetwat / scv). Without them the
      ! aquifer, surface water, wetland and pre-layer thin-snow δ values are
      ! only observable in the restart file — which rarely matches the
      ! history output cadence.
      real(r8), intent(in) :: wa, wdsrf, wetwat, scv

      integer :: itrc, j, jsnow

      IF (ntracers <= 0) RETURN

      ! Canopy water
      a_water_ldew(ipatch) = a_water_ldew(ipatch) + (ldew_rain + ldew_snow)
      DO itrc = 1, ntracers
         a_trc_ldew_mass(itrc, ipatch) = a_trc_ldew_mass(itrc, ipatch) &
            + trc_ldew_rain(itrc, ipatch) + trc_ldew_snow(itrc, ipatch)
      ENDDO

      ! Soil layers (1:nl_soil) — fixed dimension, safe for MPI gather
      DO j = 1, nl_soil
         a_water_soil(j, ipatch) = a_water_soil(j, ipatch) + wliq_soisno(j) + wice_soisno(j)
         DO itrc = 1, ntracers
            a_trc_soil_mass(itrc, j, ipatch) = a_trc_soil_mass(itrc, j, ipatch) &
               + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
         ENDDO
      ENDDO

      ! Snow layers — separate accumulator with fixed dimension abs(maxsnl).
      ! Snow layer index j ranges from snl+1 to 0 (negative indices).
      ! Map to positive index: jsnow = j - maxsnl  (1-based, from top to bottom)
      ! e.g., maxsnl=-5: j=-4 → jsnow=1, j=-3 → jsnow=2, ..., j=0 → jsnow=5
      IF (snl < 0) THEN
         DO j = snl + 1, 0
            jsnow = j - maxsnl  ! positive index (1 to abs(maxsnl))
            a_water_snow(jsnow, ipatch) = a_water_snow(jsnow, ipatch) &
               + wliq_soisno(j) + wice_soisno(j)
            DO itrc = 1, ntracers
               a_trc_snow_mass(itrc, jsnow, ipatch) = a_trc_snow_mass(itrc, jsnow, ipatch) &
                  + trc_wliq_soisno(itrc, j, ipatch) + trc_wice_soisno(itrc, j, ipatch)
            ENDDO
         ENDDO
      ENDIF

		      ! Aquifer: positive storage and signed debt are different diagnostics.
		      ! wa<=0 is a wetland/aquifer debt state, not a physical concentration.
		      ! Keep trc_wa state synchronized with wa, but filter sub-mm debt
		      ! from concentration diagnostics to avoid tiny-denominator deltas.
		      IF (wa > 1._r8) a_water_wa(ipatch) = a_water_wa(ipatch) + wa
		      IF (wa < -1._r8) a_water_wa_debt(ipatch) = a_water_wa_debt(ipatch) - wa
	      a_water_wdsrf (ipatch) = a_water_wdsrf (ipatch) + wdsrf
	      a_water_wetwat(ipatch) = a_water_wetwat(ipatch) + wetwat
         ! CoLM's `scv` is total snow water equivalent even after layered
         ! snow exists, but `trc_scv` only mirrors the pre-layer thin-snow
         ! pool. Once snl<0, tracer mass lives in trc_wliq/trc_wice snow
         ! layers and is diagnosed via f_trc_conc_soisno_*; do not pair
         ! that total scv water with the thin-snow tracer denominator.
         IF (snl == 0 .and. scv > trc_tiny) THEN
            a_water_scv(ipatch) = a_water_scv(ipatch) + scv
         ENDIF
	      DO itrc = 1, ntracers
	         IF (wa > 1._r8) THEN
	            a_trc_wa_mass(itrc, ipatch) = a_trc_wa_mass(itrc, ipatch) + trc_wa(itrc, ipatch)
	         ENDIF
		         IF (wa < -1._r8) THEN
		            a_trc_wa_debt_mass(itrc, ipatch) = a_trc_wa_debt_mass(itrc, ipatch) &
		               + max(-trc_wa(itrc, ipatch), 0._r8)
		         ENDIF
	         a_trc_wdsrf_mass (itrc, ipatch) = a_trc_wdsrf_mass (itrc, ipatch) + trc_wdsrf (itrc, ipatch)
         a_trc_wetwat_mass(itrc, ipatch) = a_trc_wetwat_mass(itrc, ipatch) + trc_wetwat(itrc, ipatch)
         IF (snl == 0 .and. scv > trc_tiny) THEN
            a_trc_scv_mass(itrc, ipatch) = a_trc_scv_mass(itrc, ipatch) + trc_scv(itrc, ipatch)
         ENDIF
      ENDDO
   END SUBROUTINE tracer_hist_accumulate

   SUBROUTINE tracer_hist_out (file_hist, itime_in_file, hist_form_in, &
      sumarea, filter, maxsnl, nl_soil, forcing_has_missing_value, forcmask_pch)
   USE MOD_Block
   USE MOD_DataType
   USE MOD_LandPatch, only: numpatch
   USE MOD_SPMD_Task, only: p_is_worker
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Namelist, only: DEF_hist_vars

   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: hist_form_in
   type(block_data_real8_2d), intent(inout) :: sumarea
   logical, intent(inout) :: filter(:)
   integer, intent(in) :: maxsnl, nl_soil
   logical, intent(in) :: forcing_has_missing_value
   logical, intent(in) :: forcmask_pch(:)
            real(r8), allocatable :: trc_mass_soisno (:,:)
            real(r8), allocatable :: water_soisno    (:,:)
            real(r8), allocatable :: trc_delta_flux  (:)
            real(r8), allocatable :: snowpack_mass_pat(:)
            real(r8), allocatable :: snowpack_water_pat(:)
            real(r8), allocatable :: ones_arr        (:)
            ! Per-patch (pseudo_mass, pseudo_water) pair for the
            ! leaf-water-weighted bulk delta. pseudo_water = leaf
            ! water moles, pseudo_mass = R_b * moles. Feeding these
            ! to write_history_tracer_delta_2d makes the gridded
            ! result mass_to_delta(sum(area * R_b * moles),
            ! sum(area * moles)) -- a true leaf-water-weighted mean
            ! that obeys the sum-then-divide invariant. Reusing
            ! trc_delta_flux for this would mix two different units
            ! across the same buffer.
            real(r8), allocatable :: leaf_mass_pat   (:)
            real(r8), allocatable :: leaf_water_pat  (:)
            character(len=64)  :: trc_varname
            character(len=128) :: trc_longname
            character(len=32)  :: trc_ratio_word, trc_ratio_units
            real(r8) :: delta_loc, leaf_R_loc
            integer :: itrc_loc, ip_loc, j_loc, jsnow_loc

      HistForm = hist_form_in

      IF ((p_is_worker) .and. (numpatch > 0)) THEN
         filter(:) = patchtype < 99
         IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
      ENDIF
      IF (HistForm == 'Gridded') THEN
         CALL mp2g_hist%get_sumarea (sumarea, filter)
      ENDIF

            IF (ntracers > 0) THEN
               IF (p_is_worker) THEN
                  allocate (trc_mass_soisno (maxsnl+1:nl_soil, numpatch))
                  allocate (water_soisno    (maxsnl+1:nl_soil, numpatch))
                  allocate (trc_delta_flux  (numpatch))
                  allocate (snowpack_mass_pat (numpatch))
                  allocate (snowpack_water_pat(numpatch))
                  allocate (ones_arr        (numpatch));  ones_arr = 1._r8
                  allocate (leaf_mass_pat   (numpatch))
                  allocate (leaf_water_pat  (numpatch))
               ELSE
                  allocate (trc_mass_soisno (0,0))
                  allocate (water_soisno    (0,0))
                  allocate (trc_delta_flux  (0))
                  allocate (snowpack_mass_pat (0))
                  allocate (snowpack_water_pat(0))
                  allocate (ones_arr        (0))
                  allocate (leaf_mass_pat   (0))
                  allocate (leaf_water_pat  (0))
               ENDIF

               DO itrc_loc = 1, ntracers
                  IF (tracer_uses_delta_diagnostics(itrc_loc)) THEN
                     trc_ratio_word  = 'ratio'
                     trc_ratio_units = 'R'
                  ELSE
                     trc_ratio_word  = 'concentration'
                     trc_ratio_units = 'tracer/water'
                  ENDIF

                  IF (tracer_uses_delta_diagnostics(itrc_loc)) THEN
	                  ! --- Evapotranspiration flux delta ---
	                  ! Pair a_trc_evap with gross evaporative water so the
	                  ! denominator never goes near zero from dew/frost
	                  ! cancellation (a_fevpa is the net flux; not valid as a
	                  ! delta denominator). write_history_tracer_delta_2d
	                  ! maps mass and water separately to the grid and divides
	                  ! once at the end -- preserving the
	                  ! mass_to_delta(sum(mass), sum(water)) invariant
	                  ! documented in MOD_Tracer_Defs.F90.
	                  write(trc_varname , '(A,A)')   'f_trc_delta_evap_', trim(tracers(itrc_loc)%name)
	                  write(trc_longname, '(A,A,A)') 'evapotranspiration tracer delta (', &
	                     trim(tracers(itrc_loc)%name), ')'
		                  CALL write_history_tracer_delta_2d (DEF_hist_vars%fevpa, &
		                     a_trc_evap(itrc_loc, :), a_water_evap_gross(itrc_loc, :), &
		                     tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
		                     itime_in_file, filter, trim(trc_longname), 'permil', &
		                     water_min_override = trc_flux_water_min_for_delta)

	                  write(trc_varname , '(A,A)')   'f_trc_delta_soilevap_', trim(tracers(itrc_loc)%name)
	                  write(trc_longname, '(A,A,A)') 'soil/surface evaporation tracer delta (', &
	                     trim(tracers(itrc_loc)%name), ')'
		                  CALL write_history_tracer_delta_2d (DEF_hist_vars%fevpa, &
		                     a_trc_soilevap(itrc_loc, :), a_water_soilevap(itrc_loc, :), &
		                     tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
		                     itime_in_file, filter, trim(trc_longname), 'permil', &
		                     water_min_override = trc_flux_water_min_for_delta)

	                  write(trc_varname , '(A,A)')   'f_trc_delta_canopyevap_', trim(tracers(itrc_loc)%name)
	                  write(trc_longname, '(A,A,A)') 'canopy evaporation tracer delta (', &
	                     trim(tracers(itrc_loc)%name), ')'
		                  CALL write_history_tracer_delta_2d (DEF_hist_vars%fevpa, &
		                     a_trc_canopyevap(itrc_loc, :), a_water_canopyevap(itrc_loc, :), &
		                     tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
		                     itime_in_file, filter, trim(trc_longname), 'permil', &
		                     water_min_override = trc_flux_water_min_for_delta)

	                  write(trc_varname , '(A,A)')   'f_trc_delta_subl_', trim(tracers(itrc_loc)%name)
	                  write(trc_longname, '(A,A,A)') 'sublimation tracer delta (', &
	                     trim(tracers(itrc_loc)%name), ')'
		                  CALL write_history_tracer_delta_2d (DEF_hist_vars%fevpa, &
		                     a_trc_subl(itrc_loc, :), a_water_subl(itrc_loc, :), &
		                     tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
		                     itime_in_file, filter, trim(trc_longname), 'permil', &
		                     water_min_override = trc_flux_water_min_for_delta)

	                  write(trc_varname , '(A,A)')   'f_trc_delta_wetland_evap_', trim(tracers(itrc_loc)%name)
	                  write(trc_longname, '(A,A,A)') 'wetland evaporation/sublimation tracer delta (', &
	                     trim(tracers(itrc_loc)%name), ')'
		                  CALL write_history_tracer_delta_2d (DEF_hist_vars%fevpa, &
		                     a_trc_wetland_evap(itrc_loc, :), a_water_wetland_evap(itrc_loc, :), &
		                     tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
		                     itime_in_file, filter, trim(trc_longname), 'permil', &
		                     water_min_override = trc_flux_water_min_for_delta)

	                  ! --- Transpiration-only flux delta ---
	                  ! write_history_tracer_delta_2d does sum-then-divide on
	                  ! the grid; no per-patch delta needed here.
	                  write(trc_varname , '(A,A)')   'f_trc_delta_transp_', trim(tracers(itrc_loc)%name)
	                  write(trc_longname, '(A,A,A)') 'transpiration tracer delta (', &
	                     trim(tracers(itrc_loc)%name), ')'
		                  CALL write_history_tracer_delta_2d (DEF_hist_vars%etr, &
		                     a_trc_transp(itrc_loc, :), a_water_transp(itrc_loc, :), &
		                     tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
		                     itime_in_file, filter, trim(trc_longname), 'permil', &
		                     water_min_override = trc_flux_water_min_for_delta)

	                  ! --- Transpiration source/xylem delta before NSS storage exchange ---
	                  write(trc_varname , '(A,A)')   'f_trc_delta_transp_src_', trim(tracers(itrc_loc)%name)
	                  write(trc_longname, '(A,A,A)') 'transpiration source/xylem tracer delta (', &
	                     trim(tracers(itrc_loc)%name), ')'
		                  CALL write_history_tracer_delta_2d (DEF_hist_vars%etr, &
		                     a_trc_transp_src(itrc_loc, :), a_water_transp(itrc_loc, :), &
		                     tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
		                     itime_in_file, filter, trim(trc_longname), 'permil', &
		                     water_min_override = trc_flux_water_min_for_delta)

	                  ! --- Leaf NSS state deltas ---
                  ! e-site: NSS evaporation-site state -- there is no
                  ! conserved per-patch water mass paired with this
                  ! delta (the Peclet number is a mixing factor, not
                  ! a mass), so area-weighted aggregation is the only
                  ! defensible choice. Longname documents that this
                  ! is a state diagnostic, not a mass-weighted bulk
                  ! pool delta.
                  IF (p_is_worker) THEN
                     DO ip_loc = 1, numpatch
                        delta_loc = trc_leaf_delta_e(itrc_loc, ip_loc)
                        IF (delta_loc /= spval .and. &
                            abs(delta_loc) <= trc_delta_sanity_max) THEN
                           trc_delta_flux(ip_loc) = delta_loc
                        ELSE
                           trc_delta_flux(ip_loc) = spval
                        ENDIF
                     ENDDO
                  ENDIF
                  write(trc_varname , '(A,A)')   'f_trc_leaf_delta_e_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(A,A,A)') 'leaf evaporation-site NSS delta, area-weighted state (', &
                     trim(tracers(itrc_loc)%name), ')'
                  CALL write_history_variable_2d (DEF_hist_vars%etr, &
                     trc_delta_flux, file_hist, trim(trc_varname), itime_in_file, &
                     sumarea, filter, trim(trc_longname), 'permil', acc_num = ones_arr)

                  ! b-site: bulk leaf water has a matching per-patch
                  ! mass scale (trc_leaf_water_moles, mol/m^2). Build
                  ! pseudo (mass, water) = (R_b * moles, moles) and
                  ! route through write_history_tracer_delta_2d so the
                  ! grid result is mass_to_delta(sum(area*R_b*moles),
                  ! sum(area*moles)) -- a leaf-water-weighted bulk
                  ! delta consistent with the sum-then-divide
                  ! invariant. Patches with delta_b == spval, out-of-
                  ! sanity, or zero leaf water are spval'd in mass
                  ! and 0'd in water so the downstream filter at
                  ! water > trc_water_min_for_delta drops them.
                  IF (p_is_worker) THEN
                     DO ip_loc = 1, numpatch
                        delta_loc = trc_leaf_delta_b(itrc_loc, ip_loc)
                        IF (delta_loc /= spval .and. &
                            abs(delta_loc) <= trc_delta_sanity_max .and. &
                            trc_leaf_water_moles(itrc_loc, ip_loc) > 0._r8) THEN
                           leaf_R_loc = (1._r8 + delta_loc * 1.0e-3_r8) * &
                                        tracers(itrc_loc)%ref_ratio
                           leaf_mass_pat (ip_loc) = leaf_R_loc * &
                                                    trc_leaf_water_moles(itrc_loc, ip_loc)
                           leaf_water_pat(ip_loc) = trc_leaf_water_moles(itrc_loc, ip_loc)
                        ELSE
                           leaf_mass_pat (ip_loc) = spval
                           leaf_water_pat(ip_loc) = 0._r8
                        ENDIF
                     ENDDO
                  ENDIF
                  write(trc_varname , '(A,A)')   'f_trc_leaf_delta_b_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(A,A,A)') 'bulk leaf-water NSS delta, leaf-water-mole weighted (', &
                     trim(tracers(itrc_loc)%name), ')'
                  CALL write_history_tracer_delta_2d (DEF_hist_vars%etr, &
                     leaf_mass_pat, leaf_water_pat, &
                     tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
                     itime_in_file, filter, trim(trc_longname), 'permil')
                  ENDIF

                  ! --- Canopy ratio ---
                  write(trc_varname , '(A,A)')   'f_trc_conc_ldew_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(5A)') 'canopy tracer ', trim(trc_ratio_word), &
                     ' (', trim(tracers(itrc_loc)%name), ')'
                  CALL write_history_tracer_ratio_2d (DEF_hist_vars%ldew, &
                     a_trc_ldew_mass(itrc_loc, :), a_water_ldew, &
                     file_hist, trim(trc_varname), itime_in_file, &
                     filter, trim(trc_longname), trim(trc_ratio_units))

                  ! --- Soil + snow combined ratio ---
                  ! Snow accumulator is keyed by jsnow = j - maxsnl (positive);
                  ! soil accumulator is keyed by j (1..nl_soil). Map tracer
                  ! mass and water separately, then divide on the history grid.
                  ! This avoids area-fraction dilution when a grid cell has
                  ! multiple patch fractions with different tracer signatures.
                  IF (p_is_worker) THEN
                     trc_mass_soisno = spval
                     water_soisno    = spval
                     DO ip_loc = 1, numpatch
                        DO j_loc = maxsnl+1, 0
                           jsnow_loc = j_loc - maxsnl
                           IF (a_water_snow(jsnow_loc, ip_loc) > trc_tiny) THEN
                              trc_mass_soisno(j_loc, ip_loc) = &
                                 a_trc_snow_mass(itrc_loc, jsnow_loc, ip_loc)
                              water_soisno(j_loc, ip_loc) = a_water_snow(jsnow_loc, ip_loc)
                           ENDIF
                        ENDDO
                        DO j_loc = 1, nl_soil
                           IF (a_water_soil(j_loc, ip_loc) > trc_tiny) THEN
                              trc_mass_soisno(j_loc, ip_loc) = &
                                 a_trc_soil_mass(itrc_loc, j_loc, ip_loc)
                              water_soisno(j_loc, ip_loc) = a_water_soil(j_loc, ip_loc)
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
                  write(trc_varname , '(A,A)')   'f_trc_conc_soisno_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(5A)') 'soil/snow layer tracer ', trim(trc_ratio_word), &
                     ' (', trim(tracers(itrc_loc)%name), ')'
                  CALL write_history_tracer_ratio_3d (DEF_hist_vars%wliq_soisno, &
                     trc_mass_soisno, water_soisno, file_hist, trim(trc_varname), itime_in_file, &
                     'soilsnow', maxsnl+1, nl_soil-maxsnl, filter, &
                     trim(trc_longname), trim(trc_ratio_units))

                  ! --- Total snowpack ratio/delta ---
                  ! CoLM's water-side scv is total SWE, while tracer_scv is
                  ! only the pre-layer thin-snow pool. Build a formal snowpack
                  ! diagnostic from the same pieces that own tracer mass:
                  ! layered snow accumulators plus thin-snow scv accumulators.
                  ! This keeps f_trc_conc_scv_* as a thin-snow diagnostic and
                  ! prevents the old thin-tracer / total-SWE denominator mix.
                  IF (p_is_worker) THEN
                     snowpack_mass_pat(:)  = a_trc_scv_mass(itrc_loc, :)
                     snowpack_water_pat(:) = a_water_scv(:)
                     DO ip_loc = 1, numpatch
                        DO jsnow_loc = 1, abs(maxsnl)
                           snowpack_mass_pat(ip_loc) = snowpack_mass_pat(ip_loc) &
                              + a_trc_snow_mass(itrc_loc, jsnow_loc, ip_loc)
                           snowpack_water_pat(ip_loc) = snowpack_water_pat(ip_loc) &
                              + a_water_snow(jsnow_loc, ip_loc)
                        ENDDO
                     ENDDO
                  ENDIF

                  write(trc_varname , '(A,A)')   'f_trc_conc_snowpack_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(5A)') 'total snowpack tracer ', trim(trc_ratio_word), &
                     ' (', trim(tracers(itrc_loc)%name), ')'
                  CALL write_history_tracer_ratio_2d (DEF_hist_vars%scv, &
                     snowpack_mass_pat, snowpack_water_pat, &
                     file_hist, trim(trc_varname), itime_in_file, &
                     filter, trim(trc_longname), trim(trc_ratio_units))

                  IF (tracer_uses_delta_diagnostics(itrc_loc)) THEN
                     write(trc_varname , '(A,A)')   'f_trc_delta_snowpack_', trim(tracers(itrc_loc)%name)
                     write(trc_longname, '(A,A,A)') 'total snowpack tracer delta (', &
                        trim(tracers(itrc_loc)%name), ')'
                     CALL write_history_tracer_delta_2d (DEF_hist_vars%scv, &
                        snowpack_mass_pat, snowpack_water_pat, &
                        tracers(itrc_loc)%ref_ratio, file_hist, trim(trc_varname), &
                        itime_in_file, filter, trim(trc_longname), 'permil')
                  ENDIF

                  ! --- Aquifer (wa) ratio ---
                  ! wa can be negative during wetland aquifer-debt episodes;
                  ! the tracer_wetland redistribution keeps signed arithmetic.
                  write(trc_varname , '(A,A)')   'f_trc_conc_wa_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(5A)') 'aquifer tracer ', trim(trc_ratio_word), &
                     ' (', trim(tracers(itrc_loc)%name), ')'
	                  CALL write_history_tracer_ratio_2d (DEF_hist_vars%wa, &
	                     a_trc_wa_mass(itrc_loc, :), a_water_wa, &
	                     file_hist, trim(trc_varname), itime_in_file, &
	                     filter, trim(trc_longname), trim(trc_ratio_units))

	                  write(trc_varname , '(A,A)')   'f_trc_conc_wa_debt_', trim(tracers(itrc_loc)%name)
	                  write(trc_longname, '(5A)') 'aquifer debt tracer ', trim(trc_ratio_word), &
	                     ' (', trim(tracers(itrc_loc)%name), ')'
	                  CALL write_history_tracer_ratio_2d (DEF_hist_vars%wa, &
	                     a_trc_wa_debt_mass(itrc_loc, :), a_water_wa_debt, &
	                     file_hist, trim(trc_varname), itime_in_file, &
	                     filter, trim(trc_longname), trim(trc_ratio_units))

                  ! --- Surface water (wdsrf) ratio ---
                  write(trc_varname , '(A,A)')   'f_trc_conc_wdsrf_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(5A)') 'surface water tracer ', trim(trc_ratio_word), &
                     ' (', trim(tracers(itrc_loc)%name), ')'
                  CALL write_history_tracer_ratio_2d (DEF_hist_vars%wdsrf, &
                     a_trc_wdsrf_mass(itrc_loc, :), a_water_wdsrf, &
                     file_hist, trim(trc_varname), itime_in_file, &
                     filter, trim(trc_longname), trim(trc_ratio_units))

                  ! --- Wetland pool (wetwat) ratio ---
                  write(trc_varname , '(A,A)')   'f_trc_conc_wetwat_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(5A)') 'wetland pool tracer ', trim(trc_ratio_word), &
                     ' (', trim(tracers(itrc_loc)%name), ')'
                  IF (p_is_worker) THEN
                     IF (numpatch > 0) THEN
                        filter(:) = patchtype == 2
                        IF (forcing_has_missing_value) THEN
                           filter = filter .and. forcmask_pch
                        ENDIF
                        filter = filter .and. patchmask
                     ENDIF
                  ENDIF
                  CALL write_history_tracer_ratio_2d (DEF_hist_vars%wetwat, &
                     a_trc_wetwat_mass(itrc_loc, :), a_water_wetwat, &
                     file_hist, trim(trc_varname), itime_in_file, &
                     filter, trim(trc_longname), trim(trc_ratio_units))

                  IF (p_is_worker) THEN
                     IF (numpatch > 0) THEN
                        filter(:) = patchtype <= 3
                        IF (forcing_has_missing_value) THEN
                           filter = filter .and. forcmask_pch
                        ENDIF
                        filter = filter .and. patchmask
                     ENDIF
                  ENDIF

                  ! --- Pre-layer thin snow (scv) ratio ---
                  ! Only meaningful when snl==0 (accumulated snow below
                  ! layer-creation threshold). CoLM's water-side scv remains
                  ! total SWE after layered snow exists, but tracer mass then
                  ! lives in trc_wliq/trc_wice snow layers; the accumulator
                  ! therefore filters out snl<0 states so this ratio cannot
                  ! divide thin-snow tracer by total layered-snow SWE.
                  write(trc_varname , '(A,A)')   'f_trc_conc_scv_', trim(tracers(itrc_loc)%name)
                  write(trc_longname, '(5A)') 'thin-snow (scv) tracer ', trim(trc_ratio_word), &
                     ' (', trim(tracers(itrc_loc)%name), ')'
                  CALL write_history_tracer_ratio_2d (DEF_hist_vars%scv, &
                     a_trc_scv_mass(itrc_loc, :), a_water_scv, &
                     file_hist, trim(trc_varname), itime_in_file, &
                     filter, trim(trc_longname), trim(trc_ratio_units))

               ENDDO

               deallocate (trc_mass_soisno, water_soisno, trc_delta_flux, &
                           snowpack_mass_pat, snowpack_water_pat, ones_arr, &
                           leaf_mass_pat, leaf_water_pat)
            ENDIF

#if defined BGC
      CALL tracer_methane_hist_out (file_hist, itime_in_file, sumarea, filter, &
         nl_soil, forcing_has_missing_value, forcmask_pch)
#endif

   END SUBROUTINE tracer_hist_out


#if defined BGC
   SUBROUTINE tracer_methane_hist_out (file_hist, itime_in_file, sumarea, filter, &
      nl_soil, forcing_has_missing_value, forcmask_pch)
   USE MOD_DataType
   USE MOD_LandPatch, only: numpatch
   USE MOD_Namelist, only: DEF_METHANE_enable_rice_paddy, DEF_METHANE_only_wetland
   USE MOD_SPMD_Task, only: p_is_worker
   USE MOD_Vars_TimeInvariants, only: patchtype, patchmask
   USE MOD_Tracer_Methane_BgcLink, only: methane_patch_active_mask
            USE MOD_Tracer_Methane_Registry, only: igas_ch4
            USE MOD_Tracer_Methane_Const,    only: mhist_on => methane_history_enabled
            USE MOD_Tracer_Methane_AccFlux,  only: &
                 a_net_methane, a_methane_prod_depth, a_o2_decomp_depth, a_co2_decomp_depth, &
                 a_methane_oxid_depth, a_o2_oxid_depth, a_co2_oxid_depth, &
                 a_methane_aere_depth, a_methane_tran_depth, a_o2_aere_depth, a_co2_aere_depth, &
                 a_methane_ebul_depth, a_o2stress, a_methane_stress, &
                 a_methane_surf_flux_tot, a_methane_surf_aere, a_methane_surf_ebul, &
	                 a_methane_surf_diff, a_methane_surf_diff_phys, a_methane_balance_residual, &
	                 a_o2_cap_loss, a_o2_cap_gain, &
	                 a_methane_ebul_tot, a_methane_prod_tot, a_methane_oxid_tot, &
                 a_co2_decomp_tot, a_co2_oxid_tot, a_co2_aere_tot, a_co2_net_tot, &
                 a_totcol_methane, a_grnd_methane_cond, a_conc_o2, a_conc_methane, &
	                 a_net_methane_unsat, a_net_methane_sat, &
	                 a_methane_prod_depth_unsat, a_methane_prod_depth_sat, &
		                 a_o2_decomp_depth_unsat, a_o2_decomp_depth_sat, &
		                 a_co2_decomp_depth_unsat, a_co2_decomp_depth_sat, &
		                 a_methane_oxid_depth_unsat, a_methane_oxid_depth_sat, &
		                 a_o2_oxid_depth_unsat, a_o2_oxid_depth_sat, &
		                 a_co2_oxid_depth_unsat, a_co2_oxid_depth_sat, &
		                 a_methane_aere_depth_unsat, a_methane_aere_depth_sat, &
		                 a_methane_tran_depth_unsat, a_methane_tran_depth_sat, &
		                 a_o2_aere_depth_unsat, a_o2_aere_depth_sat, &
		                 a_co2_aere_depth_unsat, a_co2_aere_depth_sat, &
		                 a_methane_ebul_depth_unsat, a_methane_ebul_depth_sat, &
		                 a_o2stress_unsat, a_o2stress_sat, &
		                 a_methane_stress_unsat, a_methane_stress_sat, &
		                 a_methane_surf_flux_tot_unsat, a_methane_surf_flux_tot_sat, &
		                 a_methane_surf_aere_unsat, a_methane_surf_aere_sat, &
		                 a_methane_surf_ebul_unsat, a_methane_surf_ebul_sat, &
		                 a_methane_surf_diff_unsat, a_methane_surf_diff_sat, &
		                 a_methane_ebul_tot_unsat, a_methane_ebul_tot_sat, &
		                 a_methane_prod_tot_unsat, a_methane_prod_tot_sat, &
		                 a_methane_oxid_tot_unsat, a_methane_oxid_tot_sat, &
		                 a_co2_decomp_tot_unsat, a_co2_decomp_tot_sat, &
		                 a_co2_oxid_tot_unsat, a_co2_oxid_tot_sat, &
		                 a_co2_net_tot_unsat, a_co2_net_tot_sat, &
		                 a_totcol_methane_unsat, a_totcol_methane_sat, &
		                 a_grnd_methane_cond_unsat, a_grnd_methane_cond_sat, &
	                 a_conc_o2_unsat, a_conc_o2_sat, &
	                 a_conc_methane_unsat, a_conc_methane_sat, &
                 a_methane_prod_depth_lake, a_methane_oxid_depth_lake, a_methane_ebul_depth_lake, &
                 a_co2_decomp_depth_lake, a_co2_oxid_depth_lake, &
                 a_methane_surf_ebul_lake, a_methane_surf_diff_lake, a_methane_surf_flux_tot_lake, &
                 a_methane_prod_tot_lake, a_methane_oxid_tot_lake, a_methane_ebul_tot_lake, &
                 a_co2_decomp_tot_lake, a_co2_oxid_tot_lake, a_co2_net_tot_lake, &
                 a_totcol_methane_lake, a_grnd_methane_cond_lake, &
                 a_conc_o2_lake, a_conc_methane_lake, &
	                 a_forc_pmethanem, a_layer_sat_lag, &
	                 a_annavg_agnpp, a_annavg_bgnpp, a_annavg_somhr, a_annavg_finrw, &
	                 a_methane_dfsat_tot, a_f_h2osfc, &
	                 a_methane_finundated, a_methane_soil_finundated, a_methane_soil_zwt, &
	                 a_f_inund_flood_patch, a_f_inund_flood_depth_patch, a_wetland_frac_per_patch, &
	                 a_methane_surf_flux_wetland, a_methane_surf_flux_soil, &
	                 a_methane_surf_flux_lake, a_methane_surf_flux_rice, &
	                 a_B_methanogen, a_B_methanotroph, &
                 a_B_methanogen_dormant, a_B_methanotroph_dormant, &
                 a_f_T_methanogen, a_f_S_methanogen, a_f_O2_methanogen, a_f_T_methanotroph, &
                 a_methanogen_growth_rate, a_methanotroph_growth_rate, &
                 a_microbial_prod_potential, a_microbial_oxid_potential, &
                 a_methane_acc_num, a_methane_acc_num_unsat, a_methane_acc_num_sat, &
                 a_methane_acc_num_lake, a_methane_acc_num_extra, a_methane_acc_num_microbe
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer, intent(in) :: itime_in_file
   type(block_data_real8_2d), intent(inout) :: sumarea
   logical, intent(inout) :: filter(:)
   integer, intent(in) :: nl_soil
   logical, intent(in) :: forcing_has_missing_value
   logical, intent(in) :: forcmask_pch(:)

            IF (igas_ch4 > 0) THEN
               ! Outer filter (patchtype==0) excludes wetland (patchtype==2) and
               ! lake (patchtype==4) — the only patches where the CH4 module
               ! actually runs.  Without overriding, every CH4 grid cell aggregates
               ! to zero because only soil patches are summed.  Re-derive a
               ! soil+wetland filter for the main/sat/unsat writes here, then
               ! swap to a lake-only filter before the lake section below.
	               ! Active CH4 mask mirrors the driver gate at
	               ! CoLMDRIVER.F90:381-384 — patchtype==2 and patchtype==0 are
	               ! active by default.  Urban (1) and ice (3) are NEVER active
	               ! even though "<=2" would include urban.
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  ! Use BgcLink helper so the history filter matches the
                  ! driver gate + rice paddy override exactly.  Same call
                  ! site pattern as MOD_Tracer_Methane_AccFlux.F90.
                  CALL methane_patch_active_mask (numpatch, DEF_METHANE_only_wetland, &
                                                  DEF_METHANE_enable_rice_paddy, &
                                                  patchtype, filter)
                  filter = filter .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
               IF (HistForm == 'Gridded') THEN
                  CALL mp2g_hist%get_sumarea (sumarea, filter)
               ENDIF

               CALL write_history_variable_2d (mhist_on('f_net_methane'), a_net_methane, file_hist, &
                  'f_net_methane', itime_in_file, sumarea, filter, &
                  'average net methane correction to CO2 flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth'), a_methane_prod_depth, file_hist, &
                  'f_methane_prod_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth'), a_o2_decomp_depth, file_hist, &
                  'f_o2_decomp_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 consumption during decomposition per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth'), a_co2_decomp_depth, file_hist, &
                  'f_co2_decomp_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'diagnostic CO2 decomp/methanogenesis before O2-stress scaling per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth'), a_methane_oxid_depth, file_hist, &
                  'f_methane_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth'), a_o2_oxid_depth, file_hist, &
                  'f_o2_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 consumption via oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth'), a_co2_oxid_depth, file_hist, &
                  'f_co2_oxid_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CO2 production from CH4 oxidation per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_aere_depth'), a_methane_aere_depth, file_hist, &
                  'f_methane_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via aerenchyma per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_tran_depth'), a_methane_tran_depth, file_hist, &
                  'f_methane_tran_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via transpiration per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2_aere_depth'), a_o2_aere_depth, file_hist, &
                  'f_o2_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 gain via aerenchyma per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_co2_aere_depth'), a_co2_aere_depth, file_hist, &
                  'f_co2_aere_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'diagnosed CO2 aerenchyma flux per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth'), a_methane_ebul_depth, file_hist, &
                  'f_methane_ebul_depth', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 loss via ebullition per layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_o2stress'), a_o2stress, file_hist, &
                  'f_o2stress', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 availability/demand ratio', '-', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_methane_stress'), a_methane_stress, file_hist, &
                  'f_methane_stress', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 availability/sinks ratio', '-', &
                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot'), a_methane_surf_flux_tot, file_hist, &
	                  'f_methane_surf_flux_tot', itime_in_file, sumarea, filter, &
	                  'CH4 total flux to atmosphere including diffusion, ebullition, aerenchyma, and transpiration', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_aere'), a_methane_surf_aere, file_hist, &
                  'f_methane_surf_aere', itime_in_file, sumarea, filter, &
                  'CH4 column aerenchyma flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul'), a_methane_surf_ebul, file_hist, &
                  'f_methane_surf_ebul', itime_in_file, sumarea, filter, &
                  'CH4 surface ebullition flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
		               CALL write_history_variable_2d (mhist_on('f_methane_surf_diff'), a_methane_surf_diff, file_hist, &
		                  'f_methane_surf_diff', itime_in_file, sumarea, filter, &
		                  'CH4 diffusive surface flux plus inundated-fraction and closure adjustments', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num)
		               CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_phys'), a_methane_surf_diff_phys, file_hist, &
		                  'f_methane_surf_diff_phys', itime_in_file, sumarea, filter, &
		                  'CH4 pure-physics diffusive surface flux (before negative-conc clip and closure residual)', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_methane_balance_residual'), a_methane_balance_residual, file_hist, &
	                  'f_methane_balance_residual', itime_in_file, sumarea, filter, &
	                  'CH4 numerical closure flux credited to f_methane_surf_diff', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_o2_cap_loss'), a_o2_cap_loss, file_hist, &
	                  'f_o2_cap_loss', itime_in_file, sumarea, filter, &
	                  'O2 column sink from post-solve physical concentration cap', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_o2_cap_gain'), a_o2_cap_gain, file_hist, &
	                  'f_o2_cap_gain', itime_in_file, sumarea, filter, &
	                  'O2 column source from post-solve nonnegative concentration floor', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot'), a_methane_ebul_tot, file_hist, &
                  'f_methane_ebul_tot', itime_in_file, sumarea, filter, &
                  'CH4 column ebullition', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_prod_tot'), a_methane_prod_tot, file_hist, &
                  'f_methane_prod_tot', itime_in_file, sumarea, filter, &
                  'CH4 column production', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot'), a_methane_oxid_tot, file_hist, &
                  'f_methane_oxid_tot', itime_in_file, sumarea, filter, &
                  'CH4 column oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot'), a_co2_decomp_tot, file_hist, &
                  'f_co2_decomp_tot', itime_in_file, sumarea, filter, &
	                  'diagnostic CO2 decomp/methanogenesis before O2-stress scaling', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot'), a_co2_oxid_tot, file_hist, &
                  'f_co2_oxid_tot', itime_in_file, sumarea, filter, &
                  'CO2 column source from CH4 oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_aere_tot'), a_co2_aere_tot, file_hist, &
                  'f_co2_aere_tot', itime_in_file, sumarea, filter, &
                  'diagnosed CO2 column aerenchyma flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot'), a_co2_net_tot, file_hist, &
                  'f_co2_net_tot', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source from Methane module', 'mol/m2/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_totcol_methane'), a_totcol_methane, file_hist, &
                  'f_totcol_methane', itime_in_file, sumarea, filter, &
                  'total CH4 in soil column', 'mol/m2', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond'), a_grnd_methane_cond, file_hist, &
                  'f_grnd_methane_cond', itime_in_file, sumarea, filter, &
                  'boundary-layer CH4 conductance', 'm/s', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_conc_o2'), a_conc_o2, file_hist, &
                  'f_conc_o2', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'O2 concentration per soil layer', 'mol/m3', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_3d (mhist_on('f_conc_methane'), a_conc_methane, file_hist, &
                  'f_conc_methane', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 concentration per soil layer', 'mol/m3', &
                  acc_num=a_methane_acc_num)

	               ! optional microbial-pool diagnostics
	               IF (allocated(a_B_methanogen)) THEN
	                  CALL write_history_variable_3d (mhist_on('f_methane_B_methanogen'), a_B_methanogen, file_hist, &
	                     'f_methane_B_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'active methanogen biomass', 'gC/m3', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (mhist_on('f_methane_B_methanotroph'), a_B_methanotroph, file_hist, &
	                     'f_methane_B_methanotroph', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'active methanotroph biomass', 'gC/m3', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (mhist_on('f_methane_B_methanogen_dormant'), a_B_methanogen_dormant, file_hist, &
	                     'f_methane_B_methanogen_dormant', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'dormant methanogen biomass; active only when use_microbial_dormancy is true', 'gC/m3', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (&
	                     mhist_on('f_methane_B_methanotroph_dormant'), a_B_methanotroph_dormant, file_hist, &
	                     'f_methane_B_methanotroph_dormant', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'dormant methanotroph biomass; active only when use_microbial_dormancy is true', 'gC/m3', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_T_methanogen'), a_f_T_methanogen, file_hist, &
	                     'f_methane_microbe_f_T_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'methanogen microbial temperature factor', '-', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_S_methanogen'), a_f_S_methanogen, file_hist, &
	                     'f_methane_microbe_f_S_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'methanogen substrate factor', '-', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_O2_methanogen'), a_f_O2_methanogen, file_hist, &
	                     'f_methane_microbe_f_O2_methanogen', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'methanogen oxygen inhibition factor', '-', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (mhist_on('f_methane_microbe_f_T_methanotroph'), a_f_T_methanotroph, file_hist, &
	                     'f_methane_microbe_f_T_methanotroph', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'methanotroph microbial temperature factor', '-', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (&
	                     mhist_on('f_methane_methanogen_growth_rate'), a_methanogen_growth_rate, file_hist, &
	                     'f_methane_methanogen_growth_rate', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'methanogen active-pool growth rate', 'day-1', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (&
	                     mhist_on('f_methane_methanotroph_growth_rate'), a_methanotroph_growth_rate, file_hist, &
	                     'f_methane_methanotroph_growth_rate', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'methanotroph active-pool growth rate', 'day-1', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (&
	                     mhist_on('f_methane_microbial_prod_potential'), a_microbial_prod_potential, file_hist, &
	                     'f_methane_microbial_prod_potential', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'carbon-capped microbial CH4 production potential', 'mol/m3/s', &
	                     acc_num=a_methane_acc_num_microbe)
	                  CALL write_history_variable_3d (&
	                     mhist_on('f_methane_microbial_oxid_potential'), a_microbial_oxid_potential, file_hist, &
	                     'f_methane_microbial_oxid_potential', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                     'CH4/O2-capped microbial CH4 oxidation potential', 'mol/m3/s', &
	                     acc_num=a_methane_acc_num_microbe)
	               ENDIF

               ! unsat / sat
               CALL write_history_variable_2d (mhist_on('f_net_methane_unsat'), a_net_methane_unsat, file_hist, &
                  'f_net_methane_unsat', itime_in_file, sumarea, filter, &
                  'average net methane correction (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_net_methane_sat'), a_net_methane_sat, file_hist, &
                  'f_net_methane_sat', itime_in_file, sumarea, filter, &
                  'average net methane correction (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_unsat'), a_methane_prod_depth_unsat, file_hist, &
                  'f_methane_prod_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_sat'), a_methane_prod_depth_sat, file_hist, &
                  'f_methane_prod_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 production per layer (saturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_unsat'), a_methane_oxid_depth_unsat, file_hist, &
                  'f_methane_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_sat'), a_methane_oxid_depth_sat, file_hist, &
	                  'f_methane_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 oxidation per layer (saturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth_unsat'), a_o2_decomp_depth_unsat, file_hist, &
	                  'f_o2_decomp_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'O2 consumption from decomposition/root respiration per layer (unsaturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_3d (mhist_on('f_o2_decomp_depth_sat'), a_o2_decomp_depth_sat, file_hist, &
		                  'f_o2_decomp_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'O2 consumption from decomposition/root respiration per layer (saturated)', 'mol/m3/s', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth_unsat'), a_o2_oxid_depth_unsat, file_hist, &
		                  'f_o2_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'O2 consumption from CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_3d (mhist_on('f_o2_oxid_depth_sat'), a_o2_oxid_depth_sat, file_hist, &
		                  'f_o2_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'O2 consumption from CH4 oxidation per layer (saturated)', 'mol/m3/s', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_3d (mhist_on('f_o2_aere_depth_unsat'), a_o2_aere_depth_unsat, file_hist, &
		                  'f_o2_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'O2 aerenchyma supply per layer (unsaturated)', 'mol/m3/s', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_3d (mhist_on('f_o2_aere_depth_sat'), a_o2_aere_depth_sat, file_hist, &
		                  'f_o2_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'O2 aerenchyma supply per layer (saturated)', 'mol/m3/s', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_3d (mhist_on('f_o2stress_unsat'), a_o2stress_unsat, file_hist, &
		                  'f_o2stress_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'O2 availability stress factor (unsaturated)', '-', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_3d (mhist_on('f_o2stress_sat'), a_o2stress_sat, file_hist, &
		                  'f_o2stress_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'O2 availability stress factor (saturated)', '-', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_3d (mhist_on('f_methane_stress_unsat'), a_methane_stress_unsat, file_hist, &
		                  'f_methane_stress_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'CH4 availability stress factor (unsaturated)', '-', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_3d (mhist_on('f_methane_stress_sat'), a_methane_stress_sat, file_hist, &
		                  'f_methane_stress_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'CH4 availability stress factor (saturated)', '-', &
		                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_3d (mhist_on('f_methane_aere_depth_unsat'), a_methane_aere_depth_unsat, file_hist, &
	                  'f_methane_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 aerenchyma loss per layer (unsaturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_methane_aere_depth_sat'), a_methane_aere_depth_sat, file_hist, &
	                  'f_methane_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 aerenchyma loss per layer (saturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_3d (mhist_on('f_methane_tran_depth_unsat'), a_methane_tran_depth_unsat, file_hist, &
	                  'f_methane_tran_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 transpiration loss per layer (unsaturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_methane_tran_depth_sat'), a_methane_tran_depth_sat, file_hist, &
	                  'f_methane_tran_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 transpiration loss per layer (saturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_unsat'), a_methane_ebul_depth_unsat, file_hist, &
	                  'f_methane_ebul_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 ebullition loss per layer (unsaturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_sat'), a_methane_ebul_depth_sat, file_hist, &
	                  'f_methane_ebul_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CH4 ebullition loss per layer (saturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_unsat'), a_co2_decomp_depth_unsat, file_hist, &
	                  'f_co2_decomp_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
		                  'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_sat'), a_co2_decomp_depth_sat, file_hist, &
                  'f_co2_decomp_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (saturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_unsat'), a_co2_oxid_depth_unsat, file_hist, &
                  'f_co2_oxid_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CO2 from CH4 oxidation per layer (unsaturated)', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_sat'), a_co2_oxid_depth_sat, file_hist, &
	                  'f_co2_oxid_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CO2 from CH4 oxidation per layer (saturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_3d (mhist_on('f_co2_aere_depth_unsat'), a_co2_aere_depth_unsat, file_hist, &
	                  'f_co2_aere_depth_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CO2 aerenchyma diagnostic per layer (unsaturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_co2_aere_depth_sat'), a_co2_aere_depth_sat, file_hist, &
	                  'f_co2_aere_depth_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'CO2 aerenchyma diagnostic per layer (saturated)', 'mol/m3/s', &
	                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_unsat'), a_co2_decomp_tot_unsat, file_hist, &
                  'f_co2_decomp_tot_unsat', itime_in_file, sumarea, filter, &
	                  'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_sat'), a_co2_decomp_tot_sat, file_hist, &
                  'f_co2_decomp_tot_sat', itime_in_file, sumarea, filter, &
	                  'diagnostic CO2 decomp/methanogenesis before O2-stress scaling (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_unsat'), a_co2_oxid_tot_unsat, file_hist, &
                  'f_co2_oxid_tot_unsat', itime_in_file, sumarea, filter, &
                  'CO2 oxidation column source (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_sat'), a_co2_oxid_tot_sat, file_hist, &
                  'f_co2_oxid_tot_sat', itime_in_file, sumarea, filter, &
                  'CO2 oxidation column source (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_unsat'), a_co2_net_tot_unsat, file_hist, &
                  'f_co2_net_tot_unsat', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_sat'), a_co2_net_tot_sat, file_hist, &
                  'f_co2_net_tot_sat', itime_in_file, sumarea, filter, &
                  'net diagnosed CO2 source (saturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_sat)
               CALL write_history_variable_2d (&
                  mhist_on('f_methane_surf_flux_tot_unsat'), a_methane_surf_flux_tot_unsat, file_hist, &
                  'f_methane_surf_flux_tot_unsat', itime_in_file, sumarea, filter, &
                  'CH4 surface flux (unsaturated)', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_sat'), a_methane_surf_flux_tot_sat, file_hist, &
	                  'f_methane_surf_flux_tot_sat', itime_in_file, sumarea, filter, &
	                  'CH4 surface flux (saturated)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_aere_unsat'), a_methane_surf_aere_unsat, file_hist, &
	                  'f_methane_surf_aere_unsat', itime_in_file, sumarea, filter, &
	                  'CH4 column aerenchyma flux (unsaturated)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_aere_sat'), a_methane_surf_aere_sat, file_hist, &
	                  'f_methane_surf_aere_sat', itime_in_file, sumarea, filter, &
	                  'CH4 column aerenchyma flux (saturated)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_unsat'), a_methane_surf_ebul_unsat, file_hist, &
	                  'f_methane_surf_ebul_unsat', itime_in_file, sumarea, filter, &
	                  'CH4 surface ebullition flux (unsaturated)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_sat'), a_methane_surf_ebul_sat, file_hist, &
	                  'f_methane_surf_ebul_sat', itime_in_file, sumarea, filter, &
	                  'CH4 surface ebullition flux (saturated)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_unsat'), a_methane_surf_diff_unsat, file_hist, &
	                  'f_methane_surf_diff_unsat', itime_in_file, sumarea, filter, &
	                  'CH4 diffusive surface flux (unsaturated)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_sat'), a_methane_surf_diff_sat, file_hist, &
		                  'f_methane_surf_diff_sat', itime_in_file, sumarea, filter, &
		                  'CH4 diffusive surface flux (saturated)', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_unsat'), a_methane_ebul_tot_unsat, file_hist, &
		                  'f_methane_ebul_tot_unsat', itime_in_file, sumarea, filter, &
		                  'CH4 column ebullition (unsaturated)', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_sat'), a_methane_ebul_tot_sat, file_hist, &
		                  'f_methane_ebul_tot_sat', itime_in_file, sumarea, filter, &
		                  'CH4 column ebullition (saturated)', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_unsat'), a_methane_prod_tot_unsat, file_hist, &
		                  'f_methane_prod_tot_unsat', itime_in_file, sumarea, filter, &
		                  'CH4 column production (unsaturated)', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_sat'), a_methane_prod_tot_sat, file_hist, &
		                  'f_methane_prod_tot_sat', itime_in_file, sumarea, filter, &
		                  'CH4 column production (saturated)', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_unsat'), a_methane_oxid_tot_unsat, file_hist, &
		                  'f_methane_oxid_tot_unsat', itime_in_file, sumarea, filter, &
		                  'CH4 column oxidation (unsaturated)', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_sat'), a_methane_oxid_tot_sat, file_hist, &
		                  'f_methane_oxid_tot_sat', itime_in_file, sumarea, filter, &
		                  'CH4 column oxidation (saturated)', 'mol/m2/s', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_2d (mhist_on('f_totcol_methane_unsat'), a_totcol_methane_unsat, file_hist, &
		                  'f_totcol_methane_unsat', itime_in_file, sumarea, filter, &
		                  'total CH4 in unsaturated column', 'mol/m2', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_2d (mhist_on('f_totcol_methane_sat'), a_totcol_methane_sat, file_hist, &
		                  'f_totcol_methane_sat', itime_in_file, sumarea, filter, &
		                  'total CH4 in saturated column', 'mol/m2', &
		                  acc_num=a_methane_acc_num_sat)
		               CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_unsat'), a_grnd_methane_cond_unsat, file_hist, &
		                  'f_grnd_methane_cond_unsat', itime_in_file, sumarea, filter, &
		                  'CH4 boundary-layer conductance (unsaturated)', 'm/s', &
		                  acc_num=a_methane_acc_num_unsat)
		               CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_sat'), a_grnd_methane_cond_sat, file_hist, &
		                  'f_grnd_methane_cond_sat', itime_in_file, sumarea, filter, &
		                  'CH4 boundary-layer conductance (saturated)', 'm/s', &
		                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_3d (mhist_on('f_conc_o2_unsat'), a_conc_o2_unsat, file_hist, &
	                  'f_conc_o2_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'O2 concentration per layer (unsaturated)', 'mol/m3', &
	                  acc_num=a_methane_acc_num_unsat)
	               CALL write_history_variable_3d (mhist_on('f_conc_o2_sat'), a_conc_o2_sat, file_hist, &
	                  'f_conc_o2_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
	                  'O2 concentration per layer (saturated)', 'mol/m3', &
	                  acc_num=a_methane_acc_num_sat)
	               CALL write_history_variable_3d (mhist_on('f_conc_methane_unsat'), a_conc_methane_unsat, file_hist, &
	                  'f_conc_methane_unsat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 concentration per layer (unsaturated)', 'mol/m3', &
                  acc_num=a_methane_acc_num_unsat)
               CALL write_history_variable_3d (mhist_on('f_conc_methane_sat'), a_conc_methane_sat, file_hist, &
                  'f_conc_methane_sat', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'CH4 concentration per layer (saturated)', 'mol/m3', &
                  acc_num=a_methane_acc_num_sat)


               ! lake diagnostics (CTSM lake CH4 path) — swap to lake-only filter
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  filter = (patchtype == 4) .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
               IF (HistForm == 'Gridded') THEN
                  CALL mp2g_hist%get_sumarea (sumarea, filter)
               ENDIF

               CALL write_history_variable_3d (mhist_on('f_methane_prod_depth_lake'), a_methane_prod_depth_lake, file_hist, &
                  'f_methane_prod_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 production per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_methane_oxid_depth_lake'), a_methane_oxid_depth_lake, file_hist, &
                  'f_methane_oxid_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 oxidation per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_co2_decomp_depth_lake'), a_co2_decomp_depth_lake, file_hist, &
                  'f_co2_decomp_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CO2 decomp/methanogenesis per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_co2_oxid_depth_lake'), a_co2_oxid_depth_lake, file_hist, &
                  'f_co2_oxid_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CO2 from CH4 oxidation per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_methane_ebul_depth_lake'), a_methane_ebul_depth_lake, file_hist, &
                  'f_methane_ebul_depth_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 ebullition loss per sediment layer', 'mol/m3/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_ebul_lake'), a_methane_surf_ebul_lake, file_hist, &
                  'f_methane_surf_ebul_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 surface ebullition flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_surf_diff_lake'), a_methane_surf_diff_lake, file_hist, &
                  'f_methane_surf_diff_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 surface diffusion flux', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_tot_lake'), a_methane_surf_flux_tot_lake, file_hist, &
	                  'f_methane_surf_flux_tot_lake', itime_in_file, sumarea, filter, &
	                  'lake CH4 total surface flux (ebullition + diffusion)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num_lake)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_lake'), a_methane_surf_flux_lake, file_hist, &
	                  'f_methane_surf_flux_lake', itime_in_file, sumarea, filter, &
	                  'lake contribution to CH4 total surface flux', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num_lake)
	               CALL write_history_variable_2d (mhist_on('f_methane_prod_tot_lake'), a_methane_prod_tot_lake, file_hist, &
                  'f_methane_prod_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column production', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_oxid_tot_lake'), a_methane_oxid_tot_lake, file_hist, &
                  'f_methane_oxid_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column oxidation', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_methane_ebul_tot_lake'), a_methane_ebul_tot_lake, file_hist, &
                  'f_methane_ebul_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 column ebullition', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_decomp_tot_lake'), a_co2_decomp_tot_lake, file_hist, &
                  'f_co2_decomp_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CO2 decomp/methanogenesis column source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_oxid_tot_lake'), a_co2_oxid_tot_lake, file_hist, &
                  'f_co2_oxid_tot_lake', itime_in_file, sumarea, filter, &
                  'lake CO2 oxidation column source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_co2_net_tot_lake'), a_co2_net_tot_lake, file_hist, &
                  'f_co2_net_tot_lake', itime_in_file, sumarea, filter, &
                  'lake net diagnosed CO2 source', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_totcol_methane_lake'), a_totcol_methane_lake, file_hist, &
                  'f_totcol_methane_lake', itime_in_file, sumarea, filter, &
                  'total CH4 in lake sediment column', 'mol/m2', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_2d (mhist_on('f_grnd_methane_cond_lake'), a_grnd_methane_cond_lake, file_hist, &
                  'f_grnd_methane_cond_lake', itime_in_file, sumarea, filter, &
                  'lake CH4 boundary-layer conductance', 'm/s', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_conc_o2_lake'), a_conc_o2_lake, file_hist, &
                  'f_conc_o2_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake O2 concentration per sediment layer', 'mol/m3', &
                  acc_num=a_methane_acc_num_lake)
               CALL write_history_variable_3d (mhist_on('f_conc_methane_lake'), a_conc_methane_lake, file_hist, &
                  'f_conc_methane_lake', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'lake CH4 concentration per sediment layer', 'mol/m3', &
                  acc_num=a_methane_acc_num_lake)

               ! extras + microbe diagnostics — swap back to active CH4 mask
               IF ((p_is_worker) .and. (numpatch > 0)) THEN
                  ! Use BgcLink helper so the history filter matches the
                  ! driver gate + rice paddy override exactly.  Same call
                  ! site pattern as MOD_Tracer_Methane_AccFlux.F90.
                  CALL methane_patch_active_mask (numpatch, DEF_METHANE_only_wetland, &
                                                  DEF_METHANE_enable_rice_paddy, &
                                                  patchtype, filter)
                  filter = filter .and. patchmask
                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
               ENDIF
               IF (HistForm == 'Gridded') THEN
                  CALL mp2g_hist%get_sumarea (sumarea, filter)
               ENDIF

               CALL write_history_variable_2d (mhist_on('f_forc_pmethanem'), a_forc_pmethanem, file_hist, &
                  'f_forc_pmethanem', itime_in_file, sumarea, filter, &
                  'atmospheric CH4 partial pressure forcing', 'Pa', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_3d (mhist_on('f_layer_sat_lag'), a_layer_sat_lag, file_hist, &
                  'f_layer_sat_lag', itime_in_file, 'soil', 1, nl_soil, sumarea, filter, &
                  'time-lagged saturation flag per layer', '-', &
                  acc_num=a_methane_acc_num)
               CALL write_history_variable_2d (mhist_on('f_annavg_agnpp'), a_annavg_agnpp, file_hist, &
                  'f_annavg_agnpp', itime_in_file, sumarea, filter, &
                  'annual avg above-ground NPP', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_bgnpp'), a_annavg_bgnpp, file_hist, &
                  'f_annavg_bgnpp', itime_in_file, sumarea, filter, &
                  'annual avg below-ground NPP', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_somhr'), a_annavg_somhr, file_hist, &
                  'f_annavg_somhr', itime_in_file, sumarea, filter, &
                  'annual avg SOM heterotrophic respiration', 'gC/m2/s', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_annavg_finrw'), a_annavg_finrw, file_hist, &
                  'f_annavg_finrw', itime_in_file, sumarea, filter, &
                  'respiration-weighted annual avg inundated fraction', '-', &
                  acc_num=a_methane_acc_num_extra)
               CALL write_history_variable_2d (mhist_on('f_methane_dfsat_tot'), a_methane_dfsat_tot, file_hist, &
                  'f_methane_dfsat_tot', itime_in_file, sumarea, filter, &
                  'CH4 flux from decreasing inundated fraction', 'mol/m2/s', &
                  acc_num=a_methane_acc_num_extra)
	               CALL write_history_variable_2d (mhist_on('f_f_h2osfc'), a_f_h2osfc, file_hist, &
	                  'f_f_h2osfc', itime_in_file, sumarea, filter, &
	                  'fraction of surface water (CLM5 microtopography)', '-', &
	                  acc_num=a_methane_acc_num_extra)
	               CALL write_history_variable_2d (mhist_on('f_methane_finundated'), a_methane_finundated, file_hist, &
	                  'f_methane_finundated', itime_in_file, sumarea, filter, &
	                  'actual inundated fraction used by CH4 physics', '-', &
	                  acc_num=a_methane_acc_num)
	               ! f_methane_soil_finundated moved to soil-only filter block below (avoid wetland zero dilution)
	               CALL write_history_variable_2d (mhist_on('f_inund_flood_patch'), a_f_inund_flood_patch, file_hist, &
	                  'f_inund_flood_patch', itime_in_file, sumarea, filter, &
	                  'routing-published flood inundation fraction on methane-active patches', '-', &
	                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_inund_flood_depth_patch'), a_f_inund_flood_depth_patch, file_hist, &
	                  'f_inund_flood_depth_patch', itime_in_file, sumarea, filter, &
	                  'routing-published floodplain water depth on methane-active patches', 'm', &
	                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_wetland_frac_patch'), a_wetland_frac_per_patch, file_hist, &
	                  'f_wetland_frac_patch', itime_in_file, sumarea, filter, &
	                  'static wetland fraction of parent cell for methane patch', '-', &
	                  acc_num=a_methane_acc_num)
	               ! Sub-tile CH4 flux contributions: each variable is the active-CH4 area mean
	               ! of that tile's contribution; non-source patches contribute 0.
	               ! Sum (wetland+soil+rice+lake) equals f_methane_surf_flux_tot exactly.
	               ! NOT the intrinsic per-tile flux (divide by tile area fraction for that).
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_wetland'), a_methane_surf_flux_wetland, file_hist, &
	                  'f_methane_surf_flux_wetland', itime_in_file, sumarea, filter, &
	                  'wetland-patch contribution to total CH4 surface flux (active-mean)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_soil'), a_methane_surf_flux_soil, file_hist, &
	                  'f_methane_surf_flux_soil', itime_in_file, sumarea, filter, &
	                  'non-rice soil-patch contribution to total CH4 surface flux (active-mean)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_methane_surf_flux_rice'), a_methane_surf_flux_rice, file_hist, &
	                  'f_methane_surf_flux_rice', itime_in_file, sumarea, filter, &
	                  'rice-paddy contribution to total CH4 surface flux (active-mean)', 'mol/m2/s', &
	                  acc_num=a_methane_acc_num)

	               IF ((p_is_worker) .and. (numpatch > 0)) THEN
	                  filter = (patchtype == 0) .and. patchmask
	                  IF (forcing_has_missing_value) filter = filter .and. forcmask_pch
	               ENDIF
	               IF (HistForm == 'Gridded') THEN
	                  CALL mp2g_hist%get_sumarea (sumarea, filter)
	               ENDIF
	               CALL write_history_variable_2d (mhist_on('f_methane_soil_zwt'), a_methane_soil_zwt, file_hist, &
	                  'f_methane_soil_zwt', itime_in_file, sumarea, filter, &
	                  'water-table depth on methane-active soil/rice patches', 'm', &
	                  acc_num=a_methane_acc_num)
	               CALL write_history_variable_2d (mhist_on('f_methane_soil_finundated'), a_methane_soil_finundated, file_hist, &
	                  'f_methane_soil_finundated', itime_in_file, sumarea, filter, &
	                  'CH4 inundated fraction on methane-active soil/rice patches (soil-only mean, no wetland-zero dilution)', '-', &
	                  acc_num=a_methane_acc_num)
	            ENDIF

   END SUBROUTINE tracer_methane_hist_out
#endif

#ifdef TRACER
   SUBROUTINE write_history_tracer_ratio_2d ( is_hist, &
         tracer_mass_acc, water_acc, file_hist, varname, itime_in_file, &
         filter, longname, units)

   USE MOD_Block
   USE MOD_SPMD_Task, only: p_is_worker, p_is_io
#ifdef TRACER
#endif
   USE MOD_Namelist, only: DEF_HIST_CompressLevel

   IMPLICIT NONE

   logical, intent(in) :: is_hist
   real(r8), intent(in) :: tracer_mass_acc(:)
   real(r8), intent(in) :: water_acc(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   logical, intent(in) :: filter(:)
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   real(r8), allocatable :: tracer_map_vec(:)
   real(r8), allocatable :: water_map_vec(:)
   type(block_data_real8_2d) :: tracer_xy_2d
   type(block_data_real8_2d) :: water_xy_2d
   type(block_data_real8_2d) :: ratio_xy_2d
   integer :: ip, iblkme, xblk, yblk, xloc, yloc
   integer :: compress

      IF (.not. is_hist) RETURN
      IF (HistForm /= 'Gridded') RETURN

      allocate(tracer_map_vec(size(tracer_mass_acc)))
      allocate(water_map_vec(size(water_acc)))
      tracer_map_vec = 0._r8
      water_map_vec = 0._r8

      IF (p_is_worker) THEN
         DO ip = 1, size(water_acc)
            IF (abs(water_acc(ip)) > trc_tiny .and. tracer_mass_acc(ip) /= spval) THEN
               tracer_map_vec(ip) = tracer_mass_acc(ip)
               water_map_vec(ip) = water_acc(ip)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, tracer_xy_2d)
         CALL allocate_block_data (ghist, water_xy_2d)
         CALL allocate_block_data (ghist, ratio_xy_2d)
      ENDIF

      CALL mp2g_hist%pset2grid (tracer_map_vec, tracer_xy_2d, spv = spval, msk = filter)
      CALL mp2g_hist%pset2grid (water_map_vec, water_xy_2d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)
                  IF (water_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval .and. &
                      tracer_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval .and. &
                      abs(water_xy_2d%blk(xblk,yblk)%val(xloc,yloc)) > trc_tiny) THEN
                     ratio_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = &
                        tracer_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                        / water_xy_2d%blk(xblk,yblk)%val(xloc,yloc)
                  ELSE
                     ratio_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, &
         ratio_xy_2d, compress, longname, units)

      deallocate(tracer_map_vec)
      deallocate(water_map_vec)

   END SUBROUTINE write_history_tracer_ratio_2d


   SUBROUTINE write_history_tracer_ratio_3d ( is_hist, &
         tracer_mass_acc, water_acc, file_hist, varname, itime_in_file, &
         dim1name, lb1, ndim1, filter, longname, units)

   USE MOD_Block
   USE MOD_SPMD_Task, only: p_is_worker, p_is_io
#ifdef TRACER
#endif
   USE MOD_Namelist, only: DEF_HIST_CompressLevel

   IMPLICIT NONE

   logical, intent(in) :: is_hist
   real(r8), intent(in) :: tracer_mass_acc(:,:)
   real(r8), intent(in) :: water_acc(:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name
   integer, intent(in) :: lb1, ndim1
   logical, intent(in) :: filter(:)
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   real(r8), allocatable :: tracer_map_vec(:,:)
   real(r8), allocatable :: water_map_vec(:,:)
   type(block_data_real8_3d) :: tracer_xy_3d
   type(block_data_real8_3d) :: water_xy_3d
   type(block_data_real8_3d) :: ratio_xy_3d
   integer :: ip, i1, iblkme, xblk, yblk, xloc, yloc
   integer :: lb_vec, ub_vec
   integer :: compress

      IF (.not. is_hist) RETURN
      IF (HistForm /= 'Gridded') RETURN

      lb_vec = lbound(tracer_mass_acc, 1)
      ub_vec = ubound(tracer_mass_acc, 1)
      allocate(tracer_map_vec(lb_vec:ub_vec, size(tracer_mass_acc, 2)))
      allocate(water_map_vec  (lb_vec:ub_vec, size(water_acc, 2)))
      tracer_map_vec = spval
      water_map_vec  = spval

      IF (p_is_worker) THEN
         DO ip = 1, size(water_acc, 2)
            DO i1 = lb_vec, ub_vec
               IF (water_acc(i1, ip) > trc_tiny .and. tracer_mass_acc(i1, ip) /= spval) THEN
                  tracer_map_vec(i1, ip) = tracer_mass_acc(i1, ip)
                  water_map_vec  (i1, ip) = water_acc(i1, ip)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, tracer_xy_3d, ndim1, lb1)
         CALL allocate_block_data (ghist, water_xy_3d,   ndim1, lb1)
         CALL allocate_block_data (ghist, ratio_xy_3d,   ndim1, lb1)
      ENDIF

      CALL mp2g_hist%pset2grid (tracer_map_vec, tracer_xy_3d, spv = spval, msk = filter)
      CALL mp2g_hist%pset2grid (water_map_vec,   water_xy_3d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)
                  DO i1 = ratio_xy_3d%lb1, ratio_xy_3d%ub1
                     IF (water_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval .and. &
                         tracer_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) /= spval .and. &
                         water_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) > trc_tiny) THEN
                        ratio_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) = &
                           tracer_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) &
                           / water_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc)
                     ELSE
                        ratio_xy_3d%blk(xblk,yblk)%val(i1,xloc,yloc) = spval
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_3d (file_hist, varname, dim1name, ghist, &
         itime_in_file, ratio_xy_3d, compress, longname, units)

      deallocate(tracer_map_vec)
      deallocate(water_map_vec)

   END SUBROUTINE write_history_tracer_ratio_3d


   SUBROUTINE write_history_tracer_delta_2d ( is_hist, &
         tracer_mass_acc, water_acc, ref_ratio, file_hist, varname, itime_in_file, &
         filter, longname, units, water_min_override)

   USE MOD_Block
   USE MOD_SPMD_Task, only: p_is_worker, p_is_io
#ifdef TRACER
#endif
   USE MOD_Namelist, only: DEF_HIST_CompressLevel

   IMPLICIT NONE

   logical, intent(in) :: is_hist
   real(r8), intent(in) :: tracer_mass_acc(:)
   real(r8), intent(in) :: water_acc(:)
   real(r8), intent(in) :: ref_ratio
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
	   logical, intent(in) :: filter(:)
	   character(len=*), intent(in) :: longname
	   character(len=*), intent(in) :: units
	   real(r8), intent(in), optional :: water_min_override

	   real(r8), allocatable :: tracer_map_vec(:)
	   real(r8), allocatable :: water_map_vec(:)
	   type(block_data_real8_2d) :: tracer_xy_2d
	   type(block_data_real8_2d) :: water_xy_2d
	   type(block_data_real8_2d) :: delta_xy_2d
	   integer :: ip, iblkme, xblk, yblk, xloc, yloc
	   integer :: compress
	   real(r8) :: delta_loc, water_min

      IF (.not. is_hist) RETURN
      IF (HistForm /= 'Gridded') RETURN

      allocate(tracer_map_vec(size(tracer_mass_acc)))
	      allocate(water_map_vec(size(water_acc)))
	      tracer_map_vec = 0._r8
	      water_map_vec = 0._r8
	      water_min = trc_water_min_for_delta
	      IF (present(water_min_override)) water_min = water_min_override

	      IF (p_is_worker) THEN
	         DO ip = 1, size(water_acc)
	            IF (water_acc(ip) > water_min .and. &
	                tracer_mass_acc(ip) /= spval) THEN
               tracer_map_vec(ip) = tracer_mass_acc(ip)
               water_map_vec(ip) = water_acc(ip)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, tracer_xy_2d)
         CALL allocate_block_data (ghist, water_xy_2d)
         CALL allocate_block_data (ghist, delta_xy_2d)
      ENDIF

      CALL mp2g_hist%pset2grid (tracer_map_vec, tracer_xy_2d, spv = spval, msk = filter)
      CALL mp2g_hist%pset2grid (water_map_vec, water_xy_2d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, ghist%ycnt(yblk)
               DO xloc = 1, ghist%xcnt(xblk)
                  IF (water_xy_2d%blk(xblk,yblk)%val(xloc,yloc) > trc_tiny) THEN
                     delta_loc = mass_to_delta( &
                        tracer_xy_2d%blk(xblk,yblk)%val(xloc,yloc), &
                        water_xy_2d%blk(xblk,yblk)%val(xloc,yloc), &
                        ref_ratio)
                     IF (delta_loc /= spval .and. abs(delta_loc) <= trc_delta_sanity_max) THEN
                        delta_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = delta_loc
                     ELSE
                        delta_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                     ENDIF
                  ELSE
                     delta_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel
      CALL hist_write_var_real8_2d (file_hist, varname, ghist, itime_in_file, &
         delta_xy_2d, compress, longname, units)

      deallocate(tracer_map_vec)
      deallocate(water_map_vec)

   END SUBROUTINE write_history_tracer_delta_2d
#endif

   SUBROUTINE write_history_variable_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units, acc_num, input_mode)

   USE MOD_Vars_1DAccFluxes, only: nac

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8),         intent(inout) :: acc_vec(:)
   character(len=*), intent(in)    :: file_hist
   character(len=*), intent(in)    :: varname
   integer,          intent(in)    :: itime_in_file
   character(len=*), intent(in)    :: longname
   character(len=*), intent(in)    :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   real(r8), intent(in), optional  :: acc_num(:)

   character(len=*), intent(in), optional :: input_mode

      IF (.not. is_hist) RETURN

#ifndef SinglePoint
      IF ( .not. present(acc_num) ) THEN
         IF (p_is_worker) &
            WHERE (acc_vec /= spval) acc_vec = acc_vec / nac
      ELSE
         IF (p_is_worker)  &
            WHERE (acc_vec/=spval .and. acc_num>0) acc_vec = acc_vec / acc_num
      ENDIF
#else
      IF ( .not. present(acc_num) ) THEN
         WHERE (acc_vec /= spval)  acc_vec = acc_vec / nac
      ELSE
         WHERE (acc_vec/=spval .and. acc_num>0) acc_vec = acc_vec / acc_num
      ENDIF
#endif

      select CASE (HistForm)
      CASE ('Gridded')
         IF (present(input_mode)) THEN
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units, input_mode)
         ELSE
            CALL flux_map_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
         ENDIF
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         IF (present(input_mode)) THEN
            CALL aggregate_to_vector_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, filter, longname, units, input_mode)
         ELSE
            CALL aggregate_to_vector_and_write_2d ( &
               acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
         ENDIF
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_2d


#ifdef URBAN_MODEL
   SUBROUTINE write_history_variable_urb_2d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, sumarea, filter, &
         longname, units)

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

      IF (.not. is_hist) RETURN

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_urb_2d ( &
            acc_vec, file_hist, varname, itime_in_file, sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         !TODO: currently, it is not applicable to urban variables
         CALL aggregate_to_vector_and_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_2d ( &
            acc_vec, file_hist, varname, itime_in_file, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_urb_2d
#endif


   SUBROUTINE write_history_variable_3d ( is_hist, &
         acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
         sumarea, filter, longname, units, acc_num)

   USE MOD_Vars_1DAccFluxes, only: nac

   IMPLICIT NONE

   logical, intent(in) :: is_hist

   real(r8), intent(inout) :: acc_vec(:,:)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in) :: itime_in_file
   character(len=*), intent(in) :: dim1name
   integer, intent(in) :: lb1, ndim1

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)
   character (len=*), intent(in) :: longname
   character (len=*), intent(in) :: units
   real(r8), intent(in), optional :: acc_num(:)

   ! Local variables
   integer :: iblkme, xblk, yblk, xloc, yloc, i1
   integer :: compress
   integer :: j

      IF (.not. is_hist) RETURN

      ! The lower-level 3D writers (flux_map_and_write_3d, aggregate_..., single_...)
      ! always divide acc_vec by `nac`. When a per-patch acc_num counter is supplied
      ! (e.g. methane accumulators that only tick on wetland patches), pre-scale by
      ! nac/acc_num so the lower-level division yields acc_vec/acc_num.
      IF (present(acc_num)) THEN
#ifndef SinglePoint
         IF (p_is_worker) THEN
            DO j = 1, size(acc_vec, 1)
               WHERE (acc_vec(j,:) /= spval .and. acc_num > 0._r8) &
                  acc_vec(j,:) = acc_vec(j,:) * real(nac, r8) / acc_num
            ENDDO
         ENDIF
#else
         DO j = 1, size(acc_vec, 1)
            WHERE (acc_vec(j,:) /= spval .and. acc_num > 0._r8) &
               acc_vec(j,:) = acc_vec(j,:) * real(nac, r8) / acc_num
         ENDDO
#endif
      ENDIF

      select CASE (HistForm)
      CASE ('Gridded')
         CALL flux_map_and_write_3d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            sumarea, filter, longname, units)
#if (defined UNSTRUCTURED || defined CATCHMENT)
      CASE ('Vector')
         CALL aggregate_to_vector_and_write_3d ( &
            acc_vec, file_hist, varname, itime_in_file, dim1name, lb1, ndim1, &
            filter, longname, units)
#endif
#ifdef SinglePoint
      CASE ('Single')
         CALL single_write_3d (acc_vec, file_hist, varname, itime_in_file, &
            dim1name, ndim1, longname, units)
#endif
      END select

   END SUBROUTINE write_history_variable_3d

END MODULE MOD_Tracer_Hist
#endif
