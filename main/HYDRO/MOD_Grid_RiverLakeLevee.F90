#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeLevee
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Levee (flood protection) module for grid-based river-lake routing.
!   Ported from CaMa-Flood v4 levee scheme (Zhao et al., 2025, WRR).
!
!   Splits floodplain storage into protected (behind levee) and unprotected
!   (river-side) zones. Water only enters the protected zone when the levee
!   is overtopped.
!
! Created by CoLM team, March 2026
!-------------------------------------------------------------------------------------

   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_nan
   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE
   PRIVATE

   ! ----- Levee parameters (derived in levee_init) -----
   logical,  allocatable :: has_levee     (:)   ! true where levee exists
   real(r8), allocatable :: levee_frc     (:)   ! unprotected fraction [0-1]
   real(r8), allocatable :: levee_hgt     (:)   ! levee crest height above riverbed [m]
   real(r8), allocatable :: levee_dst     (:)   ! distance from river center to levee [m]
   real(r8), allocatable :: levee_bashgt  (:)   ! ground elevation at levee position [m]
   real(r8), allocatable :: levee_bassto  (:)   ! storage at levee base [m^3]
   real(r8), allocatable :: levee_topsto  (:)   ! storage at levee crest, river-side [m^3]
   real(r8), allocatable :: levee_filsto  (:)   ! storage when both sides filled to crest [m^3]
   real(r8), allocatable :: levee_fldstomax(:,:) ! total storage at each floodplain layer top [m^3]
   real(r8), allocatable :: levee_fldgrd  (:,:) ! floodplain layer slope [m/m]

   ! ----- State variables -----
   real(r8), allocatable :: levsto        (:)   ! protected-side water storage [m^3]

   ! ----- Diagnostic -----
   real(r8), allocatable :: levdph        (:)   ! protected-side water depth [m]

   PUBLIC :: has_levee, levsto, levdph
   PUBLIC :: levee_init
   PUBLIC :: read_levee_restart
   PUBLIC :: levee_fldstg
   PUBLIC :: levee_apply_protected_flux
   PUBLIC :: levee_repartition_storage
   PUBLIC :: levee_visible_volume_from_stage
   PUBLIC :: levee_final

CONTAINS

   ! =========================================================================
   SUBROUTINE allocate_levee_arrays (ncell, max_nlfp)
   ! =========================================================================

   IMPLICIT NONE

   integer, intent(in) :: ncell
   integer, intent(in) :: max_nlfp

      allocate (has_levee      (ncell))
      allocate (levee_frc      (ncell))
      allocate (levee_hgt      (ncell))
      allocate (levee_dst      (ncell))
      allocate (levee_bashgt   (ncell))
      allocate (levee_bassto   (ncell))
      allocate (levee_topsto   (ncell))
      allocate (levee_filsto   (ncell))
      allocate (levee_fldstomax(max_nlfp, ncell))
      allocate (levee_fldgrd   (max_nlfp, ncell))
      allocate (levsto         (ncell))
      allocate (levdph         (ncell))

      has_levee       = .false.
      levee_frc       = 1._r8
      levee_hgt       = 0._r8
      levee_dst       = 0._r8
      levee_bashgt    = 0._r8
      levee_bassto    = 0._r8
      levee_topsto    = 0._r8
      levee_filsto    = 0._r8
      levee_fldstomax = 0._r8
      levee_fldgrd    = 0._r8
      levsto          = 0._r8
      levdph          = 0._r8

   END SUBROUTINE allocate_levee_arrays

   ! =========================================================================
   SUBROUTINE levee_init ()
   ! =========================================================================

   USE MOD_Namelist,              only: DEF_USE_LEVEE
   USE MOD_Grid_RiverLakeNetwork, only: numucat, topo_rivlen, topo_rivwth, &
      topo_area, topo_rivstomax, floodplain_curve, levee_frc_data, levee_hgt_data, lake_type
   IMPLICIT NONE

   ! Local variables
   integer  :: i, j, ilev, nlfp, max_nlfp
   real(r8) :: rivlen, rivwth, catarea, dwth_inc
   real(r8) :: dhgtpre, dwth_fil, dwth_add, fldgrd, dhgtdif
   real(r8) :: dhgtnow, dsto_fil, dsto_add
   integer  :: n_levee_resv_clash

      IF (allocated(has_levee)) CALL levee_final()

      IF (.not. p_is_worker) THEN
         CALL allocate_levee_arrays (0, 0)
         RETURN
      ENDIF

      IF (numucat <= 0) THEN
         CALL allocate_levee_arrays (0, 0)
         RETURN
      ENDIF

      max_nlfp = 0
      DO i = 1, numucat
         max_nlfp = MAX(max_nlfp, floodplain_curve(i)%nlfp)
      ENDDO
      CALL allocate_levee_arrays (numucat, max_nlfp)

      ! When LEVEE is disabled, leave arrays at the inert default state set by
      ! allocate_levee_arrays (has_levee=.false. everywhere). Downstream guards
      ! `IF (DEF_USE_LEVEE .and. has_levee(i) ...)` then short-circuit logically
      ! AND survive ifx -check bounds, since has_levee(i) is in-range and false.
      IF (.not. DEF_USE_LEVEE) RETURN

      ! --- Copy from Network data ---
      levee_frc = levee_frc_data
      levee_hgt = levee_hgt_data

      ! --- Validate ---
      DO i = 1, numucat
         IF (ieee_is_nan(levee_hgt(i)) .or. ieee_is_nan(levee_frc(i)) .or. &
             levee_hgt(i) <= 0._r8 .or. levee_frc(i) >= 1._r8 .or. levee_frc(i) < 0._r8) THEN
            levee_hgt(i) = 0._r8
            levee_frc(i) = 1._r8
         ENDIF
         levee_frc(i) = max(0._r8, min(1._r8, levee_frc(i)))
      ENDDO

      has_levee = (levee_frc < 1.0_r8)

      ! Bug L guard: a reservoir cell (lake_type == 2) and a leveed cell are
      ! mutually exclusive. Reservoirs use volresv as the visible-side storage
      ! and never call levee_fldstg, so any non-zero levsto on a reservoir
      ! cell would be a phantom pool that the conservation budget cannot see.
      ! If the input levee_frc dataset disagrees with the reservoir mask,
      ! force has_levee=.false. on the offending cells with a one-shot info
      ! log so the data inconsistency is visible without crashing the run.
      IF (allocated(lake_type)) THEN
         n_levee_resv_clash = 0
         DO i = 1, numucat
            IF (i > size(lake_type)) EXIT
            IF (lake_type(i) == 2 .and. has_levee(i)) THEN
               has_levee(i) = .false.
               levee_frc(i) = 1._r8
               levee_hgt(i) = 0._r8
               n_levee_resv_clash = n_levee_resv_clash + 1
            ENDIF
         ENDDO
         IF (n_levee_resv_clash > 0) THEN
            write(*,'(A,I0,A,I0,A)') &
               'INFO levee_init: forced has_levee=.false. on ', &
               n_levee_resv_clash, ' of ', numucat, &
               ' cells where lake_type==2 (reservoir) and levee_frc<1 collide; check input levee_frc / dam_seq consistency.'
         ENDIF
      ENDIF

      ! --- Build a CaMa-compatible total storage curve.
      ! Network flpstomax is floodplain-only; levee thresholds need river
      ! storage plus the river-width prism above bankfull.
      DO i = 1, numucat
         nlfp    = floodplain_curve(i)%nlfp
         rivlen  = topo_rivlen(i)
         rivwth  = topo_rivwth(i)
         catarea = topo_area(i)

         IF (rivlen <= 0._r8 .or. nlfp <= 0 .or. catarea <= 0._r8) CYCLE

         dwth_inc = catarea / (rivlen * nlfp)
         dsto_fil = floodplain_curve(i)%rivstomax
         dhgtpre  = 0._r8

         DO j = 1, nlfp
            dhgtnow = floodplain_curve(i)%flphgt(j) - dhgtpre
            IF (dhgtnow > 0._r8 .and. dwth_inc > 0._r8) THEN
               levee_fldgrd(j, i) = dhgtnow / dwth_inc
            ELSE
               levee_fldgrd(j, i) = 0._r8
               dhgtnow = 0._r8
            ENDIF

            dsto_add = (rivwth + dwth_inc * (real(j, r8) - 0.5_r8)) * dhgtnow * rivlen
            levee_fldstomax(j, i) = dsto_fil + dsto_add
            dsto_fil = levee_fldstomax(j, i)
            dhgtpre  = floodplain_curve(i)%flphgt(j)
         ENDDO
      ENDDO

      ! --- Derive parameters for levee cells ---
      DO i = 1, numucat
         IF (.not. has_levee(i)) CYCLE

         nlfp    = floodplain_curve(i)%nlfp
         rivlen  = topo_rivlen(i)
         rivwth  = topo_rivwth(i)
         catarea = topo_area(i)

         IF (rivlen <= 0. .or. nlfp <= 0 .or. catarea <= 0.) THEN
            has_levee(i) = .false.
            CYCLE
         ENDIF

         ! Distance from river center to levee
         levee_dst(i) = levee_frc(i) * (catarea / rivlen)

         ! Width of each floodplain layer
         dwth_inc = catarea / (rivlen * nlfp)

         ! Which floodplain layer contains the levee
         ilev = INT(levee_frc(i) * nlfp) + 1
         ilev = MIN(MAX(ilev, 1), nlfp)

         ! Levee base height: interpolate within floodplain layer
         IF (ilev >= 2) THEN
            dhgtpre  = floodplain_curve(i)%flphgt(ilev-1)
            dwth_fil = dwth_inc * (ilev - 1)
         ELSE
            dhgtpre  = 0.
            dwth_fil = 0.
         ENDIF

         dwth_add = levee_dst(i) - dwth_fil
         dwth_add = MAX(dwth_add, 0._r8)

         fldgrd = levee_fldgrd(ilev, i)

         levee_bashgt(i) = dhgtpre + dwth_add * fldgrd

         ! Ensure levee crest >= base
         levee_hgt(i) = MAX(levee_hgt(i), levee_bashgt(i))

         ! Storage at levee base (river + floodplain layers up to levee position)
         IF (ilev >= 2) THEN
            dsto_fil = levee_fldstomax(ilev-1, i)
         ELSE
            dsto_fil = floodplain_curve(i)%rivstomax
         ENDIF
         levee_bassto(i) = dsto_fil

         IF (dwth_add > 0._r8 .and. fldgrd > 0._r8) THEN
            levee_bassto(i) = dsto_fil &
               + (dwth_add*0.5 + dwth_fil + rivwth) * (dwth_add * fldgrd) * rivlen
         ENDIF

         ! Storage at levee crest (river-side only)
         dhgtdif = levee_hgt(i) - levee_bashgt(i)
         levee_topsto(i) = levee_bassto(i) &
            + (levee_dst(i) + rivwth) * dhgtdif * rivlen

         ! Storage when both sides are filled to the levee crest.
         levee_filsto(i) = levee_total_volume_from_depth(i, floodplain_curve(i)%rivhgt + levee_hgt(i))

         ! Ensure monotonicity: bassto <= topsto <= filsto
         levee_topsto(i) = MAX(levee_topsto(i), levee_bassto(i))
         levee_filsto(i) = MAX(levee_filsto(i), levee_topsto(i))

      ENDDO

   END SUBROUTINE levee_init


   ! =========================================================================
   SUBROUTINE levee_fldstg (i, vol_total, wdsrf, levsto_out, levdph_out, fldfrc)
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: topo_rivlen, topo_rivwth, topo_area, &
      floodplain_curve
   IMPLICIT NONE

   integer,  intent(in)  :: i           ! unit catchment index
   real(r8), intent(in)  :: vol_total   ! total water volume [m^3]
   real(r8), intent(out) :: wdsrf       ! water depth for river+unprotected side [m]
   real(r8), intent(out) :: levsto_out  ! protected-side storage [m^3]
   real(r8), intent(out) :: levdph_out  ! protected-side water depth [m]
   real(r8), intent(out) :: fldfrc      ! flooded fraction of floodplain [0-1]

   ! Local variables
   real(r8) :: rivstomax, rivlen, rivwth
   real(r8) :: flddph, rivsto, dsto_fil, dsto_add
   real(r8) :: dwth_inc, dwth_fil, dwth_add, ddph_fil, ddph_add
   real(r8) :: fldgrd, fldsto_unprot
   integer  :: ilev, j, nlfp

      rivstomax = floodplain_curve(i)%rivstomax

      ! --- No levee: standard behavior ---
      IF (.not. has_levee(i)) THEN
         wdsrf      = floodplain_curve(i)%depth(vol_total)
         levsto_out = 0.
         levdph_out = 0.
         IF (topo_area(i) > 0.) THEN
            fldfrc = floodplain_curve(i)%floodarea(wdsrf) / topo_area(i)
         ELSE
            fldfrc = 0.
         ENDIF
         RETURN
      ENDIF

      rivlen = topo_rivlen(i)
      rivwth = topo_rivwth(i)
      nlfp   = floodplain_curve(i)%nlfp
      IF (rivlen > 0._r8 .and. nlfp > 0) THEN
         dwth_inc = topo_area(i) / (rivlen * nlfp)
      ELSE
         dwth_inc = 0._r8
      ENDIF

      IF (vol_total <= rivstomax) THEN
         !--- Case 0: water only in river channel ---
         IF (floodplain_curve(i)%rivare > 0.) THEN
            wdsrf = vol_total / floodplain_curve(i)%rivare
         ELSE
            wdsrf = 0.
         ENDIF
         levsto_out = 0.
         levdph_out = 0.
         fldfrc     = 0.

      ELSE IF (vol_total < levee_bassto(i)) THEN
         !--- Case 1: water below levee base, unprotected side only ---
         j = 1
         dsto_fil = rivstomax
         dwth_fil = rivwth
         ddph_fil = 0._r8

         DO WHILE (j <= nlfp)
            IF (vol_total <= levee_fldstomax(j, i)) EXIT
            dsto_fil = levee_fldstomax(j, i)
            dwth_fil = dwth_fil + dwth_inc
            ddph_fil = floodplain_curve(i)%flphgt(j)
            j = j + 1
         ENDDO

         IF (j <= nlfp) THEN
            dsto_add = vol_total - dsto_fil
            fldgrd = levee_fldgrd(j, i)
            IF (fldgrd > 0._r8) THEN
               dwth_add = -dwth_fil + SQRT(MAX(dwth_fil**2 + 2._r8*dsto_add / rivlen / fldgrd, 0._r8))
            ELSE
               dwth_add = 0._r8
            ENDIF
            flddph = ddph_fil + fldgrd * dwth_add
         ELSE
            dsto_add = vol_total - dsto_fil
            dwth_add = 0._r8
            IF (dwth_fil > 0._r8) THEN
               flddph = ddph_fil + dsto_add / dwth_fil / rivlen
            ELSE
               flddph = ddph_fil
            ENDIF
         ENDIF

         wdsrf      = floodplain_curve(i)%rivhgt + flddph
         levsto_out = 0.
         levdph_out = 0.
         fldfrc = (-rivwth + dwth_fil + dwth_add) / (dwth_inc * nlfp)
         fldfrc = MIN(MAX(fldfrc, 0._r8), levee_frc(i))

      ELSE IF (vol_total < levee_topsto(i)) THEN
         !--- Case 2: water between levee base and crest, unprotected side ---
         IF ((levee_dst(i) + rivwth) * rivlen > 0.) THEN
            flddph = levee_bashgt(i) + (vol_total - levee_bassto(i)) &
               / ((levee_dst(i) + rivwth) * rivlen)
         ELSE
            flddph = levee_bashgt(i)
         ENDIF
         wdsrf      = floodplain_curve(i)%rivhgt + flddph
         levsto_out = 0.
         levdph_out = 0.
         fldfrc     = levee_frc(i)

      ELSE IF (vol_total < levee_filsto(i)) THEN
         !--- Case 3: levee overtopping, both sides filling ---
         wdsrf = floodplain_curve(i)%rivhgt + levee_hgt(i)
         rivsto = rivstomax + rivlen * rivwth * levee_hgt(i)
         fldsto_unprot = MAX(levee_topsto(i) - rivsto, 0._r8)
         levsto_out = MAX(vol_total - rivsto - fldsto_unprot, 0._r8)

         ilev = INT(levee_frc(i) * nlfp) + 1
         ilev = MIN(MAX(ilev, 1), nlfp)
         dsto_fil = levee_topsto(i)
         dwth_fil = 0._r8
         ddph_fil = 0._r8

         j = ilev
         DO WHILE (j <= nlfp)
            dsto_add = (levee_dst(i) + rivwth) * (levee_hgt(i) - floodplain_curve(i)%flphgt(j)) * rivlen
            IF (vol_total < levee_fldstomax(j, i) + dsto_add) EXIT
            dsto_fil = levee_fldstomax(j, i) + dsto_add
            dwth_fil = dwth_inc * j - levee_dst(i)
            ddph_fil = floodplain_curve(i)%flphgt(j) - levee_bashgt(i)
            j = j + 1
         ENDDO

         IF (j <= nlfp) THEN
            dsto_add = vol_total - dsto_fil
            fldgrd = levee_fldgrd(j, i)
            IF (fldgrd > 0._r8) THEN
               dwth_add = -dwth_fil + SQRT(MAX(dwth_fil**2 + 2._r8*dsto_add / rivlen / fldgrd, 0._r8))
            ELSE
               dwth_add = 0._r8
            ENDIF
            ddph_add = dwth_add * fldgrd
            levdph_out = levee_bashgt(i) + ddph_fil + ddph_add
            fldfrc = (dwth_fil + dwth_add + levee_dst(i)) / (dwth_inc * nlfp)
            fldfrc = MAX(fldfrc, 0._r8)
            fldfrc = MIN(fldfrc, 1._r8)
         ELSE
            dsto_add = vol_total - dsto_fil
            IF (dwth_fil > 0._r8) THEN
               ddph_add = dsto_add / dwth_fil / rivlen
            ELSE
               ddph_add = 0._r8
            ENDIF
            levdph_out = levee_bashgt(i) + ddph_fil + ddph_add
            fldfrc = 1.0_r8
         ENDIF

      ELSE
         !--- Case 4: water above levee crest, both sides fully flooded ---
         wdsrf    = levee_total_depth(i, vol_total)
         flddph   = wdsrf - floodplain_curve(i)%rivhgt
         levdph_out = flddph
         fldfrc = 1.0

         rivsto = rivstomax + rivlen * rivwth * flddph
         dsto_add = (flddph - levee_hgt(i)) * (levee_dst(i) + rivwth) * rivlen
         fldsto_unprot = MAX(levee_topsto(i) + dsto_add - rivsto, 0._r8)
         levsto_out = MAX(vol_total - rivsto - fldsto_unprot, 0._r8)
      ENDIF

   END SUBROUTINE levee_fldstg


   ! =========================================================================
   SUBROUTINE levee_apply_protected_flux (i, protected_hflux, dt, clipped_volume)
   ! =========================================================================

   IMPLICIT NONE

   integer,  intent(in)  :: i
   real(r8), intent(in)  :: protected_hflux
   real(r8), intent(in)  :: dt
   real(r8), intent(out), optional :: clipped_volume

   real(r8) :: protected_raw

      IF (present(clipped_volume)) clipped_volume = 0._r8
      IF (.not. allocated(levsto)) RETURN
      IF (i < 1 .or. i > size(levsto)) RETURN
      IF (dt <= 0._r8) RETURN

      protected_raw = levsto(i) - protected_hflux * dt
      IF (present(clipped_volume)) THEN
         IF (protected_raw < 0._r8) clipped_volume = -protected_raw
      ENDIF
      levsto(i) = max(protected_raw, 0._r8)

   END SUBROUTINE levee_apply_protected_flux


   ! =========================================================================
   SUBROUTINE levee_repartition_storage (i, visible_volume, wdsrf, fldfrc, &
      vis_vol_bef, levsto_bef)
   ! =========================================================================

   IMPLICIT NONE

   integer,  intent(in)    :: i
   real(r8), intent(inout) :: visible_volume
   real(r8), intent(out)   :: wdsrf
   real(r8), intent(out)   :: fldfrc
   real(r8), intent(out), optional :: vis_vol_bef
   real(r8), intent(out), optional :: levsto_bef

   real(r8) :: vol_total

      IF (present(vis_vol_bef)) vis_vol_bef = visible_volume
      IF (present(levsto_bef))  levsto_bef  = levsto(i)

      vol_total = visible_volume + levsto(i)
      CALL levee_fldstg(i, vol_total, wdsrf, levsto(i), levdph(i), fldfrc)
      visible_volume = vol_total - levsto(i)

   END SUBROUTINE levee_repartition_storage


   ! =========================================================================
   FUNCTION levee_total_volume_from_depth (i, wdsrf) RESULT(volume)
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: topo_rivlen, topo_rivwth, topo_area, &
      floodplain_curve
   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: wdsrf
   real(r8)             :: volume

   integer  :: j, nlfp
   real(r8) :: rivlen, rivwth, dwth_inc
   real(r8) :: flddph, dhgtpre, dhgtnow, dwth_fil, dwth_add, fldgrd
   real(r8) :: dsto_fil, dsto_add

      IF (wdsrf <= 0._r8) THEN
         volume = 0._r8
         RETURN
      ENDIF

      IF (wdsrf <= floodplain_curve(i)%rivhgt) THEN
         volume = wdsrf * floodplain_curve(i)%rivare
         RETURN
      ENDIF

      rivlen = topo_rivlen(i)
      rivwth = topo_rivwth(i)
      nlfp   = floodplain_curve(i)%nlfp
      IF (rivlen <= 0._r8 .or. topo_area(i) <= 0._r8 .or. nlfp <= 0) THEN
         volume = floodplain_curve(i)%volume(wdsrf)
         RETURN
      ENDIF

      flddph   = wdsrf - floodplain_curve(i)%rivhgt
      dwth_inc = topo_area(i) / (rivlen * nlfp)

      j         = 1
      dsto_fil  = floodplain_curve(i)%rivstomax
      dwth_fil  = rivwth
      dhgtpre   = 0._r8

      DO WHILE (j <= nlfp)
         IF (flddph <= floodplain_curve(i)%flphgt(j)) EXIT
         dsto_fil = levee_fldstomax(j, i)
         dwth_fil = rivwth + dwth_inc * j
         dhgtpre  = floodplain_curve(i)%flphgt(j)
         j = j + 1
      ENDDO

      dhgtnow = MAX(flddph - dhgtpre, 0._r8)
      IF (j <= nlfp) THEN
         fldgrd = levee_fldgrd(j, i)
         IF (fldgrd > 0._r8) THEN
            dwth_add = dhgtnow / fldgrd
         ELSE
            dwth_add = 0._r8
         ENDIF
         dsto_add = (dwth_fil + 0.5_r8 * dwth_add) * dhgtnow * rivlen
      ELSE
         dsto_add = dwth_fil * dhgtnow * rivlen
      ENDIF

      volume = MAX(dsto_fil + dsto_add, 0._r8)

   END FUNCTION levee_total_volume_from_depth


   ! =========================================================================
   FUNCTION levee_total_depth (i, volume) RESULT(wdsrf)
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: topo_rivlen, topo_rivwth, topo_area, &
      floodplain_curve
   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: volume
   real(r8)             :: wdsrf

   integer  :: j, nlfp
   real(r8) :: rivlen, rivwth, dwth_inc
   real(r8) :: dsto_fil, dsto_add, dwth_fil, dwth_add, ddph_fil, fldgrd

      IF (volume <= 0._r8) THEN
         wdsrf = 0._r8
         RETURN
      ENDIF

      IF (volume <= floodplain_curve(i)%rivstomax) THEN
         IF (floodplain_curve(i)%rivare > 0._r8) THEN
            wdsrf = volume / floodplain_curve(i)%rivare
         ELSE
            wdsrf = 0._r8
         ENDIF
         RETURN
      ENDIF

      rivlen = topo_rivlen(i)
      rivwth = topo_rivwth(i)
      nlfp   = floodplain_curve(i)%nlfp
      IF (rivlen <= 0._r8 .or. topo_area(i) <= 0._r8 .or. nlfp <= 0) THEN
         wdsrf = floodplain_curve(i)%depth(volume)
         RETURN
      ENDIF

      dwth_inc = topo_area(i) / (rivlen * nlfp)

      j         = 1
      dsto_fil  = floodplain_curve(i)%rivstomax
      dwth_fil  = rivwth
      ddph_fil  = 0._r8

      DO WHILE (j <= nlfp)
         IF (volume <= levee_fldstomax(j, i)) EXIT
         dsto_fil = levee_fldstomax(j, i)
         dwth_fil = rivwth + dwth_inc * j
         ddph_fil = floodplain_curve(i)%flphgt(j)
         j = j + 1
      ENDDO

      dsto_add = MAX(volume - dsto_fil, 0._r8)
      IF (j <= nlfp) THEN
         fldgrd = levee_fldgrd(j, i)
         IF (fldgrd > 0._r8) THEN
            dwth_add = -dwth_fil + SQRT(MAX(dwth_fil**2 + 2._r8*dsto_add / rivlen / fldgrd, 0._r8))
            wdsrf = floodplain_curve(i)%rivhgt + ddph_fil + fldgrd * dwth_add
         ELSE
            wdsrf = floodplain_curve(i)%rivhgt + ddph_fil
         ENDIF
      ELSE
         IF (dwth_fil > 0._r8 .and. rivlen > 0._r8) THEN
            wdsrf = floodplain_curve(i)%rivhgt + ddph_fil + dsto_add / dwth_fil / rivlen
         ELSE
            wdsrf = floodplain_curve(i)%rivhgt + ddph_fil
         ENDIF
      ENDIF

   END FUNCTION levee_total_depth


   ! =========================================================================
   FUNCTION levee_visible_volume_from_stage (i, wdsrf, levsto_in) RESULT(visible_volume)
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: topo_rivlen, topo_rivwth, floodplain_curve
   IMPLICIT NONE

   integer,  intent(in) :: i
   real(r8), intent(in) :: wdsrf
   real(r8), intent(in) :: levsto_in
   real(r8)             :: visible_volume

   real(r8) :: flddph

      IF (.not. allocated(has_levee)) THEN
         visible_volume = floodplain_curve(i)%volume(wdsrf)
         RETURN
      ENDIF

      IF (i < 1 .or. i > size(has_levee) .or. .not. has_levee(i)) THEN
         visible_volume = floodplain_curve(i)%volume(wdsrf)
         RETURN
      ENDIF

      flddph = wdsrf - floodplain_curve(i)%rivhgt

      IF (flddph <= levee_bashgt(i)) THEN
         visible_volume = levee_total_volume_from_depth(i, wdsrf)
      ELSEIF (flddph <= levee_hgt(i) + 1.e-6_r8) THEN
         visible_volume = levee_bassto(i) + (levee_dst(i) + topo_rivwth(i)) &
            * MAX(flddph - levee_bashgt(i), 0._r8) * topo_rivlen(i)
      ELSE
         visible_volume = levee_total_volume_from_depth(i, wdsrf) - MAX(levsto_in, 0._r8)
      ENDIF

      visible_volume = MAX(visible_volume, 0._r8)

   END FUNCTION levee_visible_volume_from_stage


   ! =========================================================================
   SUBROUTINE read_levee_restart (file_restart, fold_protected_to_visible, &
      volwater_ucat_io, volwater_ucat_valid_io)
   ! =========================================================================
   ! Read levee state from restart file. Called AFTER levee_init() so that
   ! levsto is already allocated. If runtime levee is disabled while the
   ! restart contains protected storage, fold that protected mass into the
   ! TimeVars-owned visible routing volume before zeroing levsto.
   ! =========================================================================

   USE MOD_NetCDFSerial,        only: ncio_var_exist
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart
   logical,  intent(in),    optional :: fold_protected_to_visible
   real(r8), allocatable, intent(inout), optional :: volwater_ucat_io(:)
   logical,  intent(inout), optional :: volwater_ucat_valid_io
   integer :: has_restart_var(2)
   logical :: do_fold
   real(r8) :: folded_volume

      ! NOTE: vector_read_and_scatter contains mpi_barrier(p_comm_glb),
      ! so ALL processes must participate. Never early-return for non-workers.
      has_restart_var = 0
      IF (p_is_master) THEN
         IF (ncio_var_exist(file_restart, 'levsto', readflag = .false.)) has_restart_var(1) = 1
         IF (ncio_var_exist(file_restart, 'volwater_ucat', readflag = .false.)) has_restart_var(2) = 1
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (has_restart_var, 2, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      IF (has_restart_var(1) == 1) THEN
         CALL vector_read_and_scatter (file_restart, levsto, numucat, 'levsto', ucat_data_address)
      ENDIF

      IF (has_restart_var(2) == 1) THEN
         IF (present(volwater_ucat_io)) THEN
            CALL vector_read_and_scatter (file_restart, volwater_ucat_io, numucat, 'volwater_ucat', ucat_data_address)
         ENDIF
         IF (present(volwater_ucat_valid_io)) volwater_ucat_valid_io = .true.
      ENDIF

      do_fold = .false.
      IF (present(fold_protected_to_visible)) do_fold = fold_protected_to_visible
      IF (do_fold .and. (.not. present(volwater_ucat_io))) THEN
         IF (p_is_master) THEN
            write(*,'(A)') 'ERROR read_levee_restart: fold_protected_to_visible requires volwater_ucat_io.'
         ENDIF
         do_fold = .false.
      ENDIF
      IF (do_fold .and. has_restart_var(1) == 1) THEN
         folded_volume = 0._r8
         IF (p_is_worker .and. present(volwater_ucat_io)) THEN
            IF (numucat > 0) THEN
               folded_volume = sum(max(levsto, 0._r8))
               volwater_ucat_io = volwater_ucat_io + max(levsto, 0._r8)
               levsto = 0._r8
            ENDIF
         ENDIF
         IF (present(volwater_ucat_valid_io)) volwater_ucat_valid_io = .true.
#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, folded_volume, 1, MPI_REAL8, MPI_SUM, p_comm_glb, p_err)
#endif
         IF (p_is_master .and. folded_volume > 0._r8) THEN
            write(*,'(A,ES12.4,A)') 'WARNING read_levee_restart: DEF_USE_LEVEE is false; folded ', &
               folded_volume, ' m^3 protected levee storage into volwater_ucat.'
         ENDIF
      ENDIF

   END SUBROUTINE read_levee_restart


   ! =========================================================================
   SUBROUTINE levee_final ()
   ! =========================================================================

   IMPLICIT NONE

      IF (allocated(has_levee    )) deallocate(has_levee    )
      IF (allocated(levee_frc    )) deallocate(levee_frc    )
      IF (allocated(levee_hgt    )) deallocate(levee_hgt    )
      IF (allocated(levee_dst    )) deallocate(levee_dst    )
      IF (allocated(levee_bashgt )) deallocate(levee_bashgt )
      IF (allocated(levee_bassto )) deallocate(levee_bassto )
      IF (allocated(levee_topsto )) deallocate(levee_topsto )
      IF (allocated(levee_filsto )) deallocate(levee_filsto )
      IF (allocated(levee_fldstomax)) deallocate(levee_fldstomax)
      IF (allocated(levee_fldgrd    )) deallocate(levee_fldgrd    )
      IF (allocated(levsto        )) deallocate(levsto        )
      IF (allocated(levdph        )) deallocate(levdph        )

   END SUBROUTINE levee_final

END MODULE MOD_Grid_RiverLakeLevee
#endif
