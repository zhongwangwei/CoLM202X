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

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   ! ----- Levee parameters (derived in levee_init) -----
   logical,  allocatable :: has_levee     (:)   ! true where levee exists
   real(r8), allocatable :: levee_frc     (:)   ! unprotected fraction [0-1]
   real(r8), allocatable :: levee_hgt     (:)   ! levee crest height above riverbed [m]
   real(r8), allocatable :: levee_dst     (:)   ! distance from river center to levee [m]
   real(r8), allocatable :: levee_bashgt  (:)   ! ground elevation at levee position [m]
   real(r8), allocatable :: levee_bassto  (:)   ! storage at levee base [m^3]
   real(r8), allocatable :: levee_topsto  (:)   ! storage at levee crest, river-side [m^3]
   real(r8), allocatable :: levee_filsto  (:)   ! storage when both sides filled to crest [m^3]

   ! ----- State variables -----
   real(r8), allocatable :: levsto        (:)   ! protected-side water storage [m^3]
   real(r8), allocatable :: volwater_ucat (:)   ! non-protected-side water volume [m^3]
   logical               :: volwater_ucat_valid = .false.  ! true after restore from restart or first routing call

   ! ----- Diagnostic -----
   real(r8), allocatable :: levdph        (:)   ! protected-side water depth [m]

   PUBLIC :: levee_init
   PUBLIC :: read_levee_restart
   PUBLIC :: levee_fldstg
   PUBLIC :: levee_final

CONTAINS

   ! =========================================================================
   SUBROUTINE levee_init ()
   ! =========================================================================

   USE MOD_Grid_RiverLakeNetwork, only: numucat, topo_rivlen, topo_rivwth, &
      topo_area, topo_rivstomax, floodplain_curve, levee_frc_data, levee_hgt_data
   IMPLICIT NONE

   ! Local variables
   integer  :: i, j, ilev, nlfp
   real(r8) :: rivlen, rivwth, catarea, dwth_inc
   real(r8) :: dhgtpre, dwth_fil, dwth_add, fldgrd, dhgtdif
   real(r8) :: dhgtnow, dsto_fil, dsto_add

      volwater_ucat_valid = .false.

      IF (.not. p_is_worker) RETURN
      IF (numucat <= 0) RETURN

      ! --- Allocate arrays ---
      allocate (has_levee    (numucat))
      allocate (levee_frc    (numucat))
      allocate (levee_hgt    (numucat))
      allocate (levee_dst    (numucat))
      allocate (levee_bashgt (numucat))
      allocate (levee_bassto (numucat))
      allocate (levee_topsto (numucat))
      allocate (levee_filsto (numucat))
      allocate (levsto        (numucat))
      allocate (volwater_ucat (numucat))
      allocate (levdph        (numucat))

      ! --- Copy from Network data ---
      levee_frc = levee_frc_data
      levee_hgt = levee_hgt_data

      ! --- Validate ---
      DO i = 1, numucat
         IF (levee_hgt(i) <= 0. .or. levee_frc(i) >= 1. .or. levee_frc(i) < 0.) THEN
            levee_hgt(i) = 0.
            levee_frc(i) = 1.0
         ENDIF
         levee_frc(i) = max(0._r8, min(1._r8, levee_frc(i)))
      ENDDO

      has_levee = (levee_frc < 1.0_r8)

      ! --- Initialize state ---
      levsto        = 0.
      volwater_ucat = 0.
      levdph        = 0.

      ! --- Initialize derived parameters for non-levee cells ---
      levee_dst    = 0.
      levee_bashgt = 0.
      levee_bassto = 0.
      levee_topsto = 0.
      levee_filsto = 0.

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
         ilev = MIN(ilev, nlfp)

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

         IF (floodplain_curve(i)%flphgt(ilev) > floodplain_curve(i)%flphgt(ilev-1)) THEN
            fldgrd = (floodplain_curve(i)%flphgt(ilev) - floodplain_curve(i)%flphgt(ilev-1)) / dwth_inc
         ELSE
            fldgrd = 0.
         ENDIF

         levee_bashgt(i) = dhgtpre + dwth_add * fldgrd

         ! Ensure levee crest >= base
         levee_hgt(i) = MAX(levee_hgt(i), levee_bashgt(i))

         ! Storage at levee base (river + floodplain layers up to levee position)
         levee_bassto(i) = floodplain_curve(i)%rivstomax
         DO j = 1, ilev-1
            levee_bassto(i) = levee_bassto(i) &
               + floodplain_curve(i)%flpstomax(j) - floodplain_curve(i)%flpstomax(j-1)
         ENDDO
         ! Partial layer contribution
         IF (dwth_add > 0. .and. fldgrd > 0.) THEN
            levee_bassto(i) = levee_bassto(i) &
               + (dwth_add*0.5 + dwth_fil + rivwth) * (dwth_add * fldgrd) * rivlen
         ENDIF

         ! Storage at levee crest (river-side only)
         dhgtdif = levee_hgt(i) - levee_bashgt(i)
         levee_topsto(i) = levee_bassto(i) &
            + (levee_dst(i) + rivwth) * dhgtdif * rivlen

         ! Storage when both sides are filled to the levee crest.
         ! Follow CaMa-Flood's layer-wise integration instead of using the
         ! aggregated volume(depth) curve directly.
         j = 1
         dsto_fil = floodplain_curve(i)%rivstomax
         dwth_fil = rivwth
         dhgtpre  = 0._r8

         DO WHILE (j <= nlfp)
            IF (levee_hgt(i) <= floodplain_curve(i)%flphgt(j)) EXIT
            dsto_fil = floodplain_curve(i)%rivstomax + floodplain_curve(i)%flpstomax(j)
            dwth_fil = dwth_fil + dwth_inc
            dhgtpre  = floodplain_curve(i)%flphgt(j)
            j = j + 1
         ENDDO

         IF (j <= nlfp) THEN
            dhgtnow = levee_hgt(i) - dhgtpre
            fldgrd = (floodplain_curve(i)%flphgt(j) - floodplain_curve(i)%flphgt(j-1)) / dwth_inc
            IF (fldgrd > 0._r8) THEN
               dwth_add = dhgtnow / fldgrd
            ELSE
               dwth_add = 0._r8
            ENDIF
            dsto_add = (dwth_add*0.5_r8 + dwth_fil) * dhgtnow * rivlen
            levee_filsto(i) = dsto_fil + dsto_add
         ELSE
            dhgtnow = levee_hgt(i) - dhgtpre
            dsto_add = dwth_fil * dhgtnow * rivlen
            levee_filsto(i) = dsto_fil + dsto_add
         ENDIF

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
            IF (vol_total <= floodplain_curve(i)%rivstomax + floodplain_curve(i)%flpstomax(j)) EXIT
            dsto_fil = floodplain_curve(i)%rivstomax + floodplain_curve(i)%flpstomax(j)
            dwth_fil = dwth_fil + dwth_inc
            fldgrd = (floodplain_curve(i)%flphgt(j) - floodplain_curve(i)%flphgt(j-1)) / dwth_inc
            ddph_fil = ddph_fil + fldgrd * dwth_inc
            j = j + 1
         ENDDO

         IF (j <= nlfp) THEN
            dsto_add = vol_total - dsto_fil
            fldgrd = (floodplain_curve(i)%flphgt(j) - floodplain_curve(i)%flphgt(j-1)) / dwth_inc
            IF (fldgrd > 0._r8) THEN
               dwth_add = -dwth_fil + SQRT(dwth_fil**2 + 2._r8*dsto_add / rivlen / fldgrd)
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
         dsto_fil = levee_topsto(i)
         dwth_fil = 0._r8
         ddph_fil = 0._r8

         j = ilev
         DO WHILE (j <= nlfp)
            dsto_add = (levee_dst(i) + rivwth) * (levee_hgt(i) - floodplain_curve(i)%flphgt(j)) * rivlen
            IF (vol_total < floodplain_curve(i)%rivstomax + floodplain_curve(i)%flpstomax(j) + dsto_add) EXIT
            dsto_fil = floodplain_curve(i)%rivstomax + floodplain_curve(i)%flpstomax(j) + dsto_add
            dwth_fil = dwth_inc * j - levee_dst(i)
            ddph_fil = floodplain_curve(i)%flphgt(j) - levee_bashgt(i)
            j = j + 1
         ENDDO

         IF (j <= nlfp) THEN
            dsto_add = vol_total - dsto_fil
            fldgrd = (floodplain_curve(i)%flphgt(j) - floodplain_curve(i)%flphgt(j-1)) / dwth_inc
            IF (fldgrd > 0._r8) THEN
               dwth_add = -dwth_fil + SQRT(dwth_fil**2 + 2._r8*dsto_add / rivlen / fldgrd)
            ELSE
               dwth_add = 0._r8
            ENDIF
            ddph_add = dwth_add * fldgrd
            levdph_out = levee_bashgt(i) + ddph_fil + ddph_add
            fldfrc = (dwth_fil + levee_dst(i)) / (dwth_inc * nlfp)
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
         wdsrf    = floodplain_curve(i)%depth(vol_total)
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
   SUBROUTINE read_levee_restart (file_restart)
   ! =========================================================================
   ! Read levee state from restart file. Called AFTER levee_init() so that
   ! levsto is already allocated. If the restart file doesn't contain 'levsto',
   ! the array keeps its zero initialization from levee_init().
   ! =========================================================================

   USE MOD_NetCDFSerial,        only: ncio_var_exist
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

      ! NOTE: vector_read_and_scatter contains mpi_barrier(p_comm_glb),
      ! so ALL processes must participate. Never early-return for non-workers.
      IF (ncio_var_exist(file_restart, 'levsto')) THEN
         CALL vector_read_and_scatter (file_restart, levsto, numucat, 'levsto', ucat_data_address)
      ENDIF

      IF (ncio_var_exist(file_restart, 'volwater_ucat')) THEN
         CALL vector_read_and_scatter (file_restart, volwater_ucat, numucat, 'volwater_ucat', ucat_data_address)
         volwater_ucat_valid = .true.
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
      IF (allocated(levsto        )) deallocate(levsto        )
      IF (allocated(volwater_ucat)) deallocate(volwater_ucat)
      IF (allocated(levdph        )) deallocate(levdph        )
      volwater_ucat_valid = .false.

   END SUBROUTINE levee_final

END MODULE MOD_Grid_RiverLakeLevee
#endif
