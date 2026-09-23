#include <define.h>

MODULE MOD_HtopReadin

   USE MOD_Precision
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: HTOP_readin

CONTAINS

   SUBROUTINE HTOP_readin (dir_landdata, lc_year)

! ===========================================================
! Read in the canopy tree top height
! ===========================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Const_LC
   USE MOD_Const_PFT
   USE MOD_Vars_TimeInvariants
   USE MOD_LandPatch
   USE MOD_Namelist, only: DEF_Interception_scheme
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeInvariants
   USE MOD_Vars_PFTimeVariables
#endif
   USE MOD_NetCDFVector
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year    ! which year of land cover data used
   character(len=256), intent(in) :: dir_landdata

   ! Local Variables
   character(len=256) :: c
   character(len=256) :: landdir, cstructdir, lndname, cyear
   integer :: i,j,t,p,ps,pe,m,n,npatch
   logical :: site_struct_ok
   integer :: site_required_count, site_valid_count
   integer :: canopy_counts(2), icanopy
   integer :: usgs_counts(4)

   real(r8), allocatable :: htoplc  (:)
   real(r8), allocatable :: htoppft (:)
   character(len=16) :: struct_names(3)
   logical :: struct_present(3)

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_landdata) // '/htop/' // trim(cyear)
      cstructdir = trim(dir_landdata) // '/cstructure/' // trim(cyear)

#if (defined LULC_USGS || (defined LULC_IGBP && !defined LULC_IGBP_PFT && !defined LULC_IGBP_PC))
#ifndef SinglePoint
      struct_names = (/ 'ncd_patches    ', 'ncw_patches    ', 'bcw_patches    ' /)
      lndname = trim(cstructdir)//'/ncd_patches.nc'
      struct_present(1) = ncio_vector_var_present(lndname, struct_names(1), landpatch)
      lndname = trim(cstructdir)//'/ncw_patches.nc'
      struct_present(2) = ncio_vector_var_present(lndname, struct_names(2), landpatch)
      lndname = trim(cstructdir)//'/bcw_patches.nc'
      struct_present(3) = ncio_vector_var_present(lndname, struct_names(3), landpatch)
      IF (any(struct_present) .and. .not. all(struct_present)) THEN
         IF (p_is_master) write(*,'(A)') &
            'ERROR: CoLM2024 canopy-structure input is incomplete; ' // &
            'need ncd_patches, ncw_patches, and bcw_patches.'
         CALL CoLM_stop()
      ENDIF
      IF (all(struct_present)) THEN
         lndname = trim(cstructdir)//'/ncd_patches.nc'
         CALL ncio_read_vector (lndname, 'ncd_patches', landpatch, ncd)
         lndname = trim(cstructdir)//'/ncw_patches.nc'
         CALL ncio_read_vector (lndname, 'ncw_patches', landpatch, ncw)
         lndname = trim(cstructdir)//'/bcw_patches.nc'
         CALL ncio_read_vector (lndname, 'bcw_patches', landpatch, bcw)
      ELSEIF (DEF_Interception_scheme == 8) THEN
#ifdef LULC_IGBP
         IF (p_is_master) write(*,'(A)') &
            'Warning: CoLM2024 requested but canopy-structure inputs are missing; ' // &
            'using CoLM2014 interception. Run mksrfdata with canopy_data to enable CoLM2024.'
         DEF_Interception_scheme = 1
#endif
      ENDIF
#endif
#endif


#ifdef LULC_USGS

#ifdef SinglePoint
      IF (p_is_worker) THEN
         ncd = SITE_ncd
         ncw = SITE_ncw
         bcw = SITE_bcw
      ENDIF
#endif

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            m = patchclass(npatch)

            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)

         ENDDO
      ENDIF

      IF (DEF_Interception_scheme == 8) THEN
         ! Forest classes need morphology; shrubs use only wind.  A missing
         ! forest value downgrades the entire run, never only that patch.
         usgs_counts = 0
         IF (p_is_worker) THEN
            DO icanopy = 1, numpatch
               IF (patchtype(icanopy) /= 0) CYCLE
               SELECT CASE (patchclass(icanopy))
               CASE (11,13)
                  usgs_counts(1) = usgs_counts(1) + 1
                  IF (.not. all(ieee_is_finite([bcw(icanopy), htop(icanopy)]))) CYCLE
                  IF (bcw(icanopy) <= 0._r8 .or. bcw(icanopy) >= 1000._r8 .or. &
                      htop(icanopy) <= 0._r8 .or. htop(icanopy) >= 1000._r8) CYCLE
               CASE (12,14)
                  usgs_counts(1) = usgs_counts(1) + 1
                  IF (.not. all(ieee_is_finite([ncd(icanopy), ncw(icanopy)]))) CYCLE
                  IF (ncd(icanopy) <= 0._r8 .or. ncd(icanopy) >= 1000._r8 .or. &
                      ncw(icanopy) <= 0._r8 .or. ncw(icanopy) >= 1000._r8) CYCLE
               CASE (15)
                  usgs_counts(1) = usgs_counts(1) + 1
                  IF (.not. all(ieee_is_finite([ncd(icanopy), ncw(icanopy), bcw(icanopy), htop(icanopy)]))) CYCLE
                  IF (ncd(icanopy) <= 0._r8 .or. ncd(icanopy) >= 1000._r8 .or. &
                      ncw(icanopy) <= 0._r8 .or. ncw(icanopy) >= 1000._r8 .or. &
                      bcw(icanopy) <= 0._r8 .or. bcw(icanopy) >= 1000._r8 .or. &
                      htop(icanopy) <= 0._r8 .or. htop(icanopy) >= 1000._r8) CYCLE
               CASE (8)
                  usgs_counts(3) = usgs_counts(3) + 1
                  CYCLE
               CASE DEFAULT
                  usgs_counts(4) = usgs_counts(4) + 1
                  CYCLE
               END SELECT
               usgs_counts(2) = usgs_counts(2) + 1
            ENDDO
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, usgs_counts, 4, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
         IF (usgs_counts(2) < usgs_counts(1)) THEN
            IF (p_is_master) write(*,'(A,I0,A,I0,A)') &
               'WARNING: USGS CoLM2024 canopy structure invalid for ', &
               usgs_counts(1)-usgs_counts(2), ' of ', usgs_counts(1), ' forest patches.'
            IF (p_is_master) write(*,'(A)') &
               'WARNING: Downgrading the whole run to CoLM2014 interception; no canopy structure is synthesized.'
            DEF_Interception_scheme = 1
         ELSEIF (usgs_counts(1)+usgs_counts(3) == 0) THEN
            IF (p_is_master) write(*,'(A)') &
               'WARNING: USGS CoLM2024 has no supported forest or shrub patches; using CoLM2014 interception.'
            DEF_Interception_scheme = 1
         ELSE
            IF (p_is_master) write(*,'(A,I0,A,I0,A,I0,A)') &
               'USGS CoLM2024 capacity: forest=', usgs_counts(1), ', shrub=', usgs_counts(3), &
               ', other land patches retaining CoLM2014 capacity=', usgs_counts(4), '.'
         ENDIF
      ENDIF

#endif

#ifdef LULC_IGBP
#ifdef SinglePoint
      allocate (htoplc (numpatch))
      htoplc(:) = SITE_htop
      IF (p_is_worker) THEN
         ncd(:) = SITE_ncd
         ncw(:) = SITE_ncw
         bcw(:) = SITE_bcw
      ENDIF
      IF (DEF_Interception_scheme == 8) THEN
         site_required_count = count(patchtype == 0 .and. &
                                    (patchclass == 1 .or. patchclass == 2 .or. &
                                     patchclass == 3 .or. patchclass == 4 .or. patchclass == 5))
         site_valid_count = count((patchtype == 0 .and. (patchclass == 1 .or. patchclass == 3) .and. &
                                  ieee_is_finite(SITE_ncd) .and. SITE_ncd > 0._r8 .and. SITE_ncd < 1000._r8 .and. &
                                  ieee_is_finite(SITE_ncw) .and. SITE_ncw > 0._r8 .and. SITE_ncw < 1000._r8) .or. &
                                 (patchtype == 0 .and. (patchclass == 2 .or. patchclass == 4) .and. &
                                  ieee_is_finite(SITE_bcw) .and. SITE_bcw > 0._r8 .and. SITE_bcw < 1000._r8 .and. &
                                  ieee_is_finite(SITE_htop) .and. SITE_htop > 0._r8 .and. SITE_htop < 1000._r8) .or. &
                                 (patchtype == 0 .and. patchclass == 5 .and. &
                                  ieee_is_finite(SITE_ncd) .and. SITE_ncd > 0._r8 .and. SITE_ncd < 1000._r8 .and. &
                                  ieee_is_finite(SITE_ncw) .and. SITE_ncw > 0._r8 .and. SITE_ncw < 1000._r8 .and. &
                                  ieee_is_finite(SITE_bcw) .and. SITE_bcw > 0._r8 .and. SITE_bcw < 1000._r8 .and. &
                                  ieee_is_finite(SITE_htop) .and. SITE_htop > 0._r8 .and. SITE_htop < 1000._r8))
         IF (site_required_count > 0 .and. site_valid_count == 0) THEN
            IF (p_is_master) write(*,'(A,I0,A)') &
               'Warning: CoLM2024 requested but SinglePoint has no valid canopy structure for ', &
               site_required_count, ' required tree patches; using CoLM2014 interception.'
            DEF_Interception_scheme = 1
         ELSEIF (site_valid_count < site_required_count) THEN
            IF (p_is_master) write(*,'(A,I0,A,I0,A)') &
               'ERROR: SinglePoint canopy structure is invalid for ', &
               site_required_count - site_valid_count, ' of ', site_required_count, &
               ' required tree patches.'
         CALL CoLM_stop()
         ENDIF
      ENDIF
#else
      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_read_vector (lndname, 'htop_patches', landpatch, htoplc)
#endif

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            m = patchclass(npatch)

            htop(npatch) = htop0(m)
            hbot(npatch) = hbot0(m)

            ! trees or woody savannas
            IF ( m<6 .or. m==8 ) THEN
               ! 01/06/2020, yuan: adjust htop reading
               ! 11/15/2021, yuan: adjust htop setting
               htop(npatch) = max(2., htoplc(npatch))
               hbot(npatch) = htoplc(npatch)*hbot0(m)/htop0(m)
               hbot(npatch) = max(1., hbot(npatch))
            ENDIF

         ENDDO
      ENDIF

#ifndef SinglePoint
      IF (DEF_Interception_scheme == 8 .and. all(struct_present)) THEN
         canopy_counts = 0
         IF (p_is_worker) THEN
            canopy_counts(1) = count(patchtype == 0 .and. &
                                     (patchclass == 1 .or. patchclass == 2 .or. &
                                      patchclass == 3 .or. patchclass == 4 .or. patchclass == 5))
            DO icanopy = 1, numpatch
               IF (patchtype(icanopy) /= 0) CYCLE
               SELECT CASE (patchclass(icanopy))
               CASE (1,3)
                  IF (.not. all(ieee_is_finite([ncd(icanopy), ncw(icanopy)]))) CYCLE
                  IF (ncd(icanopy) <= 0._r8 .or. ncd(icanopy) >= 1000._r8 .or. &
                      ncw(icanopy) <= 0._r8 .or. ncw(icanopy) >= 1000._r8) CYCLE
               CASE (2,4)
                  IF (.not. all(ieee_is_finite([bcw(icanopy), htop(icanopy)]))) CYCLE
                  IF (bcw(icanopy) <= 0._r8 .or. bcw(icanopy) >= 1000._r8 .or. &
                      htop(icanopy) <= 0._r8 .or. htop(icanopy) >= 1000._r8) CYCLE
               CASE (5)
                  IF (.not. all(ieee_is_finite([ncd(icanopy), ncw(icanopy), bcw(icanopy), htop(icanopy)]))) CYCLE
                  IF (ncd(icanopy) <= 0._r8 .or. ncd(icanopy) >= 1000._r8 .or. &
                      ncw(icanopy) <= 0._r8 .or. ncw(icanopy) >= 1000._r8 .or. &
                      bcw(icanopy) <= 0._r8 .or. bcw(icanopy) >= 1000._r8 .or. &
                      htop(icanopy) <= 0._r8 .or. htop(icanopy) >= 1000._r8) CYCLE
               CASE DEFAULT
                  CYCLE
               END SELECT
               canopy_counts(2) = canopy_counts(2) + 1
            ENDDO
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, canopy_counts, 2, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
         IF (canopy_counts(1) > 0 .and. canopy_counts(2) == 0) THEN
            IF (p_is_master) write(*,'(A,I0,A)') &
               'Warning: CoLM2024 requested but surface data have no valid canopy structure for ', &
               canopy_counts(1), ' required tree patches; using CoLM2014 interception.'
            DEF_Interception_scheme = 1
         ELSEIF (canopy_counts(2) < canopy_counts(1)) THEN
            IF (p_is_master) write(*,'(A,I0,A,I0,A)') &
               'WARNING: CoLM2024 surface canopy structure is invalid for ', &
               canopy_counts(1) - canopy_counts(2), ' of ', canopy_counts(1), ' required tree patches.'
            IF (p_is_master) write(*,'(A)') &
               'WARNING: Downgrading the whole run to CoLM2014 interception; no canopy structure is synthesized.'
            DEF_Interception_scheme = 1
         ENDIF
      ENDIF
#endif

      IF (allocated(htoplc))   deallocate ( htoplc )
#endif


#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
#ifdef SinglePoint
      IF (numpft > 0) THEN
         allocate(htoppft(numpft))
         htoppft = pack(SITE_htop_pfts, SITE_pctpfts > 0.)
         site_struct_ok = .false.
         IF (allocated(SITE_ncd_pfts)) THEN
            IF (allocated(SITE_ncw_pfts)) THEN
               IF (allocated(SITE_bcw_pfts)) THEN
                  IF (size(SITE_ncd_pfts) == size(SITE_pctpfts) .and. &
                      size(SITE_ncw_pfts) == size(SITE_pctpfts) .and. &
                      size(SITE_bcw_pfts) == size(SITE_pctpfts)) THEN
                     site_struct_ok = .true.
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (site_struct_ok) THEN
            ncd_p = pack(SITE_ncd_pfts, SITE_pctpfts > 0.)
            ncw_p = pack(SITE_ncw_pfts, SITE_pctpfts > 0.)
            bcw_p = pack(SITE_bcw_pfts, SITE_pctpfts > 0.)
         ELSEIF (allocated(SITE_ncd_pfts) .or. allocated(SITE_ncw_pfts) .or. allocated(SITE_bcw_pfts)) THEN
            IF (p_is_master) write(*,'(A)') &
               'ERROR: SinglePoint PFT canopy-structure arrays must all match SITE_pctpfts size.'
            CALL CoLM_stop()
         ENDIF
         IF (DEF_Interception_scheme == 8) THEN
            IF (site_struct_ok) THEN
               site_required_count = count((SITE_pfttyp >= 1 .and. SITE_pfttyp <= 3) .or. &
                                           (SITE_pfttyp >= 4 .and. SITE_pfttyp <= 8))
               site_valid_count = count((SITE_pfttyp >= 1 .and. SITE_pfttyp <= 3 .and. &
                                         ieee_is_finite(SITE_ncd_pfts) .and. &
                                         SITE_ncd_pfts > 0._r8 .and. SITE_ncd_pfts < 1000._r8 .and. &
                                         ieee_is_finite(SITE_ncw_pfts) .and. &
                                         SITE_ncw_pfts > 0._r8 .and. SITE_ncw_pfts < 1000._r8) .or. &
                                        (SITE_pfttyp >= 4 .and. SITE_pfttyp <= 8 .and. &
                                         ieee_is_finite(SITE_bcw_pfts) .and. &
                                         SITE_bcw_pfts > 0._r8 .and. SITE_bcw_pfts < 1000._r8 .and. &
                                         ieee_is_finite(SITE_htop_pfts) .and. &
                                         SITE_htop_pfts > 0._r8 .and. SITE_htop_pfts < 1000._r8))
               IF (site_required_count > 0 .and. site_valid_count == 0) THEN
                  IF (p_is_master) write(*,'(A,I0,A)') &
                     'Warning: CoLM2024 requested but SinglePoint has no valid PFT canopy structure for ', &
                     site_required_count, ' required tree PFTs; using CoLM2014 interception.'
                  DEF_Interception_scheme = 1
               ELSEIF (site_valid_count < site_required_count) THEN
                  IF (p_is_master) write(*,'(A,I0,A,I0,A)') &
                     'ERROR: SinglePoint PFT canopy structure is invalid for ', &
                     site_required_count - site_valid_count, ' of ', site_required_count, &
                     ' required tree PFTs.'
               CALL CoLM_stop()
               ENDIF
            ELSE
               IF (p_is_master) write(*,'(A)') &
                  'Warning: CoLM2024 requested but SinglePoint lacks valid ncd_pfts/ncw_pfts/bcw_pfts; ' // &
                  'using CoLM2014 interception. Add these fields to SITE_fsitedata to enable CoLM2024.'
               DEF_Interception_scheme = 1
            ENDIF
         ENDIF
      ENDIF
#else
      lndname = trim(landdir)//'/htop_pfts.nc'
      CALL ncio_read_vector (lndname, 'htop_pfts', landpft,   htoppft)

      struct_names = (/ 'ncd_pfts       ', 'ncw_pfts       ', 'bcw_pfts       ' /)
      lndname = trim(cstructdir)//'/ncd_pfts.nc'
      struct_present(1) = ncio_vector_var_present(lndname, struct_names(1), landpft)
      lndname = trim(cstructdir)//'/ncw_pfts.nc'
      struct_present(2) = ncio_vector_var_present(lndname, struct_names(2), landpft)
      lndname = trim(cstructdir)//'/bcw_pfts.nc'
      struct_present(3) = ncio_vector_var_present(lndname, struct_names(3), landpft)
      IF (any(struct_present) .and. .not. all(struct_present)) THEN
         IF (p_is_master) write(*,'(A)') &
            'ERROR: CoLM2024 PFT canopy-structure input is incomplete; ' // &
            'need ncd_pfts, ncw_pfts, and bcw_pfts.'
         CALL CoLM_stop()
      ENDIF
      IF (all(struct_present)) THEN
         lndname = trim(cstructdir)//'/ncd_pfts.nc'
         CALL ncio_read_vector (lndname, 'ncd_pfts', landpft, ncd_p)
         lndname = trim(cstructdir)//'/ncw_pfts.nc'
         CALL ncio_read_vector (lndname, 'ncw_pfts', landpft, ncw_p)
         lndname = trim(cstructdir)//'/bcw_pfts.nc'
         CALL ncio_read_vector (lndname, 'bcw_pfts', landpft, bcw_p)
      ELSEIF (DEF_Interception_scheme == 8) THEN
         IF (p_is_master) write(*,'(A)') &
            'Warning: CoLM2024 requested but PFT canopy-structure inputs are missing; ' // &
            'using CoLM2014 interception. Run mksrfdata with canopy_data to enable CoLM2024.'
         DEF_Interception_scheme = 1
      ENDIF
#endif

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            t = patchtype(npatch)
            m = patchclass(npatch)

            IF (t == 0) THEN
               ps = patch_pft_s(npatch)
               pe = patch_pft_e(npatch)

               DO p = ps, pe
                  n = pftclass(p)

                  htop_p(p) = htop0_p(n)
                  hbot_p(p) = hbot0_p(n)

                  ! for trees
                  ! 01/06/2020, yuan: adjust htop reading
                  ! 11/15/2021, yuan: adjust htop setting
                  IF ( n>0 .and. n<9 ) THEN
                     htop_p(p) = max(2., htoppft(p))
                     hbot_p(p) = htoppft(p)*hbot0_p(n)/htop0_p(n)
                     hbot_p(p) = max(1., hbot_p(p))
                  ENDIF
               ENDDO

               htop(npatch) = sum(htop_p(ps:pe)*pftfrac(ps:pe))
               hbot(npatch) = sum(hbot_p(ps:pe)*pftfrac(ps:pe))

            ELSE
               htop(npatch) = htop0(m)
               hbot(npatch) = hbot0(m)
            ENDIF

         ENDDO
      ENDIF

#ifndef SinglePoint
      IF (DEF_Interception_scheme == 8 .and. all(struct_present)) THEN
         canopy_counts = 0
         IF (p_is_worker) THEN
            canopy_counts(1) = count((pftclass >= 1 .and. pftclass <= 3) .or. &
                                     (pftclass >= 4 .and. pftclass <= 8))
            DO icanopy = 1, numpft
               SELECT CASE (pftclass(icanopy))
               CASE (1:3)
                  IF (.not. all(ieee_is_finite([ncd_p(icanopy), ncw_p(icanopy)]))) CYCLE
                  IF (ncd_p(icanopy) <= 0._r8 .or. ncd_p(icanopy) >= 1000._r8 .or. &
                      ncw_p(icanopy) <= 0._r8 .or. ncw_p(icanopy) >= 1000._r8) CYCLE
               CASE (4:8)
                  IF (.not. all(ieee_is_finite([bcw_p(icanopy), htop_p(icanopy)]))) CYCLE
                  IF (bcw_p(icanopy) <= 0._r8 .or. bcw_p(icanopy) >= 1000._r8 .or. &
                      htop_p(icanopy) <= 0._r8 .or. htop_p(icanopy) >= 1000._r8) CYCLE
               CASE DEFAULT
                  CYCLE
               END SELECT
               canopy_counts(2) = canopy_counts(2) + 1
            ENDDO
         ENDIF
#ifdef USEMPI
         CALL mpi_allreduce(MPI_IN_PLACE, canopy_counts, 2, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
         IF (canopy_counts(1) > 0 .and. canopy_counts(2) == 0) THEN
            IF (p_is_master) write(*,'(A,I0,A)') &
               'Warning: CoLM2024 requested but surface data have no valid PFT canopy structure for ', &
               canopy_counts(1), ' required tree PFTs; using CoLM2014 interception.'
            DEF_Interception_scheme = 1
         ELSEIF (canopy_counts(2) < canopy_counts(1)) THEN
            IF (p_is_master) write(*,'(A,I0,A,I0,A)') &
               'ERROR: CoLM2024 surface PFT canopy structure is invalid for ', &
               canopy_counts(1) - canopy_counts(2), ' of ', canopy_counts(1), ' required tree PFTs.'
            CALL CoLM_stop()
         ENDIF
      ENDIF
#endif

      IF (allocated(htoppft)) deallocate(htoppft)
#endif

   END SUBROUTINE HTOP_readin

END MODULE MOD_HtopReadin
