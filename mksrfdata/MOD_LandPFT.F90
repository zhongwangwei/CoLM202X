#include <define.h>

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

MODULE MOD_LandPFT

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Build pixelset "landpft" (Plant Function Type).
!
!    In CoLM, the global/regional area is divided into a hierarchical structure:
!    1. If GRIDBASED or UNSTRUCTURED is defined, it is
!       ELEMENT >>> PATCH
!    2. If CATCHMENT is defined, it is
!       ELEMENT >>> HRU >>> PATCH
!    If Plant Function Type classification is used, PATCH is further divided into PFT.
!    If Plant Community classification is used,     PATCH is further divided into PC.
!
!    "landpft" refers to pixelset PFT.
!
!  Created by Shupeng Zhang, May 2023
!    porting codes from Hua Yuan's OpenMP version to MPI parallel version.
!-----------------------------------------------------------------------

   USE MOD_Namelist
   USE MOD_Pixelset
   USE MOD_Const_LC
   USE MOD_Vars_Global
   ! SITE_wetland_class arrives through the blanket USE MOD_Namelist above.
   IMPLICIT NONE

   ! ---- Instance ----
   integer :: numpft
   type(pixelset_type) :: landpft

   integer , allocatable :: pft2patch   (:)  !patch index of a PFT
   integer , allocatable :: patch_pft_s (:)  !start PFT index of a patch
   integer , allocatable :: patch_pft_e (:)  !end PFT index of a patch

   ! ---- PUBLIC routines ----
   PUBLIC :: landpft_build

CONTAINS

   ! -------------------------------
   SUBROUTINE landpft_build (lc_year)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Namelist
   USE MOD_5x5DataReadin
   USE MOD_LandPatch
   USE MOD_Land2mWMO
   USE MOD_AggregationRequestData
   USE MOD_Const_LC
#ifdef CROP
   USE MOD_LandCrop
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year
   ! Local Variables
   character(len=256) :: dir_5x5, suffix, cyear
   type (block_data_real8_3d) :: pctpft
   real(r8), allocatable :: pctpft_patch(:,:), pctpft_one(:,:)
   real(r8), allocatable :: area_one  (:)
   logical,  allocatable :: patchmask (:)
#if (defined LULC_IGBP_WFT) && (!defined SinglePoint)
   integer, allocatable :: wft_class_patch(:)
#endif
   integer  :: ipatch, ipft, npatch, npft, npft_glb
   integer  :: wmo_src, ipft_grass
   real(r8) :: sumarea, maxgrass

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land plant function type tiles:'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#if (defined LULC_IGBP_WFT) && (!defined SinglePoint)
      ! Type the wetland tiles before they are built, the way landcrop_build
      ! reads global_CFT_surface_data.nc before typing crop patches.
      CALL wetland_class_readin (wft_class_patch)
#endif

      landpft%has_shared = .true.

      IF (p_is_io) THEN

         CALL allocate_block_data (grid_patch, pctpft, N_PFT_modis, lb1 = 0)
         CALL flush_block_data (pctpft, 1.0)

         dir_5x5 = trim(DEF_dir_rawdata) // '/plant_15s'
         ! add parameter input for time year
         write(cyear,'(i4.4)') lc_year
         suffix  = 'MOD'//trim(cyear)
         CALL read_5x5_data_pft (dir_5x5, suffix, grid_patch, 'PCT_PFT', pctpft)

#ifdef USEMPI
         CALL aggregation_data_daemon (grid_patch, data_r8_3d_in1 = pctpft, n1_r8_3d_in1 = N_PFT_modis)
#endif
      ENDIF


      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN
            allocate (pctpft_patch (0:N_PFT-1,numpatch))
            allocate (patchmask (numpatch))

            pctpft_patch(:,:) = 0
            patchmask (:) = .true.
         ENDIF

         DO ipatch = 1, numpatch

            IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN

               wmo_src  = wmo_source(landpatch%ielm(ipatch))
               maxgrass = maxval(pctpft_patch(12:14,wmo_src))

               IF (maxgrass > 0) THEN
                  ipft_grass = maxloc(pctpft_patch(12:14,wmo_src), dim=1) + 11
                  pctpft_patch(:,ipatch) = 0
                  pctpft_patch(ipft_grass,ipatch) = 1.
               ELSE
                  pctpft_patch(0,ipatch) = 1.
               ENDIF

               CYCLE
            ENDIF
#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
               CALL aggregation_request_data (landpatch, ipatch, grid_patch, zip = .false., area = area_one, &
                  data_r8_3d_in1 = pctpft, data_r8_3d_out1 = pctpft_one, n1_r8_3d_in1 = N_PFT_modis, lb1_r8_3d_in1 = 0)

               sumarea = sum(area_one * sum(pctpft_one(0:N_PFT-1,:),dim=1))

               ! in case of no PFT data, set to 100% bare when patchtype=0,
               ! be consistent with Aggregation_PercentagesPFT.F90.
               IF (sumarea <= 0.0) THEN
                  pctpft_patch(0,ipatch) = 1.
               ELSE
                  DO ipft = 0, N_PFT-1
                     pctpft_patch(ipft,ipatch) = sum(pctpft_one(ipft,:) * area_one) / sumarea
                  ENDDO
               ENDIF

            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif

         IF (numpatch > 0) THEN
            npatch = count(patchmask)
            numpft = count(pctpft_patch > 0.)
#ifdef CROP
            numpft = numpft + count(landpatch%settyp == CROPLAND)
#endif
#ifdef LULC_IGBP_WFT
            ! one wetland functional type tile per permanent-wetland patch
            numpft = numpft + count(landpatch%settyp == WETLAND)
#endif
            IF (npatch > 0) THEN
               allocate (patch_pft_s (npatch))
               allocate (patch_pft_e (npatch))
            ENDIF
         ELSE
            numpft = 0
         ENDIF

         IF (numpft > 0) THEN

            allocate (pft2patch      (numpft))

            allocate (landpft%eindex (numpft))
            allocate (landpft%settyp (numpft))
            allocate (landpft%ipxstt (numpft))
            allocate (landpft%ipxend (numpft))
            allocate (landpft%ielm   (numpft))

            allocate (landpft%pctshared (numpft))
            landpft%pctshared(:) = 1.

            npft = 0
            npatch = 0
            DO ipatch = 1, numpatch
               IF (patchmask(ipatch)) THEN
                  npatch = npatch + 1

#ifndef CROP
                  IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
                  IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
                     patch_pft_s(npatch) = npft + 1
                     patch_pft_e(npatch) = npft + count(pctpft_patch(:,ipatch) > 0)

                     DO ipft = 0, N_PFT-1
                        IF (pctpft_patch(ipft,ipatch) > 0) THEN
                           npft = npft + 1

                           landpft%ielm  (npft) = landpatch%ielm  (ipatch)
                           landpft%eindex(npft) = landpatch%eindex(ipatch)
                           landpft%ipxstt(npft) = landpatch%ipxstt(ipatch)
                           landpft%ipxend(npft) = landpatch%ipxend(ipatch)
                           landpft%settyp(npft) = ipft

                           landpft%pctshared(npft) = pctpft_patch(ipft,ipatch)

                           pft2patch(npft) = npatch
                        ENDIF
                     ENDDO
#ifdef CROP
                  ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
                     npft = npft + 1
                     patch_pft_s(npatch) = npft
                     patch_pft_e(npatch) = npft

                     landpft%ielm  (npft) = landpatch%ielm  (ipatch)
                     landpft%eindex(npft) = landpatch%eindex(ipatch)
                     landpft%ipxstt(npft) = landpatch%ipxstt(ipatch)
                     landpft%ipxend(npft) = landpatch%ipxend(ipatch)
                     landpft%settyp(npft) = cropclass(ipatch) + N_PFT - 1

                     landpft%pctshared(npft) = landpatch%pctshared(ipatch)

                     pft2patch(npft) = npatch
#endif
#ifdef LULC_IGBP_WFT
                  ELSEIF (landpatch%settyp(ipatch) == WETLAND) THEN
                     ! Permanent wetland: one WFT tile covering the whole patch,
                     ! typed by the wetland class (79..86 = npwetlandmin+class-1).
                     ! Without it patch_pft_s/e stay at -1, the CN driver never
                     ! runs, and the soil carbon pools decay with no litter input.
                     npft = npft + 1
                     patch_pft_s(npatch) = npft
                     patch_pft_e(npatch) = npft

                     landpft%ielm  (npft) = landpatch%ielm  (ipatch)
                     landpft%eindex(npft) = landpatch%eindex(ipatch)
                     landpft%ipxstt(npft) = landpatch%ipxstt(ipatch)
                     landpft%ipxend(npft) = landpatch%ipxend(ipatch)
#ifdef SinglePoint
                     landpft%settyp(npft) = npwetlandmin + SITE_wetland_class - 1
#else
                     landpft%settyp(npft) = npwetlandmin + wft_class_patch(ipatch) - 1
#endif

                     landpft%pctshared(npft) = 1.

                     pft2patch(npft) = npatch
#endif
                  ELSE
                     patch_pft_s(npatch) = -1
                     patch_pft_e(npatch) = -1
                  ENDIF
               ENDIF
            ENDDO

         ENDIF

      ENDIF

      CALL landpatch%pset_pack(patchmask, numpatch)

      landpft%nset = numpft
      CALL landpft%set_vecgs

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpft, npft_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npft_glb, ' plant function type tiles.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpft, ' plant function type tiles.'
#endif

      IF (allocated(pctpft_patch)) deallocate (pctpft_patch)
      IF (allocated(pctpft_one  )) deallocate (pctpft_one  )
      IF (allocated(area_one    )) deallocate (area_one    )
      IF (allocated(patchmask   )) deallocate (patchmask   )

   END SUBROUTINE landpft_build

#if (defined LULC_IGBP_WFT) && (!defined SinglePoint)
   SUBROUTINE wetland_class_readin (wft_class_patch)

   ! Type the wetland tiles from the eight-class rawdata map, the way
   ! landcrop_build types crop patches from global_CFT_surface_data.nc:
   ! read at build time, one class per wetland patch. The vote over the
   ! patch's 5-arcmin cells is weighted by overlap area times WET_FRAC
   ! (the GLWD wetland share of the cell), so it follows where the
   ! wetland actually is, not the dry background of the cell.

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_dir_rawdata
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_Vars_Global, only: WETLAND, N_WFT
   USE MOD_Tracer_Reactive_Methane_PHMapping, only: methane_ph_mapping_type, &
      build_methane_ph_areal_mapping
   USE netcdf
   IMPLICIT NONE

   integer, allocatable, intent(out) :: wft_class_patch(:)

   character(len=256) :: rawfile
   logical  :: raw_exists
   integer  :: ncid, vid, did, ierr, nlat, nlon
   integer  :: iset, ipart, iproc, iloc, jlat, jlon, ic
   integer  :: nfallback, nweighted_empty, tmp
   integer  :: count_loc(N_WFT), count_glb(N_WFT)
   real(r8) :: votes(N_WFT), votes_area(N_WFT), area, w
   real(r8) :: south_edge, north_edge
   real(r8),   allocatable :: lat_g(:), lon_g(:)
   integer(1), allocatable :: class_g(:,:)
   real(4),    allocatable :: wetf_g(:,:)
   type(grid_type) :: grid_wft
   type(methane_ph_mapping_type) :: map_wft

      rawfile = trim(DEF_dir_rawdata) // '/global_WFT_surface_data.nc'
      raw_exists = .false.
      IF (p_is_master) inquire (file=trim(rawfile), exist=raw_exists)
#ifdef USEMPI
      CALL mpi_bcast (raw_exists, 1, MPI_LOGICAL, p_address_master, p_comm_glb, p_err)
#endif
      IF (.not. raw_exists) THEN
         CALL CoLM_stop (' ***** ERROR: wetland class rawdata not found: '//trim(rawfile))
      ENDIF

      nlat = 0; nlon = 0
      IF (p_is_master) THEN
         ierr = nf90_open(trim(rawfile), NF90_NOWRITE, ncid)
         IF (ierr /= NF90_NOERR) CALL CoLM_stop &
            (' ***** ERROR: cannot open '//trim(rawfile))
         ierr = nf90_inq_dimid(ncid, 'lat', did)
         IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension(ncid, did, len=nlat)
         IF (ierr == NF90_NOERR) ierr = nf90_inq_dimid(ncid, 'lon', did)
         IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension(ncid, did, len=nlon)
         IF (ierr /= NF90_NOERR .or. nlat < 2 .or. nlon < 2) CALL CoLM_stop &
            (' ***** ERROR: invalid lat/lon in wetland class rawdata')
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (nlat, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (nlon, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      allocate (lat_g(nlat), lon_g(nlon))
      allocate (class_g(nlon,nlat), wetf_g(nlon,nlat))
      IF (p_is_master) THEN
         ierr = nf90_inq_varid(ncid, 'lat', vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lat_g)
         IF (ierr == NF90_NOERR) ierr = nf90_inq_varid(ncid, 'lon', vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lon_g)
         IF (ierr == NF90_NOERR) ierr = nf90_inq_varid(ncid, 'WFT_CLASS', vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, class_g)
         IF (ierr == NF90_NOERR) ierr = nf90_inq_varid(ncid, 'WET_FRAC', vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, wetf_g)
         IF (ierr /= NF90_NOERR) CALL CoLM_stop &
            (' ***** ERROR: cannot read WFT_CLASS/WET_FRAC: '//trim(nf90_strerror(ierr)))
         ierr = nf90_close(ncid)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (lat_g, nlat, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (lon_g, nlon, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (class_g, nlon*nlat, MPI_INTEGER1, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (wetf_g,  nlon*nlat, MPI_REAL4,    p_address_master, p_comm_glb, p_err)
#endif

      IF (lat_g(1) < lat_g(nlat)) THEN
         south_edge = lat_g(1)    - 0.5_r8 * (lat_g(2) - lat_g(1))
         north_edge = lat_g(nlat) + 0.5_r8 * (lat_g(nlat) - lat_g(nlat-1))
      ELSE
         south_edge = lat_g(nlat) - 0.5_r8 * (lat_g(nlat-1) - lat_g(nlat))
         north_edge = lat_g(1)    + 0.5_r8 * (lat_g(1) - lat_g(2))
      ENDIF
      south_edge = max(-90._r8, south_edge)
      north_edge = min( 90._r8, north_edge)
      CALL grid_wft%define_by_center (lat_g, lon_g, south=south_edge, north=north_edge)
      CALL build_methane_ph_areal_mapping (map_wft, grid_wft, landpatch)

      nfallback = 0
      nweighted_empty = 0
      count_loc(:) = 0
      IF (p_is_worker) THEN
         allocate (wft_class_patch(max(numpatch,1)))
         wft_class_patch(:) = 0
         DO iset = 1, numpatch
            IF (landpatch%settyp(iset) /= WETLAND) CYCLE
            votes(:) = 0._r8
            votes_area(:) = 0._r8
            DO ipart = 1, map_wft%npart(iset)
               iproc = map_wft%address(iset)%val(1,ipart)
               iloc  = map_wft%address(iset)%val(2,ipart)
               jlat  = map_wft%glist(iproc)%ilat(iloc)
               jlon  = map_wft%glist(iproc)%ilon(iloc)
               ic    = int(class_g(jlon,jlat))
               area  = map_wft%areapart(iset)%val(ipart)
               IF (ic < 1 .or. ic > N_WFT .or. area <= 0._r8) CYCLE
               w = area * max(real(wetf_g(jlon,jlat), r8), 0._r8)
               votes(ic) = votes(ic) + w
               votes_area(ic) = votes_area(ic) + area
            ENDDO
            IF (any(votes > 0._r8)) THEN
               wft_class_patch(iset) = maxloc(votes, dim=1)
            ELSEIF (any(votes_area > 0._r8)) THEN
               ! no GLWD wetland inside: fall back to plain-area majority
               nweighted_empty = nweighted_empty + 1
               wft_class_patch(iset) = maxloc(votes_area, dim=1)
            ELSE
               nfallback = nfallback + 1
               wft_class_patch(iset) = 5    ! herbaceous mineral fallback
            ENDIF
            ic = wft_class_patch(iset)
            count_loc(ic) = count_loc(ic) + 1
         ENDDO
      ELSE
         allocate (wft_class_patch(1))
         wft_class_patch(:) = 0
      ENDIF

      count_glb = count_loc
#ifdef USEMPI
      CALL mpi_allreduce (count_loc, count_glb, N_WFT, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      tmp = nfallback
      CALL mpi_allreduce (tmp, nfallback, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
      tmp = nweighted_empty
      CALL mpi_allreduce (tmp, nweighted_empty, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(A)') 'Wetland functional type classes on wetland patches:'
         DO ic = 1, N_WFT
            write(*,'(A,I2,A,I10)') '  class ', ic, ' patches: ', count_glb(ic)
         ENDDO
         write(*,'(A,I10)') '  area-majority fallbacks (no GLWD wetland in cell): ', &
            nweighted_empty
         write(*,'(A,I10)') '  hard fallbacks to class 5 (no map data at all):    ', &
            nfallback
      ENDIF

      deallocate (lat_g, lon_g, class_g, wetf_g)

   END SUBROUTINE wetland_class_readin
#endif

   ! ----------------------
   SUBROUTINE map_patch_to_pft

   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_Const_LC
   IMPLICIT NONE

   integer :: ipatch, ipft

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numpft <= 0)) RETURN

         IF (allocated(patch_pft_s)) deallocate(patch_pft_s)
         IF (allocated(patch_pft_e)) deallocate(patch_pft_e)
         IF (allocated(pft2patch  )) deallocate(pft2patch  )

         allocate (patch_pft_s (numpatch))
         allocate (patch_pft_e (numpatch))
         allocate (pft2patch   (numpft  ))

         ipft = 1
         DO ipatch = 1, numpatch
#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif

               patch_pft_s(ipatch) = ipft

               DO WHILE (ipft <= numpft)
                  IF ((landpft%eindex(ipft) == landpatch%eindex(ipatch))  &
                     .and. (landpft%ipxstt(ipft) == landpatch%ipxstt(ipatch))  &
                     .and. (landpft%settyp(ipft) < N_PFT)) THEN
                     pft2patch  (ipft  ) = ipatch
                     patch_pft_e(ipatch) = ipft
                     ipft = ipft + 1
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
#ifdef CROP
            ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
               patch_pft_s(ipatch) = ipft
               patch_pft_e(ipatch) = ipft
               pft2patch  (ipft  ) = ipatch
               ipft = ipft + 1
#endif
#ifdef LULC_IGBP_WFT
            ELSEIF (landpatch%settyp(ipatch) == WETLAND) THEN
               patch_pft_s(ipatch) = ipft
               patch_pft_e(ipatch) = ipft
               pft2patch  (ipft  ) = ipatch
               ipft = ipft + 1
#endif
            ELSE
               patch_pft_s(ipatch) = -1
               patch_pft_e(ipatch) = -1
            ENDIF

         ENDDO

      ENDIF

   END SUBROUTINE map_patch_to_pft

END MODULE MOD_LandPFT
#endif
