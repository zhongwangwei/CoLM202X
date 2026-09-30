#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_WetlandVeg

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Vegetation of the permanent-wetland tile from its GLWD v2 make-up
!   (C-13, paper V2), replacing the five-zone climate proxy when
!   DEF_METHANE%wetland_veg_glwd is true.
!
!   Per patch, from the grid cell nearest to the patch centre:
!   - wetveg_forest: forested share of the cell's wetland-tile classes,
!     (16, 18, 22, 24, 26) / (16-19, 22-27);
!   - wetveg_laicap: peak LAI of the non-forested share, the area mean of
!     wetland_lai_open_peat over the open-peatland classes (23, 25) and
!     wetland_lai_marsh over the marsh classes (17, 19, 27).
!   A cell without any tile class gives forest share 0 and no LAI cap.
!
! INPUT FILE: DEF_METHANE%wetland_veg_file, a regular lat/lon grid with
!   variables lat, lon, forested_share and area_class_NN (km2), as written
!   by v2/scripts/glwd_class_2deg.py.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: PI
   USE, INTRINSIC :: IEEE_ARITHMETIC, only: ieee_is_finite

   IMPLICIT NONE
   SAVE
   PRIVATE

   real(r8), allocatable, public :: wetveg_forest(:)   ! forested share [-]
   real(r8), allocatable, public :: wetveg_laicap(:)   ! LAI cap of the non-forested share [m2/m2]
   logical,  public :: wetveg_active = .false.

   PUBLIC :: read_methane_wetveg
   PUBLIC :: deallocate_methane_wetveg
   PUBLIC :: wetveg_cap_lai
   PUBLIC :: wetveg_bg_frac
   PUBLIC :: set_wetland_veg_glwd

CONTAINS

   SUBROUTINE deallocate_methane_wetveg ()
      IF (allocated(wetveg_forest)) deallocate(wetveg_forest)
      IF (allocated(wetveg_laicap)) deallocate(wetveg_laicap)
      wetveg_active = .false.
   END SUBROUTINE deallocate_methane_wetveg

   SUBROUTINE read_methane_wetveg (file_veg, lai_open_peat, lai_marsh, patchlatr_in, patchlonr_in, numpatch)
#ifdef USEMPI
      USE MPI
#endif
      USE netcdf
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE

      character(len=*), intent(in) :: file_veg
      real(r8), intent(in) :: lai_open_peat, lai_marsh
      real(r8), intent(in) :: patchlatr_in(:), patchlonr_in(:)   ! radians
      integer,  intent(in) :: numpatch

      integer, parameter :: nopen = 2, nmarsh = 3
      integer, parameter :: open_classes(nopen) = (/23, 25/), marsh_classes(nmarsh) = (/17, 19, 27/)
      integer :: ncid, vid, ierr, nlat, nlon, ilat, ilon, ip, k, dims(2), bad
      real(r8), allocatable :: lat_g(:), lon_g(:), forest_g(:,:), cap_g(:,:), a(:,:), aopen(:,:), amarsh(:,:)
      real(r8) :: lat_deg, lon_deg, d, dmin, dlon
      character(len=16) :: vname

      CALL deallocate_methane_wetveg ()
      bad = 0
      dims = 0

      IF (p_is_master) THEN
         ierr = nf90_open(trim(file_veg), NF90_NOWRITE, ncid)
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A,A,A)') ' ERROR: wetland vegetation file ', trim(file_veg), ': ', trim(nf90_strerror(ierr))
            bad = 1
         ELSE
            ierr = nf90_inq_dimid(ncid, 'lat', vid)
            IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension(ncid, vid, len = dims(1))
            IF (ierr == NF90_NOERR) ierr = nf90_inq_dimid(ncid, 'lon', vid)
            IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension(ncid, vid, len = dims(2))
            IF (ierr /= NF90_NOERR .or. dims(1) <= 0 .or. dims(2) <= 0) bad = 1
         ENDIF
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (bad,  1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (dims, 2, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      IF (bad /= 0) CALL CoLM_stop (' ***** ERROR: cannot read DEF_METHANE%wetland_veg_file.')
      nlat = dims(1); nlon = dims(2)
      allocate (lat_g(nlat), lon_g(nlon), forest_g(nlat,nlon), cap_g(nlat,nlon))

      IF (p_is_master) THEN
         allocate (a(nlon,nlat), aopen(nlat,nlon), amarsh(nlat,nlon))
         ! Every variable must be read; a failed read would otherwise leave the
         ! previous class's values in a and change LAI and CH4 transport silently.
         vname = 'lat'
         ierr = nf90_inq_varid(ncid, trim(vname), vid)
         IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lat_g)
         IF (ierr == NF90_NOERR) THEN
            vname = 'lon'
            ierr = nf90_inq_varid(ncid, trim(vname), vid)
            IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lon_g)
         ENDIF
         IF (ierr == NF90_NOERR) THEN
            vname = 'forested_share'
            ierr = nf90_inq_varid(ncid, trim(vname), vid)
            IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, a)
            IF (ierr == NF90_NOERR) forest_g = transpose(a)
         ENDIF
         aopen = 0._r8; amarsh = 0._r8
         DO k = 1, nopen
            IF (ierr /= NF90_NOERR) EXIT
            write(vname,'(A,I2.2)') 'area_class_', open_classes(k)
            ierr = nf90_inq_varid(ncid, trim(vname), vid)
            IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, a)
            IF (ierr == NF90_NOERR) aopen = aopen + max(transpose(a), 0._r8)
         ENDDO
         DO k = 1, nmarsh
            IF (ierr /= NF90_NOERR) EXIT
            write(vname,'(A,I2.2)') 'area_class_', marsh_classes(k)
            ierr = nf90_inq_varid(ncid, trim(vname), vid)
            IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, a)
            IF (ierr == NF90_NOERR) amarsh = amarsh + max(transpose(a), 0._r8)
         ENDDO
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A,A,A,A,A)') ' ERROR: wetland vegetation file ', trim(file_veg), &
               ', variable ', trim(vname), ': ', trim(nf90_strerror(ierr))
            bad = 1
         ENDIF
         ierr = nf90_close(ncid)
         ! netCDF fill for cells without tile classes: no forest, no cap
         WHERE (.not. ieee_is_finite(forest_g) .or. forest_g < 0._r8 .or. forest_g > 1._r8) forest_g = 0._r8
         WHERE (aopen + amarsh > 0._r8)
            cap_g = (aopen * lai_open_peat + amarsh * lai_marsh) / (aopen + amarsh)
         ELSEWHERE
            cap_g = huge(1._r8)
         END WHERE
         deallocate (a, aopen, amarsh)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (bad,      1,         MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif
      IF (bad /= 0) CALL CoLM_stop (' ***** ERROR: cannot read DEF_METHANE%wetland_veg_file.')
#ifdef USEMPI
      CALL mpi_bcast (lat_g,    nlat,      MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (lon_g,    nlon,      MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (forest_g, nlat*nlon, MPI_REAL8, p_address_master, p_comm_glb, p_err)
      CALL mpi_bcast (cap_g,    nlat*nlon, MPI_REAL8, p_address_master, p_comm_glb, p_err)
#endif

      allocate (wetveg_forest(max(numpatch,0)), wetveg_laicap(max(numpatch,0)))
      DO ip = 1, numpatch
         lat_deg = patchlatr_in(ip) * 180._r8 / PI
         lon_deg = patchlonr_in(ip) * 180._r8 / PI
         dmin = huge(1._r8); ilat = 1
         DO k = 1, nlat
            d = abs(lat_g(k) - lat_deg)
            IF (d < dmin) THEN; dmin = d; ilat = k; ENDIF
         ENDDO
         dmin = huge(1._r8); ilon = 1
         DO k = 1, nlon
            dlon = abs(modulo(lon_g(k) - lon_deg + 180._r8, 360._r8) - 180._r8)
            IF (dlon < dmin) THEN; dmin = dlon; ilon = k; ENDIF
         ENDDO
         wetveg_forest(ip) = forest_g(ilat, ilon)
         wetveg_laicap(ip) = cap_g(ilat, ilon)
      ENDDO
      deallocate (lat_g, lon_g, forest_g, cap_g)
      ! A tower sits in one wetland, not in the cell's mix: single-point runs
      ! give its own make-up through these two keys (< 0: keep the file).
      IF (DEF_METHANE%wetland_forest_share_site >= 0._r8) &
         wetveg_forest(:) = min(DEF_METHANE%wetland_forest_share_site, 1._r8)
      IF (DEF_METHANE%wetland_lai_cap_site > 0._r8) wetveg_laicap(:) = DEF_METHANE%wetland_lai_cap_site
      wetveg_active = .true.
      IF (p_is_master) write(*,'(A,A)') ' C-13 wetland vegetation read from ', trim(file_veg)

   END SUBROUTINE read_methane_wetveg

   SUBROUTINE wetveg_cap_lai ()
      ! Called after each LAI read: the non-forested share of a wetland patch
      ! keeps the remote-sensing LAI only up to its measured peak.
      USE MOD_Vars_TimeVariables, only: tlai
      USE MOD_Vars_TimeInvariants, only: patchtype
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      integer :: i
      IF (.not. (DEF_METHANE%wetland_veg_glwd .and. wetveg_active)) RETURN
      IF (.not. p_is_worker) RETURN
      IF (.not. allocated(tlai) .or. .not. allocated(patchtype)) RETURN
      DO i = 1, min(size(tlai), size(wetveg_forest))
         IF (patchtype(i) /= 2) CYCLE
         tlai(i) = wetveg_forest(i) * tlai(i) + (1._r8 - wetveg_forest(i)) * min(tlai(i), wetveg_laicap(i))
      ENDDO
   END SUBROUTINE wetveg_cap_lai

   real(r8) FUNCTION wetveg_bg_frac (ipatch)
      ! Belowground share of the wetland plant input, weighted by the
      ! forested share (C-13); wetland_bg_frac alone without C-13.
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      integer, intent(in) :: ipatch
      wetveg_bg_frac = DEF_METHANE%wetland_bg_frac
      IF (.not. (DEF_METHANE%wetland_veg_glwd .and. wetveg_active)) RETURN
      IF (ipatch < 1 .or. ipatch > size(wetveg_forest)) RETURN
      wetveg_bg_frac = wetveg_forest(ipatch) * DEF_METHANE%wetland_bg_frac_forest &
         + (1._r8 - wetveg_forest(ipatch)) * DEF_METHANE%wetland_bg_frac
   END FUNCTION wetveg_bg_frac

   SUBROUTINE set_wetland_veg_glwd (ipatch, lai_in, lai_out, annsum_npp_out, agnpp_out, bgnpp_out, rootfr_out)
      ! Replaces get_wetland_veg_proxy (C-13): NPP from the tile's own
      ! assimilation times wetland_npp_frac (as the C-12 input), split by the
      ! forest-weighted belowground share; the land class's root profile (as
      ! the C-12 input); CLM4Me grass aerenchyma defaults, the forested share
      ! at nongrassporosratio of the grass porosity (Riley et al. 2011).
      ! annsum_npp_out is 0 so the aerenchyma uses the annual means of agnpp
      ! and bgnpp accumulated by methane_annualupdate.
      USE MOD_Vars_Global, only: nl_soil
      USE MOD_Vars_1DFluxes, only: assim
      USE MOD_Vars_TimeInvariants, only: patchclass
      USE MOD_Const_LC, only: rootfr
      USE MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
      USE MOD_Tracer_Reactive_Methane_VegOverride, only: wetland_aere_poros, wetland_aere_radius, &
         wetland_aere_tillerC, wetland_aere_scale, wetland_aere_active
      integer,  intent(in)  :: ipatch
      real(r8), intent(in)  :: lai_in
      real(r8), intent(out) :: lai_out, annsum_npp_out, agnpp_out, bgnpp_out
      real(r8), intent(out) :: rootfr_out(1:nl_soil)
      real(r8) :: npp, bg, f, s

      f = 0._r8
      IF (ipatch >= 1 .and. ipatch <= size(wetveg_forest)) f = wetveg_forest(ipatch)
      lai_out = max(lai_in, 0._r8)
      npp = 0._r8
      IF (allocated(assim)) THEN
         IF (ieee_is_finite(assim(ipatch)) .and. assim(ipatch) > 0._r8 .and. assim(ipatch) < 1.e-2_r8) &
            npp = assim(ipatch) * 12.011_r8 * DEF_METHANE%wetland_npp_frac       ! [gC m-2 s-1]
      ENDIF
      bg = wetveg_bg_frac(ipatch)
      agnpp_out = (1._r8 - bg) * npp
      bgnpp_out = bg * npp
      annsum_npp_out = 0._r8

      rootfr_out(:) = max(rootfr(1:nl_soil, patchclass(ipatch)), 0._r8)
      s = sum(rootfr_out)
      IF (s > 0._r8) THEN
         rootfr_out(:) = rootfr_out(:) / s
      ELSE
         rootfr_out(:) = 0._r8
         rootfr_out(1) = 1._r8
      ENDIF

      IF (allocated(wetland_aere_active) .and. ipatch >= 1 .and. ipatch <= size(wetland_aere_active)) THEN
         wetland_aere_poros  (ipatch) = DEF_METHANE%poros_tiller * &
            (1._r8 - f + f * DEF_METHANE%nongrassporosratio)
         wetland_aere_radius (ipatch) = DEF_METHANE%aere_radius
         wetland_aere_tillerC(ipatch) = DEF_METHANE%tiller_C
         wetland_aere_scale  (ipatch) = DEF_METHANE%scale_factor_aere
         wetland_aere_active (ipatch) = .true.
      ENDIF
   END SUBROUTINE set_wetland_veg_glwd

END MODULE MOD_Tracer_Reactive_Methane_WetlandVeg
#endif
