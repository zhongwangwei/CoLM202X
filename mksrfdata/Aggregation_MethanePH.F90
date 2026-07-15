#include <define.h>

SUBROUTINE Aggregation_MethanePH (dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
! !DESCRIPTION:
!  Build a methane patch-level soil pH surface file from the native
!  PHH2O1.nc raw dataset.  This is intentionally done during mksrfdata
!  rather than at runtime: PHH2O1.nc is a 30 arc-sec global file
!  (43200 x 16800 x 4, ~2.8 GB), and broadcasting a gridded pH field
!  during colm.x initialization can deadlock or exhaust memory.
!
!  Output:
!    <dir_model_landdata>/soil/<lc_year>/methane_ph_patches.nc
!       variable methane_ph_patches(patch), pH units
!
!  Method:
!    Each worker computes a patch centroid from its local mesh pixels and
!    samples the native PHH2O1 grid around that centroid (3x3 native pixels,
!    top four depth layers).  Missing pixels are skipped; unmatched patches
!    use the methane neutral fallback pH=6.2.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global, only: PI
   USE MOD_Namelist, only: DEF_Srfdata_CompressLevel
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_Land2mWMO, only: wmo_patch, wmo_source
   USE MOD_Mesh
   USE MOD_Pixel
   USE MOD_NetCDFVector
   USE MOD_Utils, only: areaquad
   USE netcdf
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   real(r8), parameter :: methane_ph_fallback = 6.2_r8
   integer(1), parameter :: missing_byte = -100_1
   integer, parameter :: native_search_radius = 1

   character(len=256) :: landdir, lndname, cyear, rawfile
   logical :: raw_exists, ok_center
   integer :: ncid, vid, ierr
   integer :: nlat, nlon, ndepth
   integer :: ipatch, wmo_src
   real(r8), allocatable :: methane_ph_patches(:)
   real(r8), allocatable :: lat_g(:), lon_g(:)
   real(r8) :: lat_deg, lon_deg

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/soil/' // trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate methane soil pH patches ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      rawfile = trim(dir_rawdata)//'/soil/PHH2O1.nc'
      inquire(file=trim(rawfile), exist=raw_exists)

      IF (.not. raw_exists) THEN
         IF (p_is_master) THEN
            write(*,'(A,A)') '  WARNING: PHH2O raw file not found: ', trim(rawfile)
            write(*,'(A)')   '           writing fallback methane_ph_patches = 6.2.'
         ENDIF
      ENDIF

      IF (p_is_worker) THEN
         allocate(methane_ph_patches(numpatch))
         methane_ph_patches(:) = methane_ph_fallback

         IF (raw_exists) THEN
            ierr = nf90_open(trim(rawfile), NF90_NOWRITE, ncid)
            IF (ierr /= NF90_NOERR) THEN
               write(*,'(A,A)') ' WARNING: PHH2O open failed in Aggregation_MethanePH: ', &
                  trim(nf90_strerror(ierr))
            ELSE
               CALL get_dim(ncid, 'lat',   nlat)
               CALL get_dim(ncid, 'lon',   nlon)
               CALL get_dim(ncid, 'depth', ndepth)

               IF (nlat <= 0 .or. nlon <= 0 .or. ndepth <= 0) THEN
                  write(*,'(A)') ' WARNING: PHH2O missing lat/lon/depth dimensions; using fallback pH.'
               ELSE
                  allocate(lat_g(nlat), lon_g(nlon))
                  ierr = nf90_inq_varid(ncid, 'lat', vid)
                  IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lat_g)
                  IF (ierr /= NF90_NOERR) THEN
                     write(*,'(A,A)') ' WARNING: PHH2O latitude read failed: ', trim(nf90_strerror(ierr))
                  ELSE
                     ierr = nf90_inq_varid(ncid, 'lon', vid)
                     IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lon_g)
                     IF (ierr /= NF90_NOERR) THEN
                        write(*,'(A,A)') ' WARNING: PHH2O longitude read failed: ', trim(nf90_strerror(ierr))
                     ELSE
                        ierr = nf90_inq_varid(ncid, 'PHH2O', vid)
                        IF (ierr /= NF90_NOERR) THEN
                           write(*,'(A,A)') ' WARNING: PHH2O variable read failed: ', trim(nf90_strerror(ierr))
                        ELSE
                           DO ipatch = 1, numpatch
                              IF (ipatch == wmo_patch(landpatch%ielm(ipatch))) THEN
                                 wmo_src = wmo_source(landpatch%ielm(ipatch))
                                 methane_ph_patches(ipatch) = methane_ph_patches(wmo_src)
                                 CYCLE
                              ENDIF
                              CALL patch_center_deg(ipatch, lat_deg, lon_deg, ok_center)
                              IF (.not. ok_center) CYCLE
                              methane_ph_patches(ipatch) = native_patch_ph(ncid, vid, lat_g, lon_g, &
                                 nlat, nlon, min(4, ndepth), lat_deg, lon_deg)
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (allocated(lat_g)) deallocate(lat_g)
                  IF (allocated(lon_g)) deallocate(lon_g)
               ENDIF

               ierr = nf90_close(ncid)
            ENDIF
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('methane_ph_patches [pH]', methane_ph_patches)
#endif

      lndname = trim(landdir)//'/methane_ph_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'methane_ph_patches', 'patch', &
         landpatch, methane_ph_patches, DEF_Srfdata_CompressLevel)

      IF (p_is_worker) THEN
         deallocate(methane_ph_patches)
      ENDIF

CONTAINS

   SUBROUTINE patch_center_deg(ipatch, lat_deg, lon_deg, ok)
      integer, intent(in) :: ipatch
      real(r8), intent(out) :: lat_deg, lon_deg
      logical, intent(out) :: ok

      integer :: ie, ipxl, ipxstt, ipxend
      real(r8) :: area, sumarea, latc, lonc, xsum, ysum

      ok = .false.
      lat_deg = 0._r8
      lon_deg = 0._r8

      IF (ipatch < 1 .or. ipatch > numpatch) RETURN
      ie     = landpatch%ielm(ipatch)
      ipxstt = landpatch%ipxstt(ipatch)
      ipxend = landpatch%ipxend(ipatch)
      IF (ie <= 0 .or. ipxstt < 1 .or. ipxend < ipxstt) RETURN

      sumarea = 0._r8
      lat_deg = 0._r8
      xsum = 0._r8
      ysum = 0._r8

      DO ipxl = ipxstt, ipxend
         area = areaquad(pixel%lat_s(mesh(ie)%ilat(ipxl)), pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                         pixel%lon_w(mesh(ie)%ilon(ipxl)), pixel%lon_e(mesh(ie)%ilon(ipxl)))
         IF (area <= 0._r8) CYCLE
         latc = 0.5_r8 * (pixel%lat_s(mesh(ie)%ilat(ipxl)) + pixel%lat_n(mesh(ie)%ilat(ipxl)))
         lonc = 0.5_r8 * (pixel%lon_w(mesh(ie)%ilon(ipxl)) + pixel%lon_e(mesh(ie)%ilon(ipxl)))
         lat_deg = lat_deg + area * latc
         xsum = xsum + area * cos(lonc * PI / 180._r8)
         ysum = ysum + area * sin(lonc * PI / 180._r8)
         sumarea = sumarea + area
      ENDDO

      IF (sumarea <= 0._r8) RETURN
      lat_deg = lat_deg / sumarea
      lon_deg = atan2(ysum, xsum) * 180._r8 / PI
      ok = .true.
   END SUBROUTINE patch_center_deg

   real(r8) FUNCTION native_patch_ph(ncid, vid, lat_g, lon_g, nlat, nlon, ndepth, lat_deg, lon_deg)
      integer, intent(in) :: ncid, vid, nlat, nlon, ndepth
      real(r8), intent(in) :: lat_g(:), lon_g(:)
      real(r8), intent(in) :: lat_deg, lon_deg

      integer :: ilat0, ilon0, ilat, ilon, d, ioff, joff, cnt, ierr_loc
      integer :: ilat_start, ilat_end, nlat_read, nlon_read, ilon_start
      integer :: nfirst, nsecond, ix, iy
      integer :: ival
      integer(1) :: bval(2*native_search_radius+1, 2*native_search_radius+1, 4)
      integer(1) :: bsegment(2*native_search_radius+1, 2*native_search_radius+1, 4)
      real(r8) :: dlat, dlon, lon_norm, sum_x10
      real(r8) :: lat_min, lat_max
      logical :: any_read, slab_failed

      native_patch_ph = methane_ph_fallback
      IF (nlat <= 0 .or. nlon <= 0 .or. ndepth <= 0) RETURN

      lat_min = min(lat_g(1), lat_g(nlat))
      lat_max = max(lat_g(1), lat_g(nlat))
      IF (lat_deg < lat_min .or. lat_deg > lat_max) RETURN

      IF (nlat > 1) THEN
         dlat = abs(lat_g(2) - lat_g(1))
      ELSE
         dlat = 1._r8
      ENDIF
      IF (nlon > 1) THEN
         dlon = abs(lon_g(2) - lon_g(1))
      ELSE
         dlon = 1._r8
      ENDIF

      IF (lat_g(1) >= lat_g(nlat)) THEN
         ilat0 = nint((lat_g(1) - lat_deg) / max(dlat, 1.e-12_r8)) + 1
      ELSE
         ilat0 = nint((lat_deg - lat_g(1)) / max(dlat, 1.e-12_r8)) + 1
      ENDIF
      ilat0 = min(max(ilat0, 1), nlat)

      lon_norm = modulo(lon_deg - lon_g(1), 360._r8) + lon_g(1)
      ilon0 = nint((lon_norm - lon_g(1)) / max(dlon, 1.e-12_r8)) + 1
      ilon0 = modulo(ilon0 - 1, nlon) + 1

      ilat_start = max(1, ilat0 - native_search_radius)
      ilat_end = min(nlat, ilat0 + native_search_radius)
      nlat_read = ilat_end - ilat_start + 1
      nlon_read = min(nlon, 2 * native_search_radius + 1)
      bval = missing_byte
      any_read = .false.
      slab_failed = .false.

      IF (nlon <= 2 * native_search_radius) THEN
         ierr_loc = nf90_get_var(ncid, vid, bsegment(1:nlon_read,1:nlat_read,1:ndepth), &
            start=[1, ilat_start, 1], &
            count=[nlon_read, nlat_read, ndepth])
         IF (ierr_loc == NF90_NOERR) THEN
            bval(1:nlon_read,1:nlat_read,1:ndepth) = bsegment(1:nlon_read,1:nlat_read,1:ndepth)
            any_read = .true.
         ELSE
            slab_failed = .true.
         ENDIF
      ELSE
         ilon_start = modulo(ilon0 - native_search_radius - 1, nlon) + 1
         nfirst = min(nlon_read, nlon - ilon_start + 1)
         ierr_loc = nf90_get_var(ncid, vid, bsegment(1:nfirst,1:nlat_read,1:ndepth), &
            start=[ilon_start, ilat_start, 1], &
            count=[nfirst, nlat_read, ndepth])
         IF (ierr_loc == NF90_NOERR) THEN
            bval(1:nfirst,1:nlat_read,1:ndepth) = bsegment(1:nfirst,1:nlat_read,1:ndepth)
            any_read = .true.
         ELSE
            slab_failed = .true.
         ENDIF

         nsecond = nlon_read - nfirst
         IF (nsecond > 0) THEN
            ierr_loc = nf90_get_var(ncid, vid, bsegment(1:nsecond,1:nlat_read,1:ndepth), &
               start=[1, ilat_start, 1], &
               count=[nsecond, nlat_read, ndepth])
            IF (ierr_loc == NF90_NOERR) THEN
               bval(nfirst+1:nlon_read,1:nlat_read,1:ndepth) = &
                  bsegment(1:nsecond,1:nlat_read,1:ndepth)
               any_read = .true.
            ELSE
               slab_failed = .true.
            ENDIF
         ENDIF
      ENDIF

      IF (slab_failed) THEN
         bval = missing_byte
         any_read = .false.
         DO d = 1, ndepth
            DO joff = -native_search_radius, native_search_radius
               ilat = ilat0 + joff
               IF (ilat < 1 .or. ilat > nlat) CYCLE
               iy = ilat - ilat_start + 1
               DO ioff = -native_search_radius, native_search_radius
                  ilon = modulo(ilon0 + ioff - 1, nlon) + 1
                  IF (nlon <= 2 * native_search_radius) THEN
                     ix = ilon
                  ELSE
                     ix = ioff + native_search_radius + 1
                  ENDIF
                  ierr_loc = nf90_get_var(ncid, vid, bsegment(1:1,1:1,1:1), &
                     start=[ilon, ilat, d], count=[1, 1, 1])
                  IF (ierr_loc /= NF90_NOERR) CYCLE
                  bval(ix,iy,d) = bsegment(1,1,1)
                  any_read = .true.
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      IF (.not. any_read) RETURN

      sum_x10 = 0._r8
      cnt = 0
      DO d = 1, ndepth
         DO joff = -native_search_radius, native_search_radius
            ilat = ilat0 + joff
            IF (ilat < 1 .or. ilat > nlat) CYCLE
            iy = ilat - ilat_start + 1
            DO ioff = -native_search_radius, native_search_radius
               ilon = modulo(ilon0 + ioff - 1, nlon) + 1
               IF (nlon <= 2 * native_search_radius) THEN
                  ix = ilon
               ELSE
                  ix = ioff + native_search_radius + 1
               ENDIF
               IF (bval(ix,iy,d) == missing_byte) CYCLE
               IF (bval(ix,iy,d) >= 0_1) THEN
                  ival = int(bval(ix,iy,d))
               ELSE
                  ival = int(bval(ix,iy,d)) + 256
               ENDIF
               IF (ival <= 0) CYCLE
               sum_x10 = sum_x10 + real(ival, r8)
               cnt = cnt + 1
            ENDDO
         ENDDO
      ENDDO

      IF (cnt > 0) THEN
         native_patch_ph = 0.1_r8 * sum_x10 / real(cnt, r8)
         native_patch_ph = max(2._r8, min(10._r8, native_patch_ph))
      ENDIF
   END FUNCTION native_patch_ph

   SUBROUTINE get_dim(ncid, name, n)
      integer, intent(in)  :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: n
      integer :: did, ierr_loc

      n = 0
      ierr_loc = nf90_inq_dimid(ncid, trim(name), did)
      IF (ierr_loc /= NF90_NOERR) RETURN
      ierr_loc = nf90_inquire_dimension(ncid, did, len=n)
      IF (ierr_loc /= NF90_NOERR) n = 0
   END SUBROUTINE get_dim

END SUBROUTINE Aggregation_MethanePH
