#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Methane_GIEMS

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Loads GIEMS-MC v1.1 monthly inundated+saturated wetland fraction
!   (Prigent+2020) and exposes a per-patch true monthly time series to
!   the methane physics. Activated through `DEF_wetland_finundation_scheme = 5`.
!
! INPUT FILE:
!   `DEF_file_GIEMS` (namelist); required when scheme 5 is selected.
!   Expected variable: `inund_sat_wetland_frac(time, latitude, longitude)`
!   - time: months since 1992-01-01 (monthly)
!            348 months = 1992-01 to 2020-12
!   - latitude/longitude: degrees, regular grid (any resolution)
!   - Sentinel values: -999 ocean, -998 snow, -997 urban; treated as 0.
!
! INIT ALGORITHM (called once during methane init):
!   1. Master opens file, reads dims + lat/lon coords (small).
!   2. Broadcast dims + coords; allocate per-rank ts (small).
!   3. Each rank computes nearest-pixel (best_ix, best_iy) per patch.
!   4. Streaming loop over 348 months:
!      a. Master reads one (nlon, nlat) slab (~4 MB).
!      b. MPI_Bcast the slab.
!      c. Each rank fills giems_ts_wetland_frac(t, ipatch) from its
!         cached (best_ix, best_iy) indices.
!   5. Build 12-month climatology fallback from full ts.
!
! RUNTIME:
!   FUNCTION giems_finundated(year, day_of_year) returns:
!   - The actual (year, month) value when 1992 <= year <= 2020
!     (captures El Nino/La Nina anomalies, long-term trends)
!   - Else: 12-month climatology fallback
!
! Notes:
!   - Streaming read keeps master memory peak at one slab (~4 MB).
!   - Per-rank storage: 348 * numpatch_local * 8 bytes
!     (10K patches -> 28 MB; SA test trivial).
!   - Nearest-pixel mapping; for coarser CoLM grids than GIEMS, a
!     better mapping would area-average. Acceptable for now.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: PI

   IMPLICIT NONE
   SAVE
   PRIVATE

   ! Public state
   ! True monthly time series, 1992-01-01 .. 2020-12-01 (348 months)
   real(r8), allocatable, public :: giems_ts_wetland_frac(:,:)   ! (ntime, numpatch) [0-1]
   ! 12-month climatology fallback for years outside 1992-2020
   real(r8), allocatable, public :: giems_clim_wetland_frac(:,:) ! (12, numpatch) [0-1]
   integer,  public :: giems_year_start = 1992    ! GIEMS-MC v1.1 first year
   integer,  public :: giems_year_end   = 2020    ! GIEMS-MC v1.1 last year
   integer,  public :: giems_ntime      = 0       ! actual N months read
   logical,  public :: giems_active = .false.

   ! Public API
   PUBLIC :: allocate_methane_giems
   PUBLIC :: deallocate_methane_giems
   PUBLIC :: read_methane_giems
   PUBLIC :: giems_finundated

CONTAINS

   SUBROUTINE allocate_methane_giems(numpatch)
      integer, intent(in) :: numpatch
      IF (numpatch <= 0) RETURN
      IF (allocated(giems_clim_wetland_frac)) THEN
         IF (size(giems_clim_wetland_frac,1) == 12 .and. size(giems_clim_wetland_frac,2) == numpatch) RETURN
         deallocate(giems_clim_wetland_frac)
      ENDIF
      allocate(giems_clim_wetland_frac(12, numpatch))
      giems_clim_wetland_frac(:,:) = 0._r8
   END SUBROUTINE allocate_methane_giems

   SUBROUTINE deallocate_methane_giems()
      IF (allocated(giems_ts_wetland_frac))   deallocate(giems_ts_wetland_frac)
      IF (allocated(giems_clim_wetland_frac)) deallocate(giems_clim_wetland_frac)
      giems_active = .false.
   END SUBROUTINE deallocate_methane_giems

   SUBROUTINE read_methane_giems(file_giems, patchlatr_in, patchlonr_in, numpatch)
      ! Stream-read full GIEMS-MC monthly time series + build fallback
      ! climatology.  Memory-efficient: master holds only one slab at
      ! a time; each rank stores its own per-patch time series.
#ifdef USEMPI
      USE MPI
#endif
      USE netcdf

      character(len=*), intent(in) :: file_giems
      real(r8), intent(in) :: patchlatr_in(:), patchlonr_in(:)   ! radians
      integer, intent(in)  :: numpatch

      character(len=*), parameter :: varname = 'inund_sat_wetland_frac'
      character(len=*), parameter :: latname = 'latitude'
      character(len=*), parameter :: lonname = 'longitude'

      integer :: ncid, vid, ierr, ierr_bc
      integer :: ntime, nlat, nlon
      integer :: t, mo, ilat, ilon, ipatch
      integer :: ndims, vdims(3), dlen_lon
      real(r4), allocatable :: slab(:,:)
      real(r8), allocatable :: lat_g(:), lon_g(:)
      integer,  allocatable :: best_ix(:), best_iy(:)
      integer,  allocatable :: ccnt(:,:)            ! valid-count per (month, patch) for climatology
      real(r8) :: lat_deg, lon_deg, dlatd, dlond, dmin, d
      logical  :: fexists, found_match
      character(len=256) :: dname
      real(r4) :: v
      logical, save :: warned_open = .false.

      ! IMPORTANT: do NOT RETURN here on numpatch<=0.  Master rank often
      ! has no patches but must still participate in the MPI_Bcast
      ! collectives below; early return on master causes worker MPI
      ! deadlock.
      giems_active = .false.
      giems_ntime = 0

      ! ---- Step 1: master opens file, reads coords ----
      IF (p_is_master) THEN
         INQUIRE(file=trim(file_giems), exist=fexists)
         IF (.not. fexists) THEN
            IF (.not. warned_open) THEN
               write(*,'(A,A)') ' WARNING: GIEMS file not found: ', trim(file_giems)
               warned_open = .true.
            ENDIF
            ntime = 0; nlat = 0; nlon = 0
         ELSE
            ierr = nf90_open(trim(file_giems), NF90_NOWRITE, ncid)
            IF (ierr /= NF90_NOERR) THEN
               write(*,'(A,A)') ' WARNING: GIEMS open failed: ', trim(nf90_strerror(ierr))
               ntime = 0; nlat = 0; nlon = 0
            ELSE
               CALL get_dim(ncid, 'time',  ntime)
               CALL get_dim(ncid, latname, nlat)
               CALL get_dim(ncid, lonname, nlon)
               IF (ntime <= 0 .or. nlat <= 0 .or. nlon <= 0) THEN
                  write(*,'(A)') ' WARNING: GIEMS missing required dimensions; methane scheme 5 inactive.'
                  ierr = nf90_close(ncid)
               ELSE
                  write(*,'(A,I0,A,I0,A,I0)') ' GIEMS dims: time=', ntime, ' lat=', nlat, ' lon=', nlon

                  allocate(lat_g(nlat), lon_g(nlon))
                  ierr = nf90_inq_varid(ncid, latname, vid)
                  IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lat_g)
                  IF (ierr /= NF90_NOERR) THEN
                     write(*,'(A,A)') ' WARNING: GIEMS latitude read failed: ', trim(nf90_strerror(ierr))
                  ELSE
                     ierr = nf90_inq_varid(ncid, lonname, vid)
                     IF (ierr == NF90_NOERR) ierr = nf90_get_var(ncid, vid, lon_g)
                     IF (ierr /= NF90_NOERR) THEN
                        write(*,'(A,A)') ' WARNING: GIEMS longitude read failed: ', trim(nf90_strerror(ierr))
                     ELSE
                        ierr = nf90_inq_varid(ncid, varname, vid)
                        IF (ierr /= NF90_NOERR) THEN
                           write(*,'(A,A)') ' WARNING: GIEMS variable not found: ', trim(nf90_strerror(ierr))
                        ELSE
                           ierr = nf90_inquire_variable(ncid, vid, ndims=ndims, dimids=vdims)
                           IF (ierr == NF90_NOERR .and. ndims == 3) THEN
                              ierr = nf90_inquire_dimension(ncid, vdims(1), name=dname, len=dlen_lon)
                              IF (ierr == NF90_NOERR) THEN
                                 IF (index(dname,'lon') == 0 .and. dlen_lon /= nlon) THEN
                                    write(*,'(A,A,A,I0)') ' WARNING: GIEMS dim 1 (Fortran fastest) is "', &
                                       trim(dname), '" len=', dlen_lon, ' (expected longitude=', nlon, ')'
                                 ENDIF
                              ENDIF
                           ENDIF
                           giems_active = .true.
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (.not. giems_active) THEN
                     IF (allocated(lat_g)) deallocate(lat_g)
                     IF (allocated(lon_g)) deallocate(lon_g)
                     ierr = nf90_close(ncid)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL MPI_Bcast(giems_active, 1, MPI_LOGICAL, p_address_master, p_comm_glb, ierr)
      CALL MPI_Bcast(ntime,        1, MPI_INTEGER, p_address_master, p_comm_glb, ierr)
      CALL MPI_Bcast(nlat,         1, MPI_INTEGER, p_address_master, p_comm_glb, ierr)
      CALL MPI_Bcast(nlon,         1, MPI_INTEGER, p_address_master, p_comm_glb, ierr)
      IF (giems_active) THEN
         IF (.not. p_is_master) THEN
            allocate(lat_g(nlat), lon_g(nlon))
         ENDIF
         CALL MPI_Bcast(lat_g, nlat, MPI_DOUBLE_PRECISION, p_address_master, p_comm_glb, ierr)
         CALL MPI_Bcast(lon_g, nlon, MPI_DOUBLE_PRECISION, p_address_master, p_comm_glb, ierr)
      ENDIF
#endif

      IF (.not. giems_active) RETURN
      giems_ntime = ntime
      giems_year_end = giems_year_start + (ntime - 1) / 12

      ! patchlatr_in/patchlonr_in are worker-local coordinate vectors.
      ! Non-worker ranks call this collective routine with numpatch=0; worker
      ! ranks must provide arrays long enough for their local patch count.
      IF (numpatch < 0) THEN
         write(*,'(A,I0)') ' ERROR: read_methane_giems received negative numpatch=', numpatch
         CALL CoLM_stop ()
      ENDIF
      IF (numpatch > 0) THEN
         IF (size(patchlatr_in) < numpatch .or. size(patchlonr_in) < numpatch) THEN
            write(*,'(A,I0,A,I0,A,I0)') &
               ' ERROR: read_methane_giems patch-coordinate size mismatch: numpatch=', numpatch, &
               ' size(patchlatr_in)=', size(patchlatr_in), &
               ' size(patchlonr_in)=', size(patchlonr_in)
            CALL CoLM_stop ()
         ENDIF
      ENDIF

      ! ---- Step 2: each rank computes nearest-pixel per patch (cached) ----
      ! Split the lookup into independent nearest-lat and nearest-lon scans:
      ! O(npatch*(nlat+nlon)) rather than the old nested
      ! O(npatch*nlat*nlon).  This removes a global-run init hot spot while
      ! preserving the same nearest-gridcell mapping on a regular lon/lat
      ! GIEMS grid.
      allocate(best_ix(numpatch))
      allocate(best_iy(numpatch))
      best_ix(:) = -1
      best_iy(:) = -1
      DO ipatch = 1, numpatch
         lat_deg = patchlatr_in(ipatch) * 180._r8 / PI
         lon_deg = patchlonr_in(ipatch) * 180._r8 / PI
         found_match = .false.

         dmin = huge(0._r8)
         DO ilat = 1, nlat
            dlatd = abs(lat_deg - lat_g(ilat))
            IF (dlatd < dmin) THEN
               dmin = dlatd
               best_iy(ipatch) = ilat
            ENDIF
         ENDDO
         IF (dmin > 5._r8) THEN
            best_iy(ipatch) = -1
            CYCLE
         ENDIF

         dmin = huge(0._r8)
         DO ilon = 1, nlon
            dlond = modulo(abs(lon_deg - lon_g(ilon)), 360._r8)
            dlond = min(dlond, 360._r8 - dlond)
            IF (dlond < dmin) THEN
               dmin = dlond
               best_ix(ipatch) = ilon
            ENDIF
         ENDDO
         IF (dmin <= 5._r8) found_match = .true.
         IF (.not. found_match) THEN
            best_ix(ipatch) = -1
            best_iy(ipatch) = -1
         ENDIF
      ENDDO

      ! ---- Step 3: allocate per-rank time series + climatology accumulator ----
      IF (allocated(giems_ts_wetland_frac)) deallocate(giems_ts_wetland_frac)
      allocate(giems_ts_wetland_frac(ntime, numpatch))
      giems_ts_wetland_frac(:,:) = 0._r8

      IF (allocated(giems_clim_wetland_frac)) THEN
         IF (size(giems_clim_wetland_frac,1) /= 12 .or. size(giems_clim_wetland_frac,2) /= numpatch) THEN
            deallocate(giems_clim_wetland_frac)
         ENDIF
      ENDIF
      IF (.not. allocated(giems_clim_wetland_frac)) THEN
         allocate(giems_clim_wetland_frac(12, numpatch))
      ENDIF
      giems_clim_wetland_frac(:,:) = 0._r8
      allocate(ccnt(12, numpatch))
      ccnt(:,:) = 0

      ! ---- Step 4: streaming month-by-month read + broadcast + lookup ----
      ! Dimension-order note: NetCDF metadata advertises the variable as
      ! inund_sat_wetland_frac(time, latitude, longitude) in C order, but
      ! Fortran reverses dimension order for nf90_get_var so the buffer is
      ! declared slab(nlon, nlat) and start/count are [lon, lat, time].
      ! Do NOT "fix" this by swapping nlat/nlon -- silently transposes the
      ! lookup.  init_methane_giems verifies dim names against the file at
      ! Step 1 (nf90_inq_dim_ids/names), and a future maintainer should keep
      ! that check before touching this read.
      allocate(slab(nlon, nlat))

      DO t = 1, ntime
         mo = mod(t-1, 12) + 1

         IF (p_is_master) THEN
            ierr = nf90_get_var(ncid, vid, slab, &
               start=[1, 1, t], count=[nlon, nlat, 1])
         ENDIF
#ifdef USEMPI
         CALL MPI_Bcast(ierr, 1, MPI_INTEGER, p_address_master, p_comm_glb, ierr_bc)
#endif
         IF (ierr /= NF90_NOERR) THEN
            ! Silent zeroing of a failed monthly slab would look like "no
            ! inundation that month" — produces wrong CH4 budget.
            ! Fail loud instead (scheme 5 is observation-driven; partial
            ! failure is not acceptable).
            IF (p_is_master) THEN
               write(*,'(A,I0,A)') ' ERROR: GIEMS slab read failed at time index ', t, &
                  '; aborting.'
               write(*,'(A,A)') '   nf90_strerror: ', trim(nf90_strerror(ierr))
            ENDIF
            CALL CoLM_stop ()
         ENDIF

#ifdef USEMPI
         ! slab is declared real(r4); MPI_REAL4 locks the 4-byte ABI even if a
         ! future maintainer flips r4 -> r8 and forgets to update the type tag
         ! (silent half-array truncation otherwise).
         CALL MPI_Bcast(slab, nlon*nlat, MPI_REAL4, p_address_master, p_comm_glb, ierr)
#endif

         DO ipatch = 1, numpatch
            IF (best_ix(ipatch) < 1 .or. best_iy(ipatch) < 1) CYCLE
            v = slab(best_ix(ipatch), best_iy(ipatch))
            IF (v >= 0._r4 .and. v <= 1._r4) THEN
               giems_ts_wetland_frac(t, ipatch) = real(v, r8)
               giems_clim_wetland_frac(mo, ipatch) = &
                  giems_clim_wetland_frac(mo, ipatch) + real(v, r8)
               ccnt(mo, ipatch) = ccnt(mo, ipatch) + 1
            ENDIF
            ! else: ocean/snow/urban sentinel -> leave as 0
         ENDDO
      ENDDO

      ! Finalize climatology (per-patch monthly mean)
      DO ipatch = 1, numpatch
         DO mo = 1, 12
            IF (ccnt(mo, ipatch) > 0) THEN
               giems_clim_wetland_frac(mo, ipatch) = &
                  giems_clim_wetland_frac(mo, ipatch) / real(ccnt(mo, ipatch), r8)
            ENDIF
         ENDDO
      ENDDO

      IF (p_is_master) ierr = nf90_close(ncid)

      deallocate(slab, ccnt, best_ix, best_iy, lat_g, lon_g)

      IF (p_is_master) write(*,'(A,I0,A,I0,A,I0,A)') &
         ' GIEMS monthly time series loaded (', ntime, ' months, ', &
         giems_year_start, '-', giems_year_end, &
         '); years outside this range use 12-month climatology fallback.'
   END SUBROUTINE read_methane_giems

   real(r8) FUNCTION giems_finundated(ipatch, year, day_of_year)
      ! Returns per-patch finundated for the (year, day_of_year):
      ! - True monthly value when 1992 <= year <= 2020 (captures El Nino,
      !   long-term trend, inter-annual variability).
      ! - 12-month climatology fallback when year outside that range.
      USE MOD_TimeManager, only: julian2monthday
      integer, intent(in) :: ipatch, year, day_of_year
      integer :: mon, mday, t

      giems_finundated = 0._r8
      IF (.not. giems_active) RETURN
      IF (ipatch < 1) RETURN

      CALL julian2monthday(year, day_of_year, mon, mday)
      IF (mon < 1 .or. mon > 12) RETURN

      ! Preferred: true monthly time series for years in data range
      IF (allocated(giems_ts_wetland_frac) .and. &
          ipatch <= size(giems_ts_wetland_frac, 2)) THEN
         IF (year >= giems_year_start .and. year <= giems_year_end) THEN
            t = (year - giems_year_start) * 12 + mon
            IF (t >= 1 .and. t <= size(giems_ts_wetland_frac, 1)) THEN
               giems_finundated = giems_ts_wetland_frac(t, ipatch)
               RETURN
            ENDIF
         ENDIF
      ENDIF

      ! Fallback: 12-month climatology (year < 1992 or > 2020)
      IF (allocated(giems_clim_wetland_frac) .and. &
          ipatch <= size(giems_clim_wetland_frac, 2)) THEN
         giems_finundated = giems_clim_wetland_frac(mon, ipatch)
      ENDIF
   END FUNCTION giems_finundated

   SUBROUTINE get_dim(ncid, name, n)
      USE netcdf
      integer, intent(in)  :: ncid
      character(len=*), intent(in) :: name
      integer, intent(out) :: n
      integer :: did, ierr
      n = 0
      ierr = nf90_inq_dimid(ncid, trim(name), did)
      IF (ierr /= NF90_NOERR) RETURN
      ierr = nf90_inquire_dimension(ncid, did, len=n)
      IF (ierr /= NF90_NOERR) n = 0
   END SUBROUTINE get_dim

END MODULE MOD_Tracer_Methane_GIEMS
#endif
