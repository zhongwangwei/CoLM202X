#include <define.h>

#if (defined TRACER) && (defined BGC)
MODULE MOD_Tracer_Reactive_Methane_GIEMS

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
!   3. Each rank computes nearest-pixel (best_ix, best_iy) per patch and
!      deduplicates pixels shared by local patches.  Only unique pixel
!      requests are gathered to master once.
!   4. Streaming loop over 348 months:
!      a. Master reads one (nlon, nlat) slab (~4 MB).
!      b. Master packs only values requested by patches on each rank.
!      c. MPI_Scatterv sends each rank only its unique local pixel values;
!         the rank expands them back to patches using a cached map.
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
!   - The source data and true monthly series use real(r4); the climatology
!     remains real(r8) for numerically stable accumulation.  Persistent
!     per-rank storage is (4*ntime + 8*12) * numpatch_local bytes
!     (348 months, 10K patches -> about 14.9 MB).
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
   real(r4), allocatable, public :: giems_ts_wetland_frac(:,:)   ! (ntime, numpatch) [0-1]
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
      integer :: t, mo, ilat, ilon, ipatch, irequest, pixel
      integer :: ndims, vdims(3), dim_lengths(3), metadata(4)
      integer :: total_requests, n_unique, idim, comm_size
      real(r4), allocatable :: slab(:,:)
      real(r4), allocatable :: unique_values(:), requested_values(:)
      real(r8), allocatable :: lat_g(:), lon_g(:)
      integer,  allocatable :: best_ix(:), best_iy(:)
      integer,  allocatable :: all_pixel_index(:)
      integer,  allocatable :: patch_to_unique(:), pixel_index_unique(:), pixel_to_unique(:)
      integer,  allocatable :: request_counts(:), request_displs(:)
      integer,  allocatable :: ccnt(:,:)            ! valid-count per (month, patch) for climatology
      real(r8) :: lat_deg, lon_deg, dlatd, dlond, dmin
      real(r8) :: dlat_min, dlon_min
      logical  :: fexists, found_match, layout_ok
      character(len=256) :: dim_names(3)
      real(r4) :: v
      logical, save :: warned_open = .false.

      ! IMPORTANT: do NOT RETURN here on numpatch<=0.  Master rank often
      ! has no patches but must still participate in the MPI_Bcast
      ! collectives below; early return on master causes worker MPI
      ! deadlock.
      giems_active = .false.
      giems_ntime = 0
      ntime = 0
      nlat = 0
      nlon = 0
      metadata = 0

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
                           layout_ok = .false.
                           ! Query rank first: passing a three-element dimids
                           ! buffer before rejecting rank > 3 can overrun that
                           ! buffer in some NetCDF-Fortran implementations.
                           ierr = nf90_inquire_variable(ncid, vid, ndims=ndims)
                           IF (ierr /= NF90_NOERR) THEN
                              write(*,'(A,A)') ' WARNING: GIEMS variable metadata read failed: ', &
                                 trim(nf90_strerror(ierr))
                           ELSEIF (ndims /= 3) THEN
                              write(*,'(A,I0)') ' WARNING: GIEMS variable must have exactly 3 dimensions; found ', ndims
                           ELSE
                              ierr = nf90_inquire_variable(ncid, vid, dimids=vdims)
                              layout_ok = ierr == NF90_NOERR
                              IF (.not. layout_ok) THEN
                                 write(*,'(A,A)') ' WARNING: GIEMS dimension IDs read failed: ', &
                                    trim(nf90_strerror(ierr))
                              ELSE
                                 dim_names = ''
                                 dim_lengths = 0
                                 DO idim = 1, 3
                                    ierr = nf90_inquire_dimension(ncid, vdims(idim), &
                                       name=dim_names(idim), len=dim_lengths(idim))
                                    IF (ierr /= NF90_NOERR) THEN
                                       layout_ok = .false.
                                       write(*,'(A,I0,A,A)') ' WARNING: GIEMS dimension ', idim, &
                                          ' metadata read failed: ', trim(nf90_strerror(ierr))
                                       EXIT
                                    ENDIF
                                 ENDDO
                              ENDIF
                              IF (layout_ok) THEN
                                 IF (trim(dim_names(1)) /= trim(lonname) .or. &
                                     trim(dim_names(2)) /= trim(latname) .or. &
                                     trim(dim_names(3)) /= 'time' .or. &
                                     any(dim_lengths /= [nlon, nlat, ntime])) THEN
                                    layout_ok = .false.
                                    write(*,'(A)') ' WARNING: GIEMS variable dimension order/length mismatch.'
                                    write(*,'(A,3(1X,A,1X,I0))') '   Fortran dimensions:', &
                                       trim(dim_names(1)), dim_lengths(1), &
                                       trim(dim_names(2)), dim_lengths(2), &
                                       trim(dim_names(3)), dim_lengths(3)
                                    write(*,'(A,3(1X,A,1X,I0))') '   Expected:', &
                                       trim(lonname), nlon, trim(latname), nlat, 'time', ntime
                                 ENDIF
                              ENDIF
                           ENDIF
                           giems_active = layout_ok
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

      IF (p_is_master) metadata = [merge(1, 0, giems_active), ntime, nlat, nlon]
#ifdef USEMPI
      CALL MPI_Bcast(metadata, 4, MPI_INTEGER, p_address_master, p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'metadata broadcast')
      giems_active = metadata(1) == 1
      ntime = metadata(2)
      nlat = metadata(3)
      nlon = metadata(4)
      IF (giems_active) THEN
         IF (.not. p_is_master) THEN
            allocate(lat_g(nlat), lon_g(nlon))
         ENDIF
         CALL MPI_Bcast(lat_g, nlat, MPI_DOUBLE_PRECISION, p_address_master, p_comm_glb, ierr)
         CALL check_giems_mpi(ierr, 'latitude broadcast')
         CALL MPI_Bcast(lon_g, nlon, MPI_DOUBLE_PRECISION, p_address_master, p_comm_glb, ierr)
         CALL check_giems_mpi(ierr, 'longitude broadcast')
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
         dlat_min = dmin
         IF (dlat_min > 5._r8) THEN
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
         dlon_min = dmin
         IF (sqrt(dlat_min**2 + dlon_min**2) <= 5._r8) found_match = .true.
         IF (.not. found_match) THEN
            best_ix(ipatch) = -1
            best_iy(ipatch) = -1
         ENDIF
      ENDDO

      ! ---- Step 3: allocate per-rank time series + climatology accumulator ----
      IF (allocated(giems_ts_wetland_frac)) deallocate(giems_ts_wetland_frac)
      allocate(giems_ts_wetland_frac(ntime, numpatch))
      giems_ts_wetland_frac(:,:) = 0._r4

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

      ! Cache the local flattened GIEMS pixel requested by each patch.  Many
      ! land patches can map to the same GIEMS cell, so deduplicate within
      ! each rank before gathering.  The master decodes and distributes one
      ! value per rank-local unique pixel rather than one value per patch.
      allocate(patch_to_unique(max(1, numpatch)))
      allocate(pixel_index_unique(max(1, numpatch)))
      allocate(unique_values(max(1, numpatch)))
      allocate(pixel_to_unique(max(1, nlon*nlat)))
      patch_to_unique(:) = 0
      pixel_index_unique(:) = 0
      unique_values(:) = 0._r4
      pixel_to_unique(:) = 0
      n_unique = 0
      DO ipatch = 1, numpatch
         IF (best_ix(ipatch) > 0 .and. best_iy(ipatch) > 0) THEN
            pixel = (best_iy(ipatch) - 1) * nlon + best_ix(ipatch)
            IF (pixel_to_unique(pixel) == 0) THEN
               n_unique = n_unique + 1
               pixel_index_unique(n_unique) = pixel
               pixel_to_unique(pixel) = n_unique
            ENDIF
            patch_to_unique(ipatch) = pixel_to_unique(pixel)
         ENDIF
      ENDDO
      ! Coordinates and nearest-index scratch are no longer needed once the
      ! patch-to-unique map is built; release them before the monthly slabs.
      deallocate(pixel_to_unique, best_ix, best_iy)
      deallocate(lat_g, lon_g)

      comm_size = 1
#ifdef USEMPI
      CALL MPI_Comm_size(p_comm_glb, comm_size, ierr)
      CALL check_giems_mpi(ierr, 'communicator size query')
      IF (comm_size <= 0) THEN
         write(*,'(A,I0)') ' ERROR: GIEMS invalid p_comm_glb size=', comm_size
         CALL CoLM_stop ()
      ENDIF
      IF (comm_size /= p_np_glb) THEN
         write(*,'(A,I0,A,I0)') ' ERROR: GIEMS communicator/task-count mismatch: MPI size=', &
            comm_size, ' p_np_glb=', p_np_glb
         CALL CoLM_stop ()
      ENDIF
#endif
      allocate(request_counts(comm_size), request_displs(comm_size))
      request_counts(:) = 0
      request_displs(:) = 0
#ifdef USEMPI
      CALL MPI_Gather(n_unique, 1, MPI_INTEGER, request_counts, 1, MPI_INTEGER, &
         p_address_master, p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'unique request count gather')
#else
      request_counts(1) = n_unique
#endif

      total_requests = 0
      IF (p_is_master) THEN
         DO irequest = 1, comm_size
            request_displs(irequest) = total_requests
            total_requests = total_requests + request_counts(irequest)
         ENDDO
         allocate(all_pixel_index(max(1, total_requests)))
         allocate(requested_values(max(1, total_requests)))
      ELSE
         ! MPI ignores root-only receive arguments away from the root, but
         ! allocated one-element buffers keep the legacy `use mpi` interface
         ! valid on implementations that still inspect the actual argument.
         allocate(all_pixel_index(1), requested_values(1))
      ENDIF
      all_pixel_index(:) = 0
      requested_values(:) = 0._r4

#ifdef USEMPI
      CALL MPI_Gatherv(pixel_index_unique, n_unique, MPI_INTEGER, all_pixel_index, &
         request_counts, request_displs, MPI_INTEGER, p_address_master, &
         p_comm_glb, ierr)
      CALL check_giems_mpi(ierr, 'unique pixel index gather')
#else
      IF (n_unique > 0) all_pixel_index(1:n_unique) = pixel_index_unique(1:n_unique)
#endif

      ! ---- Step 4: streaming month-by-month read + directed distribution ----
      ! Dimension-order note: NetCDF metadata advertises the variable as
      ! inund_sat_wetland_frac(time, latitude, longitude) in C order, but
      ! Fortran reverses dimension order for nf90_get_var so the buffer is
      ! declared slab(nlon, nlat) and start/count are [lon, lat, time].
      ! Do NOT "fix" this by swapping nlat/nlon -- silently transposes the
      ! lookup.  init_methane_giems verifies dim names against the file at
      ! Step 1 (nf90_inq_dim_ids/names), and a future maintainer should keep
      ! that check before touching this read.
      IF (p_is_master) allocate(slab(nlon, nlat))

      DO t = 1, ntime
         mo = mod(t-1, 12) + 1

         IF (p_is_master) THEN
            ierr = nf90_get_var(ncid, vid, slab, &
               start=[1, 1, t], count=[nlon, nlat, 1])
         ENDIF
#ifdef USEMPI
         CALL MPI_Bcast(ierr, 1, MPI_INTEGER, p_address_master, p_comm_glb, ierr_bc)
         CALL check_giems_mpi(ierr_bc, 'monthly NetCDF status broadcast')
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

         IF (p_is_master) THEN
            DO irequest = 1, total_requests
               IF (all_pixel_index(irequest) <= 0) THEN
                  requested_values(irequest) = -999._r4
               ELSE
                  ilon = mod(all_pixel_index(irequest) - 1, nlon) + 1
                  ilat = (all_pixel_index(irequest) - 1) / nlon + 1
                  requested_values(irequest) = slab(ilon, ilat)
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         ! requested_values and unique_values are real(r4); MPI_REAL4 keeps
         ! the NetCDF slab ABI explicit while sending only unique requests.
         CALL MPI_Scatterv(requested_values, request_counts, request_displs, &
            MPI_REAL4, unique_values, n_unique, MPI_REAL4, p_address_master, &
            p_comm_glb, ierr)
         CALL check_giems_mpi(ierr, 'monthly unique value scatter')
#else
         IF (n_unique > 0) unique_values(1:n_unique) = requested_values(1:n_unique)
#endif

         DO ipatch = 1, numpatch
            IF (patch_to_unique(ipatch) <= 0) CYCLE
            v = unique_values(patch_to_unique(ipatch))
            IF (v >= 0._r4 .and. v <= 1._r4) THEN
               giems_ts_wetland_frac(t, ipatch) = v
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

      IF (allocated(slab)) deallocate(slab)
      deallocate(ccnt, patch_to_unique)
      deallocate(pixel_index_unique, unique_values)
      deallocate(all_pixel_index, requested_values, request_counts, request_displs)

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

#ifdef USEMPI
   SUBROUTINE check_giems_mpi(ierr, operation)
      USE MPI
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: operation
      integer :: error_length, ierr_string
      character(len=MPI_MAX_ERROR_STRING) :: error_message

      IF (ierr == MPI_SUCCESS) RETURN
      error_message = 'unknown MPI error'
      error_length = len_trim(error_message)
      CALL MPI_Error_string(ierr, error_message, error_length, ierr_string)
      IF (ierr_string == MPI_SUCCESS) THEN
         write(*,'(A,A,A,A)') ' ERROR: GIEMS MPI failure during ', trim(operation), ': ', &
            trim(error_message(1:max(1, error_length)))
      ELSE
         write(*,'(A,A,A,I0,A,I0)') ' ERROR: GIEMS MPI failure during ', trim(operation), &
            ': code=', ierr, '; MPI_Error_string code=', ierr_string
      ENDIF
      CALL CoLM_stop ()
   END SUBROUTINE check_giems_mpi
#endif

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

END MODULE MOD_Tracer_Reactive_Methane_GIEMS
#endif
