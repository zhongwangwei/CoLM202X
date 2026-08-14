PROGRAM river_hist_shard_harness

! ---------------------------------------------------------------------------
! Step 2 acceptance for .omx/plans/river-history-sharded-output.md
!
!   验收：包含零数据 worker/IO group 时 collective 次序一致、不死锁；
!         所有 shard 的 ID 合集完整且互斥。
!
! Drives the real IO-group shard writers in MOD_Vector_ReadWrite under a real
! CoLM SPMD layout, with a decomposition that deliberately contains
!
!   * a zero-length worker, and
!   * (at >= 2 groups) an entire IO group that owns nothing,
!
! then has the master read every shard back and check that the union of the
! global ids is complete, mutually exclusive and in range, and that each value
! sits with the id it belongs to.
!
! Pathway ids are dealt round-robin across workers so that consecutive global
! ids live on different ranks: a writer that reassembled in rank order would
! still produce a plausible-looking file, and this check would still fail it.
!
! Environment: RH_TOTALNUMUCAT RH_NGROUP RH_NPTHLEV RH_TOTALNPTHOUT RH_OUT
! ---------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Vars_Global, only: spval
   USE netcdf

   IMPLICIT NONE

   integer :: totalnumucat, ngroup, npthlev, totalnpthout
   integer :: nlon_g, nlat_g
   character(len=256) :: fileout, fileshard

   type(route_shard_layout_type) :: ucat_layout, path_layout, resv_layout
   integer,  allocatable :: ucat_gid(:), path_gid(:)
   real(r8), allocatable :: ucat_val(:), path_val(:,:)
   integer  :: numucat, npthout_local, numresv_local
   integer,  allocatable :: resv_gid(:)
   real(r8), allocatable :: resv_val(:)
   real(r8), allocatable :: lon_h(:), lat_h(:)
   integer  :: totalnumresv, iday, itime_rec
   character(len=64) :: seg_str
   integer  :: i, il, shard_count, failures
   real(r8) :: t_first, t_last

      CALL spmd_init ()

      CALL read_env_int ('RH_TOTALNUMUCAT', 997, totalnumucat)
      CALL read_env_int ('RH_NGROUP',         2, ngroup)
      CALL read_env_int ('RH_NPTHLEV',        3, npthlev)
      CALL read_env_int ('RH_TOTALNPTHOUT', 131, totalnpthout)
      CALL read_env_int ('RH_NLON',           40, nlon_g)
      CALL read_env_int ('RH_NLAT',           30, nlat_g)
      CALL read_env_int ('RH_TOTALNUMRESV',   37, totalnumresv)
      ! A restart writes a later day into the same period.
      CALL read_env_int ('RH_DAY',             1, iday)
      CALL read_env_str ('RH_OUT', 'river_hist_shard_test.nc', fileout)

      CALL divide_processes_into_groups (ngroup)
      itime_rec = 0
      write(seg_str,'(A,I4.4,I3.3,I6.6)') 'seg', 2000, iday, 0

      shard_count = p_np_io

      CALL build_local_ucat ()
      CALL build_local_paths ()
      CALL build_local_resv ()
      CALL build_axes ()

      ! ---- layouts: built once, reused by every variable
      IF (p_is_io .or. p_is_worker) THEN
         CALL route_shard_layout_build (ucat_layout, numucat,       ucat_gid)
         CALL route_shard_layout_build (path_layout, npthout_local, path_gid)
         CALL route_shard_layout_build (resv_layout, numresv_local, resv_gid)
      ENDIF

      ! ---- each IO rank creates its own shard and writes both fields
      IF (p_is_io) THEN
         CALL route_shard_filename (fileout, p_iam_io, fileshard, trim(seg_str))
         CALL ncio_create_file (trim(fileshard))
         CALL ncio_define_dimension (trim(fileshard), 'time', 0)
         CALL ncio_define_dimension (trim(fileshard), 'unitcat_local', &
            max(ucat_layout%ntotal, 0))
         CALL ncio_define_dimension (trim(fileshard), 'bifurcation_level', npthlev)
         CALL ncio_define_dimension (trim(fileshard), 'bifurcation_pathway_local', &
            max(path_layout%ntotal, 0))
         CALL ncio_define_dimension (trim(fileshard), 'reservoir_local', &
            max(resv_layout%ntotal, 0))
         CALL ncio_write_serial (trim(fileshard), 'resv_global_index', &
            resv_layout%gid(1:max(resv_layout%ntotal,0)), 'reservoir_local')
         ! the ids this shard owns, so the aggregator never guesses
         CALL ncio_write_serial (trim(fileshard), 'ucat_ucid', &
            ucat_layout%gid(1:max(ucat_layout%ntotal,0)), 'unitcat_local')
         CALL ncio_write_serial (trim(fileshard), 'pth_global_id', &
            path_layout%gid(1:max(path_layout%ntotal,0)), 'bifurcation_pathway_local')
         ! Grid metadata, matching what create_shard_skeleton writes, so the
         ! aggregator can be exercised end-to-end against these shards.
         CALL write_grid_metadata ()
         CALL ncio_write_time (trim(fileshard), 'time', (/2000, iday, 0/), itime_rec, 'DAILY')
      ENDIF

      ! Workers never open the file, so broadcast the record index the IO rank
      ! obtained; every rank must pass the same one to the collective writers.
#ifdef USEMPI
      CALL mpi_bcast (itime_rec, 1, MPI_INTEGER, 0, p_comm_group, p_err)
#endif

      IF (p_is_io .or. p_is_worker) THEN
         CALL route_shard_write_vector (ucat_layout, ucat_val, trim(fileshard), &
            'f_ucat_shard', 'unitcat_local', itime_rec, 'synthetic unitcat field', 'm')
         CALL route_shard_write_matrix (path_layout, path_val, npthlev, trim(fileshard), &
            'f_bifflw_lev', 'bifurcation_level', 'bifurcation_pathway_local', itime_rec, &
            'synthetic pathway flow', 'm^3/s')
         ! Production reservoir variable names, so the aggregator is exercised
         ! on the shape it will actually meet.
         CALL route_shard_write_vector (resv_layout, resv_val, trim(fileshard), &
            'volresv', 'reservoir_local', itime_rec, 'reservoir water volume', 'm^3')
         CALL route_shard_write_vector (resv_layout, resv_val, trim(fileshard), &
            'qresv_in', 'reservoir_local', itime_rec, 'reservoir inflow', 'm^3/s')
         CALL route_shard_write_vector (resv_layout, resv_val, trim(fileshard), &
            'qresv_out', 'reservoir_local', itime_rec, 'reservoir outflow', 'm^3/s')
      ENDIF

      ! Identity is stamped through the same public routine production calls,
      ! with the same argument shape (grid fingerprint from the real helper,
      ! shard_count from p_np_io), so a divergence between harness and
      ! production shows up here rather than only in a real run.
      IF (p_is_io) THEN
         CALL route_shard_write_identity (trim(fileshard), trim(fileout), &
            'harness_case', '2000-001', trim(seg_str), &
            max(p_iam_io,0), max(p_np_io,1), &
            trim(route_shard_grid_fingerprint (nlon_g, nlat_g, lon_h, lat_h)), &
            0._r8, 0._r8, 0)
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

      ! ---- master verifies the union across all shards
      failures = 0
      IF (p_is_master) CALL verify_shards ()

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_io .or. p_is_worker) THEN
         CALL route_shard_layout_free (ucat_layout)
         CALL route_shard_layout_free (path_layout)
      ENDIF

      IF (p_is_master) THEN
         IF (failures == 0) THEN
            write(*,'(A)') 'RHSHARD PASS'
         ELSE
            write(*,'(A,I0)') 'RHSHARD FAIL failures=', failures
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
      CALL spmd_exit ()

      IF (p_is_master .and. failures /= 0) STOP 1

CONTAINS

   ! Worker 0 is left empty on purpose; with >= 2 groups the whole of the
   ! last group is left empty too, so an entire IO group writes a shard that
   ! owns nothing.
   SUBROUTINE build_local_ucat ()

   integer :: iw, nactive, share, mystart, k

      numucat = 0
      IF (p_is_worker) THEN
         nactive = 0
         DO iw = 0, p_np_worker-1
            IF (worker_is_active(iw)) nactive = nactive + 1
         ENDDO
         IF (worker_is_active(p_iam_worker)) THEN
            share   = totalnumucat / nactive
            mystart = 0
            DO iw = 0, p_iam_worker-1
               IF (worker_is_active(iw)) mystart = mystart + share
            ENDDO
            numucat = share
            IF (active_rank(p_iam_worker) == nactive-1) &
               numucat = totalnumucat - mystart
            allocate (ucat_gid (max(numucat,1)))
            allocate (ucat_val (max(numucat,1)))
            DO k = 1, numucat
               ucat_gid(k) = mystart + k
               ucat_val(k) = real(ucat_gid(k), r8) + 0.5_r8 + 1000._r8*real(iday,r8)
            ENDDO
         ENDIF
      ENDIF

      IF (.not. allocated(ucat_gid)) THEN
         allocate (ucat_gid (1)); ucat_gid = 0
         allocate (ucat_val (1)); ucat_val = spval
      ENDIF

   END SUBROUTINE build_local_ucat

   ! Reservoirs are dealt round-robin like pathways, so consecutive global
   ! indices live on different shards.
   SUBROUTINE build_local_resv ()

   integer :: iw, nactive, ir, k

      numresv_local = 0
      IF (p_is_worker .and. worker_is_active(p_iam_worker)) THEN
         nactive = 0
         DO iw = 0, p_np_worker-1
            IF (worker_is_active(iw)) nactive = nactive + 1
         ENDDO
         DO ir = 1, totalnumresv
            IF (mod(ir-1, nactive) == active_rank(p_iam_worker)) &
               numresv_local = numresv_local + 1
         ENDDO
         allocate (resv_gid (max(numresv_local,1)))
         allocate (resv_val (max(numresv_local,1)))
         k = 0
         DO ir = 1, totalnumresv
            IF (mod(ir-1, nactive) == active_rank(p_iam_worker)) THEN
               k = k + 1
               resv_gid(k) = ir
               resv_val(k) = real(ir, r8) + 0.25_r8 + 1000._r8*real(iday,r8)
            ENDIF
         ENDDO
      ENDIF
      IF (.not. allocated(resv_gid)) THEN
         allocate (resv_gid (1)); resv_gid = 0
         allocate (resv_val (1)); resv_val = spval
      ENDIF

   END SUBROUTINE build_local_resv

   SUBROUTINE build_axes ()
   integer :: k
      allocate (lon_h(nlon_g), lat_h(nlat_g))
      DO k = 1, nlon_g
         lon_h(k) = -180._r8 + 360._r8 * (real(k,r8) - 0.5_r8) / real(nlon_g,r8)
      ENDDO
      DO k = 1, nlat_g
         lat_h(k) =   90._r8 - 180._r8 * (real(k,r8) - 0.5_r8) / real(nlat_g,r8)
      ENDDO
   END SUBROUTINE build_axes

   SUBROUTINE build_local_paths ()

   integer :: iw, nactive, ip, k

      npthout_local = 0
      IF (p_is_worker .and. worker_is_active(p_iam_worker)) THEN
         nactive = 0
         DO iw = 0, p_np_worker-1
            IF (worker_is_active(iw)) nactive = nactive + 1
         ENDDO
         DO ip = 1, totalnpthout
            IF (mod(ip-1, nactive) == active_rank(p_iam_worker)) &
               npthout_local = npthout_local + 1
         ENDDO
         allocate (path_gid (max(npthout_local,1)))
         allocate (path_val (npthlev, max(npthout_local,1)))
         k = 0
         DO ip = 1, totalnpthout
            IF (mod(ip-1, nactive) == active_rank(p_iam_worker)) THEN
               k = k + 1
               path_gid(k) = ip
               DO il = 1, npthlev
                  path_val(il,k) = real(ip, r8) + 0.001_r8*real(il,r8) + 1000._r8*real(iday,r8)
               ENDDO
            ENDIF
         ENDDO
      ENDIF

      IF (.not. allocated(path_gid)) THEN
         allocate (path_gid (1));          path_gid = 0
         allocate (path_val (npthlev, 1)); path_val = spval
      ENDIF

   END SUBROUTINE build_local_paths

   ! Worker 0 is left empty on purpose so the zero-length path is always
   ! exercised -- except when it is the only worker, where that would leave
   ! the whole domain unassigned and test nothing.
   logical FUNCTION worker_is_active (iw)
   integer, intent(in) :: iw
      worker_is_active = (p_np_worker <= 1) .or. (iw /= 0)
   END FUNCTION worker_is_active

   integer FUNCTION active_rank (iw)
   integer, intent(in) :: iw
   integer :: k
      active_rank = -1
      IF (.not. worker_is_active(iw)) RETURN
      active_rank = 0
      DO k = 0, iw-1
         IF (worker_is_active(k)) active_rank = active_rank + 1
      ENDDO
   END FUNCTION active_rank

   !> x/y for the ids this shard owns, plus the global axes.
   SUBROUTINE write_grid_metadata ()

   integer  :: k, n
   integer,  allocatable :: xs(:), ys(:)
   real(r8), allocatable :: lon(:), lat(:)

      n = max(ucat_layout%ntotal, 0)
      allocate (xs(max(n,1))); xs = 0
      allocate (ys(max(n,1))); ys = 0
      DO k = 1, n
         xs(k) = mod(ucat_layout%gid(k)-1, nlon_g) + 1
         ys(k) = mod((ucat_layout%gid(k)-1)/nlon_g, nlat_g) + 1
      ENDDO
      CALL ncio_write_serial (trim(fileshard), 'x_ucat', xs(1:n), 'unitcat_local')
      CALL ncio_write_serial (trim(fileshard), 'y_ucat', ys(1:n), 'unitcat_local')

      allocate (lon(nlon_g), lat(nlat_g))
      DO k = 1, nlon_g
         lon(k) = -180._r8 + 360._r8 * (real(k,r8) - 0.5_r8) / real(nlon_g,r8)
      ENDDO
      DO k = 1, nlat_g
         lat(k) =   90._r8 - 180._r8 * (real(k,r8) - 0.5_r8) / real(nlat_g,r8)
      ENDDO
      CALL ncio_define_dimension (trim(fileshard), 'lon_ucat', nlon_g)
      CALL ncio_define_dimension (trim(fileshard), 'lat_ucat', nlat_g)
      CALL ncio_write_serial (trim(fileshard), 'lon_ucat', lon, 'lon_ucat')
      CALL ncio_write_serial (trim(fileshard), 'lat_ucat', lat, 'lat_ucat')
      deallocate (xs, ys, lon, lat)

   END SUBROUTINE write_grid_metadata

   SUBROUTINE verify_shards ()

   integer :: ish, ncid, did, vid, ierr, n, k, gid
   integer,  allocatable :: ids(:)
   real(r8), allocatable :: vals(:), mat(:,:)
   integer,  allocatable :: seen_ucat(:), seen_path(:)
   character(len=256) :: fname

      allocate (seen_ucat (totalnumucat));  seen_ucat = 0
      allocate (seen_path (totalnpthout));  seen_path = 0

      DO ish = 0, shard_count-1
         CALL route_shard_filename (fileout, ish, fname, trim(seg_str))

         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         IF (ierr /= NF90_NOERR) THEN
            write(*,'(A,A)') 'RHSHARD missing shard: ', trim(fname)
            failures = failures + 1
            CYCLE
         ENDIF

         CALL check_identity (ncid, ish, fname)

         ! --- unit catchments
         CALL read_ids (ncid, 'unitcat_local', 'ucat_ucid', ids, n)
         IF (n > 0) THEN
            CALL read_vals (ncid, 'f_ucat_shard', n, vals)
            DO k = 1, n
               gid = ids(k)
               IF (gid < 1 .or. gid > totalnumucat) THEN
                  write(*,'(A,I0,A,A)') 'RHSHARD ucat id out of range: ', gid, ' in ', trim(fname)
                  failures = failures + 1
                  CYCLE
               ENDIF
               IF (seen_ucat(gid) /= 0) THEN
                  write(*,'(A,I0)') 'RHSHARD duplicate ucat id: ', gid
                  failures = failures + 1
               ENDIF
               seen_ucat(gid) = seen_ucat(gid) + 1
               IF (abs(vals(k) - (real(gid,r8) + 0.5_r8 + 1000._r8*real(iday,r8))) > 1.e-9_r8) THEN
                  write(*,'(A,I0,A,E20.10)') 'RHSHARD ucat value misplaced for id ', gid, &
                     ' got ', vals(k)
                  failures = failures + 1
               ENDIF
            ENDDO
            deallocate (vals)
         ENDIF
         IF (allocated(ids)) deallocate (ids)

         ! --- bifurcation pathways
         CALL read_ids (ncid, 'bifurcation_pathway_local', 'pth_global_id', ids, n)
         IF (n > 0) THEN
            CALL read_matrix (ncid, 'f_bifflw_lev', npthlev, n, mat)
            DO k = 1, n
               gid = ids(k)
               IF (gid < 1 .or. gid > totalnpthout) THEN
                  write(*,'(A,I0)') 'RHSHARD pathway id out of range: ', gid
                  failures = failures + 1
                  CYCLE
               ENDIF
               IF (seen_path(gid) /= 0) THEN
                  write(*,'(A,I0)') 'RHSHARD duplicate pathway id: ', gid
                  failures = failures + 1
               ENDIF
               seen_path(gid) = seen_path(gid) + 1
               DO il = 1, npthlev
                  IF (abs(mat(il,k) - (real(gid,r8) + 0.001_r8*real(il,r8) + 1000._r8*real(iday,r8))) > 1.e-9_r8) THEN
                     write(*,'(A,I0,A,I0)') 'RHSHARD pathway value misplaced for id ', gid, &
                        ' level ', il
                     failures = failures + 1
                  ENDIF
               ENDDO
            ENDDO
            deallocate (mat)
         ENDIF
         IF (allocated(ids)) deallocate (ids)

         ierr = nf90_close (ncid)
      ENDDO

      IF (any(seen_ucat /= 1)) THEN
         write(*,'(A,I0,A,I0)') 'RHSHARD ucat coverage broken: missing ', &
            count(seen_ucat == 0), ' duplicated ', count(seen_ucat > 1)
         failures = failures + 1
      ENDIF
      IF (any(seen_path /= 1)) THEN
         write(*,'(A,I0,A,I0)') 'RHSHARD pathway coverage broken: missing ', &
            count(seen_path == 0), ' duplicated ', count(seen_path > 1)
         failures = failures + 1
      ENDIF

      deallocate (seen_ucat, seen_path)

   END SUBROUTINE verify_shards

   SUBROUTINE check_identity (ncid, ish, fname)

   integer,          intent(in) :: ncid, ish
   character(len=*), intent(in) :: fname
   integer :: ierr, ival
   character(len=256) :: sval

      ierr = nf90_get_att (ncid, NF90_GLOBAL, 'river_hist_shard_schema_version', ival)
      IF (ierr /= NF90_NOERR .or. ival /= ROUTE_SHARD_SCHEMA_VERSION) THEN
         write(*,'(A,A)') 'RHSHARD bad schema version in ', trim(fname)
         failures = failures + 1
      ENDIF
      ierr = nf90_get_att (ncid, NF90_GLOBAL, 'shard_index', ival)
      IF (ierr /= NF90_NOERR .or. ival /= ish) THEN
         write(*,'(A,A)') 'RHSHARD bad shard_index in ', trim(fname)
         failures = failures + 1
      ENDIF
      ierr = nf90_get_att (ncid, NF90_GLOBAL, 'shard_count', ival)
      IF (ierr /= NF90_NOERR .or. ival /= shard_count) THEN
         write(*,'(A,A)') 'RHSHARD bad shard_count in ', trim(fname)
         failures = failures + 1
      ENDIF
      sval = ''
      ierr = nf90_get_att (ncid, NF90_GLOBAL, 'segment_id', sval)
      IF (ierr /= NF90_NOERR .or. trim(sval) /= trim(seg_str)) THEN
         write(*,'(A,A)') 'RHSHARD bad segment_id in ', trim(fname)
         failures = failures + 1
      ENDIF

   END SUBROUTINE check_identity

   SUBROUTINE read_ids (ncid, dimname, varname, ids, n)

   integer,          intent(in)  :: ncid
   character(len=*), intent(in)  :: dimname, varname
   integer, allocatable, intent(out) :: ids(:)
   integer,          intent(out) :: n
   integer :: did, vid, ierr

      n = 0
      ierr = nf90_inq_dimid (ncid, trim(dimname), did)
      IF (ierr == NF90_NOERR) ierr = nf90_inquire_dimension (ncid, did, len=n)
      IF (n <= 0) THEN
         allocate (ids(1)); ids = 0; n = 0; RETURN
      ENDIF
      allocate (ids(n))
      ierr = nf90_inq_varid (ncid, trim(varname), vid)
      IF (ierr == NF90_NOERR) ierr = nf90_get_var (ncid, vid, ids)
      IF (ierr /= NF90_NOERR) THEN
         write(*,'(A,A)') 'RHSHARD cannot read ', trim(varname)
         failures = failures + 1
         n = 0
      ENDIF

   END SUBROUTINE read_ids

   SUBROUTINE read_vals (ncid, varname, n, vals)

   integer,          intent(in)  :: ncid, n
   character(len=*), intent(in)  :: varname
   real(r8), allocatable, intent(out) :: vals(:)
   integer :: vid, ierr

      allocate (vals(n))
      ierr = nf90_inq_varid (ncid, trim(varname), vid)
      IF (ierr == NF90_NOERR) ierr = nf90_get_var (ncid, vid, vals, &
         start=[1,1], count=[n,1])
      IF (ierr /= NF90_NOERR) THEN
         write(*,'(A,A)') 'RHSHARD cannot read ', trim(varname)
         failures = failures + 1
         vals = spval
      ENDIF

   END SUBROUTINE read_vals

   SUBROUTINE read_matrix (ncid, varname, nrow, ncol, mat)

   integer,          intent(in)  :: ncid, nrow, ncol
   character(len=*), intent(in)  :: varname
   real(r8), allocatable, intent(out) :: mat(:,:)
   integer :: vid, ierr

      ! CDL order is (time, coldim, rowdim); the Fortran API reverses that to
      ! (rowdim, coldim, time), which is exactly the array the writer passed,
      ! so no transpose is involved.
      allocate (mat(nrow, ncol))
      ierr = nf90_inq_varid (ncid, trim(varname), vid)
      IF (ierr == NF90_NOERR) ierr = nf90_get_var (ncid, vid, mat, &
         start=[1,1,1], count=[nrow,ncol,1])
      IF (ierr /= NF90_NOERR) THEN
         write(*,'(A,A)') 'RHSHARD cannot read ', trim(varname)
         failures = failures + 1
         mat = spval
      ENDIF

   END SUBROUTINE read_matrix

   SUBROUTINE read_env_int (name, default_value, value)
   character(len=*), intent(in)  :: name
   integer,          intent(in)  :: default_value
   integer,          intent(out) :: value
   character(len=64) :: buf
   integer :: length, status
      value = default_value
      CALL get_environment_variable (name, buf, length, status)
      IF (status == 0 .and. length > 0) READ (buf(1:length),*,iostat=status) value
      IF (status /= 0) value = default_value
   END SUBROUTINE read_env_int

   SUBROUTINE read_env_str (name, default_value, value)
   character(len=*),   intent(in)  :: name, default_value
   character(len=256), intent(out) :: value
   character(len=256) :: buf
   integer :: length, status
      value = default_value
      CALL get_environment_variable (name, buf, length, status)
      IF (status == 0 .and. length > 0) value = buf(1:length)
   END SUBROUTINE read_env_str

END PROGRAM river_hist_shard_harness
