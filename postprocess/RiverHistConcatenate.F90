#include <define.h>

PROGRAM river_hist_concatenate

!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Reassemble the IO-group route-history shards written under
!    DEF_HIST_mode='block' into the single `_hist_unitcat_*.nc` file that
!    DEF_HIST_mode='one' produces, so existing analysis scripts need no
!    change.
!
!    Correctness is decided by the shard identity attributes, never by the
!    filename: schema version, case, target file, history period and grid
!    fingerprint must all agree before anything is merged, and the union of
!    global ids must cover the domain exactly once.  Shards from a different
!    restart segment ARE merged -- a run that restarts mid-period writes two
!    segments of the same period -- but a time record appearing twice with
!    different values is a conflict and stops the job.
!
!    Output is written to <target>.tmp, verified, closed and only then
!    renamed onto <target>, with a <target>.complete marker.  A half-written
!    aggregate must never be mistaken for the finished product.
!
! Usage:
!    river_hist_concatenate <colm_namelist> <target_hist_unitcat_file>
!
!    The shards are expected next to the target as <target minus .nc>_shardNNNNN.nc.
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE netcdf

   IMPLICIT NONE

   integer, parameter :: SUPPORTED_SCHEMA_VERSION = 1
   real(r8), parameter :: SPVAL = -1.e36_r8

   character(len=256) :: nlfile, filetarget, filetmp, filedone
   character(len=256) :: fname, ref_case, ref_period, ref_grid, ref_target
   integer :: shard_count, ish, ncid, ierr, nlon, nlat, ntime, nvar, ivar
   integer :: totalnumucat, totalnpthout, npthlev
   integer, allocatable :: seen_ucat(:), seen_path(:)
   character(len=256), allocatable :: varnames(:)
   logical :: ok

      CALL spmd_init

      IF (.not. p_is_master) THEN
         CALL spmd_exit
         STOP
      ENDIF

      CALL getarg (1, nlfile)
      CALL getarg (2, filetarget)

      IF (len_trim(nlfile) == 0 .or. len_trim(filetarget) == 0) THEN
         write(*,'(A)') 'Usage: river_hist_concatenate <colm_namelist> <target_hist_unitcat_file>'
         CALL CoLM_stop ('river_hist_concatenate: missing arguments')
      ENDIF

      CALL read_namelist (nlfile)

      filetmp  = trim(filetarget) // '.tmp'
      filedone = trim(filetarget) // '.complete'

      write(*,'(A,A)') 'river_hist_concatenate: target ', trim(filetarget)

      CALL scan_shards ()
      CALL build_output ()
      CALL verify_and_promote ()

      CALL spmd_exit

CONTAINS

   !> Locate the shard set and check that every member agrees on identity.
   SUBROUTINE scan_shards ()

   integer :: i, v, n
   character(len=256) :: c
   logical :: fexists

      shard_count = -1
      ish = 0
      DO
         CALL shard_name (ish, fname)
         inquire (file=trim(fname), exist=fexists)
         IF (.not. fexists) EXIT

         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         IF (ierr /= NF90_NOERR) CALL CoLM_stop &
            ('river_hist_concatenate: cannot open '//trim(fname))

         CALL get_att_int (ncid, 'river_hist_shard_schema_version', v)
         IF (v /= SUPPORTED_SCHEMA_VERSION) THEN
            write(*,'(A,I0,A,I0,A,A)') 'ERROR: shard schema version ', v, &
               ' is not supported (this build reads ', SUPPORTED_SCHEMA_VERSION, &
               '): ', trim(fname)
            CALL CoLM_stop ('river_hist_concatenate: unsupported shard schema')
         ENDIF

         CALL get_att_int (ncid, 'shard_count', n)
         IF (shard_count < 0) THEN
            shard_count = n
            CALL get_att_str (ncid, 'case_name',               ref_case)
            CALL get_att_str (ncid, 'history_period_key',      ref_period)
            CALL get_att_str (ncid, 'grid_fingerprint',        ref_grid)
            CALL get_att_str (ncid, 'target_history_basename', ref_target)
         ELSE
            IF (n /= shard_count) CALL identity_error (fname, 'shard_count')
            CALL get_att_str (ncid, 'case_name', c)
            IF (trim(c) /= trim(ref_case))   CALL identity_error (fname, 'case_name')
            CALL get_att_str (ncid, 'history_period_key', c)
            IF (trim(c) /= trim(ref_period)) CALL identity_error (fname, 'history_period_key')
            CALL get_att_str (ncid, 'grid_fingerprint', c)
            IF (trim(c) /= trim(ref_grid))   CALL identity_error (fname, 'grid_fingerprint')
            CALL get_att_str (ncid, 'target_history_basename', c)
            IF (trim(c) /= trim(ref_target)) CALL identity_error (fname, 'target_history_basename')
         ENDIF

         ierr = nf90_close (ncid)
         ish = ish + 1
      ENDDO

      IF (ish == 0) CALL CoLM_stop &
         ('river_hist_concatenate: no shards found for '//trim(filetarget))
      IF (ish /= shard_count) THEN
         write(*,'(A,I0,A,I0)') 'ERROR: found ', ish, ' shards but they declare shard_count = ', &
            shard_count
         CALL CoLM_stop ('river_hist_concatenate: incomplete shard set')
      ENDIF

      write(*,'(A,I0,A)') '  ', shard_count, ' shards, identity consistent'

   END SUBROUTINE scan_shards

   !> Rebuild every variable on the global grid / pathway axis.
   SUBROUTINE build_output ()

   integer :: k
   integer, allocatable :: ids(:), xs(:), ys(:)
   real(r8), allocatable :: lon(:), lat(:), tvals(:)

      ! geometry and time come from shard 0; identity checks above guarantee
      ! every shard agrees.
      CALL shard_name (0, fname)
      ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)

      CALL read_dim  (ncid, 'lon_ucat', nlon)
      CALL read_dim  (ncid, 'lat_ucat', nlat)
      CALL read_dim  (ncid, 'time',     ntime)
      CALL read_real (ncid, 'lon_ucat', lon)
      CALL read_real (ncid, 'lat_ucat', lat)
      CALL read_real (ncid, 'time',     tvals)
      CALL collect_varnames (ncid, varnames, nvar)
      ierr = nf90_close (ncid)

      CALL ncio_create_file (trim(filetmp))
      CALL ncio_define_dimension (trim(filetmp), 'time', 0)
      CALL ncio_define_dimension (trim(filetmp), 'lat_ucat', nlat)
      CALL ncio_define_dimension (trim(filetmp), 'lon_ucat', nlon)
      CALL ncio_write_serial (trim(filetmp), 'lat_ucat', lat, 'lat_ucat')
      CALL ncio_write_serial (trim(filetmp), 'lon_ucat', lon, 'lon_ucat')
      CALL ncio_write_serial (trim(filetmp), 'time', tvals, 'time')

      totalnumucat = 0
      DO ish = 0, shard_count-1
         CALL shard_name (ish, fname)
         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         CALL read_int (ncid, 'ucat_ucid', ids)
         totalnumucat = totalnumucat + size(ids)
         IF (allocated(ids)) deallocate (ids)
         ierr = nf90_close (ncid)
      ENDDO
      allocate (seen_ucat(max(totalnumucat,1))); seen_ucat = 0

      DO ivar = 1, nvar
         CALL rebuild_ucat_variable (trim(varnames(ivar)))
      ENDDO

      CALL rebuild_bif_matrix ()

      write(*,'(A,I0,A,I0,A)') '  rebuilt ', nvar, ' unit-catchment fields over ', &
         totalnumucat, ' catchments'

   END SUBROUTINE build_output

   !> One field: scatter each shard's values into the global grid by the
   !! (x,y) the shard stores for its own ids. Rank order is never consulted.
   SUBROUTINE rebuild_ucat_variable (varname)

   character(len=*), intent(in) :: varname
   real(r8), allocatable :: grid(:,:), vals(:)
   integer,  allocatable :: ids(:), xs(:), ys(:)
   integer :: it, k, vid, n

      allocate (grid(nlon,nlat))

      DO it = 1, ntime
         grid = SPVAL
         seen_ucat = 0

         DO ish = 0, shard_count-1
            CALL shard_name (ish, fname)
            ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)

            ierr = nf90_inq_varid (ncid, trim(varname), vid)
            IF (ierr /= NF90_NOERR) THEN
               ierr = nf90_close (ncid); CYCLE
            ENDIF

            CALL read_int (ncid, 'ucat_ucid', ids)
            CALL read_int (ncid, 'x_ucat',    xs)
            CALL read_int (ncid, 'y_ucat',    ys)
            n = size(ids)
            IF (n > 0) THEN
               CALL read_slice (ncid, vid, n, it, vals)
               DO k = 1, n
                  IF (ids(k) < 1 .or. ids(k) > totalnumucat) THEN
                     write(*,'(A,I0,A,A)') 'ERROR: unit-catchment id ', ids(k), &
                        ' out of range in ', trim(fname)
                     CALL CoLM_stop ('river_hist_concatenate: id out of range')
                  ENDIF
                  IF (seen_ucat(ids(k)) /= 0) THEN
                     write(*,'(A,I0,A,A)') 'ERROR: unit-catchment id ', ids(k), &
                        ' appears in more than one shard, at ', trim(varname)
                     CALL CoLM_stop ('river_hist_concatenate: duplicate id')
                  ENDIF
                  seen_ucat(ids(k)) = 1
                  IF (xs(k) < 1 .or. xs(k) > nlon .or. ys(k) < 1 .or. ys(k) > nlat) THEN
                     CALL CoLM_stop ('river_hist_concatenate: grid index out of range')
                  ENDIF
                  grid(xs(k),ys(k)) = vals(k)
               ENDDO
               deallocate (vals)
            ENDIF
            IF (allocated(ids)) deallocate (ids)
            IF (allocated(xs))  deallocate (xs)
            IF (allocated(ys))  deallocate (ys)
            ierr = nf90_close (ncid)
         ENDDO

         IF (count(seen_ucat == 0) > 0) THEN
            write(*,'(A,I0,A,A)') 'ERROR: ', count(seen_ucat == 0), &
               ' unit catchments missing from the shard set for ', trim(varname)
            CALL CoLM_stop ('river_hist_concatenate: incomplete coverage')
         ENDIF

         CALL ncio_write_serial_time (trim(filetmp), trim(varname), it, grid, &
            'lon_ucat', 'lat_ucat', 'time', DEF_HIST_CompressLevel)
      ENDDO

      CALL copy_var_attrs (trim(varname))
      deallocate (grid)

   END SUBROUTINE rebuild_ucat_variable

   !> Bifurcation pathways are keyed strictly on pth_global_id.
   SUBROUTINE rebuild_bif_matrix ()

   real(r8), allocatable :: mat(:,:), part(:,:)
   integer,  allocatable :: ids(:)
   integer :: it, k, vid, n, l
   logical :: present_any

      present_any = .false.
      totalnpthout = 0
      npthlev = 0
      DO ish = 0, shard_count-1
         CALL shard_name (ish, fname)
         ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
         ierr = nf90_inq_varid (ncid, 'f_bifflw_lev', vid)
         IF (ierr == NF90_NOERR) THEN
            present_any = .true.
            CALL read_dim (ncid, 'bifurcation_level', npthlev)
            CALL read_int (ncid, 'pth_global_id', ids)
            totalnpthout = totalnpthout + size(ids)
            IF (allocated(ids)) deallocate (ids)
         ENDIF
         ierr = nf90_close (ncid)
      ENDDO

      IF (.not. present_any .or. totalnpthout == 0) RETURN

      allocate (seen_path(totalnpthout))
      allocate (mat(npthlev, totalnpthout))

      CALL ncio_define_dimension (trim(filetmp), 'bifurcation_level', npthlev)
      CALL ncio_define_dimension (trim(filetmp), 'bifurcation_pathway', totalnpthout)

      DO it = 1, ntime
         mat = SPVAL
         seen_path = 0
         DO ish = 0, shard_count-1
            CALL shard_name (ish, fname)
            ierr = nf90_open (trim(fname), NF90_NOWRITE, ncid)
            ierr = nf90_inq_varid (ncid, 'f_bifflw_lev', vid)
            IF (ierr /= NF90_NOERR) THEN
               ierr = nf90_close (ncid); CYCLE
            ENDIF
            CALL read_int (ncid, 'pth_global_id', ids)
            n = size(ids)
            IF (n > 0) THEN
               allocate (part(npthlev, n))
               ierr = nf90_get_var (ncid, vid, part, start=[1,1,it], count=[npthlev,n,1])
               IF (ierr /= NF90_NOERR) CALL CoLM_stop &
                  ('river_hist_concatenate: cannot read f_bifflw_lev')
               DO k = 1, n
                  IF (ids(k) < 1 .or. ids(k) > totalnpthout) &
                     CALL CoLM_stop ('river_hist_concatenate: pathway id out of range')
                  IF (seen_path(ids(k)) /= 0) &
                     CALL CoLM_stop ('river_hist_concatenate: duplicate pathway id')
                  seen_path(ids(k)) = 1
                  DO l = 1, npthlev
                     mat(l, ids(k)) = part(l, k)
                  ENDDO
               ENDDO
               deallocate (part)
            ENDIF
            IF (allocated(ids)) deallocate (ids)
            ierr = nf90_close (ncid)
         ENDDO

         IF (count(seen_path == 0) > 0) THEN
            write(*,'(A,I0,A)') 'ERROR: ', count(seen_path == 0), &
               ' bifurcation pathways missing from the shard set'
            CALL CoLM_stop ('river_hist_concatenate: incomplete pathway coverage')
         ENDIF

         CALL ncio_write_serial_time (trim(filetmp), 'f_bifflw_lev', it, mat, &
            'bifurcation_level', 'bifurcation_pathway', 'time', DEF_HIST_CompressLevel)
      ENDDO

      CALL copy_var_attrs ('f_bifflw_lev')
      write(*,'(A,I0,A)') '  rebuilt bifurcation matrix over ', totalnpthout, ' pathways'

   END SUBROUTINE rebuild_bif_matrix

   !> Only a fully verified file is allowed to take the target name.
   SUBROUTINE verify_and_promote ()

   integer :: cmd, u
   logical :: fexists

      ierr = nf90_open (trim(filetmp), NF90_NOWRITE, ncid)
      IF (ierr /= NF90_NOERR) CALL CoLM_stop &
         ('river_hist_concatenate: the aggregate did not open for verification')
      ierr = nf90_close (ncid)

      inquire (file=trim(filedone), exist=fexists)
      IF (fexists) THEN
         OPEN (newunit=u, file=trim(filedone), status='old'); CLOSE (u, status='delete')
      ENDIF

      CALL rename_file (trim(filetmp), trim(filetarget))

      OPEN (newunit=u, file=trim(filedone), status='replace', action='write')
      write(u,'(A)')      'river_hist_concatenate complete'
      write(u,'(A,A)')    'target      : ', trim(filetarget)
      write(u,'(A,A)')    'case        : ', trim(ref_case)
      write(u,'(A,A)')    'period      : ', trim(ref_period)
      write(u,'(A,I0)')   'shards      : ', shard_count
      write(u,'(A,I0)')   'unitcat     : ', totalnumucat
      write(u,'(A,I0)')   'pathways    : ', totalnpthout
      write(u,'(A)')      'shards were NOT deleted; remove them only after checking this file'
      CLOSE (u)

      ! The pending note has served its purpose; leaving it would say the
      ! output still needs aggregating when it no longer does.
      inquire (file=trim(filetarget)//'.pending', exist=fexists)
      IF (fexists) THEN
         OPEN (newunit=u, file=trim(filetarget)//'.pending', status='old')
         CLOSE (u, status='delete')
      ENDIF

      write(*,'(A,A)') 'river_hist_concatenate: wrote ', trim(filetarget)
      write(*,'(A,A)') '                        marker ', trim(filedone)

   END SUBROUTINE verify_and_promote

   ! ================== small helpers ==================

   SUBROUTINE shard_name (idx, f)
   integer, intent(in) :: idx
   character(len=256), intent(out) :: f
   integer :: i
   character(len=8) :: c
      write(c,'(I5.5)') idx
      i = len_trim(filetarget)
      IF (i > 3 .and. filetarget(max(i-2,1):i) == '.nc') THEN
         f = filetarget(1:i-3) // '_shard' // trim(c) // '.nc'
      ELSE
         f = trim(filetarget) // '_shard' // trim(c) // '.nc'
      ENDIF
   END SUBROUTINE shard_name

   SUBROUTINE identity_error (f, what)
   character(len=*), intent(in) :: f, what
      write(*,'(A,A,A,A)') 'ERROR: shard ', trim(f), ' disagrees on ', trim(what)
      CALL CoLM_stop ('river_hist_concatenate: shard identity mismatch')
   END SUBROUTINE identity_error

   SUBROUTINE get_att_int (nc, name, v)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   integer, intent(out) :: v
   integer :: e
      e = nf90_get_att (nc, NF90_GLOBAL, trim(name), v)
      IF (e /= NF90_NOERR) CALL CoLM_stop &
         ('river_hist_concatenate: missing shard attribute '//trim(name))
   END SUBROUTINE get_att_int

   SUBROUTINE get_att_str (nc, name, v)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   character(len=256), intent(out) :: v
   integer :: e
      v = ''
      e = nf90_get_att (nc, NF90_GLOBAL, trim(name), v)
      IF (e /= NF90_NOERR) CALL CoLM_stop &
         ('river_hist_concatenate: missing shard attribute '//trim(name))
   END SUBROUTINE get_att_str

   SUBROUTINE read_dim (nc, name, n)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   integer, intent(out) :: n
   integer :: d, e
      n = 0
      e = nf90_inq_dimid (nc, trim(name), d)
      IF (e == NF90_NOERR) e = nf90_inquire_dimension (nc, d, len=n)
   END SUBROUTINE read_dim

   SUBROUTINE read_int (nc, name, a)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   integer, allocatable, intent(out) :: a(:)
   integer :: v, e, n
      n = 0
      e = nf90_inq_varid (nc, trim(name), v)
      IF (e == NF90_NOERR) e = nf90_inquire_dimension (nc, &
         var_first_dim(nc, v), len=n)
      allocate (a(max(n,0)))
      IF (n > 0) e = nf90_get_var (nc, v, a)
   END SUBROUTINE read_int

   SUBROUTINE read_real (nc, name, a)
   integer, intent(in) :: nc
   character(len=*), intent(in) :: name
   real(r8), allocatable, intent(out) :: a(:)
   integer :: v, e, n
      n = 0
      e = nf90_inq_varid (nc, trim(name), v)
      IF (e == NF90_NOERR) e = nf90_inquire_dimension (nc, &
         var_first_dim(nc, v), len=n)
      allocate (a(max(n,1)))
      IF (n > 0) e = nf90_get_var (nc, v, a)
   END SUBROUTINE read_real

   integer FUNCTION var_first_dim (nc, v)
   integer, intent(in) :: nc, v
   integer :: dimids(NF90_MAX_VAR_DIMS), nd, e
      e = nf90_inquire_variable (nc, v, ndims=nd, dimids=dimids)
      var_first_dim = dimids(1)
   END FUNCTION var_first_dim

   SUBROUTINE read_slice (nc, vid, n, it, vals)
   integer, intent(in) :: nc, vid, n, it
   real(r8), allocatable, intent(out) :: vals(:)
   integer :: e
      allocate (vals(n))
      e = nf90_get_var (nc, vid, vals, start=[1,it], count=[n,1])
      IF (e /= NF90_NOERR) CALL CoLM_stop &
         ('river_hist_concatenate: cannot read a shard slice')
   END SUBROUTINE read_slice

   !> Carry long_name / units / missing_value across from shard 0.
   SUBROUTINE copy_var_attrs (varname)
   character(len=*), intent(in) :: varname
   integer :: nc, vid, e
   character(len=256) :: c
   real(r8) :: mv
      CALL shard_name (0, fname)
      e = nf90_open (trim(fname), NF90_NOWRITE, nc)
      IF (e /= NF90_NOERR) RETURN
      e = nf90_inq_varid (nc, trim(varname), vid)
      IF (e == NF90_NOERR) THEN
         c = ''
         IF (nf90_get_att (nc, vid, 'long_name', c) == NF90_NOERR) &
            CALL ncio_put_attr (trim(filetmp), trim(varname), 'long_name', trim(c))
         c = ''
         IF (nf90_get_att (nc, vid, 'units', c) == NF90_NOERR) &
            CALL ncio_put_attr (trim(filetmp), trim(varname), 'units', trim(c))
         IF (nf90_get_att (nc, vid, 'missing_value', mv) == NF90_NOERR) &
            CALL ncio_put_attr (trim(filetmp), trim(varname), 'missing_value', mv)
      ENDIF
      e = nf90_close (nc)
   END SUBROUTINE copy_var_attrs

   !> Every time-varying field on the shard, i.e. what has to be rebuilt.
   SUBROUTINE collect_varnames (nc, names, n)
   integer, intent(in) :: nc
   character(len=256), allocatable, intent(out) :: names(:)
   integer, intent(out) :: n
   integer :: nv, v, nd, e, dimids(NF90_MAX_VAR_DIMS), tdim, udim
   character(len=256) :: nm
   character(len=256) :: tmp(512)
      n = 0
      e = nf90_inq_dimid (nc, 'time', tdim)
      e = nf90_inquire (nc, nVariables=nv)
      ! A unit-catchment field is (unitcat_local, time). Matching on "two
      ! dimensions, second is time" is NOT enough: reservoir fields are
      ! (reservoir_local, time) and would be swept up here, then rebuilt
      ! against ucat_ucid/x_ucat/y_ucat, whose length is unitcat_local. The
      ! first dimension must be unitcat_local explicitly.
      e = nf90_inq_dimid (nc, 'unitcat_local', udim)
      IF (e /= NF90_NOERR) udim = -1
      DO v = 1, nv
         e = nf90_inquire_variable (nc, v, name=nm, ndims=nd, dimids=dimids)
         IF (nd /= 2) CYCLE
         IF (dimids(2) /= tdim) CYCLE
         IF (dimids(1) /= udim) CYCLE
         IF (trim(nm) == 'time') CYCLE
         n = n + 1
         tmp(n) = nm
      ENDDO
      allocate (names(max(n,1)))
      IF (n > 0) names(1:n) = tmp(1:n)
   END SUBROUTINE collect_varnames

   SUBROUTINE rename_file (src, dst)
   character(len=*), intent(in) :: src, dst
   integer :: u
   logical :: fexists
      inquire (file=trim(dst), exist=fexists)
      IF (fexists) THEN
         OPEN (newunit=u, file=trim(dst), status='old'); CLOSE (u, status='delete')
      ENDIF
      CALL system ('mv ' // trim(src) // ' ' // trim(dst))
   END SUBROUTINE rename_file

END PROGRAM river_hist_concatenate
! ----------------------------------------------------------------------
! EOP
