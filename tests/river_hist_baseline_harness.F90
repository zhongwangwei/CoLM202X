PROGRAM river_hist_baseline_harness

! ---------------------------------------------------------------------------
! Step 1 of .omx/plans/river-history-sharded-output.md
!
! Drives the REAL `one`-mode river-history write path (MOD_Vector_ReadWrite)
! under a real CoLM SPMD layout, so that:
!
!   * a reference `_hist_unitcat_*.nc` file is produced from the actual writer
!     (not a reimplementation), locking variable names, dimension order,
!     attributes, missing value, time axis and values (plan 1.1, 1.2);
!   * the MPI gather stage and the NetCDF write stage are timed separately
!     (plan 1.3, "MPI communication 与 NetCDF I/O 必须分开");
!   * the layout always contains at least one zero-length worker (plan 1.1).
!
! Values encode the global unit-catchment id, so any ID misplacement in a
! later sharded implementation shows up as a value mismatch rather than as a
! silently plausible field.
!
! Configuration via environment variables (all optional, see defaults below):
!   RH_TOTALNUMUCAT RH_NLON RH_NLAT RH_NGROUP RH_NVAR RH_NPTHLEV
!   RH_TOTALNPTHOUT RH_NTIME RH_REPEAT RH_OUT
!
! The master rank prints one machine-readable `RHBASE ...` line per timed
! stage, consumed by tests/run_river_hist_baseline.sh.
! ---------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   integer :: totalnumucat, nlon, nlat, ngroup, nvar, ntime, nrepeat
   integer :: npthlev, totalnpthout
   character(len=256) :: fileout

   type(pointer_int32_1d), allocatable :: ucat_data_address (:)
   integer,  allocatable :: x_ucat(:), y_ucat(:)
   integer,  allocatable :: my_gid(:)
   real(r8), allocatable :: lon_ucat(:), lat_ucat(:)
   real(r8), allocatable :: vec_local(:)
   real(r8), allocatable :: wdata(:)
   real(r8), allocatable :: bif_local(:,:), bif_wdata(:,:)
   integer,  allocatable :: bif_gid(:)

   integer :: numucat, npthout_local
   integer :: ivar, itime, irep, i
   integer :: idate(3)
   real(r8) :: t0, t_gather, t_full, t_matrix
   real(r8) :: acc_gather, acc_full, acc_matrix
   logical  :: fexists
   integer  :: itime_dummy

      CALL spmd_init ()

      CALL read_env_int ('RH_TOTALNUMUCAT', 4000,  totalnumucat)
      CALL read_env_int ('RH_NLON',          200,  nlon)
      CALL read_env_int ('RH_NLAT',          100,  nlat)
      CALL read_env_int ('RH_NGROUP',          1,  ngroup)
      CALL read_env_int ('RH_NVAR',            8,  nvar)
      CALL read_env_int ('RH_NPTHLEV',         3,  npthlev)
      CALL read_env_int ('RH_TOTALNPTHOUT',  500,  totalnpthout)
      CALL read_env_int ('RH_NTIME',           2,  ntime)
      CALL read_env_int ('RH_REPEAT',          1,  nrepeat)
      CALL read_env_str ('RH_OUT', 'river_hist_baseline_hist_unitcat_ref.nc', fileout)

      CALL divide_processes_into_groups (ngroup)

      IF (p_np_worker < 2) THEN
         IF (p_is_master) write(*,*) 'RHBASE ERROR need at least 2 workers, got ', p_np_worker
         CALL mpi_barrier (p_comm_glb, p_err)
         CALL spmd_exit ()
         STOP 1
      ENDIF

      CALL build_decomposition ()
      CALL build_coordinates ()
      CALL build_bif_decomposition ()

      ! ---- reference file skeleton, mirroring MOD_Grid_RiverLakeHist.F90:314-321
      IF (p_is_master) THEN
         inquire (file=trim(fileout), exist=fexists)
         IF (fexists) THEN
            OPEN (unit=97, file=trim(fileout), status='old')
            CLOSE (unit=97, status='delete')
         ENDIF
         CALL ncio_create_file (trim(fileout))
         CALL ncio_define_dimension (trim(fileout), 'time', 0)
         CALL ncio_define_dimension (trim(fileout), 'lat_ucat', nlat)
         CALL ncio_define_dimension (trim(fileout), 'lon_ucat', nlon)
         CALL ncio_define_dimension (trim(fileout), 'bifurcation_level',   npthlev)
         CALL ncio_define_dimension (trim(fileout), 'bifurcation_pathway', totalnpthout)
         CALL ncio_write_serial (trim(fileout), 'lat_ucat', lat_ucat, 'lat_ucat')
         CALL ncio_write_serial (trim(fileout), 'lon_ucat', lon_ucat, 'lon_ucat')
      ENDIF
      CALL mpi_barrier (p_comm_glb, p_err)

      acc_gather = 0._r8
      acc_full   = 0._r8
      acc_matrix = 0._r8

      DO irep = 1, nrepeat
         DO itime = 1, ntime

            idate = (/2000, itime, 0/)
            IF (p_is_master) CALL ncio_write_time (trim(fileout), 'time', idate, &
               itime_dummy, 'DAILY')
            CALL mpi_barrier (p_comm_glb, p_err)

            DO ivar = 1, nvar

               CALL fill_local_vector (ivar, itime)

               ! ---- (a) gather stage only: exactly the MPI half of the path
               CALL mpi_barrier (p_comm_glb, p_err)
               t0 = wtime ()
               CALL vector_gather_to_master (vec_local, numucat, totalnumucat, &
                  ucat_data_address, wdata)
               CALL mpi_barrier (p_comm_glb, p_err)
               t_gather = wtime () - t0
               IF (p_is_master .and. allocated(wdata)) deallocate (wdata)

               ! ---- (b) full gather + map2grid + serial NetCDF write
               CALL mpi_barrier (p_comm_glb, p_err)
               t0 = wtime ()
               CALL vector_gather_map2grid_and_write ( &
                  vec_local, numucat, totalnumucat, ucat_data_address, &
                  nlon, x_ucat, nlat, y_ucat, &
                  trim(fileout), varname(ivar), 'lon_ucat', 'lat_ucat', itime, &
                  'baseline synthetic unitcat field', 'm')
               CALL mpi_barrier (p_comm_glb, p_err)
               t_full = wtime () - t0

               IF (irep == nrepeat) THEN
                  acc_gather = acc_gather + t_gather
                  acc_full   = acc_full   + t_full
               ENDIF

            ENDDO

            ! ---- (c) BIF pathway matrix gather
            CALL fill_bif_matrix (itime)
            CALL mpi_barrier (p_comm_glb, p_err)
            t0 = wtime ()
            ! Pass exact-length sections: the callee asserts
            ! size(global_id) == ncol_local, and a zero-data worker must
            ! present a zero-size section rather than its padded buffer.
            CALL vector_gather_matrix_to_master ( &
               bif_local(:,1:npthout_local), npthlev, npthout_local, &
               totalnpthout, bif_gid(1:npthout_local), bif_wdata)
            CALL mpi_barrier (p_comm_glb, p_err)
            t_matrix = wtime () - t0
            IF (irep == nrepeat) acc_matrix = acc_matrix + t_matrix

            IF (p_is_master) THEN
               CALL ncio_write_serial_time (trim(fileout), 'f_bifflw_lev', itime, &
                  bif_wdata, 'bifurcation_level', 'bifurcation_pathway', 'time', 1)
               IF (itime == 1) THEN
                  CALL ncio_put_attr (trim(fileout), 'f_bifflw_lev', 'missing_value', spval)
                  CALL ncio_put_attr (trim(fileout), 'f_bifflw_lev', 'long_name', &
                     'baseline synthetic bifurcation pathway flow')
                  CALL ncio_put_attr (trim(fileout), 'f_bifflw_lev', 'units', 'm^3/s')
               ENDIF
            ENDIF
            IF (p_is_master .and. allocated(bif_wdata)) deallocate (bif_wdata)

         ENDDO
      ENDDO

      IF (p_is_master) THEN
         write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)')                            &
            'RHBASE CONFIG ranks=', p_np_glb, ' workers=', p_np_worker,        &
            ' iogroups=', ngroup, ' totalnumucat=', totalnumucat,              &
            ' nvar=', nvar, ' ntime=', ntime
         write(*,'(A,E16.8)') 'RHBASE TIME gather_only_s ', acc_gather
         write(*,'(A,E16.8)') 'RHBASE TIME full_write_s  ', acc_full
         write(*,'(A,E16.8)') 'RHBASE TIME netcdf_est_s  ', max(acc_full - acc_gather, 0._r8)
         write(*,'(A,E16.8)') 'RHBASE TIME bif_matrix_s  ', acc_matrix
         ! Root-side lower bound, plan 1.5 formulas.
         write(*,'(A,I0)') 'RHBASE MEM ucat_bytes ', &
            8*(totalnumucat + nlon*nlat + max_local_ucat())
         write(*,'(A,I0)') 'RHBASE MEM bif_bytes  ', &
            8*npthlev*(totalnpthout + max_local_path())
         write(*,'(A)') 'RHBASE DONE'
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
      CALL spmd_exit ()

CONTAINS

   ! Deterministic block decomposition that always leaves worker 0 empty, so
   ! the zero-length-worker path (plan 1.1) is exercised on every run.
   SUBROUTINE build_decomposition ()

   integer :: iw, nave, nres, istt, iend, n

      allocate (ucat_data_address (0:p_np_worker-1))

      nave = totalnumucat / (p_np_worker - 1)
      nres = mod(totalnumucat, p_np_worker - 1)

      iend = 0
      DO iw = 0, p_np_worker-1
         IF (iw == 0) THEN
            allocate (ucat_data_address(iw)%val (0))
            CYCLE
         ENDIF
         n = nave
         IF (iw <= nres) n = n + 1
         istt = iend + 1
         iend = iend + n
         allocate (ucat_data_address(iw)%val (n))
         DO i = 1, n
            ucat_data_address(iw)%val(i) = istt + i - 1
         ENDDO
      ENDDO

      IF (p_is_worker) THEN
         numucat = size(ucat_data_address(p_iam_worker)%val)
         allocate (my_gid (max(numucat,1)))
         IF (numucat > 0) my_gid(1:numucat) = ucat_data_address(p_iam_worker)%val
         allocate (vec_local (max(numucat,1)))
      ELSE
         numucat = 0
         allocate (my_gid (1))
         allocate (vec_local (1))
      ENDIF

   END SUBROUTINE build_decomposition

   SUBROUTINE build_coordinates ()

      allocate (x_ucat (totalnumucat))
      allocate (y_ucat (totalnumucat))
      DO i = 1, totalnumucat
         x_ucat(i) = mod(i-1, nlon) + 1
         y_ucat(i) = mod((i-1)/nlon, nlat) + 1
      ENDDO

      allocate (lon_ucat (nlon))
      allocate (lat_ucat (nlat))
      DO i = 1, nlon
         lon_ucat(i) = -180._r8 + 360._r8 * (real(i,r8) - 0.5_r8) / real(nlon,r8)
      ENDDO
      DO i = 1, nlat
         lat_ucat(i) =   90._r8 - 180._r8 * (real(i,r8) - 0.5_r8) / real(nlat,r8)
      ENDDO

   END SUBROUTINE build_coordinates

   ! Pathways are dealt round-robin so consecutive global ids live on different
   ! workers; rank-order reassembly therefore cannot pass by accident.
   SUBROUTINE build_bif_decomposition ()

   integer :: ip, nloc

      nloc = 0
      IF (p_is_worker .and. p_iam_worker > 0) THEN
         DO ip = 1, totalnpthout
            IF (mod(ip-1, p_np_worker-1) == p_iam_worker-1) nloc = nloc + 1
         ENDDO
      ENDIF
      npthout_local = nloc

      allocate (bif_gid (max(nloc,1)))
      allocate (bif_local (npthlev, max(nloc,1)))

      nloc = 0
      IF (p_is_worker .and. p_iam_worker > 0) THEN
         DO ip = 1, totalnpthout
            IF (mod(ip-1, p_np_worker-1) == p_iam_worker-1) THEN
               nloc = nloc + 1
               bif_gid(nloc) = ip
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE build_bif_decomposition

   ! value = global_id + var/time signature, so a misplaced id is detectable.
   SUBROUTINE fill_local_vector (ivar_in, itime_in)

   integer, intent(in) :: ivar_in, itime_in

      IF (numucat > 0) THEN
         DO i = 1, numucat
            vec_local(i) = real(my_gid(i), r8) &
               + 1.e-3_r8 * real(ivar_in, r8) + 1.e-6_r8 * real(itime_in, r8)
         ENDDO
      ENDIF

   END SUBROUTINE fill_local_vector

   SUBROUTINE fill_bif_matrix (itime_in)

   integer, intent(in) :: itime_in
   integer :: ip, il

      DO ip = 1, npthout_local
         DO il = 1, npthlev
            bif_local(il,ip) = real(bif_gid(ip), r8) &
               + 1.e-3_r8 * real(il, r8) + 1.e-6_r8 * real(itime_in, r8)
         ENDDO
      ENDDO

   END SUBROUTINE fill_bif_matrix

   integer FUNCTION max_local_ucat ()

   integer :: iw

      max_local_ucat = 0
      DO iw = 0, p_np_worker-1
         max_local_ucat = max(max_local_ucat, size(ucat_data_address(iw)%val))
      ENDDO

   END FUNCTION max_local_ucat

   integer FUNCTION max_local_path ()

      max_local_path = totalnpthout / max(p_np_worker-1, 1) + 1

   END FUNCTION max_local_path

   character(len=32) FUNCTION varname (ivar_in)

   integer, intent(in) :: ivar_in
   character(len=8) :: cnum

      write(cnum,'(I2.2)') ivar_in
      varname = 'f_ucat_synth' // trim(cnum)

   END FUNCTION varname

   ! NOT mpi_wtime: mpif.h declares it DOUBLE PRECISION, and this build uses
   ! -fdefault-real-8 without -fdefault-double-8, so gfortran promotes that
   ! declaration to real(16) and reads the library's real(8) return as garbage
   ! (observed: a constant 0). system_clock has no such ambiguity.
   real(r8) FUNCTION wtime ()

   integer(8) :: cnt, rate

      CALL system_clock (count=cnt, count_rate=rate)
      wtime = real(cnt, r8) / real(max(rate, 1_8), r8)

   END FUNCTION wtime

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

   character(len=*),   intent(in)  :: name
   character(len=*),   intent(in)  :: default_value
   character(len=256), intent(out) :: value

   character(len=256) :: buf
   integer :: length, status

      value = default_value
      CALL get_environment_variable (name, buf, length, status)
      IF (status == 0 .and. length > 0) value = buf(1:length)

   END SUBROUTINE read_env_str

END PROGRAM river_hist_baseline_harness
