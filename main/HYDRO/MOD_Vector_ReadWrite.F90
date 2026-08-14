#include <define.h>

MODULE MOD_Vector_ReadWrite
!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Read/Write data in vector form.
!
! Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   PUBLIC :: vector_gather_and_write
   PUBLIC :: vector_gather_matrix_to_master
   PUBLIC :: vector_gather_map2grid_and_write
   PUBLIC :: vector_read_and_scatter
   PUBLIC :: vector_read_matrix_and_scatter

   ! ------------------------------------------------------------------
   ! IO-group shard output for DEF_HIST_mode = 'block' route history.
   !
   ! The routines above funnel every route-history variable through the
   ! global master: a serialized mpi_recv loop plus one global mpi_barrier
   ! per variable.  Measured on the step-1 baseline that stage grows
   ! superlinearly with rank count (30x from 4 to 8 ranks, 16x from 8 to
   ! 12), and by 12 ranks it is 87% of hist_grid_riverlake_out.
   !
   ! The routines below replace it with the same "one group of workers ->
   ! one IO rank -> one file" MPI_Gatherv pattern the regular block
   ! history already uses (share/MOD_NetCDFVector.F90).  The master takes
   ! no part: it forms a singleton group of its own, so it must not enter
   ! any p_comm_group collective.  Within a group the IO rank is p_root,
   ! because divide_processes_into_groups splits with the global rank as
   ! key and makes the lowest-ranked member the IO rank.
   !
   ! Shards are an internal on-disk format, reassembled offline into the
   ! same single file the 'one' path writes.  They therefore carry an
   ! explicit identity (schema version, case, target file, period,
   ! restart segment, shard index/count, grid fingerprint) so that the
   ! aggregator never has to infer correctness from a filename.
   ! ------------------------------------------------------------------

   ! Bump when the shard layout changes in a way an older aggregator
   ! would misread.  The aggregator refuses versions it does not know.
   integer, PUBLIC, parameter :: ROUTE_SHARD_SCHEMA_VERSION = 1

   type, PUBLIC :: route_shard_layout_type
      logical :: built  = .false.
      integer :: nlocal = 0                 ! this rank's contribution
      integer :: ntotal = 0                 ! gathered length, IO rank only
      integer, allocatable :: cnt (:)       ! IO only, per group rank
      integer, allocatable :: dsp (:)       ! IO only, displacements
      integer, allocatable :: gid (:)       ! IO only, global ids in gathered order
   END type route_shard_layout_type

   PUBLIC :: route_shard_layout_build
   PUBLIC :: route_shard_layout_free
   PUBLIC :: route_shard_write_vector
   PUBLIC :: route_shard_write_matrix
   PUBLIC :: route_shard_write_identity
   PUBLIC :: route_shard_filename
   PUBLIC :: route_shard_grid_fingerprint

CONTAINS

   ! -------
   SUBROUTINE vector_gather_to_master ( &
         vector, vlen, totalvlen, data_address, wdata)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(in) :: vector (:)
   integer,  intent(in) :: vlen
   integer,  intent(in) :: totalvlen

   type(pointer_int32_1d), intent(in) :: data_address (0:)

   real(r8), allocatable,  intent(inout) :: wdata (:)

   ! Local variables
   integer :: iwork, mesg(2), isrc, ndata
   real(r8), allocatable :: rcache(:)

      IF (totalvlen <= 0) RETURN

      IF (p_is_master) THEN
         allocate (wdata (totalvlen))
         wdata = spval
      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, vlen/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (vlen > 0) THEN
            CALL mpi_send (vector, vlen, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
      ENDIF

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate(rcache (ndata))

               CALL mpi_recv (rcache, ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               wdata(data_address(p_itis_worker(isrc))%val) = rcache

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      wdata(data_address(0)%val) = vector
#endif

   END SUBROUTINE vector_gather_to_master

   ! -------
   SUBROUTINE vector_gather_matrix_to_master ( &
         matrix, nrow, ncol_local, ncol_global, global_id, wdata)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   real(r8), intent(in) :: matrix (:,:)
   integer,  intent(in) :: nrow
   integer,  intent(in) :: ncol_local
   integer,  intent(in) :: ncol_global
   integer,  intent(in) :: global_id (:)

   real(r8), allocatable, intent(inout) :: wdata (:,:)

   integer :: iwork, isrc, ndata, ipth, mesg(2)
   integer, allocatable :: icache(:)
   real(r8), allocatable :: rcache(:)

      IF (nrow <= 0 .or. ncol_global <= 0) RETURN

      ! Current use assumes dense global IDs in 1..ncol_global, with one
      ! column per global pathway ID across all workers.
      IF (size(global_id) /= ncol_local) THEN
         CALL CoLM_stop ('vector_gather_matrix_to_master: global_id size mismatch')
      ENDIF

      IF (p_is_master) THEN
         allocate (wdata (nrow, ncol_global))
         wdata(:,:) = 0._r8
      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, ncol_local/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (ncol_local > 0) THEN
            CALL mpi_send (global_id, ncol_local, MPI_INTEGER, &
               p_address_master, mpi_tag_mesg + 1, p_comm_glb, p_err)
            CALL mpi_send (matrix, nrow*ncol_local, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)
         ENDIF
      ENDIF

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate (icache (ndata))
               allocate (rcache (nrow*ndata))

               CALL mpi_recv (icache, ndata, MPI_INTEGER, isrc, &
                  mpi_tag_mesg + 1, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (rcache, nrow*ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO ipth = 1, ndata
                  IF (icache(ipth) < 1 .or. icache(ipth) > ncol_global) THEN
                     CALL CoLM_stop ('vector_gather_matrix_to_master: global_id out of range')
                  ENDIF
                  wdata(:, icache(ipth)) = rcache((ipth-1)*nrow+1:ipth*nrow)
               ENDDO

               deallocate (icache)
               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      DO ipth = 1, ncol_local
         IF (global_id(ipth) < 1 .or. global_id(ipth) > ncol_global) THEN
            CALL CoLM_stop ('vector_gather_matrix_to_master: global_id out of range')
         ENDIF
         wdata(:, global_id(ipth)) = matrix(:, ipth)
      ENDDO
#endif

   END SUBROUTINE vector_gather_matrix_to_master

   ! -------
   SUBROUTINE vector_gather_and_write ( vector, vlen, totalvlen, data_address, &
         fileout, varname, dimname, itime_in_file, longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(in) :: vector (:)
   integer,  intent(in) :: vlen
   integer,  intent(in) :: totalvlen

   type(pointer_int32_1d), intent(in) :: data_address (0:)

   character(len=*), intent(in) :: fileout
   character(len=*), intent(in) :: varname
   character(len=*), intent(in) :: dimname

   integer,          intent(in), optional :: itime_in_file
   character(len=*), intent(in), optional :: longname
   character(len=*), intent(in), optional :: units

   ! Local variables
   real(r8), allocatable :: wdata(:)
   logical :: write_attr


      IF (totalvlen <= 0) RETURN

      CALL vector_gather_to_master (vector, vlen, totalvlen, data_address, wdata)

      IF (p_is_master) THEN

         IF (present(itime_in_file)) THEN
            CALL ncio_write_serial_time (fileout, varname, itime_in_file, wdata, &
               dimname, 'time', DEF_HIST_CompressLevel)
         ELSE
            CALL ncio_write_serial (fileout, varname, wdata, &
               dimname, DEF_REST_CompressLevel)
         ENDIF

         IF (present(itime_in_file)) THEN
            write_attr = itime_in_file <= 1
         ELSE
            write_attr = .true.
         ENDIF

         IF (write_attr) THEN
            CALL ncio_put_attr (fileout, varname, 'missing_value', spval)
            IF (present(longname)) CALL ncio_put_attr (fileout, varname, 'long_name', longname)
            IF (present(units   )) CALL ncio_put_attr (fileout, varname, 'units',     units   )
         ENDIF

         deallocate (wdata)

      ENDIF

   END SUBROUTINE vector_gather_and_write

   ! -------
   SUBROUTINE vector_gather_map2grid_and_write ( &
         vector,  vlen,    totalvlen, data_address, nlon, x_vec,   nlat, y_vec,   &
         fileout, varname, lon_name,  lat_name,     itime_in_file, longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_DataType
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   character(len=*), intent(in) :: fileout
   real(r8),         intent(in) :: vector (:)
   integer,          intent(in) :: vlen
   integer,          intent(in) :: totalvlen
   character(len=*), intent(in) :: varname
   character(len=*), intent(in) :: lon_name, lat_name

   type(pointer_int32_1d), intent(in) :: data_address (0:)

   integer, intent(in) :: nlon, x_vec (:)
   integer, intent(in) :: nlat, y_vec (:)

   integer,          intent(in), optional :: itime_in_file
   character(len=*), intent(in), optional :: longname
   character(len=*), intent(in), optional :: units

   ! Local variables
   integer :: i
   real(r8), allocatable :: wdata(:), wdata2d(:,:)
   logical :: write_attr

      IF (totalvlen <= 0) RETURN

      CALL vector_gather_to_master (vector, vlen, totalvlen, data_address, wdata)

      IF (p_is_master) THEN

         allocate (wdata2d (nlon,nlat))
         wdata2d(:,:) = spval

         DO i = 1, totalvlen
            wdata2d(x_vec(i),y_vec(i)) = wdata(i)
         ENDDO

         IF (present(itime_in_file)) THEN
            CALL ncio_write_serial_time (fileout, varname, itime_in_file, wdata2d, &
               lon_name, lat_name, 'time', DEF_HIST_CompressLevel)
         ELSE
            CALL ncio_write_serial (fileout, varname, wdata2d, &
               lon_name, lat_name, DEF_REST_CompressLevel)
         ENDIF

         IF (present(itime_in_file)) THEN
            write_attr = itime_in_file == 1
         ELSE
            write_attr = .true.
         ENDIF

         IF (write_attr) THEN
            CALL ncio_put_attr (fileout, varname, 'missing_value', spval)
            IF (present(longname)) CALL ncio_put_attr (fileout, varname, 'long_name', longname)
            IF (present(units   )) CALL ncio_put_attr (fileout, varname, 'units',     units   )
         ENDIF

         deallocate (wdata  )
         deallocate (wdata2d)

      ENDIF

   END SUBROUTINE vector_gather_map2grid_and_write

   ! -----
   SUBROUTINE vector_read_and_scatter ( &
         filein, vector, vlen, varname, data_address)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_DataType
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*),       intent(in)    :: filein
   real(r8),  allocatable, intent(inout) :: vector (:)
   integer,                intent(in)    :: vlen
   character(len=*),       intent(in)    :: varname
   type(pointer_int32_1d), intent(in)    :: data_address (0:)

   ! Local variables
   integer :: iwork, ndata, expected_length
   real(r8), allocatable :: rdata(:), rcache(:)

      IF (p_is_master) THEN
         CALL ncio_read_serial (filein, varname, rdata)

         expected_length = 0
         DO iwork = lbound(data_address, 1), ubound(data_address, 1)
            IF (allocated(data_address(iwork)%val)) &
               expected_length = expected_length + size(data_address(iwork)%val)
         ENDDO
         IF (size(rdata) /= expected_length) THEN
            CALL CoLM_stop ('vector_read_and_scatter: restart vector length mismatch')
         ENDIF

         DO iwork = lbound(data_address, 1), ubound(data_address, 1)
            IF (allocated(data_address(iwork)%val)) THEN
               IF (any(data_address(iwork)%val < 1 .or. &
                       data_address(iwork)%val > size(rdata))) THEN
                  CALL CoLM_stop ('vector_read_and_scatter: restart vector address out of range')
               ENDIF
            ENDIF
         ENDDO
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            IF (allocated(data_address(iwork)%val)) THEN

               ndata = size(data_address(iwork)%val)
               allocate(rcache (ndata))
               rcache = rdata(data_address(iwork)%val)

               CALL mpi_send (rcache, ndata, MPI_REAL8, &
                  p_address_worker(iwork), mpi_tag_data, p_comm_glb, p_err)

               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         IF (vlen > 0) THEN
            IF (.not. allocated(vector)) allocate(vector(vlen))
            CALL mpi_recv (vector, vlen, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      IF (.not. allocated(vector)) allocate(vector(vlen))
      vector = rdata(data_address(0)%val)
#endif

      IF (p_is_master) deallocate(rdata)

   END SUBROUTINE vector_read_and_scatter

   ! -----
   SUBROUTINE vector_read_matrix_and_scatter ( &
         filein, matrix, nrow, ncol_local, varname, global_id, ncol_global)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character(len=*),       intent(in)    :: filein
   real(r8),  allocatable, intent(inout) :: matrix (:,:)
   integer,                intent(in)    :: nrow
   integer,                intent(in)    :: ncol_local
   character(len=*),       intent(in)    :: varname
   integer,                intent(in)    :: global_id (:)
   integer,                intent(in)    :: ncol_global

   integer :: iwork, isrc, ndata, ipth, mesg(2)
   integer, allocatable :: icache(:)
   real(r8), allocatable :: rdata(:,:), rcache(:)

      IF (nrow <= 0 .or. ncol_global <= 0) RETURN

      IF (size(global_id) /= ncol_local) THEN
         CALL CoLM_stop ('vector_read_matrix_and_scatter: global_id size mismatch')
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_read_serial (filein, varname, rdata)

         IF (size(rdata, 1) /= nrow .or. size(rdata, 2) /= ncol_global) THEN
            CALL CoLM_stop ('vector_read_matrix_and_scatter: restart matrix shape mismatch')
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, ncol_local/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (ncol_local > 0) THEN
            CALL mpi_send (global_id, ncol_local, MPI_INTEGER, &
               p_address_master, mpi_tag_mesg + 1, p_comm_glb, p_err)

            IF (.not. allocated(matrix)) allocate (matrix(nrow, ncol_local))
            CALL mpi_recv (matrix, nrow*ncol_local, MPI_REAL8, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF

      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               allocate (icache (ndata))
               allocate (rcache (nrow*ndata))

               CALL mpi_recv (icache, ndata, MPI_INTEGER, isrc, &
                  mpi_tag_mesg + 1, p_comm_glb, p_stat, p_err)

               DO ipth = 1, ndata
                  IF (icache(ipth) < 1 .or. icache(ipth) > ncol_global) THEN
                     CALL CoLM_stop ('vector_read_matrix_and_scatter: global_id out of range')
                  ENDIF
                  rcache((ipth-1)*nrow+1:ipth*nrow) = rdata(:, icache(ipth))
               ENDDO

               CALL mpi_send (rcache, nrow*ndata, MPI_REAL8, isrc, &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (icache)
               deallocate (rcache)
            ENDIF
         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      IF (.not. allocated(matrix)) allocate (matrix(nrow, ncol_local))
      DO ipth = 1, ncol_local
         IF (global_id(ipth) < 1 .or. global_id(ipth) > ncol_global) THEN
            CALL CoLM_stop ('vector_read_matrix_and_scatter: global_id out of range')
         ENDIF
         matrix(:, ipth) = rdata(:, global_id(ipth))
      ENDDO
#endif

      IF (p_is_master) deallocate (rdata)

   END SUBROUTINE vector_read_matrix_and_scatter

   ! ==================================================================
   ! IO-group shard writers (DEF_HIST_mode = 'block')
   ! ==================================================================

   !> Collect the per-group counts, displacements and global ids once.
   !!
   !! Called by IO ranks and workers together; the master must not call it.
   !! Built when a history file is created, then reused for every variable,
   !! so no metadata is re-sent per variable.
   SUBROUTINE route_shard_layout_build (layout, nlocal, gid_local)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(route_shard_layout_type), intent(inout) :: layout
   integer, intent(in) :: nlocal
   integer, intent(in) :: gid_local (:)

   integer :: i
   integer, allocatable :: sbuf (:)

      CALL route_shard_layout_free (layout)

      IF (nlocal > 0 .and. size(gid_local) < nlocal) THEN
         CALL CoLM_stop ('route_shard_layout_build: gid_local shorter than nlocal')
      ENDIF

      layout%nlocal = nlocal

#ifdef USEMPI
      IF (p_is_io) THEN
         allocate (layout%cnt (0:p_np_group-1))
         allocate (layout%dsp (0:p_np_group-1))
      ELSE
         ! Unreferenced on non-root ranks, but must be allocated so the
         ! argument is always associated.
         allocate (layout%cnt (1))
         allocate (layout%dsp (1))
      ENDIF

      CALL mpi_gather (nlocal, 1, MPI_INTEGER, layout%cnt, 1, MPI_INTEGER, &
         p_root, p_comm_group, p_err)

      IF (p_is_io) THEN
         layout%dsp(0) = 0
         DO i = 1, p_np_group-1
            layout%dsp(i) = layout%dsp(i-1) + layout%cnt(i-1)
         ENDDO
         layout%ntotal = sum(layout%cnt)
         allocate (layout%gid (max(layout%ntotal,1)))
         layout%gid(:) = 0
      ELSE
         layout%ntotal = 0
         allocate (layout%gid (1))
         layout%gid(:) = 0
      ENDIF

      ! A zero-length worker still has to enter the collective; give it a
      ! valid but unread send buffer rather than an unallocated actual.
      IF (nlocal > 0) THEN
         allocate (sbuf (nlocal)); sbuf = gid_local(1:nlocal)
      ELSE
         allocate (sbuf (1));      sbuf = 0
      ENDIF

      IF (p_is_io) THEN
         CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
            layout%gid, layout%cnt, layout%dsp, MPI_INTEGER, &
            p_root, p_comm_group, p_err)
      ELSE
         CALL mpi_gatherv (sbuf, nlocal, MPI_INTEGER, &
            MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, &
            p_root, p_comm_group, p_err)
      ENDIF

      deallocate (sbuf)
#else
      allocate (layout%cnt (0:0)); layout%cnt(0) = nlocal
      allocate (layout%dsp (0:0)); layout%dsp(0) = 0
      layout%ntotal = nlocal
      allocate (layout%gid (max(nlocal,1))); layout%gid(:) = 0
      IF (nlocal > 0) layout%gid(1:nlocal) = gid_local(1:nlocal)
#endif

      layout%built = .true.

   END SUBROUTINE route_shard_layout_build

   ! -------
   SUBROUTINE route_shard_layout_free (layout)

   IMPLICIT NONE
   type(route_shard_layout_type), intent(inout) :: layout

      IF (allocated(layout%cnt)) deallocate (layout%cnt)
      IF (allocated(layout%dsp)) deallocate (layout%dsp)
      IF (allocated(layout%gid)) deallocate (layout%gid)
      layout%built  = .false.
      layout%nlocal = 0
      layout%ntotal = 0

   END SUBROUTINE route_shard_layout_free

   !> One unit-catchment / reservoir field, gathered to the group's IO rank
   !! and written to that group's shard. No global-master buffer is ever
   !! allocated: the IO rank holds only its own group's share.
   SUBROUTINE route_shard_write_vector (layout, vector, fileshard, varname, &
         dimname, itime_in_file, longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   type(route_shard_layout_type), intent(in) :: layout
   real(r8),         intent(in) :: vector (:)
   character(len=*), intent(in) :: fileshard, varname, dimname

   integer,          intent(in), optional :: itime_in_file
   character(len=*), intent(in), optional :: longname, units

   real(r8), allocatable :: rbuff(:), sbuff(:)
   logical :: write_attr

      IF (.not. layout%built) THEN
         CALL CoLM_stop ('route_shard_write_vector: layout not built')
      ENDIF

      IF (layout%nlocal > 0) THEN
         allocate (sbuff (layout%nlocal)); sbuff = vector(1:layout%nlocal)
      ELSE
         allocate (sbuff (1));             sbuff = spval
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         allocate (rbuff (max(layout%ntotal,1))); rbuff = spval
         CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL8, &
            rbuff, layout%cnt, layout%dsp, MPI_REAL8, &
            p_root, p_comm_group, p_err)
      ELSE
         CALL mpi_gatherv (sbuff, layout%nlocal, MPI_REAL8, &
            MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, &
            p_root, p_comm_group, p_err)
      ENDIF
#else
      allocate (rbuff (max(layout%ntotal,1))); rbuff = spval
      IF (layout%ntotal > 0) rbuff(1:layout%ntotal) = sbuff(1:layout%ntotal)
#endif

      IF (p_is_io) THEN
         ! An IO group that owns nothing still writes its shard, so the
         ! aggregator sees a complete, self-describing set of shard_count
         ! files rather than having to distinguish "empty" from "missing".
         IF (present(itime_in_file)) THEN
            CALL ncio_write_serial_time (fileshard, varname, itime_in_file, &
               rbuff(1:max(layout%ntotal,0)), dimname, 'time', DEF_HIST_CompressLevel)
            write_attr = itime_in_file <= 1
         ELSE
            CALL ncio_write_serial (fileshard, varname, &
               rbuff(1:max(layout%ntotal,0)), dimname, DEF_HIST_CompressLevel)
            write_attr = .true.
         ENDIF

         IF (write_attr) THEN
            CALL ncio_put_attr (fileshard, varname, 'missing_value', spval)
            IF (present(longname)) CALL ncio_put_attr (fileshard, varname, 'long_name', longname)
            IF (present(units   )) CALL ncio_put_attr (fileshard, varname, 'units',     units   )
         ENDIF
      ENDIF

      IF (allocated(sbuff)) deallocate (sbuff)
      IF (allocated(rbuff)) deallocate (rbuff)

   END SUBROUTINE route_shard_write_vector

   !> One bifurcation pathway matrix (nrow x local pathways) per shard.
   !! Column order follows the gathered global-id order recorded in the
   !! layout, which the shard also stores, so the aggregator never relies
   !! on rank order.
   SUBROUTINE route_shard_write_matrix (layout, matrix, nrow, fileshard, &
         varname, rowdim, coldim, itime_in_file, longname, units)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   type(route_shard_layout_type), intent(in) :: layout
   real(r8),         intent(in) :: matrix (:,:)
   integer,          intent(in) :: nrow
   character(len=*), intent(in) :: fileshard, varname, rowdim, coldim

   integer,          intent(in), optional :: itime_in_file
   character(len=*), intent(in), optional :: longname, units

   real(r8), allocatable :: rbuff(:,:), sbuff(:,:)
   integer,  allocatable :: cnt(:), dsp(:)
   integer  :: i
   logical  :: write_attr

      IF (.not. layout%built) THEN
         CALL CoLM_stop ('route_shard_write_matrix: layout not built')
      ENDIF
      IF (nrow <= 0) RETURN

      IF (layout%nlocal > 0) THEN
         allocate (sbuff (nrow, layout%nlocal))
         sbuff = matrix(1:nrow, 1:layout%nlocal)
      ELSE
         allocate (sbuff (nrow, 1)); sbuff = spval
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         ! Element counts scale by the row dimension; the column layout is
         ! exactly the one already collected for the vector case.
         allocate (cnt (0:p_np_group-1))
         allocate (dsp (0:p_np_group-1))
         DO i = 0, p_np_group-1
            cnt(i) = nrow * layout%cnt(i)
            dsp(i) = nrow * layout%dsp(i)
         ENDDO
         allocate (rbuff (nrow, max(layout%ntotal,1))); rbuff = spval
         CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL8, &
            rbuff, cnt, dsp, MPI_REAL8, p_root, p_comm_group, p_err)
         deallocate (cnt, dsp)
      ELSE
         CALL mpi_gatherv (sbuff, nrow*layout%nlocal, MPI_REAL8, &
            MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, &
            p_root, p_comm_group, p_err)
      ENDIF
#else
      allocate (rbuff (nrow, max(layout%ntotal,1))); rbuff = spval
      IF (layout%ntotal > 0) rbuff(:,1:layout%ntotal) = sbuff(:,1:layout%ntotal)
#endif

      IF (p_is_io) THEN
         IF (present(itime_in_file)) THEN
            CALL ncio_write_serial_time (fileshard, varname, itime_in_file, &
               rbuff(1:nrow,1:max(layout%ntotal,0)), rowdim, coldim, 'time', &
               DEF_HIST_CompressLevel)
            write_attr = itime_in_file <= 1
         ELSE
            CALL ncio_write_serial (fileshard, varname, &
               rbuff(1:nrow,1:max(layout%ntotal,0)), rowdim, coldim, &
               DEF_HIST_CompressLevel)
            write_attr = .true.
         ENDIF

         IF (write_attr) THEN
            CALL ncio_put_attr (fileshard, varname, 'missing_value', spval)
            IF (present(longname)) CALL ncio_put_attr (fileshard, varname, 'long_name', longname)
            IF (present(units   )) CALL ncio_put_attr (fileshard, varname, 'units',     units   )
         ENDIF
      ENDIF

      IF (allocated(sbuff)) deallocate (sbuff)
      IF (allocated(rbuff)) deallocate (rbuff)

   END SUBROUTINE route_shard_write_matrix

   !> Shard file name. Purely for locating the file: correctness of the
   !! aggregation is decided by the identity attributes, never by this.
   !!
   !! The segment is part of the name because a restart inside one history
   !! period must not reuse the previous run's shards: if the IO-group count
   !! changed, those shards partition the domain differently and stale ones
   !! would linger. Each continuous run segment writes its own complete set,
   !! and the aggregator merges the segments by time record.
   SUBROUTINE route_shard_filename (file_target, shard_index, fileshard, segment)

   IMPLICIT NONE
   character(len=*),   intent(in)  :: file_target
   integer,            intent(in)  :: shard_index
   character(len=256), intent(out) :: fileshard
   character(len=*),   intent(in), optional :: segment

   integer :: i
   character(len=8)  :: cidx
   character(len=80) :: seg

      write(cidx,'(I5.5)') shard_index
      seg = ''
      IF (present(segment)) THEN
         IF (len_trim(segment) > 0) seg = '_' // trim(segment)
      ENDIF

      i = len_trim(file_target)
      IF (i > 3) THEN
         IF (file_target(i-2:i) == '.nc') THEN
            fileshard = file_target(1:i-3) // trim(seg) // '_shard' // trim(cidx) // '.nc'
            RETURN
         ENDIF
      ENDIF
      fileshard = trim(file_target) // trim(seg) // '_shard' // trim(cidx) // '.nc'

   END SUBROUTINE route_shard_filename

   !> Deterministic description of the output grid, compared by the
   !! aggregator so shards from a different grid can never be merged.
   FUNCTION route_shard_grid_fingerprint (nlon, nlat, lon, lat) RESULT (fp)

   USE MOD_Precision
   IMPLICIT NONE
   integer,  intent(in) :: nlon, nlat
   real(r8), intent(in) :: lon (:), lat (:)
   character(len=256)   :: fp

   character(len=32) :: c1, c2, c3, c4, c5, c6

      write(c1,'(I0)')      nlon
      write(c2,'(I0)')      nlat
      write(c3,'(ES16.9)')  lon(1)
      write(c4,'(ES16.9)')  lon(max(nlon,1))
      write(c5,'(ES16.9)')  lat(1)
      write(c6,'(ES16.9)')  lat(max(nlat,1))

      fp = trim(adjustl(c1)) // 'x' // trim(adjustl(c2)) // ':' // &
           trim(adjustl(c3)) // ':' // trim(adjustl(c4)) // ':' // &
           trim(adjustl(c5)) // ':' // trim(adjustl(c6))

   END FUNCTION route_shard_grid_fingerprint

   !> Write the shard's identity in one open/close cycle.
   !!
   !! ncio_put_attr reopens the file per attribute; nine of those per shard
   !! per history file is pure overhead, and these attributes are what the
   !! aggregator keys on, so they are written together and atomically.
   SUBROUTINE route_shard_write_identity (fileshard, target_basename, case_name, &
         period_key, segment_id, shard_index, shard_count, grid_fingerprint, &
         time_first, time_last, time_count)

   USE MOD_Precision
   USE MOD_SPMD_Task, only: CoLM_stop
   USE netcdf
   IMPLICIT NONE

   character(len=*), intent(in) :: fileshard, target_basename, case_name
   character(len=*), intent(in) :: period_key, segment_id, grid_fingerprint
   integer,          intent(in) :: shard_index, shard_count, time_count
   real(r8),         intent(in) :: time_first, time_last

   integer :: ncid, ierr

      ierr = nf90_open (trim(fileshard), NF90_WRITE, ncid)
      IF (ierr /= NF90_NOERR) CALL CoLM_stop &
         ('route_shard_write_identity: cannot open '//trim(fileshard))

      ierr = nf90_redef (ncid)
      IF (ierr /= NF90_NOERR .and. ierr /= NF90_EINDEFINE) CALL CoLM_stop &
         ('route_shard_write_identity: cannot enter define mode')

      CALL put_int  ('river_hist_shard_schema_version', ROUTE_SHARD_SCHEMA_VERSION)
      CALL put_str  ('target_history_basename', target_basename)
      CALL put_str  ('case_name',               case_name)
      CALL put_str  ('history_period_key',      period_key)
      CALL put_str  ('segment_id',              segment_id)
      CALL put_int  ('shard_index',             shard_index)
      CALL put_int  ('shard_count',             shard_count)
      CALL put_str  ('grid_fingerprint',        grid_fingerprint)
      CALL put_real ('time_first',              time_first)
      CALL put_real ('time_last',               time_last)
      CALL put_int  ('time_count',              time_count)

      ierr = nf90_enddef (ncid)
      IF (ierr /= NF90_NOERR .and. ierr /= NF90_ENOTINDEFINE) CALL CoLM_stop &
         ('route_shard_write_identity: cannot leave define mode')

      ierr = nf90_close (ncid)
      IF (ierr /= NF90_NOERR) CALL CoLM_stop &
         ('route_shard_write_identity: cannot close '//trim(fileshard))

   CONTAINS

      SUBROUTINE put_str (name, value)
      character(len=*), intent(in) :: name, value
      integer :: e
         e = nf90_put_att (ncid, NF90_GLOBAL, trim(name), trim(value))
         IF (e /= NF90_NOERR) CALL CoLM_stop &
            ('route_shard_write_identity: cannot write attribute '//trim(name))
      END SUBROUTINE put_str

      SUBROUTINE put_int (name, value)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: value
      integer :: e
         e = nf90_put_att (ncid, NF90_GLOBAL, trim(name), value)
         IF (e /= NF90_NOERR) CALL CoLM_stop &
            ('route_shard_write_identity: cannot write attribute '//trim(name))
      END SUBROUTINE put_int

      SUBROUTINE put_real (name, value)
      character(len=*), intent(in) :: name
      real(r8),         intent(in) :: value
      integer :: e
         e = nf90_put_att (ncid, NF90_GLOBAL, trim(name), value)
         IF (e /= NF90_NOERR) CALL CoLM_stop &
            ('route_shard_write_identity: cannot write attribute '//trim(name))
      END SUBROUTINE put_real

   END SUBROUTINE route_shard_write_identity

END MODULE MOD_Vector_ReadWrite
