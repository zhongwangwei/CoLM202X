program netcdf_mpi_probe
   use netcdf
   implicit none
   include 'mpif.h'
   integer :: ierr, rank, nprocs, ncid, dimid, varid, scalarid
   real(kind=8) :: values(2), scalar, total
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
   if (nprocs /= 2) call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
   if (rank == 0) then
      call checked(nf90_create('mpi-probe.nc', ior(nf90_clobber, nf90_netcdf4), ncid))
      call checked(nf90_def_dim(ncid, 'sample', 2, dimid))
      call checked(nf90_def_var(ncid, 'values', nf90_double, [dimid], varid))
      call checked(nf90_def_var(ncid, 'scalar', nf90_double, scalarid))
      call checked(nf90_enddef(ncid))
      call checked(nf90_put_var(ncid, varid, [42.5d0, -1.25d0]))
      call checked(nf90_put_var(ncid, scalarid, 3.5d0))
      call checked(nf90_close(ncid))
   endif
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   call checked(nf90_open('mpi-probe.nc', nf90_nowrite, ncid))
   call checked(nf90_inq_varid(ncid, 'values', varid))
   call checked(nf90_inq_varid(ncid, 'scalar', scalarid))
   call checked(nf90_get_var(ncid, varid, values))
   call checked(nf90_get_var(ncid, scalarid, scalar))
   call checked(nf90_close(ncid))
   if (any(values /= [42.5d0, -1.25d0])) call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
   call MPI_Allreduce(scalar, total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
   if (total /= 7.d0) call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
   if (rank == 0) print *, 'PASS two-rank NetCDF-4 scalar/array I/O and MPI reduction'
   call MPI_Finalize(ierr)
contains
   subroutine checked(status)
      integer, intent(in) :: status
      if (status /= nf90_noerr) then
         print *, nf90_strerror(status)
         call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
      endif
   end subroutine
end program
