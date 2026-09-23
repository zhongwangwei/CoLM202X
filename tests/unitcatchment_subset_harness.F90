PROGRAM unitcatchment_subset_harness

! Drives mksrfdata/MOD_UnitCatchmentSubset on a network file.
!   harness  network.nc  subset.nc  closure(1|0)  cells.txt
! cells.txt: one runoff-input grid cell "x y" per line (1-based).

   USE MOD_UnitCatchmentSubset
   IMPLICIT NONE

   character(len=512) :: file_in, file_out, file_cells, arg
   integer :: nlon, nlat, nkeep, nsystem, closure, x, y, ios, unit
   logical, allocatable :: touched(:,:)

      CALL get_command_argument (1, file_in)
      CALL get_command_argument (2, file_out)
      CALL get_command_argument (3, arg)
      read (arg, *) closure
      CALL get_command_argument (4, file_cells)

      CALL unitcatchment_grid_size (trim(file_in), nlon, nlat)
      allocate (touched (nlon, nlat));  touched = .false.

      open (newunit = unit, file = trim(file_cells), status = 'old', action = 'read')
      DO
         read (unit, *, iostat = ios) x, y
         IF (ios /= 0) EXIT
         touched(x, y) = .true.
      ENDDO
      close (unit)

      CALL unitcatchment_subset_write (trim(file_in), trim(file_out), nlon, nlat, touched, closure /= 0, &
         nkeep, nsystem)

      write(*,'(A,I0,A,I0)') 'SUBSET_OK kept=', nkeep, ' systems=', nsystem

END PROGRAM unitcatchment_subset_harness
