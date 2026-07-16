PROGRAM river_levee_restart_netcdf_mutator
   USE MOD_Precision, only: r8
   USE MOD_NetCDFSerial, only: ncio_read_serial, ncio_write_serial
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   IMPLICIT NONE

   character(len=64) :: mode
   character(len=1024) :: restart_file
   real(r8), allocatable :: levsto(:)

   CALL get_command_argument(1, mode)
   CALL get_command_argument(2, restart_file)
   IF (len_trim(mode) == 0 .or. len_trim(restart_file) == 0) THEN
      ERROR STOP 'usage: river_levee_restart_netcdf_mutator MODE RESTART_FILE'
   ENDIF

   CALL ncio_read_serial(trim(restart_file), 'levsto', levsto)
   IF (size(levsto) < 1) ERROR STOP 'unexpected levsto vector shape'

   SELECT CASE (trim(mode))
   CASE ('nan')
      levsto(1) = ieee_value(levsto(1), ieee_quiet_nan)
   CASE ('negative')
      levsto(1) = -1._r8
   CASE DEFAULT
      ERROR STOP 'unknown levee mutator mode'
   END SELECT
   CALL ncio_write_serial(trim(restart_file), 'levsto', levsto)

END PROGRAM river_levee_restart_netcdf_mutator
