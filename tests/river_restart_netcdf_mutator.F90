PROGRAM river_restart_netcdf_mutator
   USE MOD_Precision, only: r8
   USE MOD_NetCDFSerial, only: ncio_read_serial, ncio_write_serial
   USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   IMPLICIT NONE

   character(len=64) :: mode
   character(len=1024) :: restart_file, second_file
   real(r8), allocatable :: identity(:,:), previous_depth(:), current_depth(:)
   real(r8), allocatable :: reservoir_identity(:,:)
   real(r8), allocatable :: base_vector(:)
   real(r8), allocatable :: pth_veloc(:,:), pth_momen(:,:), path_signature(:,:)
   real(r8), allocatable :: pth_veloc_2(:,:), pth_momen_2(:,:), path_signature_2(:,:)
   real(r8) :: base_scalar
   integer :: feature_flag

   CALL get_command_argument(1, mode)
   CALL get_command_argument(2, restart_file)
   CALL get_command_argument(3, second_file)
   IF (len_trim(mode) == 0 .or. len_trim(restart_file) == 0) THEN
      ERROR STOP 'usage: river_restart_netcdf_mutator MODE RESTART_FILE'
   ENDIF

   SELECT CASE (trim(mode))
   CASE ('identity')
      CALL ncio_read_serial(trim(restart_file), 'gridriver_ucatch_identity', identity)
      IF (size(identity, 1) /= 4 .or. size(identity, 2) < 1) ERROR STOP 'unexpected identity matrix shape'
      identity(2, 1) = identity(2, 1) + 1._r8
      CALL ncio_write_serial(trim(restart_file), 'gridriver_ucatch_identity', identity)
   CASE ('identity-nan')
      CALL ncio_read_serial(trim(restart_file), 'gridriver_ucatch_identity', identity)
      IF (size(identity, 1) /= 4 .or. size(identity, 2) < 1) ERROR STOP 'unexpected identity matrix shape'
      identity(2, 1) = ieee_value(identity(2, 1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'gridriver_ucatch_identity', identity)
   CASE ('reservoir-identity-swap')
      CALL ncio_read_serial(trim(restart_file), 'gridriver_reservoir_identity', reservoir_identity)
      IF (size(reservoir_identity, 1) /= 2 .or. size(reservoir_identity, 2) < 2) &
         ERROR STOP 'unexpected reservoir identity matrix shape'
      base_scalar = reservoir_identity(2, 1)
      reservoir_identity(2, 1) = reservoir_identity(2, 2)
      reservoir_identity(2, 2) = base_scalar
      CALL ncio_write_serial(trim(restart_file), 'gridriver_reservoir_identity', reservoir_identity)
   CASE ('reservoir-identity-nan')
      CALL ncio_read_serial(trim(restart_file), 'gridriver_reservoir_identity', reservoir_identity)
      IF (size(reservoir_identity, 1) /= 2 .or. size(reservoir_identity, 2) < 1) &
         ERROR STOP 'unexpected reservoir identity matrix shape'
      reservoir_identity(2, 1) = ieee_value(reservoir_identity(2, 1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'gridriver_reservoir_identity', reservoir_identity)
   CASE ('wdsrf-nan')
      CALL ncio_read_serial(trim(restart_file), 'wdsrf_ucat', base_vector)
      base_vector(1) = ieee_value(base_vector(1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'wdsrf_ucat', base_vector)
   CASE ('wdsrf-negative')
      CALL ncio_read_serial(trim(restart_file), 'wdsrf_ucat', base_vector)
      base_vector(1) = -1._r8
      CALL ncio_write_serial(trim(restart_file), 'wdsrf_ucat', base_vector)
   CASE ('velocity-nan')
      CALL ncio_read_serial(trim(restart_file), 'veloc_riv', base_vector)
      base_vector(1) = ieee_value(base_vector(1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'veloc_riv', base_vector)
   CASE ('velocity-out-of-range')
      CALL ncio_read_serial(trim(restart_file), 'veloc_riv', base_vector)
      base_vector(1) = 51._r8
      CALL ncio_write_serial(trim(restart_file), 'veloc_riv', base_vector)
   CASE ('acctime-nan')
      CALL ncio_read_serial(trim(restart_file), 'acctime_rnof', base_scalar)
      base_scalar = ieee_value(base_scalar, ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'acctime_rnof', base_scalar)
   CASE ('acctime-negative')
      CALL ncio_write_serial(trim(restart_file), 'acctime_rnof', -1._r8)
   CASE ('acc-runoff-nan')
      CALL ncio_read_serial(trim(restart_file), 'acc_rnof_uc', base_vector)
      base_vector(1) = ieee_value(base_vector(1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'acc_rnof_uc', base_vector)
   CASE ('volwater-nan')
      CALL ncio_read_serial(trim(restart_file), 'volwater_ucat', base_vector)
      base_vector(1) = ieee_value(base_vector(1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'volwater_ucat', base_vector)
   CASE ('volwater-negative')
      CALL ncio_read_serial(trim(restart_file), 'volwater_ucat', base_vector)
      base_vector(1) = -1._r8
      CALL ncio_write_serial(trim(restart_file), 'volwater_ucat', base_vector)
   CASE ('volresv-nan')
      CALL ncio_read_serial(trim(restart_file), 'volresv', base_vector)
      base_vector(1) = ieee_value(base_vector(1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'volresv', base_vector)
   CASE ('volresv-negative')
      CALL ncio_read_serial(trim(restart_file), 'volresv', base_vector)
      base_vector(1) = -1._r8
      CALL ncio_write_serial(trim(restart_file), 'volresv', base_vector)
   CASE ('manifest-bif-invalid')
      feature_flag = 2
      CALL ncio_write_serial(trim(restart_file), &
         'gridriver_restart_feature_bifurcation', feature_flag)
   CASE ('manifest-levee-invalid')
      feature_flag = -1
      CALL ncio_write_serial(trim(restart_file), &
         'gridriver_restart_feature_levee', feature_flag)
   CASE ('previous-depth-nan')
      CALL ncio_read_serial(trim(restart_file), 'wdsrf_ucat_prev', previous_depth)
      IF (size(previous_depth) < 1) ERROR STOP 'unexpected previous-depth vector shape'
      previous_depth(1) = ieee_value(previous_depth(1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'wdsrf_ucat_prev', previous_depth)
   CASE ('previous-depth-negative')
      CALL ncio_read_serial(trim(restart_file), 'wdsrf_ucat_prev', previous_depth)
      IF (size(previous_depth) < 1) ERROR STOP 'unexpected previous-depth vector shape'
      previous_depth(1) = -1._r8
      CALL ncio_write_serial(trim(restart_file), 'wdsrf_ucat_prev', previous_depth)
   CASE ('bif-momentum-nan')
      CALL ncio_read_serial(trim(restart_file), 'pth_momen', pth_momen)
      IF (size(pth_momen) < 1) ERROR STOP 'unexpected BIF momentum matrix shape'
      pth_momen(1,1) = ieee_value(pth_momen(1,1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'pth_momen', pth_momen)
   CASE ('bif-signature-nan')
      CALL ncio_read_serial(trim(restart_file), 'bif_path_signature', path_signature)
      IF (size(path_signature) < 1) ERROR STOP 'unexpected BIF signature matrix shape'
      path_signature(1,1) = ieee_value(path_signature(1,1), ieee_quiet_nan)
      CALL ncio_write_serial(trim(restart_file), 'bif_path_signature', path_signature)
   CASE ('check-bif-nonzero')
      CALL ncio_read_serial(trim(restart_file), 'pth_veloc', pth_veloc)
      CALL ncio_read_serial(trim(restart_file), 'pth_momen', pth_momen)
      IF (size(pth_veloc) < 1 .or. size(pth_momen) < 1) ERROR STOP 'empty BIF restart state'
      IF (any(.not. ieee_is_finite(pth_veloc))) ERROR STOP 'non-finite BIF velocity state'
      IF (any(.not. ieee_is_finite(pth_momen))) ERROR STOP 'non-finite BIF momentum state'
      IF (maxval(abs(pth_veloc)) <= 0._r8) ERROR STOP 'BIF velocity state is all zero'
      IF (maxval(abs(pth_momen)) <= 0._r8) ERROR STOP 'BIF momentum state is all zero'
   CASE ('check-bif-cold')
      CALL ncio_read_serial(trim(restart_file), 'wdsrf_ucat', current_depth)
      CALL ncio_read_serial(trim(restart_file), 'wdsrf_ucat_prev', previous_depth)
      CALL ncio_read_serial(trim(restart_file), 'pth_veloc', pth_veloc)
      CALL ncio_read_serial(trim(restart_file), 'pth_momen', pth_momen)
      IF (any(.not. ieee_is_finite(current_depth)) .or. &
          any(.not. ieee_is_finite(previous_depth))) ERROR STOP 'non-finite BIF depth state'
      IF (size(current_depth) /= size(previous_depth)) ERROR STOP 'BIF depth state shape mismatch'
      IF (any(previous_depth /= current_depth)) ERROR STOP 'BIF previous depth was not cold-started atomically'
      IF (any(.not. ieee_is_finite(pth_veloc)) .or. &
          any(.not. ieee_is_finite(pth_momen))) ERROR STOP 'non-finite cold-start BIF state'
      IF (any(pth_veloc /= 0._r8) .or. any(pth_momen /= 0._r8)) &
         ERROR STOP 'BIF pathway state was not cold-started at zero'
   CASE ('compare-bif-state')
      IF (len_trim(second_file) == 0) ERROR STOP 'compare-bif-state requires a second file'
      CALL ncio_read_serial(trim(restart_file), 'wdsrf_ucat_prev', previous_depth)
      CALL ncio_read_serial(trim(second_file), 'wdsrf_ucat_prev', current_depth)
      CALL ncio_read_serial(trim(restart_file), 'pth_veloc', pth_veloc)
      CALL ncio_read_serial(trim(second_file), 'pth_veloc', pth_veloc_2)
      CALL ncio_read_serial(trim(restart_file), 'pth_momen', pth_momen)
      CALL ncio_read_serial(trim(second_file), 'pth_momen', pth_momen_2)
      CALL ncio_read_serial(trim(restart_file), 'bif_path_signature', path_signature)
      CALL ncio_read_serial(trim(second_file), 'bif_path_signature', path_signature_2)
      IF (size(previous_depth) /= size(current_depth)) ERROR STOP 'previous-depth shape changed after repartition'
      IF (any(previous_depth /= current_depth)) ERROR STOP 'previous-depth state changed after repartition'
      CALL compare_matrix('pth_veloc', pth_veloc, pth_veloc_2)
      CALL compare_matrix('pth_momen', pth_momen, pth_momen_2)
      CALL compare_matrix('bif_path_signature', path_signature, path_signature_2)
   CASE DEFAULT
      ERROR STOP 'unknown mutator mode'
   END SELECT

CONTAINS

   SUBROUTINE compare_matrix(name, first, second)
      character(len=*), intent(in) :: name
      real(r8), intent(in) :: first(:,:), second(:,:)

      IF (size(first,1) /= size(second,1) .or. size(first,2) /= size(second,2)) THEN
         write(*,'(A,1X,A)') 'matrix shape changed after repartition:', trim(name)
         ERROR STOP 1
      ENDIF
      IF (any(.not. ieee_is_finite(first)) .or. any(.not. ieee_is_finite(second))) THEN
         write(*,'(A,1X,A)') 'non-finite matrix in repartition comparison:', trim(name)
         ERROR STOP 1
      ENDIF
      IF (any(first /= second)) THEN
         write(*,'(A,1X,A)') 'matrix changed after repartition:', trim(name)
         ERROR STOP 1
      ENDIF
   END SUBROUTINE compare_matrix
END PROGRAM river_restart_netcdf_mutator
