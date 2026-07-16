MODULE MOD_DataType
   USE MOD_Precision, only: r8
   IMPLICIT NONE

   type :: pointer_int32_1d
      integer, allocatable :: val(:)
   END type pointer_int32_1d
END MODULE MOD_DataType

MODULE MOD_TimeManager
   IMPLICIT NONE
CONTAINS
   integer FUNCTION minutes_since_1900(year, month, day)
      integer, intent(in) :: year, month, day
      minutes_since_1900 = (year - 1900) * 365 * 1440 + (month - 1) * 30 * 1440 + (day - 1) * 1440
   END FUNCTION minutes_since_1900
END MODULE MOD_TimeManager

MODULE MOD_Namelist
   IMPLICIT NONE
   integer, parameter :: DEF_REST_CompressLevel = 0
   integer, parameter :: DEF_HIST_CompressLevel = 0
   integer, parameter :: DEF_Reservoir_Method = 1
   logical, parameter :: DEF_USE_LEVEE = .false.
   logical, parameter :: DEF_USE_BIFURCATION = .false.
END MODULE MOD_Namelist

MODULE MOD_Vars_Global
   USE MOD_Precision, only: r8
   IMPLICIT NONE
   real(r8), parameter :: spval = -1.e36_r8
   integer, parameter :: nl_soil = 10
   integer, parameter :: maxsnl = -5
   integer, parameter :: nl_lake = 10
   integer, parameter :: nvegwcs = 5
END MODULE MOD_Vars_Global

MODULE MOD_Grid_RiverLakeNetwork
   USE MOD_DataType, only: pointer_int32_1d
   IMPLICIT NONE
   integer :: numucat = 0
   integer :: totalnumucat = 0
   integer :: totalnpthout = 0
   integer :: npthlev_bif = 0
   integer, allocatable :: ucat_ucid(:)
   integer, allocatable :: x_ucat(:)
   integer, allocatable :: y_ucat(:)
   integer, allocatable :: ucat_next(:)
   type(pointer_int32_1d), allocatable :: ucat_data_address(:)
END MODULE MOD_Grid_RiverLakeNetwork

MODULE MOD_Grid_Reservoir
   USE MOD_DataType, only: pointer_int32_1d
   IMPLICIT NONE
   integer :: numresv = 0
   integer :: totalnumresv = 0
   integer, allocatable :: resv_global_id(:)
   integer, allocatable :: resv_ucid(:)
   type(pointer_int32_1d), allocatable :: resv_data_address(:)
END MODULE MOD_Grid_Reservoir

MODULE MOD_Grid_RiverLakeLevee
   USE MOD_Precision, only: r8
   IMPLICIT NONE
   real(r8), allocatable :: levsto(:)
END MODULE MOD_Grid_RiverLakeLevee

MODULE MOD_Grid_RiverLakeBifurcation
   IMPLICIT NONE
CONTAINS
   SUBROUTINE write_bifurcation_restart(file_restart)
      character(len=*), intent(in) :: file_restart
   END SUBROUTINE write_bifurcation_restart
END MODULE MOD_Grid_RiverLakeBifurcation
