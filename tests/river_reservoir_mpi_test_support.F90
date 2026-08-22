MODULE MOD_DataType
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
      minutes_since_1900 = (year - 1900) * 365 * 1440 + &
         (month - 1) * 30 * 1440 + (day - 1) * 1440
   END FUNCTION minutes_since_1900
END MODULE MOD_TimeManager

MODULE MOD_Namelist
   IMPLICIT NONE
   integer, parameter :: DEF_Reservoir_Method = 1
   character(len=1024) :: DEF_ReservoirPara_file = ''
END MODULE MOD_Namelist

MODULE MOD_Vars_Global
   IMPLICIT NONE
   integer, parameter :: nl_soil = 10, maxsnl = -5, nl_lake = 10, nvegwcs = 5
END MODULE MOD_Vars_Global

MODULE MOD_Grid_RiverLakeNetwork
   IMPLICIT NONE
   integer :: numucat = 0
   integer, allocatable :: ucat_ucid(:), lake_type(:)
END MODULE MOD_Grid_RiverLakeNetwork
