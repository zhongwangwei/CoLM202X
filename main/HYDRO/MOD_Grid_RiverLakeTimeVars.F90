#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeTimeVars
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Time Variables in gridded hydrological processes.
!
! Created by Shupeng Zhang, Oct 2025
!-------------------------------------------------------------------------------------

   USE MOD_Precision
#ifdef GridRiverLakeSediment
   USE MOD_Grid_RiverLakeSediment, only: nsed, totlyrnum, sedcon, layer, seddep
#endif
   IMPLICIT NONE

   ! -- state variables --
   real(r8), allocatable :: wdsrf_ucat (:) ! river or lake water depth [m]
   real(r8), allocatable :: veloc_riv  (:) ! river velocity            [m/s]
   real(r8), allocatable :: momen_riv  (:) ! unit river momentum       [m^2/s]
   real(r8), allocatable :: volresv    (:) ! reservoir water volume    [m^3]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_GridRiverLakeTimeVars
   PUBLIC :: deallocate_GridRiverLakeTimeVars

   PUBLIC :: read_GridRiverLakeTimeVars
   PUBLIC :: write_GridRiverLakeTimeVars

CONTAINS

#ifdef GridRiverLakeSediment
   !-------------------------------------------------------------------------------------
   FUNCTION int2char(i) RESULT(str)
   ! Convert integer to character string
   !-------------------------------------------------------------------------------------
   IMPLICIT NONE
   integer, intent(in) :: i
   character(len=10) :: str
      write(str, '(I10)') i
      str = adjustl(str)
   END FUNCTION int2char
#endif

   SUBROUTINE allocate_GridRiverLakeTimeVars

   USE MOD_SPMD_Task
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir,        only: numresv
   IMPLICIT NONE

      IF (p_is_worker) THEN

         IF (numucat > 0)  allocate (wdsrf_ucat (numucat))
         IF (numucat > 0)  allocate (veloc_riv  (numucat))
         IF (numucat > 0)  allocate (momen_riv  (numucat))
         IF (numresv > 0)  allocate (volresv    (numresv))

      ENDIF

   END SUBROUTINE allocate_GridRiverLakeTimeVars


   SUBROUTINE READ_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: numresv, resv_data_address, totalnumresv
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

#ifdef GridRiverLakeSediment
   integer :: ised, ilyr
#endif

      CALL vector_read_and_scatter (file_restart, wdsrf_ucat, numucat, 'wdsrf_ucat', ucat_data_address)
      CALL vector_read_and_scatter (file_restart, veloc_riv,  numucat, 'veloc_riv',  ucat_data_address)

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN
            CALL vector_read_and_scatter (file_restart, volresv, numresv, 'volresv', resv_data_address)
         ENDIF
      ENDIF

#ifdef GridRiverLakeSediment
      IF (DEF_USE_SEDIMENT) THEN
         ! Read sediment concentration for each size class
         DO ised = 1, nsed
            CALL vector_read_and_scatter (file_restart, sedcon(ised,:), numucat, &
               'sedcon_' // trim(int2char(ised)), ucat_data_address)
         ENDDO

         ! Read active layer storage for each size class
         DO ised = 1, nsed
            CALL vector_read_and_scatter (file_restart, layer(ised,:), numucat, &
               'layer_' // trim(int2char(ised)), ucat_data_address)
         ENDDO

         ! Read deposition layer storage for each size class and layer
         DO ised = 1, nsed
            DO ilyr = 1, totlyrnum
               CALL vector_read_and_scatter (file_restart, seddep(ised,ilyr,:), numucat, &
                  'seddep_' // trim(int2char(ised)) // '_' // trim(int2char(ilyr)), ucat_data_address)
            ENDDO
         ENDDO
      ENDIF
#endif

   END SUBROUTINE READ_GridRiverLakeTimeVars


   SUBROUTINE WRITE_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: numresv, totalnumresv, resv_data_address
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

#ifdef GridRiverLakeSediment
   integer :: ised, ilyr
#endif

      IF (p_is_master) THEN
         CALL ncio_create_file (trim(file_restart))
         CALL ncio_define_dimension(file_restart, 'ucatch', totalnumucat)
#ifdef GridRiverLakeSediment
         IF (DEF_USE_SEDIMENT) THEN
            CALL ncio_define_dimension(file_restart, 'nsed', nsed)
            CALL ncio_define_dimension(file_restart, 'totlyrnum', totlyrnum)
         ENDIF
#endif
      ENDIF

      CALL vector_gather_and_write (&
         wdsrf_ucat, numucat, totalnumucat, ucat_data_address, file_restart, 'wdsrf_ucat', 'ucatch')

      CALL vector_gather_and_write (&
         veloc_riv, numucat, totalnumucat, ucat_data_address, file_restart, 'veloc_riv', 'ucatch')

      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN

            IF (p_is_master) CALL ncio_define_dimension(file_restart, 'reservoir', totalnumresv)

            CALL vector_gather_and_write (&
               volresv, numresv, totalnumresv, resv_data_address, file_restart, 'volresv', 'reservoir')
         ENDIF
      ENDIF

#ifdef GridRiverLakeSediment
      IF (DEF_USE_SEDIMENT) THEN
         ! Write sediment concentration for each size class
         DO ised = 1, nsed
            CALL vector_gather_and_write (&
               sedcon(ised,:), numucat, totalnumucat, ucat_data_address, file_restart, &
               'sedcon_' // trim(int2char(ised)), 'ucatch')
         ENDDO

         ! Write active layer storage for each size class
         DO ised = 1, nsed
            CALL vector_gather_and_write (&
               layer(ised,:), numucat, totalnumucat, ucat_data_address, file_restart, &
               'layer_' // trim(int2char(ised)), 'ucatch')
         ENDDO

         ! Write deposition layer storage for each size class and layer
         DO ised = 1, nsed
            DO ilyr = 1, totlyrnum
               CALL vector_gather_and_write (&
                  seddep(ised,ilyr,:), numucat, totalnumucat, ucat_data_address, file_restart, &
                  'seddep_' // trim(int2char(ised)) // '_' // trim(int2char(ilyr)), 'ucatch')
            ENDDO
         ENDDO
      ENDIF
#endif

   END SUBROUTINE WRITE_GridRiverLakeTimeVars

   SUBROUTINE deallocate_GridRiverLakeTimeVars

   IMPLICIT NONE

      IF (allocated (wdsrf_ucat)) deallocate (wdsrf_ucat)
      IF (allocated (veloc_riv )) deallocate (veloc_riv )
      IF (allocated (momen_riv )) deallocate (momen_riv )
      IF (allocated (volresv   )) deallocate (volresv   )

   END SUBROUTINE deallocate_GridRiverLakeTimeVars

END MODULE MOD_Grid_RiverLakeTimeVars
#endif
