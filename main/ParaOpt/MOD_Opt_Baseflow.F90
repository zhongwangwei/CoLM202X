#include <define.h>

MODULE MOD_Opt_Baseflow

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   IMPLICIT NONE

   real(r8), allocatable :: scale_baseflow (:)

   real(r8), allocatable :: scale_bf_left  (:)
   real(r8), allocatable :: scale_bf_right (:)

   integer :: ref_date(3)
   integer :: ref_lc_year
   real(r8), allocatable :: ref_zwt        (:)

   logical,  allocatable :: mask_bf_opt    (:)
   integer,  allocatable :: flag_bf_opt    (:)

   real(r8), parameter   :: tol_del_zwt = 0.01_r8
   integer               :: iter_bf_opt

CONTAINS

   ! -----
   SUBROUTINE Opt_Baseflow_init (ref_date_in, ref_lc_year_in)

   USE MOD_NetCDFVector
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_TimeInvariants, only: patchtype
   USE MOD_LandPatch,           only: numpatch, landpatch
   IMPLICIT NONE

   integer :: ref_date_in(3), ref_lc_year_in

   ! Local Variables
   character(len=256) :: file_restart


      ref_date    = ref_date_in
      ref_lc_year = ref_lc_year_in

      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(DEF_dir_restart)//'/ParaOpt')
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      file_restart = trim(DEF_dir_restart) // '/ParaOpt/' // trim(DEF_CASE_NAME) //'_baseflow.nc'
      CALL ncio_read_vector (file_restart, 'scale_baseflow', landpatch, scale_baseflow, defval = 1.)

      IF (DEF_Optimize_Baseflow) THEN

         IF (p_is_worker) THEN
            IF (numpatch > 0) THEN

               allocate (scale_bf_left  (numpatch));  scale_bf_left (:) = -1.
               allocate (scale_bf_right (numpatch));  scale_bf_right(:) = -1.
               allocate (ref_zwt        (numpatch));  ref_zwt       (:) = zwt

               allocate (mask_bf_opt (numpatch))
               mask_bf_opt (:) = patchtype <= 1

               allocate (flag_bf_opt (numpatch))
               WHERE (mask_bf_opt)
                  flag_bf_opt = 1
               ELSE WHERE
                  flag_bf_opt = 0
               END WHERE

            ENDIF
         ENDIF

         iter_bf_opt = 0
         CALL READ_TimeVariables (ref_date, ref_lc_year, DEF_CASE_NAME, DEF_dir_restart)

      ENDIF

   END SUBROUTINE Opt_Baseflow_init

   ! -----
   SUBROUTINE BaseFlow_Optimize (idate, deltim, end_of_spinup)

   USE MOD_TimeManager
   USE MOD_NetCDFVector
   USE MOD_Vars_TimeVariables
   USE MOD_LandPatch, only: numpatch, landpatch
   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   logical,  intent(in) :: end_of_spinup

   ! Local Variables
   logical :: do_opt_baseflow
   integer :: ipatch
   character(len=256) :: file_restart
   character(len=5)   :: strcyc

      IF (.not. DEF_Optimize_Baseflow) RETURN

      do_opt_baseflow = (.not. end_of_spinup) .and. isendofyear (idate, deltim)

      IF (do_opt_baseflow) THEN

         iter_bf_opt = iter_bf_opt + 1

         ! IF (p_is_worker) THEN
         !    DO ipatch = 1, numpatch
         !       IF (mask_bf_opt(ipatch)) THEN
         !          write(*,'(A,5Es10.2)') 'check scale_baseflow:', zwt(ipatch), ref_zwt(ipatch), &
         !             scale_bf_left(ipatch), scale_baseflow(ipatch), scale_bf_right(ipatch)
         !       ENDIF
         !    ENDDO
         ! ENDIF

         file_restart = trim(DEF_dir_restart) // '/ParaOpt/' // trim(DEF_CASE_NAME) //'_baseflow.nc'
         CALL ncio_create_file_vector (file_restart, landpatch)
         CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
         CALL ncio_write_vector (file_restart, 'scale_baseflow', 'patch', landpatch, scale_baseflow, 1)

         write(strcyc,'(A1,I4.4)') 'c', iter_bf_opt
         file_restart = trim(DEF_dir_restart)//'/ParaOpt/'//trim(DEF_CASE_NAME)//'_baseflow_'//strcyc//'.nc'
         CALL ncio_create_file_vector (file_restart, landpatch)
         CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
         CALL ncio_write_vector (file_restart, 'scale_baseflow', 'patch', landpatch, scale_baseflow, 1)
         CALL ncio_write_vector (file_restart, 'flag_optimize',  'patch', landpatch, flag_bf_opt,    1)
         CALL ncio_write_vector (file_restart, 'zwt', 'patch', landpatch, zwt, 1)
         CALL ncio_write_vector (file_restart, 'ref_zwt', 'patch', landpatch, ref_zwt, 1)

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               IF (mask_bf_opt(ipatch)) THEN
                  IF (zwt(ipatch) > ref_zwt(ipatch)+tol_del_zwt) THEN

                     scale_bf_right(ipatch) = scale_baseflow(ipatch)
                     IF (scale_bf_left(ipatch) > 0.) THEN
                        scale_baseflow(ipatch) = sqrt((scale_bf_left(ipatch) * scale_baseflow(ipatch)))
                     ELSE
                        scale_baseflow(ipatch) = scale_baseflow(ipatch) / 2.
                     ENDIF

                  ELSEIF (zwt(ipatch) < ref_zwt(ipatch)-tol_del_zwt) THEN

                     scale_bf_left(ipatch) = scale_baseflow(ipatch)
                     IF (scale_bf_right(ipatch) > 0.) THEN
                        scale_baseflow(ipatch) = sqrt((scale_bf_right(ipatch) * scale_baseflow(ipatch)))
                     ELSE
                        scale_baseflow(ipatch) = scale_baseflow(ipatch) * 2.
                     ENDIF

                  ELSE

                     mask_bf_opt(ipatch) = .false.
                     flag_bf_opt(ipatch) = 2

                  ENDIF
               ENDIF
            ENDDO
         ENDIF

         CALL READ_TimeVariables (ref_date, ref_lc_year, DEF_CASE_NAME, DEF_dir_restart)

      ENDIF

   END SUBROUTINE BaseFlow_Optimize

   ! -----
   SUBROUTINE Opt_Baseflow_final ()

   IMPLICIT NONE

      IF (allocated(scale_baseflow)) deallocate(scale_baseflow)

      IF (DEF_Optimize_Baseflow) THEN
         IF (allocated(scale_bf_left )) deallocate(scale_bf_left )
         IF (allocated(scale_bf_right)) deallocate(scale_bf_right)
         IF (allocated(ref_zwt       )) deallocate(ref_zwt       )
         IF (allocated(mask_bf_opt   )) deallocate(mask_bf_opt   )
         IF (allocated(flag_bf_opt   )) deallocate(flag_bf_opt   )
      ENDIF

   END SUBROUTINE Opt_Baseflow_final

END MODULE MOD_Opt_Baseflow
