#include <define.h>

MODULE MOD_Opt_Baseflow

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), allocatable :: scale_baseflow (:)

   real(r8), allocatable :: ref_zwt   (:)
   real(r8), allocatable :: prcp_year (:)
   real(r8), allocatable :: et_year   (:)
   real(r8), allocatable :: rsur_year (:)
   real(r8), allocatable :: rsub_year (:)
   integer               :: iter_bf_opt
   real(r8), parameter   :: tol_del_zwt = 0.01_r8

CONTAINS

   ! -----
   SUBROUTINE Opt_Baseflow_init ()

   USE MOD_NetCDFVector
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_TimeInvariants, only: patchtype
   USE MOD_LandPatch,           only: numpatch, landpatch
   IMPLICIT NONE

   ! Local Variables
   character(len=256) :: file_restart

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
               allocate (ref_zwt   (numpatch));  ref_zwt  (:) = zwt
               allocate (prcp_year (numpatch));  prcp_year(:) = spval
               allocate (et_year   (numpatch));  et_year  (:) = spval
               allocate (rsur_year (numpatch));  rsur_year(:) = spval
               allocate (rsub_year (numpatch));  rsub_year(:) = spval
            ENDIF
         ENDIF

         iter_bf_opt = 0

      ENDIF

   END SUBROUTINE Opt_Baseflow_init

   ! -----
   SUBROUTINE BaseFlow_Optimize (idate, deltim, is_spinup)

   USE MOD_TimeManager
   USE MOD_NetCDFVector
   USE MOD_RangeCheck
   USE MOD_Vars_TimeVariables,  only: zwt
   USE MOD_Vars_TimeInvariants, only: patchtype
   USE MOD_Vars_1DForcing,      only: forc_prc, forc_prl
   USE MOD_Vars_1DFluxes,       only: fevpa, rsur, rnof
   USE MOD_LandPatch,           only: numpatch, landpatch
   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   logical,  intent(in) :: is_spinup

   ! Local Variables
   real(r8), allocatable :: rbaseflow(:)
   real(r8) :: recharge
   logical  :: do_opt_baseflow
   integer  :: ipatch
   character(len=256) :: file_restart
   character(len=5)   :: strcyc

      IF (.not. DEF_Optimize_Baseflow) RETURN

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN

            CALL add_spv (forc_prc, prcp_year,  deltim)
            CALL add_spv (forc_prl, prcp_year,  deltim)
            CALL add_spv (fevpa,    et_year,    deltim)
            CALL add_spv (rsur,     rsur_year,  deltim)

            allocate (rbaseflow (numpatch))
            rbaseflow = spval
            WHERE ((rsur /= spval) .and. (rnof /= spval))
               rbaseflow = rnof - rsur
            END WHERE

            CALL add_spv (rbaseflow, rsub_year, deltim)

            deallocate (rbaseflow)
         ENDIF
      ENDIF

      CALL check_vector_data ('rsub_year', rsub_year)

      do_opt_baseflow = (is_spinup) .and. isendofyear (idate, deltim)

      IF (do_opt_baseflow) THEN

         iter_bf_opt = iter_bf_opt + 1

         write(strcyc,'(A1,I4.4)') 'c', iter_bf_opt
         IF (p_is_master) THEN
            CALL system('mkdir -p ' // trim(DEF_dir_restart)//'/ParaOpt/'//strcyc)
         ENDIF
#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         file_restart = trim(DEF_dir_restart)//'/ParaOpt/'//strcyc//'/'//trim(DEF_CASE_NAME)//'_baseflow.nc'
         CALL ncio_create_file_vector (file_restart, landpatch)
         CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
         CALL ncio_write_vector (file_restart, 'scale_baseflow', 'patch', landpatch, scale_baseflow, 1)
         CALL ncio_write_vector (file_restart, 'zwt', 'patch', landpatch, zwt, 1)
         CALL ncio_write_vector (file_restart, 'ref_zwt', 'patch', landpatch, ref_zwt, 1)

         IF (p_is_worker) THEN

            DO ipatch = 1, numpatch
               IF (patchtype(ipatch) <= 1) THEN

                  recharge = prcp_year(ipatch)-et_year(ipatch)-rsur_year(ipatch)

                  write(*,*) 'Check ', zwt(ipatch), ref_zwt(ipatch), recharge, rsub_year(ipatch), &
                     prcp_year(ipatch), et_year(ipatch), rsur_year(ipatch)

                  IF ((recharge > 0) .and. (rsub_year(ipatch) > 0)) THEN

                     IF ((zwt(ipatch) > ref_zwt(ipatch)+tol_del_zwt) &
                        .and. (recharge < rsub_year(ipatch))) THEN
                        scale_baseflow(ipatch) = recharge/rsub_year(ipatch) * scale_baseflow(ipatch)
                     ENDIF

                     IF ((zwt(ipatch) < ref_zwt(ipatch)-tol_del_zwt) &
                        .and. (recharge > rsub_year(ipatch))) THEN
                        scale_baseflow(ipatch) = recharge/rsub_year(ipatch) * scale_baseflow(ipatch)
                     ENDIF

                  ELSEIF (prcp_year(ipatch) < et_year(ipatch)) THEN

                     IF (scale_baseflow(ipatch) > 1.e-8) THEN
                        scale_baseflow(ipatch) = scale_baseflow(ipatch) * 0.1
                     ENDIF

                  ENDIF

                  scale_baseflow(ipatch) = max(1.e-8, min(1.e8, scale_baseflow(ipatch)))

               ENDIF
            ENDDO

            IF (numpatch > 0) THEN
               prcp_year(:) = spval
               et_year  (:) = spval
               rsur_year(:) = spval
               rsub_year(:) = spval
            ENDIF
         ENDIF

         file_restart = trim(DEF_dir_restart) // '/ParaOpt/' // trim(DEF_CASE_NAME) //'_baseflow.nc'
         CALL ncio_create_file_vector (file_restart, landpatch)
         CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')
         CALL ncio_write_vector (file_restart, 'scale_baseflow', 'patch', landpatch, scale_baseflow, 1)

      ENDIF

   END SUBROUTINE BaseFlow_Optimize

   ! -----
   SUBROUTINE Opt_Baseflow_final ()

   IMPLICIT NONE

      IF (allocated(scale_baseflow)) deallocate(scale_baseflow)

      IF (DEF_Optimize_Baseflow) THEN
         IF (allocated(ref_zwt  )) deallocate(ref_zwt  )
         IF (allocated(prcp_year)) deallocate(prcp_year)
         IF (allocated(et_year  )) deallocate(et_year  )
         IF (allocated(rsur_year)) deallocate(rsur_year)
         IF (allocated(rsub_year)) deallocate(rsub_year)
      ENDIF

   END SUBROUTINE Opt_Baseflow_final

   !-----------------------------------------------------------------------
   SUBROUTINE add_spv (var, s, dt)

   USE MOD_Precision

   IMPLICIT NONE

   real(r8), intent(in)    :: var(:)
   real(r8), intent(inout) :: s  (:)
   real(r8), intent(in), optional :: dt
   ! Local variables
   integer :: i

      IF (present(dt)) THEN
         DO i = lbound(var,1), ubound(var,1)
            IF (var(i) /= spval) THEN
               IF (s(i) /= spval) THEN
                  s(i) = s(i) + var(i)*dt
               ELSE
                  s(i) = var(i)*dt
               ENDIF
            ENDIF
         ENDDO
      ELSE
         DO i = lbound(var,1), ubound(var,1)
            IF (var(i) /= spval) THEN
               IF (s(i) /= spval) THEN
                  s(i) = s(i) + var(i)
               ELSE
                  s(i) = var(i)
               ENDIF
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE add_spv
END MODULE MOD_Opt_Baseflow
