#include <define.h>

MODULE MOD_CheckEquilibrium

!----------------------------------------------------------------------------
! !DESCRIPTION:
!
!   Check equilibrium state.
!
!  Created by Shupeng Zhang, 10/2024
!----------------------------------------------------------------------------

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Grid
   USE netcdf
   USE MOD_NetCDFSerial
   USE MOD_SpatialMapping
   USE MOD_Vars_Global, only: spval

   integer, parameter :: windowsize = 10

   ! ----- Variables -----
   integer :: nyearcheck

   integer :: timestrlen
   character(len=24) :: timeform

   real(r8), allocatable :: tws_last (:)
   real(r8), allocatable :: tws_this (:)
   real(r8), allocatable :: prcp_acc (:)

   real(r8), allocatable :: tws_prev (:,:)
   real(r8), allocatable :: prcp_prev(:,:)

   real(r8) :: num_totalck
   real(r8) :: ave_pct_dtws
   character(len=256) :: spinup_warning

#ifndef SinglePoint
   type(grid_type)            :: gridcheck
   type(grid_concat_type)     :: gcheck_concat
   type(spatial_mapping_type) :: map_check

   integer :: check_data_id = 0
#endif

   PUBLIC :: CheckEqb_init
   PUBLIC :: CheckEquilibrium
   PUBLIC :: CheckEquilibrium_EndOfSpinup
   PUBLIC :: CheckEqb_final

CONTAINS

   !-----------------------------------------------------------------------
   SUBROUTINE CheckEqb_init (n_spinupcycle)

   USE MOD_Vars_Global,        only: nl_soil
   USE MOD_Forcing,            only: gforc
   USE MOD_LandPatch,          only: numpatch, landpatch
   USE MOD_Vars_TimeVariables, only: wdsrf, ldew, scv, wetwat, wliq_soisno, wice_soisno, wa
   IMPLICIT NONE

   integer, intent(in) :: n_spinupcycle

   ! Local Variable
   integer :: ilev


      IF (.not. DEF_CheckEquilibrium) RETURN

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN

            allocate (tws_last (numpatch));   tws_last(:) = spval;
            allocate (tws_this (numpatch));   tws_this(:) = spval;
            allocate (prcp_acc (numpatch));   prcp_acc(:) = spval;

            allocate (tws_prev  (0:windowsize,numpatch));   tws_prev (:,:) = spval;
            allocate (prcp_prev (windowsize,  numpatch));   prcp_prev(:,:) = spval;

         ENDIF
      ENDIF

      IF (n_spinupcycle >= 10000) THEN
         timestrlen = 16
         timeform = "('spinup',I5.5,'-',I4.4)"
      ELSEIF (n_spinupcycle >= 1000) THEN
         timestrlen = 15
         timeform = "('spinup',I4.4,'-',I4.4)"
      ELSEIF (n_spinupcycle >= 100) THEN
         timestrlen = 14
         timeform = "('spinup',I3.3,'-',I4.4)"
      ELSEIF (n_spinupcycle >= 10) THEN
         timestrlen = 13
         timeform = "('spinup',I2.2,'-',I4.4)"
      ELSE
         timestrlen = 12
         timeform = "('spinup',I1.1,'-',I4.4)"
      ENDIF

      spinup_warning = ''

#ifndef SinglePoint
      ! grid
      CALL gridcheck%define_by_copy (gforc)
      ! grid info for output
      CALL gcheck_concat%set (gridcheck)
      ! mapping from patch to grid
      CALL map_check%build_arealweighted (gridcheck, landpatch)
#endif

      IF ((p_is_worker) .and. (numpatch > 0)) THEN
         tws_last = wdsrf                                 ! 1. surface water
         CALL add_spv (ldew, tws_last)                    ! 2. water on foliage
         CALL add_spv (scv , tws_last)                    ! 3. snow cover water equivalent
         IF (DEF_USE_VariablySaturatedFlow) THEN
            CALL add_spv (wetwat, tws_last)               ! 4. water in wetland
         ENDIF
         DO ilev = 1, nl_soil
            CALL add_spv (wliq_soisno(ilev,:), tws_last)  ! 5. liquid water in soil
            CALL add_spv (wice_soisno(ilev,:), tws_last)  ! 6. ice in soil
         ENDDO
         CALL add_spv (wa, tws_last)                      ! 7. water in aquifer
      ENDIF

      nyearcheck = 0

   END SUBROUTINE CheckEqb_init

   !-----------------------------------------------------------------------
   SUBROUTINE CheckEqb_final ()

   IMPLICIT NONE

      IF (.not. DEF_CheckEquilibrium) RETURN

      IF (allocated(tws_last)) deallocate(tws_last)
      IF (allocated(tws_this)) deallocate(tws_this)
      IF (allocated(prcp_acc)) deallocate(prcp_acc)

      IF (allocated(tws_prev )) deallocate(tws_prev )
      IF (allocated(prcp_prev)) deallocate(prcp_prev)

   END SUBROUTINE CheckEqb_final

   !-----------------------------------------------------------------------
   SUBROUTINE CheckEquilibrium (idate, deltim, i_spinupcycle, end_of_spinup, dir_out, casename)

   USE MOD_Precision
   USE MOD_TimeManager
   USE MOD_DataType
   USE MOD_LandPatch,          only: numpatch
   USE MOD_Vars_Global,        only: nl_soil
   USE MOD_Vars_1DForcing,     only: forc_prc, forc_prl
   USE MOD_Vars_TimeVariables, only: wdsrf, ldew, scv, wetwat, wliq_soisno, wice_soisno, wa

   IMPLICIT NONE

   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim
   integer,  intent(in) :: i_spinupcycle
   logical,  intent(in) :: end_of_spinup

   character(len=*), intent(in) :: dir_out
   character(len=*), intent(in) :: casename

   ! Local variables
   logical :: docheck
   integer :: ilev
   character(len=256) :: filename, timestr
   integer :: ncid, time_id, str_id, varid

   real(r8), allocatable     :: pct_dtws (:)
   logical,  allocatable     :: filter   (:)
   type(block_data_real8_2d) :: sumarea


      IF (.not. DEF_CheckEquilibrium) RETURN

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            CALL add_spv (forc_prc, prcp_acc, deltim)
            CALL add_spv (forc_prl, prcp_acc, deltim)
         ENDIF
      ENDIF

      docheck = isendofyear (idate, deltim)

      IF (docheck) THEN

         IF ((p_is_worker) .and. (numpatch > 0)) THEN
            tws_this = wdsrf                                 ! 1. surface water
            CALL add_spv (ldew, tws_this)                    ! 2. water on foliage
            CALL add_spv (scv , tws_this)                    ! 3. snow cover water equivalent
            IF (DEF_USE_VariablySaturatedFlow) THEN
               CALL add_spv (wetwat, tws_this)               ! 4. water in wetland
            ENDIF
            DO ilev = 1, nl_soil
               CALL add_spv (wliq_soisno(ilev,:), tws_this)  ! 5. liquid water in soil
               CALL add_spv (wice_soisno(ilev,:), tws_this)  ! 6. ice in soil
            ENDDO
            CALL add_spv (wa, tws_this)                      ! 7. water in aquifer
         ENDIF

         nyearcheck = nyearcheck + 1

         IF (nyearcheck >= 1) THEN

            IF ((p_is_worker) .and. (numpatch > 0)) THEN

               allocate (filter (numpatch))
               filter(:) = (tws_last /= spval) .and. (tws_this /= spval) .and. (prcp_acc > 0.)

               allocate (pct_dtws (numpatch))
               WHERE (filter)
                  pct_dtws = (tws_this - tws_last) / prcp_acc
               ELSEWHERE
                  pct_dtws = spval
               END WHERE

            ENDIF

            IF (p_is_master) THEN

               filename = trim(dir_out) // '/' // trim(casename) //'_check_equilibrium.nc'

               IF (nyearcheck == 1) THEN

                  CALL ncio_create_file (trim(filename))

                  CALL ncio_define_dimension(filename, 'iyear',   0)
                  CALL ncio_define_dimension(filename, 'timestr', timestrlen)

                  CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
                  CALL nccheck( nf90_inq_dimid(ncid, 'iyear',   time_id) )
                  CALL nccheck( nf90_inq_dimid(ncid, 'timestr', str_id ) )
                  CALL nccheck( nf90_redef(ncid) )
                  CALL nccheck( nf90_def_var(ncid, 'iyear', NF90_CHAR, (/str_id,time_id/), varid) )
                  CALL nccheck( nf90_put_att(ncid, varid, 'long_name', 'iyear in all spinup cycles') )
                  CALL nccheck( nf90_enddef(ncid) )
                  CALL nccheck( nf90_close(ncid) )

#ifndef SinglePoint
                  CALL ncio_define_dimension(filename, 'lat' , gcheck_concat%ginfo%nlat)
                  CALL ncio_define_dimension(filename, 'lon' , gcheck_concat%ginfo%nlon)

                  CALL ncio_write_serial (filename, 'lat', gcheck_concat%ginfo%lat_c, 'lat')
                  CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
                  CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

                  CALL ncio_write_serial (filename, 'lon', gcheck_concat%ginfo%lon_c, 'lon')
                  CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
                  CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')
#else
                  CALL ncio_define_dimension(filename, 'patch', numpatch)
#endif

               ENDIF


               IF (.not. end_of_spinup) THEN
                  write(timestr,timeform) i_spinupcycle, idate(1)
               ELSE
                  write(timestr,'(I4.4)') idate(1)
               ENDIF

               CALL nccheck( nf90_open(trim(filename), NF90_WRITE, ncid) )
               CALL nccheck( nf90_inq_varid(ncid, 'iyear', varid) )
               CALL nccheck( nf90_put_var(ncid, varid, timestr(1:timestrlen), &
                  (/1,nyearcheck/), (/timestrlen,1/)) )
               CALL nccheck( nf90_close(ncid) )

            ENDIF

#ifndef SinglePoint
            IF (p_is_io) CALL allocate_block_data (gridcheck, sumarea)

            CALL map_check%get_sumarea (sumarea, filter)

            CALL map_and_write_check_var ( &
               pct_dtws, filename, 'relative_tws_change', nyearcheck, sumarea, filter, &
               'The ratio of changes in terrestrial water storage to total precipitation', '-')

            CALL map_and_write_check_var ( &
               prcp_acc, filename, 'total_precipitation', nyearcheck, sumarea, filter, &
               'total precipitation in a year', 'mm')
#else
            CALL ncio_write_serial_time (filename, 'relative_tws_change', &
                                         nyearcheck, pct_dtws, 'patch', 'iyear')
            IF (nyearcheck == 1) THEN
               CALL ncio_put_attr (filename, 'relative_tws_change', 'long_name', &
                  'The ratio of changes in terrestrial water storage to total precipitation')
               CALL ncio_put_attr (filename, 'relative_tws_change', 'units', '-')
               CALL ncio_put_attr (filename, 'relative_tws_change', 'missing_value', spval)
            ENDIF

            CALL ncio_write_serial_time (filename, 'total_precipitation', &
                                         nyearcheck, prcp_acc, 'patch', 'iyear')
            IF (nyearcheck == 1) THEN
               CALL ncio_put_attr (filename, 'total_precipitation', 'long_name', &
                  'total precipitation in a year')
               CALL ncio_put_attr (filename, 'total_precipitation', 'units', 'mm')
               CALL ncio_put_attr (filename, 'total_precipitation', 'missing_value', spval)
            ENDIF
#endif
         ENDIF

         IF ((p_is_worker) .and. (numpatch > 0)) THEN
            tws_prev (0:windowsize-1,:) = tws_prev (1:windowsize,:)
            tws_prev (windowsize,:) = tws_this
            IF (nyearcheck == 1) tws_prev (windowsize-1,:) = tws_last

            prcp_prev(1:windowsize-1,:) = prcp_prev(2:windowsize,:)
            prcp_prev(windowsize,:) = prcp_acc

            prcp_acc(:) = spval
            tws_last = tws_this
         ENDIF

         IF (allocated(pct_dtws)) deallocate(pct_dtws)
         IF (allocated(filter  )) deallocate(filter  )

      ENDIF

   END SUBROUTINE CheckEquilibrium

   !-----------------------------------------------------------------------
   SUBROUTINE CheckEquilibrium_EndOfSpinup ()

   USE MOD_Precision
   USE MOD_LandPatch

   IMPLICIT NONE

   ! Local variables
   real(r8), allocatable :: prcp  (:), pct_dtws(:)
   logical,  allocatable :: filter(:)

   integer  :: actualsize, i
   real(r8) :: sum_pct_dtws

      IF (.not. DEF_CheckEquilibrium) RETURN

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN

            allocate (prcp     (numpatch))
            allocate (pct_dtws (numpatch))
            allocate (filter   (numpatch))

            actualsize = min(windowsize,nyearcheck)

            prcp = prcp_prev(windowsize,:)
            DO i = 1, actualsize-1
               CALL add_spv (prcp_prev(windowsize-i,:), prcp)
            ENDDO

            filter(:) = (tws_prev(windowsize,:) /= spval) &
               .and. (tws_prev(windowsize-actualsize,:) /= spval) .and. (prcp > 0.)

            WHERE (filter)
               pct_dtws = (tws_prev(windowsize,:) - tws_prev(windowsize-actualsize,:)) / prcp
            ELSEWHERE
               pct_dtws = spval
            END WHERE

            num_totalck = count(filter)

            IF (any(filter)) THEN
               sum_pct_dtws = sum(abs(pct_dtws), mask = filter)
            ELSE
               sum_pct_dtws = 0.
            ENDIF

            deallocate (prcp    )
            deallocate (pct_dtws)
            deallocate (filter  )

         ELSE
            num_totalck  = 0.
            sum_pct_dtws = 0.
         ENDIF

#ifdef USEMPI
         CALL mpi_allreduce (MPI_IN_PLACE, num_totalck,  1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
         CALL mpi_allreduce (MPI_IN_PLACE, sum_pct_dtws, 1, MPI_REAL8, MPI_SUM, p_comm_worker, p_err)
         IF (p_iam_worker == p_root) THEN
            CALL mpi_send (num_totalck,  1, MPI_REAL8, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
            CALL mpi_send (sum_pct_dtws, 1, MPI_REAL8, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN
#ifdef USEMPI
         CALL mpi_recv (num_totalck,  1, MPI_REAL8, p_address_worker(p_root), &
            mpi_tag_mesg, p_comm_glb, p_stat, p_err)
         CALL mpi_recv (sum_pct_dtws, 1, MPI_REAL8, p_address_worker(p_root), &
            mpi_tag_mesg, p_comm_glb, p_stat, p_err)
#endif

         IF (num_totalck > 0.) THEN
            ave_pct_dtws = sum_pct_dtws / num_totalck
         ELSE
            ave_pct_dtws = 0.
         ENDIF

         write(spinup_warning,'(A,F0.2,A)') 'Average delTWS/precipitation is ', ave_pct_dtws*100., &
            '% at the end of spinup. Check history for detail.'

      ENDIF

   END SUBROUTINE CheckEquilibrium_EndOfSpinup

   !-----------------------------------------------------------------------
#ifndef SinglePoint
   SUBROUTINE map_and_write_check_var ( &
         vector, filename, varname, itime_in_file, sumarea, filter, &
         longname, units)

   USE MOD_Block
   IMPLICIT NONE

   real(r8), intent(in) :: vector(:)

   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: varname
   integer,          intent(in) :: itime_in_file
   character(len=*), intent(in) :: longname
   character(len=*), intent(in) :: units

   type(block_data_real8_2d), intent(in) :: sumarea
   logical, intent(in) :: filter(:)

   ! Local variables
   type(block_data_real8_2d) :: data_xy_2d
   integer :: xblk, yblk, xloc, yloc, xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp
   integer :: iblkme, iblk, jblk, idata, ixseg, iyseg
   integer :: rmesg(3), smesg(3), isrc
   real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)

      IF (p_is_io) CALL allocate_block_data (gridcheck, data_xy_2d)
      CALL map_check%pset2grid (vector, data_xy_2d, spv = spval, msk = filter)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme
            xblk = gblock%xblkme(iblkme)
            yblk = gblock%yblkme(iblkme)

            DO yloc = 1, gridcheck%ycnt(yblk)
               DO xloc = 1, gridcheck%xcnt(xblk)

                  IF (sumarea%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                     IF (data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                        data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           = data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) &
                           / sumarea%blk(xblk,yblk)%val(xloc,yloc)
                     ENDIF
                  ELSE
                     data_xy_2d%blk(xblk,yblk)%val(xloc,yloc) = spval
                  ENDIF

               ENDDO
            ENDDO

         ENDDO
      ENDIF

      check_data_id = mod(check_data_id,100) + 1
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN

         allocate (vdata (gcheck_concat%ginfo%nlon, gcheck_concat%ginfo%nlat))
         vdata(:,:) = spval

#ifdef USEMPI
         DO idata = 1, gcheck_concat%ndatablk
            CALL mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
               check_data_id, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            ixseg = rmesg(2)
            iyseg = rmesg(3)

            xgdsp = gcheck_concat%xsegs(ixseg)%gdsp
            ygdsp = gcheck_concat%ysegs(iyseg)%gdsp
            xcnt  = gcheck_concat%xsegs(ixseg)%cnt
            ycnt  = gcheck_concat%ysegs(iyseg)%cnt

            allocate (rbuf(xcnt,ycnt))

            CALL mpi_recv (rbuf, xcnt*ycnt, MPI_DOUBLE, &
               isrc, check_data_id, p_comm_glb, p_stat, p_err)

            vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = rbuf
            deallocate (rbuf)

         ENDDO
#else
         DO iyseg = 1, gcheck_concat%nyseg
            DO ixseg = 1, gcheck_concat%nxseg
               iblk = gcheck_concat%xsegs(ixseg)%blk
               jblk = gcheck_concat%ysegs(iyseg)%blk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  xbdsp = gcheck_concat%xsegs(ixseg)%bdsp
                  ybdsp = gcheck_concat%ysegs(iyseg)%bdsp
                  xgdsp = gcheck_concat%xsegs(ixseg)%gdsp
                  ygdsp = gcheck_concat%ysegs(iyseg)%gdsp
                  xcnt  = gcheck_concat%xsegs(ixseg)%cnt
                  ycnt  = gcheck_concat%ysegs(iyseg)%cnt

                  vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                     data_xy_2d%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
               ENDIF
            ENDDO
         ENDDO
#endif

         CALL ncio_write_serial_time (filename, varname, itime_in_file, vdata, &
            'lon', 'lat', 'iyear', compress = 1)

         IF (itime_in_file == 1) THEN
            CALL ncio_put_attr (filename, varname, 'long_name', longname)
            CALL ncio_put_attr (filename, varname, 'units', units)
            CALL ncio_put_attr (filename, varname, 'missing_value', spval)
         ENDIF

         deallocate (vdata)

      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         DO iyseg = 1, gcheck_concat%nyseg
            DO ixseg = 1, gcheck_concat%nxseg

               iblk = gcheck_concat%xsegs(ixseg)%blk
               jblk = gcheck_concat%ysegs(iyseg)%blk

               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                  xbdsp = gcheck_concat%xsegs(ixseg)%bdsp
                  ybdsp = gcheck_concat%ysegs(iyseg)%bdsp
                  xcnt  = gcheck_concat%xsegs(ixseg)%cnt
                  ycnt  = gcheck_concat%ysegs(iyseg)%cnt

                  allocate (sbuf (xcnt,ycnt))
                  sbuf = data_xy_2d%blk(iblk,jblk)%val(xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                  smesg = (/p_iam_glb, ixseg, iyseg/)
                  CALL mpi_send (smesg, 3, MPI_INTEGER, &
                     p_address_master, check_data_id, p_comm_glb, p_err)
                  CALL mpi_send (sbuf, xcnt*ycnt, MPI_DOUBLE, &
                     p_address_master, check_data_id, p_comm_glb, p_err)

                  deallocate (sbuf)

               ENDIF
            ENDDO
         ENDDO
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE map_and_write_check_var
#endif

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

END MODULE MOD_CheckEquilibrium
