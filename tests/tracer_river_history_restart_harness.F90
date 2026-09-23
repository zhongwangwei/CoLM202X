! Real production tracer restart and history I/O; no mocked transport or I/O.
program tracer_river_history_restart_harness
 use MOD_Precision
 use MOD_SPMD_Task
 use MOD_Namelist
 use MOD_NetCDFSerial
 use MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_gdid, ucat_next, &
      ucat_data_address, x_ucat, y_ucat, griducat
 use MOD_Grid_RiverLakeHistState, only: acctime_ucat
 use MOD_Grid_RiverLakeHistRoute, only: route_hist_begin, route_hist_end
 use MOD_Tracer_RiverLake
 use MOD_Tracer_Defs, only: tracer_defs_init, tracers
 use netcdf
 implicit none
 character(len=32) :: mode
 character(len=64) :: name
 integer :: i, ncid, varid, nvars, itime
 logical :: found
 real(r8) :: expected(4), actual(1,1,1)
 call get_command_argument(1, mode)
 call spmd_init()
 call divide_processes_into_groups(1)
 if (p_np_worker /= 1) call CoLM_stop('harness requires three MPI ranks')
 totalnumucat = 1
 numucat = 0
 if (p_is_worker) numucat = 1
 allocate(ucat_gdid(numucat), ucat_next(numucat), acctime_ucat(numucat))
 ucat_gdid = 17
 ucat_next = 0
 allocate(ucat_data_address(0:0), x_ucat(1), y_ucat(1))
 allocate(ucat_data_address(0)%val(1))
 ucat_data_address(0)%val = 1
 x_ucat = 1
 y_ucat = 1
 griducat%nlon = 1
 griducat%nlat = 1
 DEF_HIST_mode = 'one'
 DEF_Reservoir_Method = 0
 DEF_USE_LEVEE = .true.
 DEF_USE_BIFURCATION = .true.
 DEF_HIST_CompressLevel = 0
 DEF_REST_CompressLevel = 0
 DEF_TRACER_NUM = 1
 DEF_TRACER_NAMES = 'dye'
 DEF_TRACER_TYPES = 'conservative'
 DEF_TRACER_MRAT = '1'
 DEF_TRACER_REF_RATIO = '1'
 DEF_TRACER_INIT_DELTA = '0'
 DEF_TRACER_REACTIVE_DECAY_RATE = '0'
 if (index(trim(mode), 'finite_') == 1) then
   call tracer_defs_init()
   tracers(1)%max_dissolved_conc = 2._r8
 endif
 call river_lake_tracer_init()
 if (index(trim(mode), 'finite_') == 1) then
   call finite_restart_check()
   call spmd_exit()
   stop
 endif
 if (p_is_worker) then
   ! Prefix has a different concentration and flux from the suffix.
   a_trc_storage_mass = 20._r8
   a_water_storage = 10._r8
   a_trc_levsto_mass = 12._r8
   a_levsto_water = 6._r8
   a_trc_out = 6._r8
   a_trc_bifout = -4._r8
   a_trc_acctime = 2._r8
   acctime_ucat = 2._r8
 endif
 if (trim(mode) /= 'continuous') then
   if (p_is_master) then
     call ncio_create_file('restart.nc')
     call ncio_define_dimension('restart.nc', 'ucatch', totalnumucat)
   endif
   call write_tracer_restart('restart.nc')
   if (p_is_master .and. trim(mode) /= 'split') then
     call nc_check(nf90_open('restart.nc', nf90_write, ncid))
     call nc_check(nf90_inquire(ncid, nVariables=nvars))
     call nc_check(nf90_redef(ncid))
     do varid = 1, nvars
       call nc_check(nf90_inquire_variable(ncid, varid, name=name))
       if (index(name, 'trc_hist_') /= 1) cycle
       if (trim(mode) == 'partial' .and. trim(name) /= 'trc_hist_out_dye') cycle
       if (trim(mode) == 'clockless' .and. trim(name) /= 'trc_hist_acctime') cycle
       call nc_check(nf90_rename_var(ncid, varid, 'old_'//trim(name)))
     enddo
     call nc_check(nf90_enddef(ncid))
     call nc_check(nf90_close(ncid))
   endif
   call mpi_barrier(p_comm_glb, p_err)
   ! Poison pre-read state: neither stale memory nor zero initialization may
   ! conceal a missed restore or an incomplete legacy-history reset.
   if (p_is_worker) then
     a_trc_storage_mass = 999._r8
     a_water_storage = 999._r8
     a_trc_levsto_mass = 999._r8
     a_levsto_water = 999._r8
     a_trc_out = 999._r8
     a_trc_bifout = 999._r8
     a_trc_acctime = 999._r8
   endif
   call read_tracer_restart('restart.nc', found)
   if (.not. found) call CoLM_stop('valid descriptor unexpectedly cold-started')
 endif
 if (p_is_worker) then
   a_trc_storage_mass = a_trc_storage_mass + 90._r8
   a_water_storage = a_water_storage + 30._r8
   a_trc_levsto_mass = a_trc_levsto_mass + 24._r8
   a_levsto_water = a_levsto_water + 12._r8
   a_trc_out = a_trc_out + 24._r8
   a_trc_bifout = a_trc_bifout + 12._r8
   a_trc_acctime = a_trc_acctime + 3._r8
   acctime_ucat = 5._r8 ! Legacy tracer history must NOT divide by this clock.
 endif
 call route_hist_begin('history.nc', [2000, 1, 0], .true., [0._r8], [0._r8], itime)
 call write_tracer_history('history.nc', itime, acctime_ucat)
 call route_hist_end()
 if (p_is_master) then
   expected = [2.75_r8, 6._r8, 7.2_r8, 1.6_r8]
   if (trim(mode) == 'legacy' .or. trim(mode) == 'partial') &
     expected = [3._r8, 8._r8, 8._r8, 4._r8]
   call nc_check(nf90_open('history.nc', nf90_nowrite, ncid))
   do i = 1, 4
     select case(i)
     case(1); name = 'f_trc_conc_dye'
     case(2); name = 'f_trc_flux_dye'
     case(3); name = 'f_trc_levsto_dye'
     case(4); name = 'f_trc_bifout_dye'
     end select
     call nc_check(nf90_inq_varid(ncid, trim(name), varid))
     call nc_check(nf90_get_var(ncid, varid, actual))
     if (abs(actual(1,1,1) - expected(i)) > 1.e-6_r8) then
       write(*,*) trim(name), actual, 'expected', expected(i)
       call CoLM_stop('restart history mismatch')
     endif
   enddo
   call nc_check(nf90_close(ncid))
   write(*,'(A)') 'TRACER_HISTORY_RESTART_OK'
 endif
 call spmd_exit()
contains
 subroutine finite_restart_check()
   if (trim(mode) /= 'finite_read_only') then
     if (p_is_worker) then
       trc_mass(1,1) = 4._r8
       trc_solid(1,1) = 6._r8
       trc_levsto(1,1) = 2._r8
       trc_levsto_solid(1,1) = 3._r8
     endif
     if (p_is_master) then
       call ncio_create_file('restart.nc')
       call ncio_define_dimension('restart.nc', 'ucatch', totalnumucat)
     endif
     call write_tracer_restart('restart.nc')
   endif
   if (trim(mode) == 'finite_write_only') return
   call mpi_barrier(p_comm_glb, p_err)
   if (p_is_worker) then
     trc_mass = 999._r8
     trc_solid = 999._r8
     trc_levsto = 999._r8
     trc_levsto_solid = 999._r8
   endif
   call read_tracer_restart('restart.nc', found)
   if (.not. found) call CoLM_stop('finite restart unexpectedly cold-started')
   if (p_is_worker) then
     if (abs(trc_mass(1,1)-4._r8) > 1.e-12_r8 .or. &
         abs(trc_levsto(1,1)-2._r8) > 1.e-12_r8) call CoLM_stop('finite restart mobile mass mismatch')
     if (trim(mode) == 'finite_read_only') then
       if (abs(trc_solid(1,1)) > 1.e-12_r8 .or. &
           abs(trc_levsto_solid(1,1)) > 1.e-12_r8) call CoLM_stop('schema1 solids not initialized')
     else
       if (abs(trc_solid(1,1)-6._r8) > 1.e-12_r8 .or. &
           abs(trc_levsto_solid(1,1)-3._r8) > 1.e-12_r8) call CoLM_stop('schema2 solids not restored')
     endif
   endif
   if (p_is_master) write(*,'(A)') 'FINITE_RESTART_OK'
 end subroutine
 subroutine nc_check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) call CoLM_stop(nf90_strerror(status))
 end subroutine
end program
