from pathlib import Path
import os
import shutil
import subprocess
import sys

import pytest

# Point at a fresh CaMa-enabled build; never accidentally test default/empty modules.
BLD = Path(os.environ['COLM_CAMA_BUILD_DIR']) if 'COLM_CAMA_BUILD_DIR' in os.environ else None

DRIVER = r'''
program cama_mpi_probe
  use MOD_Precision
  use MOD_SPMD_Task
  use MOD_Block, only: gblock
  use MOD_Grid, only: segment_type
  use MOD_DataType, only: block_data_real8_2d, allocate_block_data
  use MOD_LandPatch, only: numpatch
  use MOD_Vars_TimeInvariants, only: patchtype
  use MOD_CaMa_Vars
  implicit none
  integer :: ierr, rank, nprocs
  real(r8), allocatable :: worker(:), outworker(:), master(:,:), sinkmaster(:,:)
  type(block_data_real8_2d) :: io

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  if (nprocs /= 4) error stop 'need 4 ranks'
  call setup_spmd(rank,nprocs)
  call setup_grid()
  call setup_mapping(rank,1)

  nacc=1._r8
  if (p_is_worker) then
    if (rank == 1) then
      numpatch=1
      allocate(worker(1), outworker(1), patchtype(1))
      patchtype=1
      worker=[10._r8]
    else
      numpatch=0
      allocate(worker(0), outworker(0), patchtype(0))
    endif
  else
    allocate(worker(0), outworker(0))
  endif
  if (p_is_io) call allocate_block_data(gcama, io)
  if (p_is_master) then
    allocate(master(2,1), sinkmaster(2,1)); master=-999._r8; sinkmaster=-999._r8
  else
    allocate(master(0,0), sinkmaster(0,0))
  endif

  call colm2cama_real8(worker, io, master, integral=.true.)
  if (p_is_master) then
    if (abs(master(1,1)-20000._r8)>1.e-8_r8 .or. abs(master(2,1)-60000._r8)>1.e-8_r8) error stop 'integral conservation failed'
    sinkmaster=0._r8
  endif

  if (p_is_worker) then
    if (rank == 1) then
      if (allocated(flood_credit)) deallocate(flood_credit)
      if (allocated(flood_credit_part)) deallocate(flood_credit_part)
      if (allocated(flood_sink_part)) deallocate(flood_sink_part)
      allocate(flood_credit(1), flood_credit_part(1), flood_sink_part(1))
      allocate(flood_credit_part(1)%val(2), flood_sink_part(1)%val(2))
      flood_credit=1.625_r8; flood_credit_part(1)%val=[0.5_r8,2._r8]; worker=[10._r8]
    else
      if (allocated(flood_credit)) deallocate(flood_credit)
      if (allocated(flood_credit_part)) deallocate(flood_credit_part)
      if (allocated(flood_sink_part)) deallocate(flood_sink_part)
      allocate(flood_credit(0), flood_credit_part(0), flood_sink_part(0))
    endif
  endif
  call colm2cama_real8(worker, io, sinkmaster, integral=.true., flood_sink=.true.)
  if (p_is_master) then
    if (abs(sinkmaster(1,1)-10000._r8/1.625_r8)>1.e-8_r8 .or. &
        abs(sinkmaster(2,1)-120000._r8/1.625_r8)>1.e-8_r8) error stop 'sink provenance failed'
    master(1,1)=4._r8; master(2,1)=8._r8
  endif
  call cama2colm_real8(master, io, outworker)
  if (p_is_worker .and. rank == 1) then
    if (abs(outworker(1)-7._r8)>1.e-8_r8) error stop 'cama2colm weighted mean failed'
  endif
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (rank == 0) print *, 'mpi toy ok'
  call mpi_finalize(ierr)
contains
  subroutine setup_spmd(rank,nprocs)
    integer,intent(in):: rank,nprocs
    p_comm_glb=MPI_COMM_WORLD; p_iam_glb=rank; p_np_glb=nprocs
    p_address_master=3
    p_is_master=(rank==3); p_is_io=(rank==0); p_is_worker=(rank==1 .or. rank==2); p_is_writeback=.false.
    p_np_io=1; if (allocated(p_address_io)) deallocate(p_address_io); allocate(p_address_io(0:0)); p_address_io=[0]
    p_np_worker=2; if (allocated(p_address_worker)) deallocate(p_address_worker); allocate(p_address_worker(0:1)); p_address_worker=[1,2]
    p_iam_io=0; p_iam_worker=-1
    if(rank==1) p_iam_worker=0
    if(rank==2) p_iam_worker=1
  end subroutine
  subroutine setup_grid()
    gblock%nxblk=1; gblock%nyblk=1
    if(allocated(gblock%pio)) deallocate(gblock%pio)
    allocate(gblock%pio(1,1)); gblock%pio=0
    if(allocated(gblock%xblkme)) deallocate(gblock%xblkme,gblock%yblkme)
    if(p_is_io) then
      gblock%nblkme=1; allocate(gblock%xblkme(1),gblock%yblkme(1)); gblock%xblkme=1; gblock%yblkme=1
    else
      gblock%nblkme=0; allocate(gblock%xblkme(0),gblock%yblkme(0))
    endif
    call init_grid_fields(mp2g_cama%grid); call init_grid_fields(mg2p_cama%grid); call init_grid_fields(gcama)
    cama_gather%ndatablk=1; cama_gather%nxseg=1; cama_gather%nyseg=1
    allocate(cama_gather%xsegs(1),cama_gather%ysegs(1))
    cama_gather%xsegs(1)=segment_type(1,2,0,0); cama_gather%ysegs(1)=segment_type(1,1,0,0)
  end subroutine
  subroutine init_grid_fields(grid)
    use MOD_Grid, only: grid_type
    type(grid_type),intent(inout):: grid
    grid%nlon=2; grid%nlat=1; grid%yinc=1
    allocate(grid%xcnt(1),grid%ycnt(1),grid%xdsp(1),grid%ydsp(1))
    grid%xcnt=[2]; grid%ycnt=[1]; grid%xdsp=[0]; grid%ydsp=[0]
    allocate(grid%xblk(2),grid%yblk(1),grid%xloc(2),grid%yloc(1))
    grid%xblk=[1,1]; grid%yblk=[1]; grid%xloc=[1,2]; grid%yloc=[1]
  end subroutine
  subroutine setup_mapping(rank,patch_count)
    integer,intent(in):: rank,patch_count
    call setup_one_mapping(mp2g_cama, rank,patch_count); call setup_one_mapping(mg2p_cama, rank,patch_count)
  end subroutine
  subroutine setup_one_mapping(map, rank,patch_count)
    use MOD_SpatialMapping, only: spatial_mapping_type
    type(spatial_mapping_type), intent(inout) :: map
    integer,intent(in):: rank,patch_count
    integer :: i
    map%npset=0; allocate(map%glist(0:1))
    do i=0,1; map%glist(i)%ng=0; enddo
    if (p_is_io) then
      map%glist(0)%ng=2; allocate(map%glist(0)%ilon(2),map%glist(0)%ilat(2)); map%glist(0)%ilon=[1,2]; map%glist(0)%ilat=[1,1]
    elseif (p_is_worker .and. rank==1) then
      map%npset=patch_count
      map%glist(0)%ng=2; allocate(map%glist(0)%ilon(2),map%glist(0)%ilat(2)); map%glist(0)%ilon=[1,2]; map%glist(0)%ilat=[1,1]
      allocate(map%npart(patch_count),map%address(patch_count),map%areapart(patch_count),map%areapset(patch_count))
      map%npart=2; map%areapset=8._r8
      do i=1,patch_count
        allocate(map%address(i)%val(2,2),map%areapart(i)%val(2))
        map%address(i)%val(:,1)=[0,1]; map%address(i)%val(:,2)=[0,2]
        map%areapart(i)%val=[2._r8,6._r8]
      enddo
    elseif (p_is_worker) then
      map%npset=0; allocate(map%npart(0),map%address(0),map%areapart(0),map%areapset(0))
    endif
    if (p_is_io) then
      allocate(map%areagrid%blk(1,1)); allocate(map%areagrid%blk(1,1)%val(2,1)); map%areagrid%blk(1,1)%val(:,1)=[2._r8,6._r8]*patch_count
    endif
  end subroutine
end program
'''


def _run_probe(tmp_path: Path, driver: str) -> None:
    if BLD is None or not BLD.exists():
        pytest.skip("set COLM_CAMA_BUILD_DIR to a fresh CaMa-enabled .bld directory")
    if shutil.which("mpifort") is None or shutil.which("mpiexec") is None:
        pytest.skip("MPI compiler/runtime unavailable")
    src = tmp_path / "cama_mpi_probe.f90"
    exe = tmp_path / "cama_mpi_probe"
    src.write_text(driver)
    objects = [str(p) for p in BLD.glob("*.o") if p.name != "CoLM.o"]
    library = tmp_path / "libcama_probe.a"
    subprocess.run(["ar", "rcs", str(library), *objects], check=True, capture_output=True)
    libs = subprocess.check_output(["pkg-config", "--libs", "netcdf-fortran", "netcdf"], text=True).split()
    strip_flag = "-Wl,-dead_strip" if sys.platform == "darwin" else "-Wl,--gc-sections"
    compile_cmd = ["mpifort", "-fopenmp", "-fcheck=all", f"-I{BLD}", str(src), str(library), strip_flag, *libs, "-o", str(exe)]
    subprocess.run(compile_cmd, cwd=tmp_path, check=True, text=True, capture_output=True)
    env = os.environ.copy(); env["OMP_NUM_THREADS"] = "2"
    result = subprocess.run(["mpiexec", "-n", "4", str(exe)], cwd=tmp_path, env=env, text=True, capture_output=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "mpi toy ok" in result.stdout.lower()


MEAN_DRIVER = r"""
program cama_mpi_probe
  use MOD_Precision
  use MOD_SPMD_Task
  use MOD_Block, only: gblock
  use MOD_Grid, only: segment_type
  use MOD_DataType, only: block_data_real8_2d, allocate_block_data
  use MOD_LandPatch, only: numpatch
  use MOD_Vars_TimeInvariants, only: patchtype
  use MOD_Vars_1DFluxes, only: rnof
  use MOD_Vars_1DForcing, only: forc_rain
  use MOD_Forcing, only: forcmask_pch
  use MOD_Namelist, only: DEF_forcing
  use YOS_CMF_INPUT, only: LWEVAP,LWINFILT,LDAMIRR
  use MOD_CaMa_Vars
  implicit none
  integer :: ierr,rank,nprocs
  real(r8), allocatable :: master(:,:)
  type(block_data_real8_2d) :: io
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
  if(nprocs/=4) error stop 'need 4 ranks'
  call setup_spmd(rank,nprocs)
  call setup_grid()
  call setup_mapping(rank,2)
  if(p_is_worker) then
    numpatch=0
    if(rank==1) numpatch=2
    allocate(patchtype(numpatch),rnof(numpatch),forc_rain(numpatch),forcmask_pch(numpatch))
    patchtype=0; rnof=0._r8
  endif
  call allocate_acc_cama_fluxes()
  if(p_is_io) call allocate_block_data(gcama,io)
  if(p_is_master) then
    allocate(master(2,1))
  else
    allocate(master(0,0))
  endif
  LSEDIMENT=.true.; LWEVAP=.false.; LWINFILT=.false.; LDAMIRR=.false.
  DEF_forcing%has_missing_value=.true.

  ! One patch missing all window: its area must not dilute valid rain.
  call flush_acc_cama_fluxes()
  if(rank==1) then
    forc_rain=[10._r8,spval]; forcmask_pch=[.true.,.false.]
  endif
  call accumulate_cama_fluxes(1800._r8)
  call colm2cama_real8(a_prcp_cama,io,master,valid_time=a_prcp_time)
  if(p_is_master) then
    if(any(abs(master-10._r8)>1.e-12_r8)) error stop 'missing area diluted rain'
  endif

  ! Unequal timestep lengths and changing masks. Patch 1: 10 for 900s,
  ! patch 2: 20 for 900s and 40 for 2700s => area-time mean 30.
  ! Correct numerator is 2*(9000+18000+108000), denominator 2*(900+3600).
  call flush_acc_cama_fluxes()
  if(rank==1) then
    forc_rain=[10._r8,20._r8]; forcmask_pch=.true.; rnof=[10._r8,20._r8]
  endif
  call accumulate_cama_fluxes(900._r8)
  if(rank==1) then
    forc_rain=[spval,40._r8]; forcmask_pch=[.false.,.true.]; rnof=[spval,40._r8]
  endif
  call accumulate_cama_fluxes(2700._r8)
  if(rank==1) then
    if(any(abs(a_prcp_time-[900._r8,3600._r8])>1.e-12_r8)) error stop 'valid durations wrong'
  endif
  call colm2cama_real8(a_prcp_cama,io,master,valid_time=a_prcp_time)
  if(p_is_master) then
    if(any(abs(master-30._r8)>1.e-12_r8)) error stop 'missing time diluted rain'
  endif
  call colm2cama_real8(a_rnof_cama,io,master,integral=.true.)
  if(p_is_master) then
    if(abs(master(1,1)-75000._r8)>1.e-8_r8.or.abs(master(2,1)-225000._r8)>1.e-8_r8) &
      error stop 'volume flux was incorrectly renormalized'
  endif

  ! Raw sentinel is excluded even if forcing mask was not set.
  call flush_acc_cama_fluxes()
  DEF_forcing%has_missing_value=.false.
  if(rank==1) forc_rain=[10._r8,spval]
  call accumulate_cama_fluxes(3600._r8)
  if(rank==1) then
    if(a_prcp_time(2)/=0._r8) error stop 'sentinel counted as valid time'
  endif
  call colm2cama_real8(a_prcp_cama,io,master,valid_time=a_prcp_time)
  if(p_is_master) then
    if(any(abs(master-10._r8)>1.e-12_r8)) error stop 'sentinel area diluted rain'
  endif

  ! Generic intensive path uses the same data-valid mask in sumarea.
  call flush_acc_cama_fluxes()
  if(rank==1) forc_rain=[10._r8,spval]
  call accumulate_cama_fluxes(3600._r8)
  call colm2cama_real8(a_prcp_cama,io,master)
  if(p_is_master) then
    if(any(abs(master-10._r8)>1.e-12_r8)) error stop 'generic mean denominator mismatch'
  endif

  ! Entire grid missing: no NaN or divide by zero, and flush clears weights.
  call flush_acc_cama_fluxes()
  if(rank==1) forc_rain=spval
  call accumulate_cama_fluxes(3600._r8)
  call colm2cama_real8(a_prcp_cama,io,master,valid_time=a_prcp_time)
  if(p_is_master) then
    if(any(master/=0._r8)) error stop 'empty mean must be zero'
  endif
  call deallocate_acc_cama_fluxes()
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if(rank==0) print *, 'mpi toy ok'
  call mpi_finalize(ierr)
contains
""" + DRIVER.split("contains\n", 1)[1]


def test_actual_cama_mpi_mapping_toy(tmp_path: Path) -> None:
    _run_probe(tmp_path, DRIVER)


def test_actual_sediment_mean_excludes_missing_area_and_time(tmp_path: Path) -> None:
    _run_probe(tmp_path, MEAN_DRIVER)
