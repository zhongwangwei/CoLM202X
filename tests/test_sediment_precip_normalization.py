"""Compile shipped precipitation accumulation/yield bodies, not a Python replica."""
from pathlib import Path
import re
import subprocess

from fortran_test_support import require_runnable_fortran_compiler

ROOT = Path(__file__).resolve().parents[1]
SED = ROOT / 'main/TRACER/MOD_Tracer_Particle_Sediment.F90'
FLOW = ROOT / 'main/HYDRO/MOD_Grid_RiverLakeFlow.F90'


def routine(src, name):
    return re.search(rf'^\s*SUBROUTINE {name}\b.*?^\s*END SUBROUTINE {name}\b',
                     src, re.M | re.S | re.I).group()


def test_valid_area_time_preserves_nonlinear_yield(tmp_path):
    fc = require_runnable_fortran_compiler(tmp_path)
    src = SED.read_text()
    body = routine(src, 'sediment_forcing_put') + '\n' + routine(src, 'calc_sediment_yield')
    # The probe must use the constant the model ships, not its own copy of it.
    threshold = re.search(r'SED_PRECIP_THRESHOLD_MM_DAY\s*=\s*([0-9.]+)_r8', src).group(1)
    # Three cells: varying valid area/time, fully missing, and valid dry rain.
    code = '''
module MOD_Grid_RiverLakeNetwork
 integer :: numucat=3
end module
module probe
 use, intrinsic :: ieee_arithmetic
 implicit none
 integer,parameter :: r8=kind(1.d0), nlfp_sed=1
 logical :: p_is_worker=.true.
 real(r8) :: sed_precip(3)=0, sed_precip_yield(3)=0, sed_precip_time(3)=0
 real(r8) :: pyldpc=2, pyld=1, pyldc=1, dsylunit=1
 real(r8) :: sedinp(1,3)=0, sed_slope(1,3)=1, sed_frc(1,3)=1
 real(r8),parameter :: SED_PRECIP_THRESHOLD_MM_DAY=@THRESHOLD@_r8
 contains
 logical function sediment_particle_enabled()
 sediment_particle_enabled=.true.
 end function
 subroutine CoLM_stop(message)
 character(*), optional :: message
 error stop 'invalid forcing'
 end subroutine
''' + body + '''
end module
program check
 use MOD_Grid_RiverLakeNetwork
 use probe
 real(r8) :: p(3), a(3), nan
 nan=ieee_value(0._r8,ieee_quiet_nan)
 ! Valid rainfall 10 mm/h over half the area for 600 s.
 p=[10._r8/3600, 0._r8, 0._r8]; a=[0.5_r8,0._r8,1._r8]
 call sediment_forcing_put(p,600._r8,a)
 ! 20 mm/h over all the area for 300 s; missing cell still absent.
 p=[20._r8/3600,nan,0._r8]; a=[1._r8,0._r8,1._r8]
 call sediment_forcing_put(p,300._r8,a)
 if(any(abs(sed_precip_time-[600._r8,0._r8,900._r8])>1.e-10_r8)) stop 1
 if(abs(sed_precip(1)/sed_precip_time(1)*3600-15._r8)>1.e-10_r8) stop 2
 if(abs(sed_precip_yield(1)/sed_precip_time(1)-250._r8)>1.e-10_r8) stop 3
 call calc_sediment_yield([0._r8,0._r8,0._r8], [100._r8,100._r8,100._r8], sed_precip_time)
 if(abs(sedinp(1,1)-250._r8*100/3600)>1.e-10_r8) stop 4
 if(any(sedinp(1,2:3)/=0._r8)) stop 5
 ! Omitted coverage retains fully covered callers' historical dt behavior.
 sed_precip=0; sed_precip_yield=0; sed_precip_time=0
 call sediment_forcing_put([10._r8/3600,0._r8,0._r8],900._r8)
 if(any(sed_precip_time/=900._r8)) stop 6
 ! The threshold is per forcing step, as in CaMa-Flood's prcp_convert_sed: a
 ! step at or below it yields nothing and one above it does, whatever the
 ! window mean is.
 sed_precip=0; sed_precip_yield=0; sed_precip_time=0; sedinp=0
 p=[0.9_r8*SED_PRECIP_THRESHOLD_MM_DAY/86400,0._r8,0._r8]; a=[1._r8,1._r8,1._r8]
 call sediment_forcing_put(p,3600._r8,a)
 call calc_sediment_yield([0._r8,0._r8,0._r8], [100._r8,100._r8,100._r8], sed_precip_time)
 if(any(sedinp/=0._r8)) stop 7
 sed_precip=0; sed_precip_yield=0; sed_precip_time=0; sedinp=0
 p=[1.1_r8*SED_PRECIP_THRESHOLD_MM_DAY/86400,0._r8,0._r8]
 call sediment_forcing_put(p,3600._r8,a)
 call calc_sediment_yield([0._r8,0._r8,0._r8], [100._r8,100._r8,100._r8], sed_precip_time)
 if(sedinp(1,1)<=0._r8) stop 8
 ! A dry hour followed by a short strong burst: the window MEAN is well below
 ! the threshold, but the burst is a rain step above it and must still yield,
 ! diluted by the whole window's exposure time.
 sed_precip=0; sed_precip_yield=0; sed_precip_time=0; sedinp=0
 p=[0._r8,0._r8,0._r8]; a=[1._r8,1._r8,1._r8]
 call sediment_forcing_put(p,3600._r8,a)
 p=[1.5_r8*SED_PRECIP_THRESHOLD_MM_DAY/86400,0._r8,0._r8]
 call sediment_forcing_put(p,600._r8,a)
 if(sed_precip(1)/sed_precip_time(1)*86400 >= SED_PRECIP_THRESHOLD_MM_DAY) stop 9
 call calc_sediment_yield([0._r8,0._r8,0._r8], [100._r8,100._r8,100._r8], sed_precip_time)
 if(abs(sedinp(1,1) - ((1.5_r8*SED_PRECIP_THRESHOLD_MM_DAY/24)**2*600._r8/4200._r8)*100/3600) > 1.e-12_r8) stop 10
 ! Sub-threshold weak steps add exposure but no yield.
 sed_precip=0; sed_precip_yield=0; sed_precip_time=0; sedinp=0
 p=[0.5_r8*SED_PRECIP_THRESHOLD_MM_DAY/86400,0._r8,0._r8]
 call sediment_forcing_put(p,3600._r8,a)
 if(sed_precip_yield(1)/=0._r8 .or. sed_precip_time(1)/=3600._r8) stop 11
 ! Empty worker remains a valid no-op.
 numucat=0
 call sediment_forcing_put(p(:0),300._r8,a(:0))
 print *, 'PRECIP_OK'
end program
'''
    path = tmp_path / 'probe.f90'
    path.write_text(code.replace('@THRESHOLD@', threshold))
    exe = tmp_path / 'probe'
    result = subprocess.run([fc, '-ffree-line-length-0', '-fcheck=all',
                             '-ffpe-trap=invalid,zero,overflow', str(path), '-o', str(exe)],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    result = subprocess.run([str(exe)], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'PRECIP_OK' in result.stdout


def test_precipitation_mapping_uses_matching_mask_and_area():
    src = FLOW.read_text()
    block = src.split('IF (tracer_lifecycle_route_has_active()) THEN', 1)[1].split('#endif', 1)[0]
    assert 'filter_prcp' in block
    assert 'forcmask_pch' in block  # current forcing mask, not init-time runoff mask
    assert 'forc_prc(i) /= spval' in block and 'forc_prl(i) /= spval' in block
    assert 'ieee_is_finite' in block
    assert 'prcp_uc = prcp_uc / prcp_area_uc' in block
    # rain and valid area are exchanged together in one batched push
    assert 'worker_push_data(push_inpm2ucat, prcp_push_fields, mode = ' in block
    assert block.count('worker_push_data(') == 1
    assert 'fillvalue = 0._r8, filter = filter_prcp' in block
    assert 'prcp_area_uc / topo_area' in block
    assert 'tracer_lifecycle_route_forcing_put(prcp_uc, deltime, prcp_area_uc)' in block
    # Zero rain is valid area, not indistinguishable from missing data.
    assert 'prcp_pch = 1._r8' in block


def test_restart_preserves_each_cells_valid_exposure():
    src = SED.read_text()
    read = routine(src, 'read_sediment_restart')
    write = routine(src, 'write_sediment_restart')
    assert 'sed_precip_time(:) = buf(:)' in read
    assert 'sed_precip_time = buf(1)' not in read
    assert "'sed_precip_time_vec'" in write
    assert 'scalar_vec(:) = sed_precip_time' in write


def test_real_mpi_mapping_with_missing_dry_and_empty_workers(tmp_path):
    import os
    import shutil
    import pytest

    fc, launcher = shutil.which('mpif90'), shutil.which('mpiexec')
    if not fc or not launcher:
        pytest.skip('MPI compiler/runtime unavailable')
    # Extract the shipped mapping branch; link the real worker remap/push module.
    block = FLOW.read_text().split('IF (tracer_lifecycle_route_has_active()) THEN', 1)[1]
    block = block.split('\n         ENDIF\n#endif', 1)[0]
    driver = r'''
program mapping
 use MOD_WorkerPushData
 use, intrinsic :: ieee_arithmetic
 implicit none
 integer :: numpatch,numinpm,numucat,rank,nranks,i,scenario
 integer,allocatable :: patchtype(:),ids(:),req(:,:)
 logical,allocatable :: patchmask(:),forcmask_pch(:),filter_prcp(:)
 real(r8),allocatable :: forc_prc(:),forc_prl(:),topo_area(:),areas(:,:)
 real(r8),allocatable :: prcp_pch(:)
 real(r8),allocatable,target :: prcp_gd(:),prcp_uc(:),prcp_area_gd(:),prcp_area_uc(:)
 type(worker_push_real8_field_type) :: prcp_push_fields(2)
 real(r8),parameter :: spval=-1.e36_r8,deltime=300._r8
 type :: forcing_type
 logical :: has_missing_value=.true.
 end type
 type(forcing_type) :: DEF_forcing
 type(worker_remapdata_type) :: remap_patch2inpm
 type(worker_pushdata_type) :: push_inpm2ucat,push_ucat2inpm
 call mpi_init(p_err)
 call mpi_comm_rank(MPI_COMM_WORLD,rank,p_err)
 call mpi_comm_size(MPI_COMM_WORLD,nranks,p_err)
 p_is_worker=.true.; p_comm_worker=MPI_COMM_WORLD
 p_iam_worker=rank; p_np_worker=nranks
 numpatch=0; numinpm=0; numucat=0
 if(rank==0) then
 numpatch=2; numinpm=1
 endif
 if(rank==1) numucat=1
 allocate(patchtype(numpatch),patchmask(numpatch),forcmask_pch(numpatch))
 allocate(forc_prc(numpatch),forc_prl(numpatch),topo_area(numucat))
 allocate(ids(numinpm),req(1,numucat),areas(1,numucat))
 ids=101; req=101; areas=100._r8; topo_area=100._r8
 call build_worker_pushdata(numinpm,ids,numucat,req,areas,push_inpm2ucat)
 allocate(push_ucat2inpm%sum_area(numinpm)); push_ucat2inpm%sum_area=100._r8
 remap_patch2inpm%npset=numpatch; remap_patch2inpm%num_grid=numinpm
 allocate(remap_patch2inpm%npart(numpatch),remap_patch2inpm%part_to(numpatch), &
          remap_patch2inpm%areapart(numpatch))
 remap_patch2inpm%npart=1
 do i=1,numpatch
 allocate(remap_patch2inpm%part_to(i)%val(1),remap_patch2inpm%areapart(i)%val(1))
 remap_patch2inpm%part_to(i)%val=1; remap_patch2inpm%areapart(i)%val=50._r8
 enddo
 do scenario=1,6
 patchtype=0; patchmask=.true.; forcmask_pch=.true.; forc_prc=10._r8; forc_prl=0._r8
 DEF_forcing%has_missing_value=.true.
 if(rank==0) then
 select case(scenario)
 case(1)
 forcmask_pch(2)=.false. ! Valid half must still be 10, not 5.
 case(2)
 forc_prc(2)=0._r8 ! Dry half IS valid: catchment mean 5.
 case(3)
 forcmask_pch=.false. ! Zero exposure, not observed drought.
 case(4)
 forc_prl(2)=spval; DEF_forcing%has_missing_value=.false.
 case(5)
 forc_prc(2)=ieee_value(0._r8,ieee_quiet_nan)
 case(6)
 forc_prl(2)=-1._r8
 end select
 endif
''' + block + r'''
 enddo
 call mpi_barrier(MPI_COMM_WORLD,p_err)
 if(rank==0) print *, 'MAPPING_OK'
 call mpi_finalize(p_err)
 contains
 subroutine tracer_lifecycle_route_forcing_put(rain,dt,coverage)
 real(r8),intent(in) :: rain(:),dt,coverage(:)
 real(r8) :: expected_rain,expected_area
 if(rank/=1) return
 expected_rain=10; expected_area=0.5_r8
 if(scenario==2) then
 expected_rain=5; expected_area=1
 elseif(scenario==3) then
 expected_rain=0; expected_area=0
 endif
 if(size(rain)/=1) error stop 'wrong output size'
 if(abs(rain(1)-expected_rain)>1.e-10_r8) error stop 'rain diluted'
 if(abs(coverage(1)-expected_area)>1.e-10_r8) error stop 'bad valid area'
 end subroutine
end program
'''
    (tmp_path / 'define.h').write_text('#define USEMPI\n')
    (tmp_path / 'mapping.f90').write_text(driver)
    exe = tmp_path / 'mapping'
    result = subprocess.run([fc, '-cpp', '-fdefault-real-8', '-ffree-line-length-0',
                             '-fallow-argument-mismatch', '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                             f'-I{tmp_path}', f'-J{tmp_path}',
                             str(ROOT / 'tests/river_mpi_test_support.F90'),
                             str(ROOT / 'share/MOD_WorkerPushData.F90'),
                             str(tmp_path / 'mapping.f90'), '-o', str(exe)],
                            cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1',
               OMPI_MCA_rmaps_base_oversubscribe='1')
    result = subprocess.run([launcher, '-n', '4', str(exe)], cwd=tmp_path, env=env,
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'MAPPING_OK' in result.stdout
