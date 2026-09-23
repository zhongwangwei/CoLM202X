from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path: str) -> str:
    return (ROOT / path).read_text()


def flat(path: str) -> str:
    return " ".join(read(path).split())


def test_restart_structure_group_is_explicit_not_silent_default():
    src = read("main/MOD_Vars_TimeInvariants.F90")
    compact = flat("main/MOD_Vars_TimeInvariants.F90")

    assert "ncio_vector_group_presence(file_restart, canopy_fields, landpatch, canopy_present)" in src
    assert "ncio_vector_group_presence(file_restart, canopy_fields, landpft, canopy_present)" in src
    assert "any(canopy_present) .and. .not. all(canopy_present)" in src
    assert "DEF_Interception_scheme = 1" in src
    assert "restart lacks ncd/ncw/bcw" in src
    assert "restart lacks ncd_p/ncw_p/bcw_p" in src
    assert "CALL ncio_read_vector (file_restart, 'ncd' , landpatch, ncd, defval" not in compact
    assert "CALL ncio_read_vector (file_restart, 'ncd_p' , landpft, ncd_p, defval" not in compact


def test_mkinidata_structure_inputs_checked_before_scheme8_use():
    src = read("mkinidata/MOD_HtopReadin.F90")

    assert "ncio_vector_var_present" in src
    assert "ncd_patches.nc" in src and "ncw_patches.nc" in src and "bcw_patches.nc" in src
    assert "ncd_pfts.nc" in src and "ncw_pfts.nc" in src and "bcw_pfts.nc" in src
    assert "DEF_Interception_scheme = 1" in src
    assert "using CoLM2014 interception" in src
    assert "defval=-1.0e36_r8" not in src
    assert "surface canopy structure is invalid" in src
    assert "surface PFT canopy structure is invalid" in src
    assert "ieee_is_finite([ncd(icanopy), ncw(icanopy)])" in src
    assert "ieee_is_finite([bcw_p(icanopy), htop_p(icanopy)])" in src
    assert "Downgrading the whole run to CoLM2014 interception; no canopy structure is synthesized." in src


def test_partial_plain_igbp_structure_uses_global_non_synthetic_downgrade():
    mkinidata = read("mkinidata/MOD_HtopReadin.F90")
    restart = read("main/MOD_Vars_TimeInvariants.F90")

    warning = "Downgrading the whole run to CoLM2014 interception; no canopy structure is synthesized."
    assert warning in mkinidata
    assert warning in restart

    mkinidata_partial = mkinidata.split("ELSEIF (canopy_counts(2) < canopy_counts(1)) THEN", 1)[1]
    mkinidata_partial = mkinidata_partial.split("ENDIF", 1)[0]
    assert "DEF_Interception_scheme = 1" in mkinidata_partial
    assert "CALL CoLM_stop" not in mkinidata_partial

    patch_restart = restart.split("SUBROUTINE READ_TimeInvariants ", 1)[1]
    restart_partial = patch_restart.split("ELSEIF (canopy_valid_count < canopy_required_count) THEN", 1)[1]
    restart_partial = restart_partial.split("ENDIF", 1)[0]
    assert "DEF_Interception_scheme = 1" in restart_partial
    assert "CALL CoLM_stop" not in restart_partial


def test_singlepoint_site_can_supply_colm2024_structure():
    site = read("mksrfdata/MOD_SingleSrfdata.F90")
    htop = read("mkinidata/MOD_HtopReadin.F90")

    for name in ("SITE_ncd", "SITE_ncw", "SITE_bcw", "SITE_ncd_pfts", "SITE_ncw_pfts", "SITE_bcw_pfts"):
        assert name in site
    for name in ("'ncd'", "'ncw'", "'bcw'", "'ncd_pfts'", "'ncw_pfts'", "'bcw_pfts'"):
        assert name in site
    assert "ncd(:) = SITE_ncd" in htop
    assert "ncd_p = pack(SITE_ncd_pfts, SITE_pctpfts > 0.)" in htop
    assert "SinglePoint has no valid canopy structure" in htop
    assert "SinglePoint lacks valid ncd_pfts/ncw_pfts/bcw_pfts" in htop
    assert "SITE_ncd < 1000._r8" in htop
    assert "SITE_ncw < 1000._r8" in htop
    assert "SITE_bcw < 1000._r8" in htop
    assert "SITE_htop < 1000._r8" in htop
    assert "SITE_htop_pfts < 1000._r8" in htop


def test_restart_validation_agrees_across_mpi_and_handles_empty_ranks(tmp_path):
    """Compile the shipped validation blocks: master arrays are intentionally absent."""
    import os
    import shutil
    import subprocess
    import pytest

    fc, launcher = shutil.which('mpif90'), shutil.which('mpiexec')
    if not fc or not launcher:
        pytest.skip('MPI compiler/runtime unavailable')
    source = read('main/MOD_Vars_TimeInvariants.F90')
    blocks = []
    for name in ('READ_PFTimeInvariants', 'READ_TimeInvariants'):
        body = source.split(f'SUBROUTINE {name} ', 1)[1].split(f'END SUBROUTINE {name}', 1)[0]
        block = body[body.index('            canopy_counts = 0'):]
        block = block.split('\n         ENDIF\n', 1)[0]
        blocks.append(block)
    code = '''
program validate
 use mpi
 use, intrinsic :: ieee_arithmetic
 implicit none
 integer,parameter :: r8=kind(1.d0)
 integer :: rank,p_err,p_comm_glb,n,mode,which,DEF_Interception_scheme
 integer :: canopy_counts(2),canopy_required_count,canopy_valid_count,icanopy
 integer,allocatable :: pftclass(:),patchclass(:),patchtype(:)
 real(r8),allocatable :: ncd_p(:),ncw_p(:),bcw_p(:),htop_p(:),ncd(:),ncw(:),bcw(:),htop(:)
 logical :: p_is_worker,p_is_master
 character(10) :: arg
 call mpi_init(p_err)
 p_comm_glb=MPI_COMM_WORLD
 call mpi_comm_rank(p_comm_glb,rank,p_err)
 p_is_master=rank==0; p_is_worker=rank/=0
 call get_command_argument(1,arg); read(arg,*) mode
 call get_command_argument(2,arg); read(arg,*) which
 DEF_Interception_scheme=8
 if(p_is_worker) then
 n=1
 if(rank==2) n=0
 allocate(pftclass(n),patchclass(n),patchtype(n),ncd_p(n),ncw_p(n),bcw_p(n),htop_p(n), &
          ncd(n),ncw(n),bcw(n),htop(n))
 pftclass=1; patchclass=1; patchtype=0
 ncd=4; ncw=4; bcw=-1.e36_r8; htop=10 ! unused broadleaf metric may be missing
 if(mode==1) ncd=-1.e36_r8 ! all workers lack required structure -> globally select 2014
 if(mode==2 .and. rank==3) ncd=ieee_value(0._r8,ieee_quiet_nan) ! partial invalid
 ncd_p=ncd; ncw_p=ncw; bcw_p=bcw; htop_p=htop
 endif
 if(which==1) then
''' + blocks[0] + '\n else\n' + blocks[1] + '''
 endif
 if(mode==0 .and. DEF_Interception_scheme/=8) error stop 'valid structure downgraded'
 if(mode==1 .and. DEF_Interception_scheme/=1) error stop 'rank-divergent fallback'
 if(mode==2 .and. which==2 .and. DEF_Interception_scheme/=1) error stop 'partial patch fallback diverged'
 call mpi_barrier(p_comm_glb,p_err)
 if(rank==0) print *, 'STRUCTURE_OK'
 call mpi_finalize(p_err)
 contains
 subroutine CoLM_stop()
 print *, 'STRUCTURE_INVALID'
 call mpi_abort(p_comm_glb,7,p_err)
 end subroutine
end program
'''
    driver = tmp_path / 'validate.F90'
    driver.write_text(code)
    exe = tmp_path / 'validate'
    result = subprocess.run([fc, '-cpp', '-DUSEMPI', '-ffree-line-length-0', '-fcheck=all',
                             '-ffpe-trap=invalid,zero,overflow', str(driver), '-o', str(exe)],
                            cwd=tmp_path, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    env = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1',
               OMPI_MCA_rmaps_base_oversubscribe='1')
    for which in (1, 2):
        for mode in (0, 1, 2):
            result = subprocess.run([launcher, '-n', '4', str(exe), str(mode), str(which)],
                                    cwd=tmp_path, env=env, capture_output=True, text=True, timeout=30)
            output = result.stdout + result.stderr
            if mode == 2 and which == 1:
                assert result.returncode != 0 and 'STRUCTURE_INVALID' in output, output
                assert 'SIGFPE' not in output, output
            else:
                assert result.returncode == 0 and 'STRUCTURE_OK' in output, output
                if mode == 2:
                    assert output.count('Downgrading the whole run to CoLM2014 interception') == 1, output
