"""Executable regression for the actual conservative coupling arithmetic."""
import pathlib
import shutil
import subprocess
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]

def test_conservative_credits():
    compiler = shutil.which('gfortran')
    assert compiler, 'gfortran is required for the CaMa conservation check'
    with tempfile.TemporaryDirectory() as tmp:
        driver = pathlib.Path(tmp) / 'check.f90'
        driver.write_text('''program check
use cmf_coupling_budget_mod
implicit none
integer :: x(2,2), y(2,2)
real(8) :: w(2,2), area(2,1), v(2), fa(2), gv(2,1), gf(2,1), runoff(2,1), qr(2)
real(8) :: e(2,1), inf(2,1), er(2), ir(2), storage(2)
x(1,:)=[1,2]; x(2,:)=[2,0]; y=1; y(2,2)=0
w(1,:)=[1d0,3d0]; w(2,:)=[2d0,0d0]; area(:,1)=[10d0,20d0]
v=[8d0,12d0]; fa=[2d0,4d0]
call budget_init(x,y,w,area)
call budget_publish(v,fa,gv,gf)
if(abs(sum(gv*area)-20d0)>1d-12) stop 1
! A partially covered cell supplies volume, never an extrapolated density.
runoff(:,1)=[1d0,2d0]
call budget_runoff(runoff,qr)
if(abs(sum(qr)-3d0)>1d-12) stop 2
! All of grid 1, half of grid 2 consumed: total 2 + 9 = 11.
e(:,1)=[1d0,4d0]; inf(:,1)=[1d0,5d0]; storage=v
call budget_debit(e,inf,storage,er,ir)
if(abs(sum(storage)-9d0)>1d-12) stop 3
if(abs(sum(er)+sum(ir)-11d0)>1d-12) stop 4
if(any(storage<0)) stop 5
if(maxval(abs(storage-[3d0,6d0]))>1d-12) stop 6
! Zero-water cells must neither divide by zero nor manufacture credit.
call budget_publish([0d0,0d0],fa,gv,gf)
e=0; inf=0; storage=0
call budget_debit(e,inf,storage,er,ir)
if(any(storage/=0)) stop 7
end program
''')
        subprocess.run([compiler, '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                        str(ROOT/'extends/CaMa/src/cmf_coupling_budget_mod.F90'),
                        str(driver), '-o', str(pathlib.Path(tmp)/'check')], cwd=tmp, check=True)
        subprocess.run([str(pathlib.Path(tmp)/'check')], check=True)
