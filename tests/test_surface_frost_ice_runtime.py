"""Compiled water/enthalpy/tracer checks for late soil-frost ice relocation."""

from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
FORTRAN = r"""
module MOD_Precision
  integer, parameter :: r8 = selected_real_kind(15)
end module
module MOD_Namelist
  logical :: DEF_USE_SNICAR = .true.
  logical :: DEF_USE_VariablySaturatedFlow = .false.
end module
module MOD_Const_Physical
  use MOD_Precision
  real(r8), parameter :: denice=917._r8, denh2o=1000._r8
  real(r8), parameter :: cpice=2117.27_r8, cpliq=4188._r8, hfus=333600._r8, tfrz=273.16_r8
end module
program check_frost
  use MOD_Precision
  use MOD_Const_Physical
  use MOD_NewSnow, only: relocate_soil_frost_ice
  use MOD_SnowLayersCombineDivide, only: snowlayerscombine
  implicit none
  integer, parameter :: ms=-3
  real(r8) :: zi(ms:0), z(ms+1:1), dz(ms+1:1), t(ms+1:1)
  real(r8) :: wl(ms+1:1), wi(ms+1:1), fi(ms+1:1), fr(ms+1:0), rds(ms+1:0)
  real(r8) :: bcpho(ms+1:0),bcphi(ms+1:0),ocpho(ms+1:0),ocphi(ms+1:0)
  real(r8) :: dst1(ms+1:0),dst2(ms+1:0),dst3(ms+1:0),dst4(ms+1:0)
  real(r8) :: ti(2,ms+1:1),tl(2,ms+1:1),ts(2,ms+1:1),tc(2)
  real(r8) :: scv, dp, m0, h0, isotope0(2), cold, expected, cap
  integer :: snl, im(ms+1:1)
  cap=denice*0.4_r8*0.05_r8

  ! Thin snow + frost crosses the explicit-layer threshold.  Old thin snow
  ! isotope and the post-transport soil ice isotope must both survive.
  call reset()
  scv=2._r8; dp=0.006_r8; tc=[0.5_r8,0.05_r8]
  wi(1)=cap+5._r8; ti(:,1)=[4.668_r8,0.4668_r8]
  t(1)=280._r8
  m0=wi(1)+scv; isotope0=ti(:,1)+tc
  h0=cpice*m0*t(1)
  call relocate(.true.)
  call require(snl==-1 .and. abs(dp-(0.006_r8+5._r8/denice))<1.e-12_r8)
  call require(abs(wi(1)-cap)<1.e-11_r8 .and. abs(wi(0)-7._r8)<1.e-12_r8)
  call require(abs(wi(1)+scv-m0)<1.e-11_r8)
  call require(all(abs(ti(:,1)+ti(:,0)-isotope0)<1.e-12_r8))
  call require(all(abs(ti(:,0)-[1.5_r8,0.15_r8])<1.e-12_r8))
  call require(all(tc==0._r8) .and. abs(cpice*(wi(1)*t(1)+wi(0)*t(0))-h0)<1.e-7_r8)
  call require(im(0)==0 .and. fi(0)==1._r8 .and. fr(0)==0._r8)
  call require(zi(0)==0._r8 .and. abs(zi(-1)+dp)<1.e-12_r8)
  call require(rds(0)>0._r8 .and. bcpho(0)==0._r8 .and. dst4(0)==0._r8)
  m0=wi(1)+scv; isotope0=ti(:,1)+ti(:,0)
  call relocate(.true.) ! zero excess must not duplicate ice or isotope
  call require(abs(wi(1)+scv-m0)<1.e-12_r8)
  call require(all(abs(ti(:,1)+ti(:,0)-isotope0)<1.e-12_r8))
  ! A restart-like binary round trip of every field touched by the helper
  ! reproduces the same next-step no-excess state (no hidden helper memory).
  open(10,status='scratch',form='unformatted')
  write(10) snl,zi,z,dz,t,wl,wi,fi,im,fr,rds,scv,dp, &
       bcpho,bcphi,ocpho,ocphi,dst1,dst2,dst3,dst4,ti,tl,ts,tc
  rewind(10)
  call reset()
  read(10) snl,zi,z,dz,t,wl,wi,fi,im,fr,rds,scv,dp, &
       bcpho,bcphi,ocpho,ocphi,dst1,dst2,dst3,dst4,ti,tl,ts,tc
  close(10)
  call relocate(.true.)
  call require(snl==-1 .and. abs(wi(1)+scv-m0)<1.e-12_r8)
  call require(all(abs(ti(:,1)+ti(:,0)-isotope0)<1.e-12_r8))

  ! Subthreshold ice remains in thin-snow scv; no layer is fabricated.
  call reset()
  wi(1)=cap+2._r8; ti(:,1)=[2._r8,0.2_r8]
  call relocate(.true.)
  call require(snl==0 .and. abs(scv-2._r8)<1.e-12_r8)
  call require(abs(dp-2._r8/denice)<1.e-12_r8)
  call require(all(abs(ti(:,1)+tc-[2._r8,0.2_r8])<1.e-12_r8))

  ! Existing snow: preserve liquid, aerosol, geometry and sensible enthalpy.
  call reset()
  snl=-1; wi(0)=3._r8; wl(0)=1._r8; scv=4._r8; dp=0.02_r8
  dz(0)=dp; zi(0)=0._r8; zi(-1)=-dp; t(0)=260._r8
  wi(1)=cap+5._r8; ti(:,1)=[5._r8,0.5_r8]
  ti(:,0)=[0.6_r8,0.06_r8]; tl(:,0)=[0.1_r8,0.01_r8]
  bcpho(0)=0.25_r8; dst4(0)=0.125_r8
  t(1)=280._r8
  h0=(cpice*wi(0)+cpliq*wl(0))*t(0)+cpice*wi(1)*t(1)
  isotope0=ti(:,0)+ti(:,1)
  call relocate(.true.)
  call require(snl==-1 .and. abs(wi(0)-8._r8)<1.e-12_r8)
  call require(abs(scv-9._r8)<1.e-12_r8 .and. abs(wl(0)-1._r8)<1.e-12_r8)
  call require(all(abs(ti(:,0)+ti(:,1)-isotope0)<1.e-12_r8))
  call require(all(abs(tl(:,0)-[0.1_r8,0.01_r8])<1.e-12_r8))
  call require(abs((cpice*wi(0)+cpliq*wl(0))*t(0)+cpice*wi(1)*t(1)-h0)<1.e-7_r8)
  call require(t(0)>260._r8 .and. t(0)<280._r8)
  call require(abs(zi(-1)+dp)<1.e-12_r8)
  call require(bcpho(0)==0.25_r8 .and. dst4(0)==0.125_r8)

  ! Actual snow combine can collapse a sub-1 cm explicit layer to thin scv.
  ! The following frost then reopens the layer without losing its ice tracer.
  call reset()
  snl=-1; wi(0)=0.2_r8; dz(0)=0.005_r8; zi(0)=0._r8
  scv=0.2_r8; dp=0.005_r8; ti(:,0)=[0.04_r8,0.004_r8]
  wi(1)=cap+5._r8; ti(:,1)=[5._r8,0.5_r8]
#ifdef TRACER
  call snowlayerscombine(ms+1,snl,z,dz,zi,wl,wi,t,scv,dp, &
       trc_wliq=tl,trc_wice=ti,trc_solid=ts,trc_scv=tc)
#else
  call snowlayerscombine(ms+1,snl,z,dz,zi,wl,wi,t,scv,dp)
  tc=[0.04_r8,0.004_r8]; ti(:,0)=0._r8
#endif
  call require(snl==0 .and. abs(scv-0.2_r8)<1.e-12_r8)
  call relocate(.true.)
  call require(snl==-1 .and. abs(wi(0)-5.2_r8)<1.e-12_r8)
  call require(all(abs(ti(:,0)+ti(:,1)-[5.04_r8,0.504_r8])<1.e-12_r8))

  ! A wet thin layer returns its *liquid* to soil in the existing combine
  ! routine; this ice-only helper must neither erase nor mislabel that mass.
  call reset()
  snl=-1; wi(0)=0.2_r8; wl(0)=0.3_r8; dz(0)=0.005_r8; zi(0)=0._r8
  scv=0.5_r8; dp=0.005_r8; ti(:,0)=[0.04_r8,0.004_r8]
  tl(:,0)=[0.06_r8,0.006_r8]
  wi(1)=cap+5._r8
#ifdef TRACER
  call snowlayerscombine(ms+1,snl,z,dz,zi,wl,wi,t,scv,dp, &
       trc_wliq=tl,trc_wice=ti,trc_solid=ts,trc_scv=tc)
#else
  call snowlayerscombine(ms+1,snl,z,dz,zi,wl,wi,t,scv,dp)
  tc=[0.04_r8,0.004_r8]; tl(:,1)=[0.06_r8,0.006_r8]
#endif
  call require(snl==0 .and. abs(wl(1)-0.3_r8)<1.e-12_r8)
  call relocate(.true.)
  call require(abs(wl(1)-0.3_r8)<1.e-12_r8)
  call require(all(abs(tl(:,1)-[0.06_r8,0.006_r8])<1.e-12_r8))

  ! The helper is usable in a build without TRACER as well.
  call reset()
  wi(1)=cap+2._r8
  call relocate(.false.)
  call require(snl==0 .and. abs(wi(1)-cap)<1.e-11_r8)

contains
  subroutine reset()
    snl=0; zi=999._r8; z=0._r8; dz=0._r8; dz(1)=0.05_r8
    t=270._r8; wl=0._r8; wi=0._r8; fi=0._r8; im=0; fr=0._r8
    rds=-1._r8; bcpho=0._r8; bcphi=0._r8; ocpho=0._r8; ocphi=0._r8
    dst1=0._r8; dst2=0._r8; dst3=0._r8; dst4=0._r8
    ti=0._r8; tl=0._r8; ts=0._r8; tc=0._r8; scv=0._r8; dp=0._r8
  end subroutine
  subroutine relocate(with_tracer)
    logical, intent(in) :: with_tracer
    if (with_tracer) then
      call relocate_soil_frost_ice(ms,0.4_r8,snl,zi,z,dz,t,wl,wi,fi,im,fr,rds,scv,dp, &
        bcpho,bcphi,ocpho,ocphi,dst1,dst2,dst3,dst4,ti,tl,ts,tc)
    else
      call relocate_soil_frost_ice(ms,0.4_r8,snl,zi,z,dz,t,wl,wi,fi,im,fr,rds,scv,dp, &
        bcpho,bcphi,ocpho,ocphi,dst1,dst2,dst3,dst4)
    endif
  end subroutine
  subroutine require(ok)
    logical, intent(in) :: ok
    if (.not. ok) error stop 'surface frost ice conservation check failed'
  end subroutine
end program
"""


@pytest.mark.parametrize("tracer", [False, True])
def test_surface_frost_ice_runtime(tmp_path, tracer):
    compiler = shutil.which("gfortran")
    if compiler is None:
        pytest.skip("gfortran unavailable")
    harness = tmp_path / "check.f90"
    stubs, program = FORTRAN.split("program check_frost", 1)
    (tmp_path / "stubs.f90").write_text(stubs)
    harness.write_text("program check_frost" + program)
    binary = tmp_path / "check"
    cmd = [compiler, "-fdefault-real-8", "-ffree-line-length-0", "-fcheck=all"]
    if tracer:
        cmd.append("-DTRACER")
    (tmp_path / "define.h").write_text("\n")
    cmd += ["-I", str(tmp_path)]
    cmd += ["-cpp", str(tmp_path / "stubs.f90"), str(ROOT / "main/MOD_NewSnow.F90"),
            str(ROOT / "main/MOD_SnowLayersCombineDivide.F90"),
            str(harness), "-o", str(binary)]
    subprocess.run(cmd, cwd=tmp_path, check=True, capture_output=True, text=True)
    subprocess.run([str(binary)], cwd=tmp_path, check=True, capture_output=True, text=True)
