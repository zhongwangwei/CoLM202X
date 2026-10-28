#!/usr/bin/env python3
"""Compile the actual glacier/lake mixing blocks; host physics is a closed pool.

Optional --before PATH reproduces the repeated NSS inventory bug in an old file.
The normal check also tests refusal of orphan inventories without safe carriers.
"""
import pathlib
import re
import subprocess
import tempfile
import sys
import hashlib

root = pathlib.Path(__file__).resolve().parents[2]
before = len(sys.argv) > 1 and sys.argv[1] == '--before'
text = pathlib.Path(sys.argv[2]).read_text() if before else (root / 'main/MATERIAL/GENERIC_WATER/MOD_MT_SpecialPatches.F90').read_text()

fixture = '''module fixture
implicit none
integer, parameter :: r8=kind(1.d0), ntracers=1
real(r8), parameter :: trc_water_min_for_ratio=1.d-12, trc_tiny=1.d-20
real(r8) :: DEF_ISO_SUBL_SKIN_MM=1.d0, mobile
logical :: isotope=.true., fixed=.false.
type tracer
real(r8) :: init_delta=-20.d0
end type
type(tracer) :: tracers(1)
real(r8) :: trc_leaf_iso_storage(1,1), trc_leaf_water_moles(1,1)
real(r8) :: trc_leaf_delta_e(1,1), trc_leaf_delta_b(1,1), trc_leaf_peclet(1,1)
real(r8) :: trc_storage_beg(1,1), trc_surface_residue(1,1), trc_subsurface_residue(1,1)
real(r8) :: trc_surface_solid(1,1),trc_rnof_step(1,1),a_trc_precip(1,1),a_trc_rsur(1,1),a_trc_rnof(1,1)
real(r8) :: trc_wdsrf(1,1),trc_wliq_soisno(1,1,1)
logical, allocatable :: trc_runtime_forced(:)
integer, parameter :: TRC_EVAP_KIND_SOILEVAP=1, TRC_EVAP_KIND_SUBL=2
contains
logical function tracer_is_isotope(i)
integer :: i
tracer_is_isotope=isotope
end function
logical function tracer_can_use_fixed_signature(i)
integer :: i
tracer_can_use_fixed_signature=fixed
end function
logical function tracer_fractionation_active(i)
integer :: i
tracer_fractionation_active=.false.
end function
logical function tracer_is_nonvolatile_solute(i)
integer :: i
tracer_is_nonvolatile_solute=.false.
end function
logical function mt_has_dissolved_limit(i)
integer :: i
mt_has_dissolved_limit=.false.
end function
real(r8) function tracer_init_water_ratio(i)
integer :: i
tracer_init_water_ratio=2.d0
end function
real(r8) function solid_inventory(i,p,lb,n)
integer :: i,p,lb,n
solid_inventory=0.d0
end function
real(r8) function tracer_forcing_rain_value(i,p)
integer :: i,p
tracer_forcing_rain_value=2.d0
end function
real(r8) function tracer_forcing_snow_value(i,p)
integer :: i,p
tracer_forcing_snow_value=2.d0
end function
real(r8) function tracer_forcing_vapor_value(i,p)
integer :: i,p
tracer_forcing_vapor_value=2.d0
end function
real(r8) function tracer_equilibrium_deposition_ratio(i,r,t,ice)
integer :: i
real(r8) :: r,t
logical :: ice
tracer_equilibrium_deposition_ratio=r
end function
real(r8) function tracer_ratio_cap(i)
integer :: i
tracer_ratio_cap=0.d0
end function
subroutine glacier_evap_affine_for()
end subroutine
subroutine waterbody_evap_affine_for()
end subroutine
real(r8) function tracer_atmospheric_tracer_loss(m,w,e,t,ice,aff,tiny,cap,nv,skin_mass)
real(r8) :: m,w,e,t,tiny,cap
real(r8),optional :: skin_mass
logical :: ice,nv
external :: aff
tracer_atmospheric_tracer_loss=0.d0
if (e /= 0.d0) error stop 'fixture expects no evaporation'
end function
subroutine mt_equilibrate_dissolved(i,w,m,s)
integer :: i
real(r8) :: w,m,s
end subroutine
subroutine tracer_book_evap_loss(i,p,m,w,k)
integer :: i,p,k
real(r8) :: m,w
end subroutine
subroutine sync_tracer_patch_ratio(i,p,snl,maxsnl,nl,wl,wi,wa,ws,sc,r)
integer :: i,p,snl,maxsnl,nl
real(r8) :: wl(:),wi(:),wa,ws,sc,r
mobile=(sum(wl)+sum(wi)+wa+ws+sc)*r
end subroutine
subroutine CoLM_stop(msg)
character(*) :: msg
! The failing path must leave the signed inventory intact.
print *, 'PRESERVED_LEAF', trc_leaf_iso_storage(1,1)
print *, msg
stop 7
end subroutine
'''
for kind in ('glacier', 'waterbody'):
    routine = re.search(rf'   SUBROUTINE mt_{kind}_patch\(.*?   END SUBROUTINE mt_{kind}_patch', text, re.S).group()
    decl = routine[routine.index('      integer  :: itrc'):routine.index('      IF (ntracers <= 0) RETURN')]
    body = routine[routine.index('         trc_rnof_step'):routine.index('      ENDDO', routine.index('         trc_rnof_step'))]
    fixture += f'''subroutine step_{kind}(water)
real(r8),intent(in) :: water
integer,parameter :: ipatch=1,maxsnl=0,nl_soil=1,snl=0
real(r8) :: wliq_soisno(1),wice_soisno(1),wa,wdsrf,scv,t_grnd,forc_rain,deltim
{decl}
itrc=1
'''
    if kind == 'glacier':
        fixture += 'snl_trc=0\nrain_mass=0.d0\n'
    else:
        fixture += 'atm_precip_mass=0.d0\ndeficit_mass=0.d0\n'
    fixture += '''wliq_soisno=water
wice_soisno=0.d0
wa=0.d0
wdsrf=0.d0
scv=0.d0
t_grnd=280.d0
forc_rain=0.d0
deltim=1.d0
water_end=water
water_beg=water
water_input=0.d0
precip_mass=0.d0
rnof_mass=0.d0
evap_mass=0.d0
evap_liq_mass=0.d0
evap_ice_mass=0.d0
dep_liq_mass=0.d0
dep_ice_mass=0.d0
trc_storage_beg=mobile+trc_leaf_iso_storage
'''+body+f'end subroutine step_{kind}\n'
fixture += '''end module
program check
use fixture
implicit none
character(40) :: arg
integer :: kind, step
real(r8) :: water,leaf,expected
call get_command_argument(1,arg)
read(arg,*) kind
call get_command_argument(2,arg)
read(arg,*) leaf
call get_command_argument(3,arg)
read(arg,*) water
call get_command_argument(4,arg)
read(arg,*) fixed
call get_command_argument(5,arg)
read(arg,*) isotope
mobile=2.d0*water
trc_leaf_iso_storage=leaf
trc_leaf_water_moles=10.d0
trc_leaf_delta_e=123.d0
trc_leaf_delta_b=456.d0
trc_leaf_peclet=0.3d0
trc_surface_residue=0.d0
trc_subsurface_residue=0.d0
trc_surface_solid=0.d0
a_trc_precip=0.d0
a_trc_rsur=0.d0
a_trc_rnof=0.d0
if (leaf==0.d0 .and. isotope) then
trc_leaf_water_moles=0.d0
trc_leaf_delta_e=tracers(1)%init_delta
trc_leaf_delta_b=tracers(1)%init_delta
trc_leaf_peclet=1.d0
endif
expected=mobile+leaf
do step=1,2
 if(kind==1) then
 call step_glacier(water)
 else
 call step_waterbody(water)
 endif
 print '(I1,4ES25.16)',step,mobile+trc_leaf_iso_storage(1,1)-expected, &
 trc_leaf_iso_storage(1,1),trc_leaf_water_moles(1,1),mobile
enddo
end program
'''
with tempfile.TemporaryDirectory(prefix='leaf-inventory-') as tmp:
    src=pathlib.Path(tmp)/'check.f90'
    src.write_text(fixture)
    exe=pathlib.Path(tmp)/'check'
    subprocess.run(['gfortran','-O0','-fcheck=all','-ffree-line-length-none',str(src),'-o',str(exe)],check=True,cwd=tmp)
    count=0
    unaffected=[]
    for kind in (1,2):
        for leaf in (-0.125,0.,0.125):
            for fixed in (False,True):
                out=subprocess.run([str(exe),str(kind),str(leaf),'10',str(fixed), 'True'],text=True,capture_output=True)
                if not before and fixed and leaf:
                    assert out.returncode==7 and 'requires mixed-signature' in out.stdout
                    assert float(out.stdout.splitlines()[0].split()[-1])==leaf
                    count+=1
                    continue
                out.check_returncode()
                rows=[[float(v) for v in line.split()] for line in out.stdout.splitlines()]
                if leaf==0:
                    unaffected.append(out.stdout)
                if before:
                    if not fixed and leaf:
                        assert rows[0][1] == leaf and rows[1][1] == 2*leaf,rows
                    continue
                assert all(abs(row[1])<1e-14 and row[2]==0 and row[3]==0 for row in rows),rows
                count+=1
        if not before:
            for water,leaf in ((0.,0.125),(1e-13,0.125),(0.,-0.125),(1e-13,-0.125),(1.,-3.)):
                out=subprocess.run([str(exe),str(kind),str(leaf),str(water),'False','True'],text=True,capture_output=True)
                assert out.returncode==7 and 'PRESERVED_LEAF' in out.stdout,(out.returncode,out.stdout,out.stderr)
                assert float(out.stdout.splitlines()[0].split()[-1])==leaf
                count+=1
        out=subprocess.run([str(exe),str(kind),'0','10','True','False'],text=True,capture_output=True,check=True)
        rows=[[float(v) for v in line.split()] for line in out.stdout.splitlines()]
        assert all(row[1]==0 and row[3]==10 for row in rows),rows
        unaffected.append(out.stdout)
        count+=1
if not before:
    vars_text=(root/'main/MATERIAL/MOD_MT_Vars.F90').read_text()
    names=('compute_lulcc_snapshot_land_water_mass','compute_current_lulcc_land_water_mass',
           'accumulate_lulcc_mass_2d','accumulate_lulcc_mass_3d')
    blocks=[re.search(rf'      SUBROUTINE {name}\(.*?      END SUBROUTINE {name}',vars_text,re.S).group() for name in names]
    pools={pool:int(ndim) for ndim,pool in re.findall(r'CALL accumulate_lulcc_mass_(2|3)d\((\w+), area, total\)', ''.join(blocks[:2]))}
    declarations='\n'.join(f'real(r8),allocatable :: {pool}('+','.join(':'*ndim)+')' for pool,ndim in pools.items())
    audit='module audit\ninteger,parameter :: r8=kind(1.d0),ntracers=1\n'+declarations+'\ncontains\n'
    audit+='logical function mt_uses_generic_water(i)\ninteger :: i\nmt_uses_generic_water=.true.\nend function\n'
    audit+='\n'.join(blocks)+"""
end module
program check
use audit
real(r8) :: old(1),new(1)
allocate(lulcc_trc_leaf_iso_storage_old(1,2),trc_leaf_iso_storage(1,1))
lulcc_trc_leaf_iso_storage_old(1,:)=[0.5d0,-0.125d0]
trc_leaf_iso_storage(1,1)=0.25d0
call compute_lulcc_snapshot_land_water_mass([3.d0,2.d0],old)
call compute_current_lulcc_land_water_mass([5.d0],new)
if(old(1)/=1.25d0 .or. new(1)/=old(1)) error stop 'signed leaf audit omission'
! The audit must detect a lost anomaly after remapping, not report a false zero.
trc_leaf_iso_storage=0.d0
call compute_current_lulcc_land_water_mass([5.d0],new)
if(new(1)==old(1)) error stop 'audit missed loss'
end program
"""
    with tempfile.TemporaryDirectory(prefix='leaf-audit-') as tmp:
        src=pathlib.Path(tmp)/'audit.f90';src.write_text(audit)
        exe=pathlib.Path(tmp)/'audit'
        subprocess.run(['gfortran','-O0','-fcheck=all','-ffree-line-length-none',str(src),'-o',str(exe)],check=True,cwd=tmp)
        subprocess.run([str(exe)],check=True)
print('PASS old repeated signed residual reproduced' if before else f'PASS {count} production-block cases + signed area-weighted production LULCC audits')
print('UNAFFECTED closed-pool states SHA256', hashlib.sha256(''.join(unaffected).encode()).hexdigest())
