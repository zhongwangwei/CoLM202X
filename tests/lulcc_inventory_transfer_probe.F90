program check_full
use MOD_Precision
use MOD_Namelist, only: DEF_USE_PFT, DEF_SOLO_PFT, DEF_FAST_PC
use MOD_Const_LC, only: patchtypes
use MOD_Lulcc_Driver, only: lulcc_inventory_trace, lulcc_check_inventory_transfer, lulcc_patch_areas
use MOD_Mesh, only: mesh
use MOD_Pixel, only: pixel
use MOD_Pixelset, only: pixelset_type
use MOD_Utils, only: areaquad
use MOD_Tracer_Defs, only: ntracers, tracers, STATE_OWNER_GENERIC_WATER, FAMILY_ISOTOPE
use MOD_Tracer_Vars, only: allocate_Tracer_Vars, save_land_tracer_lulcc_state, &
 remap_land_tracer_lulcc_state,trc_wdsrf
use MOD_Tracer_Reactive_Methane_Const, only: DEF_METHANE
use MOD_Tracer_Reactive_Methane_State, only: allocate_methane_state, &
 save_methane_lulcc_state,remap_methane_lulcc_state,totcol_methane
use MOD_Tracer_Reactive_Methane_Microbes, only: allocate_methane_microbes_state, &
 save_methane_microbes_lulcc_state,remap_methane_microbes_lulcc_state,B_methanogen,B_methanogen_comp
use MOD_Vars_TimeInvariants, only: patchtype
use MOD_Vars_Global, only: dz_soi
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
implicit none
real(r8)::raw(3,0:17), oldarea(3),newarea(3),before,after
real(r8),allocatable::mapped(:,:), physical_area(:)
type(pixelset_type):: patch_geometry
integer:: oldclass(3), newclass(3)
character(len=32):: mode
integer*8:: oldelm(3), newelm(3)
allocate(mesh(1),pixel%lat_s(1),pixel%lat_n(1),pixel%lon_w(2),pixel%lon_e(2))
mesh(1)%npxl=2
allocate(mesh(1)%ilat(2),mesh(1)%ilon(2))
mesh(1)%ilat=[1,1]; mesh(1)%ilon=[1,2]
pixel%lat_s=0._r8; pixel%lat_n=1._r8
pixel%lon_w=[0._r8,1._r8]; pixel%lon_e=[1._r8,2._r8]
patch_geometry%nset=3
allocate(patch_geometry%ielm(3),patch_geometry%ipxstt(3),patch_geometry%ipxend(3), &
 patch_geometry%pctshared(3))
patch_geometry%ielm=1
patch_geometry%ipxstt=[1,2,-1]; patch_geometry%ipxend=[1,2,-1]
patch_geometry%has_shared=.true.;patch_geometry%pctshared=[1._r8,.25_r8,1._r8]
call lulcc_patch_areas(patch_geometry,physical_area)
if(abs(physical_area(1)-1.e6_r8*areaquad(0._r8,1._r8,0._r8,1._r8))>1.e-6_r8) error stop 15
if(abs(physical_area(2)-.25_r8*physical_area(1))>1.e-6_r8) error stop 16
if(physical_area(3)/=0._r8) error stop 17
print *, 'geometry',physical_area
block
 type(pixelset_type) :: empty_shared
 real(r8), allocatable :: empty_area(:)
 integer :: empty_class(0)
 integer*8 :: empty_element(0)
 real(r8) :: empty_trace(0,0:17)
 empty_shared%nset=0
 empty_shared%has_shared=.true.
 call lulcc_patch_areas(empty_shared,empty_area)
 if(size(empty_area)/=0) error stop 20
 call lulcc_check_inventory_transfer(empty_class,empty_element,empty_area, &
    empty_class,empty_element,empty_area,empty_trace)
end block
call get_command_argument(1,mode)
if (mode=='mismatch' .or. mode=='missing_donor' .or. mode=='missing_target' .or. &
    mode=='footprint' .or. mode=='bad_row' .or. mode=='nan_area') then
 block
 real(r8)::badtrace(2,0:17)
 badtrace=0._r8
 badtrace(1,1)=.2_r8;badtrace(1,12)=.8_r8
 badtrace(2,1)=.6_r8;badtrace(2,12)=.4_r8
 select case(mode)
 case('mismatch')
   call lulcc_check_inventory_transfer([1,12],[9_8,9_8],[80._r8,20._r8], &
     [1,12],[9_8,9_8],[50._r8,50._r8],badtrace)
 case('missing_donor')
   badtrace=0._r8;badtrace(:,12)=1._r8
   call lulcc_check_inventory_transfer([1],[9_8],[100._r8], &
     [12,12],[9_8,9_8],[50._r8,50._r8],badtrace)
 case('missing_target')
   badtrace=0._r8;badtrace(:,1)=1._r8
   call lulcc_check_inventory_transfer([0],[9_8],[100._r8], &
     [1,1],[9_8,9_8],[50._r8,50._r8],badtrace)
 case('footprint')
   badtrace=0._r8;badtrace(:,1)=1._r8
   call lulcc_check_inventory_transfer([1],[9_8],[100._r8], &
     [1,1],[9_8,9_8],[50._r8,60._r8],badtrace)
 case('bad_row')
   badtrace=0._r8;badtrace(:,1)=.5_r8
   call lulcc_check_inventory_transfer([1],[9_8],[100._r8], &
     [1,1],[9_8,9_8],[50._r8,50._r8],badtrace)
 case('nan_area')
   badtrace=0._r8;badtrace(:,1)=1._r8
   call lulcc_check_inventory_transfer([1],[9_8],[ieee_value(0._r8,ieee_quiet_nan)], &
     [1,1],[9_8,9_8],[50._r8,50._r8],badtrace)
 end select
 error stop 96
 end block
endif
DEF_USE_PFT=.true.
DEF_SOLO_PFT=.false.
DEF_FAST_PC=.true.
patchtypes=0
patchtypes(11)=2
patchtypes(13)=1
patchtypes(15)=3
patchtypes(17)=4
raw=0
raw(1,0)=1
raw(2,4)=.2_r8
raw(2,12)=.2_r8
raw(2,14)=.6_r8
raw(3,4)=1
call lulcc_inventory_trace(raw,mapped)
if(abs(mapped(2,1)-.2_r8)>1.e-12_r8 .or. &
   abs(mapped(2,12)-.8_r8)>1.e-12_r8) error stop 18
DEF_SOLO_PFT=.true.; DEF_FAST_PC=.false.
call lulcc_inventory_trace(raw,mapped)
if(abs(mapped(2,4)-.2_r8)>1.e-12_r8 .or. &
   abs(mapped(2,14)-.6_r8)>1.e-12_r8 .or. mapped(2,1)/=0._r8) error stop 19
DEF_SOLO_PFT=.false.; DEF_FAST_PC=.true.
call lulcc_inventory_trace(raw,mapped)
oldclass=[1,0,12]; oldelm=[9_8,4_8,9_8]; oldarea=[60._r8,20._r8,40._r8]
newclass=[0,12,1]; newelm=[4_8,9_8,9_8]; newarea=[20._r8,50._r8,50._r8]
call lulcc_check_inventory_transfer(oldclass,oldelm,oldarea,newclass,newelm,newarea,mapped)
ntracers=1
allocate(tracers(1))
tracers(1)%name='H2_18O'
tracers(1)%category='isotope'
tracers(1)%state_owner=STATE_OWNER_GENERIC_WATER
tracers(1)%family_id=FAMILY_ISOTOPE
tracers(1)%init_delta=0._r8
call allocate_Tracer_Vars(3,-1,1)
trc_wdsrf(1,:)=[10._r8,30._r8,20._r8]
before=sum(oldarea*trc_wdsrf(1,:))
call save_land_tracer_lulcc_state()
call remap_land_tracer_lulcc_state(newclass,newelm,oldclass,oldelm,mapped,newarea,oldarea)
after=sum(newarea*trc_wdsrf(1,:))
print *, 'generic',before,after,trc_wdsrf(1,:)
if (abs(after-before)>1.e-8_r8) error stop 11
if (any(abs(trc_wdsrf(1,:)-[30._r8,18._r8,10._r8])>1.e-8_r8)) error stop 12
allocate(patchtype(3)); patchtype=0
dz_soi=1._r8
DEF_METHANE%use_microbial_pools=.true.
call allocate_methane_state(3)
totcol_methane=[10._r8,30._r8,20._r8]
call save_methane_lulcc_state()
call allocate_methane_microbes_state(3)
B_methanogen=0._r8
B_methanogen(1,:)=[10._r8,30._r8,20._r8]
B_methanogen_comp=0._r8
B_methanogen_comp(1,1,:)=[10._r8,30._r8,20._r8]
B_methanogen_comp(1,2,:)=[10._r8,30._r8,20._r8]
call save_methane_microbes_lulcc_state()
call remap_methane_lulcc_state(newclass,newelm,oldclass,oldelm,mapped,newarea,oldarea)
print *, 'methane',sum(newarea*totcol_methane),totcol_methane
if(abs(sum(newarea*totcol_methane)-2000._r8)>1.e-8_r8) error stop 13
call remap_methane_microbes_lulcc_state(newclass,newelm,oldclass,oldelm,mapped,newarea,oldarea)
print *, 'microbes',sum(newarea*B_methanogen(1,:)),B_methanogen(1,:)
if(abs(sum(newarea*B_methanogen(1,:))-2000._r8)>1.e-8_r8) error stop 14
end program
