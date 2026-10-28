from pathlib import Path
import json
targets=json.loads((Path(__file__).parent/"targets.json").read_text())
assert targets, "no trace targets"
spatial=".or. &\n         ".join(f"(abs(lat-({t['lat']}_r8))<0.37501_r8.and.abs(lo-({t['lon']}_r8))<0.37501_r8)" for t in targets)
selection=".or. &\n         ".join(f"(abs(lat-({t['lat']}_r8))<0.37501_r8.and.abs(lo-({t['lon']}_r8))<0.37501_r8.and.abs(idate(2)-{t['doy']})<=1)" for t in targets)
root=Path(__file__).parent/'new'
def edit(p,old,new,count=1):
 p=root/p;s=p.read_text();assert s.count(old)>=count,(p,old);p.write_text(s.replace(old,new,count))
edit('main/MATERIAL/MOD_MT_Defs.F90','   SAVE\n','''   SAVE
   logical :: trace_active=.false.
   integer :: trace_patch=0,trace_date(3)=0,trace_class=0,trace_unit=0,trace_hist_date(3)=0
   real(r8) :: trace_lat=0._r8,trace_lon=0._r8
''')
edit('main/MATERIAL/MOD_MT_Defs.F90','CONTAINS\n','''CONTAINS
   logical FUNCTION trace_selected(lat,lon)
      real(r8),intent(in)::lat,lon
      real(r8)::lo
      lo=modulo(lon+180._r8,360._r8)-180._r8
      trace_selected=TRACE_SPATIAL
   END FUNCTION
   SUBROUTINE trace_context(ip,idate,lat,lon,pc)
      integer,intent(in)::ip,idate(3),pc
      real(r8),intent(in)::lat,lon
      real(r8)::lo
      lo=modulo(lon+180._r8,360._r8)-180._r8
      trace_patch=ip;trace_date=idate;trace_class=pc
      trace_lat=lat;trace_lon=lon
      trace_active=TRACE_SELECTION
   END SUBROUTINE
   SUBROUTINE trace_write(tag,values)
      USE MOD_SPMD_Task,only:p_iam_glb
      character(len=*),intent(in)::tag
      real(r8),intent(in)::values(:)
      character(len=80)::filename
      IF (.not.trace_active) RETURN
      IF(trace_unit==0) THEN
         WRITE(filename,'(A,I0,A)') 'trace_rank_',p_iam_glb,'.txt'
         OPEN(newunit=trace_unit,file=trim(filename),status='replace',action='write')
      ENDIF
      WRITE(trace_unit,'(A,1X,6(I0,1X),*(ES26.17E3,1X))') trim(tag),p_iam_glb,trace_patch, &
         trace_date,trace_class,trace_lat,trace_lon,values
   END SUBROUTINE
'''.replace('TRACE_SPATIAL',spatial).replace('TRACE_SELECTION',selection))
edit('main/CoLMMAIN.F90','   USE MOD_MT_Hist, only: mt_hist_accumulate','   USE MOD_MT_Defs, only: trace_context\n   USE MOD_MT_Hist, only: mt_hist_accumulate')
edit('main/CoLMMAIN.F90','      z_soisno (maxsnl+1:0)', '''#ifdef MaterialTransport
      CALL trace_context(ipatch,idate,patchlatr*180._r8/acos(-1._r8), &
         patchlonr*180._r8/acos(-1._r8),patchclass)
#endif
      z_soisno (maxsnl+1:0)''')
for p in ['main/MATERIAL/FAMILY_ISOTOPE/MOD_MT_Frac.F90','main/MATERIAL/GENERIC_WATER/MOD_MT_EvapLimit.F90','main/MATERIAL/MOD_MT_Vars.F90']:
 edit(p,'   USE MOD_Precision','   USE MOD_MT_Defs, only: trace_write\n   USE MOD_Precision')
# Trace actual affine inputs after computation, without recomputation of production quantities.
edit('main/MATERIAL/FAMILY_ISOTOPE/MOD_MT_Frac.F90','   END SUBROUTINE tracer_craig_gordon_evap_affine', '''      CALL trace_write('AFFINE',[real(itrc,r8),vapor_ratio,temp_k,relhum, &
         min(max(relhum,0._r8),max(min(DEF_ISO_CG_RELHUM_MAX,0.999999_r8),0._r8)), &
         alpha_eq,alpha_k,slope,offset])
   END SUBROUTINE tracer_craig_gordon_evap_affine''')
edit('main/MATERIAL/GENERIC_WATER/MOD_MT_EvapLimit.F90','   END FUNCTION tracer_evaporative_tracer_loss', '''      CALL trace_write('LIMIT',[pool_trc,pool_water,water_loss,temp_k, &
         source_ratio,slope,offset,r_max,residual_water,residual_ratio, &
         tracer_evaporative_tracer_loss,pool_trc-tracer_evaporative_tracer_loss])
   END FUNCTION tracer_evaporative_tracer_loss''')
edit('main/MATERIAL/MOD_MT_Vars.F90','\t      IF (present(evap_kind)) kind = evap_kind','''\t      IF (present(evap_kind)) kind = evap_kind
         CALL trace_write('BOOK',[real(itrc,r8),real(kind,r8),tracer_mass,water_mass,max(water_mass,0._r8)])''')
# Capture full history mapping and per-patch gross accumulation for the three target cells.
p='main/MATERIAL/MOD_MT_Hist.F90'
edit(p,'   USE MOD_MT_Vars','   USE MOD_MT_Defs, only: trace_context,trace_write,trace_selected,trace_hist_date\n   USE MOD_MT_Vars')
start=(root/p).read_text().index('   SUBROUTINE write_history_tracer_delta_2d')
s=(root/p).read_text();a=s[:start];b=s[start:]
b=b.replace('   USE MOD_Block','   USE MOD_Vars_TimeInvariants, only: patchlatr,patchlonr,patchclass\n   USE MOD_Block',1)
b=b.replace('      integer :: compress','      integer :: ti,tproc,tloc,tx,ty\n      real(r8)::tlat,tlon\n      integer :: compress',1)
old='      CALL mp2g_hist%pset2grid (tracer_map_vec, tracer_xy_2d, spv = spval, msk = filter)'
new='''      IF(p_is_worker.and.(index(varname,'soilevap')>0.or.index(varname,'wetland_evap')>0)) THEN
         DO ip=1,size(water_acc)
            DO ti=1,mp2g_hist%npart(ip)
               tproc=mp2g_hist%address(ip)%val(1,ti);tloc=mp2g_hist%address(ip)%val(2,ti)
               ty=mp2g_hist%glist(tproc)%ilat(tloc);tx=mp2g_hist%glist(tproc)%ilon(tloc)
               tlat=(ghist%lat_s(ty)+ghist%lat_n(ty))*0.5_r8
               tlon=(ghist%lon_w(tx)+ghist%lon_e(tx))*0.5_r8
               IF(.not.trace_selected(tlat,tlon)) CYCLE
               CALL trace_context(ip,trace_hist_date,tlat,tlon,patchclass(ip))
               CALL trace_write('PATCH_'//trim(varname),[real(tx,r8),real(ty,r8), &
                  mp2g_hist%areapart(ip)%val(ti),water_acc(ip),tracer_mass_acc(ip), &
                  water_map_vec(ip),tracer_map_vec(ip),water_min,ref_ratio,merge(1._r8,0._r8,filter(ip)), &
                  patchlatr(ip)*180._r8/acos(-1._r8),patchlonr(ip)*180._r8/acos(-1._r8)])
            ENDDO
         ENDDO
      ENDIF
'''+old
assert old in b;b=b.replace(old,new,1)
# Per-grid actual mapped numerator/denominator and final output, so mask is explicit.
old='''                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel'''
new='''                  ENDIF
                  IF(index(varname,'soilevap')>0.or.index(varname,'wetland_evap')>0) THEN
                     tx=ghist%xdsp(xblk)+xloc;ty=ghist%ydsp(yblk)+yloc
                     tlat=(ghist%lat_s(ty)+ghist%lat_n(ty))*0.5_r8
                     tlon=(ghist%lon_w(tx)+ghist%lon_e(tx))*0.5_r8
                     CALL trace_context(0,trace_hist_date,tlat,tlon,0)
                     CALL trace_write('GRID_'//trim(varname),[real(tx,r8),real(ty,r8), &
                        tracer_xy_2d%blk(xblk,yblk)%val(xloc,yloc),water_xy_2d%blk(xblk,yblk)%val(xloc,yloc), &
                        ref_ratio,delta_xy_2d%blk(xblk,yblk)%val(xloc,yloc)])
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      compress = DEF_HIST_CompressLevel'''
assert old in b;b=b.replace(old,new,1);(root/p).write_text(a+b)
with (root/'Makefile').open('a') as f:f.write('\nMOD_MT_EvapLimit.o: MOD_MT_Defs.o\n')

# Use the actual history timestamp, not the record number within each month.
edit('main/MOD_Hist.F90','   USE MOD_MT_Vars, only: flush_Tracer_Acc',
     '   USE MOD_MT_Defs, only: trace_hist_date\n   USE MOD_MT_Vars, only: flush_Tracer_Acc')
edit('main/MOD_Hist.F90','            history_window_seconds=0._r8',
     '            trace_hist_date=idate\n            history_window_seconds=0._r8')
# Raw surface-relative humidity distinguishes meteorology from both humidity caps.
edit('main/MATERIAL/FAMILY_ISOTOPE/MOD_MT_Frac.F90',
     '         tracer_surface_relhum = min(max(eair / esat, 0._r8), 0.999999_r8)',
     """         tracer_surface_relhum = min(max(eair / esat, 0._r8), 0.999999_r8)
         CALL trace_write('RH',[qair,psrf,temp_k,merge(1._r8,0._r8,over_ice), &
            eair,esat,eair/esat,tracer_surface_relhum])""")
