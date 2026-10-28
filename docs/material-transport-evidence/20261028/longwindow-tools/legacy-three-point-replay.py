"""Read only trace text. Decimal50 checks actual affine finite-pool inputs and gross-bookkeeping."""
from decimal import Decimal as D, localcontext
from pathlib import Path
from collections import defaultdict
import json,sys
root=Path(sys.argv[1]); out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
def decs(a):return [D(x) for x in a]
def target(lat,lon):
 lon=(lon+D(180))%D(360)-D(180) if lon>=D(180) else lon
 targets=[(D('-1.125'),D('-59.125'),2,'soil'),(D('-2.375'),D('-58.125'),8,'wet8'),(D('-2.375'),D('-58.375'),6,'wet6')]
 for la,lo,day,name in targets:
  if abs(lat-la)<D('.12501') and abs(lon-lo)<D('.12501'):return name,day
 return None
books=defaultdict(lambda:[D(0),D(0),0,0]); patches=[];grids=[];limits=[];errors=[];last_affine={};affines=[];pending=defaultdict(list)
with localcontext() as ctx:
 ctx.prec=50
 for file in sorted(root.glob('trace_rank_*.txt')):
  for lineno,line in enumerate(file.read_text().splitlines(),1):
   p=line.split();tag=p[0];rank,ip,yr,day,sec,pc=map(int,p[1:7]);lat,lon=decs(p[7:9]);v=decs(p[9:]);key=(rank,ip,day,sec)
   if tag=='AFFINE':
    it,rv,temp,rh,h,ae,ak,slope,offset=v
    reconstructed_slope=1/(ae*ak*(1-h));reconstructed_offset=-h*rv/(ak*(1-h))
    affines.append((abs(reconstructed_slope-slope),abs(reconstructed_offset-offset)))
    last_affine[key]=(v,lineno)
   elif tag=='LIMIT':
    m,w,loss,temp,r0,a,b,cap,w1,r1,flux,m1=v
    exactw1=w-loss
    x=(exactw1/w).ln();y=(a-1)*x;growth=y.exp()
    factor=x if y==0 else x*(growth-1)/y
    exactr=max(D(0),(m/w)*growth+b*factor)
    if cap>0:exactr=min(exactr,cap)
    exactflux=min(max(m,0)-exactr*exactw1,max(m,0))
    af=last_affine.get(key)
    if not af or af[0][-2:]!=[a,b]:errors.append(f'LIMIT without exact preceding AFFINE {file.name}:{lineno}')
    limits.append(dict(rank=rank,patch=ip,day=day,sec=sec,lat=str(lat),lon=str(lon),isotope=int(af[0][0]) if af else None,
       pool_water=str(w),water_loss=str(loss),pool_mass=str(m),ratio_pool=str(r0),ratio_vapor=str(af[0][1]) if af else None,
       rh_pre_cg=str(af[0][3]) if af else None,rh_post_cg=str(af[0][4]) if af else None,
       alpha_eq=str(af[0][5]) if af else None,alpha_k=str(af[0][6]) if af else None,
       model_loss=str(flux),decimal_loss=str(exactflux),loss_abs_error=str(abs(flux-exactflux)),
       model_after=str(m1),decimal_after=str(m-exactflux),model_mass_closure=str(m-flux-m1),
       line=lineno,file=file.name))
    if af:pending[(key,int(af[0][0]))].append(len(limits)-1)
   elif tag=='BOOK':
    it,kind,m,water,gross=v
    for ix in pending.pop((key,int(it)),[]):
     limits[ix]['next_book_kind']=int(kind)
     limits[ix]['next_book_mass']=str(m)
    actualday=day if sec>0 else day-1
    k=(rank,ip,actualday,int(it),int(kind));r=books[k];r[0]+=gross;r[1]+=m;r[2]+=1;r[3]+=m<0
   elif tag.startswith('PATCH_'):
    var=tag[6:];tx,ty,area,water,m,wm,mm,wmin,ref,mask,plat,plon=v
    name,tday=target(lat,lon)
    if day!=tday or ('soilevap' in var)!=(name=='soil'):continue
    it=1 if var.endswith('H2_18O') else 2;kind=2 if name=='soil' else 5
    bk=books[(rank,ip,day,it,kind)]
    eligible=water>wmin and mask==1
    ptarget=target(plat,plon)
    patches.append(dict(target=name,day=day,var=var,rank=rank,patch=ip,patchclass=pc,lat=str(plat),lon=str(plon),
      traced_centroid=ptarget is not None,area=str(area),water=str(water),signed_mass=str(m),mapped_water=str(wm),mapped_mass=str(mm),
      water_min=str(wmin),filter=bool(mask),eligible=eligible,ref=str(ref),book_water=str(bk[0]),book_mass=str(bk[1]),
      book_count=bk[2],negative_book_count=bk[3],water_residual=str(water-bk[0]),mass_residual=str(m-bk[1])))
   elif tag.startswith('GRID_'):
    var=tag[5:];tx,ty,m,water,ref,delta=v
    name,tday=target(lat,lon)
    if day!=tday or ('soilevap' in var)!=(name=='soil'):continue
    ps=[z for z in patches if z['target']==name and z['var']==var] # ranks can appear later; aggregate below
    grids.append(dict(target=name,day=day,var=var,lat=str(lat),lon=str(lon),signed_mass=str(m),gross_water=str(water),ref=str(ref),model_delta=str(delta)))
 for g in grids:
  ps=[p for p in patches if p['target']==g['target'] and p['var']==g['var']]
  assert ps,(g,'missing contribution patches')
  mass=sum((D(p['area'])*D(p['mapped_mass']) for p in ps if p['filter']),D(0))
  water=sum((D(p['area'])*D(p['mapped_water']) for p in ps if p['filter']),D(0))
  ratio=mass/water if water>0 else None
  delta=(ratio/D(g['ref'])-1)*1000 if ratio is not None else None
  reason='gross_water_nonpositive' if water<=0 else 'negative_signed_ratio' if ratio<0 else 'delta_sanity' if abs(delta)>2000 else 'valid'
  g.update(patches=len(ps),eligible_patches=sum(p['eligible'] for p in ps),decimal_delta=str(delta),mask_reason=reason,
    patch_reconstruction_delta_error=str(abs(delta-D(g['model_delta']))) if reason=='valid' else None)
 limerrs=[D(x['loss_abs_error']) for x in limits]
 summary=dict(trace_files=len(list(root.glob('trace_rank_*.txt'))),finite_pool_events=len(limits),negative_finite_pool_events=sum(D(x['model_loss'])<0 for x in limits),
  max_finite_pool_mass_abs_error=str(max(limerrs,default=D(0))),max_affine_slope_abs_error=str(max((x[0] for x in affines),default=D(0))),
  max_affine_offset_abs_error=str(max((x[1] for x in affines),default=D(0))),
  selected_patch_records=len(patches),all_contributing_centroids_traced=all(p['traced_centroid'] for p in patches),
  max_daily_water_book_residual=str(max((abs(D(p['water_residual'])) for p in patches),default=D(0))),
  max_daily_signed_mass_book_residual=str(max((abs(D(p['mass_residual'])) for p in patches),default=D(0))),grids=grids,errors=errors)
 for name,data in [('summary',summary),('patches',patches),('finite_pool_events',limits)]:
  (out/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n')
 daily=[]
 for p in patches:
  if not p['eligible']:continue
  isotope=1 if p['var'].endswith('H2_18O') else 2
  kind=2 if p['target']=='soil' else 5
  es=[e for e in limits if (e['rank'],e['patch'],e['day'],e['isotope'],e.get('next_book_kind')) == (p['rank'],p['patch'],p['day'],isotope,kind)]
  totalw=sum((D(e['water_loss']) for e in es),D(0))
  totalm=sum((D(e['model_loss']) for e in es),D(0))
  replaym=sum((D(e['decimal_loss']) for e in es),D(0))
  daily.append(dict(target=p['target'],isotope=isotope,rank=p['rank'],patch=p['patch'],day=p['day'],events=len(es),
    gross_water=str(totalw),signed_mass=str(totalm),decimal_signed_mass=str(replaym),
    decimal_delta=str((replaym/totalw/D(p['ref'])-1)*1000),
    max_step_mass_abs_error=str(max((D(e['loss_abs_error']) for e in es),default=D(0))),
    negative_events=sum(D(e['model_loss'])<0 for e in es)))
  assert totalw==D(p['book_water']) and totalm==D(p['book_mass']),p
 (out/'finite_daily_replay.json').write_text(json.dumps(daily,indent=2)+'\n')
 assert not errors and summary['all_contributing_centroids_traced']
 assert len(grids)==6 and len(daily)==6
 assert max(limerrs)<D('1e-12')
 assert D(summary['max_daily_water_book_residual'])<D('1e-14')
 assert D(summary['max_daily_signed_mass_book_residual'])<D('1e-18')
 assert all(D(g['patch_reconstruction_delta_error'])<D('1e-9') for g in grids if g['mask_reason']=='valid')
 assert sum(g['mask_reason']=='negative_signed_ratio' for g in grids)==1
 print(json.dumps(summary,indent=2))
