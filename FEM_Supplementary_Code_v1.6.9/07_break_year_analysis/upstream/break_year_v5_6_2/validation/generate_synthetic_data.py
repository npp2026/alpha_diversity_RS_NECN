#!/usr/bin/env python3
"""CI/build fallback synthetic-data generator.

The canonical publication pipeline generator is R/synthetic_data.R. This Python
implementation consumes the same validation/synthetic_spec.yml so that the
fixture can also be materialized in environments where R/terra is unavailable.
The stochastic stream is intentionally not claimed to be bit-identical to R;
truth classes, spatial layout, signal equations, AR(1) marginal SDs, and missing
patterns are identical by specification.
"""
from __future__ import annotations
import argparse, hashlib, json, math, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
import rasterio
from rasterio.transform import from_origin
import geopandas as gpd
from shapely.geometry import box


def load_spec(path: Path):
    with path.open('r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def class_table(spec):
    rows=[]
    for z in spec['classes']:
        rows.append({
            'code': int(z['code']), 'name': z['name'],
            'slope_before': float(z['slope_before']), 'slope_after': float(z['slope_after']),
            'break_year': np.nan if z.get('break_year') is None else int(z['break_year']),
            'noise': bool(z['noise']), 'missing_year': np.nan if z.get('missing_year') is None else int(z['missing_year'])
        })
    return pd.DataFrame(rows)


def layout_matrix(spec):
    nr,nc=int(spec['raster']['nrow']),int(spec['raster']['ncol'])
    out=np.zeros((nr,nc),dtype=np.uint8); occ=np.zeros((nr,nc),dtype=bool)
    codes={z['name']:int(z['code']) for z in spec['classes']}
    for z in spec['layout']:
        r0,r1=z['rows']; c0,c1=z['cols']; sl=(slice(r0-1,r1),slice(c0-1,c1))
        if occ[sl].any(): raise ValueError('synthetic layout overlap')
        out[sl]=codes[z['class']]; occ[sl]=True
    return out


def region_matrix(spec):
    nr,nc=int(spec['raster']['nrow']),int(spec['raster']['ncol'])
    out=np.zeros((nr,nc),dtype=np.uint8)
    for z in spec['regions']:
        r0,r1=z['rows']; c0,c1=z['cols']; sl=(slice(r0-1,r1),slice(c0-1,c1))
        if (out[sl]!=0).any(): raise ValueError('region overlap')
        out[sl]=int(z['id'])
    if (out==0).any(): raise ValueError('regions do not cover raster')
    return out


def stable_seed(master, *parts):
    h=hashlib.sha256('|'.join(map(str,(master,)+parts)).encode()).digest()
    return int.from_bytes(h[:8],'little') % (2**32-1)


def ar1(n,rho,sd,rng):
    e=np.zeros(n,float); e[0]=rng.normal(0,sd)
    innov=sd*math.sqrt(max(1-rho*rho,0))
    for i in range(1,n): e[i]=rho*e[i-1]+rng.normal(0,innov)
    return e


def write_tif(path, arr, transform, crs, nodata=-9999.0, dtype='float32'):
    path.parent.mkdir(parents=True,exist_ok=True)
    data=np.asarray(arr)
    if np.issubdtype(np.dtype(dtype),np.floating): data=np.where(np.isfinite(data),data,nodata).astype(dtype)
    else: data=data.astype(dtype)
    with rasterio.open(path,'w',driver='GTiff',height=data.shape[0],width=data.shape[1],count=1,dtype=dtype,
                       crs=crs,transform=transform,nodata=nodata if np.issubdtype(np.dtype(dtype),np.floating) else None,
                       compress='deflate') as dst: dst.write(data,1)


def generate(project: Path, spec_path: Path, out: Path, overwrite=False):
    spec=load_spec(spec_path); ct=class_table(spec); cm=layout_matrix(spec); rm=region_matrix(spec)
    if out.exists() and overwrite: shutil.rmtree(out)
    out.mkdir(parents=True,exist_ok=True)
    nr,nc=cm.shape; res=float(spec['raster']['resolution_m']); xmin=float(spec['raster']['xmin']); ymin=float(spec['raster']['ymin'])
    ymax=ymin+nr*res; transform=from_origin(xmin,ymax,res,res); crs=spec['raster']['crs']; years=np.array(spec['years'],dtype=int)
    write_tif(out/'synthetic_truth_class.tif',cm,transform,crs,dtype='uint8')
    write_tif(out/'synthetic_region_id.tif',rm,transform,crs,dtype='uint8')
    # Polygon quadrants exactly follow pixel boundaries from the spec.
    geoms=[]; names=[]; ids=[]
    for z in spec['regions']:
        r0,r1=z['rows']; c0,c1=z['cols']
        x0=xmin+(c0-1)*res; x1=xmin+c1*res; y1=ymax-(r0-1)*res; y0=ymax-r1*res
        geoms.append(box(x0,y0,x1,y1)); names.append(z['name']); ids.append(int(z['id']))
    rg=gpd.GeoDataFrame({'region_id':ids,'NAME':names},geometry=geoms,crs=crs)
    rdir=out/'regions'; rdir.mkdir(exist_ok=True); rg.to_file(rdir/'synthetic_regions.shp')
    truth=[]
    by_code=ct.set_index('code').to_dict('index'); region_names={int(z['id']):z['name'] for z in spec['regions']}
    for r in range(nr):
        for c in range(nc):
            cell=r*nc+c+1; code=int(cm[r,c]); z=by_code[code]
            x=xmin+(c+.5)*res; y=ymax-(r+.5)*res
            truth.append({'cell':cell,'row':r+1,'col':c+1,'x':x,'y':y,'region_id':int(rm[r,c]),'class_code':code,
                          'name':z['name'],'slope_before':z['slope_before'],'slope_after':z['slope_after'],
                          'break_year':z['break_year'],'noise':z['noise'],'missing_year':z['missing_year'],'region':region_names[int(rm[r,c])]})
    truth=pd.DataFrame(truth); truth.to_csv(out/'synthetic_truth.csv',index=False); ct.to_csv(out/'synthetic_classes.csv',index=False)
    master=int(spec['seed']); center=years.mean()
    files=[]
    for response,rcfg in spec['responses'].items():
        vals=np.full((nr*nc,len(years)),np.nan,float)
        for i,row in truth.iterrows():
            z=by_code[int(row.class_code)]; scale=float(rcfg.get('slope_scale',1)); sb=z['slope_before']*scale; sa=z['slope_after']*scale
            rng=np.random.default_rng(stable_seed(master,'synthetic_intercept',response,int(row['row']),int(row['col'])))
            intercept=float(rcfg['baseline'])+rng.normal(0,float(rcfg['intercept_sd']))
            y=intercept+sb*(years-center)
            tau=z['break_year']
            if np.isfinite(tau) and abs(sa-sb)>0: y=y+(sa-sb)*np.maximum(0,years-tau)
            if z['noise']:
                rng=np.random.default_rng(stable_seed(master,'synthetic_noise',response,int(row['row']),int(row['col'])))
                y=y+ar1(len(years),float(rcfg['rho']),float(rcfg['noise_sd']),rng)
            if np.isfinite(z['missing_year']): y[years==int(z['missing_year'])]=np.nan
            vals[i]=y
        rdir=out/response; rdir.mkdir(exist_ok=True)
        for j,yr in enumerate(years):
            path=rdir/f'synthetic_{response}_{yr}_1km.tif'; write_tif(path,vals[:,j].reshape(nr,nc),transform,crs,dtype='float64')
            files.append(path)
    shutil.copy2(spec_path,out/'synthetic_spec_resolved.yml')
    manifest=[]
    for p in sorted([q for q in out.rglob('*') if q.is_file()]):
        h=hashlib.md5(p.read_bytes()).hexdigest(); manifest.append({'path':str(p.relative_to(project)),'size':p.stat().st_size,'md5':h})
    pd.DataFrame(manifest).to_csv(out/'synthetic_manifest.csv',index=False)
    meta={'version':spec['version'],'seed':master,'rng_kind':spec.get('rng_kind','L\'Ecuyer-CMRG'),'generator_engine':'python-ci-fallback','years':years.tolist(),'nrow':nr,'ncol':nc,
          'responses':{k:{'rho':v['rho'],'noise_sd':v['noise_sd']} for k,v in spec['responses'].items()}}
    (out/'synthetic_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(f'generated {len(files)} annual rasters + truth/regions in {out}')

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--project',default=None); ap.add_argument('--force',action='store_true')
    ns=ap.parse_args(); project=Path(ns.project).resolve() if ns.project else Path(__file__).resolve().parents[1]
    generate(project,project/'validation'/'synthetic_spec.yml',project/'validation'/'synthetic_data',ns.force)
