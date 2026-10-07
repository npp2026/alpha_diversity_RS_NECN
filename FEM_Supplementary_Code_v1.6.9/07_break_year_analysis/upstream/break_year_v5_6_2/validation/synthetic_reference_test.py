#!/usr/bin/env python3
"""Independent synthetic-data reference test for R-less build environments.

This does not validate terra/R integration. It validates the materialized fixture
and independently re-implements the core hinge/SupF/MC/BH logic to catch design
or science regressions before the R E2E test is run.
"""
from __future__ import annotations
import json, math, re
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio, yaml

ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'validation'/'synthetic_data'

def resid_maker(X):
    return np.eye(X.shape[0])-X@np.linalg.inv(X.T@X)@X.T

def candidate_indices(n,m): return np.arange(m-1,n-m)  # R indices m..n-m -> Python m-1..n-m-1

def sse_tol(y):
    eps=np.finfo(float).eps; scale=((y-y.mean())**2).sum(); amp=max(1.0,np.abs(y).max())
    return max(100*eps*scale,len(y)*(eps*amp)**2*1000,np.finfo(float).tiny)

def fit_one(y,years,m):
    finite=np.isfinite(y)
    if finite.sum()!=len(years): return {'testable':False,'status':0,'F':np.nan,'tau':np.nan,'direction':False}
    t=years-years[0]; X0=np.c_[np.ones(len(years)),t]; b0=np.linalg.lstsq(X0,y,rcond=None)[0]; sse0=((y-X0@b0)**2).sum()
    best=np.inf; fits=[]
    for idx in candidate_indices(len(years),m):
        tau=years[idx]; X=np.c_[np.ones(len(years)),t,np.maximum(0,years-tau)]; b=np.linalg.lstsq(X,y,rcond=None)[0]; sse=((y-X@b)**2).sum()
        fits.append((sse,tau,b))
        best=min(best,sse)
    tol=max(sse_tol(y),100*np.finfo(float).eps*abs(best)); choices=[z for z in fits if abs(z[0]-best)<=tol]; sse,tau,b=choices[0]
    numer=sse0-sse; e=sse_tol(y)
    if numer<0 and abs(numer)<=e: numer=0
    if sse<=e:
        F=0.0 if numer<=e else np.inf; status=2 if numer<=e else 3
    else:
        F=(len(years)-3)*numer/sse; status=1
    return {'testable':True,'status':status,'F':F,'tau':tau,'direction':bool(b[1]<0 and b[1]+b[2]>0)}

def ar1_matrix(B,n,rho,rng):
    e=np.zeros((B,n)); e[:,0]=rng.normal(size=B); sd=math.sqrt(1-rho*rho)
    for j in range(1,n): e[:,j]=rho*e[:,j-1]+rng.normal(scale=sd,size=B)
    return e

def supf_matrix(Y,years,m):
    t=years-years[0]; M0=resid_maker(np.c_[np.ones(len(years)),t]); s0=np.einsum('ij,ij->i',Y@M0,Y)
    best=np.full(Y.shape[0],np.inf)
    for idx in candidate_indices(len(years),m):
        tau=years[idx]; M=resid_maker(np.c_[np.ones(len(years)),t,np.maximum(0,years-tau)]); s=np.einsum('ij,ij->i',Y@M,Y); best=np.minimum(best,s)
    num=s0-best; out=(len(years)-3)*num/best; return np.maximum(out,0)

def bh(p):
    p=np.asarray(p,float); q=np.full_like(p,np.nan); ok=np.isfinite(p); x=p[ok]; m=len(x); order=np.argsort(x); sx=x[order]
    adj=np.minimum.accumulate((sx*m/np.arange(1,m+1))[::-1])[::-1]; adj=np.minimum(adj,1); tmp=np.empty(m); tmp[order]=adj; q[ok]=tmp; return q

def load_stack(response,years):
    arr=[]
    for y in years:
        with rasterio.open(DATA/response/f'synthetic_{response}_{y}_1km.tif') as ds:
            a=ds.read(1).astype(float); a[a==ds.nodata]=np.nan; arr.append(a.ravel())
    return np.stack(arr,axis=1)

def hinge_rho_values(Y,years,m=3):
    good=np.isfinite(Y).all(axis=1); Y=Y[good]; n=len(years); t=years-years[0]
    best=np.full(len(Y),np.inf); Rbest=np.zeros_like(Y)
    for idx in candidate_indices(n,m):
        tau=years[idx]; M=resid_maker(np.c_[np.ones(n),t,np.maximum(0,years-tau)]); R=Y@M; ss=np.einsum('ij,ij->i',R,R); take=ss<best; best[take]=ss[take]; Rbest[take]=R[take]
    a=Rbest[:,:-1]; b=Rbest[:,1:]; den=np.sqrt(np.einsum('ij,ij->i',a,a)*np.einsum('ij,ij->i',b,b)); return np.einsum('ij,ij->i',a,b)/den

def hinge_bias_correct(raw,years,grid,B=1500):
    est=[]
    for i,rho in enumerate(grid):
        Y=ar1_matrix(B,len(years),float(rho),np.random.default_rng(9100+i)); est.append(float(np.median(hinge_rho_values(Y,years,3))))
    est=np.maximum.accumulate(np.asarray(est))
    # Collapse near-duplicate monotone plateaus before interpolation.
    keep=np.r_[True,np.diff(est)>1e-10]; return float(np.interp(raw,est[keep],np.asarray(grid)[keep],left=grid[0],right=grid[-1]))

def run():
    spec=yaml.safe_load((ROOT/'validation'/'synthetic_spec.yml').read_text()); years=np.array(spec['years'],dtype=int); truth=pd.read_csv(DATA/'synthetic_truth.csv').sort_values('cell')
    metrics=[]; hard=[]
    hard.append(("fixture:cell_count_576",len(truth)==576))
    hard.append(("fixture:four_equal_regions",truth.groupby("region").size().sort_values().tolist()==[144,144,144,144]))
    hard.append(("fixture:missing_class_count_12",int((truth["name"]=="missing_2010").sum())==12))
    scfg=yaml.safe_load((ROOT/'config'/'v5_6_synthetic.yml').read_text())
    for response in ('SR','Shannon'):
        pat=re.compile(scfg['data']['responses'][response]['pattern']); files=sorted((DATA/response).glob('*.tif')); yrs=[]
        for f in files:
            m=pat.match(f.name)
            if m: yrs.append(int(m.group(1)))
        hard.append((f'fixture:{response}_config_discovers_20_years',yrs==spec['years']))
    hard.append(('fixture:regions_shapefile_present',(DATA/'regions'/'synthetic_regions.shp').exists()))
    with rasterio.open(DATA/'synthetic_truth_class.tif') as ds:
        truth_r=ds.read(1).ravel()
    hard.append(('fixture:truth_raster_matches_csv',np.array_equal(truth_r.astype(int),truth['class_code'].to_numpy(dtype=int))))
    for response,rcfg in spec['responses'].items():
        Y=load_stack(response,years)
        # Mirror the v5.6 primary nuisance-rho strategy: best-hinge residual rho
        # plus short-series simulation inversion, then calibrate SupF at that rho.
        raw_rho=float(np.median(hinge_rho_values(Y,years,3))); grid=np.arange(-.8,.8001,.05); corrected=hinge_bias_correct(raw_rho,years,grid)
        true_rho=float(rcfg['rho']); rho_err=abs(corrected-true_rho)
        fits=[fit_one(Y[i],years,5) for i in range(len(Y))]
        F=np.array([z['F'] for z in fits]); tau=np.array([z['tau'] for z in fits]); direction=np.array([z['direction'] for z in fits]); status=np.array([z['status'] for z in fits])
        B=4999; null=np.sort(supf_matrix(ar1_matrix(B,len(years),corrected,np.random.default_rng(741+len(response))),years,5))
        p=np.full(len(Y),np.nan); ok=np.isfinite(F)|np.isposinf(F)
        for i in np.where(ok)[0]:
            k=B-np.searchsorted(null,F[i],side='left'); p[i]=(1+k)/(B+1)
        q=bh(p); recovery=(q<=.05)&direction
        def mask(name): return truth['name'].eq(name).to_numpy()
        exlin=mask('exact_linear'); exhinge=mask('exact_hinge_2010'); miss=mask('missing_2010'); rev=mask('reverse_2010'); nul=mask('null_linear')
        strong=mask('recovery_2008')|mask('recovery_2012')|mask('boundary_recovery_2005')|mask('boundary_recovery_2015')|exhinge
        exact_linear_pass=bool(np.all(status[exlin]==2)); exact_hinge_pass=bool(np.all(status[exhinge]==3)&np.all(tau[exhinge]==2010)); missing_pass=bool(np.all(~np.isfinite(p[miss])))
        rec_recall=float(recovery[strong].mean()); null_fpr=float(recovery[nul].mean()); reverse_rate=float(recovery[rev].mean())
        hard += [(f'{response}:exact_linear',exact_linear_pass),(f'{response}:exact_hinge',exact_hinge_pass),(f'{response}:missing_excluded',missing_pass),
                 (f'{response}:strong_recovery_recall>=0.60',rec_recall>=.60),(f'{response}:null_recovery_fpr<=0.05',null_fpr<=.05),(f'{response}:reverse_recovery<=0.05',reverse_rate<=.05)]
        metrics += [dict(response=response,metric='strong_recovery_recall',value=rec_recall),dict(response=response,metric='null_recovery_fpr',value=null_fpr),dict(response=response,metric='reverse_recovery_rate',value=reverse_rate)]
        # Synthetic science check: primary rho estimation now uses hinge residuals so
        # real structural breaks do not masquerade as AR(1) persistence.
        metrics += [dict(response=response,metric='hinge_rho_raw_median',value=raw_rho),dict(response=response,metric='hinge_rho_bias_corrected',value=corrected),dict(response=response,metric='hinge_rho_absolute_error',value=rho_err)]
        hard.append((f'{response}:hinge_rho_error<=0.15',rho_err<=.15))
        # Window sensitivity is tested on raw break estimates, not significance.
        tau_by_m={5:tau}; rec_by_m={5:recovery}
        for m in (5,4,3,2):
            if m==5: tt=tau
            else: tt=np.array([fit_one(Y[i],years,m)['tau'] for i in range(len(Y))])
            tau_by_m[m]=tt
            for cname,true_tau in [('early_recovery_2003',2003),('late_recovery_2017',2017)]:
                med=float(np.nanmedian(tt[mask(cname)])); metrics.append(dict(response=response,metric=f'{cname}_median_break_L{m}',value=med))
                if m==3: hard.append((f'{response}:{cname}:L3 within 1y',abs(med-true_tau)<=1))
        # Independent footprint-vs-time check for the new spatial interpretation fields.
        fits3=[fit_one(Y[i],years,3) for i in range(len(Y))]; F3=np.array([z['F'] for z in fits3]); dir3=np.array([z['direction'] for z in fits3])
        null3=np.sort(supf_matrix(ar1_matrix(B,len(years),corrected,np.random.default_rng(1741+len(response))),years,3)); p3=np.full(len(Y),np.nan); ok3=np.isfinite(F3)|np.isposinf(F3)
        for i in np.where(ok3)[0]: p3[i]=(1+(B-np.searchsorted(null3,F3[i],side='left')))/(B+1)
        rec3=(bh(p3)<=.05)&dir3; rec_by_m[3]=rec3; common=recovery&rec3; union=recovery|rec3
        iou=float(common.sum()/union.sum()) if union.sum() else np.nan; shifts=np.abs(tau_by_m[3][common]-tau_by_m[5][common])
        within1=float(np.mean(shifts<=1)) if len(shifts) else np.nan; mean_shift=float(np.mean(shifts)) if len(shifts) else np.nan
        metrics += [dict(response=response,metric='L5_vs_L3_recovery_mask_iou',value=iou),dict(response=response,metric='L5_vs_L3_within1_year_agreement',value=within1),dict(response=response,metric='L5_vs_L3_mean_abs_year_shift',value=mean_shift)]
        hard.append((f'{response}:spatial_temporal_metrics_finite',np.isfinite(iou) and np.isfinite(within1) and np.isfinite(mean_shift)))
        print(response,'recall',rec_recall,'null FPR',null_fpr,'reverse',reverse_rate)
    pd.DataFrame(metrics).to_csv(ROOT/'validation'/'synthetic_reference_metrics.csv',index=False)
    pd.DataFrame(hard,columns=['check','pass']).to_csv(ROOT/'validation'/'synthetic_reference_checks.csv',index=False)
    failed=[n for n,p in hard if not p]
    print('Synthetic reference checks:',len(hard)-len(failed),'/',len(hard),'passed')
    if failed:
        print('FAILED:',*failed,sep='\n - '); raise SystemExit(1)

if __name__=='__main__': run()
