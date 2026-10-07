from math import isinf
import numpy as np

def fit_sse(X,y):
    b=np.linalg.lstsq(X,y,rcond=None)[0];r=y-X@b;return float(r@r),b

def sse_tol(y):
    y=np.asarray(y,float); eps=np.finfo(float).eps; scale=np.sum((y-y.mean())**2); amp=max(1.0,np.max(np.abs(y)))
    return max(100*eps*scale, len(y)*(eps*amp)**2*1000, np.finfo(float).tiny)

def supf(y,years,L=5):
    years=np.asarray(years,float);y=np.asarray(y,float);t=years-years[0];s0,_=fit_sse(np.c_[np.ones(len(y)),t],y)
    ss=[];co=[]
    for idx in range(L,len(y)-L+1):
        tau=years[idx-1];X=np.c_[np.ones(len(y)),t,np.maximum(0,years-tau)];s,b=fit_sse(X,y);ss.append(s);co.append((tau,b))
    best=min(ss);tol=max(sse_tol(y),100*np.finfo(float).eps*abs(best));ties=[i for i,s in enumerate(ss) if abs(s-best)<=tol];pick=ties[0];bt,bb=co[pick]
    numer=s0-best;eps=sse_tol(y)
    if numer<0 and abs(numer)<=eps:numer=0
    if best<=eps:F=0.0 if numer<=eps else np.inf
    else:F=(len(y)-3)*numer/best
    return F,bt,bb

def discrete_lookup(k,B,method='BH'):
    k=np.asarray(k,int);counts=np.bincount(k,minlength=B+1);m=len(k);p=(1+np.arange(B+1))/(B+1);present=np.where(counts>0)[0];rank=np.cumsum(counts)[present]
    factor=sum(1/i for i in range(1,m+1)) if method=='BY' else 1;raw=factor*m*p[present]/rank;q=np.minimum.accumulate(raw[::-1])[::-1];q=np.minimum(q,1);return dict(zip(present,q))

years=np.arange(2001,2021);t=years-years[0]
F,tau,b=supf(5+2*t,years);assert abs(F)<1e-8
F,tau,b=supf(10-.5*t+1*np.maximum(0,years-2010),years);assert isinf(F) and tau==2010
rng=np.random.default_rng(4);y=2+.2*t+rng.normal(0,.03,len(t));F1,*_=supf(y,years);F2,*_=supf(y*1e-6,years);assert np.isfinite(F1) and np.isfinite(F2) and abs(F1-F2)<1e-4
null=np.array([1.,1.,2.,3.]);assert np.sum(null>=1)==4 and np.sum(null>=2)==2
k=np.array([0,0,1,1,1,3,5,9]);B=9
for method in ('BH','BY'):
    lookup=discrete_lookup(k,B,method);q=np.array([lookup[x] for x in k]);p=(1+k)/(B+1);order=np.argsort(p,kind='stable');ps=p[order];m=len(p);factor=sum(1/i for i in range(1,m+1)) if method=='BY' else 1
    raw=factor*m*ps/np.arange(1,m+1);adj=np.minimum.accumulate(raw[::-1])[::-1];adj=np.minimum(adj,1);ref=np.empty(m);ref[order]=adj;assert np.max(np.abs(q-ref))<1e-12,(method,q,ref)
print('Python reference checks passed: exact boundaries, scale invariance, >= MC ties, grouped BH/BY.')
