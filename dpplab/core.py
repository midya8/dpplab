"""Exact finite real-symmetric L-ensembles, not a sampling or ranking API."""
from fractions import Fraction as F
from itertools import combinations
import re

def subsets(n):
    return [list(s) for k in range(n+1) for s in combinations(range(n), k)]

def minor(a, s):
    return [[a[i][j] for j in s] for i in s]

def determinant(a):
    a=[row[:] for row in a]; out=F(1); n=len(a)
    for i in range(n):
        p=next((p for p in range(i,n) if a[p][i]), None)
        if p is None:return F(0)
        if p!=i:a[p],a[i]=a[i],a[p];out=-out
        pivot=a[i][i];out*=pivot
        for j in range(i+1,n):
            q=a[j][i]/pivot
            for k in range(i+1,n):a[j][k]-=q*a[i][k]
    return out

def inverse(a):
    n=len(a); b=[row[:]+[F(i==j) for j in range(n)] for i,row in enumerate(a)]
    for i in range(n):
        p=next((p for p in range(i,n) if b[p][i]),None)
        if p is None:raise ValueError('Singular matrix')
        b[p],b[i]=b[i],b[p];d=b[i][i];b[i]=[v/d for v in b[i]]
        for j in range(n):
            if j!=i:
                d=b[j][i];b[j]=[x-d*y for x,y in zip(b[j],b[i])]
    return [row[n:] for row in b]

def correlation(a):
    n=len(a);b=inverse([[a[i][j]+F(i==j) for j in range(n)] for i in range(n)])
    return [[F(i==j)-b[i][j] for j in range(n)] for i in range(n)]

def parse(spec):
    if not isinstance(spec,dict) or set(spec)!={'L','include','exclude','k'}:raise ValueError('Expected exactly L, include, exclude, k')
    a=spec['L'];n=len(a) if isinstance(a,list) else 0
    if not 1<=n<=5 or any(not isinstance(r,list) or len(r)!=n for r in a):raise ValueError('L must be square, size 1..5')
    def fraction(s):
        if not isinstance(s,str) or len(s)>12 or not re.fullmatch(r'-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?',s):raise ValueError('Matrix entries must be exact integer/fraction strings')
        v=F(s)
        if abs(v.numerator)>10000 or v.denominator>10000:raise ValueError('Reduced numerator/denominator bound 10000')
        return v
    a=[[fraction(v) for v in r] for r in a]
    if any(a[i][j]!=a[j][i] for i in range(n) for j in range(n)):raise ValueError('L must be symmetric')
    for key in ['include','exclude']:
        s=spec[key]
        if not isinstance(s,list) or any(type(i) not in (int,float) or not 0<=i<n or int(i)!=i for i in s) or len(s)!=len(set(s)):raise ValueError('Invalid '+key+' indices')
    if set(spec['include'])&set(spec['exclude']):raise ValueError('Include and exclude must be disjoint')
    k=spec['k']
    if k is not None and (type(k) not in (int,float) or not 0<=k<=n or int(k)!=k):raise ValueError('k must be null or integer 0..n')
    ss=subsets(n);weights=[determinant(minor(a,s)) for s in ss]
    if any(v<0 for v in weights):raise ValueError('L is not positive semidefinite: negative principal minor')
    return a,ss,weights

def audit(spec):
    a,ss,weights=parse(spec);n=len(a);inc=sorted(int(i) for i in spec['include']);exc=sorted(int(i) for i in spec['exclude']);k=None if spec['k'] is None else int(spec['k']);universe=list(range(n))
    z=determinant([[a[i][j]+F(i==j) for j in universe] for i in universe]);K=correlation(a)
    event=[set(inc)<=set(s) and not set(exc)&set(s) for s in ss]
    match=[e and (k is None or len(s)==k) for s,e in zip(ss,event)]
    event_z=sum((w for w,e in zip(weights,event) if e),F(0));chosen_z=sum((w for w,e in zip(weights,match) if e),F(0))
    remaining=[i for i in universe if i not in inc+exc];schur=None;ck=None
    if event_z:
        inv=inverse(minor(a,inc))
        schur=[[a[i][j]-sum((a[i][p]*inv[u][v]*a[q][j] for u,p in enumerate(inc) for v,q in enumerate(inc)),F(0)) for j in remaining] for i in remaining]
        ck=correlation(schur)
    rows=[]
    for s,w,e,m in zip(ss,weights,event,match):
        incl=sum((v for t,v in zip(ss,weights) if set(s)<=set(t)),F(0))/z
        rows.append({'subset':s,'determinant':w,'probability':w/z,'inclusion_sum':incl,'inclusion_determinant':determinant(minor(K,s)),'event_match':e,'selected_match':m,'conditional_probability':w/chosen_z if m and chosen_z else F(0) if chosen_z else None})
    probs=[r['conditional_probability'] for r in rows]
    post=[sum((p for s,p in zip(ss,probs) if i in s),F(0)) for i in universe] if chosen_z else None
    cardinality=[sum((p for s,p in zip(ss,probs) if len(s)==j),F(0)) for j in range(n+1)] if chosen_z else None
    mean=sum((F(j)*p for j,p in enumerate(cardinality)),F(0)) if chosen_z else None
    variance=sum((F(j*j)*p for j,p in enumerate(cardinality)),F(0))-mean*mean if chosen_z else None
    baseline_mean=sum((K[i][i] for i in universe),F(0));baseline_var=baseline_mean-sum((K[i][j]*K[j][i] for i in universe for j in universe),F(0))
    pairs=[{'i':i,'j':j,'joint':K[i][i]*K[j][j]-K[i][j]*K[j][i],'independent_product':K[i][i]*K[j][j],'covariance':-K[i][j]*K[j][i]} for i in universe for j in range(i+1,n)]
    best=max(weights);modes=[s for s,w in zip(ss,weights) if w==best]
    selected_best=max((w for w,m in zip(weights,match) if m),default=F(0));selected_modes=[s for s,w,m in zip(ss,weights,match) if m and w==selected_best] if chosen_z else []
    result={'spec':spec,'normalizer':z,'enumerated_normalizer':sum(weights,F(0)),'K':K,'baseline_mean':baseline_mean,'baseline_variance':baseline_var,'pairs':pairs,'rows':rows,'modes':modes,'event_weight':event_z,'event_probability':event_z/z,'selected_weight':chosen_z,'selected_probability':chosen_z/z,'remaining':remaining,'event_conditional_L':schur,'event_conditional_K':ck,'selected_inclusions':post,'selected_cardinality':cardinality,'selected_mean':mean,'selected_variance':variance,'selected_modes':selected_modes}
    def pack(x):
        if isinstance(x,F):return str(x)
        if isinstance(x,list):return [pack(v) for v in x]
        if isinstance(x,dict):return {k:pack(v) for k,v in x.items()}
        return x
    return pack(result)
