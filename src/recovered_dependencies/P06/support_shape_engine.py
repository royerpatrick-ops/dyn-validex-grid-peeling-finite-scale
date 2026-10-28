"""Qualified analytic support geometry, then integer extreme-vertex peeling.

Initial smooth shape: h(phi)=1-3/40*cos(3*phi), rho=1+3/5*cos(3*phi).
The shape is not rescaled. theta=atan2(q,p) is mathematical, evaluated
with mpmath.atan2 at 80/120 digits near integer row boundaries. Rational
translations are passed as integers tx,ty,D. Initial membership is
qualified numerical geometry, NOT global interval certification. Once
the initial points are classified, every hull/removal predicate is integer.

Run this file to execute geometry qualification only, never scientific
transfer trajectories.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.spatial import ConvexHull

ROOT=Path(__file__).resolve().parent
M=3
EPSILON=Fraction(3,5)
AMPLITUDE=Fraction(-3,40)
STATS={"mp_row_rechecks":0,"exact_boundary_recognitions":0,
       "ambiguous_oracle_shared_rechecks":0}


def geometry():
    return {"name":"support_m3","support":"h(phi)=1-3/40*cos(3*phi)",
            "radius_of_curvature":"rho(phi)=1+3/5*cos(3*phi)",
            "P":2*math.pi,"area":391*math.pi/400,
            "area_over_pi":"391/400","rho_min":.4,"rho_max":1.6,
            "rmax_bound":43/40,"rescaled":False,
            "initial_membership":"qualified floating geometry with 80/120-digit rechecks, not interval certified"}


def _parameters(R,p,q,tx,ty,D):
    if any(not isinstance(v,int) for v in (R,p,q,tx,ty,D)):
        raise TypeError("R,p,q,tx,ty,D must be Python integers")
    if R<=0 or p<=0 or q<0 or D<=0:
        raise ValueError("Require R>0,p>0,q>=0,D>0")


def support(alpha,p=1,q=0):
    theta=math.atan2(q,p)
    return 1-.075*np.cos(3*(np.asarray(alpha)-theta))


def boundary(alpha,p=1,q=0):
    theta=math.atan2(q,p)
    angle=np.asarray(alpha)
    phase=3*(angle-theta)
    h=1-.075*np.cos(phase)
    hp=.225*np.sin(phase)
    c,s=np.cos(angle),np.sin(angle)
    return h*c-hp*s,h*s+hp*c


def _aligned_rational_boundary(x,y,alpha):
    """A positive exact witness, not a general membership test."""
    c=Fraction(math.cos(float(alpha))).limit_denominator(1000000)
    s=Fraction(math.sin(float(alpha))).limit_denominator(1000000)
    if c*c+s*s!=1:return False
    cos3=c*c*c-3*c*s*s
    sin3=3*c*c*s-s*s*s
    h=1+AMPLITUDE*cos3
    hp=-3*AMPLITUDE*sin3
    return h*c-hp*s==x and h*s+hp*c==y


def _mp_section(R,p,q,tx,ty,D,row):
    STATS['mp_row_rechecks']+=1
    answers=[]
    for dps in (80,120):
        with mp.workdps(dps):
            theta=mp.atan2(q,p)
            amplitude=-mp.mpf(3)/40
            yy=(mp.mpf(row)+mp.mpf(ty)/D)/R
            ux=mp.mpf(tx)/D
            def point(alpha):
                h=1+amplitude*mp.cos(3*(alpha-theta))
                hp=-3*amplitude*mp.sin(3*(alpha-theta))
                c,s=mp.cos(alpha),mp.sin(alpha)
                return h*c-hp*s,h*s+hp*c
            # Aligned vertical tangencies have exact rational coordinates:
            # (x,y)=(9/40,+/-1) on the unit shape.
            if q==0:
                exacty=Fraction(D*row+ty,D)
                if exacty>R or exacty<-R:
                    answers.append(None);continue
                if exacty==R or exacty==-R:
                    xx=Fraction(9*R,40)-Fraction(tx,D)
                    if xx.denominator==1:
                        n=xx.numerator
                        answers.append((n,n));STATS['exact_boundary_recognitions']+=1
                    else:answers.append(None)
                    continue
            bottom,top=point(-mp.pi/2)[1],point(mp.pi/2)[1]
            if yy<bottom or yy>top:
                answers.append(None);continue
            roots=[]
            for left,right,increasing in ((-mp.pi/2,mp.pi/2,True),
                                          (mp.pi/2,3*mp.pi/2,False)):
                for _ in range(4*dps):
                    mid=(left+right)/2
                    if (point(mid)[1]<yy)==increasing:left=mid
                    else:right=mid
                alpha=(left+right)/2
                xx=R*point(alpha)[0]-ux
                n=int(mp.nint(xx))
                if abs(xx-n)<mp.power(10,-dps//2):
                    physicalx=Fraction(D*n+tx,D*R)
                    physicaly=Fraction(D*row+ty,D*R)
                    if q==0 and _aligned_rational_boundary(physicalx,physicaly,alpha):
                        xx=mp.mpf(n);STATS['exact_boundary_recognitions']+=1
                    else:
                        raise ArithmeticError("Unresolved near-exact boundary; no silent classification")
                roots.append(xx)
            # Right half-boundary first, then left half-boundary.
            lo,hi=int(mp.ceil(roots[1])),int(mp.floor(roots[0]))
            answers.append((lo,hi) if lo<=hi else None)
    if answers[0]!=answers[1]:raise ArithmeticError("80/120-digit row endpoints disagree")
    return answers[1]


def initial_rows(R,p=1,q=0,tx=0,ty=0,D=1):
    """Inclusive integer row intervals for the analytic smooth initial shape.

    dy/dalpha=rho(alpha-theta)*cos(alpha), hence each half-boundary is
    monotone. No polygonal approximation is used in the membership route.
    """
    _parameters(R,p,q,tx,ty,D)
    ux,uy=tx/D,ty/D
    bottom=R*float(boundary(-math.pi/2,p,q)[1])
    top=R*float(boundary(math.pi/2,p,q)[1])
    rows=np.arange(math.floor(bottom-uy)-1,math.ceil(top-uy)+2,dtype=np.int64)
    yy=rows.astype(float)+uy
    guard=1024*np.finfo(float).eps*(R+1)
    keep=(yy>=bottom-guard)&(yy<=top+guard)
    rows,yy=rows[keep],yy[keep]
    targets=np.clip(yy/R,bottom/R,top/R)
    roots=[]
    for l0,r0,increasing in ((-math.pi/2,math.pi/2,True),
                              (math.pi/2,3*math.pi/2,False)):
        left=np.full(targets.shape,l0);right=np.full(targets.shape,r0)
        for _ in range(58):
            mid=(left+right)/2
            ym=boundary(mid,p,q)[1]
            select=ym<targets if increasing else ym>targets
            left=np.where(select,mid,left);right=np.where(select,right,mid)
        roots.append(R*boundary((left+right)/2,p,q)[0])
    xr,xl=roots
    lo=np.ceil(xl-ux).astype(np.int64);hi=np.floor(xr-ux).astype(np.int64)
    close=((np.abs((xl-ux)-np.rint(xl-ux))<1e-8)|
           (np.abs((xr-ux)-np.rint(xr-ux))<1e-8)|
           (np.abs(yy-bottom)<guard)|(np.abs(yy-top)<guard))
    valid=lo<=hi
    for j in np.flatnonzero(close):
        ends=_mp_section(R,p,q,tx,ty,D,int(rows[j]))
        if ends is None:valid[j]=False
        else:lo[j],hi[j]=ends;valid[j]=True
    return {int(row):[int(a),int(b)] for row,a,b,ok in zip(rows,lo,hi,valid) if ok}


def support_margin_oracle(x,y,p=1,q=0):
    """Independent stationary-polynomial support test, numerically qualified.

    The extrema of h(alpha)-x*cos(alpha)-y*sin(alpha) are found from all
    unit-circle roots of a degree-six polynomial. This does not use the
    horizontal-section bisection route. Root finding is not certified.
    """
    theta=math.atan2(q,p);a=-.075;m=3
    coeff=np.zeros(2*m+1,dtype=complex)
    coeff[2*m]=-a*m*np.exp(-1j*m*theta)
    coeff[m+1]=x-1j*y
    coeff[m-1]=-x-1j*y
    coeff[0]=a*m*np.exp(1j*m*theta)
    roots=np.polynomial.polynomial.polyroots(coeff)
    angles=[float(np.angle(z)) for z in roots if abs(abs(z)-1)<2e-5]
    if not angles:raise ArithmeticError("No unit-circle stationary root")
    return min(float(support(alpha,p,q))-x*math.cos(alpha)-y*math.sin(alpha) for alpha in angles)


def full_grid_oracle(R,p,q,tx,ty,D):
    bound=math.ceil(43*R/40)+2+max(abs(tx)//D,abs(ty)//D)
    out=set()
    for row in range(-bound,bound+1):
        for i in range(-bound,bound+1):
            margin=support_margin_oracle((i+tx/D)/R,(row+ty/D)/R,p,q)
            if abs(margin)<1e-9:
                STATS['ambiguous_oracle_shared_rechecks']+=1
                ends=_mp_section(R,p,q,tx,ty,D,row)
                valid=ends is not None and ends[0]<=i<=ends[1]
            else:valid=margin>=0
            if valid:out.add((i,row))
    return out


def cross(o,a,b):
    return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])


def hull(points):
    pts=sorted({(int(p[0]),int(p[1])) for p in points})
    if len(pts)<3:raise ValueError("Degenerate lattice hull")
    lower=[];upper=[]
    for p in pts:
        while len(lower)>1 and cross(lower[-2],lower[-1],p)<=0:lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper)>1 and cross(upper[-2],upper[-1],p)<=0:upper.pop()
        upper.append(p)
    return lower[:-1]+upper[:-1]


def endpoints(rows):return [(i,j) for j,(lo,hi) in rows.items() for i in (lo,hi)]


def remove_vertices(rows,vertices):
    for i,j in vertices:
        lo,hi=rows[j]
        if i==lo:lo+=1
        elif i==hi:hi-=1
        else:raise AssertionError("Extreme vertex is not a row endpoint")
        if lo>hi:del rows[j]
        else:rows[j]=[lo,hi]


def perimeter(vertices):
    return math.fsum(math.hypot(a[0]-b[0],a[1]-b[1]) for a,b in zip(vertices,vertices[1:]+vertices[:1]))


def area(vertices):
    return abs(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(vertices,vertices[1:]+vertices[:1])))/2


def qualification():
    initial=[];layers=[]
    orientations=((1,0),(8,1),(1,1))
    translations=((0,0,1),(271,419,1000),(1,7,8))
    for R in (3,5,8):
        for p,q in orientations:
            for tx,ty,D in translations:
                rows=initial_rows(R,p,q,tx,ty,D)
                full=full_grid_oracle(R,p,q,tx,ty,D)
                reconstructed={(i,j) for j,(lo,hi) in rows.items() for i in range(lo,hi+1)}
                assert reconstructed==full,(R,p,q,tx,ty,D)
                initial.append([R,p,q,tx,ty,D,len(full)])
                for n in range(6):
                    if len(full)<3:break
                    expected=hull(full)
                    if len(expected)<3:break
                    arr=np.asarray(sorted(full),dtype=np.int64)
                    qh=[tuple(map(int,v)) for v in arr[ConvexHull(arr).vertices]]
                    got=hull(endpoints(rows))
                    assert set(got)==set(expected)==set(qh)
                    assert sum(hi-lo+1 for lo,hi in rows.values())==len(full)
                    layers.append([R,p,q,tx,ty,D,n,len(got)])
                    remove_vertices(rows,got);full.difference_update(got)
    witnesses=[]
    for c,s in ((Fraction(1),Fraction(0)),(Fraction(0),Fraction(1)),
                 (Fraction(3,5),Fraction(4,5))):
        cos3=c*c*c-3*c*s*s;sin3=3*c*c*s-s*s*s
        h=1+AMPLITUDE*cos3;hp=-3*AMPLITUDE*sin3
        x,y=h*c-hp*s,h*s+hp*c
        R=math.lcm(x.denominator,y.denominator)
        section=_mp_section(R,1,0,0,0,1,int(y*R))
        assert section is not None and int(x*R) in section
        assert _aligned_rational_boundary(x,y,math.atan2(float(s),float(c)))
        witnesses.append({'normal':[str(c),str(s)],'R':R,'point':[int(x*R),int(y*R)]})
    # Pure integer peeling retains non-vertex collinear boundary points.
    rows={-1:[-1,1],0:[-1,1],1:[-1,1]}
    vertices=hull(endpoints(rows));remove_vertices(rows,vertices)
    assert {(i,j) for j,(lo,hi) in rows.items() for i in range(lo,hi+1)}=={(0,0),(-1,0),(1,0),(0,-1),(0,1)}
    # Read-only anchor against the previous smooth-support engine.
    previous=ROOT/'baseline'/'support_shapes_initial.py'
    spec=importlib.util.spec_from_file_location('support_shapes_initial',previous)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    anchors=[]
    for R in (8,64,128,256):
        for p,q in orientations:
            tx,ty,D=271,419,1024
            got=hull(endpoints(initial_rows(R,p,q,tx,ty,D)))
            theta=math.atan2(q,p)
            reference=old.hull(old.row_points(3,.6,R,theta,tx/D,ty/D))
            assert got==reference,(R,p,q)
            anchors.append([R,p,q,len(got)])
    return {'status':'PASS_QUALIFIED_SMOOTH_SUPPORT_GEOMETRY',
        'global_interval_certified':False,'initial_full_grid_comparisons':len(initial),
        'layer_full_grid_qhull_comparisons':len(layers),'exact_rational_boundary_witnesses':witnesses,
        'historical_initial_hull_anchors':len(anchors),'collinear_nonvertices_retained':True,
        'geometry':geometry(),'initial_details':initial,'layer_details':layers,'anchor_details':anchors,
        'stats':dict(STATS),
        'oracle_scope':'stationary support polynomial independent of row bisection; ambiguous boundary shares high-precision section route',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'baseline_source_sha256':hashlib.sha256(previous.read_bytes()).hexdigest()}


if __name__=='__main__':
    result=qualification()
    (ROOT/'geometry_controls.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','global_interval_certified',
        'initial_full_grid_comparisons','layer_full_grid_qhull_comparisons',
        'historical_initial_hull_anchors','stats')}))
