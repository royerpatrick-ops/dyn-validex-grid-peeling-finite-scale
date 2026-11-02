#!/usr/bin/env python3
"""Conditional one-step integer geometry predictor, not a macroscopic flow law.

Construction reads only the current strictly convex CCW integer polygon.
Every retained integer point is represented by the lowest/highest survivor in
its integer-x column; their exact integer hull is the next strict-vertex layer.
It is one application of the peeling map, rather than a new trajectory/campaign.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, sys, time
import numpy as np

HERE=Path(__file__).resolve().parent
# Adjacent modules in a delivered package; original workspace during extraction.
ROOT=HERE.parent
PACKAGE=HERE if (HERE/'polygon_localization.py').is_file() else ROOT/'output/pi4_displacements_v01'
sys.path.insert(0,str(PACKAGE))
from polygon_localization import PolygonSupport

def cross(a,b,c):return (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0])

def strict_integer_hull(points):
    points=sorted(set(points))
    if len(points)<3:raise ValueError('Fewer than three surviving candidates')
    lo=[]
    for p in points:
        while len(lo)>=2 and cross(lo[-2],lo[-1],p)<=0:lo.pop()
        lo.append(p)
    hi=[]
    for p in reversed(points):
        while len(hi)>=2 and cross(hi[-2],hi[-1],p)<=0:hi.pop()
        hi.append(p)
    hull=lo[:-1]+hi[:-1]
    if len(hull)<3:raise ValueError('Degenerate survivor hull')
    return np.asarray(hull,dtype=np.int64)

def reconstruct_next(vertices,return_stats=False):
    """Return exact next hull from a saturated current integer polygon only.

    No reference to future vertices, recorded supports, or the peeling engine.
    All incidence, intersection, rounding and hull predicates use Python ints.
    """
    raw=np.asarray(vertices)
    if raw.ndim!=2 or raw.shape[1]!=2 or len(raw)<3 or not np.issubdtype(raw.dtype,np.integer):
        raise ValueError('Expected >=3 two-dimensional integer vertices')
    vv=[(int(x),int(y)) for x,y in raw]
    if any(cross(vv[i-1],vv[i],vv[(i+1)%len(vv)])<=0 for i in range(len(vv))):
        raise ValueError('Current vertices must be strict CCW')
    vertices_set=set(vv);edges=list(zip(vv,vv[1:]+vv[:1]))
    candidates=[];count=0;columns=0
    for x in range(min(p[0] for p in vv),max(p[0] for p in vv)+1):
        lower=upper=None
        for a,b in edges:
            ax,ay=a;bx,by=b
            if not min(ax,bx)<=x<=max(ax,bx):continue
            dx=bx-ax
            if dx==0:
                ys=[(ay,1),(by,1)]
            else:
                num=ay*dx+(by-ay)*(x-ax);den=dx
                if den<0:num,den=-num,-den
                ys=[(num,den)]
            for num,den in ys:
                if lower is None or num*lower[1]<lower[0]*den:lower=(num,den)
                if upper is None or num*upper[1]>upper[0]*den:upper=(num,den)
        if lower is None:continue
        ylo=-((-lower[0])//lower[1]);yhi=upper[0]//upper[1]
        while ylo<=yhi and (x,ylo) in vertices_set:ylo+=1
        while ylo<=yhi and (x,yhi) in vertices_set:yhi-=1
        if ylo<=yhi:
            candidates.append((x,ylo))
            if yhi!=ylo:candidates.append((x,yhi))
            count+=yhi-ylo+1;columns+=1
    hull=strict_integer_hull(candidates)
    stats={'surviving_integer_points':count,'nonempty_columns':columns,'candidate_points':len(candidates),'predicted_vertices':len(hull)}
    return (hull,stats) if return_stats else hull

def digest(vertices):
    z=np.asarray(vertices,dtype='<i8');z=z[np.lexsort((z[:,1],z[:,0]))]
    return hashlib.sha256(z.tobytes()).hexdigest()

def qualify():
    # Independent exhaustive point enumeration on small convex polygons.
    fixtures=[[(0,0),(4,0),(4,4),(0,4)],[(0,0),(5,0),(0,5)],[(0,0),(2,-1),(5,1),(4,4),(1,5),(-1,2)],[(0,0),(1,0),(3,2),(2,4),(-1,3),(-2,1)]]
    for fixture in fixtures:
        v=strict_integer_hull(fixture);vv=[tuple(map(int,p)) for p in v];vs=set(vv)
        keep=[]
        for x in range(min(p[0] for p in vv),max(p[0] for p in vv)+1):
            for y in range(min(p[1] for p in vv),max(p[1] for p in vv)+1):
                p=(x,y)
                if p not in vs and all(cross(a,b,p)>=0 for a,b in zip(vv,vv[1:]+vv[:1])):keep.append(p)
        expected=strict_integer_hull(keep);actual,stats=reconstruct_next(v,True)
        assert digest(expected)==digest(actual)
        assert len(keep)==stats['surviving_integer_points']
    return {'small_exact_exhaustive_fixtures':len(fixtures),'status':'PASS'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--production',type=Path,default=ROOT/'work/production_v04/_prepare_replication_v0_4_20261027')
    parser.add_argument('--output',type=Path,default=ROOT/'tmp/one_step_prediction_results')
    args=parser.parse_args();start=time.perf_counter();qualification=qualify()
    alpha=.6450983601211226;predictions=[];polygons=[];current_polygons=[];records=[]
    # Stage 1: build and seal all forecasts before any next-layer comparison.
    for row in (0,8,16):
        with np.load(args.production/f'campaign_v04/rows/row_{row:02d}_main.npz',allow_pickle=False) as z:
            lookup={(int(j),int(n)):i for i,(j,n) in enumerate(zip(z['layer_trajectory'],z['layer_n']))}
            flat=z['layer_vertices_ordered_flat'];offsets=z['layer_vertex_offsets']
            for j in range(6):
                R=int(z['radii'][j]);theta=float(z['orientation_theta'][j]);tr=z['translations'][j].copy();translation=tr[:2]/tr[2]
                for stage,n in [('initial',0),('t002',math.floor(.02*R**(4/3)/alpha)),('t004',math.floor(.04*R**(4/3)/alpha))]:
                    i=lookup[j,n];p,l=map(int,offsets[i]);current=flat[p:p+l].copy()
                    predicted,stats=reconstruct_next(current,True)
                    cprev=PolygonSupport(current,theta,translation,R).coefficients((0,3))
                    cpred=PolygonSupport(predicted,theta,translation,R).coefficients((0,3))
                    assert stats['surviving_integer_points']==int(z['layer_remaining_points'][i])-len(current)
                    rec={'row':row,'trajectory':j,'R':R,'theta':theta,'translation':tr.tolist(),'stage':stage,'current_n':n,'predicted_n':n+1,'current_sha256':digest(current),'predicted_sha256':digest(predicted),**stats,
                         'predicted_c0':cpred[0].real,'predicted_c3_real':cpred[3].real,'predicted_c3_imag':cpred[3].imag,'predicted_a3_real':(cpred[3]/cpred[0].real).real,'predicted_a3_imag':(cpred[3]/cpred[0].real).imag,
                         'predicted_delta_c3_real':(cpred[3]-cprev[3]).real,'predicted_delta_c3_imag':(cpred[3]-cprev[3]).imag,'predicted_delta_a3_real':(cpred[3]/cpred[0].real-cprev[3]/cprev[0].real).real,'predicted_delta_a3_imag':(cpred[3]/cpred[0].real-cprev[3]/cprev[0].real).imag}
                    records.append(rec);polygons.append(predicted);current_polygons.append(current)
    seal=hashlib.sha256(json.dumps(records,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    # Stage 2: independent reads of stored next-layer outcomes.
    max_c0=max_c3=max_a3=0.;geometry_failures=0;outcome=[]
    for row in (0,8,16):
        with np.load(args.production/f'campaign_v04/rows/row_{row:02d}_main.npz',allow_pickle=False) as z:
            lookup={(int(j),int(n)):i for i,(j,n) in enumerate(zip(z['layer_trajectory'],z['layer_n']))}
            for q in [q for q in records if q['row']==row]:
                i=lookup[q['trajectory'],q['predicted_n']];p,l=map(int,z['layer_vertex_offsets'][i]);actual=z['layer_vertices_ordered_flat'][p:p+l]
                translation=np.asarray(q['translation'][:2])/q['translation'][2]
                c=PolygonSupport(actual,q['theta'],translation,q['R']).coefficients((0,3));a=c[3]/c[0].real
                dc0=abs(c[0].real-q['predicted_c0']);dc3=abs(c[3]-complex(q['predicted_c3_real'],q['predicted_c3_imag']));da3=abs(a-complex(q['predicted_a3_real'],q['predicted_a3_imag']))
                equal=digest(actual)==q['predicted_sha256'];geometry_failures+=not equal
                max_c0=max(max_c0,dc0);max_c3=max(max_c3,dc3);max_a3=max(max_a3,da3)
                outcome.append({**q,'actual_sha256':digest(actual),'integer_geometry_equal':equal,'c0_absolute_error':dc0,'c3_absolute_error':dc3,'a3_absolute_error':da3})
    assert geometry_failures==0
    result={'status':'PASS_CONDITIONAL_ONE_STEP_INTEGER_GEOMETRY','qualification':qualification,'predictions':len(records),'forecast_seal_sha256':seal,'integer_geometry_failures':geometry_failures,'max_c0_absolute_error':max_c0,'max_c3_absolute_error':max_c3,'max_a3_absolute_error':max_a3,'new_independent_trajectory':False,'new_parameter_fits':0,'new_pde_solve':False,'interpretation':'Exact conditional one-step evaluation of the existing peeling map; not an autonomous macroscopic prediction from initial mode or 16 normals. All 54 predictions were constructed and sealed before the validation stage accessed next-layer polygons.','wall_seconds':time.perf_counter()-start,'records':outcome}
    args.output.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    co=[];po=[];c=0;p=0
    for a,b in zip(current_polygons,polygons):co.append([c,len(a)]);po.append([p,len(b)]);c+=len(a);p+=len(b)
    np.savez_compressed(args.output.with_suffix('.npz'),current_vertices=np.concatenate(current_polygons),current_offsets=np.asarray(co),predicted_vertices=np.concatenate(polygons),predicted_offsets=np.asarray(po))
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))

if __name__=='__main__':main()
