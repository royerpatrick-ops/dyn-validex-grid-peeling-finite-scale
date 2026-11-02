#!/usr/bin/env python3
"""Diagnostic only: sixteen observed directional supports do not specify a3.

No forecast is claimed: reconstruct the outer polygon determined by the exact
sixteen supports at saved marks and compare its third support mode with the
saved polygon. No peeling, new draw, PDE solve, or fitting is performed.
"""
from pathlib import Path
import sys, json, math, argparse
from collections import defaultdict
import numpy as np
from scipy.spatial import ConvexHull

ROOT=Path(__file__).resolve().parent
BASE=ROOT
sys.path.insert(0,str(BASE))
from polygon_localization import PolygonSupport
from frozen_reference_time import FrozenTimeReference

def outer_polygon(directions, levels):
    pts=[]
    for i in range(len(directions)):
        for j in range(i):
            a,b=directions[i],directions[j]
            det=float(a[0]*b[1]-a[1]*b[0])
            if det==0:continue
            x=np.array([(levels[i]*b[1]-a[1]*levels[j])/det,
                        (a[0]*levels[j]-levels[i]*b[0])/det])
            if np.max(directions@x-levels)<1e-8:pts.append(x)
    pts=np.unique(np.round(np.asarray(pts),10),axis=0)
    return pts[ConvexHull(pts).vertices]

def coefficients(vertices, theta, trans, R):
    xy=(np.asarray(vertices,dtype=float)+trans)/R
    z=(xy[:,0]+1j*xy[:,1])*np.exp(-1j*theta)
    edge=np.roll(z,-1)-z
    upper=np.mod(np.angle(-1j*edge),2*np.pi)
    lower=np.roll(upper,1)
    width=np.mod(upper-lower,2*np.pi)
    assert abs(math.fsum(width)-2*np.pi)<1e-12
    def exp_integral(k):
        if k==0:return width.astype(complex)
        return 2*np.sin(k*width/2)/k*np.exp(-1j*k*(lower+width/2))
    out={}
    for m in [0,3]:
        values=(z.conjugate()*exp_integral(m-1)+z*exp_integral(m+1))/(4*np.pi)
        out[m]=complex(math.fsum(values.real),math.fsum(values.imag))
    return out

def linear_error_modes(angles, values):
    order=np.argsort(np.mod(angles,2*np.pi))
    a=np.mod(angles[order],2*np.pi);v=values[order]
    width=np.diff(np.r_[a,a[0]+2*np.pi]);nextv=np.roll(v,-1)
    c0=float(np.sum(width*(v+nextv)/2)/(2*np.pi))
    m=3;ex=np.exp(-1j*m*width)
    j0=(1-ex)/(1j*m)
    j1=-width*ex/(1j*m)+(ex-1)/(m*m)
    term=np.exp(-1j*m*a)*(v*j0+(nextv-v)/width*j1)
    c3=complex(math.fsum(term.real),math.fsum(term.imag))/(2*np.pi)
    return c0,c3

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'sparse_replay')
    out=parser.parse_args().output
    if out.exists():raise FileExistsError('Use a fresh output directory')
    out.mkdir()
    with np.load(BASE/'inputs/projected_layers.npz',allow_pickle=False) as z:
        data={k:z[k] for k in z.files}
    with np.load(BASE/'inputs/control_polygons.npz',allow_pickle=False) as z:
        control={k:z[k] for k in z.files}
    ref=FrozenTimeReference(BASE/'inputs/reference_contraction.npz')
    alpha=ref.alpha_clock
    directions=data['directions']
    records=[];max_sample_error=0.
    for case in range(len(data['case_R'])):
        R=int(data['case_R'][case]);theta=float(data['case_theta'][case])
        t=data['case_translation'][case];trans=t[:2]/t[2]
        start,count=map(int,data['case_offsets_start_count'][case])
        for mark,stage in enumerate(control['stage']):
            n=int(control['layer_n'][case,mark]);v0,nv=map(int,control['offsets_start_count'][case,mark])
            truth=control['vertices'][v0:v0+nv]
            lev=data['levels_integer'][start+n].astype(float)
            outer=outer_polygon(directions,lev)
            max_sample_error=max(max_sample_error,float(np.max(np.abs(np.max(outer@directions.T,axis=0)-lev)))/R)
            pc=coefficients(truth,theta,trans,R);qc=coefficients(outer,theta,trans,R)
            pa=pc[3]/pc[0].real;qa=qc[3]/qc[0].real
            rp=ref.at_q(pc[0].real)
            rq=ref.at_q(qc[0].real) if qc[0].real<=ref.q_max else None
            rt=ref.at_time(alpha*n/R**(4/3))
            angles=np.arctan2(directions[:,1],directions[:,0])-theta
            hp=(lev+directions@trans)/(R*np.linalg.norm(directions,axis=1))
            href=ref.raw_support(rt['coefficients'],angles)
            de0,de3=linear_error_modes(angles,hp-href)
            linear_c0=rt['c0']+de0
            linear_c3=rt['coefficients'][ref.index3]+de3
            linear_a3=linear_c3/linear_c0
            linear_a3_with_true_c0=linear_c3/pc[0].real
            error=qa-pa
            residual=pa.real-rp['a3'].real
            records.append({'case_index':case,'row_index':int(data['case_row'][case]),'R':R,'orientation':str(data['case_orientation'][case]),'stage':str(stage),'n':n,'time':alpha*n/R**(4/3),'true_c0':pc[0].real,'outer_c0':qc[0].real,'true_a3_real':pa.real,'true_a3_imag':pa.imag,'outer_a3_real':qa.real,'outer_a3_imag':qa.imag,'error_a3_real':error.real,'error_a3_imag':error.imag,'error_a3_abs':abs(error),'true_residual_equal_c0':residual,'outer_residual_using_true_ref':qa.real-rp['a3'].real,'outer_residual_using_own_ref':qa.real-rq['a3'].real if rq is not None else None,'linear_error_a3_real':(linear_a3-pa).real,'linear_error_a3_imag':(linear_a3-pa).imag,'linear_error_a3_true_c0_real':(linear_a3_with_true_c0-pa).real,'linear_c0_error':linear_c0-pc[0].real,'true_vertices':len(truth),'outer_vertices':len(outer)})
    groups=defaultdict(list)
    for r in records:groups[(r['R'],r['orientation'],r['stage'])].append(r)
    summary=[]
    for key,rr in groups.items():
        g={'R':key[0],'orientation':key[1],'stage':key[2],'cases':len(rr)}
        for f in ['true_a3_real','outer_a3_real','error_a3_real','error_a3_imag','error_a3_abs','true_residual_equal_c0','outer_residual_using_true_ref','true_c0','outer_c0','linear_error_a3_real','linear_error_a3_imag','linear_error_a3_true_c0_real','linear_c0_error']:
            x=[r[f] for r in rr]
            g['mean_'+f]=math.fsum(x)/len(x);g['min_'+f]=min(x);g['max_'+f]=max(x)
        g['wrong_sign_residual_count']=int(sum(np.sign(r['true_residual_equal_c0'])!=np.sign(r['outer_residual_using_true_ref']) for r in rr))
        g['error_exceeds_residual_count']=sum(abs(r['error_a3_real'])>abs(r['true_residual_equal_c0']) for r in rr)
        g['linear_error_exceeds_residual_count']=int(sum(abs(r['linear_error_a3_real'])>abs(r['true_residual_equal_c0']) for r in rr))
        summary.append(g)
    result={'status':'DIAGNOSTIC_NOT_FORECAST','case_count':384,'saved_marks':4,'profiles':len(records),'fixed_support_directions':directions.tolist(),'max_sample_support_discrepancy_physical':max_sample_error,'max_abs_a3_error':max(r['error_a3_abs'] for r in records),'groups':summary,'records':records,'qualification':'Sixteen supports are observed, not predicted. This constructs a distinct convex polygon with the same directional samples, to assess spatial information loss. The outer polygon is not generally a lattice polygon.'}
    (out/'angular16_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['records','groups']},indent=2))
    print('Sparse angular diagnostics written to',out/'angular16_audit.json')

if __name__=='__main__':main()
