#!/usr/bin/env python3
"""Audit saved prospective replication geometry; no new peeling.

This program imports no delivered module. It validates stored polygons against
the exact row intervals reconstructed from their initial integer hull, and
deletes only the stored extreme vertices between successive stored layers.
It never constructs a hull for a new point set and never draws coordinates.
Initial membership in the original smooth body is outside this audit.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def digest_vertices(v):
    a=np.unique(np.asarray(v,dtype=np.int64),axis=0)
    a=a[np.lexsort((a[:,1],a[:,0]))].astype('<i8',copy=False)
    return hashlib.sha256(a.tobytes()).hexdigest()


def perimeter(v):
    return math.fsum(math.hypot(int(a[0])-int(b[0]),int(a[1])-int(b[1]))
                     for a,b in zip(v,np.roll(v,-1,axis=0)))


def rows_from_integer_hull(v):
    # Every integer point of the initial lattice hull was in the original
    # lattice set: hull subset K and the set was K intersect the lattice.
    rows={}
    for y in range(int(v[:,1].min()),int(v[:,1].max())+1):
        lower=int(v[:,0].min()); upper=int(v[:,0].max())
        for a,b in zip(v,np.roll(v,-1,axis=0)):
            ax,ay,bx,by=map(int,(*a,*b)); dx,dy=bx-ax,by-ay
            numerator=dx*(y-ay)+dy*ax
            if dy>0:upper=min(upper,numerator//dy)
            elif dy<0:lower=max(lower,-((-numerator)//dy))
            elif dx*(y-ay)<0:lower,upper=1,0;break
        if lower<=upper:rows[y]=[lower,upper]
    return rows


def coefficients(v,translation,radius,modes):
    # Independent direct integration of h(phi) on normal cones; handles m=1.
    vf=np.asarray(v,dtype=float)+translation
    edges=np.roll(v,-1,axis=0)-v
    normals=np.unwrap(np.arctan2(-edges[:,0],edges[:,1]))
    lower=np.r_[normals[-1]-2*math.pi,normals[:-1]]
    upper=normals
    z=vf[:,0]+1j*vf[:,1]
    out=[]
    for mode in modes:
        integrals=[]
        for k in [int(mode)-1,int(mode)+1]:
            integrals.append(upper-lower if k==0 else
                             (np.exp(-1j*k*upper)-np.exp(-1j*k*lower))/(-1j*k))
        value=np.sum(np.conj(z)*integrals[0]/2+z*integrals[1]/2)/(2*math.pi*radius)
        out.append(value)
    return np.array(out)


def audit_line(path,check_path):
    with np.load(path,allow_pickle=False) as archive:
        z={name:archive[name] for name in archive.files}
    check=json.loads(check_path.read_text())
    nlayer=len(z['layer_n']); maxp=0.; maxcoef=0.; maxinterpolation=0.
    assert nlayer==check['layers_compared']
    assert z['radii'].tolist()==[256,256,256,512,512,512]
    assert z['orientation_pq'].tolist()==[[1,0],[8,1],[1,1],[1,0],[8,1],[1,1]]
    digest=hashlib.sha256(); layercounts=[]; projected_jarvis_visits=0
    for ti in range(6):
        indices=np.flatnonzero(z['layer_trajectory']==ti)
        assert z['layer_n'][indices].tolist()==list(range(len(indices)))
        layercounts.append(len(indices)); rows=None; previous=None
        pvalues=[]
        for li in indices:
            start,count=map(int,z['layer_vertex_offsets'][li])
            v=z['layer_vertices_ordered_flat'][start:start+count]
            assert len(v)==z['layer_vertex_count'][li]
            assert digest_vertices(v)==str(z['layer_vertex_sha256'][li])
            edges=np.roll(v,-1,axis=0)-v
            turns=edges[:,0]*np.roll(edges[:,1],-1)-edges[:,1]*np.roll(edges[:,0],-1)
            assert (turns>0).all()
            assert max(abs(int(v.min())),abs(int(v.max())))<2048
            if rows is None:
                rows=rows_from_integer_hull(v)
                assert digest_vertices(v)==str(z['initial_vertex_sha256'][ti])
            assert sum(hi-lo+1 for lo,hi in rows.values())==int(z['layer_remaining_points'][li])
            endpoints=np.array([(x,y) for y,(lo,hi) in rows.items() for x in (lo,hi)],dtype=np.int64)
            # Coordinates are all below 2048: these integer crosses cannot overflow.
            cross=edges[:,None,0]*(endpoints[None,:,1]-v[:,None,1])-edges[:,None,1]*(endpoints[None,:,0]-v[:,None,0])
            assert (cross>=0).all(), 'stored hull fails to contain a remaining row endpoint'
            for x,y in v:
                lo,hi=rows[int(y)]
                assert int(x) in (lo,hi)
            currentp=perimeter(v);maxp=max(maxp,abs(currentp-float(z['layer_perimeter'][li])))
            pvalues.append(currentp)
            if previous is not None:
                assert currentp<perimeter(previous)
                pe=np.roll(previous,-1,axis=0)-previous
                c=pe[:,None,0]*(v[None,:,1]-previous[:,None,1])-pe[:,None,1]*(v[None,:,0]-previous[:,None,0])
                assert int(c.min())==int(z['layer_minimum_containment_cross_from_previous'][li])>=0
            record={'R':int(z['radii'][ti]),'p':int(z['orientation_pq'][ti,0]),'q':int(z['orientation_pq'][ti,1]),
                    'tx':int(z['translations'][ti,0]),'ty':int(z['translations'][ti,1]),
                    'denominator':int(z['translations'][ti,2]),'n':int(z['layer_n'][li]),
                    'remaining_points':int(z['layer_remaining_points'][li]),'vertex_count':int(count),
                    'perimeter_hex':float(z['layer_perimeter'][li]).hex(),'vertex_sha256':str(z['layer_vertex_sha256'][li]),
                    'minimum_containment_cross_from_previous':None if previous is None else int(z['layer_minimum_containment_cross_from_previous'][li])}
            digest.update((json.dumps(record,sort_keys=True,separators=(',',':'))+'\n').encode('ascii'))
            projected_jarvis_visits+=int(z['layer_remaining_points'][li])*count
            # Validate the actual saved removal sequence; no hull is recomputed.
            for x,y in v:
                x,y=int(x),int(y);lo,hi=rows[y]
                if x==lo:lo+=1
                elif x==hi:hi-=1
                else:raise AssertionError('non-endpoint deletion')
                if lo>hi:del rows[y]
                else:rows[y]=[lo,hi]
            previous=v
        for qi,target in enumerate(z['target_q_historical']):
            first=int(np.flatnonzero(np.asarray(pvalues)<=target*pvalues[0])[0])
            assert z['endpoint_layers'][ti,qi].tolist()==[first-1,first]
            pair=z['endpoint_perimeters'][ti,qi]
            assert pvalues[first-1]>target*pvalues[0]>=pvalues[first]
            assert abs(pair[0]-pvalues[first-1])<1e-11 and abs(pair[1]-pvalues[first])<1e-11
            w=(pair[0]-target*pvalues[0])/(pair[0]-pair[1])
            assert abs(w-z['weights_after'][ti,qi])<2e-12
            qabs=target*pvalues[0]/(2*math.pi*int(z['radii'][ti]))
            assert abs(qabs-z['target_q_absolute'][ti,qi])<1e-14
            cp=[]
            for end in range(2):
                start,count=map(int,z['endpoint_vertex_offsets'][ti,qi,end])
                vv=z['endpoint_vertices_ordered_flat'][start:start+count]
                li=indices[first-1+end];ls,lc=map(int,z['layer_vertex_offsets'][li])
                assert np.array_equal(vv,z['layer_vertices_ordered_flat'][ls:ls+lc])
                translation=z['translations'][ti,:2]/z['translations'][ti,2]
                coeff=coefficients(vv,translation,int(z['radii'][ti]),z['fourier_modes'])
                maxcoef=max(maxcoef,float(abs(coeff-z['endpoint_coefficients_raw'][ti,qi,end]).max()))
                cp.append(coeff)
            inter=(1-w)*cp[0]+w*cp[1]
            maxinterpolation=max(maxinterpolation,float(abs(inter-z['interpolated_coefficients_raw'][ti,qi]).max()))
            inter[0]=qabs
            global_a=inter/qabs
            shape_a=global_a*np.exp(1j*z['fourier_modes']*z['orientation_theta'][ti])
            assert np.max(abs(shape_a-z['a_shape_exp_plus_i_m_theta'][ti,qi]))<2e-13
            bound=2*w*(1-w)*(pair[0]-pair[1])/(2*math.pi*z['radii'][ti])
            assert abs(bound-z['interpolation_coefficient_bounds'][ti,qi])<1e-14
            assert abs(bound/qabs-z['epsilon_interpolation_a3'][ti,qi])<1e-14
    assert digest.hexdigest()==check['cumulative_layer_record_sha256_primary']==check['cumulative_layer_record_sha256_secondary']
    assert maxp<1e-10 and maxcoef<2e-13 and maxinterpolation<2e-13
    return {'line':path.name,'layers_checked':nlayer,'trajectory_layer_counts':layercounts,
            'uncompressed_array_bytes':sum(a.nbytes for a in z.values()),
            'max_perimeter_abs_error':maxp,'max_endpoint_fourier_abs_error':maxcoef,
            'max_interpolated_fourier_abs_error':maxinterpolation,
            'cumulative_digest_reproduced':digest.hexdigest(),
            'projected_jarvis_visits_recomputed':projected_jarvis_visits,
            'remaining_row_sets_validated_from_initial_hull':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.campaign_root.resolve();p=root/'rows'
    indices=[0,8,16,24,32,40,48,56]
    results=[]
    inputs={}
    for i in indices:
        main_path=p/f'row_{i:02d}_main.npz'
        check_path=p/f'row_{i:02d}_quickhull_check.json'
        result=audit_line(main_path,check_path)
        result['row_index']=i
        results.append(result)
        for path in [main_path,check_path]:
            inputs[str(path.relative_to(root))]=hashlib.sha256(path.read_bytes()).hexdigest()
        print(json.dumps({'row_index':i,'layers_checked':result['layers_checked'],'fourier_error':result['max_endpoint_fourier_abs_error']}),flush=True)
    out={
      'status':'SAVED_EIGHT_PRESPECIFIED_ROWS_GEOMETRY_REPRODUCED',
      'scope':'No new coordinates, seed, hull construction or peeling. All saved layer records on the eight predeclared QuickHull-check rows were validated using independent exact-integer row reconstruction and direct Fourier integration.',
      'row_indices':indices,
      'rows_checked':len(results),
      'total_layers_checked':sum(r['layers_checked'] for r in results),
      'maximum_perimeter_abs_error':max(r['max_perimeter_abs_error'] for r in results),
      'maximum_endpoint_fourier_abs_error':max(r['max_endpoint_fourier_abs_error'] for r in results),
      'maximum_interpolated_fourier_abs_error':max(r['max_interpolated_fourier_abs_error'] for r in results),
      'rows':results,
      'input_sha256':inputs,
      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'limits':[
        'Only eight predeclared rows, not all64 main rows, received this full layer geometry check.',
        'The initial lattice set is reconstructed from each initial integer hull. Membership in the original smooth body was not independently requalified.',
        'QuickHull raw secondary polygons are not separately saved; equality reports are checked by reproducing their cumulative primary and secondary digest from main saved polygons.',
        'This check does not execute either delivered geometry engine and does not create any new scientific observations.',
        'Numerical Fourier identities are checked in floating point; no interval enclosure is claimed.'
      ]
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'layers':out['total_layers_checked'],'max_fourier_error':out['maximum_endpoint_fourier_abs_error']}))

if __name__=='__main__':main()
