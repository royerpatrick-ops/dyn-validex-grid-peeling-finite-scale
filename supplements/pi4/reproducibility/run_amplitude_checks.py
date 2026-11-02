#!/usr/bin/env python3
"""Replay the causal sufficient unit-jump certificates from extracted sections.

M and r0 are saved current-polygon measurements, not future outcomes. Recomputing
these geometric sections from original polygons needs the optional extractor.
This compact replay recomputes both certificates and verifies them against the
unchanged saved integer support trajectories. No fitted parameter.
"""
from pathlib import Path
import argparse, csv, hashlib, json
import numpy as np

ROOT=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=Path('amplitude_replay'))
    args=p.parse_args();out=args.output
    if not out.is_absolute():out=ROOT/out
    if out.exists():raise FileExistsError('Use a new output directory')
    out.mkdir(parents=True)
    with np.load(ROOT/'inputs/jump_birth_certificate_strong.npz',allow_pickle=False) as z:
        A=z['predictions'];r0=z['below_line_initial_vertex_count']
    case,n,d,k,m,M,L,stored_cert,stored_jump=A.T
    assert np.array_equal(L,(m+1)//2)
    strong=M>r0+2*(L-1);weak=M>2*L
    assert np.array_equal(strong,stored_cert.astype(bool))
    assert np.array_equal(np.where(strong,1,-1),stored_jump)
    assert np.all((r0>=0)&(r0<=2)&(r0<=M))
    # Write causal amplitude predictions before loading any future levels.
    predicted=out/'amplitude_predictions.npz'
    np.savez_compressed(predicted,case=case,birth_n=n,direction=d,L=L,
                        initial_certificate=weak,strengthened_certificate=strong,
                        predicted_integer_jump=np.where(strong,1,-1))
    seal=hashlib.sha256(predicted.read_bytes()).hexdigest()
    (out/'FORECAST_SEAL.json').write_text(json.dumps({'prediction_sha256':seal},indent=2)+'\n')
    with np.load(ROOT/'inputs/projected_layers.npz',allow_pickle=False) as z:
        levels=z['levels_integer'];faces=z['face_points'];offsets=z['case_offsets_start_count']
        R=z['case_R'];ori=z['case_orientation']
    seen=n+L<offsets[case,1]
    pos=offsets[case,0]+n
    assert np.array_equal(levels[pos,d],k)
    assert np.array_equal(faces[pos,d],m)
    actual=np.full(len(A),-1,dtype=np.int64)
    actual[seen]=k[seen]-levels[pos[seen]+L[seen],d[seen]]
    groups=[]
    for radius in (256,512):
        for angle in ('theta0','atan_1_8','pi_over_4'):
            for kind,maskkind in [('informative',m>=3),('trivial',m<=2)]:
                mask=(R[case]==radius)&(ori[case]==angle)&maskkind
                for rule,cert in [('initial',weak),('strengthened',strong)]:
                    groups.append({'R':radius,'orientation':angle,'kind':kind,'rule':rule,
                      'births':int(mask.sum()),'certified':int((mask&cert).sum()),
                      'certified_complete':int((mask&cert&seen).sum()),
                      'certified_censored':int((mask&cert&~seen).sum()),
                      'abstentions':int((mask&~cert).sum()),
                      'failures':int((mask&cert&seen&(actual!=1)).sum())})
    failures=int((strong&seen&(actual!=1)).sum())
    result={'status':'PASS_AMPLITUDE_CERTIFICATES' if failures==0 else 'FAIL',
       'births':len(A),'complete_observed_jumps':int(seen.sum()),
       'nonunit_complete_observed_jumps':int((seen&(actual!=1)).sum()),
       'certified_complete':int((strong&seen).sum()),'certified_censored':int((strong&~seen).sum()),
       'certified_failures':failures,'informative_certified':int((strong&(m>=3)).sum()),
       'informative_births':int((m>=3).sum()),'informative_abstentions':int((~strong&(m>=3)).sum()),
       'fitted_parameters':0,'prediction_sha256':seal,
       'scope':'Replay from extracted current sections; geometric extraction separately qualified, not rerun here.',
       'adaptivity':'The strengthened theorem was selected after inspecting coverage of the initial certificate; no coefficient was fitted.',
       'groups':groups}
    with (out/'amplitude_groups.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(groups[0]));w.writeheader();w.writerows(groups)
    (out/'amplitude_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'complete',result['certified_complete'],'censored',result['certified_censored'],'failures',failures)
    return int(failures>0)

if __name__=='__main__':raise SystemExit(main())
