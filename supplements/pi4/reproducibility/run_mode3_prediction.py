#!/usr/bin/env python3
"""Replay conditional global-mode forecasts from current-state compact inputs.

All forecasts and their seal are written before the outcome benchmark is read.
This exact one-step microscopic map does not predict initial-to-target dynamics.
"""
from pathlib import Path
import argparse, csv, hashlib, json, math, sys, time
import numpy as np
from one_step_prediction import reconstruct_next, digest, qualify, PACKAGE
sys.path.insert(0,str(PACKAGE))
from polygon_localization import PolygonSupport,dipole_weighted_integrals
from frozen_reference_time import FrozenTimeReference

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write_csv(path,rows):
    with path.open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run(inputs,benchmark,reference,output):
    if output.exists():raise FileExistsError('Use a new output directory')
    output.mkdir(parents=True);started=time.perf_counter();qtest=qualify()
    with np.load(inputs,allow_pickle=False) as z:
        flat=z['current_vertices'];offsets=z['current_offsets'];metadata=json.loads(str(z['metadata_json']))
    ref=FrozenTimeReference(reference)
    intervals=np.column_stack((np.deg2rad(np.arange(360)),np.deg2rad(np.arange(1,361))))
    sectors=np.array([[-math.pi/36,math.pi/36],[math.pi-math.pi/36,math.pi+math.pi/36],[0,2*math.pi]])
    rows=[];resbins=[];centeredbins=[];stepbins=[];fullstepbins=[];sectorrows=[]
    for i,(m,of) in enumerate(zip(metadata,offsets)):
        pos,length=map(int,of);v=flat[pos:pos+length]
        if digest(v)!=m['current_sha256']:raise ValueError('Current polygon hash mismatch')
        pred,stats=reconstruct_next(v,True);tr=np.asarray(m['translation'][:2])/m['translation'][2]
        lower=PolygonSupport(v,m['theta'],tr,m['R']);upper=PolygonSupport(pred,m['theta'],tr,m['R'])
        cl=lower.coefficients((0,1,3));cp=upper.coefficients((0,1,3));ql=cl[0].real;qp=cp[0].real;q=m['q_absolute']
        if not qp<q<ql:raise ValueError('Target not inside predicted perimeter bracket')
        w=(q-ql)/(qp-ql);cq={k:(1-w)*cl[k]+w*cp[k] for k in (0,1,3)};a3=cq[3]/q
        rr=ref.at_q(q);refa3=complex(rr['a3']);reference_bins=ref.weighted_integrals(rr['coefficients'],intervals)
        I0=lower.weighted_integrals(intervals);I1=upper.weighted_integrals(intervals);Iq=((1-w)*I0+w*I1)/q
        residual=Iq-reference_bins;centered=residual-dipole_weighted_integrals(cq[1]/q,intervals)
        targetchange=Iq-I0/ql;fullchange=I1/qp-I0/ql
        a3lower=cl[3]/ql;a3upper=cp[3]/qp;deltaformula=(cp[3]-cl[3]-a3lower*(qp-ql))/qp
        assert abs(deltaformula-(a3upper-a3lower))<3e-16
        assert abs(math.fsum(residual)-(a3-refa3).real)<3e-15
        assert abs(math.fsum(centered)-(a3-refa3).real)<3e-15
        assert abs(math.fsum(targetchange)-(a3-a3lower).real)<3e-15
        assert np.max(np.abs(targetchange-w*qp/q*fullchange))<1e-15
        record={k:v for k,v in m.items() if k!='translation'}
        record.update(translation_x_numerator=m['translation'][0],translation_y_numerator=m['translation'][1],translation_denominator=m['translation'][2],predicted_upper_sha256=digest(pred),**stats,c0_lower=ql,predicted_c0_upper=qp,predicted_minkowski_weight=w,
          lower_c1_real=cl[1].real,lower_c1_imag=cl[1].imag,predicted_upper_c1_real=cp[1].real,predicted_upper_c1_imag=cp[1].imag,predicted_target_c1_real=cq[1].real,predicted_target_c1_imag=cq[1].imag,
          predicted_upper_c3_real=cp[3].real,predicted_upper_c3_imag=cp[3].imag,predicted_a3_real=a3.real,predicted_a3_imag=a3.imag,reference_a3_real=refa3.real,reference_a3_imag=refa3.imag,predicted_residual_real=(a3-refa3).real,predicted_residual_imag=(a3-refa3).imag,
          predicted_delta_c0=qp-ql,predicted_delta_c3_real=(cp[3]-cl[3]).real,predicted_delta_c3_imag=(cp[3]-cl[3]).imag,predicted_full_step_delta_a3_real=(a3upper-a3lower).real,predicted_full_step_delta_a3_imag=(a3upper-a3lower).imag,predicted_target_step_delta_a3_real=(a3-a3lower).real,predicted_target_step_delta_a3_imag=(a3-a3lower).imag,
          residual_angular_closure_error=abs(math.fsum(residual)-(a3-refa3).real),centered_residual_angular_closure_error=abs(math.fsum(centered)-(a3-refa3).real),target_step_angular_closure_error=abs(math.fsum(targetchange)-(a3-a3lower).real))
        rows.append(record);resbins.append(residual);centeredbins.append(centered);stepbins.append(targetchange);fullstepbins.append(fullchange)
        S0=lower.weighted_integrals(sectors);S1=upper.weighted_integrals(sectors);Sq=((1-w)*S0+w*S1)/q;Sr=ref.weighted_integrals(rr['coefficients'],sectors)
        Sd=dipole_weighted_integrals(cq[1]/q,sectors)
        for j,center in enumerate((0,180,-1)):
            sectorrows.append({'case_index':i,'row':m['row'],'R':m['R'],'orientation':('theta0','atan_1_8','pi_over_4')[m['trajectory']%3],'q_historical':m['q_historical'],'q_absolute':q,'shape_center_deg':center,'halfwidth_deg':5 if center!=-1 else 180,'predicted_residual_uncentered':float(Sq[j]-Sr[j]),'predicted_residual_centered':float(Sq[j]-Sd[j]-Sr[j]),'predicted_target_step_contribution_uncentered':float(Sq[j]-S0[j]/ql),'predicted_target_step_contribution_centered':float(Sq[j]-Sd[j]-S0[j]/ql+dipole_weighted_integrals(cl[1]/ql,sectors)[j])})
        if i%256==0:print('Forecast',i,'/',len(metadata),flush=True)
    # Separate prediction artifacts exist before the first benchmark access.
    predpath=output/'endpoint_predictions.csv';write_csv(predpath,rows)
    angularpath=output/'angular_predictions.npz'
    np.savez_compressed(angularpath,residual_bins_uncentered=np.asarray(resbins),residual_bins_centered=np.asarray(centeredbins),target_step_change_bins=np.asarray(stepbins),full_step_change_bins=np.asarray(fullstepbins))
    seal={'predictions_sha256':sha(predpath),'angular_predictions_sha256':sha(angularpath),'inputs_sha256':sha(inputs),'reference_sha256':sha(reference),'predictor_script_sha256':sha(Path(__file__)),'geometry_script_sha256':sha(Path(__file__).with_name('one_step_prediction.py'))}
    (output/'FORECAST_SEAL.json').write_text(json.dumps(seal,indent=2)+'\n')
    # The benchmark is used only here, after forecasts and their seal are saved.
    validation=[]
    with np.load(benchmark,allow_pickle=False) as z:
        for i,r in enumerate(rows):
            a3=complex(r['predicted_a3_real'],r['predicted_a3_imag']);c3=complex(r['predicted_upper_c3_real'],r['predicted_upper_c3_imag'])
            validation.append({'case_index':i,'integer_geometry_equal':r['predicted_upper_sha256']==str(z['upper_sha256'][i]),'remaining_points_equal':r['surviving_integer_points']==int(z['upper_remaining_points'][i]),'a3_absolute_error':abs(a3-z['a3_target'][i]),'upper_c3_absolute_error':abs(c3-z['upper_c3'][i]),'weight_absolute_error':abs(r['predicted_minkowski_weight']-float(z['weight'][i]))})
    write_csv(output/'endpoint_validation.csv',validation);write_csv(output/'predicted_sector_contributions.csv',sectorrows)
    groups=[];binmeans=[]
    for R in (256,512):
        for j,orientation in enumerate(('theta0','atan_1_8','pi_over_4')):
            for qh in (.95,.9,.85):
                ix=[i for i,r in enumerate(rows) if r['R']==R and r['trajectory']%3==j and r['q_historical']==qh];assert len(ix)==64
                d={'R':R,'orientation':orientation,'q_historical':qh,'count':len(ix),'mean_predicted_residual_real':math.fsum(rows[i]['predicted_residual_real'] for i in ix)/64,'mean_predicted_target_step_delta_a3_real':math.fsum(rows[i]['predicted_target_step_delta_a3_real'] for i in ix)/64,'min_predicted_residual_real':min(rows[i]['predicted_residual_real'] for i in ix),'max_predicted_residual_real':max(rows[i]['predicted_residual_real'] for i in ix),'max_a3_absolute_error':max(validation[i]['a3_absolute_error'] for i in ix)}
                groups.append(d)
                for b in range(360):
                    binmeans.append({'R':R,'orientation':orientation,'q_historical':qh,'phi_left_deg':b,'phi_right_deg':b+1,'mean_residual_uncentered':math.fsum(resbins[i][b] for i in ix)/64,'mean_residual_centered':math.fsum(centeredbins[i][b] for i in ix)/64,'mean_target_step_change':math.fsum(stepbins[i][b] for i in ix)/64})
    write_csv(output/'endpoint_groups.csv',groups);write_csv(output/'angular_prediction_groups.csv',binmeans)
    result={'status':'PASS_CONDITIONAL_ENDPOINT_MODE3_PREDICTION','count':len(rows),'qualification':qtest,'integer_geometry_failures':sum(not r['integer_geometry_equal'] for r in validation),'remaining_points_failures':sum(not r['remaining_points_equal'] for r in validation),'max_a3_absolute_error':max(r['a3_absolute_error'] for r in validation),'max_upper_c3_absolute_error':max(r['upper_c3_absolute_error'] for r in validation),'max_weight_absolute_error':max(r['weight_absolute_error'] for r in validation),'max_angular_closure_error':max(r['residual_angular_closure_error'] for r in rows),'max_centered_angular_closure_error':max(r['centered_residual_angular_closure_error'] for r in rows),'max_target_step_angular_closure_error':max(r['target_step_angular_closure_error'] for r in rows),'fitted_parameters':0,'new_independent_trajectory':False,'new_pde_solve':False,'reference_limit':'Saved affine reference shared with prior localization; numerical agreement does not certify PDE accuracy.','scope':'Conditional exact one-step microscopic geometry from full lower integer polygons. Target lower layers selected retrospectively; does not predict initial-to-target dynamics or crossing time. q_absolute preserves historical contraction fraction of each initial polygon mean support.','seal':seal,'wall_seconds':time.perf_counter()-started}
    assert result['integer_geometry_failures']==result['remaining_points_failures']==0
    (output/'endpoint_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);here=Path(__file__).resolve().parent
    p.add_argument('--inputs',type=Path,default=here/'inputs/endpoint_jump_inputs.npz');p.add_argument('--benchmark',type=Path,default=here/'inputs/endpoint_jump_benchmark.npz');p.add_argument('--reference',type=Path,default=here/'inputs/reference_contraction.npz');p.add_argument('--output',type=Path,default=here/'results_replay');a=p.parse_args()
    run(a.inputs,a.benchmark,a.reference,a.output)

if __name__=='__main__':main()
