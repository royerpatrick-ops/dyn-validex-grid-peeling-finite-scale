#!/usr/bin/env python3
"""Recompose scientific secondary summaries of the 64-row replication.

Only delivered CSV/JSON are read. No delivered analysis code is imported, and
no new grid, hull, flow reference or random sample is generated. Additional
intervals/figures below are descriptive audit summaries, not new primaries.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


def describe(values):
    z=np.asarray(values,float)
    n=len(z); mean=float(z.mean()); sd=float(z.std(ddof=1))
    half=float(student_t.ppf(.975,n-1)*sd/math.sqrt(n))
    return {'n_translation_rows':n,'mean':mean,'sd_ddof1':sd,'min':float(z.min()),
            'max':float(z.max()),'positive_count':int((z>0).sum()),
            'negative_count':int((z<0).sum()),
            'nominal_95_Student_interval_descriptive_unadjusted':[mean-half,mean+half]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preparation-root',type=Path,required=True)
    parser.add_argument('--historical-profile-csv',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args(); root=args.preparation_root.resolve()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    csv=root/'campaign_v04/campaign_cells_v04.csv'
    result_path=root/'analysis_v04/analysis_results_v04.json'
    d=pd.read_csv(csv,float_precision='round_trip');report=json.loads(result_path.read_text())
    assert len(d)==8064 and d.row_index.nunique()==64
    d['peel']=d.peel_a3_shape_real+1j*d.peel_a3_shape_imag
    d['reference']=d.reference_a3_shape_real+1j*d.reference_a3_shape_imag
    d['loss']=abs(d.peel-d.reference)**2
    bkeys=['row_index','R','orientation','q_historical']
    assert not d.duplicated(bkeys+['beta_label']).any()
    a=d[d.beta_label=='1/3'].set_index(bkeys).sort_index()
    b=d[d.beta_label=='0.32'].set_index(bkeys).sort_index()
    assert a.index.equals(b.index) and len(a)==1152
    assert (a.peel==b.peel).all() and (a.epsilon_interp_a3==b.epsilon_interp_a3).all()
    eps=a.epsilon_interp_a3.to_numpy(float)
    ea=a.epsilon_ref_a3_qualification_grid_max.to_numpy(float)
    eb=b.epsilon_ref_a3_qualification_grid_max.to_numpy(float)
    common=2*eps*abs(a.reference.to_numpy()-b.reference.to_numpy())
    refa=2*abs(a.peel.to_numpy()-a.reference.to_numpy())*ea+ea**2
    refb=2*abs(b.peel.to_numpy()-b.reference.to_numpy())*eb+eb**2
    mixed=2*eps*(ea+eb)
    contrast=pd.DataFrame({'D_raw':a.loss.to_numpy()-b.loss.to_numpy(),
                           'B_interp':common,'B_ref_affine':refa,'B_ref_032':refb,
                           'B_mixed':mixed,'B_empirical':common+refa+refb+mixed},index=a.index)
    contrast['D_minus_B']=contrast.D_raw-contrast.B_empirical
    contrast=contrast.reset_index()
    loss_by_row=d.groupby(['row_index','R','beta_label'],sort=True).loss.mean().reset_index()
    loss_summary=[]
    order=['0.30','0.31','0.32','1/3','0.34','0.35','0.36']
    max_report_error=0.
    for radius in [256,512]:
        for label in order:
            s=loss_by_row[(loss_by_row.R==radius)&(loss_by_row.beta_label==label)]
            row={'R':radius,'beta_label':label,**describe(s.loss)}
            loss_summary.append(row)
            # Reconcile raw cell-level means with every reported angle/q mean.
            for (orientation,q),g in d[(d.R==radius)&(d.beta_label==label)].groupby(['orientation','q_historical']):
                stored=next(r for r in report['predeclared_secondary']['all_seven_beta_losses']
                            if r['R']==radius and r['beta_label']==label and r['orientation']==orientation and r['q_historical']==q)
                max_report_error=max(max_report_error,abs(float(g.loss.mean())-stored['mean']))
    cells=[]
    for (radius,orientation,q),g in contrast.groupby(['R','orientation','q_historical'],sort=True):
        cells.append({'R':int(radius),'orientation':orientation,'q_historical':float(q),
                      'D_raw':describe(g.D_raw),'D_minus_B':describe(g.D_minus_B),
                      'B_empirical_mean':float(g.B_empirical.mean()),
                      'B_interp_mean':float(g.B_interp.mean()),
                      'B_reference_and_mixed_mean':float((g.B_ref_affine+g.B_ref_032+g.B_mixed).mean())})
    by_row_radius=contrast.groupby(['row_index','R'])[['D_raw','B_empirical','D_minus_B']].mean()
    rawradius=by_row_radius.D_raw.unstack('R')
    adjustedradius=by_row_radius.D_minus_B.unstack('R')
    interaction={'D512_minus_D256_raw':describe(rawradius[512]-rawradius[256]),
                 'adjusted512_minus_adjusted256':describe(adjustedradius[512]-adjustedradius[256])}
    max_report_error=max(max_report_error,abs(interaction['D512_minus_D256_raw']['mean']-
                     report['predeclared_secondary']['raw_paired_D512_minus_D256_by_translation']['mean']))
    # Primary is recomposed only to contextualize the explicitly prespecified
    # finite R512 average; it is not selected from these secondary findings.
    primary=describe(adjustedradius[512])
    lcb=primary['nominal_95_Student_interval_descriptive_unadjusted'][0]
    max_report_error=max(max_report_error,abs(primary['mean']-report['primary']['mean_W']),
                         abs(lcb-report['primary']['lower_confidence_bound_nominal']))
    assert max_report_error<5e-20

    old=pd.read_csv(args.historical_profile_csv,float_precision='round_trip')
    old=old[old.R.isin([256,512])]
    # Historical rows contain mean loss over the three q values for each
    # orientation, so the equal-angle mean gives the same nine-cell average.
    old['beta_label']=old.beta.map(lambda beta:'1/3' if abs(beta-1/3)<1e-12 else f'{beta:.2f}')
    oldrows=old.groupby(['translation_key','R','beta_label']).loss_mean_q.mean().reset_index()
    old_new=[]
    for entry in loss_summary:
        z=oldrows[(oldrows.R==entry['R'])&(oldrows.beta_label==entry['beta_label'])].loss_mean_q
        olddesc=describe(z)
        old_new.append({'R':entry['R'],'beta_label':entry['beta_label'],
                        'historical_n':len(z),'new_n':64,
                        'historical_mean_loss':olddesc['mean'],'new_mean_loss':entry['mean'],
                        'new_minus_old_mean':entry['mean']-olddesc['mean'],
                        'new_relative_to_old_minus_one':entry['mean']/olddesc['mean']-1})

    affine256=next(r['mean'] for r in loss_summary if r['R']==256 and r['beta_label']=='1/3')
    affine512=next(r['mean'] for r in loss_summary if r['R']==512 and r['beta_label']=='1/3')
    avg032256=next(r['mean'] for r in loss_summary if r['R']==256 and r['beta_label']=='0.32')
    avg032512=next(r['mean'] for r in loss_summary if r['R']==512 and r['beta_label']=='0.32')
    result={'status':'SECONDARY_SCIENTIFIC_RECOMPOSITION_OK',
            'input_sha256':{str(csv):hashlib.sha256(csv.read_bytes()).hexdigest(),
                            str(result_path):hashlib.sha256(result_path.read_bytes()).hexdigest(),
                            str(args.historical_profile_csv.resolve()):hashlib.sha256(args.historical_profile_csv.read_bytes()).hexdigest()},
            'counts':{'beta_expanded_cells':len(d),'base_cells':len(contrast),'independent_new_translation_rows':64},
            'maximum_mean_disagreement_with_report':max_report_error,
            'loss_summary_equal_nine_cells_per_translation':loss_summary,
            'cell_secondary_contrasts':cells,
            'paired_radius_interaction':interaction,
            'primary_context_R512':{**primary,'nominal_lower_bound':lcb,
                                   'interpretation':'The single prespecified R512 equal-nine-cell mean after empirical sensitivity; no uniform cellwise claim.'},
            'finite_radius_loss_changes':{'affine_loss_R256':affine256,'affine_loss_R512':affine512,
                                          'affine_loss_R512_over_R256':affine512/affine256,
                                          'beta032_loss_R256':avg032256,'beta032_loss_R512':avg032512,
                                          'beta032_loss_R512_over_R256':avg032512/avg032256},
            'historical16_vs_new64_unpaired_descriptive':old_new,
            'scientific_scope':[
                'Confirmed target: prespecified finite-scale average squared normalized mode3 mismatch, beta0.32 versus beta1/3, at R512, after empirical sensitivity B.',
                'The data use one nonelliptic C3 shape, three fixed lattice/shape orientations, three contractions, and a dyadic translation law. They do not establish all-shape or all-orientation universality.',
                'At R512 eight of nine raw mean cell contrasts favor0.32; the prespecified pi/4,q0.95 cell favors1/3 for all64 translations. Thus uniform preference is false in the observed finite cell grid.',
                'After subtracting B, six of nine R512 mean cell contrasts remain positive; all three pi/4 cell means are nonpositive. A negative D-B alone does not establish a robust preference for1/3.',
                'The affine mean loss increases from R256 to R512. These two scales neither demonstrate monotone convergence nor refute an affine asymptotic limit.',
                'The seven-beta minimum changes from1/3 at R256 to0.32 at R512; this identifies a discrete finite-grid ranking, not an exact or continuous beta optimum.',
                'New64 translations reproduce the qualitative historical16 ranking. Historical and new observations are distinct samples and must not be paired across samples or pooled posthoc into a confirmatory sample of80.',
                'Secondary intervals produced by this audit are nominal descriptive and unadjusted; they do not replace the one frozen primary decision.'
            ]}
    (out/'secondary_recomposition_v05.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    pd.DataFrame([{k:v for k,v in row.items() if k in ['R','orientation','q_historical','B_empirical_mean']}|
                  {'D_raw_mean':row['D_raw']['mean'],'D_minus_B_mean':row['D_minus_B']['mean'],
                   'n_positive_D_raw':row['D_raw']['positive_count']} for row in cells]).to_csv(out/'secondary_cells_v05.csv',index=False)
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        orientations=['theta0','atan_1_8','pi_over_4'];qs=[.95,.90,.85]
        scale=1e-7;vmax=max(abs(row[metric]['mean']) for row in cells for metric in ['D_raw','D_minus_B'])/scale
        fig,axes=plt.subplots(2,2,figsize=(10.2,6.9),layout='constrained')
        for col,radius in enumerate([256,512]):
            for ri,metric in enumerate(['D_raw','D_minus_B']):
                arr=np.array([[next(row[metric]['mean'] for row in cells if row['R']==radius and row['orientation']==o and row['q_historical']==q)/scale for q in qs] for o in orientations])
                ax=axes[ri,col];im=ax.imshow(arr,cmap='RdBu',vmin=-vmax,vmax=vmax,aspect='auto')
                for y in range(3):
                    for x in range(3):ax.text(x,y,f'{arr[y,x]:+.3f}',ha='center',va='center',color='white' if abs(arr[y,x])>vmax*.57 else 'black',fontsize=11)
                ax.set_xticks(range(3),['0,95','0,90','0,85']);ax.set_xlabel('Contraction q (convention historique)')
                ax.set_yticks(range(3),['0','atan(1/8)','pi/4']);ax.set_ylabel('Orientation')
                ax.set_title(f'R={radius}: '+('D' if metric=='D_raw' else 'D - B'))
        fig.colorbar(im,ax=axes.ravel().tolist(),label='Contraste moyen / 1e-7 ; positif : avantage à 0,32')
        fig.suptitle('64 nouvelles translations : hétérogénéité à échelle finie\nD = perte(1/3) - perte(0,32) ; B : sensibilité empirique',fontsize=13)
        fig.savefig(out/'secondary_contrasts_v05.png',dpi=190)
        plt.close(fig)
    except ImportError:
        pass
    print(json.dumps({'status':result['status'],'max_report_error':max_report_error,
                      'primary_mean':primary['mean'],'primary_lcb':lcb,
                      'affine_loss_R512_over_R256':affine512/affine256,
                      'radius_interaction_mean':interaction['D512_minus_D256_raw']['mean']},ensure_ascii=False))


if __name__=='__main__':main()
