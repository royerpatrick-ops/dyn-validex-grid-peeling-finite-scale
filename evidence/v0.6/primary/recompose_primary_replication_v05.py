"""Independent read-only primary audit of completed v0.4/v0.5 campaign.
No frozen analyzer/oracle is imported or re-executed. No geometry collection.
Reconstructs all endpoint c3 from supplied vertices, verifies targets/reference
joins and sensitivities, and recomputes the frozen primary Student rule.
"""
from pathlib import Path
from collections import defaultdict
import csv,json,math,hashlib
import numpy as np
from scipy.stats import t
ROOT=Path('/workspace/scratch/cbad571848ca')
P=ROOT/'tmp/resultats_replication_v05/_prepare_replication_v0_4_20261027/_prepare_replication_v0_4_20261027'
E=ROOT/'tmp/resultats_replication_v05/_execute_replication_v0_5_20261027/_execute_replication_v0_5_20261027'
HERE=Path(__file__).resolve().parent
C=P/'campaign_v04'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
rows=list(csv.DictReader((C/'campaign_cells_v04.csv').open(newline='',encoding='utf-8')))
reported=json.loads((P/'analysis_v04/analysis_results_v04.json').read_text())
qual=json.loads((P/'pilot/reference_qualification_v04.json').read_text())
with np.load(C/'campaign_reference_v04.npz',allow_pickle=False) as z:ref={k:z[k] for k in z.files}
ORS=['theta0','atan_1_8','pi_over_4'];RAD=[256,512];QS=[.95,.90,.85];LABELS=['0.30','0.31','0.32','1/3','0.34','0.35','0.36'];BETAS=[.30,.31,.32,1/3,.34,.35,.36]
errors=defaultdict(float);groups=defaultdict(dict);identities=defaultdict(set);allkeys=set()
for r in rows:
 u=int(r['row_index']);R=int(r['R']);o=r['orientation'];q=float(r['q_historical']);label=r['beta_label'];b=float(r['beta'])
 key=(u,R,o,q,label)
 if key in allkeys:raise ValueError('duplicate full scientific key')
 allkeys.add(key)
 vals=[float(r[k]) for k in ['q_absolute','beta','peel_a3_shape_real','peel_a3_shape_imag','reference_a3_shape_real','reference_a3_shape_imag','epsilon_interp_a3','epsilon_ref_a3_qualification_grid_max']]
 if not np.all(np.isfinite(vals)):raise ValueError('nonfinite scientific field')
 if label not in LABELS or b!=BETAS[LABELS.index(label)]:raise ValueError('beta label/value mismatch')
 identity=(int(r['tx']),int(r['ty']),int(r['denominator']),r['provenance']);identities[u].add(identity)
 p=complex(float(r['peel_a3_shape_real']),float(r['peel_a3_shape_imag']));rref=complex(float(r['reference_a3_shape_real']),float(r['reference_a3_shape_imag']))
 er=float(r['epsilon_ref_a3_qualification_grid_max']);ei=float(r['epsilon_interp_a3']);qa=float(r['q_absolute'])
 if ei<0 or er<0:raise ValueError('negative sensitivity')
 expected=qual['epsilon_ref_qualification_grid_max_final'][f'{b:.17g}']
 errors['CSV_epsilon_vs_qualification']=max(errors['CSV_epsilon_vs_qualification'],abs(er-expected))
 groups[u,R,o,q][label]={'peel':p,'ref':rref,'ei':ei,'er':er,'qa':qa,'loss':abs(p-rref)**2}
expectedkeys={(u,R,o,q,l) for u in range(64) for R in RAD for o in ORS for q in QS for l in LABELS}
if allkeys!=expectedkeys:raise ValueError('missing or extra scientific grid cells')
if len(rows)!=8064 or len(groups)!=1152:raise ValueError('row/group count mismatch')
if set(identities)!=set(range(64)) or any(len(v)!=1 for v in identities.values()):raise ValueError('inconsistent row identity')
if len(set(next(iter(v)) for v in identities.values()))!=64:raise ValueError('duplicate translation identity')
# Anchor all CSV references and exact target order to production oracle NPZ.
if ref['beta_grid'].tolist()!=BETAS or ref['modes'].tolist()!=[0,3,6,9,12]:raise ValueError('oracle array order mismatch')
errors['oracle_c0_target']=float(np.max(abs(ref['coefficients'][:,:,0]-ref['q_absolute'][None,:])))
errors['oracle_normalization']=float(np.max(abs(ref['coefficients']/ref['coefficients'][:,:,[0]].real-ref['normalized_coefficients'])))
baseorder=[(u,R,o,q) for u in range(64) for R in RAD for o in ORS for q in QS]
for j,key in enumerate(baseorder):
 g=groups[key]
 for bi,l in enumerate(LABELS):
  errors['CSV_reference_vs_oracleNPZ']=max(errors['CSV_reference_vs_oracleNPZ'],abs(g[l]['ref']-ref['normalized_coefficients'][bi,j,1]))
  errors['CSV_q_vs_oracleNPZ']=max(errors['CSV_q_vs_oracleNPZ'],abs(g[l]['qa']-ref['q_absolute'][j]))
  errors['CSV_epsilon_vs_oracleNPZ']=max(errors['CSV_epsilon_vs_oracleNPZ'],abs(g[l]['er']-ref['epsilon_ref_qualification_grid_max'][bi]))
  errors['common_peel_across_beta']=max(errors['common_peel_across_beta'],abs(g[l]['peel']-g['0.32']['peel']))
  errors['common_interpolation_epsilon_across_beta']=max(errors['common_interpolation_epsilon_across_beta'],abs(g[l]['ei']-g['0.32']['ei']))
# Geometry-to-observable check at all 1152 saved targets; no peeling rerun.
def edge(vertices,R,m):
 v=np.asarray(vertices,dtype=float);ev=np.roll(v,-1,axis=0)-v;ll=np.linalg.norm(ev,axis=1)
 phase=(ev[:,1]+1j*ev[:,0])/ll;s=(ll/R)*phase**m
 return complex(math.fsum(s.real),math.fsum(s.imag))/(2*math.pi*(1-m*m))
raw_count=0
for u in range(64):
 rawpath=C/'rows'/f'row_{u:02d}_main.npz'
 with np.load(rawpath,allow_pickle=False) as z:
  raw={k:z[k] for k in ['translations','radii','orientation_pq','orientation_theta','initial_perimeter','target_q_historical','target_q_absolute','endpoint_layers','endpoint_perimeters','endpoint_vertex_offsets','endpoint_vertices_ordered_flat','fourier_modes','endpoint_coefficients_raw','weights_after','a_shape_exp_plus_i_m_theta','epsilon_interpolation_a3']}
 ids=next(iter(identities[u]))
 if not np.all(raw['translations']==np.array(ids[:3])[None,:]):raise ValueError('raw translation mismatch')
 if raw['radii'].tolist()!=[256,256,256,512,512,512]:raise ValueError('raw radius order mismatch')
 if not np.array_equal(raw['orientation_pq'],np.array([[1,0],[8,1],[1,1]]*2)):raise ValueError('raw orientation order mismatch')
 if raw['target_q_historical'].tolist()!=QS or raw['fourier_modes'].tolist()!=[0,1,3,6,9,12]:raise ValueError('raw q/mode order mismatch')
 if not np.all(raw['endpoint_layers'][:,:,1]-raw['endpoint_layers'][:,:,0]==1):raise ValueError('nonadjacent endpoints')
 for ri,R in enumerate(RAD):
  for oi,o in enumerate(ORS):
   ti=ri*3+oi;theta=float(raw['orientation_theta'][ti])
   for qi,q in enumerate(QS):
    key=(u,R,o,q);g=groups[key]['0.32'];qa=q*float(raw['initial_perimeter'][ti])/(2*math.pi*R)
    w=float(raw['weights_after'][ti,qi]);pb,pa=map(float,raw['endpoint_perimeters'][ti,qi]);pw=(pb-q*float(raw['initial_perimeter'][ti]))/(pb-pa)
    errors['raw_weight_formula']=max(errors['raw_weight_formula'],abs(w-pw))
    errors['raw_q_formula']=max(errors['raw_q_formula'],abs(qa-g['qa']),abs(qa-float(raw['target_q_absolute'][ti,qi])))
    vals=[]
    for ep in range(2):
     start,count=map(int,raw['endpoint_vertex_offsets'][ti,qi,ep]);vertices=raw['endpoint_vertices_ordered_flat'][start:start+count]
     c0=edge(vertices,R,0);c3=edge(vertices,R,3);vals.append(c3)
     errors['raw_endpoint_c0_from_vertices']=max(errors['raw_endpoint_c0_from_vertices'],abs(c0-raw['endpoint_coefficients_raw'][ti,qi,ep,0]))
     errors['raw_endpoint_c3_from_vertices']=max(errors['raw_endpoint_c3_from_vertices'],abs(c3-raw['endpoint_coefficients_raw'][ti,qi,ep,2]))
    peel=((1-w)*vals[0]+w*vals[1])/qa*np.exp(3j*theta)
    eps=2*w*(1-w)*(pb-pa)/(2*math.pi*R*qa)
    errors['raw_vertices_to_CSV_peel']=max(errors['raw_vertices_to_CSV_peel'],abs(peel-g['peel']))
    errors['raw_CSV_peel_stored_array']=max(errors['raw_CSV_peel_stored_array'],abs(g['peel']-raw['a_shape_exp_plus_i_m_theta'][ti,qi,2]))
    errors['raw_epsilon_formula_to_CSV']=max(errors['raw_epsilon_formula_to_CSV'],abs(eps-g['ei']))
    raw_count+=1
# Recompose primary and every supplied contrast field.
cells={};byrow=defaultdict(list)
for key,g in groups.items():
 a=g['1/3'];c=g['0.32'];p=a['peel'];ra=a['ref'];rc=c['ref'];ei=a['ei'];ea=a['er'];ec=c['er']
 d=abs(p-ra)**2-abs(p-rc)**2
 parts=[2*ei*abs(ra-rc),2*abs(p-ra)*ea+ea**2,2*abs(p-rc)*ec+ec**2,2*ei*(ea+ec)]
 B=sum(parts);cells[key]={'d':d,'B':B,'Wcell':d-B,'parts':parts}
 byrow[key[0],key[1]].append(key)
contrastcsv=list(csv.DictReader((P/'analysis_v04/cell_contrasts_v04.csv').open(newline='',encoding='utf-8')))
for r in contrastcsv:
 key=(int(r['row_index']),int(r['R']),r['orientation'],float(r['q_historical']));v=cells[key]
 pub=[float(r[k]) for k in ['d_raw','B_empirical','d_minus_B','B_common_interpolation','B_reference_1_over_3','B_reference_0_32','B_mixed']]
 obs=[v['d'],v['B'],v['Wcell'],*v['parts']]
 errors['all_cell_contrasts_and_B_components']=max(errors['all_cell_contrasts_and_B_components'],float(np.max(abs(np.array(pub)-obs))))
w=np.array([math.fsum(cells[key]['Wcell'] for key in byrow[u,512])/9 for u in range(64)])
D=np.array([[math.fsum(cells[k]['d'] for k in byrow[u,R])/9 for R in RAD] for u in range(64)])
BB=np.array([[math.fsum(cells[k]['B'] for k in byrow[u,R])/9 for R in RAD] for u in range(64)])
mean=float(w.mean());sd=float(w.std(ddof=1));crit=float(t.ppf(.975,63));lcb=mean-crit*sd/8
primary={'n_complete_translation_rows':64,'cells_per_row_at_primary_R':9,'mean_W':mean,'sd_W_ddof1':sd,'student_t_0_975_df63':crit,'LCB_nominal':lcb,
 'decision_recomputed':'FINITE_R512_PREFERENCE_0_32' if lcb>0 else 'NON_DISCRIMINE','W_min':float(w.min()),'W_max':float(w.max()),'W_positive_count':int(np.sum(w>0)),
 'mean_D_R512':float(D[:,1].mean()),'mean_B_R512':float(BB[:,1].mean()),'mean_D_over_mean_B_R512':float(D[:,1].mean()/BB[:,1].mean())}
for name,obv in [('mean_W',mean),('sd_W_ddof1',sd),('student_t_quantile_0_975_df63',crit),('lower_confidence_bound_nominal',lcb)]:errors['primary_statistics_vs_receipt']=max(errors['primary_statistics_vs_receipt'],abs(obv-reported['primary'][name]))
rowreceipt={r['row_index']:r for r in reported['primary']['row_scores_with_exact_identity']}
for u in range(64):errors['all_64_W_vs_receipt']=max(errors['all_64_W_vs_receipt'],abs(w[u]-rowreceipt[u]['W']))
loo=[]
for u in range(64):
 v=np.delete(w,u);mm=float(v.mean());ll=float(mm-t.ppf(.975,62)*v.std(ddof=1)/math.sqrt(63));loo.append({'omitted_row':u,'mean_W':mm,'LCB':ll})
 rr=reported['leave_one_translation_out']['rows'][u]
 errors['all_64_LOO_statistics_vs_receipt']=max(errors['all_64_LOO_statistics_vs_receipt'],abs(mm-rr['mean_W']),abs(ll-rr['lcb']))
def desc(v):
 v=np.array(v);return {'n':len(v),'mean':float(v.mean()),'sd_ddof1':float(v.std(ddof=1)),'min':float(v.min()),'max':float(v.max())}
secondary={'D512_minus_D256':desc(D[:,1]-D[:,0]),'D_R512_by_orientation':{},'D_R512_by_q':{},'dedicated_pi4_q095':{},'loss_means_by_R_beta':{}}
for o in ORS:secondary['D_R512_by_orientation'][o]=desc([math.fsum(cells[u,512,o,q]['d'] for q in QS)/3 for u in range(64)])
for q in QS:secondary['D_R512_by_q'][str(q)]=desc([math.fsum(cells[u,512,o,q]['d'] for o in ORS)/3 for u in range(64)])
for R in RAD:
 secondary['dedicated_pi4_q095'][str(R)]=desc([cells[u,R,'pi_over_4',.95]['d'] for u in range(64)])
 secondary['loss_means_by_R_beta'][str(R)]={l:float(np.mean([groups[u,R,o,q][l]['loss'] for u in range(64) for o in ORS for q in QS])) for l in LABELS}
secondary['affine_mean_loss_R512_over_R256']=secondary['loss_means_by_R_beta']['512']['1/3']/secondary['loss_means_by_R_beta']['256']['1/3']
# Verify the predeclared descriptive summaries independently too.
sec_reported=reported['predeclared_secondary']
sec_pairs=[(secondary['D512_minus_D256'],sec_reported['raw_paired_D512_minus_D256_by_translation'])]
sec_pairs += [(v,sec_reported['raw_D_at_R512_by_orientation'][k]) for k,v in secondary['D_R512_by_orientation'].items()]
sec_pairs += [(v,sec_reported['raw_D_at_R512_by_q'][f'{float(k):.2f}']) for k,v in secondary['D_R512_by_q'].items()]
sec_pairs += [(v,sec_reported['raw_D_pi_over_4_q_0_95_by_R'][k]) for k,v in secondary['dedicated_pi4_q095'].items()]
errors['secondary_descriptive_statistics_vs_receipt']=max(abs(a[k]-b[k]) for a,b in sec_pairs for k in ['mean','sd_ddof1','min','max'])
res={'status':'INDEPENDENT_PRIMARY_RECOMPUTED_POSITIVE' if lcb>0 else 'INDEPENDENT_PRIMARY_RECOMPUTED_NONPOSITIVE',
 'no_frozen_analyzer_or_oracle_imported':True,'no_analysis_or_collection_rerun':True,
 'counts':{'CSV_beta_cells':len(rows),'base_groups':len(groups),'complete_translation_rows':len(identities),'raw_npz_targets_reconstructed_from_vertices':raw_count,'contrast_csv_rows':len(contrastcsv)},
 'maximum_absolute_errors':dict(errors),'primary':primary,'leave_one_translation_out':{'minimum_LCB':min(r['LCB'] for r in loo),'maximum_LCB':max(r['LCB'] for r in loo),'positive_count':sum(r['LCB']>0 for r in loo),'rows':loo},
 'row_scores':[{'row_index':u,'identity':list(next(iter(identities[u]))),'W':float(w[u])} for u in range(64)],
 'secondary_descriptive_only':secondary,
 'limits':['Primary recomputation validates the exact frozen finite-comparison rule, with approximate nominal Student coverage.','B is an empirical sensitivity proxy, not an interval-certified enclosure.','No uniform preference, exact beta identification or asymptotic implication.','This audit checks saved endpoint-to-coefficient consistency, not every layer or the new-seed provenance chronology.'],
 'source_hashes':{str(p.relative_to(P)):sha(p) for p in [C/'campaign_cells_v04.csv',C/'campaign_reference_v04.npz',P/'pilot/reference_qualification_v04.json',P/'analysis_v04/analysis_results_v04.json']}}
(HERE/'bilan_primary_replication_v05.json').write_text(json.dumps(res,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'counts':res['counts'],'errors':dict(errors),'primary':primary,'LOO_minimum':res['leave_one_translation_out']['minimum_LCB'],'LOO_positive':res['leave_one_translation_out']['positive_count'],'secondary_affine_loss_ratio':secondary['affine_mean_loss_R512_over_R256']},ensure_ascii=False,indent=2))
