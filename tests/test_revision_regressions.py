from pathlib import Path
import importlib.util, sys, numpy as np, pandas as pd
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'scripts'))
from model_validation import inner_selection, fit_selected_model

def load(name):
 spec=importlib.util.spec_from_file_location('stage'+name[:2],root/'scripts'/name);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
rng=np.random.default_rng(4);n=180
base=pd.DataFrame({'virus_strain_matched':np.repeat([f'v{i}' for i in range(12)],15),'serum_strain_matched':np.tile(['s0','s1','s2'],60),'serumIsolate':np.tile(['s0','s1','s2'],60),'virusYear':2000,'serumYear':2000,'temporal_distance':rng.integers(0,10,n),'change_1':rng.integers(0,2,n),'K1':rng.integers(0,2,n),'1K':rng.integers(0,2,n),'sub_1_K_N':rng.integers(0,2,n),'corrected_titer':0.0})
base['standardized_titer']=base['change_1']+0.5*base['K1']+rng.normal(size=n)*0.1
train=np.arange(150);test=np.arange(150,n)
m=load('06_train_model.py');params=m.get_lightgbm_params();params.update(n_estimators=30,min_child_samples=3,num_leaves=5,n_jobs=1)
args=dict(train_idx=train,val_idx=test,params=params,distance_cols=['temporal_distance'])
a,pa=m.run_models_for_split(df_binary=base,df_substitution=base,rng=np.random.default_rng(42),**args)
altered=base.copy();altered.loc[test,'standardized_titer']+=100
b,pb=m.run_models_for_split(df_binary=altered,df_substitution=altered,rng=np.random.default_rng(42),**args)
cols=[c for c in pa if c.startswith('predicted_') or c=='additive_only_pred']
np.testing.assert_allclose(pa[cols],pb[cols])
assert [x['best_iteration'] for x in a]==[x['best_iteration'] for x in b]
seen=[]
def baseline(a,b):
 assert set(a.virus_strain_matched).isdisjoint(b.virus_strain_matched)
 seen.append((len(a),len(b)));return np.zeros(len(a)),np.zeros(len(b))
sel=inner_selection(base.iloc[train],'virus_strain_matched',baseline)
assert len(seen)==1 and sum(seen[0])==len(train)
assert set(sel[0]).isdisjoint(sel[1])
for p in (root/'scripts').glob('*.py'):compile(p.read_text(),str(p),'exec')
print('PASS: outer target perturbation leaves all six family predictions and selected rounds unchanged.')
print('PASS: inner fit and stopping partitions contain disjoint virus groups; all Python stages compile.')

from common import fit_additive_effects, apply_additive_effects
old=fit_additive_effects(base.iloc[train], 'standardized_titer', 'virus_strain_matched', 'serumIsolate')
renamed=base.copy();renamed['virus_strain_matched']=renamed.virus_strain_matched.map(lambda x: 'reversed'+str(20-int(x[1:])))
new=fit_additive_effects(renamed.iloc[train], 'standardized_titer', 'virus_strain_matched', 'serumIsolate')
def predict_offset(frame, fit):
    return apply_additive_effects(frame, 'serumIsolate', 'virus_strain_matched', fit['intercept'], fit['serum_effects'], fit['virus_effects'])
np.testing.assert_allclose(predict_offset(base.iloc[test],old),predict_offset(renamed.iloc[test],new),atol=1e-10)
np.testing.assert_allclose(predict_offset(base.iloc[train],old),old['fitted'],atol=1e-10)
pre=load('25_preprocess_wic_model_titers.py')
assert pre.normalize_passage('Cell|Egg')=='CELL|EGG'
for value in ['EGG','MIXED','UNKNOWN']:
    assert value in {'EGG','MIXED','UNKNOWN'}
print('PASS: composite passage labels are normalized but only exact classes are excluded.')
match=load('03_match_strains.py').build_matcher(['A/NewYork/55/2004'])
assert match('A/NewYork/55/2001',100)[0] is None
assert match('A/New York/55/2004',100)[0]=='A/NewYork/55/2004'
print('PASS: unseen-virus offsets are invariant to arbitrary strain-name ordering; training fits are preserved.')
print('PASS: composite passage labels are normalized but only exact classes are excluded.')
