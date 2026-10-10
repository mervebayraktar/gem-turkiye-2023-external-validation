"""Auditable GEM Turkey validation. No source-column or taxonomy inference.

Run inventory immediately; run analyse only after reviewing an explicit adapter.
All internal accelerations are g. Synthetic tests are never validation results.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import beta as beta_distribution, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
RAW = ROOT / 'data/raw'
TABLES = ROOT / 'outputs/tables'
FIGURES = ROOT / 'outputs/figures'
SHA = '5974372ac3f4a99f25d0649eb030fbe596f23b36'
DLR = np.array([0., .05, .15, .60, 1.])
IM_FIELDS = {'PGA': 'PGA', 'SA(0.3)': 'PSA03', 'SA(0.6)': 'PSA06', 'SA(1.0)': 'PSA10'}


def dump(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')


def csv(frame, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding='utf-8-sig')


def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_model(path):
    root = ET.parse(path).getroot()
    model = root.find('{*}vulnerabilityModel')
    if model is None:
        raise ValueError('NRML vulnerabilityModel bulunamadı')
    functions = {}
    for f in model.findall('{*}vulnerabilityFunction'):
        imls = f.find('{*}imls')
        x = np.fromstring(imls.text, sep=' ')
        mean = np.fromstring(f.find('{*}meanLRs').text, sep=' ')
        cov = np.fromstring(f.find('{*}covLRs').text, sep=' ')
        if not (len(x) == len(mean) == len(cov) and len(x) >= 2):
            raise ValueError('XML dizi boyutları tutarsız')
        if not (np.isfinite(np.r_[x, mean, cov]).all() and np.all(np.diff(x) > 0)
                and np.all((mean >= 0) & (mean <= 1)) and np.all(cov >= 0)):
            raise ValueError('Geçersiz XML sayıları')
        if f.get('id') in functions:
            raise ValueError('Tekrarlı taksonomi')
        functions[f.get('id')] = dict(imt=imls.get('imt'), dist=f.get('dist'), x=x, mean=mean, cov=cov)
    return functions


def evaluate(f, im, outside='nan'):
    """Linear interpolation in IM and mean/CoV; fail closed outside XML domain.

    clamp is an explicit sensitivity option, NOT the default/OpenQuake claim.
    """
    im = np.asarray(im, dtype=float)
    low, high = im < f['x'][0], im > f['x'][-1]
    mean = np.interp(im, f['x'], f['mean'])
    cov = np.interp(im, f['x'], f['cov'])
    if outside == 'nan':
        mean = np.where(low | high, np.nan, mean)
        cov = np.where(low | high, np.nan, cov)
    elif outside != 'clamp':
        raise ValueError('outside must be nan or clamp')
    mean = np.where(im >= 0, mean, np.nan)
    return mean, cov, low | high


def read_grid(path):
    root = ET.parse(path).getroot()
    fields = sorted([e.attrib for e in root.findall('{*}grid_field')], key=lambda x: int(x['index']))
    text = root.find('{*}grid_data').text
    data = np.loadtxt(io.StringIO(text))
    if data.shape[1] != len(fields):
        raise ValueError('ShakeMap kolon sayısı tutarsız')
    df = pd.DataFrame(data, columns=[f['name'] for f in fields])
    if df.duplicated(['LAT', 'LON']).any():
        raise ValueError('Tekrarlı grid düğümü')
    return df, {f['name']: f['units'] for f in fields}, root.attrib


def sample_grid(df, field, lat, lon):
    """Bilinear interpolation; outside coverage stays NaN (no nearest extrapolation)."""
    a = df.pivot(index='LAT', columns='LON', values=field).sort_index().sort_index(axis=1)
    if a.isna().any().any():
        raise ValueError('Grid tam dikdörtgen değil')
    fun = RegularGridInterpolator((a.index.to_numpy(), a.columns.to_numpy()), a.to_numpy(),
                                  bounds_error=False, fill_value=np.nan)
    return fun(np.column_stack([lat, lon]))


def moments_interval(mean, cov):
    """Beta intervals from XML CoV, without adding 0.30 again.

    Invalid beta moments are flagged, not clipped into valid distributions.
    """
    mean, cov = np.broadcast_arrays(np.asarray(mean, float), np.asarray(cov, float))
    var = (mean * cov)**2
    valid = (mean > 0) & (mean < 1) & (var > 0) & (var < mean * (1 - mean))
    lo, hi = np.full(mean.shape, np.nan), np.full(mean.shape, np.nan)
    k = mean[valid] * (1 - mean[valid]) / var[valid] - 1
    stable = valid & (cov > 1e-6)
    ks = mean[stable] * (1 - mean[stable]) / var[stable] - 1
    lo[stable] = beta_distribution.ppf(.025, mean[stable] * ks, (1 - mean[stable]) * ks)
    hi[stable] = beta_distribution.ppf(.975, mean[stable] * ks, (1 - mean[stable]) * ks)
    narrow = valid & ~stable
    lo[narrow] = np.maximum(0, mean[narrow] - 1.959964*np.sqrt(var[narrow]))
    hi[narrow] = np.minimum(1, mean[narrow] + 1.959964*np.sqrt(var[narrow]))
    deterministic = np.isfinite(mean) & (var == 0)
    lo[deterministic] = hi[deterministic] = mean[deterministic]
    return lo, hi, valid | deterministic


def inventory():
    rows = []
    for path in sorted(RAW.glob('gem_vulnerability_*.xml')):
        model = read_model(path)
        for taxonomy, f in model.items():
            parts = taxonomy.split('/')
            rows.append(dict(component=path.stem.replace('gem_vulnerability_', ''), taxonomy=taxonomy,
                             material=parts[0], system=parts[1], design_response=parts[2],
                             height=';'.join(s for s in parts if s.startswith('H:')),
                             occupancy=parts[-1], imt=f['imt'], distribution=f['dist'],
                             n_im=len(f['x']), min_im_g=f['x'].min(), max_im_g=f['x'].max(),
                             mean_monotonic=bool(np.all(np.diff(f['mean']) >= -1e-8)),
                             invalid_beta_nodes=int(((f['mean']*f['cov'])**2 >= f['mean']*(1-f['mean'])).sum())))
    d = pd.DataFrame(rows)
    csv(d, TABLES/'taxonomy_all_components.csv')
    structural = d[d.component == 'structural']
    csv(structural, TABLES/'turkey_structural_taxonomies.csv')
    csv(structural.groupby(['material','system','design_response','imt']).size().reset_index(name='n'),
        TABLES/'taxonomy_summary.csv')
    fields = []
    for ev in ['event1','event2']:
        for kind in ['grid','uncertainty']:
            path = RAW/f'{ev}_{kind}.xml'
            root = ET.parse(path).getroot()
            for el in root.findall('{*}grid_field'):
                fields.append(dict(event=ev, kind=kind, **el.attrib))
    csv(pd.DataFrame(fields), TABLES/'shakemap_fields.csv')
    files = [dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=fingerprint(p))
             for p in sorted(RAW.iterdir()) if p.is_file()]
    dump(TABLES/'raw_checksums.json', files)
    FIGURES.mkdir(parents=True, exist_ok=True)
    counts = structural.imt.value_counts().sort_index()
    ax = counts.plot.bar(color='#226b83', rot=0, figsize=(7,4))
    ax.set(ylabel='Number of taxonomy functions', title='GEM v2026.0.0 Turkey structural model — inventory')
    plt.tight_layout(); plt.savefig(FIGURES/'model_im_inventory.png', dpi=180); plt.close()
    print(f'{len(structural)} structural taxonomies; inventory complete')


def import_buildings(source, adapter, output):
    """User-reviewed mapping from ACTUAL source names to our internal schema."""
    cfg = json.loads(Path(adapter).read_text(encoding='utf-8'))
    if not cfg.get('schema_reviewed') or not cfg.get('research_use_permitted') or not cfg.get('license_evidence'):
        raise ValueError('Şema ve araştırma kullanım koşulları doğrulanmalı')
    df = pd.read_csv(source, **cfg.get('read_csv', {}))
    required = ['building_id','latitude','longitude','damage_label','district','neighborhood',
                'structural_system','year_built','storeys','occupancy']
    mapping = cfg['columns']
    # Optional attributes are explicitly absent, never fabricated or guessed.
    absent_source = [v for v in mapping.values() if v is not None and v not in df.columns]
    if absent_source:
        raise ValueError(f'Kaynak kolonları bulunamadı: {absent_source}')
    result = pd.DataFrame({k: df[v] if v is not None else pd.Series(pd.NA, index=df.index)
                           for k,v in mapping.items()})
    for col in required:
        if col not in result:
            result[col] = pd.NA
    csv(result, output)
    dump(str(output)+'.provenance.json', dict(source_sha256=fingerprint(source), adapter_sha256=fingerprint(adapter),
                                            raw_columns=list(df.columns), license_evidence=cfg['license_evidence']))


def metrics(frame):
    y, p = frame.observed_lr.to_numpy(), frame.predicted_lr.to_numpy()
    ok = np.isfinite(y) & np.isfinite(p)
    y,p = y[ok],p[ok]
    if not len(y):
        return dict(n=0, bias_pred_obs=np.nan, rmse=np.nan, spearman=np.nan)
    rho = float(spearmanr(y,p).statistic) if len(y)>2 and len(np.unique(y))>1 and len(np.unique(p))>1 else np.nan
    return dict(n=len(y), bias_pred_obs=float(p.sum()/y.sum()) if y.sum()>0 else np.nan,
                rmse=float(np.sqrt(np.mean((p-y)**2))), spearman=rho)


def analyse(buildings, assignments, damage_map, out, outside='nan'):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    b = pd.read_csv(buildings, dtype={'building_id':str})
    required = {'building_id','latitude','longitude','damage_label','district','neighborhood',
                'structural_system','year_built','storeys'}
    if not required.issubset(b):
        raise ValueError(f'İç şema eksik: {sorted(required-set(b))}')
    csv(pd.DataFrame({'field':b.columns,'missing_fraction':b.isna().mean().values}), out/'missingness.csv')
    csv(b.damage_label.value_counts(dropna=False).rename_axis('label').reset_index(name='n'), out/'damage_counts.csv')
    if b.building_id.isna().any() or b.building_id.duplicated().any():
        raise ValueError('Bina kimlikleri dolu ve tekil olmalı; tekrar ölçümler önce uzlaştırılmalı')
    for c in ['latitude','longitude','year_built','storeys']:
        b[c] = pd.to_numeric(b[c], errors='raise')
    if not b.latitude.dropna().between(-90,90).all() or not b.longitude.dropna().between(-180,180).all():
        raise ValueError('WGS84 derece koordinatları gerekli')
    a = pd.read_csv(assignments, dtype={'building_id':str})
    needed = {'building_id','scenario','taxonomy','mapping_evidence','reviewed'}
    if not needed.issubset(a) or a.duplicated(['building_id','scenario']).any():
        raise ValueError('Eşleştirme tablosu şeması/tekilliği geçersiz')
    if not a.reviewed.astype(str).str.lower().isin(['true','1']).all() or a.mapping_evidence.isna().any():
        raise ValueError('Her taksonomi ataması kaynak/varsayım ve inceleme kaydı gerektirir')
    if not set(a.building_id).issubset(set(b.building_id)):
        raise ValueError('Atamada bilinmeyen bina kimliği')
    dm = pd.read_csv(damage_map)
    if not {'damage_label','gem_ds','mapping_evidence','reviewed'}.issubset(dm) or dm.damage_label.duplicated().any():
        raise ValueError('Hasar eşleştirmesi geçersiz')
    if not dm.reviewed.astype(str).str.lower().isin(['true','1']).all() or dm.mapping_evidence.isna().any():
        raise ValueError('Hasar tanımları incelenmeli')
    if not dm.gem_ds.dropna().isin(range(5)).all():
        raise ValueError('DS 0..4 veya açıkça dışlanmış boş değer olmalı')
    b = b.merge(dm[['damage_label','gem_ds']], how='left', on='damage_label', validate='many_to_one')
    b['observed_lr'] = b.gem_ds.map(dict(enumerate(DLR)))
    im_rows = []
    for event in ['event1','event2']:
        grid,units,metadata = read_grid(RAW/f'{event}_grid.xml')
        unc,uu,_ = read_grid(RAW/f'{event}_uncertainty.xml')
        for imt,field in IM_FIELDS.items():
            if field not in units or units[field] not in ['%g','g']:
                raise ValueError(f'{event} {field}: IM/birim bulunamadı')
            values = sample_grid(grid,field,b.latitude,b.longitude) / (100 if units[field]=='%g' else 1)
            stdfield = 'STD'+field
            sigma = sample_grid(unc,stdfield,b.latitude,b.longitude) if uu.get(stdfield)=='ln(g)' else np.full(len(b),np.nan)
            im_rows.append(pd.DataFrame(dict(building_id=b.building_id,event=event,imt=imt,im_g=values,
                                             sigma_ln_im=sigma,method='bilinear_grid',grid_sha256=fingerprint(RAW/f'{event}_grid.xml'))))
        if 'SVEL' in units and units['SVEL']=='m/s':
            b[f'vs30_{event}'] = sample_grid(grid,'SVEL',b.latitude,b.longitude)
    ims = pd.concat(im_rows,ignore_index=True)
    csv(ims,out/'im_assignments.csv')
    models = read_model(RAW/'gem_vulnerability_structural.xml')
    scenario_results=[]
    for scenario, aa in a.groupby('scenario'):
        d = b.merge(aa[['building_id','taxonomy','mapping_evidence']],how='left',on='building_id',validate='one_to_one')
        d['imt'] = d.taxonomy.map({k:v['imt'] for k,v in models.items()})
        for ev in ['event1','event2']:
            ii=ims[ims.event==ev][['building_id','imt','im_g','sigma_ln_im']]
            d=d.merge(ii.rename(columns={'im_g':f'im_{ev}','sigma_ln_im':f'sigma_{ev}'}),how='left',on=['building_id','imt'],validate='many_to_one')
            d[f'loss_{ev}']=np.nan;d[f'cov_{ev}']=np.nan
        d['max_im'] = np.maximum(d.im_event1,d.im_event2)  # missing event is NOT ignored
        d['max_loss']=np.nan;d['max_cov']=np.nan;d['out_of_model_domain']=False
        for taxonomy, idx in d.groupby('taxonomy').groups.items():
            if taxonomy not in models: continue
            f=models[taxonomy]
            for ev in ['event1','event2']:
                m,c,flag=evaluate(f,d.loc[idx,f'im_{ev}'],outside)
                d.loc[idx,f'loss_{ev}']=m;d.loc[idx,f'cov_{ev}']=c
                d.loc[idx,'out_of_model_domain'] |= flag
            m,c,flag=evaluate(f,d.loc[idx,'max_im'],outside)
            d.loc[idx,'max_loss']=m;d.loc[idx,'max_cov']=c
        d['status']=np.select([d.taxonomy.isna(),d.imt.isna(),d.gem_ds.isna(),
                               d[['im_event1','im_event2']].isna().any(axis=1),d.max_loss.isna()],
                              ['taxonomy_unassigned','taxonomy_not_in_xml','damage_unmapped','im_missing','outside_xml_domain'],default='eligible')
        # These are explicitly uncalibrated sensitivity models, not damaged-state fragilities.
        predictions={'max_im':d.max_loss,
                     'remaining_value_proxy':1-(1-d.loss_event1)*(1-d.loss_event2),
                     'capped_sum_proxy':np.minimum(1,d.loss_event1+d.loss_event2)}
        for seq,pred in predictions.items():
            z=d.copy();z['scenario']=scenario;z['sequence']=seq;z['predicted_lr']=pred
            z.loc[(z.status=='eligible') & z.predicted_lr.isna(),'status']='sequence_prediction_missing'
            z['residual']=z.observed_lr-z.predicted_lr
            if seq=='max_im':
                lo,hi,valid=moments_interval(z.max_loss,z.max_cov)
                z['xml_pi_low']=lo;z['xml_pi_high']=hi;z['xml_beta_valid']=valid
                z['within_xml_pi']=np.where(valid & z.observed_lr.notna(),(z.observed_lr>=lo)&(z.observed_lr<=hi),np.nan)
            scenario_results.append(z)
    if not scenario_results:
        raise ValueError('İncelenmiş atama yok; sonuç üretilemez')
    results=pd.concat(scenario_results,ignore_index=True)
    csv(results,out/'building_results.csv')
    csv(results.groupby(['scenario','sequence','status'],dropna=False).size().reset_index(name='n'),out/'exclusion_flow.csv')
    eligible=results[results.status=='eligible'].copy()
    stats=[]; aggregates=[]
    for keys,z in eligible.groupby(['scenario','sequence']):
        stats.append(dict(scenario=keys[0],sequence=keys[1],level='building',group='all',**metrics(z)))
        for col in ['taxonomy','structural_system','district','neighborhood']:
            groupcols=['district','neighborhood'] if col=='neighborhood' else [col]
            for name,g in z.dropna(subset=groupcols).groupby(groupcols,dropna=False):
                stats.append(dict(scenario=keys[0],sequence=keys[1],level=col,group=str(name),**metrics(g)))
            if col in ['district','neighborhood']:
                ag=z.dropna(subset=groupcols).groupby(groupcols,dropna=False).agg(observed_lr=('observed_lr','mean'),predicted_lr=('predicted_lr','mean'),n=('building_id','size')).reset_index()
                ag['scenario'],ag['sequence'],ag['level']=keys[0],keys[1],col
                aggregates.append(ag)
                stats.append(dict(scenario=keys[0],sequence=keys[1],level=col+'_means',group='all',**metrics(ag)))
    csv(pd.DataFrame(stats),out/'metrics.csv')
    if aggregates: csv(pd.concat(aggregates,ignore_index=True),out/'area_means.csv')
    residuals=[]
    for keys,z in eligible.groupby(['scenario','sequence']):
        for cov in ['vs30_event1','year_built','storeys','fault_distance_km']:
            if cov not in z: continue
            pairs=z[[cov,'residual']].dropna()
            rho=float(spearmanr(pairs[cov],pairs.residual).statistic) if len(pairs)>2 and pairs[cov].nunique()>1 and pairs.residual.nunique()>1 else np.nan
            residuals.append(dict(scenario=keys[0],sequence=keys[1],covariate=cov,n=len(pairs),spearman=rho))
    csv(pd.DataFrame(residuals),out/'residual_correlations.csv')
    calibration_plots(eligible,out)
    correction_cv(eligible,out)
    coverage=eligible[eligible.sequence=='max_im'].groupby('scenario').agg(n=('building_id','size'),
                n_valid=('within_xml_pi','count'),coverage_proxy=('within_xml_pi','mean'))
    csv(coverage.reset_index(),out/'xml_interval_diagnostic.csv')
    dump(out/'run.json',dict(inputs={str(p):fingerprint(p) for p in [buildings,assignments,damage_map]},
         model_sha=SHA,outside_policy=outside,seed=20260919,
         warning='Damage-derived proxy loss; XML interval diagnostic does not identify beta_bldg. Spatial IM correlations not integrated.'))
    print(f'{len(b)} buildings audited; {len(eligible)} eligible scenario rows')


def calibration_plots(d,out):
    bins=[0,.05,.1,.2,.4,.8,1.6,3.2,6.4,np.inf]
    rows=[]
    # Never pool different IMs on one acceleration axis.
    for idx,(keys,z) in enumerate(d.groupby(['scenario','sequence','imt'])):
        z=z.copy();z['im_bin']=pd.cut(z.max_im,bins,right=False)
        g=z.groupby('im_bin',observed=True).agg(im_g=('max_im','mean'),observed=('observed_lr','mean'),
                                              predicted=('predicted_lr','mean'),n=('building_id','size')).reset_index()
        g['scenario'],g['sequence'],g['imt']=keys
        rows.append(g)
        fig,ax=plt.subplots(figsize=(6,4))
        ax.plot(g.im_g,g.observed,'o-',label='Damage-derived loss proxy')
        ax.plot(g.im_g,g.predicted,'s--',label='GEM prediction / sequence proxy')
        ax.set(xlabel=f'{keys[2]} (g), maximum of two events',ylabel='Mean structural loss ratio',ylim=(0,1),title=' | '.join(map(str,keys)))
        ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/f'calibration_{idx:03d}.png',dpi=180);plt.close(fig)
    if rows:csv(pd.concat(rows,ignore_index=True),out/'im_bin_calibration.csv')


def correction_cv(d,out):
    """Leave-one-district-out factor; entire blocks bootstrap the descriptive factor.

    No Turkey-wide correction is claimed from a convenience sample.
    """
    rng=np.random.default_rng(20260919);rows=[]
    for keys,z in d.groupby(['scenario','sequence']):
        z=z.dropna(subset=['district','observed_lr','predicted_lr'])
        districts=z.district.unique()
        if len(districts)<3: continue
        for held in districts:
            train=z[z.district!=held];test=z[z.district==held].copy()
            if train.predicted_lr.sum()<=0:continue
            factor=train.observed_lr.sum()/train.predicted_lr.sum()
            raw=metrics(test)
            test.predicted_lr=np.minimum(1,factor*test.predicted_lr)
            rows.append(dict(scenario=keys[0],sequence=keys[1],held_out_district=held,factor=factor,
                             raw_rmse=raw['rmse'],**metrics(test)))
        block=z.groupby('district')[['observed_lr','predicted_lr']].sum().to_numpy()
        sums=block[rng.integers(0,len(block),size=(1000,len(block)))].sum(axis=1)
        ratios=sums[:,0]/np.where(sums[:,1]>0,sums[:,1],np.nan)
        q=np.nanquantile(ratios,[.025,.5,.975])
        rows.append(dict(scenario=keys[0],sequence=keys[1],held_out_district='DESCRIPTIVE_BLOCK_BOOTSTRAP',
                         factor=z.observed_lr.sum()/z.predicted_lr.sum(),factor_low=q[0],factor_median=q[1],factor_high=q[2]))
    csv(pd.DataFrame(rows),out/'correction_cross_validation.csv')


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    sub.add_parser('inventory')
    q=sub.add_parser('import');q.add_argument('source');q.add_argument('adapter');q.add_argument('output')
    q=sub.add_parser('analyse');q.add_argument('buildings');q.add_argument('assignments');q.add_argument('damage_map');q.add_argument('output');q.add_argument('--outside',choices=['nan','clamp'],default='nan')
    a=p.parse_args()
    if a.command=='inventory':inventory()
    elif a.command=='import':import_buildings(a.source,a.adapter,a.output)
    else:analyse(a.buildings,a.assignments,a.damage_map,a.output,a.outside)


if __name__=='__main__': main()
