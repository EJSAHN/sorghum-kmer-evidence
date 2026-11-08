"""Denominator-matched primary comparisons from the complete catalogue."""
from collections import Counter, defaultdict

def require(ok,message):
    if not ok: raise ValueError(message)

def primary_headlines(catalog,protocol):
    quality=protocol['primary']['quality'];threshold=protocol['primary']['minimum_support']
    calls=catalog['tables']['full_catalog_calls']
    primary=[r for r in calls if r['quality']==quality and r['minimum_support']==threshold]
    by_label=defaultdict(dict)
    for r in primary:
        key=r['pi_accession'],r['state_id']
        require(key not in by_label[r['baseline']],'Repeated primary cell in a baseline')
        by_label[r['baseline']][key]=r
    fixed=['k'+str(k) for k in protocol['k_values']]
    common=set.intersection(*(set(by_label[b]) for b in fixed))
    rows=[]
    for cohort in ('ALL_ELIGIBLE','COMMON_FIXED_K_COHORT'):
        for baseline in protocol['baseline_labels']:
            rr=list(by_label[baseline].values())
            if cohort=='COMMON_FIXED_K_COHORT':rr=[r for r in rr if (r['pi_accession'],r['state_id']) in common]
            require(rr,'Empty baseline cohort')
            result=dict(cohort=cohort,baseline=baseline,eligible_cells=len(rr),
                distinct_states=len({r['state_id'] for r in rr}),
                excluded_from_full_catalog=protocol['expected_catalog_cells']-len(rr),
                quality=quality,minimum_support=threshold,support_unit='read_pair_observations')
            for name in ('pooled','same_pair','same_read','exact'):
                result[name+'_positive_cells']=sum(bool(r[name+'_pass']) for r in rr)
            for name in ('pooled','same_pair','same_read'):
                result[name+'_positive_without_exact_cells']=sum(bool(r[name+'_pass']) and not r['exact_pass'] for r in rr)
                result['exact_positive_without_'+name+'_cells']=sum(bool(r['exact_pass']) and not r[name+'_pass'] for r in rr)
                require(result['exact_positive_without_'+name+'_cells']==0,'Call nesting failed')
            result['same_pair_count_differs_from_exact_cells']=sum(r['same_pair_pairs']!=r['exact_pairs'] for r in rr)
            result['same_read_count_differs_from_exact_cells']=sum(r['same_read_pairs']!=r['exact_pairs'] for r in rr)
            rows.append(result)
    k31=list(by_label['k31'].values())
    strata=dict(Counter(r['cohort_at_primary_k31'] for r in k31))
    require(len(k31)==protocol['expected_catalog_cells'],'Full catalog denominator differs')
    require(sum(strata.values())==len(k31),'Fixed control strata do not partition the catalog')
    return rows,dict(full_catalog_cells=len(k31),fixed_k_common_cells=len(common),
                     fixed_k_common_states=len({state for pi,state in common}),fixed_control_strata=strata)

