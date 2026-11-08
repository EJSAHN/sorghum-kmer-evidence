"""Frozen, truth-known simulation. All classifiers see the same observed pairs."""
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path
import csv, gzip, hashlib, html, io, itertools, json, os, platform, time
from .detector import Detector, require
from .simulation import build_templates, simulate_cell, stable_seed

METHODS = ('pooled', 'same_pair', 'same_read')
COUNT_KEYS = {'pooled':'minimum_word', 'same_pair':'same_pair', 'same_read':'same_read', 'exact_span':'exact'}

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save_json(path, data):
    path=Path(path); tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n'); os.replace(tmp,path)

@contextmanager
def table_writer(path, fields):
    with Path(path).open('wb') as raw, gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as gz:
        with io.TextIOWrapper(gz,encoding='utf-8',newline='') as out:
            writer=csv.DictWriter(out,fieldnames=fields); writer.writeheader(); yield writer

def write_table(path, rows, fields=None):
    rows=list(rows); fields=fields or list(rows[0])
    with table_writer(path,fields) as out: out.writerows(rows)

def cells(protocol, states):
    for target in sorted(states,key=lambda x:x['state_id']):
        for model, noise, read_length, pair_count, replicate in itertools.product(
            protocol['template_models'],protocol['noise_models'],protocol['read_lengths'],
            protocol['pair_counts'],range(protocol['replicates'])):
            # Noise variants deliberately share generating fragments. They are correlated,
            # not independent trials for biological inference.
            seed=stable_seed(protocol['master_seed'],'cell',target['state_id'],model,read_length,pair_count,replicate)
            cell_id=hashlib.sha256(('|'.join(map(str,(target['state_id'],model,noise,read_length,pair_count,replicate)))).encode()).hexdigest()[:24]
            yield target,dict(cell_id=cell_id,state_id=target['state_id'],template_model=model,noise_model=noise,
                              read_length=read_length,pair_count=pair_count,replicate=replicate,seed=seed)

def add_metric(acc,key,prediction,truth):
    v=acc[key];v['trials']+=1;v['positive_calls']+=int(prediction)
    v['true_positive']+=int(prediction and truth);v['false_negative']+=int(not prediction and truth)
    v['false_positive']+=int(prediction and not truth);v['true_negative']+=int(not prediction and not truth)

def metric_rows(acc,fields):
    for key,v in sorted(acc.items()):
        truth=bool(key[fields.index('truth_present')]); n=v['trials']; positives=v['positive_calls']
        yield dict(zip(fields,key),**v,positive_call_rate=positives/n,
                   rate_interpretation='sensitivity' if truth else 'false_positive_rate')

def run(config):
    started=time.monotonic(); root=Path(config['run_dir']); out=root/'analysis';out.mkdir(exist_ok=True)
    require(not (out/'counts_by_cell.csv.gz').exists(),'Analysis output exists; create a new run instead of overwriting it')
    protocol=json.loads(Path(config['protocol_file']).read_text()); panel=json.loads(Path(config['panel_file']).read_text())
    states=panel['states']; require(len(states)==protocol['expected_targets'],'Target count differs from frozen protocol')
    require(len({s['state_id'] for s in states})==len(states),'Duplicate target identifiers')
    require(len({s['sequence'] for s in states})==len(states),'Duplicate target sequences')
    expected=len(states)*len(protocol['template_models'])*len(protocol['noise_models'])*len(protocol['read_lengths'])*len(protocol['pair_counts'])*protocol['replicates']
    require(expected==protocol['expected_cells'],'Frozen expected cell count differs')
    counts_fields=['cell_id','state_id','template_model','noise_model','read_length','pair_count','replicate','truth_present','k','span_length','quality','distinct_canonical_words','minimum_word_pairs','minimum_word_signatures','same_pair_pairs','same_pair_signatures','same_read_pairs','same_read_signatures','exact_pairs','exact_signatures']
    audit_fields=['cell_id','state_id','template_model','noise_model','read_length','pair_count','replicate','seed','truth_present','audit_json']
    metric_fields=['template_model','noise_model','read_length','pair_count','truth_present','method','k','quality','support_unit','minimum_support']
    primary_fields=['state_id']+metric_fields
    pair_fields=['template_model','noise_model','read_length','pair_count','truth_present','comparator','quality','support_unit','minimum_support']
    metrics=defaultdict(Counter); primary_by_target=defaultdict(Counter); paired=defaultdict(Counter)
    all_templates={}; detectors={}; processed=0; total_pairs=0; row_count=0; exact_invariant_checks=0
    read_digest=hashlib.sha256(); previous_target=None
    with table_writer(out/'counts_by_cell.csv.gz',counts_fields) as count_out, table_writer(out/'generation_audit.csv.gz',audit_fields) as audit_out:
        for target,cell in cells(protocol,states):
            sid=target['state_id']
            if sid!=previous_target:
                # Keep per-cell data out of memory; only the current target matcher lives here.
                detectors={sid:Detector(target['sequence'],protocol['k_values'])}
                models=build_templates(target['sequence'],sid,protocol['master_seed'])
                require(set(models)==set(protocol['template_models']),'Generator models differ from protocol')
                all_templates[sid]=models; previous_target=sid
            model=all_templates[sid][cell['template_model']]; truth=bool(model['truth_present'])
            require(truth==cell['template_model'].startswith('present_'),'Template truth disagrees with model label')
            reads,audit=simulate_cell(model,cell['seed'],cell['pair_count'],cell['read_length'],cell['noise_model'])
            require(len(reads)==cell['pair_count'],'Read-pair count differs')
            # Preserve the entire independent generator audit, including origins and read digest.
            audit_out.writerow(dict(cell,truth_present=truth,audit_json=json.dumps(audit,sort_keys=True,separators=(',',':'))))
            read_digest.update((cell['cell_id']+'\n'+json.dumps(audit,sort_keys=True)+'\n').encode())
            results=detectors[sid].evaluate(reads); lookup={(r['quality'],r['k']):r for r in results}
            base={k:v for k,v in cell.items() if k!='seed'};base['truth_present']=truth
            for r in results:
                count_out.writerow(dict(base,**r));row_count+=1
                if r['k']==len(target['sequence']):
                    for unit in protocol['support_units']:
                        require(len({r[p+'_'+unit] for p in ('minimum_word','same_pair','same_read','exact')})==1,'Full-length/exact equivalence failed')
                        exact_invariant_checks+=1
            for quality,unit,threshold in itertools.product(protocol['quality_policies'],protocol['support_units'],protocol['support_thresholds']):
                predictions={}
                for k in protocol['k_values']:
                    r=lookup[quality,k]
                    for method in METHODS:
                        predictions[method,k]=r[COUNT_KEYS[method]+'_'+unit]>=threshold
                predictions['exact_span',0]=lookup[quality,protocol['k_values'][0]]['exact_'+unit]>=threshold
                for (method,k),prediction in predictions.items():
                    key=(cell['template_model'],cell['noise_model'],cell['read_length'],cell['pair_count'],truth,method,k,quality,unit,threshold)
                    add_metric(metrics,key,prediction,truth)
                    if (quality,unit,threshold)==(protocol['primary']['quality'],protocol['primary']['support_unit'],protocol['primary']['minimum_support']):
                        add_metric(primary_by_target,(sid,)+key,prediction,truth)
                base_prediction=predictions['pooled',31]
                for method,k in [('pooled',81),('same_pair',31),('same_read',31),('exact_span',0)]:
                    label=method if k==0 else method+'_k'+str(k)
                    key=(cell['template_model'],cell['noise_model'],cell['read_length'],cell['pair_count'],truth,label,quality,unit,threshold)
                    other=predictions[method,k]; pv=paired[key];pv['trials']+=1
                    pv['both_positive' if base_prediction and other else 'pooled31_only' if base_prediction else 'comparator_only' if other else 'both_negative']+=1
            processed+=1;total_pairs+=len(reads)
            if processed%100==0 or processed==expected:
                elapsed=time.monotonic()-started
                progress=dict(status='RUNNING',phase='simulation',completed_cells=processed,total_cells=expected,completed_read_pairs=total_pairs,
                    elapsed_seconds=elapsed,estimated_remaining_analysis_seconds=elapsed*(expected-processed)/processed,
                    estimate_scope='Current process throughput only; excludes queue and final collection.')
                save_json(root/'progress.json',progress)
                if processed%1000==0 or processed==expected: print('PROGRESS cells=%d/%d elapsed=%.1fs'%(processed,expected,elapsed),flush=True)
    require(processed==expected,'Incomplete simulation grid')
    require(total_pairs==protocol['expected_read_pairs'],'Unexpected total simulated read pairs')
    mrows=list(metric_rows(metrics,metric_fields));write_table(out/'performance_by_stratum.csv.gz',mrows)
    write_table(out/'primary_per_target.csv.gz',metric_rows(primary_by_target,primary_fields))
    paired_rows=[]
    for key,v in sorted(paired.items()):
        require(sum(v[k] for k in ('both_positive','pooled31_only','comparator_only','both_negative'))==v['trials'],'Paired comparison denominator mismatch')
        paired_rows.append(dict(zip(pair_fields,key),**{k:v[k] for k in ('trials','both_positive','pooled31_only','comparator_only','both_negative')}))
    write_table(out/'paired_method_comparisons.csv.gz',paired_rows)
    save_json(out/'generating_templates.json',all_templates)
    # Verify noise-paired geometry from generator-provided origin records/digests independently.
    # All audit fields are retained, so replay does not depend on rounded CSV summaries.
    geometry_checks=verify_geometry(out/'generation_audit.csv.gz',protocol)
    primary=[r for r in mrows if (r['quality'],r['support_unit'],r['minimum_support'])==(protocol['primary']['quality'],protocol['primary']['support_unit'],protocol['primary']['minimum_support'])]
    summary=dict(status='CODE25_COMPLETE',schema='CODE25_SIMULATION_RESULT_V1',synthetic_only=True,targets=len(states),target_length_counts=dict(Counter(len(s['sequence']) for s in states)),
        simulation_cells=processed,simulated_read_pairs=total_pairs,counts_rows=row_count,performance_rows=len(mrows),primary_per_target_rows=len(primary_by_target),
        exact_equivalence_checks=exact_invariant_checks,noise_geometry_checks=geometry_checks,master_seed=protocol['master_seed'],
        protocol_sha256=digest(config['protocol_file']),panel_sha256=digest(config['panel_file']),generation_audit_binding_sha256=read_digest.hexdigest(),
        primary=protocol['primary'],elapsed_seconds=time.monotonic()-started,python_version=platform.python_version(),hostname=platform.node(),
        biological_validation=False,GWAS=False,original_cases_reclassified=False,downloads=False,original_fastq_rescan=False,
        interpretation='Stratum-specific synthetic detection and false-positive rates on a selected target panel. No single overall accuracy or biological false-positive-rate estimate.',limitations=protocol['scope_limits'])
    save_json(out/'benchmark_summary.json',summary)
    save_json(root/'progress.json',dict(status='COMPLETE',phase='simulation',completed_cells=processed,total_cells=expected,elapsed_seconds=time.monotonic()-started))
    return summary

def verify_geometry(path,protocol):
    by_cell={}; noise_seen=defaultdict(set); checks=0
    with gzip.open(path,'rt',newline='') as f:
        for r in csv.DictReader(f):
            key=tuple(r[k] for k in ('state_id','template_model','read_length','pair_count','replicate'))
            require(r['noise_model'] in protocol['noise_models'],'Unknown noise condition in audit')
            require(r['noise_model'] not in noise_seen[key],'Duplicate noise condition in audit')
            noise_seen[key].add(r['noise_model'])
            audit=json.loads(r['audit_json'])
            origin=audit.get('origin_sha256',audit.get('origins_sha256'))
            require(origin is not None,'Generator did not provide origin SHA256')
            if key in by_cell:
                require(by_cell[key]==origin,'Noise conditions changed generating fragment geometry');checks+=1
            else:by_cell[key]=origin
    require(all(v==set(protocol['noise_models']) for v in noise_seen.values()),'Missing noise-paired audit rows')
    require(checks==protocol['expected_cells']*(len(protocol['noise_models'])-1)//len(protocol['noise_models']),'Missing noise-paired audit rows')
    return checks



def context_tables(analysis_dir):
    """Two descriptive tables computed from the completed benchmark outputs.

    dilution_sampling_ceiling: in the diluted target-present model, how many
    trials sampled zero, one, or at least two read pairs from the target-bearing
    template. Trials with fewer than two such pairs cannot meet a two-observation
    criterion from target-origin reads alone.

    support_unit_agreement: how often a call differs between the read-pair
    observation unit and the distinct ordered sequence-signature unit.
    """
    analysis_dir = Path(analysis_dir)
    ceiling = defaultdict(Counter)
    with gzip.open(analysis_dir / 'generation_audit.csv.gz', 'rt', newline='') as handle:
        for row in csv.DictReader(handle):
            if row['template_model'] != 'present_diluted10':
                continue
            audit = json.loads(row['audit_json'])
            target_pairs = int(audit['selected_target_template_pairs'])
            key = (int(row['read_length']), int(row['pair_count']), row['noise_model'])
            ceiling[key]['trials'] += 1
            ceiling[key]['target_template_pairs_0'] += target_pairs == 0
            ceiling[key]['target_template_pairs_1'] += target_pairs == 1
            ceiling[key]['target_template_pairs_at_least_2'] += target_pairs >= 2
    ceiling_rows = [dict(template_model='present_diluted10', read_length=rl, pair_count=pc, noise_model=noise,
                         **{name: counts[name] for name in ('trials', 'target_template_pairs_0',
                                                            'target_template_pairs_1',
                                                            'target_template_pairs_at_least_2')})
                    for (rl, pc, noise), counts in sorted(ceiling.items())]
    agreement = defaultdict(Counter)
    with gzip.open(analysis_dir / 'counts_by_cell.csv.gz', 'rt', newline='') as handle:
        for row in csv.DictReader(handle):
            for criterion, key in COUNT_KEYS.items():
                pairs = int(row[key + '_pairs']); signatures = int(row[key + '_signatures'])
                require(signatures <= pairs, 'Signature count exceeds observation count')
                for threshold in (1, 2, 3):
                    bucket = agreement[(int(row['read_length']), row['quality'], criterion, threshold)]
                    bucket['rows'] += 1
                    bucket['calls_differ'] += (pairs >= threshold) != (signatures >= threshold)
    agreement_rows = [dict(read_length=rl, quality=q, criterion=c, minimum_support=t,
                           count_rows=v['rows'], calls_differ_between_units=v['calls_differ'])
                      for (rl, q, c, t), v in sorted(agreement.items())]
    write_table(analysis_dir / 'dilution_sampling_ceiling.csv.gz', ceiling_rows)
    write_table(analysis_dir / 'support_unit_agreement.csv.gz', agreement_rows)
    return dict(dilution_sampling_ceiling=len(ceiling_rows), support_unit_agreement=len(agreement_rows))
