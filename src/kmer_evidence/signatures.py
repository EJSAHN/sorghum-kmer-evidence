"""Reconcile selected-case signatures without extending them to unmeasured cells.

Whole-catalogue counts are read-pair observations. Sequence signatures are a
separate complete census for selected accession-state cells only. Simulations
provide context and remain outside every real-data denominator.
"""
from collections import Counter, defaultdict
from pathlib import Path
import csv
import gzip
import json


COUNT_PREFIXES = ('minimum_word', 'same_pair', 'same_read', 'exact')
PASS_NAMES = ('pooled', 'same_pair', 'same_read', 'exact')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def number(value):
    require(not isinstance(value, bool), 'Boolean supplied as a count')
    if isinstance(value, int):
        require(value >= 0, 'Negative count')
        return value
    require(isinstance(value, str) and value.isdecimal(), 'Expected a nonnegative integer count')
    return int(value)


def boolean(value):
    if type(value) is bool:
        return value
    require(value in ('True', 'False'), 'Invalid literal boolean')
    return value == 'True'


def read_table(path):
    with gzip.open(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle)


def index_rows(rows, keys, label):
    result = {}
    for row in rows:
        key = tuple(row[k] for k in keys)
        require(key not in result, 'Duplicate '+label+' identity: '+repr(key))
        result[key] = row
    return result


def relation(counts, threshold):
    if counts['exact'] >= threshold:
        return 'CONCORDANT_POSITIVE'
    if counts['same_read'] >= threshold:
        return 'SAME_READ_WITHOUT_EXACT_SPAN'
    if counts['same_pair'] >= threshold:
        return 'MATE_UNION_ONLY'
    if counts['minimum_word'] >= threshold:
        return 'POOLED_WITHOUT_SAME_PAIR'
    return 'ALL_BELOW_THRESHOLD'


def audit_simulation(source25, primary):
    summary = json.loads((Path(source25)/'CODE25_summary.json').read_text())
    require(summary.get('status') == 'CODE25_COMPLETE', 'CODE25 did not complete')
    analysis = summary['analysis']
    require(analysis.get('status') == 'CODE25_COMPLETE' and analysis.get('synthetic_only') is True,
            'CODE25 synthetic analysis status differs')
    rows = list(read_table(Path(source25)/'analysis/performance_by_stratum.csv.gz'))
    key_fields = ('template_model', 'noise_model', 'read_length', 'pair_count',
                  'quality', 'support_unit', 'minimum_support')
    seen = set(); paired = defaultdict(dict); context = []; exact_totals = Counter()
    for raw in rows:
        row = dict(raw)
        for key in ('read_length', 'pair_count', 'k', 'minimum_support', 'trials', 'positive_calls',
                    'true_positive', 'false_negative', 'false_positive', 'true_negative'):
            row[key] = number(row[key])
        row['truth_present'] = boolean(row['truth_present'])
        require(row['trials'] > 0 and row['positive_calls'] <= row['trials'], 'Invalid simulation denominator')
        require(row['true_positive']+row['false_negative']+row['false_positive']+row['true_negative'] == row['trials'],
                'Simulation confusion counts do not sum to trials')
        require(row['true_positive']+row['false_positive'] == row['positive_calls'], 'Simulation positives differ')
        if row['truth_present']:
            require(row['false_positive'] == row['true_negative'] == 0, 'Wrong truth stratum in simulation')
        else:
            require(row['true_positive'] == row['false_negative'] == 0, 'Wrong truth stratum in simulation')
        require(abs(float(row['positive_call_rate'])-row['positive_calls']/row['trials']) < 1e-12,
                'Simulation rate differs from numerator/denominator')
        row['positive_call_rate'] = row['positive_calls']/row['trials']
        key = tuple(row[k] for k in key_fields)
        identity = key+(row['method'], row['k'])
        require(identity not in seen, 'Duplicate simulation performance row')
        seen.add(identity)
        if (row['method'], row['k']) in (('same_read', 31), ('exact_span', 0)):
            paired[key][row['method']] = row
        if (row['method'], row['k']) == ('exact_span', 0):
            exact_totals[row['quality'], row['support_unit'], row['minimum_support']] += row['trials']
        if (row['quality'], row['support_unit'], row['minimum_support']) == (primary['quality'], 'signatures', primary['minimum_support']):
            row['evidence_domain'] = 'SYNTHETIC_CODE25_ONLY'
            row['real_catalog_denominator_included'] = False
            context.append(row)
    expected = number(analysis['simulation_cells'])
    require(exact_totals and all(n == expected for n in exact_totals.values()), 'Simulation sample denominators differ')
    differences = 0
    for key, pair in paired.items():
        require(set(pair) == {'same_read', 'exact_span'}, 'Missing simulation comparator')
        first, second = pair['same_read'], pair['exact_span']
        require(first['trials'] == second['trials'] and first['truth_present'] == second['truth_present'],
                'Simulation comparator denominators/truth differ')
        differences += int(first['positive_calls'] != second['positive_calls'])
    return context, dict(simulation_stratum_comparisons=len(paired),
        simulation_same_read31_vs_exact_positive_count_differences=differences,
        simulation_comparison_scope='Aggregate positive-call counts, not a new per-cell replay or inference of per-cell identity.',
        simulation_primary_rows=len(context), simulation_samples=expected,
        synthetic_results_in_real_catalog_denominators=False)


def audit_overlay(source24, source25, native_evidence, full_catalog_calls, protocol):
    base_k = min(protocol['k_values']); qualities = tuple(protocol['qualities'])
    thresholds = tuple(protocol['support_thresholds']); primary = protocol['primary']
    native = {}
    for raw in native_evidence:
        row = dict(raw); row['k'] = number(row['k']); row['span_length'] = number(row['span_length'])
        for name in COUNT_PREFIXES:
            row[name+'_pairs'] = number(row[name+'_pairs'])
        key = (row['pi_accession'], row['state_id'], row['k'], row['quality'])
        require(key not in native, 'Duplicate native overlay identity')
        native[key] = row
    base = {(pi, sid, q):r for (pi, sid, k, q),r in native.items() if k == base_k}
    cells = {(pi, sid) for pi, sid, q in base}
    require(len(cells) == protocol['expected_catalog_cells'], 'Catalog cell coverage differs')
    require(set(base) == {(pi, sid, q) for pi,sid in cells for q in qualities}, 'Missing catalog quality/base row')
    require({(pi,sid) for pi,sid,k,q in native} == cells, 'Native cell lacks base-k coverage')
    primary_calls = {}
    for r in full_catalog_calls:
        threshold = number(r.get('minimum_support', r.get('minimum_supporting_pairs')))
        if r['baseline'] == 'k'+str(base_k) and r['quality'] == primary['quality'] and threshold == primary['minimum_support']:
            key = (r['pi_accession'], r['state_id'])
            require(key not in primary_calls, 'Duplicate primary catalog call')
            primary_calls[key] = r
    require(set(primary_calls) == cells, 'Missing primary catalog call')
    selected_expected = {cell for cell in cells if base[cell+(primary['quality'],)]['minimum_word_pairs'] >= primary['minimum_support']
                         and base[cell+(primary['quality'],)]['exact_pairs'] < primary['minimum_support']}
    require(len(selected_expected) == protocol['expected_selected_cells'], 'Selected-case count differs')
    for cell,r in primary_calls.items():
        measured = base[cell+(primary['quality'],)]
        expected_cohort = ('CONCORDANT_POSITIVE' if measured['exact_pairs'] >= primary['minimum_support'] else
                           'DISCORDANT' if measured['minimum_word_pairs'] >= primary['minimum_support'] else 'BELOW_THRESHOLD')
        require(r['cohort_at_primary_k31'] == expected_cohort, 'Primary control-stratum assignment differs')

    root24 = Path(source24)
    ledger = index_rows(read_table(root24/'analysis/case_evidence_ledger.csv.gz'), ('pi_accession','state_id'), 'CODE24 ledger')
    require(set(ledger) == selected_expected, 'CODE24 ledger selection differs from full-catalog discrepancies')
    rows = list(read_table(root24/'analysis/all_baselines.csv.gz'))
    indexed = {}; stable_counts = {}; selected_rows = []; aggregates = defaultdict(Counter); reconciled = 0
    for raw in rows:
        row = dict(raw); cell = (row['pi_accession'], row['state_id'])
        require(cell in selected_expected, 'Unexpected CODE24 selected cell')
        k = number(row['k']); span = number(row['span_length']); q = row['quality']; threshold = number(row['minimum_support'])
        require(q in qualities and threshold in thresholds, 'Unexpected CODE24 quality/support setting')
        require(k in set(protocol['k_values'])|{span}, 'Unexpected CODE24 k')
        key = cell+(k,q,threshold)
        require(key not in indexed, 'Duplicate CODE24 baseline row')
        indexed[key] = row
        base_row = base[cell+(q,)]
        require(span == base_row['span_length'] and row['pair_id'] == base_row['pair_id'], 'CODE24 target identity differs')
        require(row['group_id'] == base_row['group_id'], 'CODE24 base-group identity differs')
        require(ledger[cell]['pair_id'] == row['pair_id'], 'CODE24 ledger pair identity differs')
        if k == span:
            expected_counts = {name:base_row['exact_pairs'] for name in COUNT_PREFIXES}
            baseline = 'full_span'
        else:
            native_key = cell+(k,q)
            require(native_key in native, 'Missing native counter for selected baseline')
            expected_counts = {name:native[native_key][name+'_pairs'] for name in COUNT_PREFIXES}
            baseline = 'k'+str(k)
        values = {}
        for unit in ('pairs','signatures'):
            values[unit] = {name:number(row[name+'_'+unit]) for name in COUNT_PREFIXES}
            c = values[unit]
            require(c['exact'] <= c['same_read'] <= c['same_pair'] <= c['minimum_word'], 'CODE24 evidence nesting differs')
        require(values['pairs'] == expected_counts, 'CODE24 raw counters disagree with full-catalog native evidence: '+repr(key))
        reconciled += len(COUNT_PREFIXES)
        require(all(values['signatures'][name] <= values['pairs'][name] for name in COUNT_PREFIXES), 'Signature count exceeds raw pair count')
        if k == span:
            require(len(set(values['signatures'].values())) == 1, 'Full-span signature equivalence differs')
        stable_key = cell+(k,q)
        signature_tuple = tuple(values[unit][name] for unit in ('pairs','signatures') for name in COUNT_PREFIXES)
        require(stable_key not in stable_counts or stable_counts[stable_key] == signature_tuple, 'Counter changes across support thresholds')
        stable_counts[stable_key] = signature_tuple
        for name,prefix in zip(PASS_NAMES,COUNT_PREFIXES):
            require(boolean(row[name+'_pass_signatures']) == (values['signatures'][prefix] >= threshold), 'CODE24 signature pass flag differs')
        clean = dict(pi_accession=cell[0],state_id=cell[1],pair_id=row['pair_id'],baseline=baseline,k=k,
            span_length=span,quality=q,minimum_support=threshold,
            code24_primary_evidence_label=ledger[cell].get('primary_evidence_label',''),
            code24_placement_status=ledger[cell].get('placement_status',''),
            mechanism_scope=ledger[cell].get('mechanism_scope',''),
            signature_scope='COMPLETE_SELECTED_CELL_CENSUS_ONLY')
        for unit in ('pairs','signatures'):
            clean.update({name+'_'+unit:values[unit][name] for name in COUNT_PREFIXES})
            clean[unit+'_evidence_relation'] = relation(values[unit],threshold)
            a = aggregates[baseline,q,threshold,'read_pair_observations' if unit=='pairs' else 'ordered_sequence_signatures']
            a['cells'] += 1
            for label,name in zip(PASS_NAMES,COUNT_PREFIXES):
                a[label+'_supported_cells'] += int(values[unit][name] >= threshold)
            a['pooled_without_exact_cells'] += int(values[unit]['minimum_word'] >= threshold and values[unit]['exact'] < threshold)
            a['same_read_without_exact_cells'] += int(values[unit]['same_read'] >= threshold and values[unit]['exact'] < threshold)
        selected_rows.append(clean)
    expected_keys = {cell+(k,q,t) for cell in cells if cell in selected_expected
                     for k in set(protocol['k_values'])|{base[cell+(qualities[0],)]['span_length']}
                     for q in qualities for t in thresholds}
    require(set(indexed) == expected_keys, 'Missing or extra CODE24 selected baseline combinations')
    if 'expected_selected_states' in protocol:
        require(len({sid for pi,sid in selected_expected}) == protocol['expected_selected_states'], 'Selected target-state count differs')

    # The evidence ledger is not a substitute for numerical measurements.
    for cell,l in ledger.items():
        r = indexed[cell+(base_k,primary['quality'],primary['minimum_support'])]
        for ledger_name,source_name in [('q20_minimum_k31_signatures','minimum_word_signatures'),
                                        ('q20_same_pair_signatures','same_pair_signatures'),
                                        ('q20_same_read_signatures','same_read_signatures')]:
            if ledger_name in l:
                require(number(l[ledger_name]) == number(r[source_name]), 'CODE24 ledger signature counter differs')

    coverage = []
    for cell in sorted(cells):
        b = base[cell+(primary['quality'],)]; measured = cell in selected_expected
        r = indexed.get(cell+(base_k,primary['quality'],primary['minimum_support']))
        coverage.append(dict(pi_accession=cell[0],state_id=cell[1],pair_id=b['pair_id'],
            cohort_at_primary_k31=primary_calls[cell]['cohort_at_primary_k31'],
            signature_available=measured,signature_measurement_scope='COMPLETE_SELECTED_CELL_CENSUS' if measured else 'NOT_MEASURED',
            **{name+'_signatures':number(r[name+'_signatures']) if r else '' for name in COUNT_PREFIXES}))
    sensitivity = [dict(baseline=key[0],quality=key[1],minimum_support=key[2],support_unit=key[3],
                        selected_case_scope='CODE24_SELECTED_CELLS_ONLY',**value) for key,value in sorted(aggregates.items())]
    context = []
    for cell in sorted(cells):
        for q in qualities:
            b = base[cell+(q,)]
            if b['same_read_pairs'] <= b['exact_pairs']:
                continue
            for native_key,r in sorted(native.items()):
                if native_key[:2] != cell or native_key[3] != q:
                    continue
                match = indexed.get(cell+(r['k'],q,primary['minimum_support']))
                l = ledger.get(cell,{})
                context.append(dict(pi_accession=cell[0],state_id=cell[1],pair_id=r['pair_id'],k=r['k'],quality=q,
                    cohort_at_primary_k31=primary_calls[cell]['cohort_at_primary_k31'],
                    **{name+'_pairs':r[name+'_pairs'] for name in COUNT_PREFIXES},
                    **{name+'_signatures':number(match[name+'_signatures']) if match else '' for name in COUNT_PREFIXES},
                    signature_measurement_scope='COMPLETE_SELECTED_CELL_AND_K_CENSUS' if match else 'NOT_MEASURED',
                    code24_primary_evidence_label=l.get('primary_evidence_label','NOT_ASSESSED_IN_CODE24'),
                    code24_extra_repeat_example=l.get('code23_extra_repeat_example',''),
                    code24_placement_status=l.get('placement_status',''),
                    causal_attribution='Inherited descriptive context only; no genomic absence or genotype truth inferred.'))
    simulation, simulation_summary = audit_simulation(source25,primary)
    summary = dict(selected_cells=len(selected_expected),selected_states=len({sid for pi,sid in selected_expected}),
        catalog_cells=len(cells),signature_not_measured_cells=len(cells)-len(selected_expected),
        selected_baseline_rows=len(selected_rows),raw_counter_reconciliations=reconciled,
        selected_signature_sensitivity_rows=len(sensitivity),incremental_case_context_rows=len(context),
        incremental_count_surplus_cells=len({(r['pi_accession'],r['state_id']) for r in context}),
        signature_scope='Complete only for selected accession/state cells; unknown cells remain blank.',
        **simulation_summary)
    return {'tables':dict(signature_coverage=coverage,selected_signature_sensitivity=sensitivity,
        selected_case_baselines=selected_rows,incremental_case_context=context,simulation_context=simulation),
        'summary':summary}
