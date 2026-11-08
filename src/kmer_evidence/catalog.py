"""Recompute whole-catalogue evidence from complete read-pair counters.

Counts here are FASTQ pair observations. They are not deduplicated sequence
signatures, independent molecules, genotypes, or truth labels. Missing k>=span
baselines remain unavailable. Native groups shared by two baseline labels are
represented once in native_evidence, with explicit aliases in comparison views.
"""
from array import array
from collections import Counter, defaultdict
from pathlib import Path
import csv
import gzip
import json
import re

COUNT_BASIS = 'READ_PAIR_OBSERVATIONS'
METHOD_FIELDS = {'pooled':'minimum_word_pairs', 'same_pair':'same_pair_pairs',
                 'same_read':'same_read_pairs', 'exact':'exact_pairs'}
CLASSES = ('CONCORDANT_POSITIVE', 'SAME_READ_WITHOUT_EXACT_SPAN',
           'MATE_UNION_ONLY', 'POOLED_WITHOUT_SAME_PAIR', 'ALL_BELOW_THRESHOLD')
RC = str.maketrans('ACGT', 'TGCA')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(sequence):
    return min(sequence, sequence.translate(RC)[::-1])


def read_table(path):
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)),
                'Missing or duplicate CSV columns: '+str(path))
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()),
                    'Malformed CSV row: '+str(path))
            yield row


def as_int(value):
    require(type(value) is int or (isinstance(value, str) and re.fullmatch(r'0|[1-9][0-9]*', value)),
            'Expected nonnegative integer count: '+repr(value))
    number = int(value)
    require(0 <= number < 2**64, 'Count outside uint64 range')
    return number


def as_bool(value):
    if type(value) is bool:
        return value
    require(value in ('True', 'False'), 'Expected literal boolean: '+repr(value))
    return value == 'True'


def unique_index(rows, fields):
    result = {}
    for row in rows:
        key = tuple(row[f] for f in fields)
        require(key not in result, 'Duplicate input identity: '+repr(key))
        result[key] = row
    return result


def classification(evidence, threshold):
    if evidence['exact_pairs'] >= threshold:
        return 'CONCORDANT_POSITIVE'
    if evidence['same_read_pairs'] >= threshold:
        return 'SAME_READ_WITHOUT_EXACT_SPAN'
    if evidence['same_pair_pairs'] >= threshold:
        return 'MATE_UNION_ONLY'
    if evidence['minimum_word_pairs'] >= threshold:
        return 'POOLED_WITHOUT_SAME_PAIR'
    return 'ALL_BELOW_THRESHOLD'


def _primary(protocol):
    primary = protocol.get('primary', {})
    return primary.get('quality', 'q20'), int(primary.get('minimum_support', primary.get('minimum_supporting_pairs', 2)))


def audit_catalog(source19, protocol):
    """Return typed tables and a summary, after independent source reconciliation.

    Protocol expected-count fields are mandatory, so absent all-zero rows cannot
    silently disappear. Tests may supply a smaller complete, internally valid
    panel. Production uses the unchanged CODE19 fixed-grid eligibility k<span.
    """
    source19 = Path(source19)
    k_values = tuple(int(k) for k in protocol['k_values'])
    qualities = tuple(protocol.get('qualities', ('raw', 'q20')))
    thresholds = tuple(int(t) for t in protocol['support_thresholds'])
    require(k_values and list(k_values) == sorted(set(k_values)) and min(k_values) > 0,
            'Invalid fixed k grid')
    require(set(qualities) == {'raw', 'q20'} and len(qualities) == 2, 'Both raw and q20 are required')
    require(thresholds and len(set(thresholds)) == len(thresholds) and min(thresholds) > 0,
            'Invalid support thresholds')
    primary_q, primary_t = _primary(protocol)
    require(primary_q in qualities and primary_t in thresholds, 'Primary setting is absent from sensitivity grid')
    baseline_names = tuple('k'+str(k) for k in k_values) + ('near_full',)
    base_label, long_label = baseline_names[0], baseline_names[-2]
    common_label = 'COMMON_K31_K81' if (k_values[0], k_values[-1]) == (31, 81) else 'COMMON_FIXED_GRID'
    states = {key[0]: value for key, value in unique_index(read_table(source19/'target_sequences.csv.gz'), ('state_id',)).items()}
    require(len(states) == protocol['expected_catalog_states'], 'Catalog state count differs')
    if 'expected_catalog_cells' in protocol:
        require(len(states)*protocol['expected_accessions'] == protocol['expected_catalog_cells'], 'Expected catalog cell denominator differs')
    for sid, state in states.items():
        seq = state['sequence']
        require(seq and set(seq) <= set('ACGT'), 'Invalid target sequence: '+sid)
        require(as_int(state['span_length']) == len(seq), 'Target length mismatch: '+sid)
        if 'normalized_span' in state:
            require(state['normalized_span'] == seq, 'Normalized target sequence differs: '+sid)
        require(k_values[0] < len(seq), 'Every target must support the first fixed k')

    queries = {key[0]: value for key, value in unique_index(read_table(source19/'queries.csv.gz'), ('pattern_id',)).items()}
    require(len(queries) == protocol['expected_queries'], 'Query count differs')
    query_ids = {qid: i for i, qid in enumerate(queries)}
    for qid, query in queries.items():
        seq = query['sequence']; k = as_int(query['k'])
        require(seq and set(seq) <= set('ACGT') and k == len(seq) and seq == canonical(seq),
                'Invalid canonical query: '+qid)
        require(query['query_kind'] in ('word', 'span'), 'Invalid query kind')

    groups = {key[0]: value for key, value in unique_index(read_table(source19/'groups.csv.gz'), ('group_id',)).items()}
    require(len(groups) == protocol['expected_native_groups'], 'Native group count differs')
    group_lookup = {}; eligible = {}; group_word_indexes = {}; span_indexes = {}
    group_labels = {}; state_groups = defaultdict(list)
    for gid, group in groups.items():
        sid = group['state_id']; require(sid in states, 'Unknown group state: '+sid)
        state = states[sid]; seq = state['sequence']; k = as_int(group['k'])
        require(group['pair_id'] == state['pair_id'] and as_int(group['span_length']) == len(seq), 'Group target metadata differs')
        require(0 < k < len(seq), 'Native k must be below target length')
        require((sid, k) not in group_lookup, 'Duplicate state/k native group')
        group_lookup[sid, k] = gid; state_groups[sid].append((k, gid))
        words = json.loads(group['word_ids'])
        require(isinstance(words, list) and words and len(words) == len(set(words)) and all(w in queries for w in words),
                'Duplicate or missing group word identity')
        expected_words = {canonical(seq[i:i+k]) for i in range(len(seq)-k+1)}
        require({queries[w]['sequence'] for w in words} == expected_words and len(words) == len(expected_words),
                'Group target word set differs')
        require(all(as_int(queries[w]['k']) == k and queries[w]['query_kind'] == 'word' for w in words),
                'Group word kind or length differs')
        span = group['span_query_id']
        require(span in queries and queries[span]['sequence'] == canonical(seq) and queries[span]['query_kind'] == 'span',
                'Group exact span query differs')
        labels = []
        if as_bool(group['is_fixed_grid']):
            require(k in k_values, 'Unexpected fixed-grid k')
            labels.append('k'+str(k))
        if as_bool(group['is_near_full']):
            require(k == len(seq)-1, 'Near-full k differs from span minus one')
            labels.append('near_full')
        require(labels and json.loads(group['baseline_labels']) == labels, 'Baseline aliases differ from native group flags')
        group_labels[gid] = labels; group_word_indexes[gid] = [query_ids[w] for w in words]
        span_indexes[gid] = query_ids[span]
        for label in labels:
            require((sid, label) not in eligible, 'Duplicate native group for baseline')
            eligible[sid, label] = gid
    expected_eligible = {(sid, 'k'+str(k)) for sid, s in states.items() for k in k_values if k < len(s['sequence'])}
    expected_eligible |= {(sid, 'near_full') for sid in states}
    require(set(eligible) == expected_eligible, 'Missing or extra eligible baseline groups')
    if 'expected_comparison_rows' in protocol:
        require(len(eligible)*protocol['expected_accessions']*len(qualities)*len(thresholds) == protocol['expected_comparison_rows'],
                'Expected comparison row denominator differs')
    old_eligible = unique_index(read_table(source19/'baseline_eligibility.csv.gz'), ('state_id', 'baseline'))
    require(set(old_eligible) == expected_eligible, 'Source baseline eligibility differs')
    for key, row in old_eligible.items():
        gid = eligible[key]; group = groups[gid]
        require(as_bool(row['eligible']) and row['group_id'] == gid and as_int(row['k']) == as_int(group['k'])
                and row['pair_id'] == group['pair_id'] and as_int(row['span_length']) == as_int(group['span_length']),
                'Source eligibility metadata differs')

    # Arrays retain complete zero-valued queries without retaining 847,864 CSV
    # dictionaries. Seen bitmaps distinguish a missing row from a real zero.
    word_arrays = {}; seen_queries = {}; word_rows = 0; word_field_checks = 0
    for row in read_table(source19/'word_counts.csv.gz'):
        qid = row['pattern_id']; pi = row['pi_accession']
        require(qid in queries, 'Unknown counted query')
        query = queries[qid]; qi = query_ids[qid]
        for field in ('sequence', 'k', 'query_kind', 'catalog_owner_count'):
            if field in query:
                require(row[field] == query[field], 'Counted query metadata differs: '+qid)
        if pi not in word_arrays:
            require(len(word_arrays) < protocol['expected_accessions'], 'Too many accessions')
            word_arrays[pi] = (array('Q', [0])*len(queries), array('Q', [0])*len(queries))
            seen_queries[pi] = bytearray(len(queries))
        require(not seen_queries[pi][qi], 'Duplicate word-count identity')
        seen_queries[pi][qi] = 1
        values = {q+'_'+m: as_int(row[q+'_'+m]) for q in qualities for m in ('occurrences', 'read_support', 'pair_support')}
        for q in qualities:
            occ, reads, pairs = (values[q+'_'+m] for m in ('occurrences', 'read_support', 'pair_support'))
            require(occ >= reads >= pairs and reads <= 2*pairs, 'Query pair/read/occurrence nesting failed')
        for suffix in ('occurrences', 'read_support', 'pair_support'):
            require(values['q20_'+suffix] <= values['raw_'+suffix], 'Query quality nesting failed')
        word_arrays[pi][0][qi] = values['raw_pair_support']; word_arrays[pi][1][qi] = values['q20_pair_support']
        word_rows += 1; word_field_checks += 9
    accessions = sorted(word_arrays)
    require(len(accessions) == protocol['expected_accessions'], 'Accession count differs')
    if 'accessions' in protocol:
        require(accessions == sorted(protocol['accessions']), 'Accession identities differ')
    require(all(all(bits) for bits in seen_queries.values()), 'Missing query/accession count, including possible zero rows')
    require(word_rows == len(queries)*len(accessions), 'Incomplete word-count matrix')
    del seen_queries

    evidence = {}; native = []; group_seen = set(); group_rows = 0; group_field_checks = 0
    for row in read_table(source19/'group_counts.csv.gz'):
        pi, gid = row['pi_accession'], row['group_id']
        require(pi in word_arrays and gid in groups, 'Unknown group/accession count')
        require((pi, gid) not in group_seen, 'Duplicate group/accession count');group_seen.add((pi, gid))
        group = groups[gid]; sid = group['state_id']; k = as_int(group['k'])
        for field, value in group.items():
            require(row[field] == value, 'Group count metadata differs: '+field)
        for suffix in ('same_read_complete_reads', 'same_read_complete_pairs', 'mate_union_complete_pairs', 'exact_span_pairs'):
            require(as_int(row['q20_'+suffix]) <= as_int(row['raw_'+suffix]), 'Group quality nesting failed')
        for quality in qualities:
            nums = word_arrays[pi][0 if quality == 'raw' else 1]
            minimum = min(nums[i] for i in group_word_indexes[gid])
            same = as_int(row[quality+'_same_read_complete_pairs'])
            union = as_int(row[quality+'_mate_union_complete_pairs'])
            exact = as_int(row[quality+'_exact_span_pairs'])
            reads = as_int(row[quality+'_same_read_complete_reads'])
            require(exact <= same <= union <= minimum, 'Nested evidence criteria failed')
            require(same <= reads <= 2*same, 'Group read/pair units differ')
            require(exact == nums[span_indexes[gid]], 'Group exact count differs from literal span query')
            record = dict(pi_accession=pi, state_id=sid, pair_id=group['pair_id'], group_id=gid,
                          k=k, span_length=len(states[sid]['sequence']), quality=quality,
                          minimum_word_pairs=minimum, same_pair_pairs=union, same_read_pairs=same,
                          exact_pairs=exact, baseline_labels=json.dumps(group_labels[gid],separators=(',',':')),
                          count_basis=COUNT_BASIS)
            evidence[pi, sid, k, quality] = record; native.append(record)
            group_field_checks += 3
        group_rows += 1
    require(group_rows == len(groups)*len(accessions), 'Incomplete group/accession matrix')
    require(group_seen == {(pi,gid) for pi in accessions for gid in groups}, 'Missing group/accession count')
    del word_arrays
    monotonic_checks = 0
    for pi in accessions:
        for sid, items in state_groups.items():
            for q in qualities:
                previous = None
                for k,gid in sorted(items):
                    current = evidence[pi,sid,k,q]
                    if previous is not None:
                        for field in METHOD_FIELDS.values():
                            require(current[field] <= previous[field], 'Support increased as k increased: '+field)
                            monotonic_checks += 1
                        require(current['exact_pairs'] == previous['exact_pairs'], 'Exact count changed with baseline')
                    previous = current

    cohort = {}
    for pi in accessions:
        for sid in states:
            row = evidence[pi,sid,k_values[0],primary_q]
            cohort[pi,sid] = ('CONCORDANT_POSITIVE' if row['exact_pairs'] >= primary_t else
                              'DISCORDANT' if row['minimum_word_pairs'] >= primary_t else 'BELOW_THRESHOLD')
    selected_keys = {key for key,value in cohort.items() if value == 'DISCORDANT'}
    if 'expected_selected_cells' in protocol:
        require(len(selected_keys) == protocol['expected_selected_cells'], 'Selected discrepancy cell count differs')
    if 'expected_selected_states' in protocol:
        require(len({sid for pi,sid in selected_keys}) == protocol['expected_selected_states'], 'Selected discrepancy state count differs')
    common_states = {sid for sid in states if all((sid,'k'+str(k)) in eligible for k in k_values)}
    strata = ('ALL_CATALOG', common_label, 'K31_PRIMARY_CONCORDANT_POSITIVE',
              'K31_PRIMARY_DISCORDANT', 'K31_PRIMARY_BELOW_THRESHOLD', 'NONSELECTED_FOR_CODE20')
    def memberships(pi,sid):
        values = ['ALL_CATALOG', 'K31_PRIMARY_'+cohort[pi,sid]]
        if sid in common_states: values.append(common_label)
        if cohort[pi,sid] != 'DISCORDANT': values.append('NONSELECTED_FOR_CODE20')
        return values
    stratum_totals = Counter(s for pi in accessions for sid in states for s in memberships(pi,sid))
    coverage = []
    for sid,state in sorted(states.items()):
        for label in baseline_names:
            available = (sid,label) in eligible
            gid = eligible.get((sid,label),'')
            actual_k = len(state['sequence'])-1 if label == 'near_full' else int(label[1:])
            coverage.append(dict(state_id=sid,pair_id=state['pair_id'],baseline=label,k=actual_k,
                span_length=len(state['sequence']),eligible=available,group_id=gid,
                reason='MEASURED_NATIVE_GROUP' if available else 'K_NOT_BELOW_TARGET_SPAN_UNAVAILABLE',
                aliases_native_group_with_other_baseline=available and len(group_labels[gid])>1,
                native_baseline_labels=json.dumps(group_labels[gid],separators=(',',':')) if available else '[]',
                accession_cells=len(accessions),common_fixed_grid=sid in common_states))

    calls = []; lookup = {}; residuals = []; aggregates = defaultdict(Counter); equivalence = defaultdict(Counter)
    for row in native:
        pi,sid,q = row['pi_accession'],row['state_id'],row['quality']; member = memberships(pi,sid)
        for baseline in group_labels[row['group_id']]:
            for scope in member:
                for method,field in METHOD_FIELDS.items():
                    if method == 'exact':continue
                    eq = equivalence[baseline,q,scope,method]
                    eq['eligible_cells'] += 1;eq['equal_count_cells'] += row[field] == row['exact_pairs']
                    eq['surplus_count_cells'] += row[field] > row['exact_pairs']
                    eq['support_surplus_total'] += row[field]-row['exact_pairs']
                    for threshold in thresholds:
                        eq['call_discordant_at_support_'+str(threshold)] += (row[field]>=threshold) != (row['exact_pairs']>=threshold)
            for threshold in thresholds:
                rec = {key:value for key,value in row.items() if key!='baseline_labels'}
                rec.update(baseline=baseline,minimum_support=threshold,
                           cohort_at_primary_k31=cohort[pi,sid], selected_for_code20=cohort[pi,sid]=='DISCORDANT',
                           common_fixed_grid=sid in common_states,
                           evidence_class=classification(row,threshold),
                           **{method+'_pass':row[field]>=threshold for method,field in METHOD_FIELDS.items()})
                rec['pooled_without_exact'] = rec['pooled_pass'] and not rec['exact_pass']
                rec['same_pair_without_exact'] = rec['same_pair_pass'] and not rec['exact_pass']
                rec['same_read_without_exact'] = rec['same_read_pass'] and not rec['exact_pass']
                calls.append(rec);lookup[pi,sid,baseline,q,threshold] = rec
                if rec['pooled_without_exact']: residuals.append(rec.copy())
                for scope in member:
                    agg=aggregates[baseline,q,threshold,scope];agg['eligible_cells']+=1
                    for method in METHOD_FIELDS:agg[method+'_supported_cells']+=rec[method+'_pass']
                    for method in ('pooled','same_pair','same_read'):agg[method+'_without_exact_cells']+=rec[method+'_without_exact']
                    agg['class_'+rec['evidence_class']]+=1

    comparison_rows = 0; comparison_checks = 0; comparison_seen = set()
    old_to_new = {'minimum_word_pair_support':'minimum_word_pairs','same_read_complete_pairs':'same_read_pairs',
                  'mate_union_complete_pairs':'same_pair_pairs','exact_span_pairs':'exact_pairs'}
    old_flags = {'pooled_supported':'pooled_pass','same_read_supported':'same_read_pass',
                 'mate_union_supported':'same_pair_pass','exact_span_supported':'exact_pass','pooled_without_exact':'pooled_without_exact'}
    for old in read_table(source19/'support_comparisons.csv.gz'):
        key=(old['pi_accession'],old['state_id'],old['baseline'],old['quality'],as_int(old['minimum_supporting_pairs']))
        require(key not in comparison_seen,'Duplicate source comparison');comparison_seen.add(key)
        require(key in lookup,'Unexpected source comparison');new=lookup[key]
        for field in ('group_id','pair_id'):require(old[field]==new[field],'Source comparison metadata differs')
        for field in ('k','span_length'):require(as_int(old[field])==new[field],'Source comparison length differs')
        for oldfield,newfield in old_to_new.items():
            require(as_int(old[oldfield])==new[newfield],'Recomputed source comparison count differs: '+oldfield);comparison_checks+=1
        for oldfield,newfield in old_flags.items():
            require(as_bool(old[oldfield])==new[newfield],'Recomputed source comparison call differs: '+oldfield);comparison_checks+=1
        relation={'EXACT_SPAN_SUPPORTED':'CONCORDANT_POSITIVE','NO_COMPLETE_WORD_SUPPORT':'ALL_BELOW_THRESHOLD'}.get(old['evidence_relation'],old['evidence_relation'])
        require(relation==new['evidence_class'],'Source comparison class differs')
        for quality in ('raw','q20'):
            require(as_int(old[quality+'_exact_span_pairs'])==evidence[key[0],key[1],new['k'],quality]['exact_pairs'],'Source raw/Q20 span count differs')
        comparison_checks+=3;comparison_rows+=1
    require(comparison_seen==set(lookup),'Missing source comparison rows')
    require(comparison_rows==len(eligible)*len(accessions)*len(qualities)*len(thresholds),'Comparison denominator differs')

    summaries=[];eq_rows=[]
    for baseline in baseline_names:
        for q in qualities:
            for scope in strata:
                for threshold in thresholds:
                    agg=aggregates[baseline,q,threshold,scope];n=agg['eligible_cells']
                    rec=dict(baseline=baseline,quality=q,minimum_support=threshold,scope=scope,
                        total_scope_cells=stratum_totals[scope],eligible_cells=n,
                        unavailable_cells=stratum_totals[scope]-n,count_basis=COUNT_BASIS)
                    for method in METHOD_FIELDS:rec[method+'_supported_cells']=agg[method+'_supported_cells']
                    for method in ('pooled','same_pair','same_read'):
                        num=agg[method+'_without_exact_cells'];rec[method+'_without_exact_cells']=num
                        rec[method+'_discordance_fraction_of_eligible']=num/n if n else None
                    for name in CLASSES:rec['class_'+name]=agg['class_'+name]
                    require(sum(rec['class_'+name] for name in CLASSES)==n,'Class partition denominator failed')
                    summaries.append(rec)
                for method in ('pooled','same_pair','same_read'):
                    eq=equivalence[baseline,q,scope,method]
                    rec=dict(baseline=baseline,quality=q,scope=scope,method=method,count_basis=COUNT_BASIS,
                        total_scope_cells=stratum_totals[scope],eligible_cells=eq['eligible_cells'],
                        unavailable_cells=stratum_totals[scope]-eq['eligible_cells'],equal_count_cells=eq['equal_count_cells'],
                        surplus_count_cells=eq['surplus_count_cells'],support_surplus_total=eq['support_surplus_total'])
                    for threshold in thresholds:rec['call_discordant_at_support_'+str(threshold)]=eq['call_discordant_at_support_'+str(threshold)]
                    require(rec['equal_count_cells']+rec['surplus_count_cells']==rec['eligible_cells'],'Count equality denominator failed')
                    eq_rows.append(rec)

    simple = defaultdict(Counter)
    for pi in accessions:
        for sid in sorted(common_states):
            for q in qualities:
                for threshold in thresholds:
                    exact=lookup[pi,sid,base_label,q,threshold]['exact_pass']
                    rules={method+'_'+baseline:lookup[pi,sid,baseline,q,threshold][method+'_pass']
                           for baseline in baseline_names for method in ('pooled','same_pair','same_read')}
                    for method in ('same_pair','same_read'):
                        name='pooled_'+long_label+'_AND_'+method+'_'+base_label
                        rules[name]=lookup[pi,sid,long_label,q,threshold]['pooled_pass'] and lookup[pi,sid,base_label,q,threshold][method+'_pass']
                    for name,prediction in rules.items():
                        agg=simple[name,q,threshold];agg['eligible_cells']+=1;agg['rule_supported_cells']+=prediction;agg['exact_supported_cells']+=exact
                        agg['both_supported_cells' if prediction and exact else 'rule_only_cells' if prediction else 'exact_only_cells' if exact else 'both_below_cells']+=1
    simple_rows=[]
    for (rule,q,threshold),agg in sorted(simple.items()):
        rec=dict(rule=rule,quality=q,minimum_support=threshold,scope=common_label,count_basis=COUNT_BASIS,
                 rule_scope='CELL_LEVEL_CONJUNCTION_OF_MARGINAL_CALLS' if '_AND_' in rule else 'SINGLE_MARGINAL_THRESHOLD_CALL',
                 **{key:agg[key] for key in ('eligible_cells','rule_supported_cells','exact_supported_cells','both_supported_cells','rule_only_cells','exact_only_cells','both_below_cells')})
        require(rec['eligible_cells']==len(common_states)*len(accessions),'Simple-rule common denominator failed')
        require(sum(rec[x] for x in ('both_supported_cells','rule_only_cells','exact_only_cells','both_below_cells'))==rec['eligible_cells'],'Simple-rule contingency failed')
        simple_rows.append(rec)

    primary_class=Counter(r['evidence_class'] for r in calls if r['baseline']==base_label and r['quality']==primary_q and r['minimum_support']==primary_t)
    return dict(tables=dict(native_evidence=native,baseline_coverage=coverage,full_catalog_calls=calls,
                           baseline_summary=summaries,count_equivalence=eq_rows,residual_cases=residuals,
                           simple_rule_comparison=simple_rows),
        summary=dict(status='CODE26_CATALOG_AUDIT_COMPLETE',catalog_states=len(states),catalog_pairs=len({s['pair_id'] for s in states.values()}),
            accessions=len(accessions),accession_ids=accessions,state_accession_cells=len(states)*len(accessions),native_groups=len(groups),
            query_count=len(queries),word_count_rows=word_rows,group_count_rows=group_rows,native_evidence_rows=len(native),
            comparison_rows=len(calls),source_comparison_rows_reconciled=comparison_rows,source_comparison_field_checks=comparison_checks,
            word_counter_invariant_checks=word_field_checks,group_counter_invariant_checks=group_field_checks,k_monotonicity_checks=monotonic_checks,
            baseline_eligible_states={label:sum((sid,label) in eligible for sid in states) for label in baseline_names},
            common_fixed_grid_states=len(common_states),common_fixed_grid_cells=len(common_states)*len(accessions),
            primary_k31_cohort_counts=dict(Counter(cohort.values())),primary_k31_evidence_classes=dict(primary_class),
            count_basis=COUNT_BASIS,full_catalog_signature_deduplication=False,genomic_truth_known=False,
            baseline_alias_rows_are_independent=False,strata_form_one_partition=False,
            interpretation='Agreement and discordance of evidence criteria within a selected frozen catalog and same-study holdout. These are not biological false-positive rates or sensitivity.'))
