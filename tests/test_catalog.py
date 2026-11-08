"""Independent literal-read fixture for the complete catalog counter audit."""
import copy
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

PACKAGE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(PACKAGE/'src'))
from kmer_evidence.catalog import audit_catalog, as_int


def reverse(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def canon(s):return min(s,reverse(s))


def write(path,rows):
    with gzip.open(path,'wt',newline='') as f:
        writer=csv.DictWriter(f,list(rows[0]));writer.writeheader();writer.writerows(rows)


def read(path):
    with gzip.open(path,'rt',newline='') as f:return list(csv.DictReader(f))


def hit_count(seq,qual,word,quality):
    literals={word,reverse(word)}
    return sum(seq[i:i+len(word)] in literals and (quality=='raw' or all(ord(q)>=53 for q in qual[i:i+len(word)]))
               for i in range(len(seq)-len(word)+1))


def pair(a,b,qa=None,qb=None):return [(a,qa or 'I'*len(a)),(b,qb or 'I'*len(b))]


def fixture(root):
    states=[dict(state_id='S0',pair_id='P0',sequence='ACGTTGCA',span_length=8),
            dict(state_id='S1',pair_id='P1',sequence='AGTCGC',span_length=6),
            dict(state_id='S2',pair_id='P2',sequence='CCATA',span_length=5)]
    samples={'A':[pair('ACGTTGCAAA','CCCC'),pair('ACGTTGG','TTGCAA'),pair('ACGTT','GGGG'),
                  pair('AGTCGC','TTTT'),pair('AGTCGCAA','AAAA')],
             'B':[pair('ACG','GGGG'),pair('AGTCGC','TTTT'),pair('AGTCGCAA','AAAA'),pair('AGTCGC','CCCC'),
                  pair('CCATA','GGGG'),pair('CCATAAA','AAAA','+++++++')]}
    groups=[];queries={}
    def query(kind,word):
        key=kind,canon(word)
        if key not in queries:
            queries[key]=dict(pattern_id='Q%04d'%len(queries),sequence=canon(word),k=len(word),query_kind=kind,catalog_owner_count=1)
        return queries[key]['pattern_id']
    eligibility=[]
    for s in states:
        seq=s['sequence'];L=len(seq)
        for k in sorted({x for x in (3,5) if x<L}|{L-1}):
            fixed=k in (3,5);near=k==L-1;labels=(['k'+str(k)] if fixed else [])+(['near_full'] if near else [])
            g=dict(group_id='G%d'%len(groups),state_id=s['state_id'],pair_id=s['pair_id'],k=k,span_length=L,
                   word_ids=json.dumps([query('word',w) for w in sorted({canon(seq[i:i+k]) for i in range(L-k+1)})]),
                   span_query_id=query('span',seq),is_fixed_grid=fixed,is_near_full=near,baseline_labels=json.dumps(labels))
            groups.append(g)
            for label in labels:eligibility.append(dict(state_id=s['state_id'],pair_id=s['pair_id'],group_id=g['group_id'],baseline=label,k=k,span_length=L,eligible=True))
    by_id={q['pattern_id']:q for q in queries.values()}
    word_rows=[];word_lookup={};group_rows=[];comp=[]
    for pi,observations in samples.items():
        for q in queries.values():
            row=dict(q,pi_accession=pi)
            for quality in ('raw','q20'):
                matches=[[hit_count(seq,qual,q['sequence'],quality) for seq,qual in mates] for mates in observations]
                row.update({quality+'_occurrences':sum(sum(m) for m in matches),
                            quality+'_read_support':sum(sum(n>0 for n in m) for m in matches),
                            quality+'_pair_support':sum(any(m) for m in matches)})
            word_rows.append(row);word_lookup[pi,q['pattern_id']]=row
        for g in groups:
            row=dict(g,pi_accession=pi);target=next(s['sequence'] for s in states if s['state_id']==g['state_id'])
            ws={by_id[w]['sequence'] for w in json.loads(g['word_ids'])}
            for quality in ('raw','q20'):
                sets=[[{w for w in ws if hit_count(seq,qual,w,quality)} for seq,qual in mates] for mates in observations]
                row.update({quality+'_same_read_complete_reads':sum(sum(ws<=h for h in m) for m in sets),
                            quality+'_same_read_complete_pairs':sum(any(ws<=h for h in m) for m in sets),
                            quality+'_mate_union_complete_pairs':sum(ws<=m[0]|m[1] for m in sets),
                            quality+'_exact_span_pairs':sum(any(hit_count(seq,qual,target,quality)>0 for seq,qual in mates) for mates in observations)})
            group_rows.append(row)
            for quality in ('raw','q20'):
                minimum=min(word_lookup[pi,w][quality+'_pair_support'] for w in json.loads(g['word_ids']))
                same=row[quality+'_same_read_complete_pairs'];mate=row[quality+'_mate_union_complete_pairs'];exact=row[quality+'_exact_span_pairs']
                for threshold in (1,2,3):
                    label=('EXACT_SPAN_SUPPORTED' if exact>=threshold else 'SAME_READ_WITHOUT_EXACT_SPAN' if same>=threshold else
                           'MATE_UNION_ONLY' if mate>=threshold else 'POOLED_WITHOUT_SAME_PAIR' if minimum>=threshold else 'NO_COMPLETE_WORD_SUPPORT')
                    for baseline in json.loads(g['baseline_labels']):
                        comp.append(dict(pi_accession=pi,state_id=g['state_id'],pair_id=g['pair_id'],group_id=g['group_id'],baseline=baseline,
                            k=g['k'],span_length=g['span_length'],quality=quality,minimum_supporting_pairs=threshold,
                            minimum_word_pair_support=minimum,same_read_complete_pairs=same,mate_union_complete_pairs=mate,exact_span_pairs=exact,
                            raw_exact_span_pairs=row['raw_exact_span_pairs'],q20_exact_span_pairs=row['q20_exact_span_pairs'],
                            pooled_supported=minimum>=threshold,same_read_supported=same>=threshold,mate_union_supported=mate>=threshold,
                            exact_span_supported=exact>=threshold,pooled_without_exact=minimum>=threshold and exact<threshold,evidence_relation=label))
    for filename,rows in [('target_sequences',states),('queries',list(queries.values())),('groups',groups),('word_counts',word_rows),
                          ('group_counts',group_rows),('baseline_eligibility',eligibility),('support_comparisons',comp)]:write(root/(filename+'.csv.gz'),rows)
    return dict(k_values=[3,5],qualities=['raw','q20'],support_thresholds=[1,2,3],primary={'quality':'q20','minimum_support':2},
                expected_catalog_states=3,expected_accessions=2,accessions=['A','B'],expected_native_groups=7,expected_queries=len(queries))


class CatalogAuditTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.protocol=fixture(self.root)
    def tearDown(self):self.tmp.cleanup()
    def audit(self):return audit_catalog(self.root,self.protocol)
    def change(self,name,edit):
        path=self.root/(name+'.csv.gz');rows=read(path);edit(rows);write(path,rows)

    def test_full_literal_fixture_and_zero_denominators(self):
        result=self.audit();s=result['summary'];t=result['tables']
        self.assertEqual(s['state_accession_cells'],6);self.assertEqual(s['group_count_rows'],14)
        self.assertEqual(s['native_evidence_rows'],28);self.assertEqual(s['comparison_rows'],96)
        self.assertEqual(s['word_count_rows'],self.protocol['expected_queries']*2)
        self.assertEqual(s['common_fixed_grid_cells'],4)
        self.assertFalse(s['full_catalog_signature_deduplication']);self.assertFalse(s['genomic_truth_known'])
        self.assertEqual(s['primary_k31_cohort_counts'],{'DISCORDANT':1,'CONCORDANT_POSITIVE':2,'BELOW_THRESHOLD':3})
        self.assertTrue(all(r['count_basis']=='READ_PAIR_OBSERVATIONS' for r in t['native_evidence']))

    def test_aliases_not_native_duplicates_and_unavailable_not_zero(self):
        t=self.audit()['tables'];cov=t['baseline_coverage']
        unavailable=[r for r in cov if not r['eligible']]
        self.assertEqual(len(unavailable),1);self.assertEqual(unavailable[0]['state_id'],'S2');self.assertEqual(unavailable[0]['baseline'],'k5')
        self.assertEqual(unavailable[0]['group_id'],'')
        self.assertFalse(any(r['state_id']=='S2' and r['baseline']=='k5' for r in t['full_catalog_calls']))
        aliases=[r for r in cov if r['state_id']=='S1' and r['baseline'] in ('k5','near_full')]
        self.assertEqual(len({r['group_id'] for r in aliases}),1)
        all_k5=[r for r in t['baseline_summary'] if r['baseline']=='k5' and r['scope']=='ALL_CATALOG']
        self.assertTrue(all(r['eligible_cells']==4 and r['unavailable_cells']==2 for r in all_k5))

    def test_count_difference_and_threshold_call_difference_are_separate(self):
        t=self.audit()['tables']
        native=next(r for r in t['native_evidence'] if (r['pi_accession'],r['state_id'],r['k'],r['quality'])==('A','S0',3,'q20'))
        self.assertEqual(native['same_pair_pairs'],2);self.assertEqual(native['exact_pairs'],1)
        one=next(r for r in t['full_catalog_calls'] if (r['pi_accession'],r['state_id'],r['baseline'],r['quality'],r['minimum_support'])==('A','S0','k3','q20',1))
        two=next(r for r in t['full_catalog_calls'] if (r['pi_accession'],r['state_id'],r['baseline'],r['quality'],r['minimum_support'])==('A','S0','k3','q20',2))
        self.assertTrue(one['same_pair_pass'] and one['exact_pass']);self.assertTrue(two['same_pair_without_exact'])
        eq=next(r for r in t['count_equivalence'] if (r['baseline'],r['quality'],r['scope'],r['method'])==('k3','q20','ALL_CATALOG','same_pair'))
        self.assertEqual(eq['surplus_count_cells'],1);self.assertEqual(eq['call_discordant_at_support_1'],0);self.assertEqual(eq['call_discordant_at_support_2'],1)

    def test_cell_level_conjunctions_use_common_population_only(self):
        rows=self.audit()['tables']['simple_rule_comparison']
        composite=[r for r in rows if '_AND_' in r['rule']]
        self.assertTrue(composite)
        self.assertTrue(all(r['eligible_cells']==4 for r in rows))
        self.assertTrue(all(r['rule_scope']=='CELL_LEVEL_CONJUNCTION_OF_MARGINAL_CALLS' for r in composite))
        self.assertTrue(all('joint_support' not in r for r in rows))
        self.assertTrue(all(r['both_supported_cells']+r['rule_only_cells']+r['exact_only_cells']+r['both_below_cells']==4 for r in rows))

    def test_missing_zero_word_row_fails(self):
        self.change('word_counts',lambda rows:rows.pop(next(i for i,r in enumerate(rows) if r['raw_pair_support']=='0')))
        with self.assertRaisesRegex(ValueError,'Missing query/accession'):self.audit()

    def test_duplicate_word_row_fails(self):
        self.change('word_counts',lambda rows:rows.append(dict(rows[0])))
        with self.assertRaisesRegex(ValueError,'Duplicate word-count'):self.audit()

    def test_missing_group_and_source_comparison_fail(self):
        self.change('group_counts',lambda rows:rows.pop())
        with self.assertRaisesRegex(ValueError,'Incomplete group'):self.audit()
        fixture(self.root)
        self.change('support_comparisons',lambda rows:rows.pop())
        with self.assertRaisesRegex(ValueError,'Missing source comparison'):self.audit()

    def test_group_target_word_membership_fails(self):
        def bad(rows):
            words=json.loads(rows[0]['word_ids']);rows[0]['word_ids']=json.dumps(words[:-1])
        self.change('groups',bad)
        with self.assertRaisesRegex(ValueError,'target word set'):self.audit()

    def test_source_comparison_mismatch_fails(self):
        self.change('support_comparisons',lambda rows:rows[0].update(minimum_word_pair_support='999'))
        with self.assertRaisesRegex(ValueError,'source comparison count'):self.audit()

    def test_noninteger_count_and_quality_violation_fail(self):
        for bad in ('1.0','-1',' 1','nan',True):
            with self.assertRaises(ValueError):as_int(bad)
        self.change('word_counts',lambda rows:rows[0].update(q20_occurrences='999'))
        with self.assertRaisesRegex(ValueError,'quality nesting'):self.audit()

    def test_accession_identity_and_expected_dimensions_are_pinned(self):
        self.protocol['accessions']=['A','C']
        with self.assertRaisesRegex(ValueError,'Accession identities'):self.audit()
        self.protocol['accessions']=['A','B'];self.protocol['expected_queries']+=1
        with self.assertRaisesRegex(ValueError,'Query count'):self.audit()

    def test_explicit_expected_denominators_and_selection_are_checked(self):
        self.protocol['expected_catalog_cells']=7
        with self.assertRaisesRegex(ValueError,'catalog cell denominator'):self.audit()
        self.protocol['expected_catalog_cells']=6;self.protocol['expected_comparison_rows']=97
        with self.assertRaisesRegex(ValueError,'comparison row denominator'):self.audit()
        self.protocol['expected_comparison_rows']=96;self.protocol['expected_selected_cells']=2
        with self.assertRaisesRegex(ValueError,'discrepancy cell count'):self.audit()
        self.protocol['expected_selected_cells']=1;self.protocol['expected_selected_states']=2
        with self.assertRaisesRegex(ValueError,'discrepancy state count'):self.audit()


if __name__=='__main__':unittest.main()
