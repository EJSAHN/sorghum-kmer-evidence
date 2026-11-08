"""Selected signatures cannot become measurements for other catalog cells."""
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
from kmer_evidence.signatures import audit_overlay


def write_table(path,rows):
    with gzip.open(path,'wt',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


class OverlayTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.s24=self.root/'code24';self.s25=self.root/'code25'
        (self.s24/'analysis').mkdir(parents=True);(self.s25/'analysis').mkdir(parents=True)
        self.protocol=dict(k_values=[3,5],qualities=['raw','q20'],support_thresholds=[1,2,3],
            primary={'quality':'q20','minimum_support':2},expected_catalog_cells=4,
            expected_selected_cells=1,expected_selected_states=1)
        self.native=[]; self.calls=[]
        selected={('raw',3):(5,4,3,1),('q20',3):(4,3,2,0),('raw',5):(4,3,1,1),('q20',5):(2,1,0,0)}
        for pi,sid in [('A','S1'),('B','S1'),('A','S2'),('B','S2')]:
            for k in (3,5):
                for q in ('raw','q20'):
                    vals=selected[q,k] if (pi,sid)==('A','S1') else ((4,4,4,4) if q=='raw' else (3,3,3,3)) if (pi,sid)==('B','S1') else (1,1,1,0) if (pi,sid)==('A','S2') else (0,0,0,0)
                    self.native.append(dict(pi_accession=pi,state_id=sid,pair_id='P'+sid,group_id=sid+'K'+str(k),k=k,span_length=6,quality=q,
                        **dict(zip(('minimum_word_pairs','same_pair_pairs','same_read_pairs','exact_pairs'),vals))))
            self.calls.append(dict(pi_accession=pi,state_id=sid,baseline='k3',quality='q20',minimum_support=2,
                cohort_at_primary_k31='DISCORDANT' if (pi,sid)==('A','S1') else 'CONCORDANT_POSITIVE' if (pi,sid)==('B','S1') else 'BELOW_THRESHOLD'))
        self.rows=[]
        signatures={('raw',3):(4,3,2,1),('q20',3):(3,2,1,0),('raw',5):(3,2,1,1),('q20',5):(2,1,0,0)}
        for k in (3,5,6):
            for q in ('raw','q20'):
                pairs=selected[q,k] if k!=6 else (1,1,1,1) if q=='raw' else (0,0,0,0)
                sig=signatures[q,k] if k!=6 else pairs
                for t in (1,2,3):
                    row=dict(pi_accession='A',state_id='S1',pair_id='PS1',group_id='S1K3',k=str(k),span_length='6',quality=q,minimum_support=str(t))
                    for name,v,s in zip(('minimum_word','same_pair','same_read','exact'),pairs,sig):
                        row[name+'_pairs']=str(v);row[name+'_signatures']=str(s)
                    for name,s in zip(('pooled','same_pair','same_read','exact'),sig):row[name+'_pass_signatures']=str(s>=t)
                    self.rows.append(row)
        self.ledger=[dict(pi_accession='A',state_id='S1',pair_id='PS1',primary_evidence_label='DESCRIPTIVE_TEST_CONTEXT',placement_status='UNRESOLVED',mechanism_scope='Not causal truth',q20_minimum_k31_signatures='3',q20_same_pair_signatures='2',q20_same_read_signatures='1',code23_extra_repeat_example='False')]
        self.sim=[]
        for method,k in [('same_read',31),('exact_span',0)]:
            self.sim.append(dict(template_model='present_target',noise_model='clean',read_length='150',pair_count='8',truth_present='True',method=method,k=str(k),quality='q20',support_unit='signatures',minimum_support='2',trials='2',positive_calls='1',true_positive='1',false_negative='1',false_positive='0',true_negative='0',positive_call_rate='0.5',rate_interpretation='sensitivity'))
        (self.s25/'CODE25_summary.json').write_text(json.dumps({'status':'CODE25_COMPLETE','analysis':{'status':'CODE25_COMPLETE','synthetic_only':True,'simulation_cells':2}}))
        self.save()

    def tearDown(self):
        self.temp.cleanup()

    def save(self):
        write_table(self.s24/'analysis/all_baselines.csv.gz',self.rows)
        write_table(self.s24/'analysis/case_evidence_ledger.csv.gz',self.ledger)
        write_table(self.s25/'analysis/performance_by_stratum.csv.gz',self.sim)

    def audit(self):
        return audit_overlay(self.s24,self.s25,self.native,self.calls,self.protocol)

    def test_unknown_signatures_stay_blank_for_same_state_other_accession(self):
        out=self.audit();coverage={(r['pi_accession'],r['state_id']):r for r in out['tables']['signature_coverage']}
        self.assertTrue(coverage['A','S1']['signature_available'])
        self.assertEqual(coverage['A','S1']['minimum_word_signatures'],3)
        self.assertFalse(coverage['B','S1']['signature_available'])
        self.assertEqual(coverage['B','S1']['minimum_word_signatures'],'')
        self.assertEqual(coverage['B','S1']['signature_measurement_scope'],'NOT_MEASURED')
        self.assertEqual(out['summary']['signature_not_measured_cells'],3)
        self.assertEqual(out['summary']['raw_counter_reconciliations'],72)

    def test_support_units_keep_same_selected_denominator(self):
        rows=self.audit()['tables']['selected_signature_sensitivity']
        chosen={r['support_unit']:r for r in rows if (r['baseline'],r['quality'],r['minimum_support'])==('k3','q20',2)}
        self.assertEqual(chosen['read_pair_observations']['cells'],1)
        self.assertEqual(chosen['ordered_sequence_signatures']['cells'],1)
        self.assertEqual(chosen['read_pair_observations']['same_read_supported_cells'],1)
        self.assertEqual(chosen['ordered_sequence_signatures']['same_read_supported_cells'],0)

    def test_full_span_counts_use_exact_equivalence_and_base_group_identity(self):
        out=self.audit();rows=[r for r in out['tables']['selected_case_baselines'] if r['baseline']=='full_span']
        self.assertEqual(len(rows),6)
        for r in rows:
            self.assertEqual(len({r[k+'_pairs'] for k in ('minimum_word','same_pair','same_read','exact')}),1)

    def test_unselected_count_surplus_has_no_invented_signature_or_mechanism(self):
        rows=[r for r in self.audit()['tables']['incremental_case_context'] if (r['pi_accession'],r['state_id'])==('A','S2')]
        self.assertEqual(len(rows),4)
        self.assertTrue(all(r['signature_measurement_scope']=='NOT_MEASURED' and r['same_read_signatures']=='' for r in rows))
        self.assertTrue(all(r['code24_primary_evidence_label']=='NOT_ASSESSED_IN_CODE24' for r in rows))

    def test_missing_selected_combination_fails(self):
        self.rows.pop();self.save()
        with self.assertRaisesRegex(ValueError,'Missing or extra'):self.audit()

    def test_duplicate_selected_combination_fails(self):
        self.rows.append(dict(self.rows[0]));self.save()
        with self.assertRaisesRegex(ValueError,'Duplicate CODE24'):self.audit()

    def test_wrong_accession_selection_cannot_be_transferred(self):
        self.ledger[0]['pi_accession']='B';self.save()
        with self.assertRaisesRegex(ValueError,'ledger selection differs'):self.audit()

    def test_changed_raw_counter_fails_reconciliation(self):
        self.rows[0]['minimum_word_pairs']='6';self.save()
        with self.assertRaisesRegex(ValueError,'raw counters disagree'):self.audit()

    def test_signature_count_above_raw_fails(self):
        self.rows[0]['minimum_word_signatures']='6';self.save()
        with self.assertRaisesRegex(ValueError,'Signature count exceeds'):self.audit()

    def test_counter_cannot_change_between_thresholds(self):
        self.rows[1]['minimum_word_signatures']='3';self.save()
        with self.assertRaisesRegex(ValueError,'changes across support'):self.audit()

    def test_full_span_signature_counts_must_be_equivalent(self):
        r=next(r for r in self.rows if r['k']=='6' and r['quality']=='raw')
        r['same_read_signatures']='0';r['exact_signatures']='0';self.save()
        with self.assertRaisesRegex(ValueError,'Full-span signature equivalence'):self.audit()

    def test_native_missing_and_duplicate_identifiers_fail(self):
        self.native.append(dict(self.native[0]))
        with self.assertRaisesRegex(ValueError,'Duplicate native'):self.audit()
        self.native.pop();self.native.pop(1)
        with self.assertRaisesRegex(ValueError,'Missing catalog quality'):self.audit()

    def test_primary_stratum_mismatch_fails(self):
        self.calls[0]['cohort_at_primary_k31']='BELOW_THRESHOLD'
        with self.assertRaisesRegex(ValueError,'stratum assignment differs'):self.audit()

    def test_ledger_counter_mismatch_fails(self):
        self.ledger[0]['q20_same_read_signatures']='2';self.save()
        with self.assertRaisesRegex(ValueError,'ledger signature counter differs'):self.audit()

    def test_simulation_is_separate_aggregate_context(self):
        out=self.audit();self.assertEqual(out['summary']['simulation_stratum_comparisons'],1)
        self.assertFalse(out['summary']['synthetic_results_in_real_catalog_denominators'])
        self.assertTrue(all(not r['real_catalog_denominator_included'] for r in out['tables']['simulation_context']))
        self.sim[0].update(positive_calls='2',true_positive='2',false_negative='0',positive_call_rate='1.0');self.save()
        self.assertEqual(self.audit()['summary']['simulation_same_read31_vs_exact_positive_count_differences'],1)

    def test_simulation_missing_comparator_or_bad_rate_fails(self):
        self.sim=self.sim[1:];self.save()
        with self.assertRaisesRegex(ValueError,'Missing simulation comparator'):self.audit()
        self.setUp_simulation_again_for_rate()

    def setUp_simulation_again_for_rate(self):
        same=dict(self.sim[0],method='same_read',k='31');self.sim.append(same)
        self.sim[0]['positive_call_rate']='0.9';self.save()
        with self.assertRaisesRegex(ValueError,'Simulation rate differs'):self.audit()


if __name__=='__main__':unittest.main()
