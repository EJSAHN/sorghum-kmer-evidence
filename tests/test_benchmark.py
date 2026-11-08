"""Small complete benchmark; independent reconciliation of exported metrics."""
from collections import defaultdict
import copy
import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE / 'src'))
from kmer_evidence.benchmark import cells, run, verify_geometry
from kmer_evidence.simulation import build_templates, simulate_cell


def read_csv(path):
    with gzip.open(path, 'rt', newline='') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with gzip.open(path, 'wt', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def as_bool(value):
    if value not in ('True', 'False'):
        raise AssertionError('Invalid literal boolean in result: ' + str(value))
    return value == 'True'


class AnalysisIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.protocol = json.loads((PACKAGE / 'config/benchmark_protocol.json').read_text())
        panel = json.loads((PACKAGE / 'config/target_panel.json').read_text())
        first = next(s for s in panel['states'] if len(s['sequence']) == 82)
        second = next(s for s in panel['states'] if len(s['sequence']) == 97)
        cls.states = [first, second]
        panel['states'] = cls.states
        cls.protocol.update(replicates=2, expected_targets=2, expected_cells=288,
                            expected_read_pairs=5760)
        cls.panel_path = cls.root / 'panel.json'
        cls.protocol_path = cls.root / 'protocol.json'
        cls.panel_path.write_text(json.dumps(panel))
        cls.protocol_path.write_text(json.dumps(cls.protocol))
        cls.run_dir = cls.root / 'run'
        cls.run_dir.mkdir()
        cls.config = dict(run_dir=str(cls.run_dir), panel_file=str(cls.panel_path),
                          protocol_file=str(cls.protocol_path))
        cls.summary = run(cls.config)
        cls.out = cls.run_dir / 'analysis'
        cls.counts = read_csv(cls.out / 'counts_by_cell.csv.gz')
        cls.audits = read_csv(cls.out / 'generation_audit.csv.gz')
        cls.metrics = read_csv(cls.out / 'performance_by_stratum.csv.gz')
        cls.primary = read_csv(cls.out / 'primary_per_target.csv.gz')
        cls.paired = read_csv(cls.out / 'paired_method_comparisons.csv.gz')

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_complete_grid_and_observation_denominators(self):
        # Two targets, six templates, three noise modes, two lengths, two pair
        # counts and two replicates: 288 independent cell identifiers.
        self.assertEqual(self.summary['simulation_cells'], 288)
        self.assertEqual(self.summary['simulated_read_pairs'], 5760)
        self.assertEqual(len(self.audits), 288)
        self.assertEqual(len({a['cell_id'] for a in self.audits}), 288)
        self.assertEqual(sum(int(a['pair_count']) for a in self.audits), 5760)
        self.assertEqual(len(self.counts), 4032)  # six k values + full span, raw/Q20
        self.assertEqual(len(self.metrics), 16416)
        self.assertEqual(len(self.primary), 2736)
        self.assertEqual(len(self.paired), 3456)
        self.assertEqual(self.summary['exact_equivalence_checks'], 1152)
        self.assertEqual(self.summary['noise_geometry_checks'], 192)
        self.assertEqual(self.summary['status'], 'CODE25_COMPLETE')
        self.assertTrue(self.summary['synthetic_only'])
        self.assertFalse(self.summary['biological_validation'])

    def test_seeds_identifiers_and_noise_pairing_survive_target_order(self):
        forward = list(cells(self.protocol, self.states))
        backward = list(cells(self.protocol, list(reversed(self.states))))
        self.assertEqual(forward, backward)
        seed_groups = defaultdict(set)
        ids = set()
        for target, cell in forward:
            key = (target['state_id'], cell['template_model'], cell['read_length'],
                   cell['pair_count'], cell['replicate'])
            seed_groups[key].add(cell['seed'])
            self.assertNotIn(cell['cell_id'], ids)
            ids.add(cell['cell_id'])
        self.assertEqual(len(seed_groups), 96)
        self.assertTrue(all(len(seeds) == 1 for seeds in seed_groups.values()))
        self.assertEqual(len({next(iter(s)) for s in seed_groups.values()}), 96)

    def test_exported_replay_audits_regenerate_identical_observations(self):
        by_id = {r['cell_id']: r for r in self.audits}
        # All six templates, both target lengths and all noise settings. These
        # cells are selected by design labels, never by their detection result.
        selected = [(target, cell) for target, cell in cells(self.protocol, self.states)
                    if cell['pair_count'] == 8 and cell['read_length'] == 150 and cell['replicate'] == 0]
        self.assertEqual(len(selected), 36)
        models = {}
        for target, cell in selected:
            sid = target['state_id']
            if sid not in models:
                models[sid] = build_templates(target['sequence'], sid, self.protocol['master_seed'])
            reads, audit = simulate_cell(models[sid][cell['template_model']], cell['seed'],
                                         cell['pair_count'], cell['read_length'], cell['noise_model'])
            self.assertEqual(audit, json.loads(by_id[cell['cell_id']]['audit_json']))
            self.assertEqual(len(reads), 8)
            self.assertEqual(str(cell['seed']), by_id[cell['cell_id']]['seed'])

    def source_index(self):
        index = defaultdict(list)
        for row in self.counts:
            key = tuple(row[k] for k in ('template_model', 'noise_model', 'read_length',
                                        'pair_count', 'quality', 'k'))
            index[key].append(row)
        return index

    def test_every_summary_rate_and_primary_rate_matches_cell_evidence(self):
        index = self.source_index()
        columns = {'pooled': 'minimum_word', 'same_pair': 'same_pair',
                   'same_read': 'same_read', 'exact_span': 'exact'}
        for row in self.metrics + self.primary:
            k = '31' if row['method'] == 'exact_span' else row['k']
            key = tuple(row[x] for x in ('template_model', 'noise_model', 'read_length',
                                        'pair_count', 'quality')) + (k,)
            evidence = index[key]
            if 'state_id' in row:
                evidence = [x for x in evidence if x['state_id'] == row['state_id']]
            self.assertEqual(len(evidence), 2 if 'state_id' in row else 4)
            truth = as_bool(row['truth_present'])
            self.assertTrue(all(as_bool(x['truth_present']) == truth for x in evidence))
            positives = sum(int(x[columns[row['method']] + '_' + row['support_unit']])
                            >= int(row['minimum_support']) for x in evidence)
            self.assertEqual(int(row['trials']), len(evidence))
            self.assertEqual(int(row['positive_calls']), positives)
            self.assertAlmostEqual(float(row['positive_call_rate']), positives / len(evidence))
            self.assertEqual(int(row['true_positive']), positives if truth else 0)
            self.assertEqual(int(row['false_positive']), positives if not truth else 0)
            self.assertEqual(int(row['false_negative']), len(evidence)-positives if truth else 0)
            self.assertEqual(int(row['true_negative']), len(evidence)-positives if not truth else 0)
            self.assertEqual(row['rate_interpretation'], 'sensitivity' if truth else 'false_positive_rate')

    def test_every_paired_table_uses_identical_cells(self):
        index = self.source_index()
        mapping = {'pooled_k81': ('81', 'minimum_word'),
                   'same_pair_k31': ('31', 'same_pair'),
                   'same_read_k31': ('31', 'same_read'),
                   'exact_span': ('31', 'exact')}
        for row in self.paired:
            basekey = tuple(row[x] for x in ('template_model', 'noise_model', 'read_length',
                                            'pair_count', 'quality'))
            base = {r['cell_id']: r for r in index[basekey + ('31',)]}
            k, col = mapping[row['comparator']]
            other = {r['cell_id']: r for r in index[basekey + (k,)]}
            self.assertEqual(set(base), set(other))
            n = int(row['minimum_support']); unit = row['support_unit']
            outcomes = [(int(base[c]['minimum_word_' + unit]) >= n,
                         int(other[c][col + '_' + unit]) >= n) for c in base]
            self.assertEqual(int(row['trials']), 4)
            for field, expected in [('both_positive', (True, True)),
                                    ('pooled31_only', (True, False)),
                                    ('comparator_only', (False, True)),
                                    ('both_negative', (False, False))]:
                self.assertEqual(int(row[field]), outcomes.count(expected))

    def test_short_reads_cannot_support_long_contiguous_words(self):
        short = [r for r in self.counts if r['read_length'] == '60']
        self.assertTrue(short)
        for r in short:
            self.assertEqual(int(r['exact_pairs']), 0)
            if int(r['k']) > 60:
                self.assertEqual(int(r['minimum_word_pairs']), 0)
                self.assertEqual(int(r['same_pair_pairs']), 0)

    def test_geometry_mismatch_and_missing_audit_fail(self):
        damaged = copy.deepcopy(self.audits)
        payload = json.loads(damaged[-1]['audit_json'])
        payload['origin_sha256'] = '0' * 64
        damaged[-1]['audit_json'] = json.dumps(payload)
        bad_path = self.root / 'bad_geometry.csv.gz'
        write_csv(bad_path, damaged)
        with self.assertRaisesRegex(ValueError, 'geometry'):
            verify_geometry(bad_path, self.protocol)
        write_csv(bad_path, self.audits[:-1])
        with self.assertRaisesRegex(ValueError, 'Missing noise-paired'):
            verify_geometry(bad_path, self.protocol)

    def test_duplicate_noise_label_cannot_replace_missing_condition(self):
        damaged = copy.deepcopy(self.audits)
        target = damaged[-1]
        key_fields = ('state_id', 'template_model', 'read_length', 'pair_count', 'replicate')
        sibling = next(r for r in damaged
                       if all(r[k] == target[k] for k in key_fields)
                       and r['noise_model'] != target['noise_model'])
        # Preserve row count, origins and all numerical counts. A duplicate
        # label must still fail, rather than pretending all conditions exist.
        target['noise_model'] = sibling['noise_model']
        payload = json.loads(target['audit_json'])
        payload['noise_model'] = sibling['noise_model']
        target['audit_json'] = json.dumps(payload)
        bad_path = self.root / 'duplicate_noise.csv.gz'
        write_csv(bad_path, damaged)
        with self.assertRaisesRegex(ValueError, 'Duplicate noise condition'):
            verify_geometry(bad_path, self.protocol)

    def test_rerun_does_not_silently_overwrite_counts(self):
        with self.assertRaisesRegex(ValueError, 'output exists'):
            run(self.config)


if __name__ == '__main__':
    unittest.main()
