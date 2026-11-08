import copy
import json
from pathlib import Path
import random
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from kmer_evidence.simulation import (build_templates, simulate_cell, stable_seed, rc,
                               words, validate_template_record, _noise)


class SimulationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = random.Random(147)
        cls.target = ''.join(rng.choice('ACGT') for _ in range(82))
        cls.templates = build_templates(cls.target, 'T1', 20261107)

    def test_stable_seed_is_typed_and_not_ambiguous(self):
        self.assertEqual(stable_seed('x', 4), stable_seed('x', 4))
        self.assertNotEqual(stable_seed('x', 4), stable_seed('x', '4'))
        self.assertNotEqual(stable_seed('ab', 'c'), stable_seed('a', 'bc'))

    def test_all_model_truth_on_entire_templates_both_strands(self):
        self.assertEqual(len(self.templates), 6)
        for model, record in self.templates.items():
            observed = validate_template_record(record)
            self.assertEqual(any(observed), model.startswith('present'))
            for template, expected in zip(record['templates'], observed):
                self.assertEqual(len(template), 600)
                self.assertEqual(self.target in template or rc(self.target) in template, expected)

    def test_record_truth_rejects_reverse_complement_at_unexpected_location(self):
        record = copy.deepcopy(self.templates['absent_background'])
        record['templates'][0] = rc(self.target)+record['templates'][0][len(self.target):]
        with self.assertRaisesRegex(ValueError, 'truth metadata'):
            validate_template_record(record)

    def test_control_complementation_invariants(self):
        record = self.templates['absent_complement_internal']
        a, b = record['construction']['core_sequences']
        self.assertTrue(words(self.target, 31) <= words(a, 31) | words(b, 31))
        record = self.templates['absent_complement_endpoints']
        a, b = record['construction']['core_sequences']
        for k in range(1, len(self.target)):
            self.assertTrue(words(self.target, k) <= words(a, k) | words(b, k))
        self.assertNotIn(self.target, (a, b))

    def test_full_build_and_reads_are_reproducible(self):
        self.assertEqual(self.templates, build_templates(self.target, 'T1', 20261107))
        args = (self.templates['present_target'], 981, 32, 150, 'error_1pct')
        self.assertEqual(simulate_cell(*args), simulate_cell(*args))

    def test_fixed_sampling_center_for_even_and_odd_target_lengths(self):
        for record in self.templates.values():
            self.assertEqual(record['sampling_center0'], 300)
        odd = build_templates(self.target+'ACGTTGAGTCGATAC', 'odd97', 20261107)
        for record in odd.values():
            self.assertEqual(len(record['target_sequence']), 97)
            self.assertEqual(record['sampling_center0'], 300)
            reads, _ = simulate_cell(record, 813, 32, 60, 'clean')
            self.assertTrue(all(row['origin_fragment_start0'] <= 300 < row['origin_fragment_end0']
                                for row in reads))

    def test_noiseless_pair_geometry_and_random_strands(self):
        record = self.templates['present_target']
        reads, audit = simulate_cell(record, 193, 128, 150, 'clean')
        seen = set()
        for row in reads:
            start, end = row['origin_fragment_start0'], row['origin_fragment_end0']
            self.assertTrue(start <= record['sampling_center0'] < end)
            self.assertTrue(180 <= end-start <= 300)
            fragment = record['templates'][row['origin_template_index']][start:end]
            seen.add(row['origin_strand'])
            if row['origin_strand'] == '-':
                fragment = rc(fragment)
            self.assertEqual(row['r1_sequence'], fragment[:150])
            self.assertEqual(row['r2_sequence'], rc(fragment[-150:]))
            self.assertEqual(row['r1_quality'], 'I'*150)
        self.assertEqual(seen, {'+', '-'})
        self.assertEqual(audit['selected_target_template_pairs'], 128)
        self.assertEqual(audit['substitution_error_bases'], 0)

    def test_short_reads_do_not_change_positive_template_truth(self):
        record = self.templates['present_target']
        reads, audit = simulate_cell(record, 211, 32, 60, 'clean')
        self.assertTrue(audit['truth_present'])
        self.assertTrue(all(self.target not in row[m] and rc(self.target) not in row[m]
                            for row in reads for m in ('r1_sequence', 'r2_sequence')))

    def test_diluted_positive_can_sample_zero_target_templates(self):
        record = self.templates['present_diluted10']
        for seed in range(100):
            reads, audit = simulate_cell(record, seed, 1, 150, 'clean')
            if audit['selected_target_template_pairs'] == 0:
                self.assertTrue(audit['truth_present'])
                break
        else:
            self.fail('Deterministic test seeds unexpectedly all sampled target')

    def test_noise_does_not_change_origins_and_quality_patch_keeps_sequence(self):
        record = self.templates['present_target']
        clean, a = simulate_cell(record, 6001, 12, 150, 'clean')
        patch_reads, b = simulate_cell(record, 6001, 12, 150, 'quality_patch')
        errors, c = simulate_cell(record, 6001, 12, 150, 'error_1pct')
        self.assertEqual(a['origin_sha256'], b['origin_sha256'])
        self.assertEqual(a['origin_sha256'], c['origin_sha256'])
        self.assertEqual(b['low_quality_bases'], 12*2*10)
        self.assertEqual(b['substitution_error_bases'], 0)
        self.assertGreater(c['substitution_error_bases'], 0)
        for x, y in zip(clean, patch_reads):
            self.assertEqual(x['r1_sequence'], y['r1_sequence'])
            self.assertEqual(x['r2_sequence'], y['r2_sequence'])
        self.assertTrue(all(row['r1_quality'] == '5'*150 for row in errors))

    def test_error_created_target_is_not_filtered_or_relabelled(self):
        record = self.templates['absent_single_snv']
        # Force the noise function to emit a restored target. This adversarial
        # hook tests the essential no-post-observation-filtering contract.
        def restored(sequence, model, rng):
            if len(sequence) >= len(self.target):
                sequence = self.target + sequence[len(self.target):]
            return sequence, '5'*len(sequence), 1, 0
        with patch('kmer_evidence.simulation._noise', restored):
            reads, audit = simulate_cell(record, 788, 8, 150, 'error_1pct')
        self.assertFalse(audit['truth_present'])
        self.assertEqual(len(reads), 8)
        self.assertTrue(all(self.target in row['r1_sequence'] for row in reads))
        self.assertFalse(audit['post_observation_truth_filtering'])

    def test_invalid_inputs_fail_explicitly(self):
        with self.assertRaises(ValueError):
            build_templates('N'*82, 'bad', 1)
        with self.assertRaises(ValueError):
            simulate_cell(self.templates['present_target'], 4, 0, 150, 'clean')
        with self.assertRaises(ValueError):
            simulate_cell(self.templates['present_target'], 4, 1, 300, 'clean')


if __name__ == '__main__':
    unittest.main()
