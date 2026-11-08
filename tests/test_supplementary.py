"""Content hashing, the fraction rule, repeat profiles and simulation context tables."""
import csv
import gzip
import io
import json
from pathlib import Path
import random
import sys
import tempfile
import unittest
import zlib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from kmer_evidence.benchmark import context_tables
from kmer_evidence.io import content_digest, digest
from kmer_evidence.supplementary import (canonical, fraction_rule_table, longest_repeat,
                                         passes_fraction)


def brute_longest_repeat(sequence):
    best = 0
    for length in range(1, len(sequence)):
        words = [canonical(sequence[i:i + length]) for i in range(len(sequence) - length + 1)]
        if len(words) != len(set(words)):
            best = length
    return best


def write_gz(path, text, level):
    raw = zlib.compressobj(level, zlib.DEFLATED, 31)
    path.write_bytes(raw.compress(text.encode()) + raw.flush())


class ContentDigestTests(unittest.TestCase):
    def test_compression_level_changes_bytes_but_not_content_digest(self):
        text = ''.join('row%d,%d\n' % (i, i * 7 % 13) for i in range(20000))
        with tempfile.TemporaryDirectory() as tmp:
            fast, best = Path(tmp) / 'fast.csv.gz', Path(tmp) / 'best.csv.gz'
            write_gz(fast, text, 1)
            write_gz(best, text, 9)
            self.assertNotEqual(digest(fast), digest(best))
            self.assertEqual(content_digest(fast), content_digest(best))
            plain = Path(tmp) / 'plain.json'
            plain.write_text('{"a": 1}\n')
            self.assertEqual(content_digest(plain), digest(plain))


class FractionRuleTests(unittest.TestCase):
    def test_fraction_boundaries(self):
        supports = [2, 2, 0, 0, 0]
        self.assertTrue(passes_fraction(supports, 2, 0.4))
        self.assertFalse(passes_fraction(supports, 2, 0.5))
        self.assertTrue(passes_fraction([3, 3], 2, 1.0))
        self.assertFalse(passes_fraction([3, 1], 2, 1.0))

    def test_fraction_precision_is_not_rounded_to_tenths(self):
        self.assertFalse(passes_fraction([2, 2, 0, 0, 0], 2, 0.45))
        self.assertTrue(passes_fraction([2] * 9 + [0] * 11, 2, 0.45))
        self.assertFalse(passes_fraction([2] * 8 + [0] * 12, 2, 0.45))
        self.assertTrue(passes_fraction([2] * 29 + [0] * 71, 2, 0.29))
        self.assertFalse(passes_fraction([2] * 28 + [0] * 72, 2, 0.29))

    def test_invalid_fractions_are_rejected(self):
        for fraction in (-0.1, 1.1, float('nan'), float('inf'), float('-inf')):
            with self.subTest(fraction=fraction), self.assertRaises(ValueError):
                passes_fraction([2, 0], 2, fraction)

    def test_zero_fraction_and_empty_supports(self):
        self.assertTrue(passes_fraction([0, 0], 2, 0.0))
        with self.assertRaisesRegex(ValueError, 'at least one k-mer'):
            passes_fraction([], 2, 0.4)

    def test_full_fraction_reproduces_primary_pooled_comparison(self):
        protocol = json.loads((ROOT / 'config/catalog_protocol.json').read_text())
        rows = fraction_rule_table(ROOT / 'data/derived/catalog', protocol)
        primary = {(r['scope'], r['baseline'], r['required_kmer_fraction']): r for r in rows
                   if r['quality'] == 'q20' and r['minimum_support'] == 2}
        full = primary['ALL_CATALOG', 'k31', 1.0]
        self.assertEqual((full['eligible_cells'], full['fraction_rule_positive_cells'],
                          full['exact_positive_cells'], full['fraction_rule_positive_without_exact_cells']),
                         (5304, 2457, 2397, 60))
        self.assertEqual(primary['COMMON_K31_K81', 'k81', 1.0]['eligible_cells'], 5256)
        series = [primary['ALL_CATALOG', 'k31', f]['fraction_rule_positive_without_exact_cells']
                  for f in (0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)]
        self.assertEqual(series, sorted(series, reverse=True))


class RepeatProfileTests(unittest.TestCase):
    def test_longest_repeat_matches_brute_force_on_both_strands(self):
        rng = random.Random(7)
        for _ in range(200):
            seq = ''.join(rng.choice('ACGT') for _ in range(rng.randint(8, 40)))
            if rng.random() < .5:
                cut = rng.randint(0, len(seq) - 6)
                seq += seq[cut:cut + 6][::-1].translate(str.maketrans('ACGT', 'TGCA'))
            self.assertEqual(longest_repeat(seq), brute_longest_repeat(seq), seq)


class SimulationContextTests(unittest.TestCase):
    def test_ceiling_and_unit_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            audit_rows = [
                dict(cell_id='a', state_id='S', template_model='present_diluted10', noise_model='clean',
                     read_length=150, pair_count=8, replicate=0, seed=1, truth_present=True,
                     audit_json=json.dumps({'selected_target_template_pairs': n}))
                for n in (0, 1, 2, 5)]
            audit_rows.append(dict(audit_rows[0], template_model='present_target',
                                   audit_json=json.dumps({'selected_target_template_pairs': 8})))
            with gzip.open(out / 'generation_audit.csv.gz', 'wt', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0])); writer.writeheader()
                writer.writerows(audit_rows)
            count_row = dict(read_length=150, quality='q20', minimum_word_pairs=3, minimum_word_signatures=1,
                             same_pair_pairs=2, same_pair_signatures=2, same_read_pairs=0,
                             same_read_signatures=0, exact_pairs=0, exact_signatures=0)
            with gzip.open(out / 'counts_by_cell.csv.gz', 'wt', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(count_row)); writer.writeheader()
                writer.writerow(count_row)
            context_tables(out)
            with gzip.open(out / 'dilution_sampling_ceiling.csv.gz', 'rt') as handle:
                ceiling = list(csv.DictReader(handle))
            self.assertEqual(len(ceiling), 1)
            self.assertEqual((ceiling[0]['trials'], ceiling[0]['target_template_pairs_0'],
                              ceiling[0]['target_template_pairs_1'], ceiling[0]['target_template_pairs_at_least_2']),
                             ('4', '1', '1', '2'))
            with gzip.open(out / 'support_unit_agreement.csv.gz', 'rt') as handle:
                units = {(r['criterion'], r['minimum_support']): r['calls_differ_between_units']
                         for r in csv.DictReader(handle)}
            self.assertEqual(units['pooled', '1'], '0')
            self.assertEqual(units['pooled', '2'], '1')
            self.assertEqual(units['same_pair', '2'], '0')


if __name__ == '__main__':
    unittest.main()
