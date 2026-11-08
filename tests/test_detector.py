"""Independent substring oracle and adversarial detector semantics."""
import itertools
import random
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from kmer_evidence.detector import Detector, LiteralMatcher, canonical, rc


def rec(a, b='N', qa=None, qb=None):
    return dict(r1_sequence=a, r2_sequence=b,
                r1_quality='I' * len(a) if qa is None else qa,
                r2_quality='I' * len(b) if qb is None else qb)


def brute_hits(pattern, sequence, quality):
    """No finite-state implementation: inspect each literal substring."""
    variants = {pattern, rc(pattern)}
    return [(i, i + len(pattern), all(ord(q) >= 53 for q in quality[i:i+len(pattern)]))
            for i in range(len(sequence)-len(pattern)+1)
            if sequence[i:i+len(pattern)] in variants]


def naive_evaluate(target, reads, k_values):
    out = []
    for k in sorted({k for k in k_values if k <= len(target)} | {len(target)}):
        words = {canonical(target[i:i+k]) for i in range(len(target)-k+1)}
        for quality in ('raw', 'q20'):
            counts = {word: [] for word in words}
            criteria = {name: [] for name in ('same_pair', 'same_read', 'exact')}
            for row in reads:
                sig = (row['r1_sequence'], row['r2_sequence'])
                seen = []
                exact = False
                for mate in (1, 2):
                    seq, qual = row['r%d_sequence' % mate], row['r%d_quality' % mate]
                    seen.append({word for word in words
                                 if any(high or quality == 'raw' for _, _, high in brute_hits(word, seq, qual))})
                    exact |= any(high or quality == 'raw' for _, _, high in brute_hits(canonical(target), seq, qual))
                for word in seen[0] | seen[1]:
                    counts[word].append(sig)
                if words <= seen[0] | seen[1]:
                    criteria['same_pair'].append(sig)
                if words <= seen[0] or words <= seen[1]:
                    criteria['same_read'].append(sig)
                if exact:
                    criteria['exact'].append(sig)
            values = dict(minimum_word_pairs=min(map(len, counts.values())),
                          minimum_word_signatures=min(len(set(v)) for v in counts.values()))
            for name, evidence in criteria.items():
                values[name + '_pairs'] = len(evidence)
                values[name + '_signatures'] = len(set(evidence))
            out.append(dict(k=k, span_length=len(target), quality=quality,
                            distinct_canonical_words=len(words), **values))
    return out


class MatcherTests(unittest.TestCase):
    def test_exhaustive_short_ambiguous_reverse_palindromes_and_overlaps(self):
        patterns = {canonical(''.join(x)) for k in range(1, 4) for x in itertools.product('ACGT', repeat=k)}
        matcher = LiteralMatcher(patterns)
        for chars in itertools.product('ACNT', repeat=4):
            seq, qual = ''.join(chars), 'I5!4'
            expected = sorted((pid, start, end, high) for pid, p in enumerate(matcher.patterns)
                              for start, end, high in brute_hits(p, seq, qual))
            self.assertEqual(sorted(matcher.scan(seq, qual)), expected)

    def test_invalid_inputs(self):
        for patterns in ([], ['T'], ['N'], ['']):
            with self.assertRaises(ValueError):
                LiteralMatcher(patterns)
        matcher = LiteralMatcher(['A'])
        for seq, qual in (('', ''), ('a', 'I'), ('X', 'I'), ('A', ''), ('A', ' ')):
            with self.assertRaises(ValueError):
                matcher.scan(seq, qual)


class DetectorTests(unittest.TestCase):
    def view(self, rows, target='ACGAC', k=3):
        return {(r['k'], r['quality']): r for r in Detector(target, [k]).evaluate(rows)}

    def test_word_set_loses_order_spacing_and_multiplicity(self):
        row = self.view([rec('GACNACGNCGA')])[3, 'q20']
        self.assertEqual((row['minimum_word_signatures'], row['same_read_signatures'], row['exact_signatures']), (1, 1, 0))
        row = self.view([rec('AAA')], 'AAAAA', 3)[3, 'q20']
        self.assertEqual((row['distinct_canonical_words'], row['same_read_signatures'], row['exact_signatures']), (1, 1, 0))

    def test_mate_union_and_cross_pair_pooling_are_distinct(self):
        row = self.view([rec('ACGA', 'GAC')])[3, 'q20']
        self.assertEqual((row['minimum_word_signatures'], row['same_pair_signatures'], row['same_read_signatures']), (1, 1, 0))
        row = self.view([rec('ACGA'), rec('GAC')])[3, 'q20']
        self.assertEqual((row['minimum_word_signatures'], row['same_pair_signatures']), (1, 0))

    def test_duplicate_quality_mosaics_not_spliced_for_same_pair_or_read(self):
        rows = [rec('ACGAC', qa='III!!'), rec('ACGAC', qa='!IIII')]
        result = self.view(rows)
        row = result[3, 'q20']
        self.assertEqual((row['minimum_word_signatures'], row['same_pair_signatures'], row['same_read_signatures'], row['exact_signatures']), (1, 0, 0, 0))
        self.assertEqual((result[5, 'raw']['exact_pairs'], result[5, 'raw']['exact_signatures']), (2, 1))

    def test_ordered_signature_and_q20_boundary(self):
        result = self.view([rec('ACGAC', 'N', '55555'), rec('N', 'ACGAC', qb='55555'),
                            rec('ACGAC', 'N', '55455')])
        self.assertEqual((result[5, 'raw']['exact_pairs'], result[5, 'raw']['exact_signatures']), (3, 2))
        self.assertEqual((result[5, 'q20']['exact_pairs'], result[5, 'q20']['exact_signatures']), (2, 2))

    def test_reusable_detector_has_no_cross_cell_evidence(self):
        detector = Detector('ACGAC', [3])
        first = detector.evaluate([rec('ACGAC')])
        zero = detector.evaluate(iter(()))
        self.assertTrue(all(r['minimum_word_pairs'] == r['exact_signatures'] == 0 for r in zero))
        self.assertEqual(first, detector.evaluate(iter([rec('ACGAC')])))

    def test_full_length_exact_equivalence_and_reverse_complement(self):
        rows = [rec('N' + rc('ACGAC') + 'N'), rec('ACGAC' * 3)]
        result = self.view(rows)
        for quality in ('raw', 'q20'):
            row = result[5, quality]
            for unit in ('pairs', 'signatures'):
                self.assertEqual({row[key + '_' + unit] for key in ('minimum_word', 'same_pair', 'same_read', 'exact')}, {2})

    def test_random_independent_oracle_with_repeated_and_low_complexity_targets(self):
        rng = random.Random(240025)
        targets = ['ACGAC', 'AAAAAA', 'ACGTACGTACGT', 'ATATATAT',
                   ''.join(rng.choices('ACGT', k=37)), ''.join(rng.choices('ACGT', k=82))]
        for target in targets:
            reads = []
            for i in range(35):
                a = rng.choice([target, rc(target), target[:len(target)//2], 'N'.join(target[j:j+3] for j in range(0, len(target), 3))])
                b = rng.choice([target, rc(target), 'NNN', ''.join(rng.choices('ACGTN', k=40))])
                qa = ''.join(rng.choices('I5!4', weights=(10, 10, 1, 1), k=len(a)))
                qb = ''.join(rng.choices('I5!4', weights=(10, 10, 1, 1), k=len(b)))
                reads.append(rec(a, b, qa, qb))
                if i % 7 == 0:
                    reads.append(rec(a, b))
            ks = [2, 3, 9, 31, 41, 81]
            self.assertEqual(Detector(target, ks).evaluate(reads), naive_evaluate(target, reads, ks))

    def test_invalid_targets_and_lengths(self):
        for target in ('', 'N', 'acg'):
            with self.assertRaises(ValueError):
                Detector(target)
        for ks in ([0], [-1], [3.0], [True]):
            with self.assertRaises(ValueError):
                Detector('ACGT', ks)
        self.assertEqual(Detector('ACGT', [8]).k_values, (4,))


if __name__ == '__main__':
    unittest.main()
