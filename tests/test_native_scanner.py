"""Build the native scanner and compare tiny paired-read counters independently."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from kmer_evidence.detector import Detector, canonical


@unittest.skipUnless(shutil.which('g++'), 'Optional native checks need g++ and zlib development files')
class NativeScannerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.scanner = cls.root / 'panel_scan'
        subprocess.run([sys.executable, str(ROOT / 'build_scanner.py'), '--output', str(cls.scanner)],
                       check=True, capture_output=True, text=True)
        cls.catalogue = cls.root / 'catalogue'; cls.catalogue.mkdir()
        cls.target = 'ACGAC'
        words = sorted({canonical(cls.target[i:i+3]) for i in range(3)})
        queries = [{'pattern_id': 'Q%d' % i, 'sequence': sequence}
                   for i, sequence in enumerate(words + [canonical(cls.target)])]
        cls.words = words
        groups = [{'group_id': 'G0', 'state_id': 'S0', 'k': 3,
                   'word_ids': json.dumps([q['pattern_id'] for q in queries[:-1]]),
                   'span_query_id': queries[-1]['pattern_id']}]
        for name, rows in [('queries', queries), ('groups', groups)]:
            with gzip.open(cls.catalogue / (name + '.csv.gz'), 'wt', newline='') as handle:
                writer = csv.DictWriter(handle, list(rows[0])); writer.writeheader(); writer.writerows(rows)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def records(self, name, rows, bad_id=False):
        paths = []
        for mate in (1, 2):
            path = self.root / ('%s_%d.fastq.gz' % (name, mate))
            with gzip.open(path, 'wt', newline='') as handle:
                for i, record in enumerate(rows):
                    identifier = i + (1 if bad_id and mate == 2 else 0)
                    handle.write('@r%d/%d\n%s\n+\n%s\n' % (identifier, mate,
                        record['r%d_sequence' % mate], record['r%d_quality' % mate]))
            paths.append(path)
        return paths

    def command(self, paths, output):
        return [sys.executable, str(ROOT / 'scan_reads.py'), '--r1', str(paths[0]), '--r2', str(paths[1]),
                '--scanner', str(self.scanner), '--catalogue', str(self.catalogue), '--output', str(output)]

    def test_streamed_counts_equal_python_detector(self):
        rows = []
        for a, b, qa in [('ACGAC', 'NNNNN', 'IIIII'), ('ACGA', 'GAC', 'IIII'),
                         ('GACNACGNCGA', 'NN', 'I'*11), ('ACGAC', 'NN', 'II!II'),
                         ('GTCGT', 'NN', 'IIIII'), ('ACG', 'NN', 'III')]:
            rows.append({'r1_sequence': a, 'r2_sequence': b, 'r1_quality': qa, 'r2_quality': 'I'*len(b)})
        paths = self.records('matched', rows); out = self.root / 'matched_output'
        result = subprocess.run(self.command(paths, out), text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('ANALYSIS_STATUS=SCAN_COMPLETE', result.stdout)
        summary = json.loads((out / 'counts/summary.json').read_text())
        self.assertEqual(summary['pairs'], len(rows)); self.assertEqual(summary['oracle_mismatches'], 0)
        with (out / 'counts/group_counts.tsv').open() as handle:
            group = next(csv.DictReader(handle, delimiter='\t'))
        with (out / 'counts/word_counts.tsv').open() as handle:
            counts = list(csv.DictReader(handle, delimiter='\t'))
        predicted = {r['quality']: r for r in Detector(self.target, [3]).evaluate(rows) if r['k'] == 3}
        for quality in ('raw', 'q20'):
            native = [int(group[quality + suffix]) for suffix in
                      ('_same_read_complete_pairs', '_mate_union_complete_pairs', '_exact_span_pairs')]
            self.assertEqual(native, [predicted[quality][field] for field in ('same_read_pairs', 'same_pair_pairs', 'exact_pairs')])
            pooled = min(int(row[quality + '_pair_support']) for row in counts if row['sequence'] in self.words)
            self.assertEqual(pooled, predicted[quality]['minimum_word_pairs'])

    def test_mismatched_mates_fail_without_success_status(self):
        rows = [{'r1_sequence': 'ACGAC', 'r2_sequence': 'NNNNN', 'r1_quality': 'IIIII', 'r2_quality': 'IIIII'}]
        paths = self.records('mismatch', rows, bad_id=True); out = self.root / 'mismatch_output'
        result = subprocess.run(self.command(paths, out), text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(json.loads((out / 'status.json').read_text())['status'], 'FAILED')

    def test_corrupt_or_incomplete_gzip_is_not_accepted(self):
        rows = [{'r1_sequence': 'ACGAC', 'r2_sequence': 'NNNNN', 'r1_quality': 'IIIII', 'r2_quality': 'IIIII'}]
        for kind in ('truncated_trailer', 'bad_crc', 'trailing_bytes'):
            with self.subTest(kind=kind):
                paths = self.records(kind, rows)
                payload = bytearray(paths[0].read_bytes())
                if kind == 'truncated_trailer':
                    payload = payload[:-4]
                elif kind == 'bad_crc':
                    payload[-8] ^= 1
                else:
                    payload.extend(b'not another gzip member')
                paths[0].write_bytes(payload)
                out = self.root / (kind + '_output')
                result = subprocess.run(self.command(paths, out), text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(json.loads((out / 'status.json').read_text())['status'], 'FAILED')

    def test_unequal_mate_record_counts_fail(self):
        rows = [{'r1_sequence': 'ACGAC', 'r2_sequence': 'NNNNN', 'r1_quality': 'IIIII', 'r2_quality': 'IIIII'}]
        paths = self.records('unequal_count', rows)
        with gzip.open(paths[1], 'at', newline='') as handle:
            handle.write('@r1/2\nNNNNN\n+\nIIIII\n')
        out = self.root / 'unequal_count_output'
        result = subprocess.run(self.command(paths, out), text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(json.loads((out / 'status.json').read_text())['status'], 'FAILED')


if __name__ == '__main__':
    unittest.main()
