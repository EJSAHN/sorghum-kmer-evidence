"""Portable entry points preserve inputs and identify incomplete work."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from kmer_evidence.io import new_output, verify_manifest


class PublicInterfaceTests(unittest.TestCase):
    def test_input_manifest_rejects_changed_bytes_and_traversal(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'input.txt').write_text('original')
            row = {'file': 'input.txt', 'size_bytes': 8, 'sha256': hashlib.sha256(b'original').hexdigest()}
            manifest = root / 'input_manifest.json'
            manifest.write_text(json.dumps({'files': [row]}))
            self.assertEqual(len(verify_manifest(root, manifest)), 1)
            (root / 'input.txt').write_text('changed!')
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                verify_manifest(root, manifest)
            row['file'] = '../input.txt'
            manifest.write_text(json.dumps({'files': [row]}))
            with self.assertRaisesRegex(ValueError, 'Unsafe'):
                verify_manifest(root, manifest)

    def test_output_cannot_overwrite_inputs_or_nonempty_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'data'; source.mkdir()
            with self.assertRaisesRegex(ValueError, 'outside'):
                new_output(source / 'nested', [source])
            out = new_output(root / 'results', [source])
            (out / 'keep.txt').write_text('preserve')
            with self.assertRaisesRegex(ValueError, 'not empty'):
                new_output(out, [source])
            self.assertEqual((out / 'keep.txt').read_text(), 'preserve')

    def test_cli_help_from_unrelated_working_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            for entry in ('run_analysis.py', 'run_benchmark.py', 'scan_reads.py', 'build_scanner.py'):
                result = subprocess.run([sys.executable, str(ROOT / entry), '--help'], cwd=tmp,
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)


if __name__ == '__main__':
    unittest.main()
