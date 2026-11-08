#!/usr/bin/env python3
"""Regenerate the complete, fixed synthetic benchmark from its target panel."""
import argparse
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from kmer_evidence.benchmark import context_tables, run
from kmer_evidence.io import content_digest, json_read, json_write, new_output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/synthetic_benchmark')
    args = parser.parse_args()
    out = new_output(args.output, inputs=[ROOT / 'config', ROOT / 'data', ROOT / 'src'])
    started = time.monotonic()
    json_write(out / 'status.json', {'status': 'RUNNING'})
    try:
        result = run({'run_dir': str(out), 'protocol_file': str(ROOT / 'config/benchmark_protocol.json'),
                      'panel_file': str(ROOT / 'config/target_panel.json')})
        context_tables(out / 'analysis')
        expected = json_read(ROOT / 'config/expected_benchmark_outputs.json')
        checks = []
        for name, expected_sha in expected.items():
            actual = content_digest(out / 'analysis' / name)
            if actual != expected_sha:
                raise ValueError('Synthetic replay content differs from the study result: ' + name)
            checks.append({'file': 'analysis/' + name, 'content_sha256': actual})
        json_write(out / 'reproduction.json', {'status': 'PASS', 'scope': 'Full fixed simulation grid', 'hash_basis': 'SHA-256 of decompressed content',
                                             'elapsed_seconds': time.monotonic() - started, 'files': checks})
        json_write(out / 'status.json', {'status': 'BENCHMARK_COMPLETE', 'simulation_cells': result['simulation_cells']})
        print('ANALYSIS_STATUS=BENCHMARK_COMPLETE')
        print('OUTPUT_DIRECTORY=' + str(out))
    except BaseException as exc:
        json_write(out / 'status.json', {'status': 'FAILED', 'error': str(exc)})
        raise


if __name__ == '__main__':
    main()
