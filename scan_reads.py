#!/usr/bin/env python3
"""Scan one existing paired gzip FASTQ library against a catalogue.

This optional command reads every input record. It neither downloads reads nor
performs the complete multi-accession study orchestration or signature census.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from kmer_evidence.io import digest, json_read, json_write, new_output


def prepare_panel(catalogue, output):
    catalogue, output = Path(catalogue), Path(output)
    with gzip.open(catalogue / 'queries.csv.gz', 'rt', newline='') as handle:
        queries = list(csv.DictReader(handle))
    with gzip.open(catalogue / 'groups.csv.gz', 'rt', newline='') as handle:
        groups = list(csv.DictReader(handle))
    with (output / 'queries.tsv').open('w', encoding='ascii', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['pattern_id', 'sequence'])
        writer.writerows((row['pattern_id'], row['sequence']) for row in queries)
    with (output / 'groups.tsv').open('w', encoding='ascii', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['group_id', 'state_id', 'k', 'word_ids', 'span_query_id'])
        for row in groups:
            writer.writerow([row['group_id'], row['state_id'], row['k'],
                             ','.join(json.loads(row['word_ids'])), row['span_query_id']])
    return {'queries': len(queries), 'groups': len(groups),
            'queries_sha256': digest(catalogue / 'queries.csv.gz'),
            'groups_sha256': digest(catalogue / 'groups.csv.gz')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r1', type=Path, required=True)
    parser.add_argument('--r2', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--catalogue', type=Path, default=ROOT / 'data/derived/catalog')
    parser.add_argument('--scanner', type=Path, default=ROOT / 'bin/panel_scan')
    parser.add_argument('--validation-pairs', type=int, default=1000, choices=range(1, 10001), metavar='1..10000')
    args = parser.parse_args()
    r1, r2, scanner = (path.resolve() for path in (args.r1, args.r2, args.scanner))
    if not all(path.is_file() for path in (r1, r2, scanner)):
        parser.error('Both input FASTQs and a compiled scanner must exist')
    if r1 == r2:
        parser.error('R1 and R2 must be different files')
    out = new_output(args.output, inputs=[args.catalogue])
    before = {str(path): (path.stat().st_size, path.stat().st_mtime_ns) for path in (r1, r2)}
    json_write(out / 'status.json', {'status': 'RUNNING'})
    try:
        panel = prepare_panel(args.catalogue, out)
        command = [str(scanner), '--queries', str(out / 'queries.tsv'), '--groups', str(out / 'groups.tsv'),
                   '--r1', str(r1), '--r2', str(r2), '--out', str(out / 'counts'),
                   '--validation-pairs', str(args.validation_pairs), '--examples-per-category', '3']
        with (out / 'scanner.stdout.log').open('w') as stdout, (out / 'scanner.stderr.log').open('w') as stderr:
            subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
        summary = json_read(out / 'counts/summary.json')
        if summary.get('status') != 'PANEL_SCAN_COMPLETE' or summary.get('full_eof') is not True or summary.get('oracle_mismatches') != 0:
            raise ValueError('Scanner did not report complete EOF and successful independent prefix validation')
        after = {str(path): (path.stat().st_size, path.stat().st_mtime_ns) for path in (r1, r2)}
        if before != after:
            raise ValueError('Input FASTQ metadata changed during the scan')
        json_write(out / 'scan_provenance.json', {'status': 'SCAN_COMPLETE', 'panel': panel,
                                                'scanner_sha256': digest(scanner), 'inputs': before,
                                                'scope': 'One paired library; read-pair observations, not a signature census',
                                                'summary': summary})
        json_write(out / 'status.json', {'status': 'SCAN_COMPLETE'})
        print('ANALYSIS_STATUS=SCAN_COMPLETE')
        print('OUTPUT_DIRECTORY=' + str(out))
    except BaseException as exc:
        json_write(out / 'status.json', {'status': 'FAILED', 'error': str(exc)})
        raise


if __name__ == '__main__':
    main()
