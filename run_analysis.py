#!/usr/bin/env python3
"""Reproduce the catalogue evidence tables from the supplied derived counters."""
import argparse
import json
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from kmer_evidence import __version__
from kmer_evidence.catalog import audit_catalog
from kmer_evidence.headlines import primary_headlines
from kmer_evidence.signatures import audit_overlay
from kmer_evidence.io import content_digest, digest, json_read, json_write, new_output, verify_manifest, write_table
from kmer_evidence.supplementary import (fraction_rule_table, repeat_summary, selected_case_structure,
                                        target_repeat_profile)


def reproduce(input_root, output, protocol_file, expected_file):
    started = time.monotonic()
    records = verify_manifest(input_root, input_root / 'input_manifest.json')
    protocol = json_read(protocol_file)
    out = new_output(output, inputs=[input_root, ROOT / 'config', ROOT / 'src'])
    (out / 'tables').mkdir()
    (out / 'checks').mkdir()
    json_write(out / 'status.json', {'status': 'RUNNING'})
    try:
        print('PHASE=reconstruct_catalog', flush=True)
        catalog = audit_catalog(input_root / 'catalog', protocol)
        print('PHASE=reconcile_selected_signatures', flush=True)
        overlay = audit_overlay(input_root / 'signatures', input_root / 'simulation',
                                catalog['tables']['native_evidence'], catalog['tables']['full_catalog_calls'], protocol)
        tables = {**catalog['tables'], **overlay['tables']}
        headlines, denominators = primary_headlines(catalog, protocol)
        tables['primary_comparison'] = headlines
        print('PHASE=supplementary_catalog_analyses', flush=True)
        tables['fraction_rule'] = fraction_rule_table(input_root / 'catalog', protocol)
        profile = target_repeat_profile(input_root / 'catalog')
        tables['target_repeat_profile'] = profile
        case_structure = selected_case_structure(input_root / 'signatures/analysis/case_evidence_ledger.csv.gz')
        expected = json_read(expected_file)
        if set(expected) != {name + '.csv.gz' for name in tables}:
            raise ValueError('Expected output table identities differ from generated tables')
        checks = []
        for name, rows in tables.items():
            target = out / 'tables' / (name + '.csv.gz')
            write_table(target, rows)
            actual = content_digest(target)
            if actual != expected[target.name]:
                raise ValueError('Reproduced table content differs from the study result: ' + target.name)
            checks.append({'file': 'tables/' + target.name, 'rows': len(rows),
                           'content_sha256': actual, 'compressed_file_sha256': digest(target)})
        verify_manifest(input_root, input_root / 'input_manifest.json')
        summary = {
            'status': 'ANALYSIS_COMPLETE', 'version': __version__,
            'elapsed_seconds': time.monotonic() - started, 'python_version': platform.python_version(),
            'reproduction_scope': 'Derived native counters to all catalogue, signature-overlay and primary comparison tables',
            'raw_fastq_rescanned': False, 'synthetic_benchmark_rerun': False,
            'input_files_verified': len(records), 'scientific_output_tables_reproduced': len(checks),
            'protocol_sha256': digest(protocol_file), 'expected_output_manifest_sha256': digest(expected_file),
            'denominators': denominators, 'catalogue': catalog['summary'], 'selected_signatures': overlay['summary'],
            'target_repeat_summary': repeat_summary(profile), 'selected_case_structure': case_structure,
            'hash_basis': 'SHA-256 of decompressed table content; compressed bytes depend on the zlib implementation',
            'interpretation': 'Evidence agreement within a selected catalogue and same-study holdout accessions; no biological ground-truth classification.'}
        json_write(out / 'checks/reproduction.json', {'status': 'PASS', 'tables': checks})
        json_write(out / 'summary.json', summary)
        json_write(out / 'status.json', {'status': 'ANALYSIS_COMPLETE'})
        print('ANALYSIS_STATUS=ANALYSIS_COMPLETE')
        print('OUTPUT_DIRECTORY=' + str(out))
        return summary
    except BaseException as exc:
        json_write(out / 'status.json', {'status': 'FAILED', 'error': str(exc)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-root', type=Path, default=ROOT / 'data/derived')
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs/catalogue_audit')
    parser.add_argument('--verify-only', action='store_true', help='Check supplied input hashes without creating output')
    args = parser.parse_args()
    input_root = args.input_root.resolve()
    if args.verify_only:
        records = verify_manifest(input_root, input_root / 'input_manifest.json')
        print('INPUT_STATUS=VERIFIED')
        print('INPUT_FILES=' + str(len(records)))
        return
    reproduce(input_root, args.output, ROOT / 'config/catalog_protocol.json', ROOT / 'config/expected_catalog_outputs.json')


if __name__ == '__main__':
    main()
