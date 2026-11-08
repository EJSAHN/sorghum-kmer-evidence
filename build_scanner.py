#!/usr/bin/env python3
"""Compile the optional streaming FASTQ scanner with a local C++17 compiler."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))
from kmer_evidence.io import digest, json_write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='g++')
    parser.add_argument('--output', type=Path, default=ROOT / 'bin/panel_scan')
    args = parser.parse_args()
    compiler = shutil.which(args.compiler)
    if compiler is None:
        parser.error('Compiler not found. Use an existing C++17 compiler and zlib development headers.')
    output = args.output.resolve()
    if output.exists():
        parser.error('Build output already exists; choose a new output path')
    output.parent.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'native/panel_scan.cpp'
    command = [compiler, '-std=c++17', '-O3', '-Wall', '-Wextra', str(source), '-lz', '-o', str(output)]
    subprocess.run(command, check=True)
    version = subprocess.run([compiler, '--version'], check=True, text=True, capture_output=True).stdout.splitlines()[0]
    json_write(output.with_name(output.name + '.build.json'), {
        'compiler': version, 'flags': ['-std=c++17', '-O3', '-Wall', '-Wextra', '-lz'],
        'source_sha256': digest(source), 'binary_sha256': digest(output), 'source': 'native/panel_scan.cpp'})
    print('SCANNER=' + str(output))


if __name__ == '__main__':
    main()
