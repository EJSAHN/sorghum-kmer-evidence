"""Deterministic tables and explicit input/output integrity checks."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 ** 2), b''):
            h.update(block)
    return h.hexdigest()


def content_digest(path):
    """SHA-256 of the file content, decompressing gzip files first.

    Different zlib implementations (for example zlib and zlib-ng) can compress
    identical text to different bytes. Output checks therefore compare the
    decompressed table content rather than the compressed file.
    """
    path = Path(path)
    h = hashlib.sha256()
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 ** 2), b''):
            h.update(block)
    return h.hexdigest()


def json_read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')


def write_table(path, rows):
    rows = list(rows)
    fields = list(dict.fromkeys(key for row in rows for key in row)) or ['status']
    with Path(path).open('wb') as raw, gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
        with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as text:
            writer = csv.DictWriter(text, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)


def verify_manifest(root, manifest):
    root = Path(root).resolve()
    manifest_data = json_read(manifest)
    records = manifest_data.get('files') if isinstance(manifest_data, dict) else manifest_data
    if not isinstance(records, list) or not records:
        raise ValueError('Input manifest must be a nonempty list')
    seen = set()
    for row in records:
        name = row['file']
        relative = Path(name)
        if (not isinstance(name, str) or relative.is_absolute() or '..' in relative.parts
                or '\\' in name or name in seen or relative.as_posix() != name):
            raise ValueError('Unsafe or repeated manifest filename: ' + str(name))
        seen.add(name)
        path = root / relative
        if not path.resolve().is_relative_to(root) or not path.is_file():
            raise ValueError('Missing input or input outside root: ' + name)
        if path.stat().st_size != row['size_bytes'] or digest(path) != row['sha256']:
            raise ValueError('Input checksum mismatch: ' + name)
    return records


def new_output(path, inputs=()):
    out = Path(path).resolve()
    for root in inputs:
        root = Path(root).resolve()
        if out == root or out.is_relative_to(root):
            raise ValueError('Output directory must be outside the supplied inputs')
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError('Output directory exists and is not empty; choose a new directory')
    out.mkdir(parents=True, exist_ok=True)
    return out
