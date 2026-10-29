#!/usr/bin/env python3
"""Generate a checksum ledger for locally acquired resources without committing large source data."""
from pathlib import Path
import argparse, hashlib, csv

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--directory', default='data/external'); ap.add_argument('--out', default='data/metadata/external_file_checksums.tsv'); args=ap.parse_args()
    root=Path(args.directory); out=Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    files=[p for p in root.rglob('*') if p.is_file() and p.name != '.keepdir']
    with out.open('w', encoding='utf-8', newline='') as f:
        w=csv.writer(f, delimiter='\t'); w.writerow(['relative_path','size_bytes','sha256'])
        for p in sorted(files): w.writerow([str(p.relative_to(root)),p.stat().st_size,digest(p)])
    print(f'Checksum ledger written for {len(files)} files: {out}')
if __name__=='__main__': main()
