#!/usr/bin/env python3
"""Freeze and validate the pre-specified Project 6 NF-kB gene panel."""
from __future__ import annotations
import argparse, csv, json
from collections import Counter, defaultdict
from pathlib import Path

def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--input', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--qc', required=True)
    args=ap.parse_args()
    inp=Path(args.input); out=Path(args.output); qc=Path(args.qc)
    with inp.open(encoding='utf-8', newline='') as f:
        rows=list(csv.DictReader(f, delimiter='\t'))
    primary=[r for r in rows if r['analysis_tier']=='Primary']
    optional=[r for r in rows if r['analysis_tier']=='Optional']
    symbols=[r['gene_symbol'] for r in rows]
    duplicate=[g for g,n in Counter(symbols).items() if n>1]
    by_module=defaultdict(int)
    for r in primary: by_module[r['module_id']]+=1
    expected={'M1':5,'M2':7,'M3':5,'M4':7}
    checks={
        'exactly_24_primary_genes': len(primary)==24,
        'exactly_3_optional_genes': len(optional)==3,
        'no_duplicate_gene_symbols': len(duplicate)==0,
        'module_counts_match_protocol': dict(by_module)==expected,
        'tnfaip3_anchor_present': any(r['gene_symbol']=='TNFAIP3' for r in primary),
        'ikbkg_xlinked_qc_gene_present': any(r['gene_symbol']=='IKBKG' for r in primary),
    }
    if not all(checks.values()):
        raise SystemExit('Gene panel validation failed: '+json.dumps(checks, indent=2))
    out.parent.mkdir(parents=True, exist_ok=True); qc.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w', encoding='utf-8', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
        writer.writeheader(); writer.writerows(primary)
    payload={'status':'PASS','primary_gene_count':len(primary),'optional_gene_count':len(optional),'primary_module_counts':dict(by_module),'checks':checks,'duplicate_symbols':duplicate}
    qc.write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(payload, indent=2))
if __name__=='__main__': main()
