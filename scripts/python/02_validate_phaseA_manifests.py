#!/usr/bin/env python3
"""Validate that Phase A manifests contain required evidence categories and pending actions."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path

def read_tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--resources', required=True); ap.add_argument('--anchors', required=True); ap.add_argument('--literature', required=True); ap.add_argument('--out', required=True)
    a=ap.parse_args()
    r=read_tsv(a.resources); an=read_tsv(a.anchors); lit=read_tsv(a.literature)
    resource_text=' '.join(' '.join(str(v) for v in x.values()) for x in r).lower()
    r10_rows=[x for x in r if x.get('resource_id') == 'R10']
    r10=r10_rows[0] if len(r10_rows) == 1 else {}
    r10_text=' '.join(str(v) for v in r10.values()).lower()
    anchor_text=' '.join(' '.join(str(v) for v in x.values()) for x in an).lower()
    a1_rows=[x for x in an if x.get('anchor_id') == 'A1']
    a1=a1_rows[0] if len(a1_rows) == 1 else {}
    a1_text=' '.join(str(v) for v in a1.values()).lower()
    lit_text=' '.join(x['title'] for x in lit).lower()
    lit_all_text=' '.join(' '.join(str(v) for v in x.values()) for x in lit).lower()
    closure_path=Path('results/phaseA/formal_competitor_search_closure_decision_v1.json')
    transition_path=Path('15_logs_and_audit/competitor_search/Gate1_to_PhaseB_Technical_Extraction_Transition_Decision_2026-05-28.md')
    closure=json.loads(closure_path.read_text(encoding='utf-8')) if closure_path.exists() else {}
    transition_text=transition_path.read_text(encoding='utf-8').lower() if transition_path.exists() else ''
    checks={
      'primary_wgs_resource_recorded': any(
          x['resource_id'] == 'R01'
          and 'gnomad v3.1.2' in x['resource'].lower()
          and 'hgdp+1kgp' in x['resource'].lower()
          and 'indexed per-chromosome genotype vcfs' in x['resource'].lower()
          for x in r
      ),
      'mane_recorded': 'mane' in resource_text,
      'structural_variant_scope_resolved': (
          len(r10_rows) == 1
          and r10.get('layer') == 'Structural variants'
          and 'hgdp+1kgp-compatible public sv layer' in r10_text
          and 'negative feasibility decision' in r10_text
          and 'excluded from executable primary workflow' in r10_text
          and 'conditional sv/cnv extension only' in r10_text
      ),
      'encode_ccre_recorded': 'encode4' in resource_text or 'ccre' in resource_text,
      'response_qtl_context_recorded': 'kim et al.' in resource_text and 'fairfax' in resource_text,
      'tnfaip3_anchor_recorded': 'tnfaip3' in anchor_text and 'denisovan' in anchor_text,
      'tnfaip3_grch38_anchor_resolved': (
          len(a1_rows) == 1
          and a1.get('gene') == 'TNFAIP3'
          and a1.get('status', '').lower() == 'resolved'
          and 'rs376205580' in a1_text
          and '137874871' in a1_text
          and 'a>g' in a1_text
          and 'rs141807543' in a1_text
          and '137875820' in a1_text
          and 'a>c' in a1_text
          and 'thr108ala' in a1_text
          and 'ile207leu' in a1_text
      ),
      'primary_resource_publication_recorded': 'harmonized public resource' in lit_text,
      'recent_nfkb1_novelty_boundary_recorded': 'alternatively spliced nfkb1' in lit_text,
      'expanded_novelty_boundary_literature_recorded': (
          'mapping the gene regulatory landscape of archaic hominin introgression' in lit_text
          and 'impact of constitutional copy number variants on biological pathway evolution' in lit_text
          and 'genomic and geographic distribution of private snps and pathways' in lit_text
          and 'african-derived genetic polymorphisms in tnfaip3' in lit_text
          and 'genetic interaction between genes involved in nf-κb signaling pathway' in lit_text
      ),
      'formal_competitor_search_completed': (
          closure.get('status') == 'FORMAL_COMPETITOR_SEARCH_COMPLETE_NO_DIRECT_COMPETITOR_IDENTIFIED'
          and closure.get('gate1_competitor_search_component') == 'PASS'
          and closure.get('primary_formal_search', {}).get('direct_competitors_identified') == 0
          and closure.get('central_scope_verification', {}).get('direct_competitors_identified') == 0
          and closure.get('supplementary_direct_design_search', {}).get('direct_competitors_identified') == 0
          and closure.get('citation_chaining', {}).get('direct_competitors_identified') == 0
      ),
      'phaseB_extraction_lock_transition_recorded': (
          'first controlled phase b technical action' in transition_text
          and 'biological interpretation' in transition_text
          and 'remain prohibited until extracted-resource qc and checksum-lock validation have passed' in transition_text
      ),
    }
    all_pass=all(checks.values())
    payload={
        'status':'PASS' if all_pass else 'FAIL',
        'gate1_status':'CLOSED_PHASE_B_TECHNICAL_EXTRACTION_AUTHORISED' if all_pass else 'OPEN',
        'checks':checks,
        'resource_records':len(r),
        'anchor_records':len(an),
        'literature_records':len(lit),
        'remaining_gate1_requirements':[] if all_pass else [
            'Resolve failed Phase A checks before closing Gate 1.'
        ],
        'first_authorised_phaseB_action': (
            'Interval-restricted genotype and phased-haplotype extraction with exact URL recording, checksum locking and QC only.'
            if all_pass else 'Not authorised until Gate 1 checks pass.'
        ),
        'interpretation_boundary': (
            'Biological interpretation, population-frequency comparison, haplotype interpretation, selection analysis and functional inference remain prohibited until extracted-resource QC and checksum-lock validation pass.'
            if all_pass else 'No Phase B activity authorised.'
        ),
    }
    out=Path(a.out); out.parent.mkdir(parents=True, exist_ok=True); out.write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(payload, indent=2))
    if payload['status'] != 'PASS': raise SystemExit(1)
if __name__=='__main__': main()
