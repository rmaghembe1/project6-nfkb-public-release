#!/usr/bin/env python3
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from collections import Counter, defaultdict

GENES = ("TNFAIP3","TRAF3","IKBKB","IKBKG","TNIP1","RBCK1")
EXPECTED_OUTPUTS = (
    "phaseF_variant_binding_v1.tsv",
    "phaseF_haplotype_frequency_v1.tsv",
    "phaseF_haplotype_sample_qc_v1.tsv",
    "phaseF_anchor_pairwise_LD_v1.tsv",
    "phaseF_context_LD_v1.tsv",
    "phaseF_TNFAIP3_historical_pair_v1.tsv",
    "phaseF_execution_summary_v1.json",
)

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        r = csv.DictReader(handle, delimiter="\t")
        return r.fieldnames or [], list(r)

def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def parse_key(key):
    parts = key.split(":")
    if len(parts) != 4:
        raise ValueError(f"invalid allele key: {key}")
    chrom, pos, ref, alt = parts
    return chrom, int(pos), ref, alt

def verify_payload_manifest(root):
    manifest = root / "bundle_payload_manifest.tsv"
    fields, rows = read_tsv(manifest)
    if fields != ["relative_path","bytes","sha256"]:
        raise SystemExit("bundle payload manifest schema mismatch")
    if len(rows) != 31:
        raise SystemExit(f"bundle payload manifest member count {len(rows)} !=31")
    listed = [row["relative_path"] for row in rows]
    if len(listed) != len(set(listed)):
        raise SystemExit("bundle payload manifest contains duplicate relative paths")
    listed_set = set(listed)
    actual_set = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path != manifest
    }
    if listed_set != actual_set:
        missing = sorted(listed_set - actual_set)
        extra = sorted(actual_set - listed_set)
        raise SystemExit(
            "bundle payload manifest file-set mismatch; "
            f"missing={missing[:5]} extra={extra[:5]}"
        )
    for row in rows:
        rel = Path(row["relative_path"])
        if rel.is_absolute() or ".." in rel.parts:
            raise SystemExit(f"unsafe bundle member path: {rel}")
        path = root / rel
        if path.stat().st_size != int(row["bytes"]):
            raise SystemExit(f"bundle byte mismatch: {rel}")
        if sha256(path) != row["sha256"]:
            raise SystemExit(f"bundle hash mismatch: {rel}")
    return rows

def required_paths(root):
    return {
        "selection": root / "inputs/interfaces/Project6_PhaseF_locus_selection_v1.tsv",
        "anchors": root / "inputs/interfaces/Project6_PhaseF_anchor_variant_representation_v1.tsv",
        "contract": root / "inputs/interfaces/Project6_PhaseF_haplotype_LD_execution_contract_v1.tsv",
        "strict": root / "inputs/interfaces/phaseB_strict_phased_interoperability_technical_inventory_v1.tsv",
        "review": root / "inputs/interfaces/phaseB_technical_review_sensitivity_inventory_v1.tsv",
        "order": root / "inputs/interfaces/TNFAIP3_genotype_interoperability_sample_order_v1.txt",
        "metadata": root / "inputs/metadata/gnomad_meta_updated.tsv",
        "eligibility": root / "inputs/interfaces/phaseE_population_eligibility_used_v1.tsv",
        "historical": root / "inputs/interfaces/Project6_TNFAIP3_historical_pair_exact_alleles_v1.tsv",
        "expected_outputs": root / "contracts/Project6_PhaseF_expected_outputs_v1.tsv",
    }

def load_metadata(path, order_path):
    fields, rows = read_tsv(path)
    required = {"s","population","population_inference.pop","hgdp_tgp_meta.Genetic.region","sex_imputation.sex_karyotype"}
    if not required.issubset(fields):
        raise SystemExit("metadata schema missing required fields")
    by = {r["s"]: r for r in rows}
    order = [x.strip() for x in Path(order_path).read_text(encoding="utf-8").splitlines() if x.strip()]
    if len(order) != 4091 or len(order) != len(set(order)):
        raise SystemExit("sample order must contain 4091 unique samples")
    if any(s not in by for s in order):
        raise SystemExit("sample order contains IDs absent from metadata")
    return order, by

def load_exact_key_sets(path):
    fields, rows = read_tsv(path)
    if "locus_scope" not in fields or "allele_key" not in fields:
        raise SystemExit(f"technical inventory schema mismatch: {path}")
    out = defaultdict(set)
    for row in rows:
        out[row["locus_scope"]].add(row["allele_key"])
    return out

def load_anchor_sets(anchor_path, historical_path):
    fields, rows = read_tsv(anchor_path)
    required = {"gene_symbol","allele_key","strict_phased_mask_member","technical_review_mask_member"}
    if not required.issubset(fields):
        raise SystemExit("anchor schema mismatch")
    strict = defaultdict(list)
    review = defaultdict(list)
    for row in rows:
        gene = row["gene_symbol"]
        key = row["allele_key"]
        if gene not in GENES:
            continue
        if row["strict_phased_mask_member"] == "YES":
            strict[gene].append(key)
        if row["technical_review_mask_member"] == "YES":
            review[gene].append(key)
    hfields, hrows = read_tsv(historical_path)
    if hfields != ["gene_symbol","historical_label","allele_key","binding_method"]:
        raise SystemExit("historical pair schema mismatch")
    historical = [r["allele_key"] for r in hrows if r["gene_symbol"] == "TNFAIP3"]
    if len(historical) != 2:
        raise SystemExit("historical TNFAIP3 pair must contain exactly two alleles")
    # deterministic dedup
    for d in (strict, review):
        for gene in list(d):
            d[gene] = sorted(set(d[gene]), key=lambda x: parse_key(x)[1:])
    historical = sorted(set(historical), key=lambda x: parse_key(x)[1:])
    return strict, review, historical

def bcf_for_gene(root, gene):
    p = root / "inputs/phased" / f"{gene}.bcf"
    i = root / "inputs/phased" / f"{gene}.bcf.csi"
    if not p.is_file() or not i.is_file():
        raise SystemExit(f"bundle phased resource missing for {gene}")
    return p

def preflight(root):
    verify_payload_manifest(root)
    p = required_paths(root)
    for name, path in p.items():
        if not path.is_file():
            raise SystemExit(f"required bundle input missing: {name}")
    order, meta = load_metadata(p["metadata"], p["order"])
    k = Counter(meta[s]["sex_imputation.sex_karyotype"] for s in order)
    if k != Counter({"XX":1896,"XY":2195}):
        raise SystemExit("phased cohort karyotype counts mismatch")
    strict = load_exact_key_sets(p["strict"])
    review = load_exact_key_sets(p["review"])
    if sum(len(v) for v in strict.values()) != 20793:
        raise SystemExit("strict phased key count mismatch")
    if sum(len(v) for v in review.values()) != 2393:
        raise SystemExit("review phased key count mismatch")
    aset, rset, hist = load_anchor_sets(p["anchors"], p["historical"])
    if sum(len(v) for v in aset.values()) != 18:
        raise SystemExit("strict selected anchor count mismatch")
    if sum(len(v) for v in rset.values()) != 2:
        raise SystemExit("review selected anchor count mismatch")
    for gene in GENES:
        bcf = bcf_for_gene(root, gene)
        proc = subprocess.run(["bcftools","query","-l",str(bcf)], check=True, text=True, capture_output=True)
        samples = [x for x in proc.stdout.splitlines() if x]
        if samples != order:
            raise SystemExit(f"sample order mismatch for {gene}")
    # exact historical representation, metadata only
    bcf = bcf_for_gene(root, "TNFAIP3")
    proc = subprocess.run(
        ["bcftools","query","-r","chr6:137874871-137875820","-f","%CHROM\t%POS\t%REF\t%ALT\n",str(bcf)],
        check=True,text=True,capture_output=True
    )
    observed=set()
    for line in proc.stdout.splitlines():
        if not line:
            continue
        chrom,pos,ref,alts=line.split("\t")
        for alt in alts.split(","):
            observed.add(f"{chrom}:{pos}:{ref}:{alt}")
    if not set(hist).issubset(observed):
        raise SystemExit("historical TNFAIP3 pair not represented exactly")
    return {
        "samples": len(order),
        "strict_keys": sum(len(v) for v in strict.values()),
        "review_keys": sum(len(v) for v in review.values()),
        "strict_selected_anchors": sum(len(v) for v in aset.values()),
        "review_selected_anchors": sum(len(v) for v in rset.values()),
    }

def split_gt(gt):
    if gt in ("",".","./.",".|."):
        return None, None
    sep = "|" if "|" in gt else ("/" if "/" in gt else None)
    if sep is None:
        if gt == ".":
            return None, None
        try:
            return (int(gt),), True
        except ValueError:
            return None, None
    parts = gt.split(sep)
    if any(x == "." for x in parts):
        return None, None
    try:
        vals = tuple(int(x) for x in parts)
    except ValueError:
        return None, None
    phased = (sep == "|")
    return vals, phased

def target_bits(gt, target_alt_index):
    vals, phased = split_gt(gt)
    if vals is None:
        return None
    bits = tuple(1 if a == target_alt_index else 0 for a in vals)
    if len(bits) == 1:
        return bits
    if phased:
        return bits
    # unphased homozygous is orientation-safe
    if len(set(bits)) == 1:
        return bits
    return None

def canonical_chromosome_bits(bits, karyotype, gene):
    if bits is None:
        return None
    if gene != "IKBKG":
        if len(bits) != 2:
            return None
        return bits
    if karyotype == "XX":
        if len(bits) != 2:
            return None
        return bits
    if karyotype == "XY":
        if len(bits) == 1:
            return (bits[0],)
        if len(bits) == 2 and bits[0] == bits[1]:
            return (bits[0],)
        return None
    return None

def ld_stats(a, b):
    if len(a) != len(b):
        raise ValueError("LD vectors differ in length")
    n = len(a)
    if n == 0:
        return {"n":0,"pA":None,"pB":None,"pAB":None,"D":None,"r2":None,"Dprime_signed":None,"Dprime_abs":None}
    pa = sum(a)/n
    pb = sum(b)/n
    pab = sum(x*y for x,y in zip(a,b))/n
    D = pab - pa*pb
    denom = pa*(1-pa)*pb*(1-pb)
    r2 = None if denom <= 0 else (D*D)/denom
    if D >= 0:
        dmax = min(pa*(1-pb),(1-pa)*pb)
    else:
        dmax = min(pa*pb,(1-pa)*(1-pb))
    dp = None if dmax <= 0 else D/dmax
    return {"n":n,"pA":pa,"pB":pb,"pAB":pab,"D":D,"r2":r2,"Dprime_signed":dp,"Dprime_abs":None if dp is None else abs(dp)}

def selftest():
    s = ld_stats([0,0,1,1],[0,0,1,1])
    if not math.isclose(s["r2"],1.0) or not math.isclose(s["Dprime_abs"],1.0):
        raise SystemExit("LD perfect-correlation selftest failed")
    s = ld_stats([0,0,1,1],[0,1,0,1])
    if not math.isclose(s["r2"],0.0, abs_tol=1e-12):
        raise SystemExit("LD independence selftest failed")
    if target_bits("0|1",1) != (0,1):
        raise SystemExit("phased GT parsing selftest failed")
    if target_bits("1/0",1) is not None:
        raise SystemExit("unphased heterozygote handling selftest failed")
    if canonical_chromosome_bits((1,1),"XY","IKBKG") != (1,):
        raise SystemExit("IKBKG XY pseudo-diploid selftest failed")
    print("PHASEF_EXECUTOR_SYNTHETIC_SELFTEST=PASS")

def stream_target_records(root, gene, target_keys, order):
    """Yield exact allele key -> (target ALT index, GT vector) without retaining the full locus."""
    targets = defaultdict(list)
    target_set = set(target_keys)
    for key in sorted(target_set, key=lambda x: parse_key(x)[1:]):
        chrom,pos,ref,alt = parse_key(key)
        targets[(chrom,pos,ref)].append((alt,key))
    bcf = bcf_for_gene(root,gene)
    proc = subprocess.Popen(
        ["bcftools","query","-S",str(root/"inputs/interfaces/TNFAIP3_genotype_interoperability_sample_order_v1.txt"),
         "-f","%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",str(bcf)],
        text=True, stdout=subprocess.PIPE
    )
    seen=set()
    assert proc.stdout is not None
    for line in proc.stdout:
        parts=line.rstrip("\n").split("\t")
        if len(parts) < 5:
            continue
        chrom,pos_s,ref,alts = parts[:4]
        pos=int(pos_s)
        key3=(chrom,pos,ref)
        if key3 not in targets:
            continue
        alt_list=alts.split(",")
        gts=parts[4:]
        if len(gts) != len(order):
            raise SystemExit(f"{gene}: GT sample width mismatch at {chrom}:{pos}")
        for alt,key in targets[key3]:
            if alt not in alt_list:
                continue
            seen.add(key)
            yield key,(alt_list.index(alt)+1,gts)
    rc=proc.wait()
    if rc != 0:
        raise SystemExit(f"bcftools query failed for {gene}: rc={rc}")
    missing=target_set-seen
    if missing:
        raise SystemExit(f"{gene}: target exact alleles absent: {sorted(missing)[:5]}")

def query_target_calls(root, gene, target_keys, order):
    return dict(stream_target_records(root,gene,target_keys,order))

def sample_chrom_bits(gene, key_record, order, metadata):
    alt_index,gts=key_record
    out=[]
    for sample,gt in zip(order,gts):
        bits=target_bits(gt,alt_index)
        out.append(canonical_chromosome_bits(bits,metadata[sample]["sex_imputation.sex_karyotype"],gene))
    return out

def chromosomes_for_group(gene, key_calls, sample_indices):
    vals=[]
    for idx in sample_indices:
        bits=key_calls[idx]
        if bits is None:
            continue
        vals.extend(bits)
    return vals

def paired_chromosomes(gene, calls_a, calls_b, sample_indices):
    a=[]
    b=[]
    for idx in sample_indices:
        xa=calls_a[idx]
        xb=calls_b[idx]
        if xa is None or xb is None or len(xa) != len(xb):
            continue
        a.extend(xa)
        b.extend(xb)
    return a,b

def format_num(x):
    if x is None:
        return "NA"
    return format(x,".12g")

def execute(root, output_dir, auth_file, auth_sha):
    if output_dir.exists():
        raise SystemExit("output directory already exists")
    if not auth_file.is_file() or sha256(auth_file) != auth_sha:
        raise SystemExit("authorization file missing or hash mismatch")
    auth_text=auth_file.read_text(encoding="utf-8")
    required_tokens=(
        "PROJECT6_PHASEF_EXECUTION_AUTHORIZED=YES",
        "BUNDLE_PAYLOAD_MANIFEST_SHA256="+sha256(root/"bundle_payload_manifest.tsv"),
    )
    if any(tok not in auth_text for tok in required_tokens):
        raise SystemExit("authorization file lacks required binding token")
    pf=preflight(root)
    p=required_paths(root)
    order,metadata=load_metadata(p["metadata"],p["order"])
    strict_inventory=load_exact_key_sets(p["strict"])
    strict_anchors, review_anchors, historical_pair=load_anchor_sets(p["anchors"],p["historical"])

    # group index definitions
    all_indices=list(range(len(order)))
    primary=defaultdict(list); inferred=defaultdict(list); broad=defaultdict(list)
    for idx,s in enumerate(order):
        m=metadata[s]
        primary[m["population"]].append(idx)
        inf=m["population_inference.pop"]
        if inf not in ("","NA","N/A"):
            inferred[inf].append(idx)
        broad[m["hgdp_tgp_meta.Genetic.region"]].append(idx)

    set_defs=[]
    for gene in GENES:
        if strict_anchors.get(gene):
            set_defs.append(("PRIMARY_STRICT_ANCHORS",gene,strict_anchors[gene]))
        if review_anchors.get(gene):
            set_defs.append(("TECHNICAL_REVIEW_ANCHORS",gene,review_anchors[gene]))
    set_defs.append(("TNFAIP3_HISTORICAL_PAIR","TNFAIP3",historical_pair))

    variant_binding=[]
    hap_rows=[]
    hap_qc=[]
    anchor_ld=[]
    context_ld=[]
    tnfaip3_rows=[]
    query_cache={}

    for set_id,gene,keys in set_defs:
        keys=sorted(set(keys), key=lambda x: parse_key(x)[1:])
        if gene not in query_cache:
            all_gene_keys=set(strict_anchors.get(gene,[]))|set(review_anchors.get(gene,[]))
            if gene=="TNFAIP3":
                all_gene_keys |= set(historical_pair)
            raw=query_target_calls(root,gene,all_gene_keys,order)
            query_cache[gene]={k:sample_chrom_bits(gene,rec,order,metadata) for k,rec in raw.items()}
        calls={k:query_cache[gene][k] for k in keys}
        for key in keys:
            chrom,pos,ref,alt=parse_key(key)
            variant_binding.append({
                "analysis_set":set_id,"gene_symbol":gene,"allele_key":key,
                "chrom":chrom,"pos":str(pos),"ref":ref,"alt":alt,
                "strict_inventory_member":"YES" if key in strict_inventory.get(gene,set()) else "NO",
                "binding_status":"PASS"
            })

        group_specs=[("GLOBAL","ALL","GLOBAL",all_indices)]
        group_specs += [("PRIMARY_N20",g,"PRIMARY",idxs) for g,idxs in sorted(primary.items()) if len(idxs)>=20]
        group_specs += [("PRIMARY_N30",g,"SENSITIVITY",idxs) for g,idxs in sorted(primary.items()) if len(idxs)>=30]
        group_specs += [("INFERRED",g,"SENSITIVITY",idxs) for g,idxs in sorted(inferred.items())]
        group_specs += [("BROAD",g,"SENSITIVITY",idxs) for g,idxs in sorted(broad.items())]

        # Haplotype counts require >=2 variants.
        if len(keys) >= 2:
            for group_type,group,elig,idxs in group_specs:
                haps=Counter()
                complete_samples=0
                complete_chromosomes=0
                excluded_samples=0
                for idx in idxs:
                    per=[calls[k][idx] for k in keys]
                    if any(x is None for x in per):
                        excluded_samples += 1
                        continue
                    widths={len(x) for x in per}
                    if len(widths)!=1:
                        excluded_samples += 1
                        continue
                    width=next(iter(widths))
                    complete_samples += 1
                    complete_chromosomes += width
                    for h in range(width):
                        hap="".join(str(per[j][h]) for j in range(len(keys)))
                        haps[hap]+=1
                hap_qc.append({
                    "analysis_set":set_id,"gene_symbol":gene,"group_type":group_type,"group":group,
                    "sample_count":str(len(idxs)),"complete_samples":str(complete_samples),
                    "excluded_samples":str(excluded_samples),"complete_chromosomes":str(complete_chromosomes),
                    "variant_count":str(len(keys)),"variant_keys":";".join(keys)
                })
                for hap,count in sorted(haps.items(), key=lambda kv:(-kv[1],kv[0])):
                    hap_rows.append({
                        "analysis_set":set_id,"gene_symbol":gene,"group_type":group_type,"group":group,
                        "eligibility":elig,"sample_count":str(len(idxs)),"variant_count":str(len(keys)),
                        "variant_keys":";".join(keys),"haplotype":hap,"haplotype_count":str(count),
                        "complete_chromosomes":str(complete_chromosomes),
                        "haplotype_frequency":format_num(count/complete_chromosomes if complete_chromosomes else None)
                    })

            # Anchor pairwise LD: global, primary N>=20, and frozen N>=30 sensitivity.
            ld_group_specs=[("GLOBAL","ALL",all_indices)]
            ld_group_specs += [("PRIMARY_N20",g,idxs) for g,idxs in sorted(primary.items()) if len(idxs)>=20]
            ld_group_specs += [("PRIMARY_N30",g,idxs) for g,idxs in sorted(primary.items()) if len(idxs)>=30]
            for ka,kb in itertools.combinations(keys,2):
                for group_type,group,idxs in ld_group_specs:
                    a,b=paired_chromosomes(gene,calls[ka],calls[kb],idxs)
                    s=ld_stats(a,b)
                    anchor_ld.append({
                        "analysis_set":set_id,"gene_symbol":gene,"group_type":group_type,"group":group,
                        "allele_key_A":ka,"allele_key_B":kb,"n_chromosomes":str(s["n"]),
                        "pA":format_num(s["pA"]),"pB":format_num(s["pB"]),"pAB":format_num(s["pAB"]),
                        "D":format_num(s["D"]),"r2":format_num(s["r2"]),
                        "D_prime_signed":format_num(s["Dprime_signed"]),
                        "D_prime_abs":format_num(s["Dprime_abs"]),
                        "evaluable":"YES" if s["r2"] is not None else "NO"
                    })

    # focal TNFAIP3 rows from haplotype and LD tables
    for row in hap_rows:
        if row["analysis_set"]=="TNFAIP3_HISTORICAL_PAIR":
            tnfaip3_rows.append({
                "record_type":"HAPLOTYPE","group_type":row["group_type"],"group":row["group"],
                "metric":row["haplotype"],"value":row["haplotype_frequency"],
                "numerator":row["haplotype_count"],"denominator":row["complete_chromosomes"],
                "detail":row["variant_keys"]
            })
    for row in anchor_ld:
        if row["analysis_set"]=="TNFAIP3_HISTORICAL_PAIR":
            tnfaip3_rows.append({
                "record_type":"LD","group_type":row["group_type"],"group":row["group"],
                "metric":"r2","value":row["r2"],"numerator":"NA","denominator":row["n_chromosomes"],
                "detail":row["allele_key_A"]+";"+row["allele_key_B"]+";Dprime_abs="+row["D_prime_abs"]
            })

    # Context LD: global only, each PRIMARY strict anchor against strict inventory alleles within 50kb.
    # Historical pair is also treated as an anchor for TNFAIP3.
    for gene in GENES:
        anchors=set(strict_anchors.get(gene,[]))
        if gene=="TNFAIP3":
            anchors |= set(historical_pair)
        if not anchors:
            continue
        context=set()
        for ck in strict_inventory.get(gene,set()):
            cchrom,cpos,cref,calt=parse_key(ck)
            for ak in anchors:
                achrom,apos,aref,aalt=parse_key(ak)
                if cchrom==achrom and abs(cpos-apos)<=50000:
                    context.add(ck); break
        # Hold only the small anchor matrix in memory; stream context variants one at a time.
        anchor_raw=query_target_calls(root,gene,anchors,order)
        anchor_bits={k:sample_chrom_bits(gene,rec,order,metadata) for k,rec in anchor_raw.items()}
        del anchor_raw
        for ck,rec in stream_target_records(root,gene,context,order):
            cchrom,cpos,_,_=parse_key(ck)
            context_bits=sample_chrom_bits(gene,rec,order,metadata)
            for ak in sorted(anchors, key=lambda x:parse_key(x)[1:]):
                if ck==ak:
                    continue
                achrom,apos,_,_=parse_key(ak)
                if achrom!=cchrom or abs(cpos-apos)>50000:
                    continue
                a,b=paired_chromosomes(gene,anchor_bits[ak],context_bits,all_indices)
                s=ld_stats(a,b)
                context_ld.append({
                    "gene_symbol":gene,"anchor_allele_key":ak,"context_allele_key":ck,
                    "distance_bp":str(abs(cpos-apos)),"n_chromosomes":str(s["n"]),
                    "p_anchor":format_num(s["pA"]),"p_context":format_num(s["pB"]),
                    "D":format_num(s["D"]),"r2":format_num(s["r2"]),
                    "D_prime_signed":format_num(s["Dprime_signed"]),
                    "D_prime_abs":format_num(s["Dprime_abs"]),
                    "evaluable":"YES" if s["r2"] is not None else "NO"
                })

    output_dir.mkdir(parents=True)
    write_tsv(output_dir/"phaseF_variant_binding_v1.tsv",
              ["analysis_set","gene_symbol","allele_key","chrom","pos","ref","alt","strict_inventory_member","binding_status"],variant_binding)
    write_tsv(output_dir/"phaseF_haplotype_frequency_v1.tsv",
              ["analysis_set","gene_symbol","group_type","group","eligibility","sample_count","variant_count","variant_keys","haplotype","haplotype_count","complete_chromosomes","haplotype_frequency"],hap_rows)
    write_tsv(output_dir/"phaseF_haplotype_sample_qc_v1.tsv",
              ["analysis_set","gene_symbol","group_type","group","sample_count","complete_samples","excluded_samples","complete_chromosomes","variant_count","variant_keys"],hap_qc)
    write_tsv(output_dir/"phaseF_anchor_pairwise_LD_v1.tsv",
              ["analysis_set","gene_symbol","group_type","group","allele_key_A","allele_key_B","n_chromosomes","pA","pB","pAB","D","r2","D_prime_signed","D_prime_abs","evaluable"],anchor_ld)
    write_tsv(output_dir/"phaseF_context_LD_v1.tsv",
              ["gene_symbol","anchor_allele_key","context_allele_key","distance_bp","n_chromosomes","p_anchor","p_context","D","r2","D_prime_signed","D_prime_abs","evaluable"],context_ld)
    write_tsv(output_dir/"phaseF_TNFAIP3_historical_pair_v1.tsv",
              ["record_type","group_type","group","metric","value","numerator","denominator","detail"],tnfaip3_rows)

    summary={
        "status":"PASS","phase":"Phase_F","design":"bounded_anchor_haplotype_and_local_LD",
        "bundle_payload_manifest_sha256":sha256(root/"bundle_payload_manifest.tsv"),
        "authorization_sha256":auth_sha,
        "samples":len(order),
        "population_group_counts":{
            "primary_total":len(primary),
            "primary_n20":sum(len(idxs)>=20 for idxs in primary.values()),
            "primary_n30":sum(len(idxs)>=30 for idxs in primary.values()),
            "inferred":len(inferred),
            "broad":len(broad),
        },
        "candidate_loci":list(GENES),
        "strict_selected_anchor_count":sum(len(v) for v in strict_anchors.values()),
        "technical_review_selected_anchor_count":sum(len(v) for v in review_anchors.values()),
        "historical_TNFAIP3_pair":historical_pair,
        "counts":{
            "variant_binding_rows":len(variant_binding),
            "haplotype_frequency_rows":len(hap_rows),
            "haplotype_sample_qc_rows":len(hap_qc),
            "anchor_pairwise_LD_rows":len(anchor_ld),
            "context_LD_rows":len(context_ld),
            "TNFAIP3_historical_pair_rows":len(tnfaip3_rows),
        },
        "claim_boundaries":{
            "positive_selection":False,"adaptation":False,"causal_function":False,
            "disease_association":False,"ancestry_specific_function":False
        },
        "preflight":pf,
    }
    with (output_dir/"phaseF_execution_summary_v1.json").open("w",encoding="utf-8") as h:
        json.dump(summary,h,indent=2,sort_keys=True); h.write("\n")

    observed=sorted(p.name for p in output_dir.iterdir() if p.is_file())
    if observed != sorted(EXPECTED_OUTPUTS):
        raise SystemExit("unexpected Phase-F output set")
    print("PHASEF_EXECUTION_STATUS=PASS")
    print("PHASEF_OUTPUT_FILE_COUNT=7")
    print("PHASEF_OUTPUT_SET=PASS")

def main():
    ap=argparse.ArgumentParser()
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--selftest",action="store_true")
    mode.add_argument("--preflight",action="store_true")
    mode.add_argument("--execute",action="store_true")
    ap.add_argument("--bundle-root")
    ap.add_argument("--output-dir")
    ap.add_argument("--authorization-file")
    ap.add_argument("--authorization-sha")
    args=ap.parse_args()
    if args.selftest:
        selftest(); return
    if not args.bundle_root:
        raise SystemExit("--bundle-root is required")
    root=Path(args.bundle_root).resolve()
    if args.preflight:
        info=preflight(root)
        print("PHASEF_BUNDLE_PREFLIGHT=PASS")
        for k,v in sorted(info.items()):
            print(f"{k.upper()}={v}")
        print("GENOTYPE_VALUES_QUERIED=NO")
        print("HAPLOTYPES_CALCULATED=NO")
        print("LD_CALCULATED=NO")
        return
    if args.execute:
        if not args.output_dir or not args.authorization_file or not args.authorization_sha:
            raise SystemExit("--execute requires --output-dir, --authorization-file, --authorization-sha")
        execute(root,Path(args.output_dir),Path(args.authorization_file),args.authorization_sha)
        return

if __name__=="__main__":
    main()
