#!/usr/bin/env python3

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import re
import subprocess
import sys
from collections import Counter, defaultdict


BASE_SEED = 20261011
PERMUTATIONS = 100000

OUTPUT_NAMES = (
    "phaseE_population_eligibility_used_v1.tsv",
    "phaseE_variant_population_effects_v1.tsv",
    "phaseE_gene_architecture_v1.tsv",
    "phaseE_module_architecture_v1.tsv",
    "phaseE_global_module_permutation_v1.tsv",
    "phaseE_pairwise_module_contrasts_v1.tsv",
    "phaseE_sensitivity_summary_v1.tsv",
    "phaseE_prioritized_variant_ranking_v1.tsv",
    "phaseE_prioritized_genotype_state_v1.tsv",
    "phaseE_execution_summary_v1.json",
)


class Refusal(RuntimeError):
    pass


def sha256(path):
    path = Path(path)
    h = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


def read_tsv(path):
    with Path(path).open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )
        return reader.fieldnames or [], list(reader)


def write_tsv(path, fieldnames, rows):
    path = Path(path)

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_key_value_file(path):
    values = {}

    for line in Path(path).read_text(
        encoding="utf-8"
    ).splitlines():
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        if "=" not in line:
            continue

        key, value = line.split(
            "=",
            1,
        )

        key = key.strip()
        value = value.strip()

        if key in values:
            raise Refusal(
                "duplicate authorization key: "
                + key
            )

        values[key] = value

    return values


def read_binding(path):
    fields, rows = read_tsv(path)

    required = {
        "kind",
        "key",
        "value",
        "sha256",
        "role",
    }

    if not required.issubset(
        set(fields)
    ):
        raise Refusal(
            "implementation binding schema mismatch"
        )

    keys = {}

    for row in rows:
        key = row["key"].strip()

        if key in keys:
            raise Refusal(
                "duplicate implementation binding key "
                + key
            )

        keys[key] = row

    return rows, keys


def verify_binding_files(rows):
    count = 0

    for row in rows:
        expected = row[
            "sha256"
        ].strip()

        if not expected:
            continue

        path = Path(
            row["value"].strip()
        )

        if not path.is_file():
            raise Refusal(
                "bound file missing: "
                + row["key"]
            )

        actual = sha256(path)

        if actual != expected:
            raise Refusal(
                "bound file hash mismatch: "
                + row["key"]
            )

        count += 1

    return count


def verify_schema_binding(path):
    fields, rows = read_tsv(path)

    required = {
        "surface",
        "path",
        "analysis_layer",
        "grouping",
        "expected_rows",
        "expected_unique_keys",
        "distinct_genes",
        "gene_field",
        "allele_key_field",
        "group_field",
        "ac_field",
        "an_field",
        "af_field",
    }

    if not required.issubset(
        set(fields)
    ):
        raise Refusal(
            "D8 schema binding schema mismatch"
        )

    if len(rows) != 6:
        raise Refusal(
            "D8 schema binding rows !=6"
        )

    mapping = {}

    for row in rows:
        surface = row[
            "surface"
        ].strip()

        if surface in mapping:
            raise Refusal(
                "duplicate D8 surface "
                + surface
            )

        path = Path(
            row["path"]
        )

        if not path.is_file():
            raise Refusal(
                "D8 source missing "
                + surface
            )

        with path.open(
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.reader(
                handle,
                delimiter="\t",
            )
            header = next(reader)

        expected_fields = {
            row["gene_field"],
            row["allele_key_field"],
            row["group_field"],
            row["ac_field"],
            row["an_field"],
            row["af_field"],
        }

        if not expected_fields.issubset(
            set(header)
        ):
            raise Refusal(
                "D8 source schema drift "
                + surface
            )

        mapping[surface] = row

    required_surfaces = {
        "AUTO_POPULATION",
        "AUTO_INFERRED",
        "AUTO_BROAD",
        "CHRX_POPULATION",
        "CHRX_INFERRED",
        "CHRX_BROAD",
    }

    if set(mapping) != required_surfaces:
        raise Refusal(
            "D8 surface identity mismatch"
        )

    return mapping


def verify_source_binding(path):
    fields, rows = read_tsv(path)

    required = {
        "gene_symbol",
        "absolute_path",
        "bytes",
        "sha256",
    }

    if not required.issubset(
        set(fields)
    ):
        raise Refusal(
            "genotype source binding schema mismatch"
        )

    if len(rows) != 24:
        raise Refusal(
            "genotype source binding rows !=24"
        )

    sources = {}

    for row in rows:
        gene = row[
            "gene_symbol"
        ].strip()

        path = Path(
            row["absolute_path"]
        )

        expected_bytes = int(
            row["bytes"]
        )

        expected_sha = row[
            "sha256"
        ].strip()

        if gene in sources:
            raise Refusal(
                "duplicate genotype source gene "
                + gene
            )

        if not path.is_file():
            raise Refusal(
                "genotype source missing "
                + gene
            )

        if path.stat().st_size != expected_bytes:
            raise Refusal(
                "genotype source byte mismatch "
                + gene
            )

        if sha256(path) != expected_sha:
            raise Refusal(
                "genotype source hash mismatch "
                + gene
            )

        sources[gene] = path

    if len(sources) != 24:
        raise Refusal(
            "genotype source gene count !=24"
        )

    return sources


def median(values):
    values = sorted(
        float(value)
        for value in values
    )

    if not values:
        raise ValueError(
            "median of empty sequence"
        )

    n = len(values)
    middle = n // 2

    if n % 2:
        return values[middle]

    return (
        values[middle - 1]
        + values[middle]
    ) / 2.0


def quantile_linear(values, q):
    values = sorted(
        float(value)
        for value in values
    )

    if not values:
        raise ValueError(
            "quantile of empty sequence"
        )

    if not 0.0 <= q <= 1.0:
        raise ValueError(
            "q outside [0,1]"
        )

    if len(values) == 1:
        return values[0]

    h = (
        len(values) - 1
    ) * q

    lo = math.floor(h)
    hi = math.ceil(h)

    if lo == hi:
        return values[lo]

    return (
        values[lo]
        + (
            h - lo
        )
        * (
            values[hi]
            - values[lo]
        )
    )


def average_ranks(values):
    indexed = sorted(
        enumerate(values),
        key=lambda item: (
            item[1],
            item[0],
        ),
    )

    ranks = [
        0.0
        for _ in values
    ]

    ties = []

    i = 0

    while i < len(indexed):
        j = i + 1

        while (
            j < len(indexed)
            and indexed[j][1]
            == indexed[i][1]
        ):
            j += 1

        average = (
            (i + 1)
            + j
        ) / 2.0

        for k in range(i, j):
            ranks[
                indexed[k][0]
            ] = average

        ties.append(
            j - i
        )

        i = j

    return ranks, ties


def kruskal_h(groups):
    values = [
        float(value)
        for group in groups
        for value in group
    ]

    n = len(values)

    if n < 2:
        raise ValueError(
            "insufficient values for KW"
        )

    ranks, ties = average_ranks(
        values
    )

    offset = 0
    total = 0.0

    for group in groups:
        size = len(group)

        if size == 0:
            raise ValueError(
                "empty KW group"
            )

        rank_sum = sum(
            ranks[
                offset:
                offset + size
            ]
        )

        total += (
            rank_sum
            * rank_sum
            / size
        )

        offset += size

    h = (
        12.0
        / (
            n
            * (
                n + 1
            )
        )
        * total
        - 3.0
        * (
            n + 1
        )
    )

    tie_term = sum(
        size ** 3 - size
        for size in ties
        if size > 1
    )

    correction = (
        1.0
        - tie_term
        / (
            n ** 3 - n
        )
    )

    if correction <= 0.0:
        raise ValueError(
            "nonpositive KW tie correction"
        )

    return h / correction


def bh(pvalues):
    n = len(pvalues)

    order = sorted(
        range(n),
        key=lambda index: (
            pvalues[index],
            index,
        ),
    )

    ranked_adjusted = [
        0.0
        for _ in range(n)
    ]

    running = 1.0

    for position in range(
        n - 1,
        -1,
        -1,
    ):
        index = order[
            position
        ]

        rank = position + 1

        value = min(
            1.0,
            pvalues[index]
            * n
            / rank,
        )

        running = min(
            running,
            value,
        )

        ranked_adjusted[
            position
        ] = running

    output = [
        0.0
        for _ in range(n)
    ]

    for position, index in enumerate(
        order
    ):
        output[
            index
        ] = ranked_adjusted[
            position
        ]

    return output


def derived_seed(
    analysis_id,
    test_id,
):
    payload = (
        str(BASE_SEED)
        + "|"
        + analysis_id
        + "|"
        + test_id
    ).encode(
        "utf-8"
    )

    digest = hashlib.sha256(
        payload
    ).digest()

    return int.from_bytes(
        digest[:8],
        byteorder="big",
        signed=False,
    )


def parse_allele_key(key):
    parts = key.split(
        ":",
        3,
    )

    if len(parts) != 4:
        raise Refusal(
            "unparseable allele key "
            + key
        )

    chrom, pos, ref, alt = parts

    try:
        pos_int = int(pos)
    except ValueError as exc:
        raise Refusal(
            "noninteger allele position "
            + key
        ) from exc

    return (
        chrom,
        pos_int,
        ref,
        alt,
    )


def load_metadata(path):
    fields, rows = read_tsv(
        path
    )

    required = {
        "s",
        "population",
        "population_inference.pop",
        "hgdp_tgp_meta.Genetic.region",
    }

    if not required.issubset(
        set(fields)
    ):
        raise Refusal(
            "metadata schema mismatch"
        )

    metadata = {}

    for row in rows:
        sample = row[
            "s"
        ].strip()

        if not sample:
            continue

        if sample in metadata:
            raise Refusal(
                "duplicate metadata sample "
                + sample
            )

        metadata[
            sample
        ] = row

    if len(metadata) != 4150:
        raise Refusal(
            "metadata samples !=4150"
        )

    return metadata


def read_order(path):
    values = [
        line.strip()
        for line in Path(path).read_text(
            encoding="utf-8"
        ).splitlines()
        if line.strip()
    ]

    if len(values) != len(
        set(values)
    ):
        raise Refusal(
            "sample order contains duplicates: "
            + str(path)
        )

    return values


def load_module_map(path):
    _, rows = read_tsv(
        path
    )

    mapping = {}

    for row in rows:
        gene = row[
            "gene_symbol"
        ].strip()

        module = row[
            "module_id"
        ].strip()

        if gene in mapping:
            raise Refusal(
                "duplicate module-map gene "
                + gene
            )

        mapping[
            gene
        ] = module

    if len(mapping) != 24:
        raise Refusal(
            "module map genes !=24"
        )

    counts = Counter(
        mapping.values()
    )

    if counts != {
        "M1": 5,
        "M2": 7,
        "M3": 5,
        "M4": 7,
    }:
        raise Refusal(
            "module sizes mismatch"
        )

    return mapping


def group_counts(
    metadata,
    samples,
    field,
):
    counts = Counter()

    for sample in samples:
        if sample not in metadata:
            raise Refusal(
                "sample absent from metadata "
                + sample
            )

        value = metadata[
            sample
        ][field].strip()

        if not value:
            raise Refusal(
                "blank grouping value "
                + field
                + " "
                + sample
            )

        counts[
            value
        ] += 1

    return counts


def eligible_set(
    counts,
    threshold,
):
    return {
        group
        for group, n
        in counts.items()
        if n >= threshold
    }


def group_table_effects(
    schema,
    eligible_groups,
):
    fields, rows = read_tsv(
        schema["path"]
    )

    gene_field = schema[
        "gene_field"
    ]

    key_field = schema[
        "allele_key_field"
    ]

    group_field = schema[
        "group_field"
    ]

    an_field = schema[
        "an_field"
    ]

    af_field = schema[
        "af_field"
    ]

    by_key = {}

    expected_rows = int(
        schema[
            "expected_rows"
        ]
    )

    if len(rows) != expected_rows:
        raise Refusal(
            schema["surface"]
            + " row-count drift"
        )

    for row in rows:
        group = row[
            group_field
        ].strip()

        if group not in eligible_groups:
            continue

        key = row[
            key_field
        ].strip()

        gene = row[
            gene_field
        ].strip()

        try:
            an = int(
                row[
                    an_field
                ]
            )
        except ValueError:
            continue

        if an <= 0:
            continue

        try:
            af = float(
                row[
                    af_field
                ]
            )
        except ValueError:
            continue

        if not math.isfinite(
            af
        ):
            continue

        if key not in by_key:
            by_key[
                key
            ] = {
                "gene": gene,
                "values": [],
            }

        elif by_key[key][
            "gene"
        ] != gene:
            raise Refusal(
                "gene disagreement for "
                + key
            )

        by_key[
            key
        ][
            "values"
        ].append(
            (
                group,
                af,
            )
        )

    effects = {}

    for key, payload in by_key.items():
        values = payload[
            "values"
        ]

        if len(values) < 2:
            continue

        minimum = min(
            values,
            key=lambda item: (
                item[1],
                item[0],
            ),
        )

        maximum = max(
            values,
            key=lambda item: (
                item[1],
                item[0],
            ),
        )

        effects[
            key
        ] = {
            "gene": payload[
                "gene"
            ],
            "delta_AF":
                maximum[1]
                - minimum[1],
            "valid_group_count":
                len(values),
            "min_group":
                minimum[0],
            "min_AF":
                minimum[1],
            "max_group":
                maximum[0],
            "max_AF":
                maximum[1],
        }

    return effects


def combine_effects(*effect_maps):
    combined = {}

    for effect_map in effect_maps:
        for key, payload in effect_map.items():
            if key in combined:
                raise Refusal(
                    "duplicate allele key across analysis layers "
                    + key
                )

            combined[
                key
            ] = payload

    return combined


def gene_scores_from_effects(
    effects,
    genes,
    method="median",
):
    by_gene = defaultdict(
        list
    )

    for payload in effects.values():
        by_gene[
            payload["gene"]
        ].append(
            payload[
                "delta_AF"
            ]
        )

    output = {}

    for gene in sorted(
        genes
    ):
        values = by_gene.get(
            gene,
            [],
        )

        if not values:
            output[
                gene
            ] = {
                "score": None,
                "valid_variant_count": 0,
                "evaluable": False,
            }
            continue

        if method == "median":
            score = median(
                values
            )
        elif method == "p90":
            score = quantile_linear(
                values,
                0.9,
            )
        else:
            raise ValueError(
                "unknown gene score method "
                + method
            )

        output[
            gene
        ] = {
            "score": score,
            "valid_variant_count":
                len(values),
            "evaluable": True,
        }

    return output


def infer_modules(
    analysis_id,
    gene_scores,
    module_map,
    remove_genes=None,
):
    remove_genes = set(
        remove_genes or []
    )

    required_genes = sorted(
        set(module_map)
        - remove_genes
    )

    missing = [
        gene
        for gene in required_genes
        if (
            gene not in gene_scores
            or not gene_scores[
                gene
            ][
                "evaluable"
            ]
        )
    ]

    module_rows = []
    global_rows = []
    pairwise_rows = []

    if missing:
        for module in (
            "M1",
            "M2",
            "M3",
            "M4",
        ):
            genes = [
                gene
                for gene in required_genes
                if module_map[
                    gene
                ] == module
            ]

            evaluable = [
                gene
                for gene in genes
                if gene not in missing
            ]

            values = [
                gene_scores[
                    gene
                ][
                    "score"
                ]
                for gene in evaluable
            ]

            module_rows.append(
                {
                    "analysis_id":
                        analysis_id,
                    "module_id":
                        module,
                    "required_gene_count":
                        len(genes),
                    "evaluable_gene_count":
                        len(evaluable),
                    "module_median_gene_score":
                        (
                            ""
                            if not values
                            else format(
                                median(
                                    values
                                ),
                                ".17g",
                            )
                        ),
                    "inference_status":
                        "NON_EVALUABLE",
                }
            )

        global_rows.append(
            {
                "analysis_id":
                    analysis_id,
                "observed_H":
                    "",
                "permutations":
                    PERMUTATIONS,
                "seed":
                    "",
                "extreme_or_equal":
                    "",
                "permutation_p":
                    "",
                "status":
                    "SKIPPED_NON_EVALUABLE_GENES",
                "missing_genes":
                    ",".join(
                        missing
                    ),
            }
        )

        for left, right in itertools.combinations(
            (
                "M1",
                "M2",
                "M3",
                "M4",
            ),
            2,
        ):
            pairwise_rows.append(
                {
                    "analysis_id":
                        analysis_id,
                    "module_left":
                        left,
                    "module_right":
                        right,
                    "observed_absolute_median_difference":
                        "",
                    "permutations":
                        PERMUTATIONS,
                    "seed":
                        "",
                    "extreme_or_equal":
                        "",
                    "raw_permutation_p":
                        "",
                    "BH_adjusted_p":
                        "",
                    "status":
                        "SKIPPED_NON_EVALUABLE_GENES",
                }
            )

        return (
            module_rows,
            global_rows,
            pairwise_rows,
            False,
        )

    module_values = {}

    for module in (
        "M1",
        "M2",
        "M3",
        "M4",
    ):
        genes = [
            gene
            for gene in required_genes
            if module_map[
                gene
            ] == module
        ]

        values = [
            gene_scores[
                gene
            ][
                "score"
            ]
            for gene in genes
        ]

        module_values[
            module
        ] = values

        module_rows.append(
            {
                "analysis_id":
                    analysis_id,
                "module_id":
                    module,
                "required_gene_count":
                    len(genes),
                "evaluable_gene_count":
                    len(genes),
                "module_median_gene_score":
                    format(
                        median(
                            values
                        ),
                        ".17g",
                    ),
                "inference_status":
                    "EVALUABLE",
            }
        )

    ordered_modules = (
        "M1",
        "M2",
        "M3",
        "M4",
    )

    groups = [
        module_values[
            module
        ]
        for module in ordered_modules
    ]

    observed_h = kruskal_h(
        groups
    )

    pooled = [
        value
        for group in groups
        for value in group
    ]

    sizes = [
        len(group)
        for group in groups
    ]

    seed = derived_seed(
        analysis_id,
        "GLOBAL",
    )

    rng = random.Random(
        seed
    )

    extreme = 0

    work = list(
        pooled
    )

    for _ in range(
        PERMUTATIONS
    ):
        rng.shuffle(
            work
        )

        offset = 0
        perm_groups = []

        for size in sizes:
            perm_groups.append(
                list(
                    work[
                        offset:
                        offset + size
                    ]
                )
            )
            offset += size

        statistic = kruskal_h(
            perm_groups
        )

        if statistic >= (
            observed_h
            - 1e-15
        ):
            extreme += 1

    global_p = (
        1
        + extreme
    ) / (
        PERMUTATIONS
        + 1
    )

    global_rows.append(
        {
            "analysis_id":
                analysis_id,
            "observed_H":
                format(
                    observed_h,
                    ".17g",
                ),
            "permutations":
                PERMUTATIONS,
            "seed":
                seed,
            "extreme_or_equal":
                extreme,
            "permutation_p":
                format(
                    global_p,
                    ".17g",
                ),
            "status":
                "EVALUABLE",
            "missing_genes":
                "",
        }
    )

    pair_temp = []

    for left, right in itertools.combinations(
        ordered_modules,
        2,
    ):
        left_values = module_values[
            left
        ]

        right_values = module_values[
            right
        ]

        observed = abs(
            median(
                left_values
            )
            - median(
                right_values
            )
        )

        pooled_pair = (
            list(
                left_values
            )
            + list(
                right_values
            )
        )

        left_n = len(
            left_values
        )

        seed_pair = derived_seed(
            analysis_id,
            left
            + "_VS_"
            + right,
        )

        rng_pair = random.Random(
            seed_pair
        )

        extreme_pair = 0

        work_pair = list(
            pooled_pair
        )

        for _ in range(
            PERMUTATIONS
        ):
            rng_pair.shuffle(
                work_pair
            )

            statistic = abs(
                median(
                    work_pair[
                        :left_n
                    ]
                )
                - median(
                    work_pair[
                        left_n:
                    ]
                )
            )

            if statistic >= (
                observed
                - 1e-15
            ):
                extreme_pair += 1

        raw_p = (
            1
            + extreme_pair
        ) / (
            PERMUTATIONS
            + 1
        )

        pair_temp.append(
            {
                "analysis_id":
                    analysis_id,
                "module_left":
                    left,
                "module_right":
                    right,
                "observed_absolute_median_difference":
                    format(
                        observed,
                        ".17g",
                    ),
                "permutations":
                    PERMUTATIONS,
                "seed":
                    seed_pair,
                "extreme_or_equal":
                    extreme_pair,
                "raw_permutation_p":
                    raw_p,
                "status":
                    "EVALUABLE",
            }
        )

    adjusted = bh(
        [
            row[
                "raw_permutation_p"
            ]
            for row in pair_temp
        ]
    )

    for row, adjusted_p in zip(
        pair_temp,
        adjusted,
    ):
        row[
            "raw_permutation_p"
        ] = format(
            row[
                "raw_permutation_p"
            ],
            ".17g",
        )

        row[
            "BH_adjusted_p"
        ] = format(
            adjusted_p,
            ".17g",
        )

        pairwise_rows.append(
            row
        )

    return (
        module_rows,
        global_rows,
        pairwise_rows,
        True,
    )


def parse_gt_tokens(gt):
    gt = gt.strip()

    if (
        not gt
        or gt == "."
        or gt == "./."
        or gt == ".|."
    ):
        return None

    tokens = re.split(
        r"[|/]",
        gt,
    )

    if (
        not tokens
        or any(
            token == "."
            or token == ""
            for token in tokens
        )
    ):
        return None

    try:
        return [
            int(token)
            for token in tokens
        ]
    except ValueError:
        return None


def ploidy_adjusted_alleles(
    gt,
    karyotype,
    chrX=False,
):
    tokens = parse_gt_tokens(
        gt
    )

    if tokens is None:
        return None

    if not chrX:
        if len(tokens) != 2:
            return None

        return tokens

    if karyotype == "XX":
        if len(tokens) != 2:
            return None

        return tokens

    if karyotype == "XY":
        if len(tokens) == 1:
            return tokens

        if (
            len(tokens) == 2
            and tokens[0]
                == tokens[1]
        ):
            return [
                tokens[0]
            ]

        return None

    return None


def read_phase_mask(path):
    _, rows = read_tsv(
        path
    )

    by_gene = defaultdict(
        set
    )

    all_keys = set()

    for row in rows:
        gene = row[
            "locus_scope"
        ].strip()

        key = row[
            "allele_key"
        ].strip()

        by_gene[
            gene
        ].add(
            key
        )

        all_keys.add(
            key
        )

    return by_gene, all_keys


def phased_effects(
    mask_path,
    sources,
    order_path,
    metadata,
    population_eligible,
):
    by_gene, mask_keys = read_phase_mask(
        mask_path
    )

    sample_order = read_order(
        order_path
    )

    if len(
        sample_order
    ) != 4091:
        raise Refusal(
            "phased execution order !=4091"
        )

    sample_groups = [
        metadata[
            sample
        ][
            "population"
        ].strip()
        for sample in sample_order
    ]

    sample_karyotypes = [
        metadata[
            sample
        ].get(
            "sex_imputation.sex_karyotype",
            "",
        ).strip()
        for sample in sample_order
    ]

    effects = {}

    for gene in sorted(
        by_gene
    ):
        if gene not in sources:
            raise Refusal(
                "no genotype source for "
                + gene
            )

        target_keys = by_gene[
            gene
        ]

        positions = sorted({
            (
                parse_allele_key(
                    key
                )[0],
                parse_allele_key(
                    key
                )[1],
            )
            for key in target_keys
        })

        region_text = ",".join(
            chrom
            + ":"
            + str(pos)
            + "-"
            + str(pos)
            for chrom, pos
            in positions
        )

        command = [
            "bcftools",
            "query",
            "-S",
            str(order_path),
            "-r",
            region_text,
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT"
            "[\t%GT]\n",
            str(
                sources[
                    gene
                ]
            ),
        ]

        proc = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        found = set()

        for line in proc.stdout.splitlines():
            if not line:
                continue

            fields = line.split(
                "\t"
            )

            if len(fields) != (
                4
                + len(
                    sample_order
                )
            ):
                raise Refusal(
                    gene
                    + " GT query width mismatch"
                )

            chrom = fields[0]
            pos = fields[1]
            ref = fields[2]
            alts = fields[3].split(
                ","
            )

            gts = fields[
                4:
            ]

            chrX = (
                gene == "IKBKG"
            )

            per_alt = {}

            for alt_index, alt in enumerate(
                alts,
                1,
            ):
                key = (
                    chrom
                    + ":"
                    + pos
                    + ":"
                    + ref
                    + ":"
                    + alt
                )

                if key not in target_keys:
                    continue

                per_alt[
                    alt_index
                ] = {
                    "key": key,
                    "AC":
                        Counter(),
                    "AN":
                        Counter(),
                }

            if not per_alt:
                continue

            for index, gt in enumerate(
                gts
            ):
                group = sample_groups[
                    index
                ]

                if group not in population_eligible:
                    continue

                alleles = ploidy_adjusted_alleles(
                    gt,
                    sample_karyotypes[
                        index
                    ],
                    chrX=chrX,
                )

                if alleles is None:
                    continue

                for alt_index, payload in per_alt.items():
                    payload[
                        "AN"
                    ][
                        group
                    ] += len(
                        alleles
                    )

                    payload[
                        "AC"
                    ][
                        group
                    ] += sum(
                        allele
                        == alt_index
                        for allele in alleles
                    )

            for payload in per_alt.values():
                key = payload[
                    "key"
                ]

                values = []

                for group in sorted(
                    population_eligible
                ):
                    an = payload[
                        "AN"
                    ][
                        group
                    ]

                    if an <= 0:
                        continue

                    ac = payload[
                        "AC"
                    ][
                        group
                    ]

                    af = ac / an

                    values.append(
                        (
                            group,
                            af,
                        )
                    )

                found.add(
                    key
                )

                if len(values) < 2:
                    continue

                minimum = min(
                    values,
                    key=lambda item: (
                        item[1],
                        item[0],
                    ),
                )

                maximum = max(
                    values,
                    key=lambda item: (
                        item[1],
                        item[0],
                    ),
                )

                effects[
                    key
                ] = {
                    "gene": gene,
                    "delta_AF":
                        maximum[1]
                        - minimum[1],
                    "valid_group_count":
                        len(values),
                    "min_group":
                        minimum[0],
                    "min_AF":
                        minimum[1],
                    "max_group":
                        maximum[0],
                    "max_AF":
                        maximum[1],
                }

        missing_raw = (
            target_keys
            - found
        )

        if missing_raw:
            raise Refusal(
                gene
                + " phased mask keys absent from queried genotype source; example="
                + sorted(
                    missing_raw
                )[0]
            )

    if set(
        effects
    ) - mask_keys:
        raise Refusal(
            "phased effect key outside mask"
        )

    return effects


def load_global_af(
    auto_path,
    chrx_path,
):
    result = {}

    for path in (
        auto_path,
        chrx_path,
    ):
        _, rows = read_tsv(
            path
        )

        for row in rows:
            key = row[
                "allele_key"
            ].strip()

            if key in result:
                raise Refusal(
                    "duplicate global-AF key "
                    + key
                )

            result[
                key
            ] = float(
                row[
                    "AF"
                ]
            )

    if len(result) != 39821:
        raise Refusal(
            "global AF key universe !=39821"
        )

    return result


def parse_source_spans(path):
    _, rows = read_tsv(
        path
    )

    spans = {}

    for row in rows:
        gene = row[
            "gene_symbol"
        ].strip()

        region = row[
            "region_1based_inclusive"
        ].strip()

        match = re.search(
            r"(\d+)\s*-\s*(\d+)$",
            region,
        )

        if not match:
            raise Refusal(
                "unparseable extraction span "
                + gene
            )

        start = int(
            match.group(1)
        )

        end = int(
            match.group(2)
        )

        if end < start:
            raise Refusal(
                "invalid extraction span "
                + gene
            )

        spans[
            gene
        ] = (
            end
            - start
            + 1
        )

    if len(
        spans
    ) != 24:
        raise Refusal(
            "extraction spans !=24"
        )

    return spans


def genotype_state_rows(
    prioritized,
    sources,
    metadata,
    auto_order,
    chrx_order,
):
    output = []

    for gene in sorted(
        prioritized
    ):
        payload = prioritized[
            gene
        ]

        key = payload[
            "allele_key"
        ]

        chrom, pos, ref, alt = parse_allele_key(
            key
        )

        source = sources[
            gene
        ]

        chrX = (
            gene == "IKBKG"
        )

        order_path = (
            chrx_order
            if chrX
            else auto_order
        )

        sample_order = read_order(
            order_path
        )

        proc = subprocess.run(
            [
                "bcftools",
                "query",
                "-S",
                str(
                    order_path
                ),
                "-r",
                chrom
                + ":"
                + str(pos)
                + "-"
                + str(pos),
                "-f",
                "%CHROM\t%POS\t%REF\t%ALT"
                "[\t%GT]\n",
                str(
                    source
                ),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        target_gts = None
        target_index = None

        for line in proc.stdout.splitlines():
            fields = line.split(
                "\t"
            )

            if len(fields) < 4:
                continue

            if (
                fields[0] != chrom
                or int(
                    fields[1]
                ) != pos
                or fields[2] != ref
            ):
                continue

            alts = fields[
                3
            ].split(
                ","
            )

            if alt not in alts:
                continue

            target_index = (
                alts.index(
                    alt
                )
                + 1
            )

            target_gts = fields[
                4:
            ]

            break

        if target_gts is None:
            raise Refusal(
                "prioritized genotype-state key not recovered "
                + key
            )

        if len(
            target_gts
        ) != len(
            sample_order
        ):
            raise Refusal(
                "prioritized genotype-state query width mismatch "
                + gene
            )

        strata_counts = defaultdict(
            Counter
        )

        for sample, gt in zip(
            sample_order,
            target_gts,
        ):
            karyotype = metadata[
                sample
            ].get(
                "sex_imputation.sex_karyotype",
                "",
            ).strip()

            alleles = ploidy_adjusted_alleles(
                gt,
                karyotype,
                chrX=chrX,
            )

            if not chrX:
                stratum = "AUTOSOMAL"

                if alleles is None:
                    state = "MISSING_OR_UNAVAILABLE"

                elif alleles == [
                    0,
                    0,
                ]:
                    state = "REFERENCE_HOMOZYGOUS"

                elif alleles == [
                    target_index,
                    target_index,
                ]:
                    state = "TARGET_ALT_HOMOZYGOUS"

                elif (
                    len(alleles) == 2
                    and alleles.count(
                        target_index
                    ) == 1
                    and alleles.count(
                        0
                    ) == 1
                ):
                    state = "TARGET_HETEROZYGOUS"

                else:
                    state = "OTHER_OR_MULTIALLELIC"

            else:
                stratum = (
                    karyotype
                    if karyotype
                    else "UNKNOWN"
                )

                if alleles is None:
                    state = "MISSING_OR_UNAVAILABLE"

                elif karyotype == "XX":
                    if alleles == [
                        0,
                        0,
                    ]:
                        state = "REFERENCE_HOMOZYGOUS"

                    elif alleles == [
                        target_index,
                        target_index,
                    ]:
                        state = "TARGET_ALT_HOMOZYGOUS"

                    elif (
                        len(alleles) == 2
                        and alleles.count(
                            target_index
                        ) == 1
                        and alleles.count(
                            0
                        ) == 1
                    ):
                        state = "TARGET_HETEROZYGOUS"

                    else:
                        state = "OTHER_OR_MULTIALLELIC"

                elif karyotype == "XY":
                    if alleles == [
                        0
                    ]:
                        state = "REFERENCE_HEMIZYGOUS"

                    elif alleles == [
                        target_index
                    ]:
                        state = "TARGET_ALT_HEMIZYGOUS"

                    else:
                        state = "OTHER_OR_MULTIALLELIC"

                else:
                    state = "MISSING_OR_UNAVAILABLE"

            strata_counts[
                stratum
            ][
                state
            ] += 1

        for stratum in sorted(
            strata_counts
        ):
            total = sum(
                strata_counts[
                    stratum
                ].values()
            )

            for state in sorted(
                strata_counts[
                    stratum
                ]
            ):
                count = strata_counts[
                    stratum
                ][
                    state
                ]

                output.append(
                    {
                        "gene_symbol":
                            gene,
                        "allele_key":
                            key,
                        "stratum":
                            stratum,
                        "genotype_state":
                            state,
                        "count":
                            count,
                        "stratum_total":
                            total,
                        "frequency":
                            format(
                                count
                                / total,
                                ".17g",
                            ),
                    }
                )

    return output


def synthetic_selftest():
    assert median(
        [1, 2, 3, 4]
    ) == 2.5

    assert abs(
        quantile_linear(
            [0, 10],
            0.9,
        )
        - 9.0
    ) < 1e-12

    groups = (
        [1, 2, 3],
        [2, 3, 4],
        [4, 5, 6],
        [7, 8, 9],
    )

    h = kruskal_h(
        groups
    )

    if not math.isfinite(
        h
    ):
        raise RuntimeError(
            "KW self-test failed"
        )

    adjusted = bh(
        [
            0.001,
            0.02,
            0.03,
            0.2,
            0.5,
            0.8,
        ]
    )

    if any(
        not (
            0.0
            <= value
            <= 1.0
        )
        for value in adjusted
    ):
        raise RuntimeError(
            "BH self-test failed"
        )

    seed_a = derived_seed(
        "PRIMARY",
        "GLOBAL",
    )

    seed_b = derived_seed(
        "PRIMARY",
        "GLOBAL",
    )

    if seed_a != seed_b:
        raise RuntimeError(
            "seed derivation not deterministic"
        )

    return {
        "status": "PASS",
        "median": True,
        "P90_linear": True,
        "Kruskal_Wallis_tie_corrected": True,
        "Benjamini_Hochberg": True,
        "deterministic_seed": seed_a,
    }


def verify_authorization(
    path,
    executor_path,
    project,
    results_dir,
    logs_dir,
):
    path = Path(
        path
    )

    if not path.is_file():
        raise Refusal(
            "REFUSE_PHASE_E_AUTHORIZATION_FILE_MISSING"
        )

    auth = parse_key_value_file(
        path
    )

    required = {
        "PHASE_E_EXECUTION_AUTHORIZED":
            "YES",
        "AUTHORIZATION_SINGLE_USE":
            "YES",
        "EXECUTOR_SHA256":
            sha256(
                executor_path
            ),
        "RESULTS_DIR":
            str(
                results_dir
            ),
        "LOGS_DIR":
            str(
                logs_dir
            ),
    }

    for key, expected in required.items():
        if auth.get(
            key
        ) != expected:
            raise Refusal(
                "authorization mismatch "
                + key
            )

    head = subprocess.check_output(
        [
            "git",
            "-C",
            str(
                project
            ),
            "rev-parse",
            "HEAD",
        ],
        text=True,
    ).strip()

    if auth.get(
        "AUTHORIZED_REPOSITORY_HEAD"
    ) != head:
        raise Refusal(
            "authorization repository head mismatch"
        )

    return auth


def build_eligibility_rows(
    metadata,
    auto_samples,
    chrx_samples,
    phased_samples,
):
    specifications = (
        (
            "population",
            "population",
        ),
        (
            "inferred",
            "population_inference.pop",
        ),
        (
            "broad",
            "hgdp_tgp_meta.Genetic.region",
        ),
    )

    rows = []

    for grouping, field in specifications:
        auto_counts = group_counts(
            metadata,
            auto_samples,
            field,
        )

        chrx_counts = group_counts(
            metadata,
            chrx_samples,
            field,
        )

        phased_counts = group_counts(
            metadata,
            phased_samples,
            field,
        )

        groups = sorted(
            set(auto_counts)
            | set(chrx_counts)
            | set(phased_counts)
        )

        for group in groups:
            rows.append(
                {
                    "grouping":
                        grouping,
                    "group":
                        group,
                    "autosomal_n":
                        auto_counts[
                            group
                        ],
                    "chrX_n":
                        chrx_counts[
                            group
                        ],
                    "phased_n":
                        phased_counts[
                            group
                        ],
                    "autosomal_N20":
                        int(
                            auto_counts[
                                group
                            ] >= 20
                        ),
                    "chrX_N20":
                        int(
                            chrx_counts[
                                group
                            ] >= 20
                        ),
                    "phased_N20":
                        int(
                            phased_counts[
                                group
                            ] >= 20
                        ),
                    "autosomal_N30":
                        int(
                            auto_counts[
                                group
                            ] >= 30
                        ),
                    "chrX_N30":
                        int(
                            chrx_counts[
                                group
                            ] >= 30
                        ),
                    "phased_N30":
                        int(
                            phased_counts[
                                group
                            ] >= 30
                        ),
                }
            )

    return rows


def execute_analysis(
    args,
    binding_rows,
    binding,
    schema,
    sources,
):
    results_dir = Path(
        args.results_dir
    )

    logs_dir = Path(
        args.logs_dir
    )

    if results_dir.exists():
        raise Refusal(
            "results directory already exists"
        )

    results_dir.mkdir(
        parents=True,
        exist_ok=False,
    )

    logs_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    metadata = load_metadata(
        binding[
            "metadata"
        ][
            "value"
        ]
    )

    auto_order_path = Path(
        binding[
            "autosomal_sample_order"
        ][
            "value"
        ]
    )

    chrx_order_path = Path(
        binding[
            "chrX_sample_order"
        ][
            "value"
        ]
    )

    phased_order_path = Path(
        binding[
            "interoperability_order"
        ][
            "value"
        ]
    )

    auto_samples = read_order(
        auto_order_path
    )

    chrx_samples = read_order(
        chrx_order_path
    )

    phased_samples = read_order(
        phased_order_path
    )

    if len(auto_samples) != 4150:
        raise Refusal(
            "autosomal sample order !=4150"
        )

    if len(chrx_samples) != 4120:
        raise Refusal(
            "chrX sample order !=4120"
        )

    if len(phased_samples) != 4091:
        raise Refusal(
            "phased sample order !=4091"
        )

    module_map = load_module_map(
        binding[
            "module_map"
        ][
            "value"
        ]
    )

    genes = set(
        module_map
    )

    auto_pop_counts = group_counts(
        metadata,
        auto_samples,
        "population",
    )

    chrx_pop_counts = group_counts(
        metadata,
        chrx_samples,
        "population",
    )

    auto_inf_counts = group_counts(
        metadata,
        auto_samples,
        "population_inference.pop",
    )

    chrx_inf_counts = group_counts(
        metadata,
        chrx_samples,
        "population_inference.pop",
    )

    auto_broad_counts = group_counts(
        metadata,
        auto_samples,
        "hgdp_tgp_meta.Genetic.region",
    )

    chrx_broad_counts = group_counts(
        metadata,
        chrx_samples,
        "hgdp_tgp_meta.Genetic.region",
    )

    phased_pop_counts = group_counts(
        metadata,
        phased_samples,
        "population",
    )

    primary_auto_groups = eligible_set(
        auto_pop_counts,
        20,
    )

    primary_chrx_groups = eligible_set(
        chrx_pop_counts,
        20,
    )

    N30_auto_groups = eligible_set(
        auto_pop_counts,
        30,
    )

    N30_chrx_groups = eligible_set(
        chrx_pop_counts,
        30,
    )

    inferred_auto_groups = eligible_set(
        auto_inf_counts,
        20,
    )

    inferred_chrx_groups = eligible_set(
        chrx_inf_counts,
        20,
    )

    broad_auto_groups = eligible_set(
        auto_broad_counts,
        20,
    )

    broad_chrx_groups = eligible_set(
        chrx_broad_counts,
        20,
    )

    phased_groups = eligible_set(
        phased_pop_counts,
        20,
    )

    eligibility_rows = build_eligibility_rows(
        metadata,
        auto_samples,
        chrx_samples,
        phased_samples,
    )

    write_tsv(
        results_dir
        / OUTPUT_NAMES[0],
        (
            "grouping",
            "group",
            "autosomal_n",
            "chrX_n",
            "phased_n",
            "autosomal_N20",
            "chrX_N20",
            "phased_N20",
            "autosomal_N30",
            "chrX_N30",
            "phased_N30",
        ),
        eligibility_rows,
    )

    # --------------------------------------------------------
    # Primary / grouping sensitivities from exact D8 tables.
    # --------------------------------------------------------

    primary_effects = combine_effects(
        group_table_effects(
            schema[
                "AUTO_POPULATION"
            ],
            primary_auto_groups,
        ),
        group_table_effects(
            schema[
                "CHRX_POPULATION"
            ],
            primary_chrx_groups,
        ),
    )

    N30_effects = combine_effects(
        group_table_effects(
            schema[
                "AUTO_POPULATION"
            ],
            N30_auto_groups,
        ),
        group_table_effects(
            schema[
                "CHRX_POPULATION"
            ],
            N30_chrx_groups,
        ),
    )

    inferred_effects = combine_effects(
        group_table_effects(
            schema[
                "AUTO_INFERRED"
            ],
            inferred_auto_groups,
        ),
        group_table_effects(
            schema[
                "CHRX_INFERRED"
            ],
            inferred_chrx_groups,
        ),
    )

    broad_effects = combine_effects(
        group_table_effects(
            schema[
                "AUTO_BROAD"
            ],
            broad_auto_groups,
        ),
        group_table_effects(
            schema[
                "CHRX_BROAD"
            ],
            broad_chrx_groups,
        ),
    )

    # --------------------------------------------------------
    # Full Phase-B phased masks on 4091 genotype-side cohort.
    # --------------------------------------------------------

    strict_phased_effects = phased_effects(
        binding[
            "strict_phased_inventory"
        ][
            "value"
        ],
        sources,
        phased_order_path,
        metadata,
        phased_groups,
    )

    review_phased_effects = phased_effects(
        binding[
            "review_phased_inventory"
        ][
            "value"
        ],
        sources,
        phased_order_path,
        metadata,
        phased_groups,
    )

    effect_sets = {
        "PRIMARY_N20":
            primary_effects,
        "SENS_N30":
            N30_effects,
        "SENS_INFERRED":
            inferred_effects,
        "SENS_BROAD":
            broad_effects,
        "SENS_PHASED_STRICT":
            strict_phased_effects,
        "SENS_PHASED_REVIEW":
            review_phased_effects,
    }

    # --------------------------------------------------------
    # Variant-effect output.
    # --------------------------------------------------------

    variant_rows = []

    for analysis_id in sorted(
        effect_sets
    ):
        for key in sorted(
            effect_sets[
                analysis_id
            ]
        ):
            payload = effect_sets[
                analysis_id
            ][
                key
            ]

            gene = payload[
                "gene"
            ]

            variant_rows.append(
                {
                    "analysis_id":
                        analysis_id,
                    "gene_symbol":
                        gene,
                    "module_id":
                        module_map[
                            gene
                        ],
                    "allele_key":
                        key,
                    "delta_AF":
                        format(
                            payload[
                                "delta_AF"
                            ],
                            ".17g",
                        ),
                    "valid_group_count":
                        payload[
                            "valid_group_count"
                        ],
                    "min_group":
                        payload[
                            "min_group"
                        ],
                    "min_AF":
                        format(
                            payload[
                                "min_AF"
                            ],
                            ".17g",
                        ),
                    "max_group":
                        payload[
                            "max_group"
                        ],
                    "max_AF":
                        format(
                            payload[
                                "max_AF"
                            ],
                            ".17g",
                        ),
                }
            )

    write_tsv(
        results_dir
        / OUTPUT_NAMES[1],
        (
            "analysis_id",
            "gene_symbol",
            "module_id",
            "allele_key",
            "delta_AF",
            "valid_group_count",
            "min_group",
            "min_AF",
            "max_group",
            "max_AF",
        ),
        variant_rows,
    )

    # --------------------------------------------------------
    # Inferential gene scores.
    # --------------------------------------------------------

    gene_score_sets = {
        "PRIMARY_N20":
            gene_scores_from_effects(
                primary_effects,
                genes,
                "median",
            ),
        "SENS_N30":
            gene_scores_from_effects(
                N30_effects,
                genes,
                "median",
            ),
        "SENS_INFERRED":
            gene_scores_from_effects(
                inferred_effects,
                genes,
                "median",
            ),
        "SENS_BROAD":
            gene_scores_from_effects(
                broad_effects,
                genes,
                "median",
            ),
        "SENS_P90":
            gene_scores_from_effects(
                primary_effects,
                genes,
                "p90",
            ),
        "SENS_PHASED_STRICT":
            gene_scores_from_effects(
                strict_phased_effects,
                genes,
                "median",
            ),
        "SENS_PHASED_REVIEW":
            gene_scores_from_effects(
                review_phased_effects,
                genes,
                "median",
            ),
    }

    if any(
        not payload[
            "evaluable"
        ]
        for payload in gene_score_sets[
            "PRIMARY_N20"
        ].values()
    ):
        missing = [
            gene
            for gene, payload
            in gene_score_sets[
                "PRIMARY_N20"
            ].items()
            if not payload[
                "evaluable"
            ]
        ]

        raise Refusal(
            "primary zero-valid-variant gene(s): "
            + ",".join(
                missing
            )
        )

    # --------------------------------------------------------
    # Descriptive sensitivity gene scores.
    # --------------------------------------------------------

    global_af = load_global_af(
        binding[
            "auto_global_AF"
        ][
            "value"
        ],
        binding[
            "chrx_global_AF"
        ][
            "value"
        ],
    )

    rare_effect_sets = {}

    for label, threshold in (
        (
            "DESC_RARE_LT_0.01",
            0.01,
        ),
        (
            "DESC_RARE_LT_0.001",
            0.001,
        ),
    ):
        rare_effect_sets[
            label
        ] = {
            key: payload
            for key, payload
            in primary_effects.items()
            if (
                key in global_af
                and global_af[
                    key
                ] < threshold
            )
        }

    descriptive_gene_scores = {}

    for label, effects in rare_effect_sets.items():
        descriptive_gene_scores[
            label
        ] = gene_scores_from_effects(
            effects,
            genes,
            "median",
        )

    by_gene_position = defaultdict(
        lambda: defaultdict(
            list
        )
    )

    for key, payload in primary_effects.items():
        chrom, pos, _, _ = parse_allele_key(
            key
        )

        by_gene_position[
            payload[
                "gene"
            ]
        ][
            (
                chrom,
                pos,
            )
        ].append(
            payload[
                "delta_AF"
            ]
        )

    unique_position_scores = {}

    for gene in sorted(
        genes
    ):
        positions = by_gene_position.get(
            gene,
            {}
        )

        collapsed = [
            median(
                values
            )
            for values in positions.values()
        ]

        if collapsed:
            unique_position_scores[
                gene
            ] = {
                "score":
                    median(
                        collapsed
                    ),
                "valid_variant_count":
                    len(
                        collapsed
                    ),
                "evaluable":
                    True,
            }
        else:
            unique_position_scores[
                gene
            ] = {
                "score":
                    None,
                "valid_variant_count":
                    0,
                "evaluable":
                    False,
            }

    descriptive_gene_scores[
        "DESC_UNIQUE_POSITION"
    ] = unique_position_scores

    spans = parse_source_spans(
        binding[
            "source_manifest"
        ][
            "value"
        ]
    )

    _, strict_rows = read_tsv(
        binding[
            "phaseC_strict_records"
        ][
            "value"
        ]
    )

    strict_counts = Counter(
        row[
            "gene_symbol"
        ].strip()
        for row in strict_rows
    )

    density_scores = {}

    for gene in sorted(
        genes
    ):
        span_kb = (
            spans[
                gene
            ]
            / 1000.0
        )

        density_scores[
            gene
        ] = {
            "score":
                strict_counts[
                    gene
                ]
                / span_kb,
            "valid_variant_count":
                strict_counts[
                    gene
                ],
            "evaluable":
                True,
        }

    descriptive_gene_scores[
        "DESC_VARIANT_DENSITY_PER_KB"
    ] = density_scores

    # --------------------------------------------------------
    # Gene architecture.
    # --------------------------------------------------------

    gene_rows = []

    for analysis_id in sorted(
        gene_score_sets
    ):
        method = (
            "P90_LINEAR_DELTA_AF"
            if analysis_id
                == "SENS_P90"
            else "MEDIAN_DELTA_AF"
        )

        for gene in sorted(
            genes
        ):
            payload = gene_score_sets[
                analysis_id
            ][
                gene
            ]

            gene_rows.append(
                {
                    "analysis_id":
                        analysis_id,
                    "analysis_type":
                        "INFERENTIAL",
                    "gene_symbol":
                        gene,
                    "module_id":
                        module_map[
                            gene
                        ],
                    "score_method":
                        method,
                    "score":
                        (
                            ""
                            if not payload[
                                "evaluable"
                            ]
                            else format(
                                payload[
                                    "score"
                                ],
                                ".17g",
                            )
                        ),
                    "valid_unit_count":
                        payload[
                            "valid_variant_count"
                        ],
                    "evaluable":
                        int(
                            payload[
                                "evaluable"
                            ]
                        ),
                }
            )

    for analysis_id in sorted(
        descriptive_gene_scores
    ):
        for gene in sorted(
            genes
        ):
            payload = descriptive_gene_scores[
                analysis_id
            ][
                gene
            ]

            gene_rows.append(
                {
                    "analysis_id":
                        analysis_id,
                    "analysis_type":
                        "DESCRIPTIVE",
                    "gene_symbol":
                        gene,
                    "module_id":
                        module_map[
                            gene
                        ],
                    "score_method":
                        analysis_id,
                    "score":
                        (
                            ""
                            if not payload[
                                "evaluable"
                            ]
                            else format(
                                payload[
                                    "score"
                                ],
                                ".17g",
                            )
                        ),
                    "valid_unit_count":
                        payload[
                            "valid_variant_count"
                        ],
                    "evaluable":
                        int(
                            payload[
                                "evaluable"
                            ]
                        ),
                }
            )

    write_tsv(
        results_dir
        / OUTPUT_NAMES[2],
        (
            "analysis_id",
            "analysis_type",
            "gene_symbol",
            "module_id",
            "score_method",
            "score",
            "valid_unit_count",
            "evaluable",
        ),
        gene_rows,
    )

    # --------------------------------------------------------
    # Module inference.
    # --------------------------------------------------------

    module_rows = []
    global_rows = []
    pairwise_rows = []

    inferential_order = (
        "PRIMARY_N20",
        "SENS_N30",
        "SENS_INFERRED",
        "SENS_BROAD",
        "SENS_P90",
        "SENS_PHASED_STRICT",
        "SENS_PHASED_REVIEW",
    )

    inference_status = {}

    for analysis_id in inferential_order:
        (
            module_part,
            global_part,
            pairwise_part,
            evaluable,
        ) = infer_modules(
            analysis_id,
            gene_score_sets[
                analysis_id
            ],
            module_map,
        )

        module_rows.extend(
            module_part
        )

        global_rows.extend(
            global_part
        )

        pairwise_rows.extend(
            pairwise_part
        )

        inference_status[
            analysis_id
        ] = evaluable

    # M2 excluding IKBKG sensitivity.
    (
        module_part,
        global_part,
        pairwise_part,
        evaluable,
    ) = infer_modules(
        "SENS_M2_NO_IKBKG",
        gene_score_sets[
            "PRIMARY_N20"
        ],
        module_map,
        remove_genes={
            "IKBKG"
        },
    )

    module_rows.extend(
        module_part
    )

    global_rows.extend(
        global_part
    )

    pairwise_rows.extend(
        pairwise_part
    )

    inference_status[
        "SENS_M2_NO_IKBKG"
    ] = evaluable

    write_tsv(
        results_dir
        / OUTPUT_NAMES[3],
        (
            "analysis_id",
            "module_id",
            "required_gene_count",
            "evaluable_gene_count",
            "module_median_gene_score",
            "inference_status",
        ),
        module_rows,
    )

    write_tsv(
        results_dir
        / OUTPUT_NAMES[4],
        (
            "analysis_id",
            "observed_H",
            "permutations",
            "seed",
            "extreme_or_equal",
            "permutation_p",
            "status",
            "missing_genes",
        ),
        global_rows,
    )

    write_tsv(
        results_dir
        / OUTPUT_NAMES[5],
        (
            "analysis_id",
            "module_left",
            "module_right",
            "observed_absolute_median_difference",
            "permutations",
            "seed",
            "extreme_or_equal",
            "raw_permutation_p",
            "BH_adjusted_p",
            "status",
        ),
        pairwise_rows,
    )

    # --------------------------------------------------------
    # Sensitivity summary.
    # --------------------------------------------------------

    sensitivity_rows = []

    for analysis_id in inferential_order:
        scores = gene_score_sets[
            analysis_id
        ]

        sensitivity_rows.append(
            {
                "analysis_id":
                    analysis_id,
                "analysis_type":
                    "INFERENTIAL",
                "evaluable_genes":
                    sum(
                        payload[
                            "evaluable"
                        ]
                        for payload in scores.values()
                    ),
                "required_genes":
                    24,
                "module_inference_status":
                    (
                        "EVALUABLE"
                        if inference_status[
                            analysis_id
                        ]
                        else "SKIPPED_NON_EVALUABLE"
                    ),
                "notes":
                    "",
            }
        )

    sensitivity_rows.append(
        {
            "analysis_id":
                "SENS_M2_NO_IKBKG",
            "analysis_type":
                "INFERENTIAL",
            "evaluable_genes":
                23,
            "required_genes":
                23,
            "module_inference_status":
                (
                    "EVALUABLE"
                    if inference_status[
                        "SENS_M2_NO_IKBKG"
                    ]
                    else "SKIPPED_NON_EVALUABLE"
                ),
            "notes":
                "IKBKG excluded from M2",
        }
    )

    for analysis_id in sorted(
        descriptive_gene_scores
    ):
        scores = descriptive_gene_scores[
            analysis_id
        ]

        sensitivity_rows.append(
            {
                "analysis_id":
                    analysis_id,
                "analysis_type":
                    "DESCRIPTIVE",
                "evaluable_genes":
                    sum(
                        payload[
                            "evaluable"
                        ]
                        for payload in scores.values()
                    ),
                "required_genes":
                    24,
                "module_inference_status":
                    "NOT_AUTHORIZED",
                "notes":
                    (
                        "No forced permutation inference"
                    ),
            }
        )

    write_tsv(
        results_dir
        / OUTPUT_NAMES[6],
        (
            "analysis_id",
            "analysis_type",
            "evaluable_genes",
            "required_genes",
            "module_inference_status",
            "notes",
        ),
        sensitivity_rows,
    )

    # --------------------------------------------------------
    # Full primary variant ranking.
    # --------------------------------------------------------

    ranking_rows = []

    global_sorted = sorted(
        primary_effects.items(),
        key=lambda item: (
            -item[
                1
            ][
                "delta_AF"
            ],
            item[
                0
            ],
        ),
    )

    within_gene_rank = {}

    by_gene_primary = defaultdict(
        list
    )

    for key, payload in primary_effects.items():
        by_gene_primary[
            payload[
                "gene"
            ]
        ].append(
            (
                key,
                payload,
            )
        )

    prioritized = {}

    for gene in sorted(
        genes
    ):
        ordered = sorted(
            by_gene_primary[
                gene
            ],
            key=lambda item: (
                -item[
                    1
                ][
                    "delta_AF"
                ],
                item[
                    0
                ],
            ),
        )

        if not ordered:
            raise Refusal(
                "primary ranking gene has no valid variant "
                + gene
            )

        prioritized[
            gene
        ] = {
            "allele_key":
                ordered[
                    0
                ][
                    0
                ],
            "delta_AF":
                ordered[
                    0
                ][
                    1
                ][
                    "delta_AF"
                ],
        }

        for rank, (
            key,
            payload,
        ) in enumerate(
            ordered,
            1,
        ):
            within_gene_rank[
                key
            ] = rank

    for global_rank, (
        key,
        payload,
    ) in enumerate(
        global_sorted,
        1,
    ):
        gene = payload[
            "gene"
        ]

        ranking_rows.append(
            {
                "global_rank":
                    global_rank,
                "within_gene_rank":
                    within_gene_rank[
                        key
                    ],
                "gene_symbol":
                    gene,
                "module_id":
                    module_map[
                        gene
                    ],
                "allele_key":
                    key,
                "delta_AF":
                    format(
                        payload[
                            "delta_AF"
                        ],
                        ".17g",
                    ),
                "valid_group_count":
                    payload[
                        "valid_group_count"
                    ],
                "min_group":
                    payload[
                        "min_group"
                    ],
                "min_AF":
                    format(
                        payload[
                            "min_AF"
                        ],
                        ".17g",
                    ),
                "max_group":
                    payload[
                        "max_group"
                    ],
                "max_AF":
                    format(
                        payload[
                            "max_AF"
                        ],
                        ".17g",
                    ),
            }
        )

    write_tsv(
        results_dir
        / OUTPUT_NAMES[7],
        (
            "global_rank",
            "within_gene_rank",
            "gene_symbol",
            "module_id",
            "allele_key",
            "delta_AF",
            "valid_group_count",
            "min_group",
            "min_AF",
            "max_group",
            "max_AF",
        ),
        ranking_rows,
    )

    # --------------------------------------------------------
    # Descriptive genotype states for top variant per gene.
    # --------------------------------------------------------

    genotype_rows = genotype_state_rows(
        prioritized,
        sources,
        metadata,
        auto_order_path,
        chrx_order_path,
    )

    write_tsv(
        results_dir
        / OUTPUT_NAMES[8],
        (
            "gene_symbol",
            "allele_key",
            "stratum",
            "genotype_state",
            "count",
            "stratum_total",
            "frequency",
        ),
        genotype_rows,
    )

    # --------------------------------------------------------
    # Machine-readable execution summary.
    # --------------------------------------------------------

    output_hashes = {}

    for name in OUTPUT_NAMES[
        :9
    ]:
        path = results_dir / name

        if not path.is_file():
            raise Refusal(
                "required scientific output missing "
                + name
            )

        output_hashes[
            name
        ] = {
            "bytes":
                path.stat().st_size,
            "sha256":
                sha256(
                    path
                ),
        }

    summary = {
        "status": "PASS",
        "failures": [],
        "phase": "Phase_E",
        "design": (
            "reviewer_driven_population_architecture"
        ),
        "scientific_boundaries": {
            "haplotype_inference": False,
            "LD_inference": False,
            "selection_inference": False,
            "adaptation_claims": False,
            "causal_functional_claims": False,
            "disease_association": False,
        },
        "counts": {
            "genes": 24,
            "modules": 4,
            "phaseC_strict_keys": 39821,
            "strict_phased_mask": 20793,
            "technical_review_phased_mask": 2393,
            "phased_union": 23186,
            "phased_samples": 4091,
            "primary_valid_variant_effects":
                len(
                    primary_effects
                ),
            "strict_phased_valid_variant_effects":
                len(
                    strict_phased_effects
                ),
            "review_phased_valid_variant_effects":
                len(
                    review_phased_effects
                ),
        },
        "statistics": {
            "permutations":
                PERMUTATIONS,
            "base_seed":
                BASE_SEED,
            "P90_method":
                "linear",
            "global_test":
                "tie_corrected_Kruskal_Wallis_H",
            "pairwise_multiple_testing":
                "Benjamini_Hochberg",
        },
        "output_hashes":
            output_hashes,
    }

    summary_path = (
        results_dir
        / OUTPUT_NAMES[9]
    )

    summary_path.write_text(
        json.dumps(
            summary,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    if len(
        list(
            results_dir.iterdir()
        )
    ) != 10:
        raise Refusal(
            "scientific result directory does not contain exactly 10 outputs"
        )

    return summary


def preflight(
    args,
    binding_rows,
    binding,
    schema,
    sources,
):
    if Path(
        args.results_dir
    ).exists():
        raise Refusal(
            "future results surface not pristine"
        )

    if Path(
        args.logs_dir
    ).exists():
        raise Refusal(
            "future logs surface not pristine"
        )

    _, planned = read_tsv(
        args.planned_output_contract
    )

    if len(
        planned
    ) != 10:
        raise Refusal(
            "planned output contract rows !=10"
        )

    planned_names = {
        row[
            "planned_filename"
        ].strip()
        for row in planned
    }

    if planned_names != set(
        OUTPUT_NAMES
    ):
        raise Refusal(
            "planned output filenames mismatch executor contract"
        )

    report = {
        "status": "PASS",
        "failures": [],
        "mode": "preflight",
        "phaseE_executed": False,
        "phaseE_authorized": False,
        "scientific_values_computed": False,
        "genotype_values_queried": False,
        "network_accessed": False,
        "bound_hash_files_verified":
            verify_binding_files(
                binding_rows
            ),
        "D8_schema_surfaces":
            len(
                schema
            ),
        "genotype_sources":
            len(
                sources
            ),
        "planned_outputs":
            len(
                planned
            ),
        "future_results_pristine":
            True,
        "future_logs_pristine":
            True,
        "synthetic_selftest":
            synthetic_selftest(),
    }

    return report


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=(
            "selftest",
            "preflight",
            "execute",
        ),
        required=True,
    )

    parser.add_argument(
        "--project",
    )

    parser.add_argument(
        "--binding",
    )

    parser.add_argument(
        "--schema-binding",
    )

    parser.add_argument(
        "--source-binding",
    )

    parser.add_argument(
        "--planned-output-contract",
    )

    parser.add_argument(
        "--results-dir",
    )

    parser.add_argument(
        "--logs-dir",
    )

    parser.add_argument(
        "--authorization-file",
    )

    args = parser.parse_args()

    if args.mode == "selftest":
        print(
            json.dumps(
                synthetic_selftest(),
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    required_args = (
        "project",
        "binding",
        "schema_binding",
        "source_binding",
        "planned_output_contract",
        "results_dir",
        "logs_dir",
    )

    for name in required_args:
        if getattr(
            args,
            name
        ) is None:
            raise Refusal(
                "required argument missing --"
                + name.replace(
                    "_",
                    "-",
                )
            )

    binding_rows, binding = read_binding(
        args.binding
    )

    verify_binding_files(
        binding_rows
    )

    schema = verify_schema_binding(
        args.schema_binding
    )

    sources = verify_source_binding(
        args.source_binding
    )

    if args.mode == "preflight":
        report = preflight(
            args,
            binding_rows,
            binding,
            schema,
            sources,
        )

        print(
            json.dumps(
                report,
                indent=2,
                sort_keys=True,
            )
        )

        return 0

    if not args.authorization_file:
        print(
            "REFUSE_PHASE_E_EXECUTION_WITHOUT_EXPLICIT_AUTHORIZATION_FILE"
        )
        print(
            "PHASE_E_EXECUTED=NO"
        )
        return 96

    project = Path(
        args.project
    )

    results_dir = Path(
        args.results_dir
    )

    logs_dir = Path(
        args.logs_dir
    )

    try:
        verify_authorization(
            args.authorization_file,
            Path(
                __file__
            ).resolve(),
            project,
            results_dir,
            logs_dir,
        )
    except Refusal as exc:
        print(
            "REFUSE_PHASE_E_EXECUTION_AUTHORIZATION_VALIDATION_FAILED="
            + str(
                exc
            )
        )
        print(
            "PHASE_E_EXECUTED=NO"
        )
        return 96

    summary = execute_analysis(
        args,
        binding_rows,
        binding,
        schema,
        sources,
    )

    print(
        json.dumps(
            {
                "status": "PASS",
                "phaseE_executed":
                    True,
                "result_output_count":
                    10,
                "execution_summary_sha256":
                    sha256(
                        Path(
                            args.results_dir
                        )
                        / OUTPUT_NAMES[9]
                    ),
                "summary":
                    summary,
            },
            indent=2,
            sort_keys=True,
        )
    )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(
            main()
        )
    except Refusal as exc:
        print(
            "FINAL_STATUS=REFUSE"
        )
        print(
            "REASON="
            + str(
                exc
            )
        )
        print(
            "PHASE_E_EXECUTED=NO"
        )
        raise SystemExit(
            97
        )
