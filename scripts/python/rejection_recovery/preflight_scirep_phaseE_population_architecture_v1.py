#!/usr/bin/env python3

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import sys


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


def read_binding(path):
    rows = []

    with Path(path).open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        required = {
            "kind",
            "key",
            "value",
            "sha256",
            "role",
        }

        missing = required - set(
            reader.fieldnames or []
        )

        if missing:
            raise RuntimeError(
                "binding fields missing: "
                + ",".join(sorted(missing))
            )

        rows = list(reader)

    return rows


def verify_files(rows):
    verified = 0

    for row in rows:
        if row["kind"] not in {
            "FILE",
            "BOUNDARY",
        }:
            continue

        expected = row[
            "sha256"
        ].strip()

        value = row[
            "value"
        ].strip()

        if not expected:
            continue

        path = Path(value)

        if not path.is_file():
            raise RuntimeError(
                "bound file missing: "
                + row["key"]
                + " "
                + value
            )

        actual = sha256(path)

        if actual != expected:
            raise RuntimeError(
                "bound hash mismatch: "
                + row["key"]
            )

        verified += 1

    return verified


def verify_schema(path):
    with Path(path).open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        rows = list(reader)

    if len(rows) != 6:
        raise RuntimeError(
            "D8 schema binding rows !=6"
        )

    expected_surfaces = {
        "AUTO_POPULATION",
        "AUTO_INFERRED",
        "AUTO_BROAD",
        "CHRX_POPULATION",
        "CHRX_INFERRED",
        "CHRX_BROAD",
    }

    observed = {
        row["surface"]
        for row in rows
    }

    if observed != expected_surfaces:
        raise RuntimeError(
            "D8 schema surfaces mismatch"
        )

    for row in rows:
        source = Path(
            row["path"]
        )

        if not source.is_file():
            raise RuntimeError(
                "D8 source missing "
                + row["surface"]
            )

        with source.open(
            encoding="utf-8",
            newline="",
        ) as handle:
            reader = csv.reader(
                handle,
                delimiter="\t",
            )

            header = next(
                reader
            )

        required_fields = {
            row["gene_field"],
            row["allele_key_field"],
            row["group_field"],
            row["ac_field"],
            row["an_field"],
            row["af_field"],
        }

        if not required_fields.issubset(
            set(header)
        ):
            raise RuntimeError(
                "D8 schema drift "
                + row["surface"]
            )

    return len(rows)


def verify_sources(path):
    with Path(path).open(
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        rows = list(reader)

    if len(rows) != 24:
        raise RuntimeError(
            "genotype source binding rows !=24"
        )

    genes = set()

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

        if gene in genes:
            raise RuntimeError(
                "duplicate genotype source gene "
                + gene
            )

        genes.add(gene)

        if not path.is_file():
            raise RuntimeError(
                "genotype source missing "
                + gene
            )

        if path.stat().st_size != expected_bytes:
            raise RuntimeError(
                "genotype source size mismatch "
                + gene
            )

        if sha256(path) != expected_sha:
            raise RuntimeError(
                "genotype source SHA mismatch "
                + gene
            )

    if len(genes) != 24:
        raise RuntimeError(
            "genotype source gene count !=24"
        )

    return len(rows)


def median(values):
    values = sorted(
        float(value)
        for value in values
    )

    if not values:
        raise ValueError(
            "empty median"
        )

    n = len(values)
    middle = n // 2

    if n % 2:
        return values[
            middle
        ]

    return (
        values[
            middle - 1
        ]
        + values[
            middle
        ]
    ) / 2.0


def quantile_linear(values, q):
    values = sorted(
        float(value)
        for value in values
    )

    if not values:
        raise ValueError(
            "empty quantile"
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


def bh(pvalues):
    n = len(pvalues)

    order = sorted(
        range(n),
        key=lambda i: (
            pvalues[i],
            i,
        ),
    )

    adjusted_sorted = [
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

        adjusted_sorted[
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
        ] = adjusted_sorted[
            position
        ]

    return output


def synthetic_selftest():
    if median(
        [1, 2, 3, 4]
    ) != 2.5:
        raise RuntimeError(
            "median self-test failed"
        )

    if abs(
        quantile_linear(
            [0, 10],
            0.9,
        )
        - 9.0
    ) > 1e-12:
        raise RuntimeError(
            "P90 self-test failed"
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
        value < 0.0
        or value > 1.0
        for value in adjusted
    ):
        raise RuntimeError(
            "BH self-test failed"
        )

    r1 = random.Random(
        20261011
    )

    r2 = random.Random(
        20261011
    )

    a = list(
        range(24)
    )

    b = list(
        range(24)
    )

    r1.shuffle(a)
    r2.shuffle(b)

    if a != b:
        raise RuntimeError(
            "random-seed determinism failed"
        )

    return hashlib.sha256(
        ",".join(
            str(value)
            for value in a
        ).encode("ascii")
    ).hexdigest()


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=(
            "preflight",
            "execute",
        ),
        required=True,
    )

    parser.add_argument(
        "--binding",
        required=True,
    )

    parser.add_argument(
        "--schema-binding",
        required=True,
    )

    parser.add_argument(
        "--source-binding",
        required=True,
    )

    parser.add_argument(
        "--future-results",
        required=True,
    )

    parser.add_argument(
        "--future-logs",
        required=True,
    )

    args = parser.parse_args()

    if args.mode == "execute":
        print(
            "REFUSE_EXECUTION_FP1R3_R1_PREFLIGHT_DRIVER_HAS_NO_SCIENTIFIC_EXECUTION_ENGINE"
        )
        print(
            "PHASE_E_EXECUTED=NO"
        )
        return 95

    if Path(
        args.future_results
    ).exists():
        raise RuntimeError(
            "future Phase-E result surface is not pristine"
        )

    if Path(
        args.future_logs
    ).exists():
        raise RuntimeError(
            "future Phase-E log surface is not pristine"
        )

    rows = read_binding(
        args.binding
    )

    verified_files = verify_files(
        rows
    )

    schema_rows = verify_schema(
        args.schema_binding
    )

    source_rows = verify_sources(
        args.source_binding
    )

    random_digest = synthetic_selftest()

    try:
        bcftools_version = subprocess.check_output(
            [
                "bcftools",
                "--version",
            ],
            text=True,
        ).splitlines()[0]
    except Exception as exc:
        raise RuntimeError(
            "bcftools unavailable"
        ) from exc

    report = {
        "gate": "FP1R3_R1",
        "status": "PASS",
        "failures": [],
        "mode": "preflight",
        "phaseE_executed": False,
        "phaseE_authorized": False,
        "final_phaseE_executor_present": False,
        "scientific_values_computed": False,
        "genotype_values_queried": False,
        "network_accessed": False,
        "bound_files_verified": verified_files,
        "D8_schema_surfaces_verified": schema_rows,
        "genotype_sources_verified": source_rows,
        "python_version": sys.version.split()[0],
        "bcftools_version": bcftools_version,
        "permutation_base_seed": 20261011,
        "permutation_count": 100000,
        "p90_quantile_method": "linear",
        "synthetic_random_permutation_sha256": random_digest,
        "future_results_pristine": True,
        "future_logs_pristine": True,
        "next_gate": (
            "FP1R3_R2_FREEZE_FINAL_PHASE_E_EXECUTOR_CODE_"
            "AND_SINGLE_USE_EXECUTION_PREFLIGHT_NO_ANALYSIS"
        ),
    }

    print(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
        )
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(
        main()
    )
