#!/usr/bin/env python3
"""Annotate exact probe matches and summarize Ensembl isoform coverage."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict

from probe_ensembl_common import (
    DEFAULT_PROBES, SCRIPT_DIR, ensure_ensembl_file, parse_ensembl_gtf, read_probes,
    write_tsv,
)


GENE_ALIASES = {"C9orf16": {"BBLN"}, "DDX58": {"RIGI"}}


MATCH_FIELDS = [
    "target_gene", "probe_id", "probe_sequence", "transcript_accession", "ensembl_transcript_id",
    "matched_gene", "ensembl_gene_id", "is_intended_gene", "genomic_seqid", "match_orientation",
    "transcript_start_0based", "transcript_end_0based", "transcript_name", "transcript_biotype",
    "transcript_source", "transcript_support_level", "tags",
]
PROBE_FIELDS = [
    "target_gene", "probe_id", "probe_sequence", "status", "intended_isoform_count",
    "intended_isoforms", "other_gene_count", "other_genes", "all_transcript_match_count",
]
GENE_FIELDS = [
    "target_gene", "probe_count", "annotated_ensembl_isoform_count", "matched_probe_count",
    "unmatched_probe_count", "hit_isoform_count", "hit_isoforms", "common_to_all_probes_count",
    "common_to_all_probes", "distinct_probe_isoform_sets", "classification",
]


def classify(signatures):
    nonempty = [s for s in signatures if s]
    hit_union = set().union(*nonempty) if nonempty else set()
    matched_count = sum(bool(s) for s in signatures)
    if not nonempty:
        label = "no_exact_intended_match"
    elif matched_count < len(signatures):
        label = "incomplete_probe_mapping"
    elif len(hit_union) == 1:
        label = "single_isoform"
    elif len(set(signatures)) == 1:
        label = "multiple_isoforms_same_probe_coverage"
    else:
        label = "multiple_isoforms_differential_probe_coverage"
    return label, hit_union, matched_count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probes", default=DEFAULT_PROBES)
    parser.add_argument("--alignments", default=SCRIPT_DIR / "results_ensembl/probe_transcript_alignments.tsv")
    parser.add_argument("--gtf", help="Ensembl GTF; release 116 is downloaded by default")
    parser.add_argument("--matches", default=SCRIPT_DIR / "results_ensembl/probe_isoform_matches.tsv")
    parser.add_argument("--probe-summary", default=SCRIPT_DIR / "results_ensembl/probe_isoform_summary.tsv")
    parser.add_argument("--gene-summary", default=SCRIPT_DIR / "results_ensembl/gene_isoform_summary.tsv")
    args = parser.parse_args()

    probes = read_probes(args.probes)
    gtf = args.gtf or ensure_ensembl_file("gtf")
    metadata, _ = parse_ensembl_gtf(gtf)
    aliases = GENE_ALIASES
    annotated_by_gene = defaultdict(set)
    for accession, meta in metadata.items():
        if meta["gene"]:
            annotated_by_gene[meta["gene"]].add(accession)
            for input_gene, annotation_names in aliases.items():
                if meta["gene"] in annotation_names:
                    annotated_by_gene[input_gene].add(accession)

    with open(args.alignments, encoding="utf-8", newline="") as handle:
        raw = list(csv.DictReader(handle, delimiter="\t"))

    match_rows = []
    intended_by_probe = defaultdict(set)
    all_by_probe = defaultdict(set)
    other_genes_by_probe = defaultdict(set)
    for row in raw:
        accession = row["transcript_accession"]
        meta = metadata.get(accession, {
            "gene": "", "gene_id": "", "transcript_name": "", "biotype": "",
            "source": "", "tags": "", "transcript_support_level": "",
        })
        intended_names = {row["target_gene"], *aliases.get(row["target_gene"], set())}
        intended = meta["gene"].casefold() in {name.casefold() for name in intended_names}
        all_by_probe[row["probe_id"]].add(accession)
        if intended:
            intended_by_probe[row["probe_id"]].add(accession)
        elif meta["gene"]:
            other_genes_by_probe[row["probe_id"]].add(meta["gene"])
        match_rows.append({
            **row, "ensembl_transcript_id": accession.split(".", 1)[0],
            "matched_gene": meta["gene"], "ensembl_gene_id": meta["gene_id"],
            "is_intended_gene": "yes" if intended else "no",
            "transcript_name": meta["transcript_name"], "transcript_biotype": meta["biotype"],
            "transcript_source": meta["source"],
            "transcript_support_level": meta["transcript_support_level"], "tags": meta["tags"],
        })

    write_tsv(args.matches, MATCH_FIELDS, match_rows)

    probes_by_gene = defaultdict(list)
    probe_rows = []
    for probe in probes:
        probes_by_gene[probe.target].append(probe)
        intended = intended_by_probe[probe.probe_id]
        all_matches = all_by_probe[probe.probe_id]
        others = other_genes_by_probe[probe.probe_id]
        status = "matches_intended_gene" if intended else ("off_target_only" if all_matches else "no_exact_match")
        probe_rows.append({
            "target_gene": probe.target, "probe_id": probe.probe_id, "probe_sequence": probe.sequence,
            "status": status, "intended_isoform_count": len(intended),
            "intended_isoforms": ",".join(sorted(intended)), "other_gene_count": len(others),
            "other_genes": ",".join(sorted(others)), "all_transcript_match_count": len(all_matches),
        })
    write_tsv(args.probe_summary, PROBE_FIELDS, probe_rows)

    gene_rows = []
    for gene in sorted(probes_by_gene):
        gene_probes = probes_by_gene[gene]
        signatures = [frozenset(intended_by_probe[p.probe_id]) for p in gene_probes]
        label, hit_union, matched_count = classify(signatures)
        common_all = set.intersection(*(set(s) for s in signatures)) if signatures else set()
        gene_rows.append({
            "target_gene": gene, "probe_count": len(gene_probes),
            "annotated_ensembl_isoform_count": len(annotated_by_gene.get(gene, set())),
            "matched_probe_count": matched_count, "unmatched_probe_count": len(gene_probes) - matched_count,
            "hit_isoform_count": len(hit_union), "hit_isoforms": ",".join(sorted(hit_union)),
            "common_to_all_probes_count": len(common_all),
            "common_to_all_probes": ",".join(sorted(common_all)),
            "distinct_probe_isoform_sets": len(set(signatures)), "classification": label,
        })
    write_tsv(args.gene_summary, GENE_FIELDS, gene_rows)
    print(f"Annotated {len(match_rows)} exact probe/transcript matches; "
          f"wrote {len(probe_rows)} probe summaries and {len(gene_rows)} gene summaries")


if __name__ == "__main__":
    main()
