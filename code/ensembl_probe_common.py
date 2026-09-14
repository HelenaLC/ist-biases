#!/usr/bin/env python3
"""Utilities for exact probe mapping to Ensembl transcript annotations."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

from probe_refseq_common import (
    Probe, TranscriptLocus, fasta_records, reverse_complement,
)


ENSEMBL_TRANSCRIPT_RE = re.compile(r"^(ENST[0-9]+(?:\.[0-9]+)?)\s")
GTF_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)";')


def parse_ensembl_header(header: str):
    match = ENSEMBL_TRANSCRIPT_RE.search(header)
    if not match:
        return None
    accession = match.group(1)
    genomic_seqid = ""
    fields = header.split()
    if len(fields) >= 3:
        location = fields[2].split(":")
        if len(location) >= 6:
            genomic_seqid = location[2]
    return accession, genomic_seqid


def map_probes_to_ensembl(probes: list[Probe], fasta_paths: list[str | Path], seed_length: int = 15):
    seed_index: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for i, probe in enumerate(probes):
        for orientation, pattern in (("as_provided", probe.sequence),
                                     ("reverse_complement", reverse_complement(probe.sequence))):
            seed_index[pattern[:seed_length]].append((i, orientation, pattern))

    seen: set[tuple[int, str, str, str, int]] = set()
    for fasta_path in fasta_paths:
        for header, sequence in fasta_records(fasta_path):
            parsed = parse_ensembl_header(header)
            if not parsed:
                continue
            accession, genomic_seqid = parsed
            for pos in range(max(0, len(sequence) - seed_length + 1)):
                candidates = seed_index.get(sequence[pos:pos + seed_length])
                if not candidates:
                    continue
                for probe_index, orientation, pattern in candidates:
                    if sequence.startswith(pattern, pos):
                        key = (probe_index, accession, genomic_seqid, orientation, pos)
                        if key not in seen:
                            seen.add(key)
                            yield key + (len(pattern),)


def parse_gtf_attributes(text: str) -> dict[str, list[str]]:
    attrs: dict[str, list[str]] = defaultdict(list)
    for key, value in GTF_ATTRIBUTE_RE.findall(text):
        attrs[key].append(value)
    return attrs


def first(attrs: dict[str, list[str]], key: str, default: str = "") -> str:
    values = attrs.get(key)
    return values[0] if values else default


def versioned_id(attrs: dict[str, list[str]], id_key: str, version_key: str) -> str:
    stable_id = first(attrs, id_key)
    version = first(attrs, version_key)
    return f"{stable_id}.{version}" if stable_id and version else stable_id


def parse_ensembl_gtf(path: str | Path):
    metadata: dict[str, dict[str, str]] = {}
    exon_groups: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    locus_attrs: dict[tuple[str, str, str], tuple[str, str, str]] = {}

    from probe_refseq_common import open_text
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon"}:
                continue
            seqid, source, feature, start, end, _, strand, _, attr_text = fields
            attrs = parse_gtf_attributes(attr_text)
            accession = versioned_id(attrs, "transcript_id", "transcript_version")
            if not accession:
                continue
            gene = first(attrs, "gene_name")
            gene_id = versioned_id(attrs, "gene_id", "gene_version")
            biotype = first(attrs, "transcript_biotype", first(attrs, "gene_biotype"))
            metadata.setdefault(accession, {
                "gene": gene,
                "gene_id": gene_id,
                "transcript_name": first(attrs, "transcript_name"),
                "biotype": biotype,
                "source": first(attrs, "transcript_source", source),
                "tags": ",".join(attrs.get("tag", [])),
                "transcript_support_level": first(attrs, "transcript_support_level"),
            })
            if feature == "exon":
                key = (accession, seqid, strand)
                exon_groups[key].append((int(start) - 1, int(end)))
                locus_attrs[key] = (gene, first(attrs, "transcript_name"), biotype)

    loci: dict[str, list[TranscriptLocus]] = defaultdict(list)
    for (accession, seqid, strand), exons in exon_groups.items():
        gene, transcript_name, biotype = locus_attrs[(accession, seqid, strand)]
        ordered = sorted(exons, reverse=(strand == "-"))
        loci[accession].append(TranscriptLocus(
            accession, gene, transcript_name, biotype, seqid, strand, ordered
        ))
    return metadata, loci


def ucsc_chrom_name(seqid: str, assembly_names: dict[str, str]) -> str:
    if seqid in assembly_names:
        return assembly_names[seqid]
    if seqid == "MT":
        return "chrM"
    if seqid in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        return "chr" + seqid
    return seqid
