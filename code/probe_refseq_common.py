#!/usr/bin/env python3
"""Shared exact-mapping utilities for the RefSeq probe scripts."""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
TRANSCRIPT_RE = re.compile(r"\[transcript_id=([^\]]+)\]")
GENOMIC_SEQ_RE = re.compile(r"^lcl\|((?:NC|NW|NT|NZ)_[0-9]+\.[0-9]+)_")


@dataclass(frozen=True)
class Probe:
    target: str
    probe_id: str
    sequence: str
    barcode: str


@dataclass
class TranscriptLocus:
    accession: str
    gene: str
    product: str
    biotype: str
    seqid: str
    strand: str
    exons: list[tuple[int, int]]


def open_text(path: str | Path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, encoding="utf-8")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def read_probes(path: str | Path) -> list[Probe]:
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        required = {"Target Name", "Probe ID", "Probe Sequence", "Barcode"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Expected semicolon-delimited columns: {sorted(required)}")
        probes = []
        for line_no, row in enumerate(reader, 2):
            seq = row["Probe Sequence"].strip().upper()
            if not seq or set(seq) - set("ACGT"):
                raise ValueError(f"Invalid DNA sequence at input line {line_no}: {seq!r}")
            probes.append(Probe(row["Target Name"].strip(), row["Probe ID"].strip(), seq,
                                row["Barcode"].strip()))
    if len({p.probe_id for p in probes}) != len(probes):
        raise ValueError("Probe ID values must be unique")
    return probes


def fasta_records(path: str | Path):
    with open_text(path) as handle:
        header = None
        chunks: list[str] = []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].rstrip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks).upper()


def map_probes_to_transcripts(probes: list[Probe], fasta_path: str | Path, seed_length: int = 15):
    """Yield exact matches using a short-seed index and full-length verification."""
    seed_index: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for i, probe in enumerate(probes):
        for orientation, pattern in (("as_provided", probe.sequence),
                                     ("reverse_complement", reverse_complement(probe.sequence))):
            seed_index[pattern[:seed_length]].append((i, orientation, pattern))

    seen: set[tuple[int, str, str, str, int]] = set()
    for header, sequence in fasta_records(fasta_path):
        match = TRANSCRIPT_RE.search(header)
        if not match:
            continue
        accession = match.group(1)
        seq_match = GENOMIC_SEQ_RE.search(header)
        genomic_seqid = seq_match.group(1) if seq_match else ""
        limit = len(sequence) - seed_length + 1
        for pos in range(max(0, limit)):
            candidates = seed_index.get(sequence[pos:pos + seed_length])
            if not candidates:
                continue
            for probe_index, orientation, pattern in candidates:
                if sequence.startswith(pattern, pos):
                    key = (probe_index, accession, genomic_seqid, orientation, pos)
                    if key not in seen:
                        seen.add(key)
                        yield key + (len(pattern),)


def parse_attributes(text: str) -> dict[str, str]:
    result = {}
    for field in text.rstrip().split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = unquote(value)
    return result


def parse_refseq_gff(path: str | Path):
    """Return transcript metadata and exon loci keyed by accession."""
    metadata: dict[str, dict[str, str]] = {}
    aliases_by_gene: dict[str, set[str]] = defaultdict(set)
    exon_groups: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    locus_attrs: dict[tuple[str, str, str], tuple[str, str, str]] = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"gene", "transcript", "mRNA", "exon"}:
                continue
            seqid, source, feature, start, end, _, strand, _, attr_text = fields
            attrs = parse_attributes(attr_text)
            if feature == "gene":
                gene = attrs.get("gene", "")
                if gene:
                    aliases_by_gene[gene].update(filter(None, attrs.get("gene_synonym", "").split(",")))
                continue
            accession = attrs.get("transcript_id")
            if not accession:
                continue
            gene = attrs.get("gene", "")
            product = attrs.get("product", "")
            biotype = attrs.get("gbkey", feature)
            metadata.setdefault(accession, {
                "gene": gene, "aliases": ",".join(sorted(aliases_by_gene.get(gene, set()))),
                "product": product, "biotype": biotype, "source": source
            })
            if feature == "exon":
                key = (accession, seqid, strand)
                exon_groups[key].append((int(start) - 1, int(end)))
                locus_attrs[key] = (gene, product, biotype)

    loci: dict[str, list[TranscriptLocus]] = defaultdict(list)
    for (accession, seqid, strand), exons in exon_groups.items():
        gene, product, biotype = locus_attrs[(accession, seqid, strand)]
        ordered = sorted(exons, reverse=(strand == "-"))
        loci[accession].append(TranscriptLocus(accession, gene, product, biotype,
                                               seqid, strand, ordered))
    return metadata, loci


def read_ucsc_names(assembly_report: str | Path) -> dict[str, str]:
    names = {}
    with open(assembly_report, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) >= 10:
                ucsc = fields[9] if fields[9] != "na" else fields[6]
                for alias in (fields[0], fields[4], fields[6]):
                    if alias != "na":
                        names[alias] = ucsc
    return names


def transcript_interval_to_blocks(locus: TranscriptLocus, start: int, end: int):
    """Convert a 0-based half-open transcript interval to genomic BED blocks."""
    blocks = []
    transcript_offset = 0
    for exon_start, exon_end in locus.exons:
        exon_length = exon_end - exon_start
        overlap_start = max(start, transcript_offset)
        overlap_end = min(end, transcript_offset + exon_length)
        if overlap_start < overlap_end:
            within_start = overlap_start - transcript_offset
            within_end = overlap_end - transcript_offset
            if locus.strand == "+":
                blocks.append((exon_start + within_start, exon_start + within_end))
            else:
                blocks.append((exon_end - within_end, exon_end - within_start))
        transcript_offset += exon_length
        if transcript_offset >= end:
            break
    return sorted(blocks)


def write_tsv(path: str | Path, fieldnames: list[str], rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
