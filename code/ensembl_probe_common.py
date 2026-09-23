#!/usr/bin/env python3
"""Shared utilities for reproducible exact probe mapping to Ensembl 116."""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
import shutil
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from urllib.request import urlopen


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
REFERENCE_DIR = SCRIPT_DIR / "ensembl" / "116"
DEFAULT_PROBES = DATA_DIR / "cosmx-lung-probes.csv"

# SHA-256 digests pin the exact release-116 files used for the analysis.
ENSEMBL_FILES = {
    "cdna": (
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cdna/"
        "Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "683eb19310c40bf1396e4718f45afa2ce86755717c0990f47a171f535d248ea1",
        183898799,
    ),
    "ncrna": (
        "Homo_sapiens.GRCh38.ncrna.fa.gz",
        "https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/ncrna/"
        "Homo_sapiens.GRCh38.ncrna.fa.gz",
        "7f03bb303e939517b322f7887a74c1ee15bbb82a6affbb54f0be84e82f89cff1",
        41013640,
    ),
    "gtf": (
        "Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz",
        "https://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/"
        "Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz",
        "ef38b04cde03949d3b6a965575f1cad7861098ef029e7e8ed35c84debcaf13b8",
        145778890,
    ),
}

UCSC_CHROM_ALIAS = (
    "hg38.chromAlias.txt.gz",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz",
    "43aac315a93939c54d8b168ea8118fb7760b5e7ff8653f59dd712cfc7467be56",
    13393,
)

DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
ENSEMBL_TRANSCRIPT_RE = re.compile(r"^(ENST[0-9]+(?:\.[0-9]+)?)\s")
GTF_ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)";')


@dataclass(frozen=True)
class Probe:
    target: str
    probe_id: str
    sequence: str


@dataclass
class TranscriptLocus:
    accession: str
    gene: str
    transcript_name: str
    biotype: str
    seqid: str
    strand: str
    exons: list[tuple[int, int]]


def open_text(path: str | Path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, encoding="utf-8")


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_download(reference, directory: Path = REFERENCE_DIR) -> Path:
    """Download a pinned reference if absent and verify its size and SHA-256."""
    filename, url, expected_sha256, expected_size = reference
    path = directory / filename
    if path.exists() and path.stat().st_size == expected_size and sha256(path) == expected_sha256:
        return path

    directory.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".part")
    print(f"Downloading {url}", file=sys.stderr)
    try:
        with urlopen(url) as response, open(partial, "wb") as output:
            shutil.copyfileobj(response, output)
        observed_size = partial.stat().st_size
        observed_sha256 = sha256(partial)
        if observed_size != expected_size or observed_sha256 != expected_sha256:
            raise RuntimeError(
                f"Verification failure for {filename}: expected {expected_size} bytes and "
                f"SHA-256 {expected_sha256}, observed {observed_size} bytes and "
                f"SHA-256 {observed_sha256}"
            )
        partial.replace(path)
    finally:
        if partial.exists():
            partial.unlink()
    return path


def ensure_ensembl_file(kind: str) -> Path:
    return ensure_download(ENSEMBL_FILES[kind])


def ensure_ucsc_chrom_alias() -> Path:
    return ensure_download(UCSC_CHROM_ALIAS)


def read_probes(path: str | Path) -> list[Probe]:
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        required = {"Target Name", "Probe ID", "Probe Sequence", "Barcode"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Expected semicolon-delimited columns: {sorted(required)}")
        probes = []
        for line_no, row in enumerate(reader, 2):
            sequence = row["Probe Sequence"].strip().upper()
            if len(sequence) < 15 or set(sequence) - set("ACGT"):
                raise ValueError(
                    f"Expected an A/C/G/T probe of at least 15 bases at line {line_no}: "
                    f"{sequence!r}"
                )
            probes.append(Probe(
                row["Target Name"].strip(), row["Probe ID"].strip(), sequence
            ))
    if len({probe.probe_id for probe in probes}) != len(probes):
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


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


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


def map_probes_to_ensembl(
    probes: list[Probe], fasta_paths: list[str | Path], seed_length: int = 15
):
    seed_index: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for index, probe in enumerate(probes):
        for orientation, pattern in (
            ("as_provided", probe.sequence),
            ("reverse_complement", reverse_complement(probe.sequence)),
        ):
            seed_index[pattern[:seed_length]].append((index, orientation, pattern))

    seen: set[tuple[int, str, str, str, int]] = set()
    for fasta_path in fasta_paths:
        for header, sequence in fasta_records(fasta_path):
            parsed = parse_ensembl_header(header)
            if not parsed:
                continue
            accession, genomic_seqid = parsed
            for position in range(max(0, len(sequence) - seed_length + 1)):
                candidates = seed_index.get(sequence[position:position + seed_length])
                if not candidates:
                    continue
                for probe_index, orientation, pattern in candidates:
                    if sequence.startswith(pattern, position):
                        key = (probe_index, accession, genomic_seqid, orientation, position)
                        if key not in seen:
                            seen.add(key)
                            yield key + (len(pattern),)


def parse_gtf_attributes(text: str) -> dict[str, list[str]]:
    attributes: dict[str, list[str]] = defaultdict(list)
    for key, value in GTF_ATTRIBUTE_RE.findall(text):
        attributes[key].append(value)
    return attributes


def first(attributes: dict[str, list[str]], key: str, default: str = "") -> str:
    values = attributes.get(key)
    return values[0] if values else default


def versioned_id(attributes: dict[str, list[str]], id_key: str, version_key: str) -> str:
    stable_id = first(attributes, id_key)
    version = first(attributes, version_key)
    return f"{stable_id}.{version}" if stable_id and version else stable_id


def parse_ensembl_gtf(path: str | Path):
    metadata: dict[str, dict[str, str]] = {}
    exon_groups: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    locus_attributes: dict[tuple[str, str, str], tuple[str, str, str]] = {}

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon"}:
                continue
            seqid, source, feature, start, end, _, strand, _, attribute_text = fields
            attributes = parse_gtf_attributes(attribute_text)
            accession = versioned_id(attributes, "transcript_id", "transcript_version")
            if not accession:
                continue
            gene = first(attributes, "gene_name")
            gene_id = versioned_id(attributes, "gene_id", "gene_version")
            biotype = first(
                attributes, "transcript_biotype", first(attributes, "gene_biotype")
            )
            transcript_name = first(attributes, "transcript_name")
            metadata.setdefault(accession, {
                "gene": gene,
                "gene_id": gene_id,
                "transcript_name": transcript_name,
                "biotype": biotype,
                "source": first(attributes, "transcript_source", source),
                "tags": ",".join(attributes.get("tag", [])),
                "transcript_support_level": first(attributes, "transcript_support_level"),
            })
            if feature == "exon":
                key = (accession, seqid, strand)
                exon_groups[key].append((int(start) - 1, int(end)))
                locus_attributes[key] = (gene, transcript_name, biotype)

    loci: dict[str, list[TranscriptLocus]] = defaultdict(list)
    for (accession, seqid, strand), exons in exon_groups.items():
        gene, transcript_name, biotype = locus_attributes[(accession, seqid, strand)]
        loci[accession].append(TranscriptLocus(
            accession=accession,
            gene=gene,
            transcript_name=transcript_name,
            biotype=biotype,
            seqid=seqid,
            strand=strand,
            exons=sorted(exons, reverse=(strand == "-")),
        ))
    return metadata, loci


def read_ucsc_chrom_aliases(path: str | Path) -> dict[str, str]:
    aliases = {}
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) >= 2:
                aliases[fields[0]] = fields[1]
    return aliases


def transcript_interval_to_blocks(locus: TranscriptLocus, start: int, end: int):
    """Convert a zero-based, half-open transcript interval to genomic BED blocks."""
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
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
