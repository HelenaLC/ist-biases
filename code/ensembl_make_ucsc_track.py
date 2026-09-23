#!/usr/bin/env python3
"""Map probes to Ensembl transcripts and create a UCSC hg38 BED12 track."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from ensembl_probe_common import (
    DEFAULT_PROBES, SCRIPT_DIR, ensure_ensembl_file, ensure_ucsc_chrom_alias,
    map_probes_to_ensembl, parse_ensembl_gtf, read_probes, read_ucsc_chrom_aliases,
    transcript_interval_to_blocks, write_tsv,
)


ALIGNMENT_FIELDS = [
    "target_gene", "probe_id", "probe_sequence", "transcript_accession",
    "genomic_seqid", "match_orientation", "transcript_start_0based", "transcript_end_0based",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probes", default=DEFAULT_PROBES)
    parser.add_argument("--transcripts", nargs="+",
                        help="Ensembl cDNA and ncRNA FASTA files; release 116 is downloaded by default")
    parser.add_argument("--gtf", help="Ensembl GTF; release 116 is downloaded by default")
    parser.add_argument("--bed", default=SCRIPT_DIR / "ensembl/probes_hg38_ensembl116.bed")
    parser.add_argument("--alignments", default=SCRIPT_DIR / "ensembl/probe_transcript_alignments.tsv")
    args = parser.parse_args()

    transcripts = args.transcripts or [ensure_ensembl_file("cdna"), ensure_ensembl_file("ncrna")]
    gtf = args.gtf or ensure_ensembl_file("gtf")
    probes = read_probes(args.probes)
    _, loci = parse_ensembl_gtf(gtf)
    chromosome_aliases = read_ucsc_chrom_aliases(ensure_ucsc_chrom_alias())
    matches = sorted(map_probes_to_ensembl(probes, transcripts),
                     key=lambda x: (probes[x[0]].target, probes[x[0]].probe_id, x[1], x[2], x[4]))

    Path(args.bed).parent.mkdir(parents=True, exist_ok=True)
    alignment_rows = []
    bed_records = {}
    missing_locus = set()
    for probe_index, accession, genomic_seqid, orientation, tx_start, length in matches:
        probe = probes[probe_index]
        tx_end = tx_start + length
        alignment_rows.append({
            "target_gene": probe.target, "probe_id": probe.probe_id,
            "probe_sequence": probe.sequence, "transcript_accession": accession,
            "genomic_seqid": genomic_seqid, "match_orientation": orientation,
            "transcript_start_0based": tx_start, "transcript_end_0based": tx_end,
        })
        matched_locus = False
        for locus in loci.get(accession, []):
            if genomic_seqid and locus.seqid != genomic_seqid:
                continue
            blocks = transcript_interval_to_blocks(locus, tx_start, tx_end)
            if not blocks or sum(b - a for a, b in blocks) != length:
                continue
            matched_locus = True
            chrom = chromosome_aliases.get(locus.seqid, locus.seqid)
            chrom_start, chrom_end = blocks[0][0], blocks[-1][1]
            key = (probe.probe_id, chrom, locus.strand, tuple(blocks))
            bed_records[key] = (
                chrom, chrom_start, chrom_end, f"{probe.target}|{probe.probe_id}", 0,
                locus.strand, chrom_start, chrom_end, "166,54,3", len(blocks),
                ",".join(str(b - a) for a, b in blocks) + ",",
                ",".join(str(a - chrom_start) for a, _ in blocks) + ",",
            )
        if not matched_locus:
            missing_locus.add((accession, genomic_seqid))

    write_tsv(args.alignments, ALIGNMENT_FIELDS, alignment_rows)
    with open(args.bed, "w", encoding="utf-8", newline="") as handle:
        handle.write('track name="Probe mappings Ensembl 116" '
                     'description="Exact Ensembl 116 probe mappings" visibility=pack itemRgb="On"\n')
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(sorted(bed_records.values(), key=lambda r: (r[0], r[1], r[2], r[3])))

    matched_probes = len({probes[i].probe_id for i, *_ in matches})
    print(f"Probes: {len(probes)}; exact transcript matches: {len(matches)}; "
          f"matched probes: {matched_probes}; BED loci: {len(bed_records)}")
    if missing_locus:
        print(f"Warning: {len(missing_locus)} matched transcript/locus pairs could not be projected to GTF exons")


if __name__ == "__main__":
    main()
