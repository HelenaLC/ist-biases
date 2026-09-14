# Probe-to-Ensembl mapping

This analysis uses the latest Ensembl human release available on 2026-07-22:

- Ensembl release: 116 (June 2026)
- assembly: GRCh38.p14 / `GCA_000001405.29`
- coding and pseudogene transcripts: `Homo_sapiens.GRCh38.cdna.all.fa.gz`
- non-coding transcripts: `Homo_sapiens.GRCh38.ncrna.fa.gz`
- annotation: `Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz`

The three downloads in `ensembl/116/` were checked successfully with the
official Ensembl `CHECKSUMS` manifests stored alongside them.

## Reproduce the outputs

```sh
python3 make_ensembl_ucsc_track.py \
  --transcripts \
    ensembl/116/Homo_sapiens.GRCh38.cdna.all.fa.gz \
    ensembl/116/Homo_sapiens.GRCh38.ncrna.fa.gz \
  --gtf ensembl/116/Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz \
  --assembly-report refseq/GCF_000001405.40_GRCh38.p14_assembly_report.txt

python3 match_probes_to_ensembl.py \
  --gtf ensembl/116/Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz
```

No non-standard Python packages are required. `gene_aliases.tsv` reconciles the
input symbols `C9orf16` and `DDX58` with the current annotation symbols `BBLN`
and `RIGI`. This affects intended-gene labeling only, not sequence matches.

## Outputs

- `results_ensembl/probes_hg38_ensembl116.bed`: UCSC hg38 BED12 custom track.
- `results_ensembl/probe_transcript_alignments.tsv`: raw exact transcript hits.
- `results_ensembl/probe_isoform_matches.tsv`: annotated long match table with
  versioned Ensembl transcript and gene identifiers, biotypes, tags, transcript
  support levels, source assembly regions, and orientations.
- `results_ensembl/probe_isoform_summary.tsv`: one row per input probe.
- `results_ensembl/gene_isoform_summary.tsv`: one row per target gene.

The gene-level classifications have the same definitions as the RefSeq run:

- `single_isoform`: collectively, the probes match one intended-gene transcript.
- `multiple_isoforms_same_probe_coverage`: every probe matches the same set of
  multiple transcripts.
- `multiple_isoforms_differential_probe_coverage`: different probes match
  different transcript sets.
- `incomplete_probe_mapping`: at least one probe lacks an exact intended-gene hit.
- `no_exact_intended_match`: no probe has an exact intended-gene hit.

Both orientations are searched; all intended-gene matches are reverse-complement
hits, consistent with antisense hybridization probes. Matching is exact and uses
the complete Ensembl gene set rather than restricting to GENCODE Basic, canonical,
protein-coding, or well-supported transcripts. Those attributes are retained in
the long table for optional filtering.
