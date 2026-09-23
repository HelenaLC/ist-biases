# Probe-to-Ensembl mapping

This analysis uses these pinned human references:

- Ensembl release: 116 (June 2026)
- assembly: GRCh38.p14 / `GCA_000001405.29`
- coding and pseudogene transcripts: `Homo_sapiens.GRCh38.cdna.all.fa.gz`
- non-coding transcripts: `Homo_sapiens.GRCh38.ncrna.fa.gz`
- annotation: `Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz`

The scripts download missing references directly from the Ensembl release-116
archive into `ensembl/116/`. Each download is checked against its exact byte
count and a pinned SHA-256 digest (and was independently checked against the
corresponding Ensembl `CHECKSUMS` entry). GTF is the GFF2-derived annotation
format consumed by the coordinate and metadata parser.

```text
683eb19310c40bf1396e4718f45afa2ce86755717c0990f47a171f535d248ea1  Homo_sapiens.GRCh38.cdna.all.fa.gz
7f03bb303e939517b322f7887a74c1ee15bbb82a6affbb54f0be84e82f89cff1  Homo_sapiens.GRCh38.ncrna.fa.gz
ef38b04cde03949d3b6a965575f1cad7861098ef029e7e8ed35c84debcaf13b8  Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz
43aac315a93939c54d8b168ea8118fb7760b5e7ff8653f59dd712cfc7467be56  hg38.chromAlias.txt.gz
```

The small pinned UCSC chromosome-alias table converts Ensembl primary,
patch, and haplotype sequence names to valid hg38 custom-track names.

Python 3.8 or later is required. There are no third-party Python dependencies.
The probe input is `../data/cosmx-lung-probes.csv`, a semicolon-delimited CosMx
probe export. Its SHA-256 digest is:

```text
4983f4c247e7c2bf7923c6cffebc6de0d8ab44bb81d61870319e6ec5ca06d107
```

## Reproduce the outputs

```sh
cd code
python3 ensembl_make_ucsc_track.py
python3 ensembl_match_probes.py
```

The aliases `C9orf16` -> `BBLN` and `DDX58` -> `RIGI` are recorded in 
`ensembl_match_probes.py`. They affect intended-gene labeling only, 
not sequence matches.

## Outputs

- `ensembl/probes_hg38_ensembl116.bed`: UCSC hg38 BED12 custom track.
- `ensembl/probe_transcript_alignments.tsv`: raw exact transcript hits.
- `ensembl/probe_isoform_matches.tsv`: annotated long match table with
  versioned Ensembl transcript and gene identifiers, biotypes, tags, transcript
  support levels, source assembly regions, and orientations.
- `ensembl/probe_isoform_summary.tsv`: one row per input probe.
- `ensembl/gene_isoform_summary.tsv`: one row per target gene.

A successful reproduction creates a gene summary identical to
`../data/cosmx-lung.tsv` (961 lines including its header), with SHA-256 digest:

```text
62be89cc4e4197cd4bb589d77a48c4937a1544237589876357c6e33ef30a8d12
```

Verify it with:

```sh
cmp ensembl/gene_isoform_summary.tsv ../data/cosmx-lung.tsv
shasum -a 256 ensembl/gene_isoform_summary.tsv
```

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
