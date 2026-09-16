#!/usr/bin/env Rscript

# Summarize probe/isoform coverage in data/cosmx-lung.tsv.
#
# Usage:
#   Rscript summarize_cosmx_lung_probe_coverage.R [input.tsv] [output.tsv] [annotation.gtf.gz]
#
# The Ensembl GTF is needed only for the MANE Select metric because the input
# gene summary contains transcript accessions, but not transcript tags.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
    getwd()
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) {
    stop(
        "Usage: Rscript summarize_cosmx_lung_probe_coverage.R ",
        "[input.tsv] [output.tsv] [annotation.gtf.gz]",
        call. = FALSE
    )
}

input_tsv <- if (length(args) >= 1L) args[[1L]] else
    file.path(script_dir, "..", "data", "cosmx-lung.tsv")
output_tsv <- if (length(args) >= 2L) args[[2L]] else
    file.path(script_dir, "results_ensembl", "probe_coverage_summary.tsv")
gtf_file <- if (length(args) >= 3L) args[[3L]] else
    file.path(
        script_dir, "ensembl", "116",
        "Homo_sapiens.GRCh38.116.chr_patch_hapl_scaff.gtf.gz"
    )

for (path in c(input_tsv, gtf_file)) {
    if (!file.exists(path)) {
        stop("Required input does not exist: ", path, call. = FALSE)
    }
}

coverage <- read.delim(
    input_tsv,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_columns <- c(
    "target_gene",
    "probe_count",
    "annotated_ensembl_isoform_count",
    "hit_isoform_count",
    "common_to_all_probes_count",
    "common_to_all_probes",
    "classification"
)
missing_columns <- setdiff(required_columns, names(coverage))
if (length(missing_columns) > 0L) {
    stop(
        "Input is missing required column(s): ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
    )
}
if (anyDuplicated(coverage$target_gene)) {
    stop("Input must contain exactly one row per target gene.", call. = FALSE)
}

count_columns <- c(
    "probe_count",
    "annotated_ensembl_isoform_count",
    "hit_isoform_count",
    "common_to_all_probes_count"
)
if (anyNA(coverage[count_columns]) ||
        any(vapply(coverage[count_columns], function(x) any(x < 0), logical(1)))) {
    stop("Coverage count columns must contain non-negative, non-missing values.", call. = FALSE)
}
if (any(coverage$hit_isoform_count > coverage$annotated_ensembl_isoform_count) ||
        any(coverage$common_to_all_probes_count > coverage$hit_isoform_count)) {
    stop("Input contains inconsistent transcript coverage counts.", call. = FALSE)
}

read_mane_accessions <- function(path, chunk_size = 100000L) {
    connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    } else {
        file(path, open = "rt")
    }
    on.exit(close(connection), add = TRUE)

    accessions <- character()
    repeat {
        lines <- readLines(connection, n = chunk_size, warn = FALSE)
        if (length(lines) == 0L) break

        is_mane_transcript <-
            grepl("\ttranscript\t", lines, fixed = TRUE) &
            grepl('tag "MANE_Select"', lines, fixed = TRUE)
        lines <- lines[is_mane_transcript]
        if (length(lines) == 0L) next

        transcript_id <- sub(
            '.*transcript_id "([^"]+)".*', "\\1", lines, perl = TRUE
        )
        has_version <- grepl('transcript_version "', lines, fixed = TRUE)
        transcript_version <- sub(
            '.*transcript_version "([^"]+)".*', "\\1", lines, perl = TRUE
        )
        accessions <- c(
            accessions,
            ifelse(
                has_version,
                paste0(transcript_id, ".", transcript_version),
                transcript_id
            )
        )
    }
    unique(accessions)
}

mane_accessions <- read_mane_accessions(gtf_file)
if (length(mane_accessions) == 0L) {
    stop("No MANE Select transcripts were found in the GTF.", call. = FALSE)
}

has_annotated_isoforms <- coverage$annotated_ensembl_isoform_count > 0L
has_five_probes <- coverage$probe_count == 5L
all_isoforms_hit <-
    has_annotated_isoforms &
    coverage$hit_isoform_count == coverage$annotated_ensembl_isoform_count
all_isoforms_hit_by_five <-
    has_annotated_isoforms &
    has_five_probes &
    coverage$common_to_all_probes_count ==
        coverage$annotated_ensembl_isoform_count

common_accessions <- strsplit(
    coverage$common_to_all_probes,
    split = ",",
    fixed = TRUE
)
mane_hit_by_all_five <- vapply(
    common_accessions,
    function(x) any(x[nzchar(x)] %in% mane_accessions),
    logical(1)
) & has_five_probes

summary_table <- data.frame(
    genes_with_differential_probe_coverage = sum(
        coverage$classification ==
            "multiple_isoforms_differential_probe_coverage"
    ),
    genes_with_all_isoforms_covered = sum(all_isoforms_hit),
    genes_with_every_isoform_covered_by_all_five_probes = sum(
        all_isoforms_hit_by_five
    ),
    annotated_transcripts_covered_by_panel = sum(coverage$hit_isoform_count),
    transcripts_covered_by_all_five_probes = sum(
        coverage$common_to_all_probes_count[has_five_probes]
    ),
    transcripts_with_no_exact_probe_match = sum(
        coverage$annotated_ensembl_isoform_count -
            coverage$hit_isoform_count
    ),
    genes_with_mane_select_covered_by_all_five_probes = sum(
        mane_hit_by_all_five
    ),
    check.names = FALSE
)

output_dir <- dirname(output_tsv)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}
write.table(
    summary_table,
    file = output_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
)

message("Wrote probe coverage summary to ", normalizePath(output_tsv))
