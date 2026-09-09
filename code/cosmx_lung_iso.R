# data available from ... were pre-downloaded b/c retrieval requires logging in.
# https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/

# dependencies
suppressPackageStartupMessages({
    library(scater)
    library(Matrix)
    library(tidyverse)
    library(SingleCellExperiment)
})

# loading
dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
tsv <- file.path(dir, "data", "cosmx-lung.tsv")
iso <- read_tsv(tsv, show_col_types=FALSE)

# assay
sub <- list.dirs(
    file.path(dir, "data", "cosmx-lung"), 
    full.names=TRUE, recursive=FALSE)
names(sub) <- basename(sub)
sce <- lapply(names(sub), \(.) {
    csv <- list.files(sub[.], "csv$", full.names=TRUE)
    cd <- read.csv(grepv("metadata", csv)); cd$sid <- .
    mx <- read.csv(grepv("exprMat", csv), check.names=FALSE)
    mx <- mx[match(cd$cell_ID, mx$cell_ID, nomatch=0), ]
    mx <- as(as.matrix(mx[, -c(1, 2)]), "dgCMatrix")
    se <- SingleCellExperiment(assays=list(counts=t(mx)), colData=cd)
}) |> do.call(what=cbind)

# quality
ids <- split(seq(ncol(sce)), sce$sid)
qc <- perFeatureQCMetrics(sce, subsets=ids)

# wrangling
pat <- "subsets_(.*)_detected"
df <- data.frame(qc) |>
    select(-c(mean, detected)) |>
    select(all_of(matches(pat))) |>
    mutate(target=rownames(sce)) |>
    pivot_longer(all_of(matches(pat))) |>
    mutate(name=gsub(pat, "\\1", name)) |>
    inner_join(iso, by=c("target"="target_gene")) |>
    mutate(iso=factor(distinct_probe_isoform_sets))

# summarization
fd <- df |> 
    select(target, name, value, iso) |>
    group_by(name, iso) |>
    summarise_at("value", median)
ns <- df |>
    distinct(target, iso) |>
    with(table(iso))

# plotting
pal <- c(stepped2()[c(2, 4)], stepped2()[5:7], stepped2()[11], stepped2()[c(14, 16)])
gg <- ggplot(fd, aes(iso, value, col=name, fill=name)) +
    labs(x="distinct isoforms", y="median detection") +
    geom_path(aes(group=name), show.legend=FALSE) + 
    geom_point(shape=21, size=0.8, stroke=0.2, col="white") + 
    scale_y_continuous(limits=c(0, 15), breaks=seq(0, 15, 5)) +
    guides(fill=guide_legend(NULL, override.aes=list(stroke=NA))) +
    annotate("text", levels(fd$iso), 0, size=1.2, label=sprintf("(%s)", ns)) +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    coord_equal(1/2.25) +
    theme_bw(4) + theme(
        legend.key.size=unit(0, "pt"),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1))

# saving
pdf <- file.path(dir, "figs", "cosmx_iso.pdf")
ggsave(pdf, gg, units="cm", width=4, height=3)
