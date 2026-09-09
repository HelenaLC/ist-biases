# data available from ... were pre-downloaded b/c retrieval requires logging in.
# https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/

# dependencies
suppressPackageStartupMessages({
    library(pals)
    library(dplyr)
    library(readxl)
    library(scater)
    library(ggplot2)
    library(tidyverse)
    library(SingleCellExperiment)
})

# loading
dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))

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

# panel
md <- file.path(dir, "data", "cosmx-lung", "metadata.xlsx")
suppressWarnings(md <- read_xlsx(md, sheet=3))
names(md) <- c("rna", "id", "seq", "bc")
md$gc <- 100*str_count(md$seq, "G|C")/nchar(md$seq)

# quality
ids <- split(seq(ncol(sce)), sce$sid)
qc <- perFeatureQCMetrics(sce, subsets=ids)

# wrangling
df <- data.frame(rna=rownames(qc), qc)
df <- filter(df, !grepl("^NegPrb", rna))
df <- md |>
    group_by(rna) |>
    summarise_at("gc", mean) |>
    right_join(df, by="rna")
xs <- with(df, cut(gc, 
    seq(min(gc), max(gc), l=4),
    right=TRUE, include.lowest=TRUE))

# aggregation
pat <- "subsets_(.*)_detected"
fd <- df |>
    mutate(xs) |>
    pivot_longer(all_of(matches(pat))) |>
    mutate(name=gsub(pat, "\\1", name)) |>
    group_by(xs, name) |>
    summarise_at("value", median) |>
    arrange(value)

# plotting
pal <- c(stepped2()[c(2, 4)], stepped2()[5:7], stepped2()[11], stepped2()[c(14, 16)])
gg <- ggplot(fd, aes(xs, value, fill=name, col=name)) + 
    labs(x="mean %GC", y="median detection") +
    geom_path(aes(group=name), show.legend=FALSE) + 
    geom_point(shape=21, size=0.8, stroke=0.2, col="white") + 
    #scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    scale_y_continuous(limits=c(0, 15), breaks=seq(0, 15, 5)) +
    guides(fill=guide_legend(NULL, override.aes=list(stroke=NA))) +
    annotate("text", levels(fd$xs), 0, size=1.2, label=sprintf("(%s)", table(xs))) +
    scale_color_manual(values=pal) +
    scale_fill_manual(values=pal) +
    coord_equal(1/2.25) +
    theme_bw(4) + theme(
        legend.key.size=unit(0, "pt"),
        panel.grid.minor=element_blank(),
        axis.text.x=element_text(angle=30, hjust=1))

# saving
pdf <- file.path(dir, "figs", "cosmx_gc.pdf")
ggsave(pdf, gg, units="cm", width=3, height=3)
