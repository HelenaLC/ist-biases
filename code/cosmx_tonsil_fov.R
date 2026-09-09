# dependencies
suppressPackageStartupMessages({
    library(zoo)
    library(ggplot2)
    library(ggrastr)
    library(SpatialExperiment)
})

# retrieval
options(timeout=1e3)
url <- "https://zenodo.org/records/16782547/files/C1.rds"
dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
rds <- file.path(dir, "data", basename(url))
download.file(url, rds)

# loading
pa <- dirname(tf)
unzip(tf, exdir=pa)
# importing
library(alabaster.sce)
nm <- file.path(pa, "C1")
sce <- readObject(nm)
# coersion
xy <- sprintf("Center%s_global_mm", c("X", "Y"))
xy <- as.matrix(colData(sce)[xy])
(spe <- toSpatialExperiment(sce, spatialCoords=xy))

# wrangling
df <- data.frame(colData(sce))
df <- df |>
    arrange(CenterX_global_mm) |>
    mutate(mx=zoo::rollmedian(val, 5e3+1, fill=NA, align="center")) |>
    arrange(CenterY_global_mm) |>
    mutate(my=zoo::rollmedian(val, 5e3+1, fill=NA, align="center"))

# plotting
gg <- ggplot(df, aes(CenterX_global_mm, val)) + 
    geom_point_rast(shape=16, stroke=0, size=0.02, alpha=0.02, col="blue") +
    geom_line(aes(y=mx), col="gold", linewidth=0.2) +
    ggplot(df, aes(CenterY_global_mm, val)) + 
    geom_point_rast(shape=16, stroke=0, size=0.02, alpha=0.02, col="blue") +
    geom_line(aes(y=my), col="gold", linewidth=0.2) +
    plot_layout(ncol=1) &
    labs(y="counts per area") &
    scale_y_continuous(limits=c(0, 25), n.breaks=3) &
    theme_bw(4) & theme(aspect.ratio=1/3, panel.grid.minor=element_blank())

# saving
pdf <- file.path(dir, "figs", "fov_cosmx.pdf")
ggsave(pdf, gg, units="cm", width=4, height=3)