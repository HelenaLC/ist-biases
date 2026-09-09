# dependencies
suppressPackageStartupMessages({
    library(zoo)
    library(ggplot2)
    library(ggrastr)
    library(SpatialExperiment)
})

# loading
options(timeout=1e3)
url <- "https://zenodo.org/records/16782547/files/C1.rds?download=1"
dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
rds <- file.path(dir, "data", "cosmx-tonsil.rds")
download.file(url, rds)
sce <- readRDS(rds)

# wrangling
df <- data.frame(colData(sce))
.x <- df[[x <- "CenterX_local_mm"]]
.y <- df[[y <- "CenterY_local_mm"]]
xs <- cut(.x, seq(min(.x), max(.x), l=100), include.lowest=TRUE, right=TRUE)
ys <- cut(.y, seq(min(.y), max(.y), l=100), include.lowest=TRUE, right=TRUE)
df <- mutate(df, xs, ys, val=nCount_RNA/Area.um2)

# smoothing
df <- df |>
    arrange(.data[[x]]) |> 
    mutate(mx=rollmedian(val, 5e3+1, fill=NA, align="center")) |>
    arrange(.data[[y]]) |> 
    mutate(my=rollmedian(val, 5e3+1, fill=NA, align="center"))

# plotting
px <- ggplot(df, aes(CenterX_global_mm, val)) + 
    geom_point_rast(shape=16, stroke=0, size=0.02, alpha=0.02, col="blue") +
    geom_line(aes(y=mx), col="gold", linewidth=0.2) +
    labs(x="global X coordinate (mm)") 
py <- ggplot(df, aes(CenterY_global_mm, val)) + 
    geom_point_rast(shape=16, stroke=0, size=0.02, alpha=0.02, col="blue") +
    geom_line(aes(y=my), col="gold", linewidth=0.2) +
    labs(x="global Y coordinate (mm)") 
gg <- wrap_plots(px, py, ncol=1) &
    labs(y="counts per area") &
    scale_y_continuous(limits=c(0, 25), n.breaks=3) &
    theme_bw(4) & theme(aspect.ratio=1/3, panel.grid.minor=element_blank())

# saving
pdf <- file.path(dir, "figs", "cosmx_fov.pdf")
ggsave(pdf, gg, units="cm", width=4, height=3)
