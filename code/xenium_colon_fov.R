# dependencies
suppressPackageStartupMessages({
    library(zoo)
    library(dplyr)
    library(ggrastr)
    library(ggplot2)
    library(patchwork)
    library(OSTA.data)
    library(SpatialExperiment)
    library(SpatialExperimentIO)
})

# loading
id <- "Xenium_HumanColon_Oliveira"
pa <- OSTA.data_load(id, mol=FALSE)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
se <- readXeniumSXE(td, addTx=FALSE)

# wrangling
df <- data.frame(
    colData(se), 
    spatialCoords(se)) |>
    mutate(across(matches("centroid"), \(.) ./1e3))
.x <- df[[x <- "x_centroid"]]
.y <- df[[y <- "y_centroid"]]
xs <- cut(.x, seq(min(.x), max(.x), l=100), include.lowest=TRUE, right=TRUE)
ys <- cut(.y, seq(min(.y), max(.y), l=100), include.lowest=TRUE, right=TRUE)
df <- mutate(df, xs, ys, val=total_counts/cell_area)

# smoothing
df <- df |>
    arrange(.data[[x]]) |> 
    mutate(mx=rollmedian(val, 5e3+1, fill=NA, align="center")) |>
    arrange(.data[[y]]) |> 
    mutate(my=rollmedian(val, 5e3+1, fill=NA, align="center"))

# fields of view
sf  <- 0.2125/1e3 # px to mm
wum <- 3520*sf    # width
hum <- 2960*sf    # height
oum <-  125*sf    # overlap
xl <- seq(0, max(df[[y]]), hum-oum)
yl <- seq(0, max(df[[x]]), wum-oum)

# plotting
px <- ggplot(df, aes(.data[[x]], val)) + 
    labs(x="global X coordinate (mm)") +
    geom_point_rast(shape=16, stroke=0, size=0.02, alpha=0.02, col="blue") +
    geom_line(aes(y=mx), col="gold", linewidth=0.2) +
    geom_vline(lty=2, lwd=0.2, xintercept=xl, col="red") 
py <- ggplot(df, aes(.data[[y]], val)) + 
    labs(x="global Y coordinate (mm)") +
    geom_point_rast(shape=16, stroke=0, size=0.02, alpha=0.02, col="blue") +
    geom_line(aes(y=my), col="gold", linewidth=0.2) +
    geom_vline(lty=2, lwd=0.2, xintercept=yl, col="red") 
gg <- wrap_plots(px, py, ncol=1) &
    labs(y="counts per area") &
    scale_y_continuous(limits=c(0, 5), breaks=seq(0, 4, 2)) &
    theme_bw(4) & theme(aspect.ratio=1/3, panel.grid.minor=element_blank())

# saving
dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
pdf <- file.path(dir, "figs", "fov_xenium.pdf")
ggsave(pdf, gg, units="cm", width=4, height=3)
