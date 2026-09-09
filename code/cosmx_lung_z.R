# dependencies
suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(KernSmooth)
})

# loading
dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
sub <- file.path(dir, "data", "cosmx-lung")
sub <- list.dirs(sub, recursive=FALSE)
csv <- list.files(sub, "tx", full.names=TRUE)
csv <- grepv("Lung(5|9)_Rep2", csv)
tx <- lapply(csv, read_csv_arrow)
names(tx) <- basename(dirname(csv))

# lines ----

ps <- lapply(tx, \(df) {
    # count transcripts by FOV
    fd <- df |>
        group_by(fov, z) |>
        summarise(n=n(), .groups="drop")
    
    # keep z-stacks containing all FOVs
    nf <- n_distinct(fd$fov)
    fd <- fd |>
        group_by(z) |>
        filter(n() == nf) |>
        arrange(desc(z)) |>
        mutate(z=factor(z, seq(0, 8)))
    
    # aesthetics
    nz <- nlevels(fd$z)
    pal <- rev(hcl.colors(nz, "Spectral"))
    names(pal) <- levels(fd$z)
    dy <- with(fd, c(
        10^floor(log10(min(n+1))),
        10^ceiling(log10(max(n+1)))))
    
    # plotting
    ggplot(fd, aes(
        x=factor(fov), y=n+1, 
        group=factor(z), color=factor(z))) + # avoid zeros
        geom_line(linewidth=1/6, alpha=2/3, show.legend=FALSE) +
        geom_point(size=1/4, alpha=2/3, show.legend=TRUE) +
        scale_x_discrete("field of view (FOV)",
            breaks=c(seq(5, 50, 5), range(fd$fov))) +
        scale_y_log10(
            expression(log[10]*"(counts+1)"),
            labels=scales::label_log(base=10),
            breaks=10^seq(1, 6), limits=dy) +
        guides(col=guide_legend(override.aes=list(size=2/3, alpha=1))) +
        scale_color_manual("z", values=pal, 
            drop=FALSE, limits=levels(fd$z)) +
        ggtitle(names(tx_list)[i]) +
        coord_cartesian(expand=FALSE, clip="off") +
        theme_minimal(4) + theme(
            legend.key.size=unit(0, "pt"),
            panel.grid.minor=element_blank())
})

# saving
gg <- wrap_plots(ps, ncol=1, guides="collect") 
pdf <- file.path(dir, "figs", "cosmx_z.pdf")
ggsave(pdf, gg, units="cm", width=4, height=3.5)

# space ----

.kde <- \(tx, zi, nx, ny) {
    df <- filter(tx, z == zi)
    xy <- c("x_global_px", "y_global_px")
    kde <- bkde2D(
        x=as.matrix(df[, xy]),
        bandwidth=c(100, 100),
        gridsize=c(nx, ny))
    kde_df <- expand.grid(x=kde$x1, y=kde$x2) |>
        mutate(density=as.vector(kde$fhat), z=zi) |>
        # scale by total number of transcripts
        mutate(intensity=density*nrow(df)) 
}

dfs <- lapply(tx, \(df) {
    # estimate grid size
    bw <- 100
    xr <- range(df$x_global_px)
    yr <- range(df$y_global_px)
    nx <- ceiling(diff(xr)/bw)
    ny <- ceiling(diff(yr)/bw)
    
    # calculate density
    zs <- seq(min(df$z), max(df$z))
    ds <- lapply(zs, \(z) .kde(df, z, nx, ny))
    df <- do.call(rbind, ds)
})

for (. in seq_along(dfs)) {
    # plotting
    zs <- seq(0, ifelse(. == 1, 8, 7))
    df <- filter(dfs[[.]], z %in% zs)
    gg <- ggplot(df, aes(x, y)) +
        facet_wrap(~z, nrow=1) +
        geom_raster(aes(fill=intensity)) +
        scale_fill_gradientn(colors=pals::brewer.rdpu(100)) +
        coord_equal(expand=FALSE) +
        theme_void(4) + theme(
            plot.margin=margin(), 
            legend.key.width=unit(0.4, "lines"),
            legend.key.height=unit(0.8, "lines"))
    
    # saving
    pdf <- "cosmx_z_%s.pdf"
    pdf <- file.path(dir, "figs", pdf)
    pdf <- sprintf(pdf, names(tx)[.])
    ggsave(pdf, gg, units="cm", width=diff(range(df$x))/1e3, height=5)
}
