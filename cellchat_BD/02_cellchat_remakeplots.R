library(CellChat)
library(Seurat)
library(patchwork)
library(ComplexHeatmap)
library(DT)
library(knitr)

netVisual_circle_edit <- function (net, color.use = NULL, title.name = NULL, sources.use = NULL,
                                   targets.use = NULL, idents.use = NULL, remove.isolate = FALSE,
                                   top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL,
                                   vertex.size.max = NULL, vertex.label.cex = 1, vertex.label.color = "black",
                                   edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6,
                                   label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8,
                                   edge.curved = 0.2, shape = "circle", layout = in_circle(),
                                   margin = 0.2, vertex.size = NULL, arrow.width = 1, arrow.size = 0.2,
                                   text.x = 0, text.y = 1.5, labels=NULL, label.cex=1, label.dist.weight=3)
{
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    }
    else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1 - top)
  net[net < thresh] <- 0
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use))) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) |
                         (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]],
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  g <- graph_from_adjacency_matrix(net, mode = "directed",
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max +
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g),
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g),
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max *
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,
                                                                             1]], alpha.edge)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[,
                                                                1])] <- loop.angle[edge.start[which(edge.start[,
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)),
                               direction = -1, start = 0)
  label.dist <- (vertex.weight/max(vertex.weight) + 2)*label.dist.weight
  plot(g, edge.curved = edge.curved, vertex.shape = shape,
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist,
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica",
       edge.label.family = "Helvetica", vertex.label = labels, vertex.label.cex = label.cex)
  if (!is.null(title.name)) {
    text(text.x, text.y, title.name, cex = label.cex + 0.1)
  }
  gg <- recordPlot()
  return(gg)
}



cellchat_path <- "/data/NIDDK_IDSS/projects/NIDDK-15/cellchat_BD/trimester1/"
savedir <- '/data/NIDDK_IDSS/projects/NIDDK-15/cellchat_BD/RemakePlots/trimester1/'
dir.create(savedir, recursive=T)


cellchat.without <- readRDS(paste0(cellchat_path, "Analysis1_without_cellchat.rds"))

test <- levels(cellchat.without@idents)
i <- test %>% grep(pattern = 'FCRL')

mat <- cellchat.without@net$weight
mat2[i, ] <- mat[i, ]
netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[4])

test[test %>% grep(pattern = 'FCRL')] <- expression(atop("FCRL"^"+" * " B cell"))

png(file.path(savedir, 'WeightSplitCellType_FCRL1_Without.png'), res=300, width=7, height=7, units='in')
netVisual_circle_edit(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = expression(atop("B"^"without")), labels=test, margin = 0.5, text.y=1.5, label.cex=1.5, label.dist.weight=2.25)
dev.off()

test[test %>% grep(pattern = 'FCRL')] <- expression(atop("B"^"Without"))
png(file.path(savedir, 'WeightSplitCellType_FCRL1_Without_BWithout.png'), res=300, width=7, height=7, units='in')
netVisual_circle_edit(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = expression(atop("B"^"without")), labels=test, margin = 0.5, text.y=1.5, label.cex=1.5, label.dist.weight=2.25)
dev.off()
netVisual_circle_edit(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = expression(atop("B"^"without")), labels=test, margin = 0.5, text.y=1.75, label.cex=1.5, label.dist.weight=2.25)




cellchat_path <- "/data/NIDDK_IDSS/projects/NIDDK-15/cellchat_BD/trimester2/"
savedir <- '/data/NIDDK_IDSS/projects/NIDDK-15/cellchat_BD/RemakePlots/trimester2/'
dir.create(savedir, recursive=T)


cellchat.without <- readRDS(paste0(cellchat_path, "trimester2_without_cellchat.rds"))

test <- levels(cellchat.without@idents)
i <- test %>% grep(pattern = 'FCRL')

mat <- cellchat.without@net$weight
mat2[i, ] <- mat[i, ]

test[test %>% grep(pattern = 'FCRL')] <- expression(atop("FCRL"^"+" * " B cell"))

png(file.path(savedir, 'WeightSplitCellType_FCRL1_Without.png'), res=300, width=7, height=7, units='in')
netVisual_circle_edit(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = expression(atop("B"^"without")), labels=test, margin = 0.5, text.y=1.5, label.cex=1.5, label.dist.weight=2.25)
dev.off()

test[test %>% grep(pattern = 'FCRL')] <- expression(atop("B"^"Without"))
png(file.path(savedir, 'WeightSplitCellType_FCRL1_Without_BWithout.png'), res=300, width=7, height=7, units='in')
netVisual_circle_edit(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = expression(atop("B"^"without")), labels=test, margin = 0.5, text.y=1.5, label.cex=1.5, label.dist.weight=2.25)
dev.off()




cellchat_path <- "/data/NIDDK_IDSS/projects/NIDDK-15/cellchat_BD/trimester3/"
savedir <- '/data/NIDDK_IDSS/projects/NIDDK-15/cellchat_BD/RemakePlots/trimester3/'
dir.create(savedir, recursive=T)


cellchat.without <- readRDS(paste0(cellchat_path, "trimester3_without_cellchat.rds"))

test <- levels(cellchat.without@idents)
i <- test %>% grep(pattern = 'FCRL')

mat <- cellchat.without@net$weight
mat2[i, ] <- mat[i, ]

test[test %>% grep(pattern = 'FCRL')] <- expression(atop("FCRL"^"+" * " B cell"))

png(file.path(savedir, 'WeightSplitCellType_FCRL1_Without.png'), res=300, width=7, height=7, units='in')
netVisual_circle_edit(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = expression(atop("B"^"without")), labels=test, margin = 0.5, text.y=1.5, label.cex=1.5, label.dist.weight=2.25)
dev.off()

test[test %>% grep(pattern = 'FCRL')] <- expression(atop("B"^"Without"))
png(file.path(savedir, 'WeightSplitCellType_FCRL1_Without_BWithout.png'), res=300, width=7, height=7, units='in')
netVisual_circle_edit(mat2, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = expression(atop("B"^"without")), labels=test, margin = 0.5, text.y=1.5, label.cex=1.5, label.dist.weight=2.25)
dev.off()


