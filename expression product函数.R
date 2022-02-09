library(Seurat)  
library(openxlsx)
library(ComplexHeatmap)
cc=readRDS("/home/wjyan/integrated.rds")  #cc--Colorectal cancer
seurat_c=cc[,cc@seurat_c$cellmeta.data$TumorSite=="C"] #seurat_c--cancer
seurat_n=cc[,cc@meta.data$TumorSite=="N"] #seurat_n--normal
df_c=as.data.frame(seurat_c@assays$integrated@data)+5 #矩阵有负值
df_n=as.data.frame(seurat_n@assays$integrated@data)+5
seurat_c$celltype=seurat_c@meta.data$CellType
seurat_n$celltype=seurat_n@meta.data$CellType
intercell_interactions <- read.xlsx("/home/wjyan/immune_lr2.xlsx")
###################
mtx=df_n/
mtx=df_c
LR <- FindLR(DB=intercell_interactions,mtx)
mtx_lr <- mtx[rownames(mtx) %in% c(LR$ligand,LR$receptor),] ##取出配体受体对所对应的表达矩阵
LR_mean <- exp.mean(mtx_lr,celltypes=seurat_n$celltype,LR=LR)
interactions <- interaction_product(LR_mean)
interactions=interactions[interactions$res_product>50,]
write.csv(interactions,file="c_interactions.csv",row.names = T)

####################单独热图-相互作用强度
res_c <- read.csv("C:/Users/YWJ/Desktop/c_interactions.csv", header=T)
res_n <- read.csv("C:/Users/YWJ/Desktop/n_interactions.csv", header=T)
dt_c <- for_heatmap(res_c)
dt_n <- for_heatmap(res_n)
Visual_heatmap(dt_c,net.diff=dt_c,color.heatmap = "Reds")
Visual_heatmap(dt_n,net.diff=dt_n,color.heatmap = "Reds")

##################### 两组间差异热图，差异circos
net.diff=dt_c-dt_n
Visual_diffInteraction(object=net.diff, weight.scale = T,measure = "weight")
Visual_heatmap_diff(object=net.diff, measure = "weight")


###########################
interaction_product <- function(LR_mean){
  #function description:
  #interaction_product,To calculate the product of ligands from sender celltype interacts with the receptors to receiver celltype
  #arguments:
  #param LR_mean -- the mean of ligand and receptor genes in cell type
  #output:
  #interaction_product,the product of the ligands from sender celltype interacts with the receptors to receiver celltype
  mean.ligand=LR_mean[["ligand"]]
  mean.receptor=LR_mean[["receptor"]]

  res=data.frame()
  for(cell1 in colnames(mean.ligand)){
    for(cell2 in colnames(mean.receptor)){
      col_res=as.data.frame(cbind(mean.ligand[,cell1],mean.receptor[,cell2]))
    #  print(head(col_res))
      col_res$ligand=rownames(mean.ligand)
      col_res$receptor=rownames(mean.receptor)
      print(dim(col_res))
      col_res$from=cell1
      col_res$to=cell2
      res=rbind(res,col_res)}}
  res_product <- res$V1*res$V2
  res <- cbind(res,res_product)
  res <- res[order(res$res_product,decreasing = T),]
  interaction_product =res
  rownames(interaction_product ) <- NULL
  return(interaction_product )
}
#########
for_heatmap <- function(res){
  #function description:
  #for_heatmap,Calculate the sum of the product of ligand receptor interactions in each cell type
  #arguments:
  #param res -- result of ligand and receptor interactions between cell types
  #output:
  #dt,the sum of the product of ligand receptor interactions in each cell type
  dt <- aggregate(res$res_product,by=list(res$from,res$to),FUN = sum)[,1:3]
  dt<- reshape2::dcast(dt,Group.1~Group.2,fun.aggregate = mean,value.var = "x")
  rownames(dt) <- dt$Group.1
  dt$Group.1 <- NULL
  dt[is.na(dt)]<-0
  dt <- as.matrix(dt)
  return(dt)}


#from CellChat--show number of interactions or interaction strength in the cell-cell communication network in one dataset;
Visual_heatmap <- function (object, comparison = c(1, 2), measure = c("weight", "count"), 
                            color.use = NULL,title.name = NULL, width = NULL, height = NULL, 
                            font.size = 8, font.size.title = 10, cluster.rows = FALSE, net.diff=NULL,
                            cluster.cols = FALSE, sources.use = NULL, targets.use = NULL, 
                            remove.isolate = FALSE, row.show = NULL, col.show = NULL,color.heatmap = "Reds") 
{
  if (measure == "count") {
    if (is.null(title.name)) {
      title.name = "Number of interactions"
    }
  }
  else if (measure == "weight") {
    if (is.null(title.name)) {
      title.name = "Interaction strength"
    }
  }
  legend.name <- title.name
  net <- net.diff
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat, 
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                                                                                                     "\\1", max(mat, na.rm = T))) + 1))
  }
  else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), round(max(mat, 
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                                                                                                  "\\1", max(mat, na.rm = T))) + 1))
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = legend.name, bottom_annotation = col_annotation, 
                left_annotation = row_annotation, top_annotation = ha2, 
                right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), column_title = title.name, 
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title), 
                row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                                fontface = "plain"), title_position = "leftcenter-rot", 
                                                                border = NA, legend_height = unit(20, "mm"), 
                                                                labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                                                                                                                  "mm")))
  return(ht1)
}

#from CellChat --Circle plot showing differential cell-cell communication network between two datasets
Visual_diffInteraction <- function (object,measure = c("count","weight"), 
          color.use = NULL, color.edge = c("#b2182b", "#2166ac"),
          title.name = NULL, top = 1, weight.scale = FALSE, vertex.weight = 20,
          vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex = 1, 
          vertex.label.color = "black", edge.weight.max = NULL, 
          edge.width.max = 8,alpha.edge = 0.6, label.edge = FALSE, 
          edge.label.color = "black",edge.label.cex = 0.8, edge.curved = 0.2, 
          shape = "circle", layout = in_circle(), margin = 0.2, 
          arrow.width = 1,arrow.size = 0.2) 
  {
  options(warn = -1)
  measure <- match.arg(measure)
  if (measure %in% c("count")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- object
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top)] <- 0
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
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
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
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}

############This heatmap can be used to show differential number of interactions or interaction strength in the cell-cell communication network between two datasets;
Visual_heatmap_diff <- function (object, measure = c("count", 
                              "weight"),color.use = NULL, color.heatmap = c("#2166ac","#b2182b"), 
          title.name = NULL, width = NULL, height = NULL, 
          font.size = 8, font.size.title = 10, cluster.rows = FALSE, 
          cluster.cols = FALSE, sources.use = NULL, targets.use = NULL, 
          remove.isolate = FALSE, row.show = NULL, col.show = NULL) 
{  
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  net <- object
  if (measure == "count") {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure == "weight") {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  legend.name = "Relative values"
  net[is.na(net)] <- 0
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat, 
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                                                                                                     "\\1", max(mat, na.rm = T))) + 1))
  }
  else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), round(max(mat, 
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                                                                                                  "\\1", max(mat, na.rm = T))) + 1))
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = legend.name, bottom_annotation = col_annotation, 
                left_annotation = row_annotation, top_annotation = ha2, 
                right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), column_title = title.name, 
                column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                row_title = "Sources (Sender)", row_title_gp = gpar(fontsize = font.size.title), 
                row_title_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                                fontface = "plain"), title_position = "leftcenter-rot", 
                                                                border = NA, legend_height = unit(20, "mm"), 
                                                                labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                                                                                                                  "mm")))
  return(ht1)
}
  
  
  
  