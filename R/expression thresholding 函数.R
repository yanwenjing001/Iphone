library(openxlsx)
library(bipartite)
library(igraph)
library(CellChat)    
pbmc <- readRDS("output/pbmc3k_final.rds")
mtx <- as.data.frame(pbmc@assays$RNA@data)
pbmc$celltype <- pbmc@active.ident
intercell_interactions <- read.xlsx("C:/Users/YWJ/Desktop/Omnipath/immune_lr2.xlsx")

LR <- FindLR(DB=intercell_interactions,mtx) ## 1.从表达矩阵中寻找配体和受体基因
mtx_lr <- mtx[rownames(mtx) %in% c(LR$ligand,LR$receptor),] ##取出配体受体对所对应的表达矩阵
LR_mean <- exp.mean(mtx_lr,LR,celltype=pbmc$celltype) ## 2.计算不同细胞类型中配体和受体基因的平均值
interactions <- interaction(LR_mean,threshold=1) ## 3.寻找细胞间相互作用
write.csv(SP1vsD2,file="C:/Users/YWJ/Desktop/SP1vsD2.csv",row.names = T)

##可视化
#plotPAC
df_ligand <- ligand_pac(interactions)
df_receptor <- receptor_pac(interactions)
plotPAC(df_ligand,scaling=0.8,plot.scale=1, 
        fill.col="black", arrow.col="black", circles=TRUE, radius=1)
plotPAC(df_receptor,scaling=0.8,plot.scale=1, 
        fill.col="black", arrow.col="black", circles=TRUE, radius=1)
#circle
res <- interactions
dt <- table(res$from,res$to)
groupSize <-rowSums(dt)
netVisual_circle(dt, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
mat <- dt
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


##############################
FindLR <- function(DB,mtx){
  #function description:
  #FindLR,find ligand and receptor genes from expression matrix 
  #arguments:
  #param DB -- database of ligand-receptor pairs from priori knowledge 
  #param mtx -- expression matrix,the rows are the genes, the columns are the cells
  #output:
  #LR,ligand-receptor pairs from expression matrix 
                   LR <- DB[,c("ligand","receptor")]
                   LR=LR[!duplicated(paste0(LR$ligand,LR$receptor)),]
                   LR_gene=intersect(rownames(mtx), c(LR$ligand,LR$receptor))##取出数据中属于配体或受体的基因
                   LR=LR[LR$ligand %in% LR_gene & LR$receptor %in% LR_gene,] ##取出LR_gene对应的配体受体
                   return(LR)
                   }

exp.mean <- function(mtx_lr,celltypes,LR=LR){
  #function description:
  #exp.mean,Calculate the mean of ligand and receptor genes in cell type
  #arguments:
  #param mtx_lr -- expression matrix about ligand-receptor genes
  #param celltypes -- expression matrix clustering results
  #param LR -- ligand-receptor pairs from expression matrix 
  #output:
  #LR_mean,the mean of ligand and receptor genes in cell type
                    res=data.frame()
                    for(celltype in unique(celltypes)){    ##求每种细胞类型对应的基因表达均值
                     print(celltype)
                     cells=names(celltypes)[celltypes==celltype]
                     mean_lr=apply(mtx_lr[,colnames(mtx_lr) %in% cells],1,mean)
                     #mean_receptor=apply(mtx_receptor[,colnames(mtx_receptor) %in% cells],1,mean)
                     res=rbind(res,mean_lr)
                   }   
                   mean_lr=as.data.frame(t(res))
                   colnames(mean_lr)=unique(celltypes)
                   rownames(mean_lr)=rownames(mtx_lr)
                   
                   mean_ligand=list()
                   mean_receptor=list()
                   for(col in colnames(mean_lr)){
                     #col=colnames(mean_lr)[1]
                     print(col)
                     col_expr=mean_lr[,col]
                     names(col_expr)=rownames(mean_lr)
                     col=gsub(" ","_",col)#替换细胞类型中的空格
                     print(col)
                     mean_ligand[[col]]=col_expr[LR$ligand]
                     mean_receptor[[col]]=col_expr[LR$receptor]
                   }
                   mean.ligand=do.call(cbind,mean_ligand)
                   mean.receptor=do.call(cbind,mean_receptor)
                   LR_mean<- list(ligand=mean.ligand,receptor=mean.receptor)
                   return(LR_mean)
}

interaction <- function(LR_mean,threshold=1){
  #function description:
  #interaction,Find the ligands from sender celltype interacts with the receptors to receiver celltype
  #arguments:
  #param LR_mean -- the mean of ligand and receptor genes in cell type
  #param threshold -- the mean expression value of screening ligands and receptors is greater than 1
  #output:
  #interactions,the ligands from sender celltype interacts with the receptors to receiver celltype
                 mean.ligand=LR_mean[["ligand"]]
                 mean.receptor=LR_mean[["receptor"]]
                 
                 binary_ligand=mean.ligand
                 binary_ligand[binary_ligand >=threshold] <- 1
                 binary_ligand[binary_ligand < threshold] <- 0
                 binary_receptor=mean.receptor
                 binary_receptor[binary_receptor >= threshold] <- 1
                 binary_receptor[binary_receptor <  threshold] <- 0
                 
                 #过滤行，保证L/R在至少一种细胞中大于阈值
                 judge=rowSums(binary_ligand)>0 & rowSums(binary_receptor)>0
                 binary_ligand=binary_ligand[judge,]
                 binary_receptor=binary_receptor[judge,]
                 
                 res=data.frame()
                 for(cell1 in colnames(binary_ligand)){
                 for(cell2 in colnames(binary_receptor)){
                 col_res=as.data.frame(cbind(binary_ligand[,cell1],binary_receptor[,cell2]))
                 print(head(col_res))
                 col_res$ligand=rownames(binary_ligand)
                 col_res$receptor=rownames(binary_receptor)
                 print(dim(col_res))
                 col_res$from=cell1
                 col_res$to=cell2
                 res=rbind(res,col_res)}}
                 interactions=res[rowSums(res[,1:2])==2,]
                 interactions=interactions[c("ligand","receptor","from","to")]
                 rownames(interactions) <- NULL
                 return(interactions)
                 }

ligand_pac <- function(interactions){
  #function description:
  #ligand_pac,Count the number of interactions between all receptors with ligand genes for different cell types
  #arguments:
  #param interactions -- the ligands from sender celltype interacts with the receptors to receiver celltype
  #output:
  #df_ligand,a data.frame of the number of interactions between all receptors with ligand genes for different cell types
                df_ligand=data.frame()
                ligand=(unique(interactions$ligand))
                
                for(celltype in unique(interactions$to)){  
                vec=rep(0,length(ligand))
                res_celltype=interactions[interactions$to==celltype,]
                tab=table(res_celltype$ligand)
                names(vec)=ligand
                vec[names(tab)]=tab
                df_ligand=rbind(df_ligand,vec)
                }
                rownames(df_ligand)=unique(interactions$to)
                colnames(df_ligand)=ligand
                return(df_ligand)}
                
receptor_pac <- function(interactions){
  #function description:
  #receptor_pac,Count the number of interactions between all ligands with receptor genes for different cell types
  #arguments:
  #param interactions -- the ligands from sender celltype interacts with the receptors to receiver celltype
  #output:
  #df_receptor,a data.frame of the number of interactions between all ligands with receptor genes for different cell types
              df_receptor=data.frame()
              receptor=(unique(interactions$receptor))

              for(celltype in unique(interactions$from)){
              vec=rep(0,length(receptor))
              res_celltype=interactions[interactions$from==celltype,]
              tab=table(res_celltype$receptor)
              names(vec)=receptor
              vec[names(tab)]=tab
              df_receptor=rbind(df_receptor,vec)
              }
              rownames(df_receptor)=unique(interactions$from)
              colnames(df_receptor)=receptor
              return(df_receptor)}

### from CellChat 
scPalette <- function (n) 
{
  colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", 
                  "#984EA3", "#F29403", "#F781BF", "#BC9DCC", 
                  "#A65628", "#54B0E4", "#222F75", "#1B9E77", 
                  "#B2DF8A", "#E3BE00", "#FB9A99", "#E7298A", 
                  "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", 
                  "#B3DE69", "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}




netVisual_circle <- function (net, color.use = NULL, title.name = NULL, sources.use = NULL, 
          targets.use = NULL, remove.isolate = FALSE, top = 1, weight.scale = FALSE, 
          vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, 
          vertex.label.cex = 1, vertex.label.color = "black", 
          edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
          label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
          edge.curved = 0.2, shape = "circle", layout = in_circle(), 
          margin = 0.2, vertex.size = NULL, arrow.width = 1, arrow.size = 0.2) 
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
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
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
