library(openxlsx)
library(Seurat)
library(CCInx)  
pbmc <- readRDS("output/pbmc3k_final.rds")
mtx <- as.data.frame(pbmc@assays$RNA@data)  
pbmc$celltype <- pbmc@active.ident
intercell_interactions <- read.xlsx("C:/Users/YWJ/Desktop/Omnipath/immune_lr2.xlsx")

LR <- FindLR(DB=intercell_interactions,mtx)
LR$key <- paste0(LR$ligand,"_",LR$receptor)
pbmc_l <- pbmc[rownames(pbmc@assays$RNA@data)%in% LR$ligand]
pbmc_r <- pbmc[rownames(pbmc@assays$RNA@data)%in% LR$receptor]
pbmc.l_markers <- FindAllMarkers(pbmc_l, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.01)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc_l, features = top10$gene) + NoLegend()

pbmc.r_markers <- FindAllMarkers(pbmc_r, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.01)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc_r, features = top10$gene) + NoLegend()

pbmc.markers <- rbind(pbmc.l_markers,pbmc.r_markers)
pbmc.markers <- pbmc.markers[!duplicated(pbmc.markers$gene),]
dt <- lapply(unique(pbmc.markers$cluster),function(x){pbmc.markers[pbmc.markers$cluster==x,]})
names(dt) <-  unique(pbmc.markers$cluster)

inx <- BuildCCI(GeneStatList=dt,
                  GeneMagnitude="avg_log2FC",
                  GeneStatistic="p_val_adj",
                  Species="hsapiens")

PlotCCInx(INX=inx,
          cellTypeA="CD14+ Mono",cellTypeB="CD8 T",
          proteinTypeA="Ligand",proteinTypeB="Receptor")

######from CCInx --Build cell-cell interaction predictions between cell types
##改用自己的配体受体库
BuildCCI <- function (GeneStatList, GeneMagnitude = "MeanNormGeneExpr", 
          GeneStatistic, Species = "hsapiens") 
{
  if (length(names(GeneStatList)) != length(GeneStatList)) {
    stop("GeneStatList must be a named list where names are cell types.")
  }
  if (any(duplicated(names(GeneStatList)))) {
    stop("GeneStatList names must be unique.")
  }
  if (any(grepl("_", names(GeneStatList)))) {
    stop("GeneStatList names cannot contain '_' due to internal naming conventions.")
  }
  if (any(grepl("~", names(GeneStatList)))) {
    stop("GeneStatList names cannot contain '~' due to internal naming conventions.")
  }
  message("Scaling node weights per cell type...")
  if (missing(GeneStatistic)) {
    temp_scaled <- pbapply::pbsapply(X = GeneStatList, FUN = CalcExprScaled, 
                                     expr = GeneMagnitude, simplify = F)
  }
  else {
    temp_scaled <- pbapply::pbsapply(X = GeneStatList, FUN = CalcDiffExprScaled, 
                                     DEmagn = GeneMagnitude, DEstat = GeneStatistic, simplify = F)
  }
  inx <- list()
  message("Building node metadata...")
  temp_cellNames <- names(GeneStatList)
  inx$nodes <- do.call(rbind, GeneStatList)
  temp_gene <- inx$nodes$gene
  temp_cellType <- unlist(mapply(function(N, X) rep(N, X), 
                                 N = temp_cellNames, X = sapply(GeneStatList, nrow), SIMPLIFY = F), 
                          use.names = F)
  switch(Species, hsapiens = load(system.file("LigRecDB_RData/BaderCCIeditedbyBI_human.RData", 
                                              package = "CCInx")), mmusculus = load(system.file("LigRecDB_RData/BaderCCIeditedbyBI_mouse.RData", 
                                                                                                package = "CCInx")), MillerKaplan = load(system.file("LigRecDB_RData/MillerKaplan_mouse.RData", 
                                                                                                                                                     package = "CCInx")), FANTOM5 = load(system.file("LigRecDB_RData/FANTOM5_human.RData", 
                                                                                                                                                                                                     package = "CCInx")), stop("Species must be one of 'hsapiens' or 'mmusculus'."))
  if (sum(rownames(geneInfo) %in% temp_gene) < 20) {
    warning(paste("Less than 20 genes from GeneStatList were detected in the CCInx database.", 
                  "  Please ensure that you've set the Species argument correctly.", 
                  "  Rownames of each entry in GeneStatList must be official gene symbols.", 
                  sep = "\n"))
  }
  temp_proteinType <- geneInfo[temp_gene, "protein_type"]
  inx$nodes <- cbind(data.frame(node = paste(temp_gene, temp_cellType, 
                                             sep = "_"), gene = temp_gene, cellType = temp_cellType, 
                                proteinType = temp_proteinType, nodeWeight = unlist(temp_scaled), 
                                stringsAsFactors = F), inx$nodes)
  rownames(inx$nodes) <- inx$nodes$node
  inx$nodes <- inx$nodes[!is.na(inx$nodes$proteinType), ]
  tempCN <- c()
  for (a in temp_cellNames) {
    for (b in temp_cellNames) {
      temp <- paste(sort(c(a, b)), collapse = "~")
      if (!temp %in% tempCN) {
        tempCN <- append(tempCN, temp)
      }
    }
  }
  rm(a, b)
  tempComp <- strsplit(tempCN, "~")
  names(tempComp) <- tempCN
  message("Building edge list...")
  inx$edges <- pbapply::pbsapply(tempComp, function(Z) {
    a <- Z[1]
    b <- Z[2]
    if (sum(inx$nodes$cellType == a) < 1 | sum(inx$nodes$cellType == 
                                               b) < 1) {
      return(NULL)
    }
    keysAB <- LR$key[LR$ligand %in% inx$nodes$gene[inx$nodes$cellType == 
                                                          a] & LR$receptor %in% inx$nodes$gene[inx$nodes$cellType == 
                                                                                                 b]]
    if (length(keysAB) < 1) {
      edgesAB <- data.frame(row.names = c("key", 
                                          "nodeA", "nodeB"))
    }
    else {
      edgesAB <- data.frame(sapply(strsplit(keysAB, "_"), 
                                   function(X) paste(paste(X[1], a, sep = "_"), 
                                                     paste(X[2], b, sep = "_"), sep = "~")), 
                            t(sapply(strsplit(keysAB, "_"), function(X) c(paste(X[1], 
                                                                                a, sep = "_"), paste(X[2], b, sep = "_")))), 
                            stringsAsFactors = F)
      colnames(edgesAB) <- c("key", "nodeA", 
                             "nodeB")
      rownames(edgesAB) <- edgesAB$key
    }
    keysBA <- LR$key[LR$ligand %in% inx$nodes$gene[inx$nodes$cellType == 
                                                          b] & LR$receptor %in% inx$nodes$gene[inx$nodes$cellType == 
                                                                                                 a]]
    if (length(keysBA) < 1) {
      edgesBA <- data.frame(row.names = c("key", 
                                          "nodeA", "nodeB"))
    }
    else {
      edgesBA <- data.frame(sapply(strsplit(keysBA, "_"), 
                                   function(X) paste(paste(X[2], a, sep = "_"), 
                                                     paste(X[1], b, sep = "_"), sep = "~")), 
                            t(sapply(strsplit(keysBA, "_"), function(X) c(paste(X[2], 
                                                                                a, sep = "_"), paste(X[1], b, sep = "_")))), 
                            stringsAsFactors = F)
      colnames(edgesBA) <- c("key", "nodeA", 
                             "nodeB")
      rownames(edgesBA) <- edgesBA$key
    }
    temp <- rbind(edgesAB, edgesBA)
    if (nrow(temp) < 1) {
      return(NULL)
    }
    else {
      return(temp)
    }
  }, simplify = F)
  inx$edges <- do.call(rbind, inx$edges)
  rownames(inx$edges) <- inx$edges$key
  inx$edges <- inx$edges[, 2:3]
  inx$edges$edgeWeight <- rowMeans(cbind(inx$nodes[inx$edges$nodeA, 
                                                   "nodeWeight"], inx$nodes[inx$edges$nodeB, "nodeWeight"]))
  attr(inx, "GeneMagnitude") <- GeneMagnitude
  if (!missing(GeneStatistic)) {
    attr(inx, "GeneStatistic") <- GeneStatistic
  }
  return(inx)
}

#######CCInx
CalcDiffExprScaled <- function (gdb, DEmagn, DEstat) 
{
  if (any(is.na(gdb[[DEmagn]]))) {
    stop(paste("This function doesn't tolerate missing", 
               DEmagn, "values."))
  }
  if (any(is.na(gdb[[DEstat]]))) {
    stop(paste("This function doesn't tolerate missing", 
               DEstat, "values."))
  }
  return(sapply(gdb[[DEmagn]], function(X) {
    if (X == Inf) {
      1.1
    } else if (X == -Inf) {
      -1.1
    } else {
      X/switch(as.character(X >= 0), `TRUE` = max(gdb[[DEmagn]][!is.infinite(gdb[[DEmagn]])]), 
               `FALSE` = min(gdb[[DEmagn]][!is.infinite(gdb[[DEmagn]])]) * 
                 -1)
    }
  }))
}


#####
######## CCInx
PlotCCInx <- function (INX, cellTypeA, cellTypeB, proteinTypeA, proteinTypeB, 
          GeneMagnitudeThreshold, GeneStatisticThreshold, TopEdges, 
          GeneNames, YSpacing = "relative") 
{
  if (missing(INX) | missing(cellTypeA) | missing(cellTypeB) | 
      missing(proteinTypeA) | missing(proteinTypeB)) {
    stop("The following arguments are required: INX, cellTypeA, cellTypeB, proteinTypeA, proteinTypeB.")
  }
  INX <- FilterInx_step1(INX, cellTypeA = cellTypeA, cellTypeB = cellTypeB, 
                         proteinTypeA = proteinTypeA, proteinTypeB = proteinTypeB)
  if (!missing(GeneMagnitudeThreshold)) {
    INX <- FilterInx_GeneMagnitude(INX, GeneMagnitudeThreshold)
  }
  else if (!missing(GeneStatisticThreshold)) {
    INX <- FilterInx_GeneStatistic(INX, GeneStatisticThreshold)
  }
  else if (!missing(TopEdges)) {
    INX <- FilterInx_topN(INX, TopEdges)
  }
  else if (!missing(GeneNames)) {
    INX <- FilterInx_genenames(INX, GeneNames)
  }
  DoPlotInx(INX, YSpacing)
}
