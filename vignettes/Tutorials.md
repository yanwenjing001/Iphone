# Tutorials

## Step 0: Load required packages

`library(Iphone)`

`library(openxlsx)`

## Step 1: Data input & processing 

Here we load PBMC 3k scRNA-seq seurat object

`pbmc <- readRDS("output/pbmc3k_final.rds")
 
 mtx <- as.data.frame(pbmc@assays$RNA@data)

pbmc$celltype <- pbmc@active.ident

intercell_interactions <- read.xlsx("C:/Users/YWJ/Desktop/Omnipath/immune_lr2.xlsx")`

## Step 2: Using three different methods to calculate interactions between cell types

### expression threshold
`LR <- FindLR(DB=intercell_interactions,mtx)

 mtx_lr <- mtx[rownames(mtx) %in% c(LR$ligand,LR$receptor),]
 
 LR_mean <- exp.mean(mtx_lr,LR,celltype=pbmc$celltype)
 
 threshold_interactions <- interaction(LR_mean,threshold=0.5)
 
 #write.csv(threshold_interactions,file="C:/Users/YWJ/Desktop/threshold_interactions.csv",row.names = T)`
