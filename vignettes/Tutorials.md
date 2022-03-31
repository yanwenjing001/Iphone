# Tutorials

## Step 0: Load required packages

`library(Iphone)`   
`library(openxlsx)`

## Step 1: Data input & processing 

Here we load PBMC 3k scRNA-seq seurat object

`pbmc <- readRDS("output/pbmc3k_final.rds")`    
`mtx <- as.matrix(pbmc@assays$RNA@data)`          
`pbmc$celltype <- pbmc@active.ident`                    
`intercell_interactions <- read.xlsx("C:/Users/YWJ/Desktop/Omnipath/immune_lr2.xlsx")`

## Step 2: Using three different methods to calculate interactions between cell types

### 1.expression threshold
`LR <- findLR(DB=intercell_interactions,mtx)`            
`mtx_lr <- mtx[rownames(mtx) %in% c(LR$ligand,LR$receptor),]`            
`LR_mean <- LRmean(mtx_lr,LR,celltype=pbmc$celltype)`                 
`threshold.interactions <- interaction_threshold(LR_mean,threshold=0.5)`            
`#write.csv(threshold.interactions,file="C:/Users/YWJ/Desktop/threshold.interactions.csv",row.names = T)`            

### Visualization     
#### plotPAC potential apparent competition
`df_ligand <- ligand_pac(interactions)`        
`df_receptor <- receptor_pac(interactions)`            
`plotPAC(df_ligand,scaling=0.8,plot.scale=1,fill.col="black", arrow.col="black", circles=TRUE, radius=1)`            
`plotPAC(df_receptor,scaling=0.8,plot.scale=1,fill.col="black", arrow.col="black", circles=TRUE, radius=1)`           
#### Chord plot
`dt_ threshold <- table(intercell_interactions$from,intercell_interactions$to)`            
`groupSize <-rowSums(dt_threshold)`         
`netVisual_circle(dt_threshold, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")`             
`par(mfrow = c(3,3), xpd=TRUE)`    
`for (i in 1:nrow(dt)) {`       
  `dt2 <- matrix(0, nrow = nrow(dt_threshold), ncol = ncol(dt_threshold), dimnames = dimnames(dt_threshold))`           
  `dt2[i,] <- dt_threshold[i,]`          
  `netVisual_circle(dt2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = dt_threshold(dt_threshold), title.name = rownames(dt_threshold)[i])`          
`}`         

### 2.expression product
`LR <- findLR(DB=intercell_interactions,mtx)`            
`mtx_lr <- mtx[rownames(mtx) %in% c(LR$ligand,LR$receptor),]`            
`LR_mean <- LRmean(mtx_lr,LR,celltype=pbmc$celltype)`                 
`product.interactions <- interaction_product(LR_mean)`     

### Visualization      
`dt_product <- for_heatmap(product.interactions)`        
`Visual_heatmap(net,color.heatmap = "Reds")`        

### 3.Differential combinations
`LR <- FindLR(DB=intercell_interactions,mtx)`         
`expr.markers <- find_markers(object=pbmc,ident="celltype")`    
`combination.interactions <- interaction_combination(expr.markers,LR)`       

### Visualization     






































