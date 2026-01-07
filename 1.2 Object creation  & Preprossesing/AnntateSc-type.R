library (Seurat)#Load Seurat
library(dplyr) #Load dlypr
library(readr)#Load en 
library(tidyverse)
library(ggplot2)
#Analisis y pre-procesamiento
#![load object]
SRRXL <- JoinLayers(SRRXL) #if merge object Seurat v5>
SRRXL <- UpdateSeuratObject(SRRXL) #if <v4
#Standard Pipeline Seurat v5 
gc() #clear unsused memory to avoid saturation
#Sctype Annotation
pbmc <- SRRXL #ease of integration with pipeline for some reason xd
DefaultAssay(pbmc)<- "RNA"
#sctype
# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openSRRXLsx"), library, character.only = T)
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.SRRXLsx";
tissue <- "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 
# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(pbmc[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))
# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(pbmc[["RNA"]]$scale.data) else as.matrix(pbmc[["RNA"]]@scale.data)
# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(pbmc@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])
#add to metada
pbmc@meta.data$Ctype = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  pbmc@meta.data$Ctype[pbmc@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
print(sctype_scores[,1:3])
SRRXL <- pbmc
#Graphs
DimPlot(SRRXL,reduction = "umap")
DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Ctype')
#Heatmap excerpt (Heatmap works better with seuratv4 objects as its currently not widely updated to accept layer instead of slot)
# Read file with genes 
gnames <- readLines("Radlgi.txt")
gnames <- gnames[gnames != ""] # Remove empty lines and whitespace
gnames <- trimws(gnames)
#Subset by needed group
#Radplx <- subset(SRRXL, Ctype == "Radial glial cells") #celltype
#Cells8w <- subset(SRRXL, Age == "56") #Age
DoHeatmap(Radplx,features = names,group.by = "orig.ident")
#Tables for proportions
WT8w <- subset(Cells8w, Condition == "WT") 
RTT8w <- subset(Cells8w, Condition == "RTT") 
WT10w <- subset(Cells10w, Condition == "WT") 
RTT10w <- subset(Cells10w, Condition == "RTT") 
WT14w <- subset(Cells14w, Condition == "WT") 
RTT14w <- subset(Cells14w, Condition == "RTT") 
table(WT8w$Ctype)
table(RTT8w$Ctype)
table(WT10w$Ctype)
table(RTT10w$Ctype)
table(WT14w$Ctype)
table(RTT14w$Ctype)
