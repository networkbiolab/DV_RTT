setwd("[home_dir]")
library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(SingleCellExperiment)
set.seed(42)
#Prep Sample
XLraw<- readRDS("[home_dir]/all objects/SRRXL.rds")
DefaultAssay(XLraw)<- "RNA"
XL <- XLraw
XL <- JoinLayers(XL)
#sigh
XL[["percent.mito"]] <- PercentageFeatureSet(XL, pattern = "^MT-")
#XL[["RNA"]] <- as(object = XL[["RNA"]], Class = "Assay")
XL <- subset(XL, subset = nFeature_RNA > 500 & nFeature_RNA < 3147 &percent.mito < 10)
# <- split(XL[["RNA"]],f= XL$orig.ident)
XL<- SplitObject(XL, split.by = "orig.ident")
for (i in seq_along(XL)) {
  XL[[i]] <- NormalizeData(XL[[i]]) %>% FindVariableFeatures()
}
features <- SelectIntegrationFeatures(XL)
for (i in seq_along(along.with = XL)) {
  XL[[i]] <- ScaleData(XL[[i]], features = features) %>% RunPCA(features = features)
}
#originally 300 dims but lets check xd)
anchors <- FindIntegrationAnchors(XL, reference = c(1, 2), reduction = "rpca", dims = 1:30)
integrated <- IntegrateData(anchors, dims = 1:30)
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated)
cds <- as.cell_data_set(integrated)
cds <- cluster_cells(cds,resolution = 1.05e-4)
integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
table(integrated.sub$monocle3_clusters)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
seur <- as.Seurat(cds,assay = NULL)

#cds <- readRDS("MonoJumbosize.rds")

