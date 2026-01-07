setwd("/mnt/610a9597-0f86-4e53-9413-6c6512763c35/sofia/PP")
library (Seurat)#Cargar Seurat
library(dplyr) #Carga dlypr para cálculos
library(cowplot)#Cargar cowplot to gráficos
library(readr)#Cargar en 
WT00.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/WT/SRR13551700/filtered_feature_bc_matrix") #Incorporar la matriz de origen obtenida en Cellranger de una de los especímenes de la condición deseada a un objeto R.
WT00 <- CreateSeuratObject(counts = WT00.data, min.cells = 3, min.features  = 120, project = "WT_700", assay = "RNA") # Crea el objeto Seurat con el archivo SRRWT data. Manteniendo los genes expresados en > = 3 células (~0,1% de los datos) y todas las células con al menos 120 genes detectados por célula.
WT00@meta.data$Age <- "14w"
WT01.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/WT/SRR13551701/filtered_feature_bc_matrix") 
WT01 <- CreateSeuratObject(counts = WT01.data, min.cells = 3, min.features  = 120, project = "WT_701", assay = "RNA")
WT01@meta.data$Age <- "14w"
WT02.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/WT/SRR13551702/filtered_feature_bc_matrix") 
WT02 <- CreateSeuratObject(counts = WT02.data, min.cells = 3, min.features  = 120, project = "WT_702", assay = "RNA")
WT02@meta.data$Age <- "14w"
WT03.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/WT/SRR13551703/filtered_feature_bc_matrix") 
WT03 <- CreateSeuratObject(counts = WT03.data, min.cells = 3, min.features  = 120, project = "WT_703", assay = "RNA")
WT03@meta.data$Age <- "14w"
WT92.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/WT/SRR13551692/filtered_feature_bc_matrix") 
WT92 <- CreateSeuratObject(counts = WT92.data, min.cells = 3, min.features  = 120, project = "WT_692", assay = "RNA")
WT92@meta.data$Age <- "8w"
WT93.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/WT/SRR13551693/filtered_feature_bc_matrix") 
WT93 <- CreateSeuratObject(counts = WT93.data, min.cells = 3, min.features  = 120, project = "WT_693", assay = "RNA")
WT93@meta.data$Age <- "8w"
WT96.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/WT/SRR13551696/filtered_feature_bc_matrix") 
WT96 <- CreateSeuratObject(counts = WT96.data, min.cells = 3, min.features  = 120, project = "WT_696", assay = "RNA")
WT96@meta.data$Age <- "10w"
WT97.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/WT/SRR13551697/filtered_feature_bc_matrix") 
WT97 <- CreateSeuratObject(counts = WT97.data, min.cells = 3, min.features  = 120, project = "WT_697", assay = "RNA")
WT97@meta.data$Age <- "10w"
SRRWT <-merge(WT00,y= c(WT01,WT02,WT03,WT92,WT93,WT96,WT97),project = "WT")
SRRWT@meta.data$Condition <-"WT"
RTT04.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/RTT/SRR13551704/filtered_feature_bc_matrix") #Incorporar la matriz de origen obtenida en Cellranger de una de los especímenes de la condición deseada a un objeto R.
RTT04 <- CreateSeuratObject(counts = RTT04.data, min.cells = 3, min.features  = 120, project = "RTT_704", assay = "RNA") # Crea el objeto Seurat con el archivo SRRRTT data. Manteniendo los genes expresados en > = 3 células (~0,1% de los datos) y todas las células con al menos 120 genes detectados por célula.
RTT04@meta.data$Age <- "14w"
RTT05.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/RTT/SRR13551705/filtered_feature_bc_matrix") 
RTT05 <- CreateSeuratObject(counts = RTT05.data, min.cells = 3, min.features  = 120, project = "RTT_705", assay = "RNA")
RTT05@meta.data$Age <- "14w"
RTT06.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/RTT/SRR13551706/filtered_feature_bc_matrix") 
RTT06 <- CreateSeuratObject(counts = RTT06.data, min.cells = 3, min.features  = 120, project = "RTT_706", assay = "RNA")
RTT06@meta.data$Age <- "14w"
RTT07.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/RTT/SRR13551707/filtered_feature_bc_matrix") 
RTT07 <- CreateSeuratObject(counts = RTT07.data, min.cells = 3, min.features  = 120, project = "RTT_707", assay = "RNA")
RTT07@meta.data$Age <- "14w"
RTT94.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/RTT/SRR13551694/filtered_feature_bc_matrix") 
RTT94 <- CreateSeuratObject(counts = RTT94.data, min.cells = 3, min.features  = 120, project = "RTT_694", assay = "RNA")
RTT94@meta.data$Age <- "8w"
RTT95.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/RTT/SRR13551695/filtered_feature_bc_matrix") 
RTT95 <- CreateSeuratObject(counts = RTT95.data, min.cells = 3, min.features  = 120, project = "RTT_695", assay = "RNA")
RTT95@meta.data$Age <- "8w"
RTT98.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/RTT/SRR13551698/filtered_feature_bc_matrix") 
RTT98 <- CreateSeuratObject(counts = RTT98.data, min.cells = 3, min.features  = 120, project = "RTT_698", assay = "RNA")
RTT98@meta.data$Age <- "10w"
RTT99.data <-Read10X (data.dir = "[home_dir]l/[Matrices]/RTT/SRR13551699/filtered_feature_bc_matrix") 
RTT99 <- CreateSeuratObject(counts = RTT99.data, min.cells = 3, min.features  = 120, project = "RTT_699", assay = "RNA")
RTT99@meta.data$Age <- "10w"
SRRRTT <-merge(RTT04,y= c(RTT05,RTT06,RTT07,RTT94,RTT95,RTT98,RTT99),project = "RTT")
SRRRTT@meta.data$Condition <-"RTT"
#merge WT and RTT seurat objects
SRRXL<- merge(SRRWTump,SRRRTTump)
#Save Raw un-processed objects
saveRDS(SRRWT, file = "SRRWT.rds")
saveRDS(SRRRTT, file = "SRRRTT.rds")
saveRDS(SRRXL, file = "SRRXL.rds")
#Pre-procesing
#Change to object of interest [WT][RTT][XL]
SRRXL <- JoinLayers(SRRXL) #v5
SRRXL[["percent.mt"]] <- PercentageFeatureSet(SRRXL, pattern = "^MT-")
SRRXL<- NormalizeData(SRRXL, normalization.method = "LogNormalize", scale.factor = 10000)
SRRXL<- FindVariableFeatures(SRRXL, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(SRRXL)
SRRXL<- ScaleData(SRRXL, features = all.genes)
#memory clean up to save ram
gc()
SRRXL<- RunPCA(SRRXL, features = VariableFeatures(object = SRRXL))
#memeory clean up to save ram
gc()
SRRXL<- FindNeighbors(SRRXL, dims = 1:10)
SRRXL<- FindClusters(SRRXL, resolution = 0.5)
SRRXL<- RunUMAP(SRRXL, dims = 1:10)
#save processed object
saveRDS(SRRXL, file = "SRRXLumap.rds")