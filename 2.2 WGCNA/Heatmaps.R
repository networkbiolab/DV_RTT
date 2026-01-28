#Get makers for the groups 
#Set working Dir
setwd("[Dir]")
print("SetWD")
#Load object
SRRXLV4xd <- readRDS("[Dir]/SRRXLv4.rds")
print("Obj Loded")
#Load Libraries
library (Seurat)#Cargar Seurat
library(dplyr) #Carga dlypr para cálculos
library(cowplot)#Cargar cowplot to gráficos
library(readr)#Cargar en
library(ggplot2)
print("Libraries loaded")
#If needed Upd8
#SRRXL<- UpdateSeuratObject(SRRXLV4xd)
#print("update Sv4 -> Sv5")
table(SRRXL$Age)
#Scale
all.genes <- rownames(SRRXL)
SRRXL<- ScaleData(SRRXL, features = all.genes)
#Create compound label of Age+Condition
SRRXL$CA<-paste(SRRXL$Condition, SRRXL$Age, sep = "_")
table(SRRXL$CA)
#RTT_14w 24578,RTT_8w 16023,RTT_10w 11616, WT_14w 24736,WT_8w 17646, WT_10w 19673.
#Check cell types labels
table(SRRXL$Ctype)
#*Cancer cells 2172,Dopaminergic neurons 12130,GABAergic neurons 8765,Immature neurons 11415,Immune system cells 3492,Mature neurons 25539,Neural Progenitor cells  11156 ,Neuroblasts 6380,Radial glial cells 20328,Unknown 12781                                                                           
#Load Full list of Eigengenes and Create Heatmap for Radial Glial Cells
Rad <- c("CENPF","UBE2C","MKI67","NUSAP1","TPX2","TOP2A","ASPM","CENPE","PTTG1","HMGB2","FAM162A","BNIP3","IGFBP2","HILPDA","FTL","VEGFA","PGK1","NRN1","ΕΝO1","ADM")
RadGXL<- subset(SRRXL, Ctype=="Radial glial cells")
DoHeatmap(RadGXL,features = Rad,group.by = "CA")
print("Rad Done")
#Load Full list of Eigengenes and Create Heatmap for Dopaminergic cells
Dp <- c("MT-CO3","MT-CO2","MT-CO1","MT-CYB","MT-ND3","MT-ATP6","TMSB4X","MT-ND2","MT-ND3","MT-ND4","MT-ND5","B2M","PAX6","AP1S2","SYNE2","ASCL1","SOX2","MCUB","ZBTB20","RPS6","RASD1")
DpXL <-subset(SRRXL, Ctype=="Dopaminergic neurons")
DoHeatmap(DpXL,features = Dp,group.by = "CA")
ggsave(filename = "DPheat.svg", plot = Dop, width = 15, height = 12, units = "in")
print("Dp Done")