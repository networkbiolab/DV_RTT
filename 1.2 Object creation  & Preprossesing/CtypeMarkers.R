#Ctype Markers
setwd("[Home-dir]")
library (Seurat)#Cargar Seurat
library(SeuratObject)
SRRXLV4xd <- readRDS("[Home-dir]/SRRXLV4xd.rds")
CancerCll <- subset(SRRXLV4xd,ctype =="Cancer cells")
NeuroStemark <- c("ASCL1","BCL11B","CD24","CD44","CDH1","CDH2","CUX1","CUX2","CXCR4","EMX2","EOMES","EPCAM","FEZF2","FUT4","GFAP","GLI1","GLI2","GLI3","HES1","HES3","HES5","HEY1","HEY2","ID4","ITGA4","ITGB1","L1CAM","LMO4","MAP2","MCAM","MKI67","NES","NEUROG2","NKX2-2","NKX6-2","NOTCH1","OLIG1","OLIG2","PAX6","PDGFA","PROM1","SATB2","SOX1","SOX10","SOX2","SOX3","SOX5","SOX9","Svet1","TBR1","TUBB3","VIM")
VlnPlot(CancerCll, features = NeuroStemark, slot = "counts", log = TRUE)
table(SRRXLV4xd$ctype)
#save plot (38 x 28)
DopaCll<- subset(SRRXLV4xd,Ctype =="Dopaminergic neurons")
Dopamark <- c("TH","SLC6A3","FOXA2","KCNJ6","NR4A2","LMX1B","CACNA2D2","CADPS2","CHRNA6","MAPK8IP2","NTN1","PRKCG","SCN2A","TENM1","ZIM3","DDC","SLC18A2","SLC18A3","NEUROD6","SMAD3","PITX3")
VlnPlot(DopaCll, features = Dopamark, group.by = monocle3_clusters, slot = "counts", log = TRUE)
#interneur
Unkw4 <- subset(SRRXLV4xd,monocle3_clusters =="4")
Unkw19 <- subset(SRRXLV4xd,monocle3_clusters =="19")
Intermark <- c("CALB1","NPY","NOS1","SST","GNG4","ABAT","EML5","ISCA1","NAT14","SLC12A5","FRYL","PDE4A","SURF1","KCNC1","RAB3B","ABCC5","VLDLR","PIANP","ABR","NRN1","NXPH4","CDH4","EOMES","CASK","CACNA1E","ARHGAP21","CAMKV","FRRS1L","ACSL4","NPTXR","LINGO1","HPCAL4","LMO3","CHN1","SLC17A7","NCALD","FBXW11","ADCY1","FBXO11","IGSF21","SRGAP3","SCAMP1","PSD3","RORA","DCAF7","EID2","RGMA","POU2F2","NBEA","SRGAP1","ATG12","AGAP1","TMEM163","ZDHHC17","NEFM","L1CAM","ID5","DNM3","DISP2","PCSK2","AAK1","CDKL5","CCDC136","KCNMA1","PTPN4","PDP1","GABBR2","AKAP8L","SPTY2D1","PCF11","TCEAL3","KCNQ5","PAX6","PDE1A","DLX2","DLX5","ARX")
VlnPlot(Unkw4, features = Intermark, slot = "counts",log = TRUE)
VlnPlot(Unkw19, features = Intermark, slot = "counts",log = TRUE)
 
