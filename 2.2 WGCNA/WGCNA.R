setwd([Obj_dir]) # Set the working directory
.library(Seurat) # Load the required libraries
library(dplyr)
library(cowplot)
library(readr)
library(tidyverse)
library(patchwork)
library(WGCNA)
library(hdWGCNA) #dev.tools
library(harmony)
library(igraph)
cds <- readRDS("MonoJumboCTYPE.rds")
table(cds$Ctype)
#select Appropiate ctype
cdsi <- subset(cds, Ctype== "Neuroblasts")

# Prepare the object for processing, selecting the "fraction" parameter to indicate that 
# clustering should be done using genes only expressed in a fraction of the cells
seuratobj <- SetupForWGCNA(cdsi,gene_select = "fraction",fraction = 0.05,wgcna_name = "CtypeGab") 

# Calculate metacells based on provided metadata using the K-Nearest Neighbors (KNN) 
# algorithm, using average expression to determine groups
seuratobj <- MetacellsByGroups(seurat_obj = seuratobj,reduction = 'PCA', group.by = c("Age", "Condition"),k = 25, ident.group = 'Age') 

seuratobj <- NormalizeMetacells(seuratobj) # Normalize the metacells

# Extract the expression matrix in the optimal format for WGCNA processing
seuratobj <- SetDatExpr(seuratobj,group_name ="14w",group.by='Age') 

# Calculate a correlation adjacency matrix, determining a "power" value for interactions 
# in a scale-free topology process
seuratobj <- TestSoftPowers(seuratobj) 

powertable <- GetPowerTable(seuratobj) # Create a table with these scores
write.table(powertable, file='PowerNP100.txt',col.names = TRUE) # Export the generated table
seuratobj <- ConstructNetwork(seuratobj, soft_power=9) # Construct the interaction network

# Calculate the modules for each cell group according to the variable of interest
seuratobj <- ModuleEigengenes(seuratobj,group.by.vars="orig.ident") 

hMEs <- GetMEs(seuratobj) # Create an object containing the harmonized modules
MEs <-