#Monocle 3 Split
 SRRXLV4xd$Ctype[SRRXLV4xd$Ctype == "Cancer cells"] <- "Neuroepithelial cell*" 
 # Createstable of counts
 ftable(SRRXLV4xd$monocle3_clusters, 
        SRRXLV4xd$Condition, 
        SRRXLV4xd$Ctype)
 
 # 1. Create the summary table
 library(dplyr)
 summary_data <- SRRXLV4xd@meta.data %>%
   # removes NA rows
   filter(!is.na(monocle3_clusters), !is.na(Condition), !is.na(Ctype)) %>%
   group_by(monocle3_clusters, Condition, Ctype) %>%
   summarise(Cell_Count = n(), .groups = 'drop')
 
 # Save the csv
 write.csv(summary_data, "Seurat_Cluster_Summary_Clean.csv", row.names = FALSE)
 plot <- ggplot(summary_data, aes(x = Condition, y = Cell_Count, fill = Ctype)) +
   # Create the bars
   geom_bar(stat = "identity", position = "stack") +
   # Split the view 
   facet_wrap(~monocle3_clusters, scales = "free") +
   theme_bw() +
   labs(title = "Cell Type Composition per Cluster and Condition",
        y = "Number of Cells",
        x = "Condition",
        fill = "Cell Type") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Tilt x labels if needed
 
 # View the plot
 print(plot)
 
