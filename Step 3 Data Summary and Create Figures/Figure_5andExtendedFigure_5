
library(ComplexHeatmap)
library(circlize)
cd <- read.csv("heatmap data input.csv",header = TRUE,row.names = 1)
cd <- as.matrix(cd)

mycol<- colorRampPalette(c("blue","dodgerblue","deepskyblue","white","orange","tomato","red2"))(41)
pdf(file='heatmap for all variants.pdf',width=9, height=9)
Heatmap(cd,col = mycol,cluster_rows = T,cluster_columns = T,
        column_dend_height = unit(2.5, "cm"), 
        row_dend_width = unit(2.5, "cm"),
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        row_names_side = "left",
        row_title_side = "left",
        row_dend_side = "right",
        #row_title_rot = 0,
        row_title ="Regions",   
        column_title = "Lineages",
        column_title_side = "bottom",
        column_names_rot = 45,  
        heatmap_legend_param = list(
          title= "Genetic distance", #title_position = "topleft", 
          legend_height=unit(8,"cm"), legend_direction="vertical"),
        heatmap_width = unit(18, "cm"), heatmap_height = unit(16, "cm"))
while (!is.null(dev.list()))  dev.off()

pdf(file='heatmap for Omicron variant.pdf',width=9, height=9)
Heatmap(cd,col = mycol,cluster_rows = T,cluster_columns = T,
        column_dend_height = unit(2.5, "cm"), 
        row_dend_width = unit(2.5, "cm"),
        column_names_side = "bottom",
        #column_names_gp = gpar(fontsize = 10),
        row_names_side = "left",
        row_title_side = "left",
        row_dend_side = "right",
        #row_title_rot = 0,
        row_title ="Regions",   
        column_title = "Lineages",
        column_title_side = "bottom",
        column_names_rot = 45,  
        heatmap_legend_param = list(
          title= "Genetic distance", #title_position = "topleft", 
          legend_height=unit(8,"cm"), legend_direction="vertical"),
        heatmap_width = unit(18, "cm"), heatmap_height = unit(16, "cm"))
while (!is.null(dev.list()))  dev.off()
