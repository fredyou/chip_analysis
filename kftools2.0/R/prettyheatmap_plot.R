prettyheatmap_plot <-
function(datafile = c("CS1_VS_M1_p005fc2.txt", "CS1_VS_M1_p001fc2.txt"), list = "list.tsv", 
    anno = FALSE) {
    cellheight = 12
    cellwidth = 50
    
    require(pheatmap)
    for (f_name in datafile) {
        data <- read.table(paste(f_name, sep = ""), sep = "\t", quote = "\"", header = T, fill = T, 
            stringsAsFactors = F)
        if (is.na(data[1, 1] == TRUE)) {
            print(paste(f_name, " is actually empty!!"))
            next
        }
        
        test = as.matrix(data[, grep("NS", colnames(data))])
        rownames(test) = data[, grep("probe|gene|systematic", ignore.case = T, colnames(data))[1]]
        
        if(nrow(test)>60){
          cellheight=1
          rownames(test)=NULL
        }
        
        if (!missing(anno) && anno == "TRUE") {
            # annotation
            if (missing(list)) {
                stop("You cannot add annotation with list.tsv!!")
            } else if (class(list) == "character") {
                list <- read.delim(list, sep = "\t", header = T, stringsAsFactors = F)
            } else if (class(list) == "data.frame") {
                list = list
            }
            rownames(list) = list$name
            datacols = colnames(test)[grep("_NS", colnames(test))]
            datacols = gsub("_NS", "", datacols)
            groups = list[datacols, 3]
            annotation = data.frame(Group = groups)
            rownames(annotation) = colnames(test)
            color = colorRampPalette(c("blue", "green", "yellow", "red"))(length(unique(groups)))
            names(color) = c(as.character(unique(annotation$Group)))
            color = list(Group = color)
            pheatmap(test, scale = "row", clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", 
                method = "average", color = colorRampPalette(c("green", "black", "red"))(50), cellwidth = cellwidth, 
                cellheight = cellheight, main = "heatmap", fontsize = 12, cluster_rows = T, cluster_cols = T, filename = paste(f_name, 
                  "_heatmap_plot.pdf", sep = ""), annotation = annotation, annotation_colors = color)
        } else {
            pheatmap(test, scale = "row", clustering_distance_cols = "correlation", clustering_distance_rows = "correlation", 
                method = "average", color = colorRampPalette(c("green", "black", "firebrick3"))(50), cellwidth = cellwidth, 
                cellheight = cellheight, filename = paste(f_name, "_heatmap_plot.pdf", sep = ""))
        }
    }  #end of for
}
