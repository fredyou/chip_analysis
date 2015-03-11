plot_go_keggEnrich <-
function(file, sample1, sample2, type = "GO") {
    colTerm = ""
    if (!missing("type")) {
        if (type == "GO") {
            colTerm = "GO_term"
        } else {
            colTerm = "KEGG_name"
        }
    } else {
        colTerm = "GO_term"
    }
    
    if (grepl("_VS_", file)) {
        if (missing(sample1) && missing(sample2)) {
            sample1 = gsub("([^_]*)_VS.*", "\\1", file, perl = T)
            sample2 = gsub(paste("[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_.*)", type, ".*.(txt|xls)", sep = ""), "\\1", 
                file, perl = T, ignore.case = T)  #sample2=gsub('[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_)?GO.txt','\\1',files[i],perl=T,ignore.case=T)
        } else if (missing(sample1) || missing(sample2)) {
            warning("You assigned only one of sample1 or sample2,\n so we decide to deduct those two from filename\nOtherwise you can just assign two of them or none of them instead!\n")
            sample1 = gsub("([^_]*)_VS.*", "\\1", file, perl = T)
            sample2 = gsub(paste("[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_.*)", type, ".*.(txt|xls)", sep = ""), "\\1", 
                file, perl = T, ignore.case = T)  #sample2=gsub('[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_)?GO.txt','\\1',files[i],perl=T,ignore.case=T)
        } else {
            sample1 = sample1
            sample2 = sample2
        }
    } else if (missing(sample1) || missing(sample2)) {
        stop("Your file name doesnt contain _VS_ pattern or you just lack sample1 sample2 parameter!\n")
    }
    
    
    
    pdf(paste(sample1, "_VS_", sample2, "_enrichment_", type, ".pdf", sep = ""), width = 12, height = 10)
    par(mai = c(0, 5, 1, 0.3))
    go_data <- read.delim(file, sep = "\t", quote = "", header = T)
    go_data <- go_data[order(go_data[, 13], decreasing = F), 1:13]
    go_part <- go_data[40:1, ]
    score <- as.numeric(-1 * log10(go_part["FDR"][, 1]))  #go_part['FDR'][,1] == go_part[,'FDR']
    
    barplot(score, names.arg = go_part[colTerm][, 1], horiz = T, col = 1, las = 2, axes = F, cex.names = 0.5)
    axis(3)
    mtext("-Lg(FDR)", side = 3, line = 2.5, font = 2, cex = 1.5)
    mtext(colTerm, side = 2, line = 22, font = 2, cex = 1.5)
    dev.off()
}
