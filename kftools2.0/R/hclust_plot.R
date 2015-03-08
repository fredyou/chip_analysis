hclust_plot <-
function(datafile = "Book1.txt", contentfile = "list.tsv", suffixes2 = "_NS", filepathin = "./", 
    filepathout = "./", picName = picName) {
    group_tb = ""
    if (class(contentfile) == "character") {
        group_tb <- read.table(paste(filepathin, contentfile, sep = ""), sep = "\t", quote = "\"", header = T, 
            fill = T)
    } else {
        group_tb = contentfile
        contentfile = paste(unique(contentfile$Group)[1], "_VS_", unique(contentfile$Group)[2])
    }
    tb1 = ""
    if (class(datafile) == "character") {
        tb1 <- read.delim(paste(filepathin, datafile, sep = ""), sep = "\t", quote = "", header = T, fill = T, 
            stringsAsFactors = F)
    } else {
        tb1 = datafile
        datafile = contentfile
    }
    expr_tb <- data.matrix(tb1[, paste(as.character(group_tb[, 2]), "_NS", sep = "")])
    colnames(expr_tb) <- gsub(suffixes2, "", colnames(expr_tb))
    # create Hierarchical cluster analysis dendogram
    d <- as.dist(dist(t(expr_tb), method = "euclidean"))
    fit <- hclust(d, method = "ward")
    # save Hierarchical cluster analysis dendogram
    sam_names <- colnames(expr_tb)
    sam_cols <- as.character(group_tb[match(fit$labels, group_tb[, "name"]), "Group"])
    names(sam_cols) <- fit$labels
    # Generate color platte
    groups = as.character(unique(group_tb$Group))
    colorCodes = colorRampPalette(c("Purple", "darkblue", "lightblue", "darkgreen", "lightgreen", "dark orange", 
        "darkred", "pink"))(length(unique(groups)))
    names(colorCodes) = unique(groups)
    d <- dendrapply(as.dendrogram(fit), labelCol, sam_cols, colorCodes)
    # dfit=as.dendrogram(fit) :contains branches as vector. dfit[1] is the leaves contained in this branch
    # d=dendrapply(dfit,colLab); pdf(file=paste(filepathout,datafile,'_hclust_plot.pdf',sep=''));
    # par(cex=min(1,par('din')[1]/(length(sam_names)*par('csi'))); #plot(fit, hang=-1) plot(d); dev.off();
    if (missing(picName)) {
        picName = datafile
    } else {
        picName = picName
    }
    png(file = paste(filepathout, picName, "_hclust_plot.png", sep = ""), height = 14, width = 14, units = "cm", 
        res = 600)
    par(cex = min(1, par("din")[1]/(length(sam_names) * par("csi"))))
    # plot(fit, hang=-1)
    plot(d)
    dev.off()
    
    pdf(file = paste(filepathout, picName, "_hclust_plot.pdf", sep = ""), height = 14, width = 14)
    par(cex = min(1, par("din")[1]/(length(sam_names) * par("csi"))))
    # plot(fit, hang=-1)
    plot(d)
    dev.off()
    
}
