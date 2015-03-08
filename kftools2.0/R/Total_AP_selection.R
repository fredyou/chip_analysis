Total_AP_selection <-
function(datafile = "Book1.txt", contentfile = "list.tsv", suffixes2 = "_flags", filepathin = "", 
    filepathout = "./") {
    # group_tb <- read.delim(file=paste(filepathin,contentfile,sep=''),sep='\t',quote='\'',header=T,fill=T);
    group_tb = ""
    if (class(contentfile) == "data.frame") {
        group_tb = contentfile
    } else {
        group_tb <- read.delim(contentfile, sep = "\t", header = T, fill = T, check.names = FALSE)
        
    }
    
    outname = ""
    # tb1 <- read.csv(paste(filepathin,datafile,sep=''),sep='\t',header=T,stringsAsFactors=F,check.names=FALSE)
    if (class(datafile) == "character") {
        tb1 <- read.delim(datafile, sep = "\t", quote = "\"", header = T, fill = T)
        outname = datafile
    } else {
        tb1 = datafile
        rm(datafile)
        outname = "Book1.txt"
    }
    
    flags_tb <- tb1[, grep(suffixes2, colnames(tb1))]
    colnames(flags_tb) <- gsub(suffixes2, "", colnames(flags_tb))
    
    total_group <- as.character(unique(group_tb[, "Group"]))
    all_flags <- rep(FALSE, nrow(flags_tb))
    for (group_i in total_group) {
        sample_i <- as.character(group_tb[grep(group_i, group_tb[, "Group"], fixed = T), "name"])
        print(sample_i)
        print(group_i)
        tmp_tb <- as.data.frame(flags_tb[, match(sample_i, colnames(flags_tb))])
        tmp_flags <- apply(tmp_tb, 1, function(x) {
            length(grep("A", x)) == 0
        })
        all_flags <- all_flags | tmp_flags
    }
    AP_selection_tb <- tb1[all_flags, ]
    ## 
    
    write.table(AP_selection_tb, file = paste(basename(outname), "_AP_selection.xls", sep = ""), sep = "\t", quote = FALSE, 
        col.names = TRUE, row.names = FALSE)
}
