Batch_ttest <-
function(datafile = "Book1.txt", contentfile = "list.tsv", Batch_infofile = "Combination.txt", paired_ttest_fname = "paired.txt", 
    filepathin = "", filepathout = "", colname1 = "ProbeName", colname2 = "name", prefixes = "", suffixes = "_NS", 
    suffixes2 = "_flags", T.type = 2, is.log = 1, cutoff = 1, is.adjustP = 0) {
    # -------------Read data--------------# read expression file
    data = ""
    if (class(datafile) == "character") {
        data <- read.delim(datafile, na.strings = "", quote = "", sep = "\t", stringsAsFactors = F)
    } else {
        data = datafile
        rm(datafile)
    }
    print(dim(data))
    
    probeset <- as.vector(data[[colname1]])
    
    Conditions = ""
    if (class(Batch_infofile) == "data.frame") {
        Conditions = Batch_infofile
    } else {
        Conditions <- read.delim(file.path(filepathin, Batch_infofile), sep = "\t", header = T, check.names = FALSE)
        
    }
    
    
    
    Batch_info <- subset(Conditions, Selected == "Y")
    # print
    print(paste("The number of iterations is ", as.character(nrow(Batch_info))))
    
    # read sample file
    
    sampleinfo = ""
    if (class(contentfile) == "data.frame") {
        sampleinfo = contentfile
    } else {
        sampleinfo <- read.delim(file.path(filepathin, contentfile), sep = "\t", header = T, check.names = FALSE, 
            quote = "\"", stringsAsFactors = F)
        
    }
    data_annot <- data[, -grep(suffixes, colnames(data))]
    data_annot <- data_annot[, -grep(suffixes2, colnames(data_annot))]
    if (ncol(as.matrix(data_annot)) == 1) {
        is.annot = "N"
    } else {
        is.annot = "Y"
    }
    for (i in 1:nrow(Batch_info)) {
        label1_raw <- subset(sampleinfo, Group == as.character(Batch_info[i, 1]), select = colname2)
        label1_raw <- as.matrix(label1_raw)
        label1 <- paste(prefixes, label1_raw, suffixes, sep = "")
        
        label2_raw <- subset(sampleinfo, Group == as.character(Batch_info[i, 2]), select = colname2)
        label2_raw <- as.matrix(label2_raw)
        label2 <- paste(prefixes, label2_raw, suffixes, sep = "")
        
        # read Flag
        label1_Flag <- paste(prefixes, label1_raw, suffixes2, sep = "")
        label2_Flag <- paste(prefixes, label2_raw, suffixes2, sep = "")
        
        # print group
        print(paste("the NO is", as.character(i)))
        print(as.character(Batch_info[i, 1]))
        print(label1_raw)
        print(as.character(Batch_info[i, 2]))
        print(label2_raw)
        
        # ------------t-test------------#
        
        cat("label1: ", label1, "\\ label1Flag: ", label1_Flag, "\\\n")
        cat("label2: ", label2, "\\ label2Flag: ", label2_Flag, "\\\n")
        
        data1 <- subset(data, select = c(label1, label1_Flag))
        data2 <- subset(data, select = c(label2, label2_Flag))
        result <- T_test(data1, data2, T_type = T.type, islog = is.log, suff = suffixes2, isadjustP = is.adjustP)
        if (is.annot == "N") 
            result <- cbind(probeset, result) else result <- cbind(probeset, result, data_annot)
        
        result <- subset(result, pvalues <= cutoff)
        colnames(result)[1] = colname1
        
        # -------------output result-----------------# create output result directory with annotation
        # outputdir<-paste(filepathout,Batch_info[i,1],'_VS_',Batch_info[i,2],sep='') dir.create(outputdir)
        outputdir <- filepathout
        assign(paste(Batch_info[i, 1], "_VS_", Batch_info[i, 2], sep = ""), result, envir = .GlobalEnv)
        write.table(result, file = paste(outputdir, "/", Batch_info[i, 1], "_VS_", Batch_info[i, 2], ".txt", sep = ""), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    }
    
}
