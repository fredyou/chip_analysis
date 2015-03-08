BatchFoldchangeHTA <-
function(DataFile = "data.txt", ContentFile = "Content.txt", BatchInfoFile = "Condition.txt", 
    OutfileSuffixes = "_fc2", FilePath = getwd(), colname1 = "ProbeName", colname2 = "filename", prefixes = "", 
    suffixes = "_NS", suffixes2 = "_flags", Paired = FALSE, is.log = 1, is.FC = 1, envelope = 1, FC = 2) {
    # -------------Read data--------------# read Batch information file
    Conditions = ""
    if (class(BatchInfoFile) == "data.frame") {
        Conditions = BatchInfoFile
    } else {
        Conditions <- read.delim(file.path(FilePath, BatchInfoFile), sep = "\t", header = T, fill = T, check.names = FALSE)
        
    }
    # Conditions <- read.delim(file.path(FilePath, BatchInfoFile), sep = '\t', header = T, fill =
    # T,check.names=FALSE)
    BatchInfo <- subset(Conditions, Selected == "Y")
    # print
    print(paste("The number of iterations is ", as.character(nrow(BatchInfo))))
    
    
    # test colname match
    if (class(DataFile) == "character") {
        tmp = read.delim(file.path(FilePath, DataFile), nrow = 1)
    } else if (class(DataFile) == "data.frame") {
        tmp = DataFile[1, ]
    }
    coltmp = colnames(tmp)[grep("_NS", colnames(tmp))]
    listtmp = ""
    if (class(ContentFile) == "data.frame") {
        listtmp = ContentFile
        
    } else {
        listtmp = read.delim(ContentFile)
    }
    collist = paste(listtmp$name, "_NS", sep = "")
    if (!identical(sort(intersect(coltmp, collist)), sort(collist))) {
        cat("coltmp: ", sort(coltmp), "\n")
        cat("collist: ", sort(collist), "\n")
        stop("Colnames not matching listnames!\n")
    }
    
    
    
    if (!exists("dataG")) {
        print("dataG not exist!")
        if (class(DataFile) == "character") {
            data <- read.delim(file.path(FilePath, DataFile), sep = "\t", quote = "", na.strings = "", header = T, 
                fill = T, check.names = FALSE)
        } else if (class(DataFile) == "data.frame") {
            data = DataFile
        }
        assign("dataG", data, envir = .GlobalEnv)
    } else {
        print("dataG exsits!")
        coltmp = colnames(dataG)[grep("_NS", colnames(dataG))]
        if (identical(sort(intersect(coltmp, collist)), sort(collist))) {
            data = dataG
        } else {
            if (class(DataFile) == "character") {
                data <- read.delim(file.path(FilePath, DataFile), sep = "\t", quote = "", na.strings = "", header = T, 
                  fill = T, check.names = FALSE)
            } else if (class(DataFile) == "data.frame") {
                data = DataFile
            }
            assign("dataG", data, envir = .GlobalEnv)
        }
    }
    
    
    # read expression data
    probeset <- as.vector(data[[colname1]])
    
    
    # Extract the gene annotation
    data_annot <- data[, -grep(suffixes, colnames(data))]
    # data_annot <- data_annot[ , -grep(suffixes2, colnames(data_annot))]
    if (ncol(as.matrix(data_annot)) == 1) {
        is.annot = "N"
    } else {
        is.annot = "Y"
        data_annot <- as.matrix(data_annot[, -1])
    }
    
    
    for (i in 1:nrow(BatchInfo)) {
        
        if (Paired) {
            # read label file
            cat("UnPaired: ", Paired, "\n")
            label1 <- paste(prefixes, BatchInfo[i, 1], suffixes, sep = "")
            label2 <- paste(prefixes, BatchInfo[i, 2], suffixes, sep = "")
            
            # #read Flag label1Flag <- paste(prefixes, BatchInfo[i, 1], suffixes2, sep = '') label2Flag <- paste(prefixes,
            # BatchInfo[i, 2], suffixes2, sep = '')
            
            # print group
            print(paste("the NO is", as.character(i)))
            print(as.character(BatchInfo[i, 1]))
            print(as.character(BatchInfo[i, 2]))
        } else {
            print("Paired")
            # read sample file
            SampleInfo = ""
            if (class(ContentFile) == "data.frame") {
                SampleInfo = ContentFile
            } else {
                SampleInfo <- read.delim(file.path(FilePath, ContentFile), sep = "\t", header = T, fill = T, check.names = FALSE)
                
            }
            # read label file
            label1Raw <- subset(SampleInfo, Group == as.character(BatchInfo[i, 1]), select = colname2)
            label1Raw <- as.matrix(label1Raw)
            label1 <- paste(prefixes, label1Raw, suffixes, sep = "")
            
            label2Raw <- subset(SampleInfo, Group == as.character(BatchInfo[i, 2]), select = colname2)
            label2Raw <- as.matrix(label2Raw)
            label2 <- paste(prefixes, label2Raw, suffixes, sep = "")
            
            # read Flag label1Flag <- paste(prefixes, label1Raw, suffixes2, sep = ''); label2Flag <- paste(prefixes,
            # label2Raw, suffixes2, sep = '');
            
            # print group
            print(paste("the NO is", as.character(i)))
            print(as.character(BatchInfo[i, 1]))
            print(label1Raw)
            print(as.character(BatchInfo[i, 2]))
            print(label2Raw)
        }
        
        
        # ------------the Group of sample file-------------------#
        data1 <- subset(data, select = c(label1))
        data2 <- subset(data, select = c(label2))
        
        # ------------compute fold change----------#
        if (is.log == 1) {
            mean1 = apply(2^data1[label1], 1, mean, na.rm = T)
            mean2 = apply(2^data2[label2], 1, mean, na.rm = T)
        } else {
            mean1 = apply(data1[label1], 1, mean, na.rm = T)
            mean2 = apply(data2[label2], 1, mean, na.rm = T)
        }
        foldchange = mean1/mean2
        
        if (is.annot == "N") {
            Groupdata <- cbind(probeset, data1, data2, foldchange)
        } else {
            Groupdata <- cbind(probeset, data1, data2, foldchange, data_annot)
        }
        
        
        # ------------Filter procedure------------# data1Flag <- subset(Groupdata, select = label1Flag) data2Flag <-
        # subset(Groupdata, select = label2Flag) #Stat Flag condition of two group Flag1 <- stat_Flag(data1Flag) Flag2
        # <- stat_Flag(data2Flag) result<-cbind(Groupdata,Flag1,Flag2)
        result = Groupdata
        
        # Whether Filter
        if (is.FC == 1) {
            # Filter the Flag Filter <- Filter_Flag(result, Propotion = envelope) Filter the foldchange
            Filter <- subset(result, (foldchange >= FC | foldchange <= 1/FC))
            # Filter <- subset(Filter, select = -c(Flag1,Flag2))
        } else {
            Filter <- result
        }
        
        
        # -------------output Filter result-----------------#
        colnames(Filter)[1] = colname1
        outputdir <- paste(BatchInfo[i, 1], "_VS_", BatchInfo[i, 2], sep = "")
        outputdir <- file.path(FilePath, outputdir)
        dir.create(outputdir, showWarnings = FALSE)
        OutName = paste(BatchInfo[i, 1], "_VS_", BatchInfo[i, 2], OutfileSuffixes, ".txt", sep = "")
        # if(i=nrow(BatchInfo)){ diff_up_down_gene_count(datafile=Filter,filepathin=outputdir,filepathout='./') }
        write.table(Filter, file = file.path(outputdir, OutName), sep = "\t", quote = FALSE, col.names = TRUE, 
            row.names = FALSE)
    }
    
}
