overlappedGroupFunc <-
function(datafile = "Book1.txt", group = "groups.txt", BatchInfoFile = "Combination.txt", 
    func = functionA) {
    
    if (missing(func)) {
        stop("You should assign a function!!")
    } else if (class(func) != "character") {
        stop("You`ve assigned the func par, but it should be a string!")
    }
    print("Processing pars!")
    if (missing(group)) {
        group = read.delim("groups.txt", na.strings = "", stringsAsFactors = F, quote = "")
    } else if (class(group) == "character") {
        group = read.delim(group)
    } else {
        group = group
    }
    group = as.data.frame(t(apply(as.matrix(group), 1, function(x) {
        gsub("\\s+", "", x)
    })))  #remove spaces
    
    if (missing(BatchInfoFile)) {
        BatchInfoFile = read.delim("Combination.txt", stringsAsFactors = F)
    } else if (class(BatchInfoFile) == "character") {
        BatchInfoFile = read.delim(BatchInfoFile, stringsAsFactors = F)
    } else {
        BatchInfoFile = BatchInfoFile
    }
    
    
    if (missing(datafile)) {
        datafile = read.delim("Book1.txt", na.strings = "", quote = "", sep = "\t", stringsAsFactors = F)
    } else if (class(datafile) == "character") {
        datafile = read.delim(datafile, na.strings = "", quote = "", sep = "\t", stringsAsFactors = F)
    } else {
        datafile = datafile
    }
    
    
    picName = "Book1.txt"
    
    
    for (i in 1:dim(BatchInfoFile)[1]) {
        con1 = BatchInfoFile[i, 1]
        con2 = BatchInfoFile[i, 2]
        print(paste(con1, con2, sep = " VS "))
        group1 = group[, con1]
        group1 = group1[!is.na(group1)]
        group1 = gsub("(^[0-9])", "s\\1", group1)
        # group1=paste(group1,'_NS',sep='')
        
        group2 = group[, con2]
        group2 = group2[!is.na(group2)]
        group2 = gsub("(^[0-9])", "s\\1", group2)
        # group2=paste(group2,'_NS',sep='')
        group1=group1[group1!=""]
        group2=group2[group2!=""]
        
        list = data.frame(Samples = c(group1, group2), name = c(group1, group2), Group = c(rep(con1, length(group1)), 
            rep(con2, length(group2))))
        combination = BatchInfoFile[i, ]
        picName = paste(con1, con2, sep = "_")
        # datafile=data hclust_plot(datafile=datafile,contentfile =
        # list,suffixes2='_NS',filepathin='./',filepathout='./',picName=picName) PCA_plot(datafile=datafile,contentfile
        # = list,suffixes2='_NS',filepathin='./',filepathout='./',picName=picName)
        if (exists("list") && class(list) == "data.frame") {
            print("Got list !")
        } else {
            stop("Something wrong with list generation!!")
        }
        
        if (exists("combination") && class(combination) == "data.frame") {
            print("Got combination !")
        } else {
            stop("Something wrong with combination generation!!")
        }
        eval(parse(text = func))
        
        
    }
}
