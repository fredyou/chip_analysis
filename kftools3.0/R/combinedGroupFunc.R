combinedGroupFunc <-
function(datafile = "Book1.txt", group = "group.txt", combinedGroup = TRUE, func = functionA) {
    if (missing(func)) {
        stop("You should assign a function!!")
    } else if (class(func) != "character") {
        stop("You`ve assigned the func par, but it should be a string!")
    }
    print("Processing pars!")
    if (missing(group)) {
        group = read.delim("group.txt")
    } else if (class(group) == "character") {
        group = read.delim(group)
    } else {
        group = group
    }
    group = as.data.frame(t(apply(as.matrix(group), 1, function(x) {
        gsub("\\s+", "", x)
    })))  #remove spaces
    
    if (missing(combinedGroup)) {
        combinedGroup = TRUE
    } else if (class(combinedGroup) == "logical") {
        combinedGroup = combinedGroup
    } else {
        stop("combinedGroup should be logical!")
    }
    
    list_tb = ""
    if (!missing(combinedGroup) && combinedGroup == TRUE) {
        list_tb = read.delim("list.tsv")
    }
    if (missing(datafile)) {
        datafile = read.delim("Book1.txt", na.strings = "", quote = "", sep = "\t", stringsAsFactors = F)
    } else if (class(datafile) == "character") {
        datafile = read.delim(datafile, na.strings = "", quote = "", sep = "\t", stringsAsFactors = F)
        
    } else {
        datafile = datafile
    }
    
    picName = ""  #to change the pic name of PCA COR etc...
    
    
    print("looping the Comparisions!!!!")
    for (i in 1:dim(group)[1]) {
        group1 = unlist(strsplit(as.character(group[i, 1]), ",|;"))
        group2 = unlist(strsplit(as.character(group[i, 2]), ",|;"))
        combination = data.frame(Condition1 = paste(group1, collapse = "_"), Condition2 = paste(group2, collapse = "_"), 
            Selected = "Y")
        list = ""
        # combinedGroup means some groups will be combined into a larger group
        list_tbtemp = ""
        if (combinedGroup) {
            list_tbtemp = list_tb[list_tb$Group %in% c(group1, group2), ]
            list_tbtemp$Group = as.character(list_tbtemp$Group)
            for (k in 1:2) {
                groupTemp = eval(parse(text = paste("group", k, sep = "")))
                altername = paste(groupTemp, sep = "", collapse = "_")
                list_tbtemp[list_tbtemp$Group %in% groupTemp, "Group"] = altername
            }
            list = list_tbtemp
        } else {
            list = data.frame(Samples = c(group1, group2), name = c(group1, group2), Group = c(rep(paste(group1, 
                collapse = "_"), length(group1)), rep(paste(group2, collapse = "_"), length(group2))))
        }
        picName = paste(paste(group1, collapse = "_"), "_VS_", paste(group2, collapse = "_"), sep = "")
        
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
        # Function Bit
        # batch_ttest_ALL(datafile=datafile,list_tb=list,combination=combination,chipType='mirna',picName=picName)
        eval(parse(text = func))
        
        
        # 
    }  #end of for
}
