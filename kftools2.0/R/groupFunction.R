groupFunction <-
function(file = "Book1.txt", list = "list.tsv", group = "NS", func = "func") {
    data = read.delim(file)
    gt = data[, c(1, grep(group, colnames(data)))]
    list = read.delim(list)
    list = list[, 2:3]
    rownames(list) = list$name
    list[, 2] = as.character(list[, 2])
    list[, 1] = as.character(list[, 1])
    
    require(reshape2)
    gt.m = melt(gt)
    gt.m[, 2] = as.character(gt.m[, 2])
    gt.m[, 2] = gsub("_NS", "", gt.m[, 2])
    gt.m[, 2] = list[gt.m[, 2], 2]
    
    method = ""
    if (!missing(func)) {
        method = func
    } else {
        method = function(x) {
            log2(mean(2^x))
        }
    }
    if (class(func) != "function") {
        stop("func parameter should be set to a function!!")
    }
    
    ag.gt = as.data.frame(aggregate(gt.m[, 3], list(probe = gt.m[, 1], sample = gt.m$variable), method))
    ag.gt = as.data.frame(acast(ag.gt, probe ~ sample, mean, margins = TRUE))
    # acast convert back aggregated data to dataframe, probe is the rows, sample is the columns, mean is an added
    # column margins =TRUE will add colnames to df
    ag.gt = ag.gt[-dim(ag.gt)[1], -grep("(all)", colnames(ag.gt))]
    ag.gt$Probes = rownames(ag.gt)
    ag.gt = ag.gt[, sort(colnames(ag.gt))]
    write.table(ag.gt, file = paste(group, ".txt", sep = ""), sep = "\t", row.names = F)
}
