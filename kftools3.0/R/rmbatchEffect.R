rmbatchEffect <-
function(datafile = "Book1.txt", phenodata = "pheno.tsv") {
    # pheno.tsv header: Samples sample outcome batch syndrome
    pheno1 = read.delim(phenodata, stringsAsFactors = FALSE, quote = "")
    row.names(pheno1) = pheno1[, 1]
    pheno1 = pheno1[, -1]
    batch1 = as.numeric(pheno1$batch)
    modcombat1 = model.matrix(~1, data = pheno1)
    
    # read data
    temp = read.delim(datafile, nrow = 2, stringsAsFactors = FALSE, quote = "")
    # colsNotNeeded=grep('_NS|Probe|systematic|Gene(_)?symbol',colnames(temp),ignore.case=T,invert=T)
    colsNotNeeded = grep("_NS|Probe|systematic", colnames(temp), ignore.case = T, invert = T)
    classes = sapply(temp, class)
    classes[colsNotNeeded] = "NULL"
    data = read.delim("Book1.txt", stringsAsFactors = FALSE, quote = "", colClasses = classes, comment.char = "")
    row.names(data) = data[, 1]
    data = as.matrix(data[, -1])
    
    # combat!!
    combat_edata = ComBat(dat = data, batch = batch1, mod = modcombat1, numCovs = NULL, par.prior = TRUE, prior.plots = TRUE)
    combat_edata = as.data.frame(combat_edata)
    combat_edata = data.frame(ProbeName = row.names(combat_edata), combat_edata)
    write.table(combat_edata, file = "nordata_after_rmBatch.txt", sep = "\t", row.names = F)
}
