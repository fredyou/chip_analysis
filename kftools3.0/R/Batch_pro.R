Batch_pro <-
function(contentfile, filepathin, filepathout, Num, headline) {
    content <- read.table(contentfile, header = T)
    
    for (i in 1:Num) {
        d1 <- readLines(paste(filepathin, content[[2]][i], ".txt", sep = ""))
        d2 = d1[headline:length(d1)]
        writeLines(d2, paste(filepathout, content[[2]][i], ".txt", sep = ""), sep = "\n")
    }
}
