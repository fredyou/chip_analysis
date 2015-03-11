batch_plotGO_KEGG <-
function(files = list.files(pattern = ".*VS.*GO.txt"), path = ".", type = "GO") {
    if (missing(path)) {
        path = "."
    } else {
        path = path
    }
    setwd(path)
    colTerm = ""
    if (!missing("type")) {
        type = type
    } else {
        type = "GO"
    }
    stopifnot(grepl("GO*|KEGG*", type, ignore.case = T))
    
    
    if (missing(files)) {
        if (type == "GO") {
            files = list.files(pattern = paste(".*VS.*", type, ".txt$", sep = ""), ignore.case = T)
        } else {
            files = list.files(pattern = paste(".*VS.*", type, ".*.xls$", sep = ""), ignore.case = T)
        }
    } else {
        files = files
    }
    stopifnot(!is.element(TRUE, file.info(files)$isdir))
    
    for (i in 1:length(files)) {
        sample1 = gsub("([^_]*)_VS.*", "\\1", files[i], perl = T)
        sample2 = gsub(paste("[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_.*)", type, ".*.(txt|xls)", sep = ""), "\\1", 
            files[i], perl = T, ignore.case = T)  #sample2=gsub('[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_)?GO.txt','\\1',files[i],perl=T,ignore.case=T)
        plot_go_keggEnrich(files[i], sample1, sample2, type = type)
        # plot_kegg(files[i],sample1,sample2)
        
    }
    
}
