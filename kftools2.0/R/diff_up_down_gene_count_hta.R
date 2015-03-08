diff_up_down_gene_count_hta <-
function(datafile = c("CS1_VS_M1_p005fc2.txt", "CS1_VS_M1_p001fc2.txt"), filepathin = "./", 
    filepathout = "./") {
    tb_count <- NULL
    for (f_name in datafile) {
        tb1 <- read.csv(paste(filepathin, f_name, sep = ""), sep = "\t", header = T, stringsAsFactors = F)
        diff_count <- nrow(tb1)
        up_count <- nrow(tb1[tb1[, "foldchange"] > 1, ])
        down_count <- nrow(tb1[tb1[, "foldchange"] < 1, ])
        tmp_li <- cbind(file_name = f_name, diff_count = diff_count, up_count = up_count, down_count = down_count)
        tb_count <- rbind(tb_count, tmp_li)
    }
    ## comparison
    comparison <- paste(sapply(tb_count[, "file_name"], function(x) {
        strsplit(x, "_")[[1]][1]
    }), sapply(tb_count[, "file_name"], function(x) {
        strsplit(x, "_")[[1]][2]
    }), sapply(tb_count[, "file_name"], function(x) {
        strsplit(x, "_")[[1]][3]
    }))
    comparison <- gsub(".txt", "", comparison)
    comparison[is.na(comparison)] <- ""
    ## p_value
    p_value <- sapply(tb_count[, "file_name"], function(x) {
        strsplit(x, "_p")[[1]][2]
    })
    p_value <- sapply(p_value, function(x) {
        strsplit(x, "fc")[[1]][1]
    })
    p_value <- sapply(p_value, function(x) {
        strsplit(x, ".txt")[[1]][1]
    })
    p_value <- gsub("005", "0.05", p_value)
    p_value <- gsub("001", "0.01", p_value)
    p_value[is.na(p_value)] <- ""
    ## fc
    fc <- sapply(tb_count[, "file_name"], function(x) {
        strsplit(x, "fc")[[1]][2]
    })
    fc <- gsub(".txt", "", fc)
    fc[is.na(fc)] <- ""
    ## 
    tb_count_2 <- cbind(file_name = tb_count[, "file_name"], comparison = comparison, p_value = p_value, fc = fc, 
        diff_count = tb_count[, "diff_count"], up_count = tb_count[, "up_count"], down_count = tb_count[, "down_count"])
    colnames(tb_count_2) <- c("file name", "comparison", "p-value threshold", "fold change threshold", "number of differentially expression genes", 
        "number of up-regulated genes", "number of down-regulated genes")
    
    write.table(tb_count_2, file = "diff_up_down_gene_count.xls", append = T, sep = "\t", col.names = T, row.names = F, 
        quote = F)
}
