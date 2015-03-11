scatterplot_NoFlag <-
function(combination = "Combination.txt", list = "list.tsv", datafile = "Book1.txt", identify = "hsa-miR-362-5p") {
    if (missing(combination)) {
        combination <- read.delim("Combination.txt", sep = "\t", header = T, stringsAsFactors = F)
    } else if (class(combination) == "character") {
        combination <- read.delim(combination, sep = "\t", header = T, stringsAsFactors = F)
    } else if (class(combination) == "data.frame") {
        combination = combination
    }
    if (missing(list)) {
        list <- read.delim("list.tsv", sep = "\t", header = T, stringsAsFactors = F)
    } else if (class(list) == "character") {
        list <- read.delim(list, sep = "\t", header = T, stringsAsFactors = F)
    } else if (class(list) == "data.frame") {
        list = list
    }
    if (missing(datafile)) {
        tb1 <- read.delim("Book1.txt", sep = "\t", header = T, stringsAsFactors = F)
    } else if (class(datafile) == "character") {
        tb1 <- read.delim(datafile, sep = "\t", header = T, stringsAsFactors = F)
    } else if (class(datafile) == "data.frame") {
        tb1 = datafile
    }
    # list <- read.delim(list,sep='\t',header=T,stringsAsFactors=F) tb1 <-
    # read.delim(datafile,sep='\t',header=T,stringsAsFactors=F)
    rownames(tb1) = tb1[, 1]
    ## for combination
    for (i in 1:nrow(combination)) {
        if (combination[i, "Selected"] == "Y") {
            
            list = within(list, {
                Group = as.character(Group)
                name = as.character(name)
            })
            combination = within(combination, {
                Condition1 = as.character(Condition1)
                Condition2 = as.character(Condition2)
            })
            g1_names <- list[list[, "Group"] == combination[i, "Condition1"], "name"]
            g2_names <- list[list[, "Group"] == combination[i, "Condition2"], "name"]
            g1_names = paste(g1_names, "_NS", sep = "")
            g2_names = paste(g2_names, "_NS", sep = "")
            
            g1_tb <- tb1[, g1_names, drop = F]
            g2_tb <- tb1[, g2_names, drop = F]
            
            
            g1_values <- apply(g1_tb, 1, mean)
            g2_values <- apply(g2_tb, 1, mean)
            
            p_cols <- rep("gray", length(g1_values))
            p_cols[2^g1_values/2^g2_values >= 2] <- "red"
            p_cols[2^g1_values/2^g2_values <= 0.5] <- "green"
            
            p_sizes <- rep(1, length(g1_values))
            p_sizes[!(g2_values > g1_values - 1 & g2_values < g1_values + 1)] <- 1.2
            p_sizes[(g2_values > g1_values - 1 & g2_values < g1_values + 1)] <- 1
            
            
            pdf(file = paste("ScatterPlot_", combination[i, "Condition1"], "_VS_", combination[i, "Condition2"], 
                ".pdf", sep = ""))
            plot(g2_values, g1_values, col = p_cols, pch = 16, cex = p_sizes, xlab = combination[i, "Condition2"], 
                ylab = combination[i, "Condition1"])
            abline(a = 1, b = 1, col = "red")
            abline(a = -1, b = 1, col = "green")
            abline(a = 0, b = 1, col = "black")
            if (!missing(identify)) {
                text(g1_values[identify] + 0.5, g2_values[identify] + (floor(max(g2_values) - g2_values[identify] - 
                  2)), label = identify)
                arrows(x0 = g1_values[identify] + 0.5, y0 = g2_values[identify] - 0.5 + (floor(max(g2_values) - 
                  g2_values[identify] - 2)), x1 = g1_values[identify], y1 = g2_values[identify], col = 2)
            }
            ## text()
            dev.off()
            ## if
        }
        ## for combination
    }
    
}
