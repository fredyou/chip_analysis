Choose_group <-
function(contentfile = "content.txt", filepath = "", colname = "GroupID", fileoutname = "Condition.txt") {
    
    d1 <- read.table(paste(filepath, contentfile, sep = ""), header = TRUE)
    Group <- levels(d1[[colname]])
    
    maxiter <- length(Group) - 1
    for (i in 1:maxiter) {
        temp <- matrix("", length(Group) - i, 3)
        for (j in 1:nrow(temp)) {
            temp[j, 1] <- Group[i]
            temp[j, 2] <- Group[i + j]
            temp[j, 3] <- "Y"
        }
        if (i == 1) {
            Condition <- temp
        } else {
            Condition <- rbind(Condition, temp)
        }
    }
    Condition.names <- c("Condition1", "Condition2", "Selected")
    write.table(Condition, file = paste(filepath, fileoutname, sep = ""), sep = "\t", quote = FALSE, col.names = Condition.names, 
        row.names = FALSE)
}
