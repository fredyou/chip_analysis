stat_Flag <-
function(Flagdata) {
    
    stat_Flag = numeric(nrow(Flagdata))
    
    for (i in 1:nrow(Flagdata)) {
        sum_PM = 0
        for (j in 1:ncol(Flagdata)) {
            if (Flagdata[i, j] != "A") 
                sum_PM = sum_PM + 1
        }
        stat_Flag[i] = sum_PM/ncol(Flagdata)
    }
    
    # -------------output stat_Flag-----------------#
    stat_Flag
    
}
