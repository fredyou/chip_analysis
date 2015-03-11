Filter_Flag <-
function(data, Propotion = 1) {
    
    Is = numeric(nrow(data))
    
    for (i in 1:nrow(data)) {
        if (data$foldchange[i] >= 1 & data$Flag1[i] >= Propotion) 
            Is[i] = 1 else if (data$foldchange[i] < 1 & data$Flag2[i] >= Propotion) 
            Is[i] = 1
    }
    data <- cbind(data, Is)
    new <- subset(data, Is == 1, select = -Is)
}
