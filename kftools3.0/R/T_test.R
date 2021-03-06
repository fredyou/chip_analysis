T_test <-
function(data1, data2, T_type = 1, islog = 1, suff = "_flags", isadjustP = 0) {
    
    if (nrow(data1) != nrow(data2)) {
        stop("The dim of two groups are not the same!\n")
    }
    
    t_statistic <- numeric(nrow(data1))
    pvalues <- numeric(nrow(data1))
    mean1 <- numeric(nrow(data1))
    mean2 <- numeric(nrow(data1))
    median1 <- numeric(nrow(data1))
    median2 <- numeric(nrow(data1))
    stdvar1 <- numeric(nrow(data1))
    stdvar2 <- numeric(nrow(data1))
    foldchange <- numeric(nrow(data1))
    
    # seperate the gene normalized data
    data1_flag <- data1[, grep(suff, colnames(data1))]
    data1 <- as.data.frame(data1[, -grep(suff, colnames(data1))])
    data2_flag <- data2[, grep(suff, colnames(data2))]
    data2 <- as.data.frame(data2[, -grep(suff, colnames(data2))])
    
    for (i in 1:nrow(data1)) {
        mean1[i] <- mean(as.numeric(data1[i, ]))
        mean2[i] <- mean(as.numeric(data2[i, ]))
        median1[i] <- median(as.numeric(data1[i, ]))
        median2[i] <- median(as.numeric(data2[i, ]))
        stdvar1[i] <- sd(as.numeric(data1[i, ]))
        stdvar2[i] <- sd(as.numeric(data2[i, ]))
        if (stdvar1[i] == 0 & stdvar2[i] == 0) {
            t_statistic[i] <- NA
            pvalues[i] <- 1
        } else {
            if (T_type == 0) {
                t.result <- t.test(data1[i, ], data2[i, ], var.equal = T)
            }
            if (T_type == 1) {
                t.result <- t.test(data1[i, ], data2[i, ], var.equal = F)
            }
            ## paired ttest
            if (T_type == 2) {
                t.result <- t.test(as.numeric(data1[i, ]), as.numeric(data2[i, ]), var.equal = F, paired = T)
            }
            ## 
            t_statistic[i] <- t.result$statistic
            pvalues[i] <- t.result$p.value
        }
        
        if (islog == 1) {
            x <- 2^data1[i, ]
            y <- 2^data2[i, ]
            foldchange[i] <- mean(x)/mean(y)
        } else {
            foldchange[i] <- mean2[i]/mean1[i]
        }
    }
    
    # -------------output result-----------------# result <-
    # cbind(data1,data1_flag,data2,data2_flag,t_statistic,pvalues,foldchange,mean1,mean2,median1,median2,stdvar1,stdvar2);
    # whether computes adjusted p-values
    if (isadjustP == 1) {
        fdr <- p.adjust(pvalues, "bonferroni")
        result <- cbind(pvalues, foldchange, fdr, data1, data1_flag, data2, data2_flag)
    } else {
        result <- cbind(pvalues, foldchange, data1, data1_flag, data2, data2_flag)
    }
    
}
