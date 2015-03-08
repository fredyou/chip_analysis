Batch_volcanoplot <-
function(Batch_infofile = "Combination.txt", pvalue_cut = 0.05, fc_cut = 2) {
    
    if (class(Batch_infofile) == "character") {
        batchinfo = read.delim(Batch_infofile)
    } else if (class(Batch_infofile) == "data.frame") {
        batchinfo = Batch_infofile
    } else {
        stop("Hell! Where is the combination!?")
    }
    
    if (missing(pvalue_cut)) {
        pvalue_cut = 0.05
    } else {
        pvalue_cut = pvalue_cut
    }
    cat("pvalue_cut: ", pvalue_cut, "\n")
    if (missing(fc_cut)) {
        fc_cut = 2
    } else {
        fc_cut = fc_cut
    }
    cat("fc_cut: ", fc_cut, "\n")
    
    for (i in 1:dim(batchinfo)[1]) {
        filename = paste(batchinfo$Condition1[i], "_VS_", batchinfo$Condition2[i], ".txt", sep = "")
        data = read.delim(filename, na.strings = "", quote = "", sep = "\t", stringsAsFactors = F)
        fileout = gsub(".txt", ".png", filename)
        
        pvalues = data$pvalues
        fc = data$foldchange
        df = data.frame(pvalues, fc)
        df$threshold = rep("grey", dim(df)[1])
        df$threshold[abs(df$fc) > fc_cut & df$pvalues < pvalue_cut] = "red"
        df$threshold[abs(df$fc) < 1/fc_cut & df$pvalues < pvalue_cut] = "green"
        color = df$threshold
        xaxisMax = ceiling(max(abs(log2(df$fc))))
        # df$threshold = as.factor((abs(df$fc) > fc_cut|abs(df$fc)<1/fc_cut) & df$pvalues < pvalue_cut)
        # df$threshold[abs(df$fc) > fc_cut] df$threshold[grep('FALSE',df$threshold)]='grey' label range
        x <- log2(df$fc)
        y <- -log10(df$pvalues)
        # label range a<-seq(from=floor(min(x)),to=ceiling(max(x)),by=0.5)
        b <- seq(from = 0, to = ceiling(max(y)))
        # b<-seq(from=0,to=3.5,by=0.5)
        
        require(ggplot2)
        g = ggplot(data = df, aes(x = log2(fc), y = -log10(pvalues), colour = threshold)) + geom_point(alpha = 0.8, 
            size = 3.5, colour = color) + theme(legend.position = "none") + # xlim(c(-xaxisMax, xaxisMax)) + ylim(c(0, ceiling(max(-log10(pvalues))))) +
        xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle(fileout) + geom_vline(xintercept = c(log2(fc_cut), 
            log2(1/fc_cut)), linetype = "dotted") + geom_hline(yintercept = (-log10(pvalue_cut)), linetype = "dotted") + 
            theme_classic() + # scale_x_continuous(breaks=a)+
        scale_y_continuous(breaks = b, limits = c(0, max(b)))
        
        ggsave(file = fileout, plot = g)
    }
    
}
