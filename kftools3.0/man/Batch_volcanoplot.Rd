\name{Batch_volcanoplot}
\alias{Batch_volcanoplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch_volcanoplot
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Batch_volcanoplot(Batch_infofile = "Combination.txt", pvalue_cut = 0.05, fc_cut = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Batch_infofile}{
%%     ~~Describe \code{Batch_infofile} here~~
}
  \item{pvalue_cut}{
%%     ~~Describe \code{pvalue_cut} here~~
}
  \item{fc_cut}{
%%     ~~Describe \code{fc_cut} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Batch_infofile = "Combination.txt", pvalue_cut = 0.05, 
    fc_cut = 2) 
{
    if (class(Batch_infofile) == "character") {
        batchinfo = read.delim(Batch_infofile)
    }
    else if (class(Batch_infofile) == "data.frame") {
        batchinfo = Batch_infofile
    }
    else {
        stop("Hell! Where is the combination!?")
    }
    if (missing(pvalue_cut)) {
        pvalue_cut = 0.05
    }
    else {
        pvalue_cut = pvalue_cut
    }
    cat("pvalue_cut: ", pvalue_cut, "\n")
    if (missing(fc_cut)) {
        fc_cut = 2
    }
    else {
        fc_cut = fc_cut
    }
    cat("fc_cut: ", fc_cut, "\n")
    for (i in 1:dim(batchinfo)[1]) {
        filename = paste(batchinfo$Condition1[i], "_VS_", batchinfo$Condition2[i], 
            ".txt", sep = "")
        data = read.delim(filename, na.strings = "", quote = "", 
            sep = "\t", stringsAsFactors = F)
        fileout = gsub(".txt", ".png", filename)
        pvalues = data$pvalues
        fc = data$foldchange
        df = data.frame(pvalues, fc)
        df$threshold = rep("grey", dim(df)[1])
        df$threshold[abs(df$fc) > fc_cut & df$pvalues < pvalue_cut] = "red"
        df$threshold[abs(df$fc) < 1/fc_cut & df$pvalues < pvalue_cut] = "green"
        color = df$threshold
        xaxisMax = ceiling(max(abs(log2(df$fc))))
        x <- log2(df$fc)
        y <- -log10(df$pvalues)
        b <- seq(from = 0, to = ceiling(max(y)))
        require(ggplot2)
        g = ggplot(data = df, aes(x = log2(fc), y = -log10(pvalues), 
            colour = threshold)) + geom_point(alpha = 0.8, size = 3.5, 
            colour = color) + theme(legend.position = "none") + 
            xlab("log2 fold change") + ylab("-log10 p-value") + 
            ggtitle(fileout) + geom_vline(xintercept = c(log2(fc_cut), 
            log2(1/fc_cut)), linetype = "dotted") + geom_hline(yintercept = (-log10(pvalue_cut)), 
            linetype = "dotted") + theme_classic() + scale_y_continuous(breaks = b, 
            limits = c(0, max(b)))
        ggsave(file = fileout, plot = g)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
