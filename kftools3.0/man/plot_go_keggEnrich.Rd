\name{plot_go_keggEnrich}
\alias{plot_go_keggEnrich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot_go_keggEnrich
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_go_keggEnrich(file, sample1, sample2, type = "GO")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{sample1}{
%%     ~~Describe \code{sample1} here~~
}
  \item{sample2}{
%%     ~~Describe \code{sample2} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (file, sample1, sample2, type = "GO") 
{
    colTerm = ""
    if (!missing("type")) {
        if (type == "GO") {
            colTerm = "GO_term"
        }
        else {
            colTerm = "KEGG_name"
        }
    }
    else {
        colTerm = "GO_term"
    }
    if (grepl("_VS_", file)) {
        if (missing(sample1) && missing(sample2)) {
            sample1 = gsub("([^_]*)_VS.*", "\\1", file, perl = T)
            sample2 = gsub(paste("[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_.*)", 
                type, ".*.(txt|xls)", sep = ""), "\\1", file, 
                perl = T, ignore.case = T)
        }
        else if (missing(sample1) || missing(sample2)) {
            warning("You assigned only one of sample1 or sample2,\n so we decide to deduct those two from filename\nOtherwise you can just assign two of them or none of them instead!\n")
            sample1 = gsub("([^_]*)_VS.*", "\\1", file, perl = T)
            sample2 = gsub(paste("[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_.*)", 
                type, ".*.(txt|xls)", sep = ""), "\\1", file, 
                perl = T, ignore.case = T)
        }
        else {
            sample1 = sample1
            sample2 = sample2
        }
    }
    else if (missing(sample1) || missing(sample2)) {
        stop("Your file name doesnt contain _VS_ pattern or you just lack sample1 sample2 parameter!\n")
    }
    pdf(paste(sample1, "_VS_", sample2, "_enrichment_", type, 
        ".pdf", sep = ""), width = 12, height = 10)
    par(mai = c(0, 5, 1, 0.3))
    go_data <- read.delim(file, sep = "\t", quote = "", header = T)
    go_data <- go_data[order(go_data[, 13], decreasing = F), 
        1:13]
    go_part <- go_data[40:1, ]
    score <- as.numeric(-1 * log10(go_part["FDR"][, 1]))
    barplot(score, names.arg = go_part[colTerm][, 1], horiz = T, 
        col = 1, las = 2, axes = F, cex.names = 0.5)
    axis(3)
    mtext("-Lg(FDR)", side = 3, line = 2.5, font = 2, cex = 1.5)
    mtext(colTerm, side = 2, line = 22, font = 2, cex = 1.5)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
