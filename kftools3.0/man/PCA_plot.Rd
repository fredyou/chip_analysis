\name{PCA_plot}
\alias{PCA_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PCA_plot
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PCA_plot(datafile = "Book1.txt", contentfile = "list.tsv", suffixes2 = "_NS", filepathin = "./", filepathout = "./", picName = "g1_vs_g2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{contentfile}{
%%     ~~Describe \code{contentfile} here~~
}
  \item{suffixes2}{
%%     ~~Describe \code{suffixes2} here~~
}
  \item{filepathin}{
%%     ~~Describe \code{filepathin} here~~
}
  \item{filepathout}{
%%     ~~Describe \code{filepathout} here~~
}
  \item{picName}{
%%     ~~Describe \code{picName} here~~
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
function (datafile = "Book1.txt", contentfile = "list.tsv", suffixes2 = "_NS", 
    filepathin = "./", filepathout = "./", picName = "g1_vs_g2") 
{
    group_tb = ""
    if (class(contentfile) == "character") {
        group_tb <- read.delim(paste(filepathin, contentfile, 
            sep = ""), sep = "\t", quote = "", header = T, fill = T, 
            stringsAsFactors = F)
    }
    else {
        group_tb = contentfile
        contentfile = paste(contentfile[1, 3], "_VS_", contentfile[2, 
            3])
    }
    tb1 = ""
    if (class(datafile) == "character") {
        tb1 <- read.delim(paste(filepathin, datafile, sep = ""), 
            sep = "\t", quote = "", header = T, fill = T, stringsAsFactors = F)
    }
    else {
        tb1 = datafile
        datafile = contentfile
    }
    tb2 <- tb1[, grep(suffixes2, colnames(tb1))]
    colnames(tb2) <- gsub(suffixes2, "", colnames(tb2))
    fit.pca <- princomp(tb2)
    var_prop <- fit.pca$sdev^2/sum(fit.pca$sdev^2)
    sam_names <- rownames(fit.pca$loadings)
    sam_colors <- as.factor(group_tb[match(rownames(fit.pca$loadings), 
        group_tb[, "name"]), "Group"])
    pch_pca <- as.character(sam_colors)
    for (i in unique(pch_pca)) {
        pch_pca[!is.na(match(pch_pca, i))] <- 1:length(which(!is.na(match(pch_pca, 
            i))))
    }
    pch_pca <- as.integer(pch_pca)
    if (!missing(picName)) {
        datafile = picName
    }
    else {
        if (class(datafile) != "character") {
            datafile = "Book1.txt"
        }
    }
    pdf(file = paste(filepathout, datafile, "_PCA_plot_2D.pdf", 
        sep = ""))
    par(xpd = NA, oma = c(0, 0, 0, 5))
    plot(fit.pca$loadings[, 1:2], xlab = paste("Comp 1", round(var_prop[1], 
        3), sep = " : "), ylab = paste("Comp 2", round(var_prop[2], 
        3), sep = " : "), col = sam_colors, cex = 1, pch = pch_pca, 
        main = gsub(".txt", "", datafile))
    legend(par("usr")[2] + (max(fit.pca$loadings[, 1]) - min(fit.pca$loadings[, 
        1]))/50, par("usr")[4], legend = sam_names, col = sam_colors, 
        pch = pch_pca, xjust = 0, yjust = 1, ncol = 1, cex = min(1, 
            8/max(nchar(sam_names))))
    dev.off()
    png(file = paste(filepathout, datafile, "_PCA_plot_2D.png", 
        sep = ""), height = 14, width = 14, units = "cm", res = 600)
    par(xpd = NA, oma = c(0, 0, 0, 5))
    plot(fit.pca$loadings[, 1:2], xlab = paste("Comp 1", round(var_prop[1], 
        3), sep = " : "), ylab = paste("Comp 2", round(var_prop[2], 
        3), sep = " : "), col = sam_colors, cex = 1, pch = pch_pca, 
        main = gsub(".txt", "", datafile))
    legend(par("usr")[2] + (max(fit.pca$loadings[, 1]) - min(fit.pca$loadings[, 
        1]))/50, par("usr")[4], legend = sam_names, col = sam_colors, 
        pch = pch_pca, xjust = 0, yjust = 1, ncol = 1, cex = min(1, 
            8/max(nchar(sam_names))))
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
