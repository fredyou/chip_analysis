\name{heatmap_plot}
\alias{heatmap_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
heatmap_plot
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heatmap_plot(datafile = c("CS1_VS_M1_p005fc2.txt", "CS1_VS_M1_p001fc2.txt"), contentfile = "list.tsv", colname1 = "ProbeName", suffixes2 = "_NS", heatmap_range = 1, filepathin = "./", filepathout = "./", heatmap_only_Sample = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{contentfile}{
%%     ~~Describe \code{contentfile} here~~
}
  \item{colname1}{
%%     ~~Describe \code{colname1} here~~
}
  \item{suffixes2}{
%%     ~~Describe \code{suffixes2} here~~
}
  \item{heatmap_range}{
%%     ~~Describe \code{heatmap_range} here~~
}
  \item{filepathin}{
%%     ~~Describe \code{filepathin} here~~
}
  \item{filepathout}{
%%     ~~Describe \code{filepathout} here~~
}
  \item{heatmap_only_Sample}{
%%     ~~Describe \code{heatmap_only_Sample} here~~
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
function (datafile = c("CS1_VS_M1_p005fc2.txt", "CS1_VS_M1_p001fc2.txt"), 
    contentfile = "list.tsv", colname1 = "ProbeName", suffixes2 = "_NS", 
    heatmap_range = 1, filepathin = "./", filepathout = "./", 
    heatmap_only_Sample = FALSE) 
{
    group_tb = ""
    if (class(contentfile) == "data.frame") {
        group_tb = contentfile
    }
    else {
        group_tb <- read.delim(file.path(filepathin, contentfile), 
            sep = "\t", header = T, fill = T, check.names = FALSE)
    }
    for (f_name in datafile) {
        tb1 <- read.table(paste(filepathin, f_name, sep = ""), 
            sep = "\t", quote = "\"", header = T, fill = T, stringsAsFactors = F)
        if (is.na(tb1[1, 1] == TRUE)) {
            print(paste(f_name, " is actually empty!!"))
            next
        }
        rownames(tb1) <- tb1[, colname1]
        tb2 <- tb1[, grep(suffixes2, colnames(tb1))]
        colnames(tb2) <- gsub(suffixes2, "", colnames(tb2))
        tb2_heatmap_tmp <- apply(tb2, 1, function(x) {
            x - median(x)
        })
        tb2_heatmap <- t(tb2_heatmap_tmp)
        if (nrow(tb2_heatmap) < 2 | ncol(tb2_heatmap) < 2) {
            print(paste(f_name, " has less than 2 columns or 2 rows! Heatmap cannot be plotted!"))
            next
        }
        if (nrow(tb2_heatmap) > 60) {
            pdf(file = paste(filepathout, f_name, "_heatmap_plot.pdf", 
                sep = ""))
            if (heatmap_only_Sample == FALSE) {
                heatmap.2(as.matrix(tb2_heatmap), col = greenred, 
                  scale = "none", breaks = c(seq(-heatmap_range, 
                    heatmap_range, 0.1)), trace = c("none"), 
                  cexRow = 0.7, dendrogram = "both", main = gsub(".txt", 
                    "", f_name), labRow = rep("", nrow(tb2_heatmap)), 
                  cexCol = min(1, 8/max(nchar(colnames(tb2_heatmap)))), 
                  margins = c(6, 2))
            }
            if (heatmap_only_Sample == TRUE) {
                heatmap.2(as.matrix(tb2_heatmap), col = greenred, 
                  scale = "none", breaks = c(seq(-heatmap_range, 
                    heatmap_range, 0.1)), trace = c("none"), 
                  cexRow = 0.7, dendrogram = "column", main = gsub(".txt", 
                    "", f_name), labRow = rep("", nrow(tb2_heatmap)), 
                  cexCol = min(1, 8/max(nchar(colnames(tb2_heatmap)))), 
                  margins = c(6, 2))
            }
            dev.off()
            png(file = paste(filepathout, f_name, "_heatmap_plot.png", 
                sep = ""), height = 14 * 8/6, width = 14, units = "cm", 
                res = 600)
            if (heatmap_only_Sample == FALSE) {
                heatmap.2(as.matrix(tb2_heatmap), col = greenred, 
                  scale = "none", breaks = c(seq(-heatmap_range, 
                    heatmap_range, 0.1)), trace = c("none"), 
                  cexRow = 0.7, dendrogram = "both", main = gsub(".txt", 
                    "", f_name), labRow = rep("", nrow(tb2_heatmap)), 
                  cexCol = min(1, 8/max(nchar(colnames(tb2_heatmap)))), 
                  margins = c(6, 2))
            }
            if (heatmap_only_Sample == TRUE) {
                heatmap.2(as.matrix(tb2_heatmap), col = greenred, 
                  scale = "none", breaks = c(seq(-heatmap_range, 
                    heatmap_range, 0.1)), trace = c("none"), 
                  cexRow = 0.7, dendrogram = "column", main = gsub(".txt", 
                    "", f_name), labRow = rep("", nrow(tb2_heatmap)), 
                  cexCol = min(1, 8/max(nchar(colnames(tb2_heatmap)))), 
                  margins = c(6, 2))
            }
            dev.off()
        }
        if (nrow(tb2_heatmap) <= 60) {
            pdf(file = paste(filepathout, f_name, "_heatmap_plot.pdf", 
                sep = ""))
            if (heatmap_only_Sample == FALSE) {
                heatmap.2(as.matrix(tb2_heatmap), col = greenred, 
                  scale = "none", breaks = c(seq(-heatmap_range, 
                    heatmap_range, 0.1)), trace = c("none"), 
                  cexRow = 0.7, dendrogram = "both", main = gsub(".txt", 
                    "", f_name), cexCol = min(1, 8/max(nchar(colnames(tb2_heatmap)))), 
                  margins = c(6, 8))
            }
            if (heatmap_only_Sample == TRUE) {
                heatmap.2(as.matrix(tb2_heatmap), col = greenred, 
                  scale = "none", breaks = c(seq(-heatmap_range, 
                    heatmap_range, 0.1)), trace = c("none"), 
                  cexRow = 0.7, dendrogram = "column", main = gsub(".txt", 
                    "", f_name), cexCol = min(1, 8/max(nchar(colnames(tb2_heatmap)))), 
                  margins = c(6, 8))
            }
            dev.off()
            png(file = paste(filepathout, f_name, "_heatmap_plot.png", 
                sep = ""), height = 14 * 8/6, width = 14, units = "cm", 
                res = 600)
            if (heatmap_only_Sample == FALSE) {
                heatmap.2(as.matrix(tb2_heatmap), col = greenred, 
                  scale = "none", breaks = c(seq(-heatmap_range, 
                    heatmap_range, 0.1)), trace = c("none"), 
                  cexRow = 0.7, dendrogram = "both", main = gsub(".txt", 
                    "", f_name), cexCol = min(1, 8/max(nchar(colnames(tb2_heatmap)))), 
                  margins = c(6, 8))
            }
            if (heatmap_only_Sample == TRUE) {
                heatmap.2(as.matrix(tb2_heatmap), col = greenred, 
                  scale = "none", breaks = c(seq(-heatmap_range, 
                    heatmap_range, 0.1)), trace = c("none"), 
                  cexRow = 0.7, dendrogram = "column", main = gsub(".txt", 
                    "", f_name), cexCol = min(1, 8/max(nchar(colnames(tb2_heatmap)))), 
                  margins = c(6, 8))
            }
            dev.off()
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
