\name{hclust_plot}
\alias{hclust_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
hclust_plot
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hclust_plot(datafile = "Book1.txt", contentfile = "list.tsv", suffixes2 = "_NS", filepathin = "./", filepathout = "./", picName = picName)
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
    filepathin = "./", filepathout = "./", picName = picName) 
{
    group_tb = ""
    if (class(contentfile) == "character") {
        group_tb <- read.table(paste(filepathin, contentfile, 
            sep = ""), sep = "\t", quote = "\"", header = T, 
            fill = T)
    }
    else {
        group_tb = contentfile
        contentfile = paste(unique(contentfile$Group)[1], "_VS_", 
            unique(contentfile$Group)[2])
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
    expr_tb <- data.matrix(tb1[, paste(as.character(group_tb[, 
        2]), "_NS", sep = "")])
    colnames(expr_tb) <- gsub(suffixes2, "", colnames(expr_tb))
    d <- as.dist(dist(t(expr_tb), method = "euclidean"))
    fit <- hclust(d, method = "ward")
    sam_names <- colnames(expr_tb)
    sam_cols <- as.character(group_tb[match(fit$labels, group_tb[, 
        "name"]), "Group"])
    names(sam_cols) <- fit$labels
    groups = as.character(unique(group_tb$Group))
    colorCodes = colorRampPalette(c("Purple", "darkblue", "lightblue", 
        "darkgreen", "lightgreen", "dark orange", "darkred", 
        "pink"))(length(unique(groups)))
    names(colorCodes) = unique(groups)
    d <- dendrapply(as.dendrogram(fit), labelCol, sam_cols, colorCodes)
    if (missing(picName)) {
        picName = datafile
    }
    else {
        picName = picName
    }
    png(file = paste(filepathout, picName, "_hclust_plot.png", 
        sep = ""), height = 14, width = 14, units = "cm", res = 600)
    par(cex = min(1, par("din")[1]/(length(sam_names) * par("csi"))))
    plot(d)
    dev.off()
    pdf(file = paste(filepathout, picName, "_hclust_plot.pdf", 
        sep = ""), height = 14, width = 14)
    par(cex = min(1, par("din")[1]/(length(sam_names) * par("csi"))))
    plot(d)
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
