\name{cor_plot}
\alias{cor_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cor_plot
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cor_plot(datafile = "CS1_VS_M1_p001fc2.txt", contentfile = "list.tsv", colname1 = "ProbeName", suffixes2 = "_NS", filepathin = "./", filepathout = "./", picName = "g1_vs_g2")
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
function (datafile = "CS1_VS_M1_p001fc2.txt", contentfile = "list.tsv", 
    colname1 = "ProbeName", suffixes2 = "_NS", filepathin = "./", 
    filepathout = "./", picName = "g1_vs_g2") 
{
    if (class(datafile) == "character") {
        tb1 <- read.table(paste(filepathin, datafile, sep = ""), 
            sep = "\t", quote = "\"", header = T, fill = T, stringsAsFactors = F)
        f_name = datafile
    }
    else if (class(datafile) == "data.frame") {
        tb1 = datafile
        f_name = "Book1_AP_selection"
    }
    else {
        stop("datafile is neither character or dataframe!!")
    }
    rownames(tb1) <- tb1[, colname1]
    tb2 <- tb1[, grep(suffixes2, colnames(tb1))]
    colnames(tb2) <- gsub(suffixes2, "", colnames(tb2))
    if (nrow(tb2) < 2 | ncol(tb2) < 2) {
        print(paste(f_name, " has less than 2 columns or 2 rows! correlation cannot be plotted!"))
    }
    t_cor = cor(tb2)
    if (exists(picName)) {
        f_name = picName
    }
    write.table(t_cor, file = paste(filepathout, f_name, "_correlation_coefficient.xls", 
        sep = ""), sep = "\t", col.names = T, row.names = T, 
        quote = F)
    min <- min(t_cor)
    max <- max(t_cor)
    ColorRamp <- colorRampPalette(c("white", "blue"))(100)
    ColorLevels <- seq(min, max, length = length(ColorRamp))
    pdf(file = paste(filepathout, f_name, "_correlation_plot.pdf", 
        sep = ""))
    layout(matrix(data = c(1, 2), nrow = 1, ncol = 2), widths = c(10, 
        2), heights = c(1, 1))
    par(mar = c(5, 5, 2.5, 2))
    image(1:ncol(t_cor), 1:nrow(t_cor), t_cor, col = ColorRamp, 
        xlab = "", ylab = "", axes = FALSE)
    axis(BELOW <- 1, at = 1:ncol(t_cor), labels = colnames(t_cor), 
        las = HORIZONTAL <- 2, cex.axis = 0.5)
    axis(LEFT <- 2, at = 1:nrow(t_cor), labels = rownames(t_cor), 
        las = HORIZONTAL <- 1, cex.axis = 0.5)
    par(mar = c(5, 2, 2.5, 2))
    image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
        nrow = 1), col = ColorRamp, xlab = "", ylab = "", xaxt = "n")
    dev.off()
    png(file = paste(filepathout, f_name, "_correlation_plot.png", 
        sep = ""), height = 14 * 7/6, width = 14, units = "cm", 
        res = 600)
    layout(matrix(data = c(1, 2), nrow = 1, ncol = 2), widths = c(10, 
        2), heights = c(1, 1))
    par(mar = c(5, 5, 2.5, 2))
    image(1:ncol(t_cor), 1:nrow(t_cor), t_cor, col = ColorRamp, 
        xlab = "", ylab = "", axes = FALSE)
    axis(BELOW <- 1, at = 1:ncol(t_cor), labels = colnames(t_cor), 
        las = HORIZONTAL <- 2, cex.axis = 0.5)
    axis(LEFT <- 2, at = 1:nrow(t_cor), labels = rownames(t_cor), 
        las = HORIZONTAL <- 1, cex.axis = 0.5)
    par(mar = c(5, 2, 2.5, 2))
    image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), 
        nrow = 1), col = ColorRamp, xlab = "", ylab = "", xaxt = "n")
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
