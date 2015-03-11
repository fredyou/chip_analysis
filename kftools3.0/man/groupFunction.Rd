\name{groupFunction}
\alias{groupFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
groupFunction
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
groupFunction(file = "Book1.txt", list = "list.tsv", group = "NS", func = "func")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
%%     ~~Describe \code{file} here~~
}
  \item{list}{
%%     ~~Describe \code{list} here~~
}
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{func}{
%%     ~~Describe \code{func} here~~
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
function (file = "Book1.txt", list = "list.tsv", group = "NS", 
    func = "func") 
{
    data = read.delim(file)
    gt = data[, c(1, grep(group, colnames(data)))]
    list = read.delim(list)
    list = list[, 2:3]
    rownames(list) = list$name
    list[, 2] = as.character(list[, 2])
    list[, 1] = as.character(list[, 1])
    require(reshape2)
    gt.m = melt(gt)
    gt.m[, 2] = as.character(gt.m[, 2])
    gt.m[, 2] = gsub("_NS", "", gt.m[, 2])
    gt.m[, 2] = list[gt.m[, 2], 2]
    method = ""
    if (!missing(func)) {
        method = func
    }
    else {
        method = function(x) {
            log2(mean(2^x))
        }
    }
    if (class(func) != "function") {
        stop("func parameter should be set to a function!!")
    }
    ag.gt = as.data.frame(aggregate(gt.m[, 3], list(probe = gt.m[, 
        1], sample = gt.m$variable), method))
    ag.gt = as.data.frame(acast(ag.gt, probe ~ sample, mean, 
        margins = TRUE))
    ag.gt = ag.gt[-dim(ag.gt)[1], -grep("(all)", colnames(ag.gt))]
    ag.gt$Probes = rownames(ag.gt)
    ag.gt = ag.gt[, sort(colnames(ag.gt))]
    write.table(ag.gt, file = paste(group, ".txt", sep = ""), 
        sep = "\t", row.names = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
