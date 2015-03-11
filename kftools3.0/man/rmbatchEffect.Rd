\name{rmbatchEffect}
\alias{rmbatchEffect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
rmbatchEffect
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rmbatchEffect(datafile = "Book1.txt", phenodata = "pheno.tsv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{phenodata}{
%%     ~~Describe \code{phenodata} here~~
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
function (datafile = "Book1.txt", phenodata = "pheno.tsv") 
{
    pheno1 = read.delim(phenodata, stringsAsFactors = FALSE, 
        quote = "")
    row.names(pheno1) = pheno1[, 1]
    pheno1 = pheno1[, -1]
    batch1 = as.numeric(pheno1$batch)
    modcombat1 = model.matrix(~1, data = pheno1)
    temp = read.delim(datafile, nrow = 2, stringsAsFactors = FALSE, 
        quote = "")
    colsNotNeeded = grep("_NS|Probe|systematic", colnames(temp), 
        ignore.case = T, invert = T)
    classes = sapply(temp, class)
    classes[colsNotNeeded] = "NULL"
    data = read.delim("Book1.txt", stringsAsFactors = FALSE, 
        quote = "", colClasses = classes, comment.char = "")
    row.names(data) = data[, 1]
    data = as.matrix(data[, -1])
    combat_edata = ComBat(dat = data, batch = batch1, mod = modcombat1, 
        numCovs = NULL, par.prior = TRUE, prior.plots = TRUE)
    combat_edata = as.data.frame(combat_edata)
    combat_edata = data.frame(ProbeName = row.names(combat_edata), 
        combat_edata)
    write.table(combat_edata, file = "nordata_after_rmBatch.txt", 
        sep = "\t", row.names = F)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
