\name{batch_plotGO_KEGG}
\alias{batch_plotGO_KEGG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
batch_plotGO_KEGG
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
batch_plotGO_KEGG(files = list.files(pattern = ".*VS.*GO.txt"), path = ".", type = "GO")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
%%     ~~Describe \code{files} here~~
}
  \item{path}{
%%     ~~Describe \code{path} here~~
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
function (files = list.files(pattern = ".*VS.*GO.txt"), path = ".", 
    type = "GO") 
{
    if (missing(path)) {
        path = "."
    }
    else {
        path = path
    }
    setwd(path)
    colTerm = ""
    if (!missing("type")) {
        type = type
    }
    else {
        type = "GO"
    }
    stopifnot(grepl("GO*|KEGG*", type, ignore.case = T))
    if (missing(files)) {
        if (type == "GO") {
            files = list.files(pattern = paste(".*VS.*", type, 
                ".txt$", sep = ""), ignore.case = T)
        }
        else {
            files = list.files(pattern = paste(".*VS.*", type, 
                ".*.xls$", sep = ""), ignore.case = T)
        }
    }
    else {
        files = files
    }
    stopifnot(!is.element(TRUE, file.info(files)$isdir))
    for (i in 1:length(files)) {
        sample1 = gsub("([^_]*)_VS.*", "\\1", files[i], perl = T)
        sample2 = gsub(paste("[^_]*_VS_(.*)_(p005fc[0-9]_|fc[0-9]_.*)", 
            type, ".*.(txt|xls)", sep = ""), "\\1", files[i], 
            perl = T, ignore.case = T)
        plot_go_keggEnrich(files[i], sample1, sample2, type = type)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
