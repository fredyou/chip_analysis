\name{labelCol}
\alias{labelCol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
labelCol
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
labelCol(x, sam_cols, colorCodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{sam_cols}{
%%     ~~Describe \code{sam_cols} here~~
}
  \item{colorCodes}{
%%     ~~Describe \code{colorCodes} here~~
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
function (x, sam_cols, colorCodes) 
{
    if (is.leaf(x)) {
        label <- attr(x, "label")
        code <- substr(label, 1, 1)
        attr(x, "nodePar") <- list(lab.col = colorCodes[sam_cols[match(label, 
            names(sam_cols))]])
        attr(x, "edgePar") <- list(col = colorCodes[sam_cols[match(label, 
            names(sam_cols))]])
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
