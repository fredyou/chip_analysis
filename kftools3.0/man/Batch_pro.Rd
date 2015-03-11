\name{Batch_pro}
\alias{Batch_pro}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch_pro
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Batch_pro(contentfile, filepathin, filepathout, Num, headline)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contentfile}{
%%     ~~Describe \code{contentfile} here~~
}
  \item{filepathin}{
%%     ~~Describe \code{filepathin} here~~
}
  \item{filepathout}{
%%     ~~Describe \code{filepathout} here~~
}
  \item{Num}{
%%     ~~Describe \code{Num} here~~
}
  \item{headline}{
%%     ~~Describe \code{headline} here~~
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
function (contentfile, filepathin, filepathout, Num, headline) 
{
    content <- read.table(contentfile, header = T)
    for (i in 1:Num) {
        d1 <- readLines(paste(filepathin, content[[2]][i], ".txt", 
            sep = ""))
        d2 = d1[headline:length(d1)]
        writeLines(d2, paste(filepathout, content[[2]][i], ".txt", 
            sep = ""), sep = "\n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
