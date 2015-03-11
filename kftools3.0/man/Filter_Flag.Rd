\name{Filter_Flag}
\alias{Filter_Flag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter_Flag
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Filter_Flag(data, Propotion = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{Propotion}{
%%     ~~Describe \code{Propotion} here~~
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
function (data, Propotion = 1) 
{
    Is = numeric(nrow(data))
    for (i in 1:nrow(data)) {
        if (data$foldchange[i] >= 1 & data$Flag1[i] >= Propotion) 
            Is[i] = 1
        else if (data$foldchange[i] < 1 & data$Flag2[i] >= Propotion) 
            Is[i] = 1
    }
    data <- cbind(data, Is)
    new <- subset(data, Is == 1, select = -Is)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
