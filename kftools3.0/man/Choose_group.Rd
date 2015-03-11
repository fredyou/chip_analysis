\name{Choose_group}
\alias{Choose_group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Choose_group
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Choose_group(contentfile = "content.txt", filepath = "", colname = "GroupID", fileoutname = "Condition.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contentfile}{
%%     ~~Describe \code{contentfile} here~~
}
  \item{filepath}{
%%     ~~Describe \code{filepath} here~~
}
  \item{colname}{
%%     ~~Describe \code{colname} here~~
}
  \item{fileoutname}{
%%     ~~Describe \code{fileoutname} here~~
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
function (contentfile = "content.txt", filepath = "", colname = "GroupID", 
    fileoutname = "Condition.txt") 
{
    d1 <- read.table(paste(filepath, contentfile, sep = ""), 
        header = TRUE)
    Group <- levels(d1[[colname]])
    maxiter <- length(Group) - 1
    for (i in 1:maxiter) {
        temp <- matrix("", length(Group) - i, 3)
        for (j in 1:nrow(temp)) {
            temp[j, 1] <- Group[i]
            temp[j, 2] <- Group[i + j]
            temp[j, 3] <- "Y"
        }
        if (i == 1) {
            Condition <- temp
        }
        else {
            Condition <- rbind(Condition, temp)
        }
    }
    Condition.names <- c("Condition1", "Condition2", "Selected")
    write.table(Condition, file = paste(filepath, fileoutname, 
        sep = ""), sep = "\t", quote = FALSE, col.names = Condition.names, 
        row.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
