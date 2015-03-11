\name{T_testHTA}
\alias{T_testHTA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
T_testHTA
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
T_testHTA(data1, data2, T_type = 1, islog = 1, suff = "_flags", isadjustP = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
%%     ~~Describe \code{data1} here~~
}
  \item{data2}{
%%     ~~Describe \code{data2} here~~
}
  \item{T_type}{
%%     ~~Describe \code{T_type} here~~
}
  \item{islog}{
%%     ~~Describe \code{islog} here~~
}
  \item{suff}{
%%     ~~Describe \code{suff} here~~
}
  \item{isadjustP}{
%%     ~~Describe \code{isadjustP} here~~
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
function (data1, data2, T_type = 1, islog = 1, suff = "_flags", 
    isadjustP = 0) 
{
    if (nrow(data1) != nrow(data2)) {
        stop("The dim of two groups are not the same!\n")
    }
    t_statistic <- numeric(nrow(data1))
    pvalues <- numeric(nrow(data1))
    mean1 <- numeric(nrow(data1))
    mean2 <- numeric(nrow(data1))
    median1 <- numeric(nrow(data1))
    median2 <- numeric(nrow(data1))
    stdvar1 <- numeric(nrow(data1))
    stdvar2 <- numeric(nrow(data1))
    foldchange <- numeric(nrow(data1))
    for (i in 1:nrow(data1)) {
        mean1[i] <- mean(as.numeric(data1[i, ]))
        mean2[i] <- mean(as.numeric(data2[i, ]))
        median1[i] <- median(as.numeric(data1[i, ]))
        median2[i] <- median(as.numeric(data2[i, ]))
        stdvar1[i] <- sd(as.numeric(data1[i, ]))
        stdvar2[i] <- sd(as.numeric(data2[i, ]))
        if (stdvar1[i] == 0 & stdvar2[i] == 0) {
            t_statistic[i] <- NA
            pvalues[i] <- 1
        }
        else {
            if (T_type == 0) {
                t.result <- t.test(data1[i, ], data2[i, ], var.equal = T)
            }
            if (T_type == 1) {
                t.result <- t.test(data1[i, ], data2[i, ], var.equal = F)
            }
            if (T_type == 2) {
                t.result <- t.test(as.numeric(data1[i, ]), as.numeric(data2[i, 
                  ]), var.equal = F, paired = T)
            }
            t_statistic[i] <- t.result$statistic
            pvalues[i] <- t.result$p.value
        }
        if (islog == 1) {
            x <- 2^data1[i, ]
            y <- 2^data2[i, ]
            foldchange[i] <- mean(x)/mean(y)
        }
        else {
            foldchange[i] <- mean2[i]/mean1[i]
        }
    }
    if (isadjustP == 1) {
        fdr <- p.adjust(pvalues, "BH")
        result <- cbind(pvalues, foldchange, fdr, data1, data2)
    }
    else {
        result <- cbind(pvalues, foldchange, data1, data2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
