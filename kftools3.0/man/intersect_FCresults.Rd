\name{intersect_FCresults}
\alias{intersect_FCresults}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
intersect_FCresults
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intersect_FCresults(files = files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{
%%     ~~Describe \code{files} here~~
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
function (files = files) 
{
    for (i in 1:length(files)) {
        print(files[i])
        dataname = paste("data", i, sep = "_")
        print(dataname)
        input = read.delim(file = files[i], header = T)
        assign(dataname, input)
    }
    temp = read.delim(files[1], nrow = 1)
    probelabel = colnames(temp)[grep("probe|systematic", colnames(temp), 
        ignore.case = T)][1]
    flag = colnames(temp)[grep("call|flags", colnames(temp))[1]]
    flagsSuffix = regmatches(flag, regexpr("_[^_]*$", flag))
    loop = 1
    for (i in 1:(length(files) - 1)) {
        if (loop == 1) {
            data1n = get(paste("data_", i, sep = ""))
            data2n = get(paste("data_", i + 1, sep = ""))
            overlap = intersect(unique(data1n[, probelabel[1]]), 
                unique(data2n[, probelabel[1]]))
            overlap = as.data.frame(overlap)
            rownames(overlap) = overlap[, 1]
            data1n_NS = grep("_NS", colnames(data1n))
            end1 = data1n_NS[length(data1n_NS)]
            data1n_NS = colnames(data1n[, grep("_NS", colnames(data1n))])
            data2n_NS = colnames(data2n[, grep("_NS", colnames(data2n))])
            boolean_2in1 = !(data2n_NS \%in\% data1n_NS)
            data2n_NS = data2n_NS[boolean_2in1]
            fcIndex = grep("foldchange", colnames(data1n))
            flagdata2n = unlist(lapply(gsub("_NS", flagsSuffix[1], 
                data2n_NS), function(x) {
                grep(x, colnames(data2n))
            }))
            colnames(data2n)[fcIndex] = paste(gsub("fc2.*.txt", 
                "", files[i + 1]), colnames(get(paste("data_", 
                i + 1, sep = "")))[fcIndex], sep = "", collapse = "")
            colnames(data1n)[fcIndex] = paste(gsub("fc2.*.txt", 
                "", files[i]), colnames(get(paste("data_", i, 
                sep = "")))[fcIndex], sep = "", collapse = "")
            data1n = data1n[data1n[, probelabel[1]] \%in\% rownames(overlap), 
                ]
            data2n = data2n[data2n[, probelabel[1]] \%in\% rownames(overlap), 
                ]
            newdata = cbind(data1n[, 1:end1], data2n[, c(colnames(data2n)[fcIndex], 
                data2n_NS, colnames(data2n)[flagdata2n])], data1n[, 
                (end1 + 1):(ncol(data1n))])
        }
        if (loop == 2) {
            data3n = get(paste("data_", i + 1, sep = ""))
            overlap = intersect(unique(newdata[, probelabel[1]]), 
                unique(data3n[, probelabel[1]]))
            overlap = as.data.frame(overlap)
            rownames(overlap) = overlap[, 1]
            newdata = newdata[newdata[, probelabel[1]] \%in\% rownames(overlap), 
                ]
            newdata_NS = grep("_NS", colnames(newdata))
            end2 = newdata_NS[length(newdata_NS)]
            newdata_NS = colnames(newdata[, grep("_NS", colnames(newdata))])
            data3n_NS = colnames(data3n[, grep("_NS", colnames(data3n))])
            boolean_3innew = !(data3n_NS \%in\% newdata_NS)
            data3n_NS = data3n_NS[boolean_3innew]
            fcIndex = grep("foldchange", colnames(data3n))
            flagdata3n = unlist(lapply(gsub("_NS", flagsSuffix[1], 
                data3n_NS), function(x) {
                grep(x, colnames(data3n))
            }))
            colnames(data3n)[fcIndex] = paste(gsub("fc2.txt", 
                "", files[i + 1]), colnames(data3n)[fcIndex], 
                sep = "", collapse = "")
            data3n = data3n[data3n[, probelabel[1]] \%in\% rownames(overlap), 
                ]
            if (length(data3n_NS) != 0) {
                newdata = cbind(newdata[, 1:end2], data3n[, c(colnames(data3n)[fcIndex], 
                  data3n_NS, colnames(data3n)[flagdata3n])], 
                  newdata[, (end2 + 1):ncol(newdata)])
            }
            else {
                newdata = cbind(newdata[, 1:end2], data3n[, colnames(data3n)[fcIndex]], 
                  newdata[, (end2 + 1):ncol(newdata)])
            }
        }
        loop = 2
    }
    write.table(newdata, file = "Overlap.txt", sep = "\t", quote = FALSE, 
        col.names = TRUE, row.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
