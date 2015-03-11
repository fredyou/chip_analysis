\name{intersect_Ttestfiles}
\alias{intersect_Ttestfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
intersect_Ttestfiles
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intersect_Ttestfiles(files = files)
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
    temp = read.delim(files[1], nrow = 1)
    probe = colnames(temp)[grep("probe|systematic", colnames(temp), 
        ignore.case = T)][1]
    flag = colnames(temp)[grep("_call|_flags", colnames(temp))[1]]
    flags = regmatches(flag, regexpr("_[^_]*$", flag))
    for (i in 1:length(files)) {
        print(files[i])
        dataname = paste("data", i, sep = "_")
        print(dataname)
        input = read.delim(file = files[i], header = T)
        assign(dataname, input)
    }
    loop = 1
    for (i in 1:(length(files) - 1)) {
        if (loop == 1) {
            data1n = get(paste("data_", i, sep = ""))
            data2n = get(paste("data_", i + 1, sep = ""))
            overlap = intersect(unique(data1n[, probe]), unique(data2n[, 
                probe]))
            overlap = as.data.frame(overlap)
            rownames(overlap) = overlap[, 1]
            data1n_NS = grep("_NS", colnames(data1n))
            end1 = rev(data1n_NS)[1]
            data1n_NS = colnames(data1n[, grep("_NS", colnames(data1n))])
            data2n_NS = colnames(data2n[, grep("_NS", colnames(data2n))])
            boolean_2in1 = !(data2n_NS \%in\% data1n_NS)
            data2n_NS = data2n_NS[boolean_2in1]
            colnames(data2n)[2] = paste(gsub("p005.*.txt", "", 
                files[i + 1]), colnames(get(paste("data_", i + 
                1, sep = "")))[2], sep = "", collapse = "")
            colnames(data2n)[3] = paste(gsub("p005*.txt", "", 
                files[i + 1]), colnames(get(paste("data_", i + 
                1, sep = "")))[3], sep = "", collapse = "")
            colnames(data1n)[2] = paste(gsub("p005.*.txt", "", 
                files[i]), colnames(get(paste("data_", i, sep = "")))[2], 
                sep = "", collapse = "")
            colnames(data1n)[3] = paste(gsub("p005.*.txt", "", 
                files[i]), colnames(get(paste("data_", i, sep = "")))[3], 
                sep = "", collapse = "")
            data1n = data1n[data1n[, probe] \%in\% rownames(overlap), 
                ]
            data2n = data2n[data2n[, probe] \%in\% rownames(overlap), 
                ]
            newdata = cbind(data1n[, 1:end1], data2n[, c(colnames(data2n)[2:3], 
                data2n_NS)], data1n[, (end1 + 1):(ncol(data1n))])
        }
        if (loop == 2) {
            data3n = get(paste("data_", i + 1, sep = ""))
            overlap = intersect(unique(newdata[, probe]), unique(data3n[, 
                probe]))
            overlap = as.data.frame(overlap)
            rownames(overlap) = overlap[, 1]
            newdata = newdata[newdata[, probe] \%in\% rownames(overlap), 
                ]
            newdata_NS = grep("_NS", colnames(newdata))
            end2 = rev(newdata_NS)[1]
            newdata_NS = colnames(newdata[, grep("_NS", colnames(newdata))])
            data3n_NS = colnames(data3n[, grep("_NS", colnames(data3n))])
            boolean_3innew = !(data3n_NS \%in\% newdata_NS)
            data3n_NS = data3n_NS[boolean_3innew]
            colnames(data3n)[2] = paste(gsub("p005.*.txt", "", 
                files[i + 1]), colnames(data3n)[2], sep = "", 
                collapse = "")
            colnames(data3n)[3] = paste(gsub("p005.*.txt", "", 
                files[i + 1]), colnames(data3n)[3], sep = "", 
                collapse = "")
            data3n = data3n[data3n[, probe] \%in\% rownames(overlap), 
                ]
            if (length(data3n_NS) != 0) {
                newdata = cbind(newdata[, 1:end2], data3n[, c(colnames(data3n)[2:3], 
                  data3n_NS)], newdata[, (end2 + 1):ncol(newdata)])
            }
            else {
                newdata = cbind(newdata[, 1:end2], data3n[, colnames(data3n)[2:3]], 
                  newdata[, (end2 + 1):ncol(newdata)])
            }
        }
        loop = 2
    }
    flags = grep(flags, colnames(newdata))
    newdata = newdata[, -flags]
    write.table(newdata, file = "Overlap.txt", sep = "\t", quote = FALSE, 
        col.names = TRUE, row.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
