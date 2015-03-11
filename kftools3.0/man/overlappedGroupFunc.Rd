\name{overlappedGroupFunc}
\alias{overlappedGroupFunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
overlappedGroupFunc
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overlappedGroupFunc(datafile = "Book1.txt", group = "groups.txt", BatchInfoFile = "Combination.txt", func = functionA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{BatchInfoFile}{
%%     ~~Describe \code{BatchInfoFile} here~~
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
function (datafile = "Book1.txt", group = "groups.txt", BatchInfoFile = "Combination.txt", 
    func = functionA) 
{
    if (missing(func)) {
        stop("You should assign a function!!")
    }
    else if (class(func) != "character") {
        stop("You`ve assigned the func par, but it should be a string!")
    }
    print("Processing pars!")
    if (missing(group)) {
        group = read.delim("groups.txt", na.strings = "", stringsAsFactors = F, 
            quote = "")
    }
    else if (class(group) == "character") {
        group = read.delim(group)
    }
    else {
        group = group
    }
    group = as.data.frame(t(apply(as.matrix(group), 1, function(x) {
        gsub("\\s+", "", x)
    })))
    if (missing(BatchInfoFile)) {
        BatchInfoFile = read.delim("Combination.txt", stringsAsFactors = F)
    }
    else if (class(BatchInfoFile) == "character") {
        BatchInfoFile = read.delim(BatchInfoFile, stringsAsFactors = F)
    }
    else {
        BatchInfoFile = BatchInfoFile
    }
    if (missing(datafile)) {
        datafile = read.delim("Book1.txt", na.strings = "", quote = "", 
            sep = "\t", stringsAsFactors = F)
    }
    else if (class(datafile) == "character") {
        datafile = read.delim(datafile, na.strings = "", quote = "", 
            sep = "\t", stringsAsFactors = F)
    }
    else {
        datafile = datafile
    }
    picName = "Book1.txt"
    for (i in 1:dim(BatchInfoFile)[1]) {
        con1 = BatchInfoFile[i, 1]
        con2 = BatchInfoFile[i, 2]
        print(paste(con1, con2, sep = " VS "))
        group1 = group[, con1]
        group1 = group1[!is.na(group1)]
        group1 = gsub("(^[0-9])", "s\\1", group1)
        group2 = group[, con2]
        group2 = group2[!is.na(group2)]
        group2 = gsub("(^[0-9])", "s\\1", group2)
        list = data.frame(Samples = c(group1, group2), name = c(group1, 
            group2), Group = c(rep(con1, length(group1)), rep(con2, 
            length(group2))))
        combination = BatchInfoFile[i, ]
        picName = paste(con1, con2, sep = "_")
        if (exists("list") && class(list) == "data.frame") {
            print("Got list !")
        }
        else {
            stop("Something wrong with list generation!!")
        }
        if (exists("combination") && class(combination) == "data.frame") {
            print("Got combination !")
        }
        else {
            stop("Something wrong with combination generation!!")
        }
        eval(parse(text = func))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
