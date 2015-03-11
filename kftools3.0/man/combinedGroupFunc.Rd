\name{combinedGroupFunc}
\alias{combinedGroupFunc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
combinedGroupFunc
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
combinedGroupFunc(datafile = "Book1.txt", group = "group.txt", combinedGroup = TRUE, func = functionA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{combinedGroup}{
%%     ~~Describe \code{combinedGroup} here~~
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
function (datafile = "Book1.txt", group = "group.txt", combinedGroup = TRUE, 
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
        group = read.delim("group.txt")
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
    if (missing(combinedGroup)) {
        combinedGroup = TRUE
    }
    else if (class(combinedGroup) == "logical") {
        combinedGroup = combinedGroup
    }
    else {
        stop("combinedGroup should be logical!")
    }
    list_tb = ""
    if (!missing(combinedGroup) && combinedGroup == TRUE) {
        list_tb = read.delim("list.tsv")
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
    picName = ""
    print("looping the Comparisions!!!!")
    for (i in 1:dim(group)[1]) {
        group1 = unlist(strsplit(as.character(group[i, 1]), ",|;"))
        group2 = unlist(strsplit(as.character(group[i, 2]), ",|;"))
        combination = data.frame(Condition1 = paste(group1, collapse = "_"), 
            Condition2 = paste(group2, collapse = "_"), Selected = "Y")
        list = ""
        list_tbtemp = ""
        if (combinedGroup) {
            list_tbtemp = list_tb[list_tb$Group \%in\% c(group1, 
                group2), ]
            list_tbtemp$Group = as.character(list_tbtemp$Group)
            for (k in 1:2) {
                groupTemp = eval(parse(text = paste("group", 
                  k, sep = "")))
                altername = paste(groupTemp, sep = "", collapse = "_")
                list_tbtemp[list_tbtemp$Group \%in\% groupTemp, 
                  "Group"] = altername
            }
            list = list_tbtemp
        }
        else {
            list = data.frame(Samples = c(group1, group2), name = c(group1, 
                group2), Group = c(rep(paste(group1, collapse = "_"), 
                length(group1)), rep(paste(group2, collapse = "_"), 
                length(group2))))
        }
        picName = paste(paste(group1, collapse = "_"), "_VS_", 
            paste(group2, collapse = "_"), sep = "")
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
