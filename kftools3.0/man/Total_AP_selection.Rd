\name{Total_AP_selection}
\alias{Total_AP_selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Total_AP_selection
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Total_AP_selection(datafile = "Book1.txt", contentfile = "list.tsv", suffixes2 = "_flags", filepathin = "", filepathout = "./")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{contentfile}{
%%     ~~Describe \code{contentfile} here~~
}
  \item{suffixes2}{
%%     ~~Describe \code{suffixes2} here~~
}
  \item{filepathin}{
%%     ~~Describe \code{filepathin} here~~
}
  \item{filepathout}{
%%     ~~Describe \code{filepathout} here~~
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
function (datafile = "Book1.txt", contentfile = "list.tsv", suffixes2 = "_flags", 
    filepathin = "", filepathout = "./") 
{
    group_tb = ""
    if (class(contentfile) == "data.frame") {
        group_tb = contentfile
    }
    else {
        group_tb <- read.delim(contentfile, sep = "\t", header = T, 
            fill = T, check.names = FALSE)
    }
    outname = ""
    if (class(datafile) == "character") {
        tb1 <- read.delim(datafile, sep = "\t", quote = "\"", 
            header = T, fill = T)
        outname = datafile
    }
    else {
        tb1 = datafile
        rm(datafile)
        outname = "Book1.txt"
    }
    flags_tb <- tb1[, grep(suffixes2, colnames(tb1))]
    colnames(flags_tb) <- gsub(suffixes2, "", colnames(flags_tb))
    total_group <- as.character(unique(group_tb[, "Group"]))
    all_flags <- rep(FALSE, nrow(flags_tb))
    for (group_i in total_group) {
        sample_i <- as.character(group_tb[grep(group_i, group_tb[, 
            "Group"], fixed = T), "name"])
        print(sample_i)
        print(group_i)
        tmp_tb <- as.data.frame(flags_tb[, match(sample_i, colnames(flags_tb))])
        tmp_flags <- apply(tmp_tb, 1, function(x) {
            length(grep("A", x)) == 0
        })
        all_flags <- all_flags | tmp_flags
    }
    AP_selection_tb <- tb1[all_flags, ]
    write.table(AP_selection_tb, file = paste(basename(outname), 
        "_AP_selection.xls", sep = ""), sep = "\t", quote = FALSE, 
        col.names = TRUE, row.names = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
