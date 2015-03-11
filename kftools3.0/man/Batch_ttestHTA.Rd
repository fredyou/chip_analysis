\name{Batch_ttestHTA}
\alias{Batch_ttestHTA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch_ttestHTA
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Batch_ttestHTA(datafile = "Book1.txt", contentfile = "list.tsv", Batch_infofile = "Combination.txt", paired_ttest_fname = "paired.txt", filepathin = "", filepathout = "", colname1 = "ProbeName", colname2 = "filename", prefixes = "", suffixes = "_NS", T.type = 2, is.log = 1, cutoff = 1, is.adjustP = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{contentfile}{
%%     ~~Describe \code{contentfile} here~~
}
  \item{Batch_infofile}{
%%     ~~Describe \code{Batch_infofile} here~~
}
  \item{paired_ttest_fname}{
%%     ~~Describe \code{paired_ttest_fname} here~~
}
  \item{filepathin}{
%%     ~~Describe \code{filepathin} here~~
}
  \item{filepathout}{
%%     ~~Describe \code{filepathout} here~~
}
  \item{colname1}{
%%     ~~Describe \code{colname1} here~~
}
  \item{colname2}{
%%     ~~Describe \code{colname2} here~~
}
  \item{prefixes}{
%%     ~~Describe \code{prefixes} here~~
}
  \item{suffixes}{
%%     ~~Describe \code{suffixes} here~~
}
  \item{T.type}{
%%     ~~Describe \code{T.type} here~~
}
  \item{is.log}{
%%     ~~Describe \code{is.log} here~~
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{is.adjustP}{
%%     ~~Describe \code{is.adjustP} here~~
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
function (datafile = "Book1.txt", contentfile = "list.tsv", Batch_infofile = "Combination.txt", 
    paired_ttest_fname = "paired.txt", filepathin = "", filepathout = "", 
    colname1 = "ProbeName", colname2 = "filename", prefixes = "", 
    suffixes = "_NS", T.type = 2, is.log = 1, cutoff = 1, is.adjustP = 0) 
{
    data = ""
    if (class(datafile) == "character") {
        data <- read.delim(paste(filepathin, datafile, sep = ""), 
            sep = "\t", quote = "\"", header = T, fill = T)
    }
    else {
        data = datafile
        rm(datafile)
    }
    probeset <- as.vector(data[[colname1]])
    Conditions = ""
    if (class(Batch_infofile) == "data.frame") {
        Conditions = Batch_infofile
    }
    else {
        Conditions <- read.delim(file.path(filepathin, Batch_infofile), 
            sep = "\t", header = T, fill = T, check.names = FALSE)
    }
    Batch_info <- subset(Conditions, Selected == "Y")
    print(paste("The number of iterations is ", as.character(nrow(Batch_info))))
    sampleinfo = ""
    if (class(contentfile) == "data.frame") {
        sampleinfo = contentfile
    }
    else {
        sampleinfo <- read.delim(file.path(filepathin, contentfile), 
            sep = "\t", header = T, fill = T, check.names = FALSE)
    }
    data_annot <- data[, -grep(suffixes, colnames(data))]
    if (ncol(as.matrix(data_annot)) == 1) {
        is.annot = "N"
    }
    else {
        is.annot = "Y"
        data_annot <- as.matrix(data_annot[, -1])
    }
    for (i in 1:nrow(Batch_info)) {
        label1_raw <- subset(sampleinfo, Group == as.character(Batch_info[i, 
            1]), select = colname2)
        label1_raw <- as.matrix(label1_raw)
        label1 <- paste(prefixes, label1_raw, suffixes, sep = "")
        label2_raw <- subset(sampleinfo, Group == as.character(Batch_info[i, 
            2]), select = colname2)
        label2_raw <- as.matrix(label2_raw)
        label2 <- paste(prefixes, label2_raw, suffixes, sep = "")
        print(paste("the NO is", as.character(i)))
        print(as.character(Batch_info[i, 1]))
        print(label1_raw)
        print(as.character(Batch_info[i, 2]))
        print(label2_raw)
        data1 <- subset(data, select = c(label1))
        data2 <- subset(data, select = c(label2))
        result <- T_testHTA(data1, data2, T_type = T.type, islog = is.log, 
            isadjustP = is.adjustP)
        if (is.annot == "N") 
            result <- cbind(probeset, result)
        else result <- cbind(probeset, result, data_annot)
        result <- subset(result, pvalues <= cutoff)
        colnames(result)[1] = colname1
        outputdir <- filepathout
        assign(paste(Batch_info[i, 1], "_VS_", Batch_info[i, 
            2], sep = ""), result, envir = .GlobalEnv)
        write.table(result, file = paste(outputdir, "/", Batch_info[i, 
            1], "_VS_", Batch_info[i, 2], ".txt", sep = ""), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
