\name{BatchFoldchange}
\alias{BatchFoldchange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BatchFoldchange
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BatchFoldchange(DataFile = "data.txt", ContentFile = "Content.txt", BatchInfoFile = "Condition.txt", OutfileSuffixes = "_fc2", FilePath = getwd(), colname1 = "ProbeName", colname2 = "filename", prefixes = "", suffixes = "_NS", suffixes2 = "_flags", Paired = FALSE, is.log = 1, is.FC = 1, envelope = 1, FC = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataFile}{
%%     ~~Describe \code{DataFile} here~~
}
  \item{ContentFile}{
%%     ~~Describe \code{ContentFile} here~~
}
  \item{BatchInfoFile}{
%%     ~~Describe \code{BatchInfoFile} here~~
}
  \item{OutfileSuffixes}{
%%     ~~Describe \code{OutfileSuffixes} here~~
}
  \item{FilePath}{
%%     ~~Describe \code{FilePath} here~~
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
  \item{suffixes2}{
%%     ~~Describe \code{suffixes2} here~~
}
  \item{Paired}{
%%     ~~Describe \code{Paired} here~~
}
  \item{is.log}{
%%     ~~Describe \code{is.log} here~~
}
  \item{is.FC}{
%%     ~~Describe \code{is.FC} here~~
}
  \item{envelope}{
%%     ~~Describe \code{envelope} here~~
}
  \item{FC}{
%%     ~~Describe \code{FC} here~~
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
function (DataFile = "data.txt", ContentFile = "Content.txt", 
    BatchInfoFile = "Condition.txt", OutfileSuffixes = "_fc2", 
    FilePath = getwd(), colname1 = "ProbeName", colname2 = "filename", 
    prefixes = "", suffixes = "_NS", suffixes2 = "_flags", Paired = FALSE, 
    is.log = 1, is.FC = 1, envelope = 1, FC = 2) 
{
    Conditions = ""
    if (class(BatchInfoFile) == "data.frame") {
        Conditions = BatchInfoFile
    }
    else {
        Conditions <- read.delim(file.path(FilePath, BatchInfoFile), 
            sep = "\t", header = T, fill = T, check.names = FALSE)
    }
    BatchInfo <- subset(Conditions, Selected == "Y")
    print(paste("The number of iterations is ", as.character(nrow(BatchInfo))))
    if (class(DataFile) == "character") {
        tmp = read.delim(file.path(FilePath, DataFile), nrow = 1)
    }
    else if (class(DataFile) == "data.frame") {
        tmp = DataFile[1, ]
    }
    coltmp = colnames(tmp)[grep("_NS", colnames(tmp))]
    listtmp = ""
    if (class(ContentFile) == "data.frame") {
        listtmp = ContentFile
    }
    else {
        listtmp = read.delim(ContentFile)
    }
    collist = paste(listtmp$name, "_NS", sep = "")
    if (!identical(sort(intersect(coltmp, collist)), sort(collist))) {
        cat("coltmp: ", sort(coltmp), "\n")
        cat("collist: ", sort(collist), "\n")
        stop("Colnames not matching listnames!\n")
    }
    if (!exists("dataG")) {
        print("dataG not exist!")
        if (class(DataFile) == "character") {
            data <- read.delim(file.path(FilePath, DataFile), 
                sep = "\t", quote = "", na.strings = "", header = T, 
                fill = T, check.names = FALSE)
        }
        else if (class(DataFile) == "data.frame") {
            data = DataFile
        }
        assign("dataG", data, envir = .GlobalEnv)
    }
    else {
        print("dataG exsits!")
        coltmp = colnames(dataG)[grep("_NS", colnames(dataG))]
        if (identical(sort(intersect(coltmp, collist)), sort(collist))) {
            data = dataG
        }
        else {
            if (class(DataFile) == "character") {
                data <- read.delim(file.path(FilePath, DataFile), 
                  sep = "\t", quote = "", na.strings = "", header = T, 
                  fill = T, check.names = FALSE)
            }
            else if (class(DataFile) == "data.frame") {
                data = DataFile
            }
            assign("dataG", data, envir = .GlobalEnv)
        }
    }
    print(dim(data))
    probeset <- as.vector(data[[colname1]])
    data_annot <- data[, -grep(suffixes, colnames(data))]
    data_annot <- data_annot[, -grep(suffixes2, colnames(data_annot))]
    if (ncol(as.matrix(data_annot)) == 1) {
        is.annot = "N"
    }
    else {
        is.annot = "Y"
        data_annot <- as.matrix(data_annot[, -1])
    }
    for (i in 1:nrow(BatchInfo)) {
        if (Paired) {
            cat("UnPaired: ", Paired, "\n")
            label1 <- paste(prefixes, BatchInfo[i, 1], suffixes, 
                sep = "")
            label2 <- paste(prefixes, BatchInfo[i, 2], suffixes, 
                sep = "")
            label1Flag <- paste(prefixes, BatchInfo[i, 1], suffixes2, 
                sep = "")
            label2Flag <- paste(prefixes, BatchInfo[i, 2], suffixes2, 
                sep = "")
            print(paste("the NO is", as.character(i)))
            print(as.character(BatchInfo[i, 1]))
            print(as.character(BatchInfo[i, 2]))
        }
        else {
            print("Paired")
            SampleInfo = ""
            if (class(ContentFile) == "data.frame") {
                SampleInfo = ContentFile
            }
            else {
                SampleInfo <- read.delim(file.path(FilePath, 
                  ContentFile), sep = "\t", header = T, fill = T, 
                  check.names = FALSE)
            }
            label1Raw <- subset(SampleInfo, Group == as.character(BatchInfo[i, 
                1]), select = colname2)
            label1Raw <- as.matrix(label1Raw)
            label1 <- paste(prefixes, label1Raw, suffixes, sep = "")
            label2Raw <- subset(SampleInfo, Group == as.character(BatchInfo[i, 
                2]), select = colname2)
            label2Raw <- as.matrix(label2Raw)
            label2 <- paste(prefixes, label2Raw, suffixes, sep = "")
            label1Flag <- paste(prefixes, label1Raw, suffixes2, 
                sep = "")
            label2Flag <- paste(prefixes, label2Raw, suffixes2, 
                sep = "")
            print(paste("the NO is", as.character(i)))
            print(as.character(BatchInfo[i, 1]))
            print(label1Raw)
            print(as.character(BatchInfo[i, 2]))
            print(label2Raw)
        }
        cat("label1: ", label1, " label1Flag: ", label1Flag, 
            "\n")
        cat("label2: ", label2, " label2Flag: ", label2Flag, 
            "\n")
        data1 <- subset(data, select = c(label1, label1Flag))
        data2 <- subset(data, select = c(label2, label2Flag))
        if (is.log == 1) {
            mean1 = apply(2^data1[label1], 1, mean, na.rm = T)
            mean2 = apply(2^data2[label2], 1, mean, na.rm = T)
        }
        else {
            mean1 = apply(data1[label1], 1, mean, na.rm = T)
            mean2 = apply(data2[label2], 1, mean, na.rm = T)
        }
        foldchange = mean1/mean2
        if (is.annot == "N") {
            Groupdata <- cbind(probeset, data1, data2, foldchange)
        }
        else {
            Groupdata <- cbind(probeset, data1, data2, foldchange, 
                data_annot)
        }
        data1Flag <- subset(Groupdata, select = label1Flag)
        data2Flag <- subset(Groupdata, select = label2Flag)
        Flag1 <- stat_Flag(data1Flag)
        Flag2 <- stat_Flag(data2Flag)
        result <- cbind(Groupdata, Flag1, Flag2)
        if (is.FC == 1) {
            Filter <- Filter_Flag(result, Propotion = envelope)
            Filter <- subset(Filter, (foldchange >= FC | foldchange <= 
                1/FC))
            Filter <- subset(Filter, select = -c(Flag1, Flag2))
        }
        else {
            Filter <- subset(result, select = -c(Flag1, Flag2))
        }
        colnames(Filter)[1] = colname1
        outputdir <- paste(BatchInfo[i, 1], "_VS_", BatchInfo[i, 
            2], sep = "")
        outputdir <- file.path(FilePath, outputdir)
        dir.create(outputdir, showWarnings = FALSE)
        OutName = paste(BatchInfo[i, 1], "_VS_", BatchInfo[i, 
            2], OutfileSuffixes, ".txt", sep = "")
        write.table(Filter, file = file.path(outputdir, OutName), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
