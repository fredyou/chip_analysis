\name{batch_diff_ALL}
\alias{batch_diff_ALL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
batch_diff_ALL
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
batch_diff_ALL(datafile = datafile, list = list, combination = combination, chipType = "ge", picName = "picName", pairedFC, setScriptPath = "path")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{list}{
%%     ~~Describe \code{list} here~~
}
  \item{combination}{
%%     ~~Describe \code{combination} here~~
}
  \item{chipType}{
%%     ~~Describe \code{chipType} here~~
}
  \item{picName}{
%%     ~~Describe \code{picName} here~~
}
  \item{pairedFC}{
%%     ~~Describe \code{pairedFC} here~~
}
  \item{setScriptPath}{
%%     ~~Describe \code{setScriptPath} here~~
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
function (datafile = datafile, list = list, combination = combination, 
    chipType = "ge", picName = "picName", pairedFC, setScriptPath = "path") 
{
    print("Sourcing files start...")
    if (!missing(setScriptPath)) {
        sourcefiles = list.files(path = setScriptPath, pattern = "*.r", 
            ignore.case = T, full.names = T)
        if (length(sourcefiles) != 0 && sourcefiles != "") {
            sapply(sourcefiles, source, .GlobalEnv)
        }
        else {
            stop("No r files in the set Path!!")
        }
    }
    else {
        if (grepl("linux", Sys.info()["sysname"], ignore.case = T)) {
            source("/home/corona/gaoy/programs/work/stat_Flag.r")
            source("/home/corona/gaoy/programs/work/Filter_Flag.r")
            source("/home/corona/gaoy/programs/work/BatchFoldchange.r")
            source("/home/corona/gaoy/programs/work/BatchFoldchange_HTA.r")
            source("/home/corona/gaoy/programs/work/diff_up_down_gene_count_hta.r")
            source("/home/corona/gaoy/programs/work/PCA_plot.r")
            source("/home/corona/gaoy/programs/work/hclust_plot.r")
        }
        else if (grepl("windows", Sys.info()["sysname"], ignore.case = T)) {
            source("D:/Expression_Scripts/work/stat_Flag.r")
            source("D:/Expression_Scripts/work/Filter_Flag.r")
            source("D:/Expression_Scripts/work/BatchFoldchange.r")
            source("D:/Expression_Scripts/work/BatchFoldchange_HTA.r")
            source("D:/Expression_Scripts/work/diff_up_down_gene_count_hta.r")
            source("D:/Expression_Scripts/work/PCA_plot.r")
            source("D:/Expression_Scripts/work/hclust_plot.r")
        }
        else {
            writeLines("Something is wrong with your system type.\nIf your system is neither linux nor windows, please comment this sentence")
            stop()
        }
    }
    print("Sourcing files End...")
    suffixes = "_NS"
    colname2 = "name"
    colname1 = ""
    suffixes2 = ""
    if (grepl("ge", chipType, ignore.case = T) || grepl("lnc.*|linc.*", 
        chipType, ignore.case = T)) {
        colname1 = "ProbeName"
        suffixes2 = "_flags"
    }
    else if (grepl(paste("affy", "|", "hta", sep = ""), chipType, 
        ignore.case = T)) {
        colname1 = "Probe_Set_ID"
        suffixes2 = "_call"
    }
    else if (grepl("miRNA", chipType, ignore.case = T)) {
        colname1 = "SystematicName"
        suffixes2 = "_flags"
    }
    else {
        stop("You ve assigned a wrong chipType!!ChipType could be ge,affy,hta,miRNA")
    }
    if (!missing(picName)) {
        picName = "Book1.txt_AP_selection"
    }
    else {
        if (picName == "") {
            stop("If u assign a picName, it mustn`t be blank!!")
        }
    }
    if (missing(list)) {
        list = read.delim("list.tsv", colClasses = rep("character", 
            3))
    }
    else if (class(list) == "character") {
        list = read.delim(list)
    }
    else if (class(list) == "data.frame") {
        list = list
    }
    else {
        stop("list should be filename or dataframe!!")
    }
    if (missing(combination)) {
        combination = read.delim("Combination.txt", colClasses = rep("character", 
            3))
    }
    else if (class(combination) == "character") {
        combination = read.delim(combination)
    }
    else if (class(combination) == "data.frame") {
        combination = combination
    }
    else {
        stop("combination should be filename or dataframe!!")
    }
    if (missing(datafile)) {
        datafile = read.delim("Book1.txt")
    }
    if (missing(pairedFC)) {
        pairedFC = FALSE
    }
    if (length(unique(list$Group)) < dim(list)[1]) {
        print("From list I can see there are multi samples in a group")
        print("Let`s do paired foldchange calculation!")
        pairedFC = TRUE
    }
    else {
        pairedFC = FALSE
    }
    if (!grepl("hta", chipType, ignore.case = T)) {
        print(colname1)
        print(colname2)
        print(suffixes2)
        diffAllSub = function(datafile = datafile) {
            if (pairedFC == FALSE) {
                print("UnpairedFC!")
                BatchFoldchange(DataFile = datafile, ContentFile = list, 
                  BatchInfoFile = combination, OutfileSuffixes = "", 
                  FilePath = "./", colname1 = colname1, colname2 = colname2, 
                  prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                  Paired = TRUE, is.log = 1, is.FC = 0)
                BatchFoldchange(DataFile = datafile, ContentFile = list, 
                  BatchInfoFile = combination, OutfileSuffixes = "_fc2", 
                  FilePath = "./", colname1 = colname1, colname2 = colname2, 
                  prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                  Paired = TRUE, is.log = 1, is.FC = 1, envelope = 1, 
                  FC = 2)
                BatchFoldchange(DataFile = datafile, ContentFile = list, 
                  BatchInfoFile = combination, OutfileSuffixes = "_fc3", 
                  FilePath = "./", colname1 = colname1, colname2 = colname2, 
                  prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                  Paired = TRUE, is.log = 1, is.FC = 1, envelope = 1, 
                  FC = 3)
            }
            else {
                print("use pairedFC")
                BatchFoldchange(DataFile = datafile, ContentFile = list, 
                  BatchInfoFile = combination, OutfileSuffixes = "", 
                  FilePath = "./", colname1 = colname1, colname2 = colname2, 
                  prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                  is.log = 1, is.FC = 0)
                BatchFoldchange(DataFile = datafile, ContentFile = list, 
                  BatchInfoFile = combination, OutfileSuffixes = "_fc2", 
                  FilePath = "./", colname1 = colname1, colname2 = colname2, 
                  prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                  is.log = 1, is.FC = 1, envelope = 1, FC = 2)
                BatchFoldchange(DataFile = datafile, ContentFile = list, 
                  BatchInfoFile = combination, OutfileSuffixes = "_fc3", 
                  FilePath = "./", colname1 = colname1, colname2 = colname2, 
                  prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                  is.log = 1, is.FC = 1, envelope = 1, FC = 3)
            }
        }
        if (grepl("lnc.*|linc.*", chipType, ignore.case = T)) {
            if (exists("dataG")) {
                rm(dataG, envir = .GlobalEnv)
            }
            if (length(Sys.glob("Book*mRNA*.txt")) != 0 && length(Sys.glob("Book*lnc*.txt")) != 
                0) {
                dir.create("mRNA", showWarnings = F)
                dir.create("lncRNA", showWarnings = F)
                currentDir = getwd()
                setwd(paste(currentDir, "/mRNA", sep = ""))
                diffAllSub("../Book_mRNA.txt")
                setwd(paste(currentDir, "/lncRNA", sep = ""))
                if (exists("dataG")) {
                  rm(dataG, envir = .GlobalEnv)
                }
                diffAllSub("../Book_lncRNA.txt")
                setwd(currentDir)
            }
            else {
                stop("If you r analysing lncRNA, you should have both Book_mRNA.txt and Book_lnc.txt right here!")
            }
        }
        else {
            diffAllSub(datafile = datafile)
        }
    }
    else {
        if (pairedFC == FALSE) {
            BatchFoldchangeHTA(DataFile = datafile, ContentFile = list, 
                BatchInfoFile = combination, OutfileSuffixes = "", 
                FilePath = "./", colname1 = colname1, colname2 = colname2, 
                prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                Paired = TRUE, is.log = 1, is.FC = 0)
            BatchFoldchangeHTA(DataFile = datafile, ContentFile = list, 
                BatchInfoFile = combination, OutfileSuffixes = "_fc2", 
                FilePath = "./", colname1 = colname1, colname2 = colname2, 
                prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                Paired = TRUE, is.log = 1, is.FC = 1, envelope = 1, 
                FC = 2)
            BatchFoldchangeHTA(DataFile = datafile, ContentFile = list, 
                BatchInfoFile = combination, OutfileSuffixes = "_fc3", 
                FilePath = "./", colname1 = colname1, colname2 = colname2, 
                prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                Paired = TRUE, is.log = 1, is.FC = 1, envelope = 1, 
                FC = 3)
        }
        else {
            BatchFoldchangeHTA(DataFile = datafile, ContentFile = list, 
                BatchInfoFile = combination, OutfileSuffixes = "", 
                FilePath = "./", colname1 = colname1, colname2 = colname2, 
                prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                is.log = 1, is.FC = 0)
            BatchFoldchangeHTA(DataFile = datafile, ContentFile = list, 
                BatchInfoFile = combination, OutfileSuffixes = "_fc2", 
                FilePath = "./", colname1 = colname1, colname2 = colname2, 
                prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                is.log = 1, is.FC = 1, envelope = 1, FC = 2)
            BatchFoldchangeHTA(DataFile = datafile, ContentFile = list, 
                BatchInfoFile = combination, OutfileSuffixes = "_fc3", 
                FilePath = "./", colname1 = colname1, colname2 = colname2, 
                prefixes = "", suffixes = "_NS", suffixes2 = suffixes2, 
                is.log = 1, is.FC = 1, envelope = 1, FC = 3)
        }
        dirs = dir(pattern = "VS")[file.info(dir(pattern = "VS"))$isdir]
        for (i in 1:length(dirs)) {
            files = list.files(path = dirs[i], pattern = "fc")
            diff_up_down_gene_count(datafile = files, filepathin = paste(dirs[i], 
                "/", sep = ""), filepathout = "./")
        }
        scatterplot_NoFlag(combination = combination, list = list, 
            datafile = datafile)
        PCA_plot(datafile = datafile, contentfile = list, suffixes2 = "_NS", 
            filepathin = "./", filepathout = "./")
        hclust_plot(datafile = datafile, contentfile = list, 
            suffixes2 = "_NS", filepathin = "./", filepathout = "./")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
