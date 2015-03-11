\name{batch_ttest_ALL}
\alias{batch_ttest_ALL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
batch_ttest_ALL
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
batch_ttest_ALL(datafile = datafile, list = list, combination = combination, paired, chipType = "chipType", picName = "picName", setScriptPath = "scriptPath")
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
  \item{paired}{
%%     ~~Describe \code{paired} here~~
}
  \item{chipType}{
%%     ~~Describe \code{chipType} here~~
}
  \item{picName}{
%%     ~~Describe \code{picName} here~~
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
    paired, chipType = "chipType", picName = "picName", setScriptPath = "scriptPath") 
{
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
            source("/home/corona/gaoy/programs/work/T_test.r")
            source("/home/corona/gaoy/programs/work/Batch_ttest.r")
            source("/home/corona/gaoy/programs/work/stat_Flag.r")
            source("/home/corona/gaoy/programs/work/Filter_Flag.r")
            source("/home/corona/gaoy/programs/work/Batch_filter.r")
            source("/home/corona/gaoy/programs/work/prettyVolcanoPlot.R")
            source("/home/corona/gaoy/programs/work/heatmap_plot.r")
            source("/home/corona/gaoy/programs/work/prettyheatmap_plot.r")
            source("/home/corona/gaoy/programs/work/PCA_plot.r")
            source("/home/corona/gaoy/programs/work/hclust_plot.r")
            source("/home/corona/gaoy/programs/work/Total_AP_selection.r")
            source("/home/corona/gaoy/programs/work/diff_up_down_gene_count.r")
            source("/home/corona/gaoy/programs/work/cor_plot.r")
            source("/home/corona/gaoy/programs/work/T_test_HTA.r")
            source("/home/corona/gaoy/programs/work/Batch_ttest_HTA.r")
            source("/home/corona/gaoy/programs/work/Batch_filter_HTA.r")
        }
        else if (grepl("windows", Sys.info()["sysname"], ignore.case = T)) {
            source("D://Expression_Scripts//work//T_test.r")
            source("D://Expression_Scripts//work//Batch_ttest.r")
            source("D://Expression_Scripts//work//stat_Flag.r")
            source("D://Expression_Scripts//work//Filter_Flag.r")
            source("D://Expression_Scripts//work//Batch_filter.r")
            source("D://Expression_Scripts//work//prettyVolcanoPlot.R")
            source("D://Expression_Scripts//work//heatmap_plot.r")
            source("/home/corona/gaoy/programs/work/prettyheatmap_plot.r")
            source("D://Expression_Scripts//work//PCA_plot.r")
            source("D://Expression_Scripts//work//hclust_plot.r")
            source("D://Expression_Scripts//work//Total_AP_selection.r")
            source("D://Expression_Scripts//work//diff_up_down_gene_count.r")
            source("D://Expression_Scripts//work//cor_plot.r")
            source("D://Expression_Scripts//work//T_test_HTA.r")
            source("D://Expression_Scripts//work//Batch_ttest_HTA.r")
            source("D://Expression_Scripts//work//Batch_filter_HTA.r")
        }
        else {
            writeLines("Something is wrong with your system type.\nIf your system is neither linux nor windows, please comment this sentence")
            stop()
        }
    }
    suffixes = "_NS"
    colname2 = "name"
    colname1 = ""
    suffixes2 = ""
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
    if (missing(datafile) && !grepl("lnc.*|linc.*", chipType, 
        ignore.case = T)) {
        datafile = read.delim("Book1.txt", na.strings = "", quote = "", 
            sep = "\t", stringsAsFactors = F)
    }
    ttype = 1
    if (!missing(paired)) {
        ttype = 2
    }
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
    if (missing(picName)) {
        picName = "Book1.txt_AP_selection"
    }
    else {
        if (picName == "") {
            stop("If u assign a picName, it mustn`t be blank!!")
        }
    }
    if (!grepl("hta", chipType, ignore.case = T)) {
        ttestSub = function(datafile = datafile) {
            Batch_ttest(datafile = datafile, contentfile = list, 
                Batch_infofile = combination, paired_ttest_fname = "paired.txt", 
                filepathin = "./", filepathout = "./", colname1 = colname1, 
                colname2 = "name", prefixes = "", suffixes = suffixes, 
                suffixes2 = suffixes2, T.type = ttype, is.log = 1, 
                cutoff = 1, is.adjustP = 0)
            Batch_filter(Batch_infofile = combination, contentfile = list, 
                outfile_suffixes = "_p005", filepath = "./", 
                colname2 = "name", prefixes = "", suffixes1 = suffixes, 
                suffixes2 = suffixes2, P_value = 0.05, is.FC = 0, 
                FC = 2, envelope = 1, is.adjustP = 0, is.mean_sel = 0, 
                mean_sel = 6)
            Batch_filter(Batch_infofile = combination, contentfile = list, 
                outfile_suffixes = "_p005fc2", filepath = "./", 
                colname2 = "name", prefixes = "", suffixes1 = suffixes, 
                suffixes2 = suffixes2, P_value = 0.05, is.FC = 1, 
                FC = 2, envelope = 1, is.adjustP = 0, is.mean_sel = 0, 
                mean_sel = 6)
            if (grepl("miRNA", chipType, ignore.case = T)) {
                Batch_filter(Batch_infofile = combination, contentfile = list, 
                  outfile_suffixes = "_p005fc2mean3", filepath = "./", 
                  colname2 = "name", prefixes = "", suffixes1 = suffixes, 
                  suffixes2 = suffixes2, P_value = 0.05, is.FC = 1, 
                  FC = 2, envelope = 1, is.adjustP = 0, is.mean_sel = 1, 
                  mean_sel = 3)
            }
            else {
                Batch_filter(Batch_infofile = combination, contentfile = list, 
                  outfile_suffixes = "_p005fc2mean7", filepath = "./", 
                  colname2 = "name", prefixes = "", suffixes1 = suffixes, 
                  suffixes2 = suffixes2, P_value = 0.05, is.FC = 1, 
                  FC = 2, envelope = 1, is.adjustP = 0, is.mean_sel = 1, 
                  mean_sel = 7)
            }
            Batch_filter(Batch_infofile = combination, contentfile = list, 
                outfile_suffixes = "_p001", filepath = "./", 
                colname2 = "name", prefixes = "", suffixes1 = suffixes, 
                suffixes2 = suffixes2, P_value = 0.01, is.FC = 0, 
                FC = 2, envelope = 1, is.adjustP = 0, is.mean_sel = 0, 
                mean_sel = 6)
            Batch_filter(Batch_infofile = combination, contentfile = list, 
                outfile_suffixes = "_p001fc2", filepath = "./", 
                colname2 = "name", prefixes = "", suffixes1 = suffixes, 
                suffixes2 = suffixes2, P_value = 0.01, is.FC = 1, 
                FC = 2, envelope = 1, is.adjustP = 0, is.mean_sel = 0, 
                mean_sel = 6)
            Batch_volcanoplot(Batch_infofile = combination, pvalue_cut = 0.05, 
                fc_cut = 2)
            scatterplot_NoFlag(combination = combination, list = list, 
                datafile = datafile)
            Total_AP_selection(datafile = datafile, contentfile = list, 
                suffixes2 = suffixes2, filepathin = "./", filepathout = "./")
            file = list.files(pattern = "AP_selection.xls")[1]
            print("hclustering....!")
            hclust_plot(datafile = file, contentfile = list, 
                suffixes2 = suffixes, filepathin = "./", filepathout = "./", 
                picName = picName)
            print("PCA plotting...!")
            PCA_plot(datafile = file, contentfile = list, suffixes2 = suffixes, 
                filepathin = "./", filepathout = "./", picName = picName)
            print("Cor plotting...!")
            cor_plot(datafile = file, contentfile = list, colname1 = colname1, 
                suffixes2 = "_NS", filepathin = "./", filepathout = "./")
            f_names <- list.files(pattern = "VS_.+_p.+txt$")
            print("heatmapping...!")
            prettyheatmap_plot(datafile = f_names, list = list)
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
                ttestSub("../Book_mRNA.txt")
                setwd(paste(currentDir, "/lncRNA", sep = ""))
                if (exists("dataG")) {
                  rm(dataG, envir = .GlobalEnv)
                }
                ttestSub("../Book_lncRNA.txt")
                setwd(currentDir)
            }
            else {
                stop("If you r analysing lncRNA, you should have both Book_mRNA.txt and Book_lnc.txt right here!")
            }
        }
        else {
            print(chipType)
            cat("Class of datafile: ", class(datafile), "\n")
            ttestSub(datafile)
        }
    }
    else {
        Batch_ttestHTA(datafile = datafile, contentfile = list, 
            Batch_infofile = combination, filepathin = "./", 
            filepathout = "./", colname1 = "Probe_Set_ID", colname2 = "name", 
            prefixes = "", suffixes = "_NS", T.type = ttype, 
            is.log = 1, cutoff = 1, is.adjustP = 0)
        Batch_filterHTA(Batch_infofile = combination, contentfile = list, 
            outfile_suffixes = "_p005", filepath = "./", colname2 = "name", 
            prefixes = "", suffixes1 = "_NS", P_value = 0.05, 
            is.FC = 0, FC = 2, envelope = 1, is.adjustP = 0, 
            is.mean_sel = 0, mean_sel = 6)
        Batch_filterHTA(Batch_infofile = combination, contentfile = list, 
            outfile_suffixes = "_p005fc2", filepath = "./", colname2 = "name", 
            prefixes = "", suffixes1 = "_NS", P_value = 0.05, 
            is.FC = 1, FC = 2, envelope = 1, is.adjustP = 0, 
            is.mean_sel = 0, mean_sel = 6)
        Batch_filterHTA(Batch_infofile = combination, contentfile = list, 
            outfile_suffixes = "_p005fc2mean7", filepath = "./", 
            colname2 = "name", prefixes = "", suffixes1 = "_NS", 
            P_value = 0.05, is.FC = 1, FC = 2, envelope = 1, 
            is.adjustP = 0, is.mean_sel = 1, mean_sel = 7)
        Batch_filterHTA(Batch_infofile = combination, contentfile = list, 
            outfile_suffixes = "_p001", filepath = "./", colname2 = "name", 
            prefixes = "", suffixes1 = "_NS", P_value = 0.01, 
            is.FC = 0, FC = 2, envelope = 1, is.adjustP = 0, 
            is.mean_sel = 0, mean_sel = 6)
        Batch_filterHTA(Batch_infofile = combination, contentfile = list, 
            outfile_suffixes = "_p001fc2", filepath = "./", colname2 = "name", 
            prefixes = "", suffixes1 = "_NS", P_value = 0.01, 
            is.FC = 1, FC = 2, envelope = 1, is.adjustP = 0, 
            is.mean_sel = 0, mean_sel = 6)
        Batch_volcanoplot(Batch_infofile = combination, pvalue_cut = 0.05, 
            fc_cut = 2)
        hclust_plot(datafile = datafile, contentfile = list, 
            suffixes2 = "_NS", filepathin = "./", filepathout = "./", 
            picName = picName)
        PCA_plot(datafile = datafile, contentfile = list, suffixes2 = "_NS", 
            filepathin = "./", filepathout = "./", picName = picName)
        cor_plot(datafile = datafile, contentfile = list, colname1 = "Probe_Set_ID", 
            suffixes2 = "_NS", filepathin = "./", filepathout = "./")
        f_names <- list.files(pattern = "VS.+txt$")
        diff_up_down_gene_count(datafile = f_names, filepathin = "./", 
            filepathout = "./")
        f_names <- list.files(pattern = "VS_.+_p.+txt$")
        print("heatmapping...!")
        prettyheatmap_plot(datafile = f_names, list = list)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
