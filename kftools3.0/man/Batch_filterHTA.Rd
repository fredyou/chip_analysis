\name{Batch_filterHTA}
\alias{Batch_filterHTA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Batch_filterHTA
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Batch_filterHTA(Batch_infofile = "Condition.txt", contentfile = "Content.txt", outfile_suffixes = "_p005", filepath = "", colname2 = "filename", prefixes = "", suffixes1 = "_NS", P_value = 0.05, is.FC = 1, FC = 2, envelope = 1, is.adjustP = 0, is.mean_sel = 1, mean_sel = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Batch_infofile}{
%%     ~~Describe \code{Batch_infofile} here~~
}
  \item{contentfile}{
%%     ~~Describe \code{contentfile} here~~
}
  \item{outfile_suffixes}{
%%     ~~Describe \code{outfile_suffixes} here~~
}
  \item{filepath}{
%%     ~~Describe \code{filepath} here~~
}
  \item{colname2}{
%%     ~~Describe \code{colname2} here~~
}
  \item{prefixes}{
%%     ~~Describe \code{prefixes} here~~
}
  \item{suffixes1}{
%%     ~~Describe \code{suffixes1} here~~
}
  \item{P_value}{
%%     ~~Describe \code{P_value} here~~
}
  \item{is.FC}{
%%     ~~Describe \code{is.FC} here~~
}
  \item{FC}{
%%     ~~Describe \code{FC} here~~
}
  \item{envelope}{
%%     ~~Describe \code{envelope} here~~
}
  \item{is.adjustP}{
%%     ~~Describe \code{is.adjustP} here~~
}
  \item{is.mean_sel}{
%%     ~~Describe \code{is.mean_sel} here~~
}
  \item{mean_sel}{
%%     ~~Describe \code{mean_sel} here~~
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
function (Batch_infofile = "Condition.txt", contentfile = "Content.txt", 
    outfile_suffixes = "_p005", filepath = "", colname2 = "filename", 
    prefixes = "", suffixes1 = "_NS", P_value = 0.05, is.FC = 1, 
    FC = 2, envelope = 1, is.adjustP = 0, is.mean_sel = 1, mean_sel = 6) 
{
    Conditions = ""
    if (class(Batch_infofile) == "data.frame") {
        Conditions = Batch_infofile
    }
    else {
        Conditions <- read.delim(Batch_infofile, sep = "\t", 
            header = T, fill = T, check.names = FALSE)
    }
    Batch_info <- subset(Conditions, Selected == "Y")
    print(paste("The number of iterations is ", as.character(nrow(Batch_info))))
    sampleinfo = ""
    if (class(contentfile) == "data.frame") {
        sampleinfo = contentfile
    }
    else {
        sampleinfo <- read.delim(contentfile, sep = "\t", header = T, 
            fill = T, check.names = FALSE)
    }
    for (i in 1:nrow(Batch_info)) {
        outputdir <- filepath
        filename <- paste(Batch_info[i, 1], "_VS_", Batch_info[i, 
            2], sep = "")
        data = ""
        if (exists(paste(Batch_info[i, 1], "_VS_", Batch_info[i, 
            2], sep = ""))) {
            data = get(filename)
        }
        else {
            data <- read.table(paste(outputdir, "/", filename, 
                ".txt", sep = ""), sep = "\t", quote = "", header = T, 
                fill = T)
            assign(paste(Batch_info[i, 1], "_VS_", Batch_info[i, 
                2], sep = ""), data, envir = .GlobalEnv)
        }
        label1_raw <- subset(sampleinfo, Group == as.character(Batch_info[i, 
            1]), select = colname2)
        label1_raw <- as.matrix(label1_raw)
        label2_raw <- subset(sampleinfo, Group == as.character(Batch_info[i, 
            2]), select = colname2)
        label2_raw <- as.matrix(label2_raw)
        print(paste("the NO is", as.character(i)))
        print(as.character(Batch_info[i, 1]))
        print(label1_raw)
        print(as.character(Batch_info[i, 2]))
        print(label2_raw)
        Filter <- data
        if (is.adjustP == 0) {
            Filter <- subset(Filter, pvalues <= P_value)
        }
        else {
            if ("fdr" \%in\% colnames(result)) {
                Filter <- subset(Filter, fdr <= P_value)
            }
            else {
                stop("The data don't compute djusted p-values for simple multiple testing procedures")
            }
        }
        if (is.FC == 1) {
            Filter <- subset(Filter, (foldchange >= FC | foldchange <= 
                1/FC))
        }
        if (is.mean_sel == 1) {
            flags_tb <- Filter[, grep(suffixes1, colnames(Filter))]
            colnames(flags_tb) <- gsub(suffixes1, "", colnames(flags_tb))
            group_tb_total <- sampleinfo
            group_tb <- group_tb_total[match(colnames(flags_tb), 
                group_tb_total[, "name"]), ]
            total_group <- as.character(unique(group_tb[, "Group"]))
            all_flags <- rep(FALSE, nrow(flags_tb))
            for (group_i in total_group) {
                sample_i <- as.character(group_tb[grep(paste("^", 
                  group_i, "$", sep = ""), group_tb[, "Group"], 
                  perl = T), "name"])
                tmp_tb <- flags_tb[, match(sample_i, colnames(flags_tb))]
                tmp_flags <- apply(tmp_tb, 1, function(x) {
                  mean(x) >= mean_sel
                })
                all_flags <- all_flags | tmp_flags
            }
            Filter <- subset(Filter, all_flags)
        }
        write.table(Filter, file = paste(outputdir, "/", Batch_info[i, 
            1], "_VS_", Batch_info[i, 2], outfile_suffixes, ".xls", 
            sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, 
            row.names = FALSE)
        write.table(Filter, file = paste(outputdir, "/", Batch_info[i, 
            1], "_VS_", Batch_info[i, 2], outfile_suffixes, ".txt", 
            sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, 
            row.names = FALSE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
