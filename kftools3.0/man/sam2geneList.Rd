\name{sam2geneList}
\alias{sam2geneList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sam2geneList
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sam2geneList(datafile = "Book1.txt", combination = "Combination.txt", list = "list.tsv", chipType = "affy", filterFlags = FALSE, fc = 2, qvalue = "0.05", anno = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datafile}{
%%     ~~Describe \code{datafile} here~~
}
  \item{combination}{
%%     ~~Describe \code{combination} here~~
}
  \item{list}{
%%     ~~Describe \code{list} here~~
}
  \item{chipType}{
%%     ~~Describe \code{chipType} here~~
}
  \item{filterFlags}{
%%     ~~Describe \code{filterFlags} here~~
}
  \item{fc}{
%%     ~~Describe \code{fc} here~~
}
  \item{qvalue}{
%%     ~~Describe \code{qvalue} here~~
}
  \item{anno}{
%%     ~~Describe \code{anno} here~~
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
function (datafile = "Book1.txt", combination = "Combination.txt", 
    list = "list.tsv", chipType = "affy", filterFlags = FALSE, 
    fc = 2, qvalue = "0.05", anno = FALSE) 
{
    require(samr)
    suffixes = "_NS"
    colname2 = "name"
    colname1 = ""
    suffixes2 = ""
    print("Reading files......")
    if (missing(list)) {
        list = read.delim("list.tsv", stringsAsFactors = F, colClasses = rep("character", 
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
    if (missing(datafile) && !grepl("lnc.*|linc.*", chipType, 
        ignore.case = T)) {
        temp = read.delim("Book1.txt", na.strings = "", quote = "", 
            sep = "\t", stringsAsFactors = F, nrow = 2)
        classes = rep("NULL", ncol(temp))
        signals = grep(suffixes, colnames(temp))
        flags = grep(suffixes2, colnames(temp))
        probes = grep(colname1, colnames(temp))[1]
        if (filterFlags == TRUE) {
            classes[c(probes, flags)] = "character"
        }
        else {
            classes[probes] = "character"
        }
        classes[signals] = "numeric"
        datafile = read.delim("Book1.txt", na.strings = "", quote = "", 
            sep = "\t", stringsAsFactors = F, colClasses = classes)
        if (!missing(anno) && anno == TRUE) {
            classes[-c(signals, flags)] = "character"
            classes[c(signals, flags)] = "NULL"
            annotation = read.delim("Book1.txt", na.strings = "", 
                quote = "", sep = "\t", stringsAsFactors = F, 
                colClasses = classes)
        }
    }
    else {
        stop("If you r analyzing lnc chip, you have to assign datafile")
    }
    rownames(datafile) = datafile[, colname1]
    for (i in 1:nrow(combination)) {
        print("Computing SAM")
        con1 = combination[i, 1]
        con2 = combination[i, 2]
        sublist = paste(subset(list, Group == con1 | Group == 
            con2)[, 2], suffixes, sep = "")
        data = datafile[, sublist]
        x = as.matrix(data)
        y = c(rep(1, length(subset(list, Group == c(con1))[, 
            2])), rep(2, length(subset(list, Group == c(con2))[, 
            2])))
        data = list(x = x, y = y, geneid = as.character(1:nrow(x)), 
            genenames = rownames(x), logged2 = TRUE)
        samr.obj <- samr(data, resp.type = "Two class unpaired", 
            nperms = 200)
        delta = 0.4
        delta.table <- samr.compute.delta.table(samr.obj)
        siggenes.table <- samr.compute.siggenes.table(samr.obj, 
            delta, data, delta.table, all.genes = T)
        qvalue.up <- as.numeric(siggenes.table$genes.up[, ncol(siggenes.table$genes.up)])/100
        qvalue.lo <- as.numeric(siggenes.table$genes.lo[, ncol(siggenes.table$genes.lo)])/100
        foldchange.up <- as.numeric(siggenes.table$genes.up[, 
            ncol(siggenes.table$genes.up) - 1])
        foldchange.lo <- as.numeric(siggenes.table$genes.lo[, 
            ncol(siggenes.table$genes.lo) - 1])
        row.up <- as.numeric(siggenes.table$genes.up[, 3])
        row.lo <- as.numeric(siggenes.table$genes.lo[, 3])
        qvalues <- c(qvalue.up, qvalue.lo)
        foldchanges <- c(foldchange.up, foldchange.lo)
        rows <- c(row.up, row.lo)
        qvalues[rows] = qvalues[1:length(rows)]
        foldchanges[rows] <- foldchanges[1:length(rows)]
        result <- cbind(rownames(x), foldchanges, qvalues)
        result = as.data.frame(result, stringsAsFactors = F)
        colnames(result)[1] = "Probe"
        if (!missing(qvalue)) {
            result <- result[which(result$qvalues < as.numeric(qvalue)), 
                ]
        }
        if (!missing(fc)) {
            result <- result[which(result$foldchanges < as.numeric(1/fc) | 
                result$foldchanges > as.numeric(fc)), ]
        }
        if (!missing(anno) && anno == TRUE) {
            print("Combining annotation....")
            rownames(annotation) = annotation[, colname1]
            result = transform(result, Probe = as.character(Probe))
            result = cbind(result, annotation[result$Probe, ])
        }
        write.table(result, file = paste(con1, "_vs_", con2, 
            "SAM.txt", sep = ""), sep = "\t", quote = FALSE, 
            row.names = FALSE)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
