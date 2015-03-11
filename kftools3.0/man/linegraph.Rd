\name{linegraph}
\alias{linegraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
linegraph
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
linegraph(data = "Book1.txt", xlab = "xaxis", ylab = "yaxis", legendTitle = "Samples")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{legendTitle}{
%%     ~~Describe \code{legendTitle} here~~
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
function (data = "Book1.txt", xlab = "xaxis", ylab = "yaxis", 
    legendTitle = "Samples") 
{
    require(ggplot2)
    require(reshape2)
    if (class(data) == "data.frame") {
        data = data
    }
    else if (class(data) == "character") {
        data = read.delim(data, quote = "", stringsAsFactors = F, 
            check.names = F)
    }
    data = setNames(data.frame(t(data[, -1])), data[, 1])
    data$num = rownames(data)
    data.m = melt(data)
    colours = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", 
        "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", 
        "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", 
        "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", 
        "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", 
        "#599861")
    if (dim(data)[1] > 26) {
        colours = colorRampPalette(c("red", "orange", "green", 
            "purple", "blue"))(dim(data)[1])
    }
    colours = data.frame(colour = sample(colours, dim(data)[1]), 
        group = row.names(data))
    row.names(colours) = colours[, 2]
    data.m$colours = colours[data.m$num, 1]
    colours = data.frame(colour = colours[sample(nrow(colours), 
        dim(data)[1]), ], group = row.names(data))
    row.names(colours) = colours[, 2]
    data.m$colours = colours[data.m$num, 1]
    ggplot(data = data.m, aes(x = variable, y = value, group = num, 
        colour = colours)) + geom_line() + geom_point() + xlab(xlab) + 
        ylab(ylab) + scale_color_discrete(name = legendTitle, 
        breaks = data.m$colours, labels = data.m$num)
    ggsave(filename = "result_line.png", plot = last_plot(), 
        height = par("din")[2] - 3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
