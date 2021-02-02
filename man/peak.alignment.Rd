\name{peak.alignment}
\alias{peak.alignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
peak.alignment(xset, bw = 5, mzwid = 0.015, max = 100, quantitative.method = "maxo")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xset}{
%%     ~~Describe \code{xset} here~~
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
}
  \item{mzwid}{
%%     ~~Describe \code{mzwid} here~~
}
  \item{max}{
%%     ~~Describe \code{max} here~~
}
  \item{quantitative.method}{
%%     ~~Describe \code{quantitative.method} here~~
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
function (xset, bw = 5, mzwid = 0.015, max = 100, quantitative.method = "maxo") 
{
    xset <- group(xset, bw = bw, minfrac = 0.5, mzwid = mzwid, 
        minsamp = 1, max = max)
    xset <- retcor(xset, method = "obiwarp", profStep = 1)
    xset <- group(xset, bw = bw, minfrac = 0.5, mzwid = mzwid, 
        minsamp = 1, max = max)
    xset <- fillPeaks(xset)
    XCMt <- data.frame(xset@groups)
    xcmI <- groupval(xset, value = quantitative.method)
    featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, 
        XCMt$rtmax, xcmI)
    colnames(featureTable)[1:4] <- c("mz", "rt", "rtmin", "rtmax")
    featureTable <- as.data.frame(featureTable)
    featureTable <- featureTable[order(featureTable$mz, decreasing = F), 
        ]
    rownames(featureTable) <- paste("F", 1:nrow(featureTable), 
        sep = "")
    return(featureTable)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
