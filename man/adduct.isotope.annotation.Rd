\name{adduct.isotope.annotation}
\alias{adduct.isotope.annotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
adduct.isotope.annotation(featureTable, xset, polarity = "positive")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{featureTable}{
%%     ~~Describe \code{featureTable} here~~
}
  \item{xset}{
%%     ~~Describe \code{xset} here~~
}
  \item{polarity}{
%%     ~~Describe \code{polarity} here~~
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
function (featureTable, xset, polarity = "positive") 
{
    xsa <- xsAnnotate(xset)
    anF <- groupFWHM(xsa, perfwhm = 0.6)
    anI <- findIsotopes(anF, mzabs = 0.01)
    anIC <- groupCorr(anI, cor_eic_th = 0.75)
    anFA <- findAdducts(anIC, polarity = polarity)
    peaklist <- getPeaklist(anFA)
    peaklist <- peaklist[order(peaklist$mz), ]
    featureTable <- cbind(featureTable, peaklist$isotopes)
    colnames(featureTable)[ncol(featureTable)] <- "Isotopes"
    featureTable <- cbind(featureTable, peaklist$adduct)
    colnames(featureTable)[ncol(featureTable)] <- "Adduct"
    featureTable <- cbind(featureTable, as.numeric(peaklist$pcgroup))
    colnames(featureTable)[ncol(featureTable)] <- "pcgroup"
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
