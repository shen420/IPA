\name{plot.features}
\alias{plot.features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.features(dir, featureTable, input.files, plot.type, plotmz.tol = 0.01, plotrt.tol = 60, smooth = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{featureTable}{
%%     ~~Describe \code{featureTable} here~~
}
  \item{input.files}{
%%     ~~Describe \code{input.files} here~~
}
  \item{plot.type}{
%%     ~~Describe \code{plot.type} here~~
}
  \item{plotmz.tol}{
%%     ~~Describe \code{plotmz.tol} here~~
}
  \item{plotrt.tol}{
%%     ~~Describe \code{plotrt.tol} here~~
}
  \item{smooth}{
%%     ~~Describe \code{smooth} here~~
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
function (dir, featureTable, input.files, plot.type, plotmz.tol = 0.01, 
    plotrt.tol = 60, smooth = 2) 
{
    peak_smooth <- function(x, level = smooth) {
        n <- level
        if (length(x) < 2 * n) {
            return(x)
        }
        else if (length(unique(x)) == 1) {
            return(x)
        }
        else {
            y <- vector(length = length(x))
            for (i in 1:n) {
                y[i] <- sum(c((n - i + 2):(n + 1), n:1) * x[1:(i + 
                  n)])/sum(c((n - i + 2):(n + 1), n:1))
            }
            for (i in (n + 1):(length(y) - n)) {
                y[i] <- sum(c(1:(n + 1), n:1) * x[(i - n):(i + 
                  n)])/sum(c(1:(n + 1), n:1))
            }
            for (i in (length(y) - n + 1):length(y)) {
                y[i] <- sum(c(1:n, (n + 1):(n + i - length(x) + 
                  1)) * x[(i - n):length(x)])/sum(c(1:n, (n + 
                  1):(n + i - length(x) + 1)))
            }
            return(y)
        }
    }
    setwd(dir)
    if (plot.type <= 3) {
        xraw <- xcmsRaw(input.files[1], profstep = 0, mslevel = 1)
        dir.create(paste0("level_", plot.type))
        setwd(paste0("level_", plot.type))
        plot.matrix <- featureTable[featureTable$level == plot.type, 
            ]
        if (nrow(plot.matrix) != 0) {
            for (k in 1:nrow(plot.matrix)) {
                rt.lower.limit <- plot.matrix$rt[k] - plotrt.tol
                rt.upper.limit <- plot.matrix$rt[k] + plotrt.tol
                mass.lower.limit <- plot.matrix$mz[k] - plotmz.tol
                mass.upper.limit <- plot.matrix$mz[k] + plotmz.tol
                mzRange <- as.double(cbind(mass.lower.limit, 
                  mass.upper.limit))
                RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
                eeic <- rawEIC(xraw, mzrange = mzRange, rtrange = RTRange)
                points <- cbind(xraw@scantime[eeic$scan], peak_smooth(eeic$intensity))
                png(file = paste0(rownames(plot.matrix)[k], "_", 
                  round(plot.matrix$mz[k], digits = 4), "_", 
                  round(plot.matrix$rt[k], digits = 0), ".png"), 
                  width = 480, height = 480)
                eic <- plot(points, type = "l", main = paste("Extracted Ion Chromatogram  m/z  ", 
                  mzRange[1], " - ", mzRange[2], sep = ""), xlab = "Seconds", 
                  ylab = "Intensity", xlim = RTRange)
                dev.off()
            }
        }
    }
    else {
        dir.create("IPA_EIC")
        setwd("IPA_EIC")
        plot.matrix <- featureTable[, 5:(ncol(featureTable))]
        xrawList <- list()
        for (n in 1:length(input.files)) {
            xrawList[n] <- xcmsRaw(filepaths(xset)[n], profstep = 0)
        }
        for (k in 1:nrow(plot.matrix)) {
            rt.lower.limit <- featureTable$rt[k] - plotrt.tol
            rt.upper.limit <- featureTable$rt[k] + plotrt.tol
            mass.lower.limit <- featureTable$mz[k] - plotmz.tol
            mass.upper.limit <- featureTable$mz[k] + plotmz.tol
            maxIndex <- as.numeric(which.max(plot.matrix[k, ]))
            mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
            RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
            eeic <- rawEIC(xrawList[[maxIndex]], mzrange = mzRange, 
                rtrange = RTRange)
            points <- cbind(xrawList[[maxIndex]]@scantime[eeic$scan], 
                peak_smooth(eeic$intensity))
            png(file = paste0(rownames(featureTable)[k], "_", 
                round(featureTable$mz[k], digits = 4), "_", round(featureTable$rt[k], 
                  digits = 0), ".png"), width = 480, height = 480)
            eic <- plot(points, type = "l", main = paste("Extracted Ion Chromatogram  m/z  ", 
                mzRange[1], " - ", mzRange[2], sep = ""), xlab = "Seconds", 
                ylab = "Intensity", xlim = RTRange)
            dev.off()
        }
    }
    setwd(dir)
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
