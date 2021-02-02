\name{peak.picking}
\alias{peak.picking}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
peak.picking(dir, mz.tol = 10, ppm = 10, peakwidth = c(5, 20), mzdiff = 0.01, snthresh = 6, integrate = 1, prefilter = c(3, 100), noise = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{mz.tol}{
%%     ~~Describe \code{mz.tol} here~~
}
  \item{ppm}{
%%     ~~Describe \code{ppm} here~~
}
  \item{peakwidth}{
%%     ~~Describe \code{peakwidth} here~~
}
  \item{mzdiff}{
%%     ~~Describe \code{mzdiff} here~~
}
  \item{snthresh}{
%%     ~~Describe \code{snthresh} here~~
}
  \item{integrate}{
%%     ~~Describe \code{integrate} here~~
}
  \item{prefilter}{
%%     ~~Describe \code{prefilter} here~~
}
  \item{noise}{
%%     ~~Describe \code{noise} here~~
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
function (dir, mz.tol = 10, ppm = 10, peakwidth = c(5, 20), mzdiff = 0.01, 
    snthresh = 6, integrate = 1, prefilter = c(3, 100), noise = 100) 
{
    cwp <- CentWaveParam(ppm = ppm, peakwidth = peakwidth, mzdiff = mzdiff, 
        snthresh = snthresh, integrate = integrate, prefilter = prefilter, 
        noise = noise)
    no_cores <- detectCores() - 1
    print("Using cores:")
    print(no_cores)
    registerDoParallel(no_cores)
    setwd(dir)
    input.files <- list.files(pattern = ".mzXML")
    data <- readMSData(input.files, centroided = TRUE, mode = "onDisk")
    data <- findChromPeaks(data, param = cwp)
    data_filtered <- filterMsLevel(data, msLevel = 1L)
    xset <- as(data_filtered, "xcmsSet")
    featureTable <- as.data.frame(xset@peaks)
    featureTable <- cbind(featureTable, 2)
    colnames(featureTable)[ncol(featureTable)] <- "level"
    if (length(input.files) == 1) {
        premass.matrix <- data@featureData@data[data@featureData@data$msLevel == 
            2 & data@featureData@data$basePeakIntensity > 0, 
            c(10, 18, 20)]
        colnames(premass.matrix) <- c("rt", "mz", "int")
        is.level12 <- logical(length = nrow(premass.matrix))
        for (i in 1:nrow(premass.matrix)) {
            mass.lower.limit <- premass.matrix$mz[i] * (1 - mz.tol * 
                1e-06)
            mass.upper.limit <- premass.matrix$mz[i] * (1 + mz.tol * 
                1e-06)
            short.list <- featureTable[featureTable$mz >= mass.lower.limit & 
                featureTable$mz <= mass.upper.limit, ]
            short.list <- short.list[short.list$rtmin <= premass.matrix$rt[i] & 
                short.list$rtmax >= premass.matrix$rt[i], ]
            if (nrow(short.list) > 0) {
                featureTable$level[featureTable$mz >= mass.lower.limit & 
                  featureTable$mz <= mass.upper.limit & featureTable$rtmin <= 
                  premass.matrix$rt[i] & featureTable$rtmax >= 
                  premass.matrix$rt[i]] <- 1
                is.level12[i] <- TRUE
            }
        }
    }
    else {
        width <- floor(log10(length(input.files))) + 1
        label12 <- foreach(n = 1:length(input.files)) \%dopar\% 
            {
                currFeatureTable <- featureTable[featureTable$sample == 
                  n, ]
                premass.matrix <- data@featureData@data[data@featureData@data$msLevel == 
                  2 & data@featureData@data$basePeakIntensity > 
                  0, c(10, 18, 20)]
                sample.extract <- grep(paste("F", formatC(n, 
                  width = width, flag = "0"), sep = ""), rownames(premass.matrix))
                premass.matrix <- premass.matrix[sample.extract, 
                  ]
                colnames(premass.matrix) <- c("rt", "mz", "int")
                level12Vector <- logical(length = nrow(premass.matrix))
                for (i in 1:nrow(premass.matrix)) {
                  if (is.na(premass.matrix$mz[i])) 
                    next
                  mass.lower.limit <- premass.matrix$mz[i] * 
                    (1 - mz.tol * 1e-06)
                  mass.upper.limit <- premass.matrix$mz[i] * 
                    (1 + mz.tol * 1e-06)
                  short.list <- featureTable[featureTable$mz >= 
                    mass.lower.limit & featureTable$mz <= mass.upper.limit & 
                    featureTable$sample == n, ]
                  short.list <- short.list[short.list$rtmin <= 
                    premass.matrix$rt[i] & short.list$rtmax >= 
                    premass.matrix$rt[i], ]
                  if (nrow(short.list) > 0) {
                    level12Vector[i] <- TRUE
                    currFeatureTable[, 12][currFeatureTable[, 
                      1] >= mass.lower.limit & currFeatureTable[, 
                      1] <= mass.upper.limit & currFeatureTable[, 
                      5] <= premass.matrix$rt[i] & currFeatureTable[, 
                      6] >= premass.matrix$rt[i]] <- 1
                  }
                }
                return(currFeatureTable)
            }
        featureTable <- as.data.frame(bind_rows(label12))
    }
    stopImplicitCluster()
    featureTable <- featureTable[order(featureTable$mz, decreasing = F), 
        ]
    rownames(featureTable) <- paste("F", 1:nrow(featureTable), 
        sep = "")
    xset@peaks <- as.matrix(featureTable)
    featureTable <- featureTable[, c("mz", "rt", "rtmin", "rtmax", 
        "maxo", "sample", "level")]
    input.files <<- input.files
    xset <<- xset
    MSdata <<- data
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
