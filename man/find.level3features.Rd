\name{find.level3features}
\alias{find.level3features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.level3features(data, input.files, xset, mz.tol = 10, mass.tol = 0.05, rt.tol = 60, level3.threshold = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{input.files}{
%%     ~~Describe \code{input.files} here~~
}
  \item{xset}{
%%     ~~Describe \code{xset} here~~
}
  \item{mz.tol}{
%%     ~~Describe \code{mz.tol} here~~
}
  \item{mass.tol}{
%%     ~~Describe \code{mass.tol} here~~
}
  \item{rt.tol}{
%%     ~~Describe \code{rt.tol} here~~
}
  \item{level3.threshold}{
%%     ~~Describe \code{level3.threshold} here~~
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
function (data, input.files, xset, mz.tol = 10, mass.tol = 0.05, 
    rt.tol = 60, level3.threshold = 2) 
{
    no_cores <- detectCores() - 1
    print("Using cores:")
    print(no_cores)
    registerDoParallel(no_cores)
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
            tmpFT <- as.data.frame(xset@peaks)
            short.list <- tmpFT[tmpFT$mz >= mass.lower.limit & 
                tmpFT$mz <= mass.upper.limit, ]
            short.list <- short.list[short.list$rtmin <= premass.matrix$rt[i] & 
                short.list$rtmax >= premass.matrix$rt[i], ]
            if (nrow(short.list) > 0) {
                is.level12[i] <- TRUE
            }
        }
        xraw <- xcmsRaw(input.files[1], profstep = 0, mslevel = 1)
        putative.level3 <- premass.matrix[is.level12 == FALSE, 
            ]
        level3.matrix <- data.frame(matrix(nrow = nrow(putative.level3), 
            ncol = ncol(tmpFT)))
        colnames(level3.matrix) <- colnames(tmpFT)
        if (nrow(putative.level3 != 0)) {
            for (j in 1:nrow(putative.level3)) {
                if (putative.level3$mz[j] > xraw@mzrange[2] | 
                  putative.level3$mz[j] < xraw@mzrange[1]) 
                  next
                if (putative.level3$rt[j] > tail(xraw@scantime, 
                  n = 1) | putative.level3$rt[j] < xraw@scantime[1]) 
                  next
                rt.lower.limit <- putative.level3$rt[j] - rt.tol
                rt.upper.limit <- putative.level3$rt[j] + rt.tol
                mass.lower.limit <- putative.level3$mz[j] - mass.tol
                mass.upper.limit <- putative.level3$mz[j] + mass.tol
                mzRange <- as.double(cbind(mass.lower.limit, 
                  mass.upper.limit))
                RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
                eeic <- rawEIC(xraw, mzrange = mzRange, rtrange = RTRange)
                eic.matrix <- eeic[["intensity"]]
                avg.int <- (sum(eic.matrix) - putative.level3$int[j])/(length(eic.matrix) - 
                  1)
                if (is.na(avg.int)) 
                  next
                if (putative.level3$int[j] >= level3.threshold * 
                  avg.int) {
                  level3.matrix[j, 1] <- putative.level3$mz[j]
                  level3.matrix[j, 2] <- putative.level3$mz[j]
                  level3.matrix[j, 3] <- putative.level3$mz[j]
                  level3.matrix[j, 4] <- putative.level3$rt[j]
                  level3.matrix[j, 5] <- putative.level3$rt[j]
                  level3.matrix[j, 6] <- putative.level3$rt[j]
                  level3.matrix[j, 7] <- 0
                  level3.matrix[j, 8] <- 0
                  level3.matrix[j, 9] <- putative.level3$int[j]
                  level3.matrix[j, 10] <- 0
                  level3.matrix[j, 11] <- 1
                  level3.matrix[j, 12] <- 3
                }
            }
            level3.matrix <- level3.matrix[is.na(level3.matrix$mz) == 
                FALSE, ]
            dereplicate.level3 <- data.frame(matrix(ncol = ncol(level3.matrix), 
                nrow = 1))
            colnames(dereplicate.level3) <- colnames(level3.matrix)
            for (q in 1:nrow(level3.matrix)) {
                mass.lower.limit <- level3.matrix$mz[q] - 0.01
                mass.upper.limit <- level3.matrix$mz[q] + 0.01
                rt.lower.limit <- level3.matrix$rt[q] - 30
                rt.upper.limit <- level3.matrix$rt[q] + 30
                temp <- dereplicate.level3[(dereplicate.level3$mz >= 
                  mass.lower.limit & dereplicate.level3$mz <= 
                  mass.upper.limit & dereplicate.level3$rt >= 
                  rt.lower.limit & dereplicate.level3$rt <= rt.upper.limit), 
                  ]
                temp <- temp[complete.cases(temp), ]
                if (nrow(temp) == 0) {
                  dereplicate.level3 <- rbind(dereplicate.level3, 
                    level3.matrix[q, ])
                }
            }
            dereplicate.level3 <- dereplicate.level3[complete.cases(dereplicate.level3), 
                ]
        }
    }
    else {
        width <- floor(log10(length(input.files))) + 1
        featureT <- as.data.frame(xset@peaks)
        label12 <- foreach(n = 1:length(input.files)) \%dopar\% 
            {
                currFeatureTable <- xset@peaks[xset@peaks[, 11] == 
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
                  short.list <- featureT[featureT$mz >= mass.lower.limit & 
                    featureT$mz <= mass.upper.limit & featureT$sample == 
                    n, ]
                  short.list <- short.list[short.list$rtmin <= 
                    premass.matrix$rt[i] & short.list$rtmax >= 
                    premass.matrix$rt[i], ]
                  if (nrow(short.list) > 0) {
                    level12Vector[i] <- TRUE
                  }
                }
                outList <- list(level12Vector, premass.matrix)
                return(outList)
            }
        is.level12List <- list()
        premass.matrixList <- list()
        for (n in 1:length(input.files)) {
            is.level12List[[n]] <- label12[[n]][[1]]
            premass.matrixList[[n]] <- label12[[n]][[2]]
        }
        level3List <- foreach(n = (1:length(input.files)), .packages = c("xcms", 
            "dplyr")) \%dopar\% {
            xraw <- xcmsRaw(input.files[n], profstep = 0, mslevel = 1)
            premass.matrix <- premass.matrixList[[n]]
            is.level12 <- is.level12List[[n]]
            putative.level3 <- premass.matrix[is.level12 == FALSE, 
                ]
            putative.level3 <- putative.level3[is.na(putative.level3$mz) == 
                FALSE, ]
            is.level3 <- logical(length = nrow(putative.level3))
            level3.matrix <- data.frame(matrix(nrow = nrow(putative.level3), 
                ncol = ncol(xset@peaks)))
            colnames(level3.matrix) <- colnames(xset@peaks)
            for (j in 1:nrow(putative.level3)) {
                if (putative.level3$mz[j] > xraw@mzrange[2] | 
                  putative.level3$mz[j] < xraw@mzrange[1]) 
                  next
                if (putative.level3$rt[j] > tail(xraw@scantime, 
                  n = 1) | putative.level3$rt[j] < xraw@scantime[1]) 
                  next
                mass.lower.limit <- putative.level3$mz[j] - mass.tol
                mass.upper.limit <- putative.level3$mz[j] + mass.tol
                rt.lower.limit <- putative.level3$rt[j] - rt.tol
                rt.upper.limit <- putative.level3$rt[j] + rt.tol
                mzRange <- as.double(cbind(mass.lower.limit, 
                  mass.upper.limit))
                RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
                eeic <- rawEIC(xraw, mzrange = mzRange, rtrange = RTRange)
                eic.matrix <- eeic[["intensity"]]
                avg.int <- (sum(eic.matrix) - putative.level3$int[j])/(length(eic.matrix) - 
                  1)
                if (is.na(avg.int)) 
                  next
                if (putative.level3$int[j] >= level3.threshold * 
                  avg.int) {
                  is.level3 <- TRUE
                  level3.matrix[j, 1] <- putative.level3$mz[j]
                  level3.matrix[j, 2] <- putative.level3$mz[j]
                  level3.matrix[j, 3] <- putative.level3$mz[j]
                  level3.matrix[j, 4] <- putative.level3$rt[j]
                  level3.matrix[j, 5] <- putative.level3$rt[j]
                  level3.matrix[j, 6] <- putative.level3$rt[j]
                  level3.matrix[j, 7] <- 0
                  level3.matrix[j, 8] <- 0
                  level3.matrix[j, 9] <- putative.level3$int[j]
                  level3.matrix[j, 10] <- 0
                  level3.matrix[j, 11] <- n
                  level3.matrix[j, 12] <- 3
                }
            }
            level3.matrix <- level3.matrix[is.na(level3.matrix$mz) == 
                FALSE, ]
            dereplicate.level3 <- data.frame(matrix(ncol = ncol(level3.matrix), 
                nrow = 1))
            colnames(dereplicate.level3) <- colnames(level3.matrix)
            for (q in 1:nrow(level3.matrix)) {
                mass.lower.limit <- level3.matrix$mz[q] - 0.01
                mass.upper.limit <- level3.matrix$mz[q] + 0.01
                rt.lower.limit <- level3.matrix$rt[q] - 30
                rt.upper.limit <- level3.matrix$rt[q] + 30
                temp <- dereplicate.level3[(dereplicate.level3$mz >= 
                  mass.lower.limit & dereplicate.level3$mz <= 
                  mass.upper.limit & dereplicate.level3$rt >= 
                  rt.lower.limit & dereplicate.level3$rt <= rt.upper.limit), 
                  ]
                temp <- temp[complete.cases(temp), ]
                if (nrow(temp) == 0) {
                  dereplicate.level3 <- rbind(dereplicate.level3, 
                    level3.matrix[q, ])
                }
            }
            dereplicate.level3 <- dereplicate.level3[complete.cases(dereplicate.level3), 
                ]
            return(dereplicate.level3)
        }
        dereplicate.level3 <- as.data.frame(bind_rows(level3List))
    }
    stopImplicitCluster()
    featureTable <- rbind(xset@peaks, as.matrix(dereplicate.level3))
    xset@peaks <- as.matrix(featureTable)
    featureTable <- as.data.frame(featureTable)
    featureTable <- featureTable[order(featureTable$mz, decreasing = F), 
        ]
    rownames(featureTable) <- paste("F", 1:nrow(featureTable), 
        sep = "")
    featureTable <- featureTable[, c("mz", "rt", "rtmin", "rtmax", 
        "maxo", "sample", "level")]
    xset <<- xset
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
