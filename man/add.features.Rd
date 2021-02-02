\name{add.features}
\alias{add.features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
add.features(input.files, xset, FTdir, FTname, mass.tol = 0.01, rt.tol = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input.files}{
%%     ~~Describe \code{input.files} here~~
}
  \item{xset}{
%%     ~~Describe \code{xset} here~~
}
  \item{FTdir}{
%%     ~~Describe \code{FTdir} here~~
}
  \item{FTname}{
%%     ~~Describe \code{FTname} here~~
}
  \item{mass.tol}{
%%     ~~Describe \code{mass.tol} here~~
}
  \item{rt.tol}{
%%     ~~Describe \code{rt.tol} here~~
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
function (input.files, xset, FTdir, FTname, mass.tol = 0.01, 
    rt.tol = 30) 
{
    no_cores <- detectCores() - 1
    print("Using cores:")
    print(no_cores)
    registerDoParallel(no_cores)
    setwd(FTdir)
    addFT <- read.csv(file = FTname, header = T, stringsAsFactors = F)
    featureTable <- as.data.frame(xset@peaks)
    if (length(input.files) == 1) {
        newFT <- data.frame(matrix(nrow = nrow(addFT), ncol = ncol(featureTable)))
        colnames(newFT) <- colnames(featureTable)
        xraw <- xcmsRaw(filepaths(xset)[1], profstep = 0, mslevel = 1)
        for (i in 1:nrow(addFT)) {
            if (addFT$mz[i] > xraw@mzrange[2] | addFT$mz[i] < 
                xraw@mzrange[1]) 
                next
            if (addFT$rt[i] > tail(xraw@scantime, n = 1) | addFT$rt[i] < 
                xraw@scantime[1]) 
                next
            rt.lower.limit <- addFT$rt[i] - rt.tol
            rt.upper.limit <- addFT$rt[i] + rt.tol
            mass.lower.limit <- addFT$mz[i] - mass.tol
            mass.upper.limit <- addFT$mz[i] + mass.tol
            mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
            RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
            eeic <- rawEIC(xraw, mzrange = mzRange, rtrange = RTRange)
            eic.matrix <- eeic[["intensity"]]
            max.int <- max(eic.matrix)
            if (is.na(max.int)) 
                (next)()
            newFT[i, 1] <- addFT$mz[i]
            newFT[i, 2] <- addFT$mz[i]
            newFT[i, 3] <- addFT$mz[i]
            newFT[i, 4] <- addFT$rt[i]
            newFT[i, 5] <- addFT$rt[i]
            newFT[i, 6] <- addFT$rt[i]
            newFT[i, 7] <- 0
            newFT[i, 8] <- 0
            newFT[i, 9] <- max.int
            newFT[i, 10] <- 0
            newFT[i, 11] <- 1
            newFT[i, 12] <- 4
        }
        featureTable <- rbind(featureTable, newFT)
        dereplicatedFT <- data.frame(matrix(ncol = ncol(featureTable), 
            nrow = 0))
        colnames(dereplicatedFT) <- colnames(featureTable)
        for (m in (1:nrow(featureTable))) {
            mass.lower.limit <- featureTable$mz[m] - 0.01
            mass.upper.limit <- featureTable$mz[m] + 0.01
            rt.lower.limit <- featureTable$rt[m] - 30
            rt.upper.limit <- featureTable$rt[m] + 30
            temp <- dereplicatedFT[dereplicatedFT$mz >= mass.lower.limit & 
                dereplicatedFT$mz <= mass.upper.limit, ]
            temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= 
                rt.upper.limit, ]
            if (nrow(temp) == 0) {
                dereplicatedFT[nrow(dereplicatedFT) + 1, ] = featureTable[m, 
                  ]
            }
            else {
                index <- which(dereplicatedFT$mz == temp$mz[1])[1]
                if (sum(dereplicatedFT$maxo[index]) < temp$maxo[1]) {
                  dereplicatedFT[index, ] <- temp[1, ]
                }
            }
        }
    }
    else {
        newFTlist <- foreach(n = (1:length(input.files)), .packages = c("xcms", 
            "dplyr")) \%dopar\% {
            newFT <- data.frame(matrix(nrow = nrow(addFT), ncol = ncol(featureTable)))
            colnames(newFT) <- colnames(featureTable)
            xraw <- xcmsRaw(filepaths(xset)[n], profstep = 0, 
                mslevel = 1)
            for (i in 1:nrow(addFT)) {
                if (addFT$mz[i] > xraw@mzrange[2] | addFT$mz[i] < 
                  xraw@mzrange[1]) 
                  next
                if (addFT$rt[i] > tail(xraw@scantime, n = 1) | 
                  addFT$rt[i] < xraw@scantime[1]) 
                  next
                mass.lower.limit <- addFT$mz[i] - mass.tol
                mass.upper.limit <- addFT$mz[i] + mass.tol
                rt.lower.limit <- addFT$rt[i] - rt.tol
                rt.upper.limit <- addFT$rt[i] + rt.tol
                mzRange <- as.double(cbind(mass.lower.limit, 
                  mass.upper.limit))
                RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
                eeic <- rawEIC(xraw, mzrange = mzRange, rtrange = RTRange)
                eic.matrix <- eeic[["intensity"]]
                max.int <- max(eic.matrix)
                if (is.na(max.int)) 
                  (next)()
                newFT[i, 1] <- addFT$mz[i]
                newFT[i, 2] <- addFT$mz[i]
                newFT[i, 3] <- addFT$mz[i]
                newFT[i, 4] <- addFT$rt[i]
                newFT[i, 5] <- addFT$rt[i]
                newFT[i, 6] <- addFT$rt[i]
                newFT[i, 7] <- 0
                newFT[i, 8] <- 0
                newFT[i, 9] <- max.int
                newFT[i, 10] <- 0
                newFT[i, 11] <- n
                newFT[i, 12] <- 4
            }
            FT <- rbind(featureTable[featureTable$sample == n, 
                ], newFT)
            dereplicatedFT <- data.frame(matrix(ncol = ncol(FT), 
                nrow = 0))
            colnames(dereplicatedFT) <- colnames(FT)
            for (m in (1:nrow(FT))) {
                mass.lower.limit <- FT$mz[m] - 0.01
                mass.upper.limit <- FT$mz[m] + 0.01
                rt.lower.limit <- FT$rt[m] - 30
                rt.upper.limit <- FT$rt[m] + 30
                temp <- dereplicatedFT[dereplicatedFT$mz >= mass.lower.limit & 
                  dereplicatedFT$mz <= mass.upper.limit, ]
                temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= 
                  rt.upper.limit, ]
                if (nrow(temp) == 0) {
                  dereplicatedFT[nrow(dereplicatedFT) + 1, ] = FT[m, 
                    ]
                }
                else {
                  index <- which(dereplicatedFT$mz == temp$mz[1])[1]
                  if (sum(dereplicatedFT$maxo[index]) < temp$maxo[1]) {
                    dereplicatedFT[index, ] <- temp[1, ]
                  }
                }
            }
            dereplicatedFT <- dereplicatedFT[complete.cases(dereplicatedFT), 
                ]
            return(dereplicatedFT)
        }
        dereplicatedFT <- as.data.frame(bind_rows(newFTlist))
    }
    stopImplicitCluster()
    featureTable <- dereplicatedFT
    featureTable <- featureTable[order(featureTable$mz, decreasing = F), 
        ]
    xset@peaks <- as.matrix(featureTable)
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
