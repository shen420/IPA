\name{feature.annotation}
\alias{feature.annotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
feature.annotation(featureTable, lib_directory, lib_name, dp = 0.7, ms1.tol = 0.01, ms2.tol = 0.02)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{featureTable}{
%%     ~~Describe \code{featureTable} here~~
}
  \item{lib_directory}{
%%     ~~Describe \code{lib_directory} here~~
}
  \item{lib_name}{
%%     ~~Describe \code{lib_name} here~~
}
  \item{dp}{
%%     ~~Describe \code{dp} here~~
}
  \item{ms1.tol}{
%%     ~~Describe \code{ms1.tol} here~~
}
  \item{ms2.tol}{
%%     ~~Describe \code{ms2.tol} here~~
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
function (featureTable, lib_directory, lib_name, dp = 0.7, ms1.tol = 0.01, 
    ms2.tol = 0.02) 
{
    dp.score <- function(x, y) {
        if (nrow(x) == 0 | nrow(y) == 0) {
            return(0)
        }
        x[, 2] <- 100 * x[, 2]/max(x[, 2])
        y[, 2] <- 100 * y[, 2]/max(y[, 2])
        alignment <- data.frame(matrix(nrow = nrow(x), ncol = 3))
        alignment[, 1:2] <- x[, 1:2]
        y1 <- y
        for (i in 1:nrow(x)) {
            mass.diff <- abs(y1[, 1] - x[i, 1])
            if (min(mass.diff) <= ms2.tol) {
                alignment[i, 3] <- y1[mass.diff == min(mass.diff), 
                  2][1]
                y1[mass.diff == min(mass.diff), 1][1] <- NA
                y1 <- y1[complete.cases(y1), ]
                if (is.null(nrow(y1)) == TRUE) 
                  break
                if (nrow(y1) == 0) 
                  break
            }
        }
        alignment <- alignment[complete.cases(alignment), ]
        if (nrow(alignment) == 0) {
            score <- 0
        }
        if (nrow(alignment) > 0) {
            AB <- sum(alignment[, 2] * alignment[, 3])
            A <- sum(x[, 2]^2)
            B <- sum(y[, 2]^2)
            dp.score <- AB/sqrt(A * B)
            score <- as.numeric(dp.score)
        }
        match_No <- nrow(alignment)
        return <- c(score, match_No)
        return(return)
    }
    no_cores <- detectCores() - 1
    print("Using cores:")
    print(no_cores)
    registerDoParallel(no_cores)
    setwd(lib_directory)
    database <- read.msp(lib_name, only.org = FALSE)
    featureTable <- cbind(featureTable, 0)
    colnames(featureTable)[ncol(featureTable)] <- "Annotation"
    featureTable <- cbind(featureTable, 0)
    colnames(featureTable)[ncol(featureTable)] <- "DPscore"
    rez <- foreach(x = 1:nrow(featureTable)) \%dopar\% {
        premass.Q <- featureTable$mz[x]
        if (featureTable$MS2mz[x] == 0) {
            return(c("unknown", 0))
        }
        ms2.Q <- data.frame(m.z = strsplit(featureTable$MS2mz[x], 
            ";")[[1]], int = strsplit(featureTable$MS2int[x], 
            ";")[[1]])
        ms2.Q$m.z <- as.numeric(as.character(ms2.Q$m.z))
        ms2.Q$int <- as.numeric(as.character(ms2.Q$int))
        output <- data.frame(matrix(ncol = 3))
        colnames(output) <- c("std.name", "DP.score", "match_No")
        h <- 1
        for (i in 1:length(database)) {
            if (is.null(database[[i]]$PrecursorMZ) == TRUE) 
                next
            premass.L <- database[[i]]$PrecursorMZ
            if (abs(premass.L - premass.Q) > ms1.tol) 
                next
            ms2.L <- as.data.frame(database[[i]]$pspectrum)
            name.L <- database[[i]]$Name
            output[h, 1] <- name.L
            output[h, 2] <- dp.score(ms2.Q, ms2.L)[1]
            output[h, 3] <- dp.score(ms2.Q, ms2.L)[2]
            h <- h + 1
        }
        output <- output[complete.cases(output), ]
        output <- output[output[, 2] >= dp, ]
        if (nrow(output) == 0) {
            return(c("unknown", 0))
        }
        else {
            output <- output[order(-output[, 2]), ]
            return(c(output[1, 1], output[1, 2]))
        }
    }
    for (x in 1:nrow(featureTable)) {
        featureTable[x, c("Annotation", "DPscore")] <- rez[[x]]
    }
    featureTable$DPscore <- as.numeric(featureTable$DPscore)
    stopImplicitCluster()
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
