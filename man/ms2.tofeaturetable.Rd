\name{ms2.tofeaturetable}
\alias{ms2.tofeaturetable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ms2.tofeaturetable(data, featureTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{featureTable}{
%%     ~~Describe \code{featureTable} here~~
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
function (data, featureTable) 
{
    matchMS2 <- function(x, featuretable, expandRt = 0, expandMz = 0, 
        ppm = 0) {
        pks <- featuretable
        pks <- cbind(pks, 0)
        colnames(pks)[ncol(pks)] <- "mzmin"
        pks <- cbind(pks, 0)
        colnames(pks)[ncol(pks)] <- "mzmax"
        if (ppm != 0) {
            mz.diff <- pks[, "mz"] * ppm/1e+06
            expandMz <- rep(expandMz, nrow(pks))
        }
        else {
            mz.diff <- rep(0, nrow(pks))
            expandMz <- rep(expandMz, nrow(pks))
        }
        pks[, "mzmin"] <- pks[, "mz"] - expandMz - mz.diff
        pks[, "mzmax"] <- pks[, "mz"] + expandMz + mz.diff
        pks[, "rtmin"] <- pks[, "rt"] - expandRt
        pks[, "rtmax"] <- pks[, "rt"] + expandRt
        peak_ids <- rownames(pks)
        sps <- spectra(x)
        pmz <- precursorMz(x)
        rtm <- rtime(x)
        res <- vector(mode = "list", nrow(pks))
        for (i in 1:nrow(pks)) {
            if (is.na(pks[i, "mz"])) 
                (next)()
            idx <- which(pmz >= pks[i, "mzmin"] & pmz <= pks[i, 
                "mzmax"] & rtm >= pks[i, "rtmin"] & rtm <= pks[i, 
                "rtmax"])
            if (length(idx)) {
                res[[i]] <- lapply(sps[idx], function(z) {
                  z
                })
            }
        }
        names(res) <- peak_ids
        return(res)
    }
    MS2spectra <- matchMS2(data, featureTable, expandRt = 10, 
        expandMz = 0.01, ppm = 0)
    featureTable <- cbind(featureTable, F, 0, 0, 0, 0)
    colnames(featureTable)[(ncol(featureTable) - 4):ncol(featureTable)] <- c("MS2_match", 
        "MS2mz", "MS2int", "PeaksCount", "fromFile")
    for (i in 1:nrow(featureTable)) {
        if (!is.null(MS2spectra[[i]])) {
            tmpSpectra <- MS2spectra[[i]]
            for (j in 1:length(tmpSpectra)) {
                if (tmpSpectra[[j]]@peaksCount == 0) {
                  tmpSpectra[[j]] <- NA
                }
            }
            tmpSpectra <- tmpSpectra[is.na(tmpSpectra) == FALSE]
            if (length(tmpSpectra) > 0) {
                currInt = tmpSpectra[[1]]@precursorIntensity
                currIdx = 1
                for (k in 1:length(tmpSpectra)) {
                  if (tmpSpectra[[k]]@precursorIntensity > currInt) {
                    currIdx = k
                    currInt = tmpSpectra[[k]]@precursorIntensity
                  }
                }
                finalSpectra = tmpSpectra[[currIdx]]
                indices <- finalSpectra@intensity >= 0
                finalSpectra@intensity <- finalSpectra@intensity[indices]
                finalSpectra@mz <- finalSpectra@mz[indices]
                featureTable$MS2_match[i] <- TRUE
                featureTable$MS2mz[i] <- paste(round(finalSpectra@mz, 
                  4), collapse = ";")
                featureTable$MS2int[i] <- paste(finalSpectra@intensity, 
                  collapse = ";")
                featureTable$PeaksCount[i] <- finalSpectra@peaksCount
                featureTable$fromFile[i] <- finalSpectra@fromFile
            }
        }
    }
    featureTable[featureTable == ""] <- 0
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
