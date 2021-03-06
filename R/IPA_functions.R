###############################################################
#This is the main script to perform Integrated data Processing and automated metabolite Annotation (IPA)
#Tao Huan, Sam Shen, Jian Guo 2020-07-28
#Copyright @ University of British Columbia
###############################################################

#XCMS Feature Extraction
peak.picking <- function(dir, mz.tol = 10, ppm=10, peakwidth=c(5,20), mzdiff = 0.01, snthresh = 6, integrate = 1,
                         prefilter = c(3,100), noise = 100){
  #XCMS feature detection
  cwp <- CentWaveParam(ppm = ppm,
                       peakwidth = peakwidth,
                       mzdiff = mzdiff,
                       snthresh = snthresh,
                       integrate = integrate,
                       prefilter = prefilter,
                       noise = noise)

  # Calculate the number of cores
  no_cores <- detectCores() - 1
  print("Using cores:")
  print(no_cores)
  # Initiate cluster
  registerDoParallel(no_cores)

  setwd(dir)
  input.files <- list.files(pattern = ".mzXML")
  data <- readMSData(input.files, centroided = TRUE, mode = "onDisk")
  #PICK LEVEL 1 AND 2 FEATURES
  data <- findChromPeaks(data, param = cwp)
  data_filtered <- filterMsLevel(data, msLevel = 1L)
  xset <- as(data_filtered, 'xcmsSet')
  featureTable <- as.data.frame(xset@peaks)
  featureTable <- cbind(featureTable, 2)
  colnames(featureTable)[ncol(featureTable)] <- "level"

  if(length(input.files) == 1){
    #Label Level 1,2 features
    premass.matrix <- data@featureData@data[data@featureData@data$msLevel == 2 &
                                              data@featureData@data$basePeakIntensity > 0,c(10, 18, 20)]
    colnames(premass.matrix) <- c("rt", "mz", "int")

    is.level12 <- logical(length = nrow(premass.matrix))
    for(i in 1:nrow(premass.matrix)){
      mass.lower.limit <- premass.matrix$mz[i] * (1 - mz.tol * 1e-6)
      mass.upper.limit <- premass.matrix$mz[i] * (1 + mz.tol * 1e-6)
      short.list <- featureTable[featureTable$mz >= mass.lower.limit & featureTable$mz <= mass.upper.limit,]
      short.list <- short.list[short.list$rtmin <= premass.matrix$rt[i] & short.list$rtmax >= premass.matrix$rt[i],]
      if(nrow (short.list) > 0){
        featureTable$level[featureTable$mz >= mass.lower.limit &
                             featureTable$mz <= mass.upper.limit &
                             featureTable$rtmin <= premass.matrix$rt[i] &
                             featureTable$rtmax >= premass.matrix$rt[i]] <- 1
        is.level12[i] <- TRUE
      }
    }
  }else{
    width <- floor(log10(length(input.files))) + 1

    #Label Level 1,2 features
    label12 <- foreach(n = 1:length(input.files)) %dopar% {
      currFeatureTable <- featureTable[featureTable$sample == n,]
      premass.matrix <- data@featureData@data[data@featureData@data$msLevel == 2 &
                                                data@featureData@data$basePeakIntensity > 0,c(10, 18, 20)]
      sample.extract <- grep(paste("F", formatC(n, width=width, flag="0"), sep=""), rownames(premass.matrix))
      premass.matrix <- premass.matrix[sample.extract,]
      colnames(premass.matrix) <- c("rt", "mz", "int")
      level12Vector <- logical(length = nrow(premass.matrix))
      for(i in 1:nrow(premass.matrix)){
        if(is.na(premass.matrix$mz[i])) next
        mass.lower.limit <- premass.matrix$mz[i] * (1 - mz.tol * 1e-6)
        mass.upper.limit <- premass.matrix$mz[i] * (1 + mz.tol * 1e-6)
        short.list <- featureTable[featureTable$mz >= mass.lower.limit &
                                     featureTable$mz <= mass.upper.limit &
                                     featureTable$sample == n,]
        short.list <- short.list[short.list$rtmin <= premass.matrix$rt[i] & short.list$rtmax >= premass.matrix$rt[i],]
        if(nrow (short.list) > 0){
          level12Vector[i] <- TRUE
          currFeatureTable[,12][currFeatureTable[,1] >= mass.lower.limit &
                                  currFeatureTable[,1] <= mass.upper.limit &
                                  currFeatureTable[,5] <= premass.matrix$rt[i] &
                                  currFeatureTable[,6] >= premass.matrix$rt[i]] <- 1
        }
      }
      return(currFeatureTable)
    }

    featureTable <- as.data.frame(bind_rows(label12))
  }

  #clean up the cluster
  stopImplicitCluster()

  featureTable <- featureTable[order(featureTable$mz, decreasing = F ),]
  rownames(featureTable) <- paste("F", 1:nrow(featureTable), sep="")
  xset@peaks <- as.matrix(featureTable)
  featureTable <- featureTable[, c("mz", "rt", "rtmin", "rtmax", "maxo", "sample", "level")]

  input.files <<- input.files
  xset <<- xset
  MSdata <<- data
  return(featureTable)
}

#Performs feature alignment for multi-sample analysis
peak.alignment <- function(xset, bw = 5, mzwid = 0.015, max = 100, quantitative.method = "maxo"){
  #ALIGNMENT
  xset <- group(xset, bw = bw, minfrac = 0.5, mzwid = mzwid, minsamp = 1, max = max)
  xset <- retcor(xset, method = "obiwarp", profStep = 1)
  xset <- group(xset, bw = bw, minfrac = 0.5, mzwid = mzwid, minsamp = 1, max = max)
  xset <- fillPeaks(xset)
  XCMt <- data.frame(xset@groups)
  xcmI <- groupval(xset, value = quantitative.method)
  featureTable <- cbind(XCMt$mzmed, XCMt$rtmed, XCMt$rtmin, XCMt$rtmax, xcmI)
  colnames(featureTable)[1:4] <- c("mz", "rt", "rtmin", "rtmax")
  #Output
  featureTable <- as.data.frame(featureTable)
  featureTable <- featureTable[order(featureTable$mz, decreasing = F ),]
  rownames(featureTable) <- paste("F", 1:nrow(featureTable), sep="")
  return(featureTable)
}

#Add additional level 3 features to featureTable
find.level3features <- function(data, input.files, xset, mz.tol = 10, mass.tol = 0.05, rt.tol = 60,
                                level3.threshold = 2){

  # Calculate the number of cores
  no_cores <- detectCores() - 1
  print("Using cores:")
  print(no_cores)
  # Initiate cluster
  registerDoParallel(no_cores)

  if(length(input.files) == 1){
    #Label Level 1,2 features
    premass.matrix <- data@featureData@data[data@featureData@data$msLevel == 2 &
                                              data@featureData@data$basePeakIntensity > 0,c(10, 18, 20)]
    colnames(premass.matrix) <- c("rt", "mz", "int")

    is.level12 <- logical(length = nrow(premass.matrix))
    for(i in 1:nrow(premass.matrix)){
      mass.lower.limit <- premass.matrix$mz[i] * (1 - mz.tol * 1e-6)
      mass.upper.limit <- premass.matrix$mz[i] * (1 + mz.tol * 1e-6)
      tmpFT <- as.data.frame(xset@peaks)
      short.list <- tmpFT[tmpFT$mz >= mass.lower.limit & tmpFT$mz <= mass.upper.limit,]
      short.list <- short.list[short.list$rtmin <= premass.matrix$rt[i] & short.list$rtmax >= premass.matrix$rt[i],]
      if(nrow (short.list) > 0){
        is.level12[i] <- TRUE
      }
    }
    xraw <- xcmsRaw(input.files[1],profstep=0, mslevel = 1)

    # Confirm level 3 features
    putative.level3 <- premass.matrix[is.level12 == FALSE,]
    level3.matrix <- data.frame(matrix(nrow = nrow(putative.level3), ncol = ncol(tmpFT)))
    colnames(level3.matrix) <- colnames(tmpFT)
    if(nrow(putative.level3 != 0)){
      for (j in 1:nrow(putative.level3)){
        if(putative.level3$mz[j] > xraw@mzrange[2] | putative.level3$mz[j] < xraw@mzrange[1]) next
        if(putative.level3$rt[j] > tail(xraw@scantime, n=1) | putative.level3$rt[j] < xraw@scantime[1]) next
        rt.lower.limit <- putative.level3$rt[j] - rt.tol
        rt.upper.limit <- putative.level3$rt[j] + rt.tol
        mass.lower.limit <- putative.level3$mz[j] - mass.tol
        mass.upper.limit <- putative.level3$mz[j] + mass.tol

        mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
        RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
        eeic <- rawEIC(xraw, mzrange=mzRange, rtrange=RTRange)
        eic.matrix <- eeic[["intensity"]]
        avg.int <- (sum(eic.matrix) - putative.level3$int[j]) / (length(eic.matrix) - 1)
        if(is.na(avg.int)) next
        if(putative.level3$int[j] >= level3.threshold * avg.int){
          #Put level 3 features in featureTable
          level3.matrix[j,1]  <- putative.level3$mz[j]
          level3.matrix[j,2]  <- putative.level3$mz[j]
          level3.matrix[j,3]  <- putative.level3$mz[j]
          level3.matrix[j,4]  <- putative.level3$rt[j]
          level3.matrix[j,5]  <- putative.level3$rt[j]
          level3.matrix[j,6]  <- putative.level3$rt[j]
          level3.matrix[j,7]  <- 0
          level3.matrix[j,8]  <- 0
          level3.matrix[j,9]  <- putative.level3$int[j]
          level3.matrix[j,10]  <- 0
          level3.matrix[j,11] <- 1
          level3.matrix[j,12] <- 3
        }
      }
      level3.matrix <- level3.matrix[is.na(level3.matrix$mz)==FALSE,]

      # Dereplication of level 3 features
      dereplicate.level3 <- data.frame(matrix(ncol = ncol(level3.matrix), nrow = 1))
      colnames(dereplicate.level3) <- colnames(level3.matrix)
      for(q in 1:nrow(level3.matrix)){
        mass.lower.limit <- level3.matrix$mz[q] - 0.01
        mass.upper.limit <- level3.matrix$mz[q] + 0.01
        rt.lower.limit <- level3.matrix$rt[q] - 30
        rt.upper.limit <- level3.matrix$rt[q] + 30
        temp <- dereplicate.level3[(dereplicate.level3$mz >= mass.lower.limit &
                                      dereplicate.level3$mz <= mass.upper.limit &
                                      dereplicate.level3$rt >= rt.lower.limit &
                                      dereplicate.level3$rt <= rt.upper.limit),]
        temp <- temp[complete.cases(temp),]
        if(nrow(temp) == 0) {
          dereplicate.level3 <- rbind(dereplicate.level3, level3.matrix[q,])
        }
      }
      dereplicate.level3 <- dereplicate.level3[complete.cases(dereplicate.level3),]
    }

  }else{
    width <- floor(log10(length(input.files))) + 1
    featureT <- as.data.frame(xset@peaks)

    label12 <- foreach(n = 1:length(input.files)) %dopar% {
      currFeatureTable <- xset@peaks[xset@peaks[,11] == n,]
      premass.matrix <- data@featureData@data[data@featureData@data$msLevel == 2 &
                                                data@featureData@data$basePeakIntensity > 0,c(10, 18, 20)]
      sample.extract <- grep(paste("F", formatC(n, width=width, flag="0"), sep=""), rownames(premass.matrix))
      premass.matrix <- premass.matrix[sample.extract,]
      colnames(premass.matrix) <- c("rt", "mz", "int")
      level12Vector <- logical(length = nrow(premass.matrix))
      for(i in 1:nrow(premass.matrix)){
        if(is.na(premass.matrix$mz[i])) next
        mass.lower.limit <- premass.matrix$mz[i] * (1 - mz.tol * 1e-6)
        mass.upper.limit <- premass.matrix$mz[i] * (1 + mz.tol * 1e-6)
        short.list <- featureT[featureT$mz >= mass.lower.limit &
                                 featureT$mz <= mass.upper.limit &
                                 featureT$sample == n,]
        short.list <- short.list[short.list$rtmin <= premass.matrix$rt[i] & short.list$rtmax >= premass.matrix$rt[i],]
        if(nrow (short.list) > 0){
          level12Vector[i] <- TRUE
        }
      }
      outList <- list(level12Vector, premass.matrix)
      return(outList)
    }
    is.level12List <- list()
    premass.matrixList <- list()

    for(n in 1:length(input.files)){
      is.level12List[[n]] <- label12[[n]][[1]]
      premass.matrixList[[n]] <- label12[[n]][[2]]
    }

    #PICK LEVEL 3 FEATURES AND ADD TO XCMS RESULTS
    level3List <- foreach(n = (1:length(input.files)), .packages = c("xcms", "dplyr")) %dopar% {
      # Find potential level 3 features
      xraw <- xcmsRaw(input.files[n],profstep=0, mslevel = 1)
      premass.matrix <- premass.matrixList[[n]]
      is.level12 <- is.level12List[[n]]
      # Confirm level 3 features
      putative.level3 <- premass.matrix[is.level12 == FALSE,]
      putative.level3 <- putative.level3[is.na(putative.level3$mz) == FALSE,]
      is.level3 <- logical(length = nrow(putative.level3))
      level3.matrix <- data.frame(matrix(nrow = nrow(putative.level3), ncol = ncol(xset@peaks)))
      colnames(level3.matrix) <- colnames(xset@peaks)

      for (j in 1:nrow(putative.level3)){
        if(putative.level3$mz[j] > xraw@mzrange[2] | putative.level3$mz[j] < xraw@mzrange[1]) next
        if(putative.level3$rt[j] > tail(xraw@scantime, n=1) | putative.level3$rt[j] < xraw@scantime[1]) next
        mass.lower.limit <- putative.level3$mz[j] - mass.tol
        mass.upper.limit <- putative.level3$mz[j] + mass.tol
        rt.lower.limit <- putative.level3$rt[j] - rt.tol
        rt.upper.limit <- putative.level3$rt[j] + rt.tol

        mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
        RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
        eeic <- rawEIC(xraw, mzrange=mzRange, rtrange=RTRange)
        eic.matrix <- eeic[["intensity"]]
        avg.int <- (sum(eic.matrix) - putative.level3$int[j]) / (length(eic.matrix) - 1)
        if(is.na(avg.int)) next
        if(putative.level3$int[j] >= level3.threshold * avg.int){
          is.level3 <- TRUE
          #Put level 3 features in featureTable
          level3.matrix[j,1]  <- putative.level3$mz[j]
          level3.matrix[j,2]  <- putative.level3$mz[j]
          level3.matrix[j,3]  <- putative.level3$mz[j]
          level3.matrix[j,4]  <- putative.level3$rt[j]
          level3.matrix[j,5]  <- putative.level3$rt[j]
          level3.matrix[j,6]  <- putative.level3$rt[j]
          level3.matrix[j,7]  <- 0
          level3.matrix[j,8]  <- 0
          level3.matrix[j,9]  <- putative.level3$int[j]
          level3.matrix[j,10]  <- 0
          level3.matrix[j,11] <- n
          level3.matrix[j,12] <- 3
        }
      }
      level3.matrix <- level3.matrix[is.na(level3.matrix$mz)==FALSE,]

      # Dereplication of level 3 features
      dereplicate.level3 <- data.frame(matrix(ncol = ncol(level3.matrix), nrow = 1))
      colnames(dereplicate.level3) <- colnames(level3.matrix)
      for(q in 1:nrow(level3.matrix)){
        mass.lower.limit <- level3.matrix$mz[q] - 0.01
        mass.upper.limit <- level3.matrix$mz[q] + 0.01
        rt.lower.limit <- level3.matrix$rt[q] - 30
        rt.upper.limit <- level3.matrix$rt[q] + 30
        temp <- dereplicate.level3[(dereplicate.level3$mz >= mass.lower.limit &
                                      dereplicate.level3$mz <= mass.upper.limit &
                                      dereplicate.level3$rt >= rt.lower.limit &
                                      dereplicate.level3$rt <= rt.upper.limit),]
        temp <- temp[complete.cases(temp),]
        if(nrow(temp) == 0) {
          dereplicate.level3 <- rbind(dereplicate.level3, level3.matrix[q,])
        }
      }
      dereplicate.level3 <- dereplicate.level3[complete.cases(dereplicate.level3),]
      return(dereplicate.level3)
    }
    dereplicate.level3 <- as.data.frame(bind_rows(level3List))
  }


  #clean up the cluster
  stopImplicitCluster()

  featureTable <- rbind(xset@peaks, as.matrix(dereplicate.level3))
  xset@peaks <- as.matrix(featureTable)
  featureTable <- as.data.frame(featureTable)
  featureTable <- featureTable[order(featureTable$mz, decreasing = F ),]
  rownames(featureTable) <- paste("F", 1:nrow(featureTable), sep="")
  featureTable <- featureTable[, c("mz", "rt", "rtmin", "rtmax", "maxo", "sample", "level")]
  xset <<- xset
  return(featureTable)
}

#Add additional features
add.features <- function(input.files, xset, FTdir, FTname, mass.tol = 0.01, rt.tol = 30){
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  print("Using cores:")
  print(no_cores)
  # Initiate cluster
  registerDoParallel(no_cores)

  setwd(FTdir)
  addFT <- read.csv(file = FTname, header = T, stringsAsFactors = F)
  featureTable <- as.data.frame(xset@peaks)

  if(length(input.files) == 1){
    newFT <- data.frame(matrix(nrow = nrow(addFT), ncol = ncol(featureTable)))
    colnames(newFT) <- colnames(featureTable)
    xraw <- xcmsRaw(filepaths(xset)[1],profstep=0, mslevel = 1)
    for(i in 1:nrow(addFT)){
      if(addFT$mz[i] > xraw@mzrange[2] | addFT$mz[i] < xraw@mzrange[1]) next
      if(addFT$rt[i] > tail(xraw@scantime, n=1) | addFT$rt[i] < xraw@scantime[1]) next
      rt.lower.limit <- addFT$rt[i] - rt.tol
      rt.upper.limit <- addFT$rt[i] + rt.tol
      mass.lower.limit <- addFT$mz[i] - mass.tol
      mass.upper.limit <- addFT$mz[i] + mass.tol

      mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
      RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
      eeic <- rawEIC(xraw, mzrange=mzRange, rtrange=RTRange)
      eic.matrix <- eeic[["intensity"]]
      max.int <- max(eic.matrix)
      if(is.na(max.int)) next()

      newFT[i,1]  <- addFT$mz[i]
      newFT[i,2]  <- addFT$mz[i]
      newFT[i,3]  <- addFT$mz[i]
      newFT[i,4]  <- addFT$rt[i]
      newFT[i,5]  <- addFT$rt[i]
      newFT[i,6]  <- addFT$rt[i]
      newFT[i,7]  <- 0
      newFT[i,8]  <- 0
      newFT[i,9]  <- max.int
      newFT[i,10]  <- 0
      newFT[i,11] <- 1
      newFT[i,12] <- 4
    }

    featureTable <- rbind(featureTable, newFT)
    dereplicatedFT <- data.frame(matrix(ncol = ncol(featureTable), nrow = 0)) #generate data frame with dereplicated features
    colnames(dereplicatedFT) <- colnames(featureTable)
    for(m in (1:nrow(featureTable))) {
      mass.lower.limit <- featureTable$mz[m] - 0.01
      mass.upper.limit <- featureTable$mz[m] + 0.01
      rt.lower.limit <- featureTable$rt[m] - 30
      rt.upper.limit <- featureTable$rt[m] + 30
      temp <- dereplicatedFT[dereplicatedFT$mz >= mass.lower.limit & dereplicatedFT$mz <= mass.upper.limit,]
      temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= rt.upper.limit,]
      if(nrow(temp) == 0) {
        dereplicatedFT[nrow(dereplicatedFT) + 1,] = featureTable[m,]
      }else{
        index <- which(dereplicatedFT$mz == temp$mz[1])[1]
        if(sum(dereplicatedFT$maxo[index]) < temp$maxo[1]){
          dereplicatedFT[index, ] <- temp[1, ]
        }
      }
    }

  }else{

    newFTlist <- foreach(n = (1:length(input.files)), .packages = c("xcms", "dplyr")) %dopar% {
      newFT <- data.frame(matrix(nrow = nrow(addFT), ncol = ncol(featureTable)))
      colnames(newFT) <- colnames(featureTable)
      xraw <- xcmsRaw(filepaths(xset)[n],profstep=0, mslevel = 1)

      for (i in 1:nrow(addFT)){
        if(addFT$mz[i] > xraw@mzrange[2] | addFT$mz[i] < xraw@mzrange[1]) next
        if(addFT$rt[i] > tail(xraw@scantime, n=1) | addFT$rt[i] < xraw@scantime[1]) next
        mass.lower.limit <- addFT$mz[i] - mass.tol
        mass.upper.limit <- addFT$mz[i] + mass.tol
        rt.lower.limit <- addFT$rt[i] - rt.tol
        rt.upper.limit <- addFT$rt[i] + rt.tol

        mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
        RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
        eeic <- rawEIC(xraw, mzrange=mzRange, rtrange=RTRange)
        eic.matrix <- eeic[["intensity"]]
        max.int <- max(eic.matrix)
        if(is.na(max.int)) next()

        newFT[i,1]  <- addFT$mz[i]
        newFT[i,2]  <- addFT$mz[i]
        newFT[i,3]  <- addFT$mz[i]
        newFT[i,4]  <- addFT$rt[i]
        newFT[i,5]  <- addFT$rt[i]
        newFT[i,6]  <- addFT$rt[i]
        newFT[i,7]  <- 0
        newFT[i,8]  <- 0
        newFT[i,9]  <- max.int
        newFT[i,10]  <- 0
        newFT[i,11] <- n
        newFT[i,12] <- 4
      }

      FT <- rbind(featureTable[featureTable$sample == n, ], newFT)
      dereplicatedFT <- data.frame(matrix(ncol = ncol(FT), nrow = 0)) #generate data frame with dereplicated features
      colnames(dereplicatedFT) <- colnames(FT)
      for(m in (1:nrow(FT))) {
        mass.lower.limit <- FT$mz[m] - 0.01
        mass.upper.limit <- FT$mz[m] + 0.01
        rt.lower.limit <- FT$rt[m] - 30
        rt.upper.limit <- FT$rt[m] + 30
        temp <- dereplicatedFT[dereplicatedFT$mz >= mass.lower.limit & dereplicatedFT$mz <= mass.upper.limit,]
        temp <- temp[temp$rt >= rt.lower.limit & temp$rt <= rt.upper.limit,]
        if(nrow(temp) == 0) {
          dereplicatedFT[nrow(dereplicatedFT) + 1,] = FT[m,]
        }else{
          index <- which(dereplicatedFT$mz == temp$mz[1])[1]
          if(sum(dereplicatedFT$maxo[index]) < temp$maxo[1]){
            dereplicatedFT[index, ] <- temp[1, ]
          }
        }
      }
      dereplicatedFT <- dereplicatedFT[complete.cases(dereplicatedFT),]
      return(dereplicatedFT)
    }
    dereplicatedFT <- as.data.frame(bind_rows(newFTlist))
  }

  #clean up the cluster
  stopImplicitCluster()

  featureTable <- dereplicatedFT
  featureTable <- featureTable[order(featureTable$mz, decreasing = F ),]
  xset@peaks <- as.matrix(featureTable)
  rownames(featureTable) <- paste("F", 1:nrow(featureTable), sep="")
  featureTable <- featureTable[, c("mz", "rt", "rtmin", "rtmax", "maxo", "sample", "level")]
  xset <<- xset
  return(featureTable)
}

#Plot level 1,2,3,IPA features
plot.features <- function(dir, featureTable, input.files, plot.type, plotmz.tol = 0.01, plotrt.tol = 60, smooth = 2){
  #peak smooth function
  peak_smooth <- function(x,level=smooth){
    n <- level
    if(length(x) < 2*n){
      return(x)
    } else if(length(unique(x))==1){
      return(x)
    } else{
      y <- vector(length=length(x))
      for(i in 1:n){
        y[i] <- sum(c((n-i+2):(n+1),n:1)*x[1:(i+n)])/sum(c((n-i+2):(n+1),n:1))
      }
      for(i in (n+1):(length(y)-n)){
        y[i] <-  sum(c(1:(n+1),n:1)*x[(i-n):(i+n)])/sum(c(1:(n+1),n:1))
      }
      for(i in (length(y)-n+1):length(y)){
        y[i] <- sum(c(1:n,(n+1):(n+i-length(x)+1))*x[(i-n):length(x)])/sum(c(1:n,(n+1):(n+i-length(x)+1)))
      }
      return(y)
    }
  }

  setwd(dir)

  if(plot.type <= 3){
    xraw <- xcmsRaw(input.files[1],profstep=0, mslevel = 1)
    dir.create(paste0("level_", plot.type))
    setwd(paste0("level_", plot.type))
    plot.matrix <- featureTable[featureTable$level == plot.type,]
    if(nrow(plot.matrix)!=0){
      for(k in 1:nrow(plot.matrix)){
        rt.lower.limit <- plot.matrix$rt[k] - plotrt.tol
        rt.upper.limit <- plot.matrix$rt[k] + plotrt.tol
        mass.lower.limit <- plot.matrix$mz[k] - plotmz.tol
        mass.upper.limit <- plot.matrix$mz[k] + plotmz.tol
        mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
        RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
        eeic <- rawEIC(xraw, mzrange=mzRange, rtrange=RTRange) #extracted EIC object
        points <- cbind(xraw@scantime[eeic$scan], peak_smooth(eeic$intensity))
        png(file = paste0(rownames(plot.matrix)[k], "_",
                          round(plot.matrix$mz[k], digits = 4), "_",
                          round(plot.matrix$rt[k], digits = 0), ".png"),
            width = 480, height = 480)
        eic <- plot(points, type="l", main=paste("Extracted Ion Chromatogram  m/z  ",mzRange[1]," - ",mzRange[2],sep=""), xlab="Seconds",
                    ylab="Intensity", xlim=RTRange)
        dev.off()
      }
    }
  }else{
    dir.create("IPA_EIC")
    setwd("IPA_EIC")
    plot.matrix <- featureTable[,5:(ncol(featureTable))]
    xrawList <- list()
    for(n in 1:length(input.files)){
      xrawList[n] <- xcmsRaw(filepaths(xset)[n],profstep=0)
    }
    for(k in 1:nrow(plot.matrix)){
      rt.lower.limit <- featureTable$rt[k] - plotrt.tol
      rt.upper.limit <- featureTable$rt[k] + plotrt.tol
      mass.lower.limit <- featureTable$mz[k] - plotmz.tol
      mass.upper.limit <- featureTable$mz[k] + plotmz.tol
      maxIndex <- as.numeric(which.max(plot.matrix[k,]))

      mzRange <- as.double(cbind(mass.lower.limit, mass.upper.limit))
      RTRange <- as.integer(cbind(rt.lower.limit, rt.upper.limit))
      eeic <- rawEIC(xrawList[[maxIndex]], mzrange=mzRange, rtrange=RTRange) #extracted EIC object
      points <- cbind(xrawList[[maxIndex]]@scantime[eeic$scan], peak_smooth(eeic$intensity))
      png(file = paste0(rownames(featureTable)[k], "_",
                        round(featureTable$mz[k], digits = 4), "_",
                        round(featureTable$rt[k], digits = 0), ".png"),
          width = 480, height = 480)
      eic <- plot(points, type="l", main=paste("Extracted Ion Chromatogram  m/z  ",mzRange[1]," - ",mzRange[2],sep=""), xlab="Seconds",
                  ylab="Intensity", xlim=RTRange)
      dev.off()
    }
  }

  setwd(dir)
}

#MS2 Assignment
ms2.tofeaturetable <- function(data, featureTable){
  #Multi sample MS2 matching
  matchMS2 <- function(x,featuretable, expandRt = 0, expandMz = 0, ppm = 0) {
    pks <- featuretable
    pks <- cbind(pks, 0)
    colnames(pks)[ncol(pks)] <- "mzmin"
    pks <- cbind(pks, 0)
    colnames(pks)[ncol(pks)] <- "mzmax"
    if (ppm != 0){
      mz.diff <- pks[, "mz"] * ppm / 1e6
      expandMz <- rep(expandMz, nrow(pks))
    }else{
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
      if (is.na(pks[i, "mz"])) next()
      idx <- which(pmz >= pks[i, "mzmin"] & pmz <= pks[i, "mzmax"] &
                     rtm >= pks[i, "rtmin"] & rtm <= pks[i, "rtmax"])
      if (length(idx)) {
        res[[i]] <- lapply(sps[idx], function(z) {
          z
        })
      }
    }
    names(res) <- peak_ids
    return(res)
  }

  MS2spectra <- matchMS2(data, featureTable, expandRt = 10, expandMz = 0.01, ppm = 0)
  featureTable <- cbind(featureTable,F,0,0,0,0)
  colnames(featureTable)[(ncol(featureTable)-4):ncol(featureTable)] <- c("MS2_match", "MS2mz", "MS2int",
                                                                         "PeaksCount", "fromFile")

  for (i in 1:nrow(featureTable)) {
    if(!is.null(MS2spectra[[i]])){
      tmpSpectra <- MS2spectra[[i]]
      for (j in 1:length(tmpSpectra)){
        if(tmpSpectra[[j]]@peaksCount == 0){
          tmpSpectra[[j]] <- NA
        }
      }
      tmpSpectra <- tmpSpectra[is.na(tmpSpectra)==FALSE]
      if(length(tmpSpectra) > 0){
        currInt = tmpSpectra[[1]]@precursorIntensity
        currIdx = 1
        for(k in 1:length(tmpSpectra)){
          if(tmpSpectra[[k]]@precursorIntensity > currInt){
            currIdx = k
            currInt = tmpSpectra[[k]]@precursorIntensity
          }
        }
        finalSpectra = tmpSpectra[[currIdx]]
        indices <- finalSpectra@intensity >= 0;
        finalSpectra@intensity <- finalSpectra@intensity[indices];
        finalSpectra@mz <- finalSpectra@mz[indices];
        featureTable$MS2_match[i] <- TRUE
        featureTable$MS2mz[i] <- paste(round(finalSpectra@mz,4),collapse = ";")
        featureTable$MS2int[i] <- paste(finalSpectra@intensity, collapse = ";")
        featureTable$PeaksCount[i] <- finalSpectra@peaksCount
        featureTable$fromFile[i] <- finalSpectra@fromFile
      }
    }
  }
  featureTable[featureTable == ""] <- 0
  return(featureTable)
}

#Feature Annotation
feature.annotation <- function(featureTable, lib_directory, lib_name, dp = 0.7, ms1.tol = 0.01, ms2.tol = 0.02){
  # Dot product function
  dp.score <- function(x,y){
    if(nrow(x)==0 | nrow(y)==0){return(0)}
    x[,2] <- 100*x[,2]/max(x[,2])
    y[,2] <- 100*y[,2]/max(y[,2])
    alignment <- data.frame(matrix(nrow=nrow(x), ncol=3))
    alignment[,1:2] <- x[,1:2]
    y1 <- y  ##in case one row in y can be selected multiple times
    for(i in 1:nrow(x)){
      mass.diff <- abs(y1[,1] - x[i,1])
      if(min(mass.diff) <= ms2.tol){
        alignment[i,3] <- y1[mass.diff==min(mass.diff),2][1]
        y1[mass.diff==min(mass.diff),1][1] <- NA   # after matched, NA assigned
        y1 <- y1[complete.cases(y1),]
        if(is.null(nrow(y1)) ==TRUE) break
        if(nrow(y1)==0) break
      }
    }
    alignment <- alignment[complete.cases(alignment),]
    if(nrow(alignment)==0){score <- 0}
    if(nrow(alignment)>0){
      #dot product calculation
      AB <- sum(alignment[,2]*alignment[,3])
      A <- sum(x[,2]^2)
      B <- sum(y[,2]^2)
      dp.score <- AB/sqrt(A*B)
      score <- as.numeric(dp.score)
    }
    match_No <- nrow(alignment)
    return  <- c(score,match_No)
    return(return)
  }

  # Calculate the number of cores
  no_cores <- detectCores() - 1
  print("Using cores:")
  print(no_cores)
  # Initiate cluster
  registerDoParallel(no_cores)

  setwd(lib_directory)
  database <- read.msp(lib_name, only.org = FALSE)
  featureTable <- cbind(featureTable, 0)
  colnames(featureTable)[ncol(featureTable)] <- "Annotation"
  featureTable <- cbind(featureTable, 0)
  colnames(featureTable)[ncol(featureTable)] <- "DPscore"

  rez <- foreach(x = 1:nrow(featureTable)) %dopar% {
    premass.Q <- featureTable$mz[x]     ###query precursor ion mass

    if(featureTable$MS2mz[x] == 0){
      # featureTable$Annotation[x] <- "unknown"
      return(c("unknown", 0))
    }

    ms2.Q <- data.frame(m.z = strsplit(featureTable$MS2mz[x], ";")[[1]],
                        int = strsplit(featureTable$MS2int[x], ";")[[1]])  ###query ms2 input, ncol = 2, m.z & int
    ms2.Q$m.z <- as.numeric(as.character(ms2.Q$m.z))
    ms2.Q$int <- as.numeric(as.character(ms2.Q$int))

    output <- data.frame(matrix(ncol=3))
    colnames(output) <- c('std.name','DP.score','match_No')
    h <- 1
    for(i in 1:length(database)){
      if(is.null(database[[i]]$PrecursorMZ)==TRUE) next # no precursor mass

      premass.L <- database[[i]]$PrecursorMZ # database precursor
      if(abs(premass.L-premass.Q) > ms1.tol) next # precursor filter

      ms2.L <- as.data.frame(database[[i]]$pspectrum) # database spectrum
      name.L <- database[[i]]$Name

      output[h,1] <- name.L
      output[h,2] <- dp.score(ms2.Q,ms2.L)[1]
      output[h,3] <- dp.score(ms2.Q,ms2.L)[2]

      h <- h + 1
    }
    output <- output[complete.cases(output),]

    # Dp score threshold, Dp score >= 0.7 , match_No >= 6 (used in GNPS identification)
    output <- output[output[,2] >= dp,]
    # output <- output[output[,3] >= match.number.threshold,]

    if(nrow(output)==0) {
      # featureTable$Annotation[x] <- "unknown"
      return(c("unknown", 0))
    } else {
      output <- output[order(-output[,2]),] # sort by scores
      # featureTable$Annotation[x] <- output[1,1]
      # featureTable$DPscore[x] <- output[1,2]
      return(c(output[1,1], output[1,2]))
    }
  }

  for(x in 1:nrow(featureTable)){
    featureTable[x, c("Annotation", "DPscore")] <- rez[[x]]
  }
  featureTable$DPscore <- as.numeric(featureTable$DPscore)
  #clean up the cluster
  stopImplicitCluster()

  return(featureTable)
}

#CAMERA Annotation
adduct.isotope.annotation <- function(featureTable, xset, polarity = "positive"){
  xsa <- xsAnnotate(xset)
  anF <- groupFWHM(xsa, perfwhm = 0.6)
  anI <- findIsotopes(anF, mzabs = 0.01)
  anIC <- groupCorr(anI, cor_eic_th = 0.75)
  anFA <- findAdducts(anIC, polarity=polarity)
  peaklist <- getPeaklist(anFA)
  peaklist <- peaklist[order(peaklist$mz),]
  featureTable <- cbind(featureTable, peaklist$isotopes)
  colnames(featureTable)[ncol(featureTable)] <- "Isotopes"
  featureTable <- cbind(featureTable, peaklist$adduct)
  colnames(featureTable)[ncol(featureTable)] <- "Adduct"
  featureTable <- cbind(featureTable, as.numeric(peaklist$pcgroup))
  colnames(featureTable)[ncol(featureTable)] <- "pcgroup"
  return(featureTable)
}



