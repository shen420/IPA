library(IPA)

dir <- "X:/Users/Sam_Shen/IPApackageTest20210129/multi"
lib_directory <- "E:/SAM"
lib_name <- "convertedLibraryNeg.msp"
FTdir <- "X:/Users/Sam_Shen/IPApackageTest20210129"
FTname <- "updated301standardsdatabaseforIPAreadcsvTest_20210201.csv"

featureTable <- peak.picking(dir = dir)
featureTable <- find.level3features(data = MSdata, input.files = input.files, xset = xset)
featureTable <- IPA::add.features(input.files = input.files, xset = xset, FTdir = FTdir, FTname = FTname)

featureTable <- peak.alignment(xset = xset)

plot.features(dir = dir, featureTable = featureTable, input.files = input.files, plot.type = 4)

featureTable <- ms2.tofeaturetable(data = MSdata, featureTable = featureTable)
featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.1)

featureTable <- adduct.isotope.annotation(featureTable = featureTable, xset = xset, polarity = "negative")
