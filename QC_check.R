library(ngsReports)

# arrange input data dir
rawDir <- "./QCFiles/1_FastQC/Raw/"
trimmedDir <- "./QCFiles/1_FastQC/Trimmed/"
mapLogDir <- "./QCFiles/barleyMapping" 

rawData <- list.files(rawDir, pattern = "fastqc.zip$", full.names = TRUE) %>% 
  FastqcDataList()
trimmedData <- list.files(trimmedDir, pattern = "fastqc.zip$", full.names = TRUE) %>% 
  FastqcDataList()

# row data 
plotSummary(rawData)
plotSummary(trimmedData)

# GC content
plotGcContent(rawData, plotType = "line", theoreticalGC = FALSE) # kont S2
plotGcContent(trimmedData, plotType = "line", theoreticalGC = FALSE) # kont S2

# Sequence Duplication Levels
plotDupLevels(rawData)
plotDupLevels(trimmedData)

# overrepresent sequnence
plotOverrep(rawData)
plotOverrep(trimmedData)

# STAR mapping
logfiles <- list.files(mapLogDir, pattern = ".out", full.names = TRUE)
starLog <- importNgsLogs(logfiles, type = "star") 
