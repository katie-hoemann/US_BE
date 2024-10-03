library(psych)
library(MASS)
library("FactoMineR")
library("factoextra")
library(data.table)


#SET PARAMETERS
dataSet <- "Dutch transcripts" # data set corresponding to folder name
fileDate <- "2024-07-26" # date MEH file was generated
situationType <- "EP" # EP, DP, EN, DN

setwd(paste0("C:/Users/Katie/Documents/MEH/",dataSet,"/MEH/",fileDate," analyses/",situationType))


#LOAD DATA
DF <- read.csv(paste0(fileDate,"_MEH_DTM_Binary_",situationType,".csv"), fileEncoding = 'UTF-8-BOM')
beginningColumn <- 5 # the first column including variables to be included in the PCA
DF[beginningColumn:length(DF)] <- apply(DF[beginningColumn:length(DF)], 2, as.character)
DF[beginningColumn:length(DF)] <- apply(DF[beginningColumn:length(DF)], 2, as.numeric)
colnames(DF)[1] = 'Filename'


#EXPLORATORY PCA & PARALLEL ANALYSIS
explorePCA <- PCA(DF[beginningColumn:length(DF)],graph=FALSE)
exploreEVs <- get_eigenvalue(explorePCA)
exploreEVs <- as.data.frame(exploreEVs)
EVgt1 <- nrow(subset(exploreEVs, exploreEVs$eigenvalue >= 1))
EVgt2 <- nrow(subset(exploreEVs, exploreEVs$eigenvalue >= 2))
fviz_eig(explorePCA, choice="eigenvalue", geom="line", linecolor = "red", ncp=EVgt1)
  