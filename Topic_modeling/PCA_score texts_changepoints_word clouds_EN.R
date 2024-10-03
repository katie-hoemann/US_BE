library(psych)
library(GPArotation)
library(MASS)
library("FactoMineR")
library("factoextra")
library(data.table)
library(changepoint)
library(dplyr)
library(tidyr)
library(tibble)
library(wordcloud)


#SET PARAMETERS
dataSet <- "English transcripts" # data set corresponding to folder name
fileDate <- "2024-06-20" # date MEH file was generated
situationType <- "DN" # EP, DP, EN, DN
numComponents <- 15 # number of components/factors to extract in PCA
numScores <- 3 # number of high scores to keep when selecting example sentences/participants per situation
numWords <- 10 # number of words to keep when illustrating contents of each component

setwd(paste0("C:/Users/Katie/Documents/MEH/",dataSet,"/MEH/",fileDate," analyses/",situationType))
dirName <- paste0("PCA results_",numComponents)
dir.create(dirName, showWarnings = FALSE)


#LOAD DATA
DF_Binary <- read.csv(paste0(fileDate,"_MEH_DTM_Binary_",situationType,".csv"), fileEncoding = 'UTF-8-BOM')
beginningColumn <- 5 # the first column including variables to be included in the PCA
DF_Binary[beginningColumn:length(DF_Binary)] <- apply(DF_Binary[beginningColumn:length(DF_Binary)], 2, as.character)
DF_Binary[beginningColumn:length(DF_Binary)] <- apply(DF_Binary[beginningColumn:length(DF_Binary)], 2, as.numeric)
colnames(DF_Binary)[1] = "Filename"


#DEFINE FUNCTIONS
kmo = function( data ) {
  X <- cor(as.matrix(data)) 
  iX <- ginv(X) 
  S2 <- diag(diag((iX^-1)))
  AIS <- S2%*%iX%*%S2                      # anti-image covariance matrix
  IS <- X+AIS-2*S2                         # image covariance matrix
  Dai <- sqrt(diag(diag(AIS)))
  IR <- ginv(Dai)%*%IS%*%ginv(Dai)         # image correlation matrix
  AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       # anti-image correlation matrix
  a <- apply((AIR - diag(diag(AIR)))^2, 2, sum)
  AA <- sum(a) 
  b <- apply((X - diag(nrow(X)))^2, 2, sum)
  BB <- sum(b)
  MSA <- b/(b+a)                        # indiv. measures of sampling adequacy
  AIR <- AIR-diag(nrow(AIR))+diag(MSA)  # Examine the anti-image of the correlation matrix. That is the  negative of the partial correlations, partialling out all other variables.
  kmo <- BB/(AA+BB)                     # overall KMO statistic
  # Reporting the conclusion 
  if (kmo >= 0.00 && kmo < 0.50){test <- 'The KMO test yields a POOR degree of common variance.'} 
  else if (kmo >= 0.50 && kmo < 0.60){test <- 'The KMO test yields a SOMEWHAT POOR degree of common variance.'} 
  else if (kmo >= 0.60 && kmo < 0.70){test <- 'The KMO test yields a DECENT degree of common variance.'} 
  else if (kmo >= 0.70 && kmo < 0.80){test <- 'The KMO test yields a GOOD degree of common variance.' } 
  else if (kmo >= 0.80 && kmo < 0.90){test <- 'The KMO test yields a VERY GOOD degree of common variance.' }
  else { test <- 'The KMO test yields a FANTASTIC degree of common variance.' }
  
  ans <- list( overall = kmo,
               report = test,
               individual = MSA,
               AIS = AIS,
               AIR = AIR )
  return(ans)
} 

Bartlett.sphericity.test <- function(x) {
  method <- "Bartlett's test of sphericity"
  data.name <- deparse(substitute(x))
  x <- subset(x, complete.cases(x))
  n <- nrow(x)
  p <- ncol(x)
  chisq <- (1-n+(2*p+5)/6)*log(det(cor(x)))
  df <- p*(p-1)/2
  p.value <- pchisq(chisq, df, lower.tail=FALSE)
  names(chisq) <- "X-squared"
  names(df) <- "df"
  return(structure(list(statistic=chisq, parameter=df, p.value=p.value, method=method, data.name=data.name), class="htest"))
}

remove_zero_var <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(out > 1)
  return(dat[, c(want)])
}

PCAFunction <- function(InputData, NumberOfFactors) {
  InputData = remove_zero_var(InputData) #remove columns with no variance

  options(width=10000)
  options(max.print=2000000000)
  
  cat("\n\nKMO_TEST:\n\n") #run the KMO test and print results
  KMO_Results <- kmo(InputData)
  cat(paste("KMO_METRIC: ", KMO_Results$overall, "\n", KMO_Results$report, sep=""))
  
  cat("\n\nBARTLETT_SPHERICITY_TEST:\n\n") #same for the Bartlett test of sphericity
  print(Bartlett.sphericity.test(InputData))
  
  cat("\n\nFACTOR_ANALYSIS_RESULTS:\n\n") #run the PCA and print results
  PCA <- principal(InputData, nfactors=NumberOfFactors, residuals=FALSE, rotate="varimax", method="regression")
  
  print(PCA)
  return(PCA)
}


#RUN PCA
sink(paste0(dirName, "/", Sys.Date(), "_-_PCA_Results_",situationType,".txt")) #This is where you call the PCA. This saves all results to a file.
PCA_Results = PCAFunction(DF_Binary[beginningColumn:length(DF_Binary)], NumberOfFactors=numComponents)
sink()

DF_PCA_Scores = data.frame(PCA_Results$scores)
DF_PCA_Scores$Filename = as.character(DF_Binary$Filename)

write.csv(PCA_Results$weights, paste0(dirName, "/", Sys.Date(), "_-_PCA_Results - Weights_",situationType,".csv"), fileEncoding = 'UTF-8')
write.csv(PCA_Results$loadings, paste0(dirName, "/", Sys.Date(), "_-_PCA_Results - Loadings_",situationType,".csv"), fileEncoding = 'UTF-8')
write.csv(DF_PCA_Scores, paste0(dirName, "/", Sys.Date(), "_-_PCA_Results - Scores_",situationType,".csv"), row.names=F, fileEncoding = 'UTF-8')


#SCORE TEXTS FOR PCA COMPONENTS
DF_Verbose = data.frame(read.csv(paste0(fileDate,"_MEH_DTM_Verbose_",situationType,".csv"), fileEncoding = "UTF-8-BOM")) #read in the VERBOSE (i.e., "relative frequency") document term matrix
setwd(paste0("C:/Users/Katie/Documents/MEH/",dataSet,"/MEH/",fileDate," analyses/",situationType,"/PCA results_",numComponents)) #jump to PCA results directory

DF_Loadings <- matrix(NA, nrow = length(DF_Binary[beginningColumn:length(DF_Binary)]), ncol = ncol(PCA_Results$loadings)+1) #create an empty matrix for loadings
DF_Loadings[, -1] <- PCA_Results$loadings #fill the loadings matrix with variable names and loadings
DF_Loadings[, 1] <- colnames(DF_Binary[beginningColumn:length(DF_Binary)])
colnames(DF_Loadings) = c('Term',colnames(PCA_Results$loadings))
DF_Loadings <- as.data.frame(DF_Loadings) %>%
  mutate_at(vars(-Term), as.numeric)
DF_Loadings[,c(2:ncol(DF_Loadings))][ abs(DF_Loadings[,c(2:ncol(DF_Loadings))]) < .10] = NA #suppress small values
DF_Scored = data.frame(DF_Verbose$Filename) #create an empty data frame to fill out
colnames(DF_Scored)[1] = 'Filename'
DF_Scored$Filename = as.character(DF_Scored$Filename)

#this loop goes through and calculates everything that we need
for(i in c(2:ncol(DF_Loadings))){
  retain_terms = as.character(DF_Loadings[ , c("Term")][!is.na(DF_Loadings[ , i])])
  term_loadings = DF_Loadings[ , c(i)][!is.na(DF_Loadings[ , i])]
  
  DF_Scored[ ,c(colnames(DF_Loadings)[i])] = 0 #go in and actually score the texts
  
  for(j in c(1:length(retain_terms))){
    if(term_loadings[j] > 0){
      DF_Scored[ ,c(colnames(DF_Loadings)[i])] = DF_Scored[ ,c(colnames(DF_Loadings)[i])] + DF_Verbose[ , c(retain_terms[j])]
    }
    else
    {
      DF_Scored[ ,c(colnames(DF_Loadings)[i])] = DF_Scored[ ,c(colnames(DF_Loadings)[i])] - DF_Verbose[ , c(retain_terms[j])]
    }
  }
  DF_Scored[ ,c(colnames(DF_Loadings)[i])] = DF_Scored[ ,c(colnames(DF_Loadings)[i])] #put on a percentage scale analogous to LIWC
}

#insert separate columns for participant and sentence ID
original_filename <- DF_Verbose$Filename
DF_Verbose <- separate(DF_Verbose, "Filename", c("PPID", "Rest"), "-")
DF_Verbose <- separate(DF_Verbose, "Rest", c("Situation", "Text", "Sentence"), ";")
DF_Verbose <- subset(DF_Verbose, select = -c(Situation, Text))
DF_Verbose <- add_column(DF_Verbose, Filename = original_filename, .before = "PPID")

#load original text, index into it to get top-loading sentences
setwd(paste0("C:/Users/Katie/Documents/MEH/",dataSet,"/For analysis_split by sentence_397 files_corrected_3 Mar 24_C101023 C101033 C101118 C101262 C102189 removed"))
DF_Text = data.frame(readr::read_csv(paste0("ERC English_split by sentence_397 files_corrected_3 Mar 24_",situationType,"_outliers removed.csv"))) #https://stackoverflow.com/questions/19610966/invalid-input-causes-read-csv-to-cut-off-data
setwd(paste0("C:/Users/Katie/Documents/MEH/",dataSet,"/MEH/",fileDate," analyses/",situationType,"/PCA results_",numComponents)) #jump to PCA results directory
DF_Text$Filename <- paste(DF_Text$TextID, DF_Text$Segment, sep=";")

#get highest scoring documents per component
DF_Sentences <- data.frame(matrix(nrow = numComponents, ncol = numScores*2))
DF_Scores <- subset(DF_Scored, select = -Filename)
rownames(DF_Sentences) <- colnames(DF_Scores)
for (col in colnames(DF_Scores)) {
  top_indices <- head(order(DF_Scored[[col]], decreasing = TRUE), numScores)
  top_filenames <- DF_Scored$Filename[top_indices]
  raw_text <- DF_Text$Text[match(top_filenames, DF_Text$Filename)]
  DF_Sentences[col,] <- c(top_filenames, raw_text)
}
colnames(DF_Sentences) <- c(paste0("ID", seq(1, numScores, by=1)),paste0("Text", seq(1, numScores, by=1)))
write.csv(DF_Sentences, paste0("MEM Themes - ",dataSet," - Top Sentences by Relative Frequency_",situationType,".csv"), row.names=T, na='', fileEncoding = 'UTF-8')

#take the original data set, keep only the filenames, participant and sentence IDs, and Word Count and keep the theme scores
DF_ID = DF_Verbose[ , c("Filename", "PPID", "Sentence", "WC")] 
DF_Scored = merge(DF_ID, DF_Scored, by="Filename")
write.csv(DF_Scored, paste0("MEM Themes - ",dataSet," - Scored via Relative Frequencies_",situationType,".csv"), row.names=F, na='', fileEncoding = 'UTF-8') 

#calculate descriptive statistics and save them as an output file
descriptives = describe(DF_Scored)
write.csv(descriptives, paste0("MEM Themes - ",dataSet," - Descriptives_",situationType,".csv"), fileEncoding = 'UTF-8')

#aggregate across sentences
DF_Agg <- DF_Scored %>%
  group_by(DF_Scored$PPID) %>%
  summarise(across(starts_with("RC"), mean))
colnames(DF_Agg) <- c("PPID", paste0("mean_", colnames(DF_Agg)[-1]))
write.csv(DF_Agg, paste0("MEM Themes - ",dataSet," - Scored via Relative Frequencies - Aggregated_",situationType,".csv"), row.names=F, na='', fileEncoding = 'UTF-8') 

#get highest scoring participants per component
DF_Participants <- data.frame(matrix(nrow = numComponents, ncol = numScores))
DF_Scores <- subset(DF_Agg, select = -PPID)
rownames(DF_Participants) <- colnames(DF_Scores)
for (col in colnames(DF_Scores)) {
  top_indices <- head(order(DF_Agg[[col]], decreasing = TRUE), numScores)
  top_filenames <- DF_Agg$PPID[top_indices]
  DF_Participants[col,] <- top_filenames
}
colnames(DF_Participants) <- c(paste0("PP", seq(1, numScores, by=1)))
write.csv(DF_Participants, paste0("MEM Themes - ",dataSet," - Top Participants by Relative Frequency_",situationType,".csv"), row.names=T, na='', fileEncoding = 'UTF-8')


#CHANGEPOINTS ANALYSIS
DF_Changepoints = DF_Scored
for (i_comp in 1:numComponents) {
  colName <- paste0("RC",i_comp)
  colData <- DF_Changepoints[,colName] 
  cpoints = cpt.mean(sort(colData[!is.na(colData)]), method="AMOC", Q=1); cp_low = cpoints@param.est$mean[1]; cp_high = cpoints@param.est$mean[length(cpoints@param.est$mean)];
  DF_Changepoints[,colName] = 0; DF_Changepoints[colData >= cp_high,colName] = 1;
}
write.csv(DF_Changepoints, paste0("MEM Themes - ",dataSet," - Scored via Relative Frequencies - Changepoints High_",situationType,".csv"), row.names = F, na='', fileEncoding = 'UTF-8')

#aggregate across sentences
DF_Agg <- DF_Changepoints %>%
  group_by(DF_Changepoints$PPID) %>%
  summarise(across(starts_with("RC"), mean))
colnames(DF_Agg) <- c("PPID", paste0("mean_", colnames(DF_Agg)[-1]))
write.csv(DF_Agg, paste0("MEM Themes - ",dataSet," - Scored via Relative Frequencies - Changepoints High - Aggregated_",situationType,".csv"), row.names=F, na='', fileEncoding = 'UTF-8')


#WORD CLOUDS (https://towardsdatascience.com/create-a-word-cloud-with-r-bde3e7422e8a)
DF_Word_Clouds <- matrix(NA, nrow = length(DF_Binary[beginningColumn:length(DF_Binary)]), ncol = ncol(PCA_Results$loadings)+1) #create an empty matrix for loadings
DF_Word_Clouds[, -1] <- PCA_Results$loadings #fill the loadings matrix with variable names and loadings
DF_Word_Clouds[, 1] <- colnames(DF_Binary[beginningColumn:length(DF_Binary)])
colnames(DF_Word_Clouds) = c('Term',colnames(PCA_Results$loadings))
DF_Word_Clouds <- as.data.frame(DF_Word_Clouds) %>%
  mutate_at(vars(-Term), as.numeric)
DF_Word_Clouds[,c(2:ncol(DF_Word_Clouds))][DF_Word_Clouds[,c(2:ncol(DF_Word_Clouds))] < 0] = 0 #suppress negative values
dirName <- paste0(Sys.Date(),"_Word clouds")
dir.create(dirName, showWarnings = FALSE)
for (col in colnames(DF_Word_Clouds)[-1]) {
  words <- DF_Word_Clouds$Term[DF_Word_Clouds[[col]] > 0]
  loadings <- DF_Word_Clouds[[col]][DF_Word_Clouds[[col]] > 0]*100
  output <- file.path(dirName, paste("wordcloud", col, ".png", sep = "_"))
  png(output)
  wordcloud(words=words, freq=loadings, scale=c(3, 0.25), min.freq = 10, random.order=FALSE, rot.per=0, colors=brewer.pal(4, "Dark2"))
  title(main = col)
  dev.off()
}


#SUMMARY FILE
DF_Summary <- DF_Sentences
top_words_string <- character()
for (row in rownames(DF_Summary)) {
  top_indices <- head(order(DF_Word_Clouds[[row]], decreasing = TRUE), numWords)
  top_words <- DF_Word_Clouds$Term[top_indices]
  top_words_string <- append(top_words_string, paste(top_words, collapse = ", "))
}
DF_Summary <- add_column(DF_Summary, TopWords = top_words_string)
write.csv(DF_Summary, paste0("MEM Themes - ",dataSet," - Summary_",situationType,".csv"), row.names=T, na='', fileEncoding = 'UTF-8')
