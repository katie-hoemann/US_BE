# set directories
setwd("C:/Users/Katie/Documents/R/LIWC_Comparisons/")
base_dir <- 'C:/Users/Katie/Documents/R/LIWC_Comparisons/'

# load required libraries (install if necessary)
library(readr)
library(lme4) 
library(lmerTest)
library(haven)
library(tidyverse)
library(psych)
library(emmeans)
library(lsr)


### [BE Dutch vs US English]
# read in data
df <- read.csv("LIWC2015 Results_BE 25 Jan 2022 (400) vs. US 19 May 2022 (397)_AGG.csv",fileEncoding='UTF-8-BOM')

## T-TESTS
# run a series of independent samples t-tests
t <- list()
p <- list()
d <- list()
for (i in names(df[,6:ncol(df)])){
  f <- formula(paste(i,"~Culture"))
  tt <- t.test(f,data=df)
  ef <- cohensD(f,data=df)
  t[i] <- round(tt$statistic,3)
  p[i] <- round(tt$p.value,3)
  d[i] <- round(ef,3)
}
all_results <- cbind(t,p,d)
write.csv(all_results,"results_target_tt.csv")


### [BE English vs US English]
df_wl <- read.csv("LIWC2015 Results_BE 7 Jul 2022 pseudonym translated (400) vs. US 19 May 2022 (397)_AGG.csv",fileEncoding='UTF-8-BOM')

## T-TESTS
# run a series of independent samples t-tests
t <- list()
p <- list()
d <- list()
for (i in names(df_wl[,4:ncol(df_wl)])){
  f <- formula(paste(i,"~Culture"))
  tt <- t.test(f,data=df_wl)
  ef <- cohensD(f,data=df_wl)
  t[i] <- round(tt$statistic,3)
  p[i] <- round(tt$p.value,3)
  d[i] <- round(ef,3)
}
all_results <- cbind(t,p,d)
write.csv(all_results,"results_withinLang_tt.csv")


### [BE English vs BE Dutch]
df_wc <- read.csv("LIWC2015 Results_BE 16 Jun 2022 pseudonym vs. 7 Jul 2022 pseudonym translated (400)_AGG_long.csv",fileEncoding='UTF-8-BOM')

## T-TESTS
# run a series of paired samples t-tests
t <- list()
p <- list()
d <- list()
for (i in names(df_wc[,4:ncol(df_wc)])){
  f <- formula(paste(i,"~Language"))
  tt <- t.test(f,data=df_wc,paired=TRUE)
  ef <- cohensD(f,data=df_wc,method="paired")
  t[i] <- round(tt$statistic,3)
  p[i] <- round(tt$p.value,3)
  d[i] <- round(ef,3)
}
all_results <- cbind(t,p,d)
write.csv(all_results,"results_withinCult_tt.csv")
