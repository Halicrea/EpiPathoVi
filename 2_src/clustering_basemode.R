#!/usr/bin/env Rscript
# Date:
# ********************
# Author: Elyna Bouchereau
# Modification:
# ********************
# **********************************************************************
#                           Set arguments                            ----
# **********************************************************************
#args = commandArgs(trailingOnly=TRUE)

# **********************************************************************
#                           Load libraries                          ----
# **********************************************************************
list.of.packages <- c(
  "ggplot2"
)

# Load the packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(library(
    package.i,
    character.only = T
  ))
}



# **********************************************************************
#                       Define global functions                      ----
# **********************************************************************
#source("~/Documents/00-R_personal_package/functions.R")
#'%ni%' <- Negate('%in%')

# **********************************************************************
#                    _________MAIN SCRIPT__________                 ----
# **********************************************************************
setwd("C:/Users/eboucher/Dropbox/M2 Ifremer Elyna Bouchereau/bioinfo_epivib_gitlab_m2-elynab/Epigenetic")
#### - I. Loading data ----
data.temoin <- read.csv("basemodif.Temoin.tsv", header = T, sep="\t")
data.old <- read.csv("Methylation_V_aesturianus.txt", header = T, sep="\t")
plot(data.temoin$coverage ~ data.temoin$score)
plot(data.old$coverage ~ data.old$score)
abline(a = 0, b = 0.51899, col="blue", lwd = 4 )
abline(a = 0, b = 1.2, col="red", lwd = 4 )
data.old.lm <- lm(coverage ~ score, data = data.old)
summary(data.old.lm)
