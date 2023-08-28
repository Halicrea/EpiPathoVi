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
  "ggplot2",
  "dplyr"
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
setwd("C:/Users/ebouc/Dropbox/M2 Ifremer Elyna Bouchereau/bioinfo_epivib_gitlab_m2-elynab/Epigenetic")
#### - I. Loading data ----
data.temoin <- read.csv("Temoins_Methylation_6mA-4mC.tsv", header = T, sep = "\t")
data.cuivre <- read.csv("Cuivre_Methylation_6mA-4mC.tsv", header = T, sep = "\t")

merge(data.temoin, data.cuivre, by.x = "start",
      by.y = "start", all.x =T, all.y = T)


data.merged <- merge(data.temoin, data.cuivre, by = c("start","seqid"), all = T)
data.merged <- data.merged[order(data.merged$seqid),]

data.merged$statut[is.na(data.merged$type.x)] <- "differential"
data.merged$statut[is.na(data.merged$type.y)] <- "differential"
data.merged$type.x[is.na(data.merged$type.x)] <- data.merged$type.y[is.na(data.merged$type.x)]
data.merged$end.x[is.na(data.merged$end.x)] <- data.merged$end.y[is.na(data.merged$end.x)]
data.merged$Strand[is.na(data.merged$Strand)] <- data.merged$strand[is.na(data.merged$Strand)]
data.merged$Colonne1.x[is.na(data.merged$Colonne1.x)] <- data.merged$Colonne1.y[is.na(data.merged$Colonne1.x)]
data.merged$context.x[is.na(data.merged$context.x)] <- data.merged$context.y[is.na(data.merged$context.x)]

rows_patch(data.cuivre, data.temoin, by = c("start","seqid"))


write.table(data.merged, file = "Diff_epigenetic_Cuivre-Temoin.tsv",quote=F, sep = "\t")
