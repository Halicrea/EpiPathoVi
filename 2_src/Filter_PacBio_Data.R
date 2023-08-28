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
  "foreach",
  "doParallel",
  "ranger",
  "tidyverse",
  "kableExtra",
  "stringr"
)

# Load the packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(library(
    package.i,
    character.only = T
  ))
}
n.cores <- parallel::detectCores() - 1
# Create the cluster
my.cluster <- parallel::makeCluster(n.cores, type ="PSOCK")
# And register them to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Check if registered
foreach::getDoParRegistered()



# **********************************************************************
#                       Define global functions                      ----
# **********************************************************************
source("~/Documents/00-R_personal_package/functions.R")
'%ni%' <- Negate('%in%')

# **********************************************************************
#                    _________MAIN SCRIPT__________                 ----
# **********************************************************************
#### - I. Loading data ----
setwd("/home/ntmo83/Dropbox/M2 Ifremer Elyna Bouchereau/bioinfo_epivib_gitlab_m2-elynab/Differences_epi_12-016_02-041_07-115")
df_PACBIO <- read.csv(file = "Methylation_Vaestu_alt.csv", header = T)
df_12_016 <- df_PACBIO[!is.na(df_PACBIO$Score_016),1:11]
df_07_115 <- df_PACBIO[!is.na(df_PACBIO$Score_115),c(1:7,12:15)]
df_02_041 <- df_PACBIO[!is.na(df_PACBIO$Score_041),c(1:7,16:19)]

#### - 1.1. Separate in multiple DATA-FRAME the motifs for a strain
head(df_12_016)

motif_GATC_12_016          <- df_12_016[df_12_016$motif == "GATC" | df_12_016$motif == "GATC (GARANNNNNNCTG,GATC)", 1:3]
motif_CAGNNNNNNTYTC_12_016 <- df_12_016[df_12_016$motif == "CAGNNNNNNTYTC", 1:3]
motif_GARANNNNNNCTG_12_016 <- df_12_016[df_12_016$motif == "GARANNNNNNCTG" | df_12_016$motif == "GARANNNNNNCTG (ACCNNNNNNNTTCY,GARANNNNNNCTG)", 1:3]
motif_GGWCC_12_016         <- df_12_016[df_12_016$motif == "GGWCC", 1:3]
motif_GTAYNNNNGTTA_12_016  <- df_12_016[df_12_016$motif == "GTAYNNNNGTTA", 1:3]
motif_TAACNNNNRTAC_12_016  <- df_12_016[df_12_016$motif == "TAACNNNNRTAC", 1:3]
motif_RGAANNNNNNNGGT_12_016<- df_12_016[df_12_016$motif == "RGAANNNNNNNGGT", 1:3]
motif_ACCNNNNNNNTTCY_12_016<- df_12_016[df_12_016$motif == "ACCNNNNNNNTTCY" | df_12_016$motif == "ACCNNNNNNNTTCY,RGAANNNNNNNGGT", 1:3]

## We calculate and add the end position of the motif
motif_GATC_12_016["end_pos"] <- motif_GATC_12_016["position"]+4
motif_CAGNNNNNNTYTC_12_016["end_pos"] <- motif_CAGNNNNNNTYTC_12_016["position"]+13
motif_GARANNNNNNCTG_12_016["end_pos"] <- motif_GARANNNNNNCTG_12_016["position"]+13
motif_GGWCC_12_016["end_pos"] <- motif_GGWCC_12_016["position"]+5
motif_GTAYNNNNGTTA_12_016["end_pos"] <- motif_GTAYNNNNGTTA_12_016["position"]+12
motif_GATC_12_016["end_pos"] <- motif_GATC_12_016["position"]+12
motif_GATC_12_016["end_pos"] <- motif_GATC_12_016["position"]+14
motif_GATC_12_016["end_pos"] <- motif_GATC_12_016["position"]+14

#### - 1.2. Creating a function to automatise the processus of exporting the coordinates data for each motifs
export_motifs <- function(df, strain_name){
  head(df)
  motif_GATC          <- df[df$motif == "GATC" | df$motif == "GATC (GARANNNNNNCTG,GATC)", 1:3]
  motif_CAGNNNNNNTYTC <- df[df$motif == "CAGNNNNNNTYTC", 1:3]
  motif_GARANNNNNNCTG <- df[df$motif == "GARANNNNNNCTG" | df$motif == "GARANNNNNNCTG (ACCNNNNNNNTTCY,GARANNNNNNCTG)", 1:3]
  motif_GGWCC         <- df[df$motif == "GGWCC", 1:3]
  motif_GTAYNNNNGTTA  <- df[df$motif == "GTAYNNNNGTTA", 1:3]
  motif_TAACNNNNRTAC  <- df[df$motif == "TAACNNNNRTAC", 1:3]
  motif_RGAANNNNNNNGGT<- df[df$motif == "RGAANNNNNNNGGT", 1:3]
  motif_ACCNNNNNNNTTCY<- df[df$motif == "ACCNNNNNNNTTCY" | df$motif == "ACCNNNNNNNTTCY,RGAANNNNNNNGGT", 1:3]
  
  ## We calculate and add the end position of the motif
  motif_GATC["end_pos"] <- motif_GATC["position"]+4
  motif_CAGNNNNNNTYTC["end_pos"] <- motif_CAGNNNNNNTYTC["position"]+13
  motif_GARANNNNNNCTG["end_pos"] <- motif_GARANNNNNNCTG["position"]+13
  motif_GGWCC["end_pos"] <- motif_GGWCC["position"]+5
  motif_GTAYNNNNGTTA["end_pos"] <- motif_GTAYNNNNGTTA["position"]+12
  motif_TAACNNNNRTAC["end_pos"] <- motif_TAACNNNNRTAC["position"]+12
  motif_RGAANNNNNNNGGT["end_pos"] <- motif_RGAANNNNNNNGGT["position"]+14
  motif_ACCNNNNNNNTTCY["end_pos"] <- motif_ACCNNNNNNNTTCY["position"]+14
  
  ## We can then exports those DATAFRAMES
  write.table(motif_GATC,paste("./MOTIFS/","motif_GATC_",strain_name,".txt",sep = ""),
              quote = F, row.names = F, sep = " ", col.names=FALSE)
  write.table(motif_CAGNNNNNNTYTC,paste("./MOTIFS/","motif_CAGNNNNNNTYTC_",strain_name,".txt",sep = ""),
              quote = F, row.names = F, sep = " ", col.names=FALSE)
  write.table(motif_GARANNNNNNCTG,paste("./MOTIFS/","motif_GARANNNNNNCTG_",strain_name,".txt",sep = ""),
              quote = F, row.names = F, sep = " ", col.names=FALSE)
  write.table(motif_GGWCC,paste("./MOTIFS/","motif_GGWCC_",strain_name,".txt",sep = ""),
              quote = F, row.names = F, sep = " ",  col.names=FALSE)
  write.table(motif_GTAYNNNNGTTA,paste("./MOTIFS/","motif_GTAYNNNNGTTA_",strain_name,".txt",sep = ""),
              quote = F, row.names = F, sep = " ",  col.names=FALSE)
  write.table(motif_TAACNNNNRTAC,paste("./MOTIFS/","motif_TAACNNNNRTAC_",strain_name,".txt",sep = ""),
              quote = F, row.names = F, sep = " ", col.names=FALSE)
  write.table(motif_RGAANNNNNNNGGT,paste("./MOTIFS/","motif_RGAANNNNNNNGGT_",strain_name,".txt",sep = ""),
              quote = F, row.names = F, sep = " ", col.names=FALSE)
  write.table(motif_ACCNNNNNNNTTCY,paste("./MOTIFS/","motif_ACCNNNNNNNTTCY_",strain_name,".txt",sep = ""),
              quote = F, row.names = F, sep = " ",  col.names=FALSE)
  print("DONE")
}

#### - 1.3. Executing the function

finding_motifs <- function(genome.path, motif){
  nb.motifs <- dim(motif)[1]
  df.motifs <- data.frame(motif = rep("",nb.motifs),
                          context = rep("",nb.motifs),
                          pos_start = rep(0,nb.motifs),
                          pos_end = rep(0,nb.motifs),
                          contig = rep(2,nb.motifs))
  for(i.motif in 1:nb.motifs){
    motif.pos.tmp = system(paste("../0_bin/Finding_context_of_motif",
                                 genome.path,
                                 motif$context[i.motif],sep=" "),intern = T)  
    motif.pos_contig = strsplit(motif.pos.tmp, " ")
    #print( strsplit(motif.pos.tmp, " "))
    if(length(motif.pos.tmp) != 0){
      df.motifs$pos_start[i.motif] = motif.pos_contig[[1]][1]
      df.motifs$contig[i.motif] = motif.pos_contig[[1]][2]    
    }
  }
  
  return(df.motifs)
}
finding_motifs("../Vaestuarianus_GENOME/Vaestlr12016-contigs.fasta",
               motif_CAGNNNNNNTYTC_12_016)


nb.motifs <- dim(motif_CAGNNNNNNTYTC_12_016)[1]
df.motifs <- data.frame(motif = rep("",nb.motifs),
                       context = rep("",nb.motifs),
                       pos_start = rep(0,nb.motifs),
                       pos_end = rep(0,nb.motifs),
                       contig = rep(2,nb.motifs))

motif.pos.tmp <- {}
for(i.motif in 1:nb.motifs){
  motif.pos.tmp = system(paste("../0_bin/Finding_context_of_motif",
                               "../Vaestuarianus_GENOME/Vaestlr12016-contigs.fasta",
                                motif_CAGNNNNNNTYTC_12_016$context[i.motif],sep=" "),intern = T)  
  motif.pos_contig = strsplit(motif.pos.tmp, " ")
  #print( strsplit(motif.pos.tmp, " "))
  if(length(motif.pos.tmp) != 0){
    df.motifs$pos_start[i.motif] = motif.pos_contig[[1]][1]
    df.motifs$contig[i.motif] = motif.pos_contig[[1]][2]    
  }
}

search_motif_coord <- function(motif,nb.motifs, genome_path){
  df.motifs <- data.frame(motif = rep("",nb.motifs),
                          context = rep("",nb.motifs),
                          pos_start = rep(0,nb.motifs),
                          pos_end = rep(0,nb.motifs),
                          contig = rep(2,nb.motifs))
  
  motif.pos.tmp <- {}
  for(i.motif in 1:nb.motifs){
    motif.pos.tmp = system(paste("../0_bin/Finding_context_of_motif",
                                 genome_path,
                                 motif$context[i.motif],sep=" "),intern = T)  
    motif.pos_contig = strsplit(motif.pos.tmp, " ")
    #print( strsplit(motif.pos.tmp, " "))
    if(length(motif.pos.tmp) != 0){
      df.motifs$pos_start[i.motif] = motif.pos_contig[[1]][1]
      df.motifs$contig[i.motif] = motif.pos_contig[[1]][2]    
    }
  }
  return(df.motifs)
}

search_motif_coord(motif_CAGNNNNNNTYTC_12_016,
                   dim(motif_CAGNNNNNNTYTC_12_016)[1],
                   "../Vaestuarianus_GENOME/Vaestlr12016-contigs.fasta")


i.motifs = 4
motif.pos.tmp =system(paste("../0_bin/Finding_context_of_motif",
             "../Vaestuarianus_GENOME/Vaestlr12016-contigs.fasta",
             motif_CAGNNNNNNTYTC_12_016$context[i.motifs],sep=" "), intern = T)
motif.pos_contig <- strsplit(motif.pos.tmp, " ")
df.motifs$pos_start[i.motif] <- as.numeric(motif.pos_contig[[1]][1])
df.motifs$contig[i.motif] <- as.numeric(motif.pos_contig[[1]][2])


#export_motifs(df_12_016,"12_016")
#export_motifs(df_07_115,"07_115")
#export_motifs(df_02_041,"02_041")
