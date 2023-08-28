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
library(ggplot2)
library(GenomicRanges)
library(ggbio)
# **********************************************************************
#                       Define global functions                      ----
# **********************************************************************
source("~/Documents/00-R_personal_package/functions.R")
'%ni%' <- Negate('%in%')

# **********************************************************************
#                    _________MAIN SCRIPT__________                 ----
# **********************************************************************
set.seed(1); N <- 100; gr <- GRanges(seqnames = sample(c("chr1", "chr2", "chr3"), size = N, replace = TRUE), IRanges(start = sample(1:300, size = N, replace = TRUE), width = sample(70:75, size = N,replace = TRUE)), strand = sample(c("+", "-"), size = N, replace = TRUE), value = rnorm(N, 10, 3), score = rnorm(N, 100, 30), sample = sample(c("Normal", "Tumor"), size = N, replace = TRUE), pair = sample(letters, size = N, replace = TRUE))
autoplot(gr, aes(color = strand, fill = strand), facets = strand ~ seqnames, stat = "coverage")
ggplot(gr) + 
  layout_circle(aes(fill = seqnames), geom = "rect")
seqlengths(gr) <- c(400, 500, 700)
values(gr)$to.gr <- gr[sample(1:length(gr), size = length(gr))]
idx <- sample(1:length(gr), size = 50)
gr <- gr[idx]
ggplot() + 
  layout_circle(gr, geom = "ideo", fill = "gray70", radius = 7, trackWidth = 3) +
  layout_circle(gr, geom = "bar", radius = 10, trackWidth = 4, aes(fill = score, y = score)) +
  layout_circle(gr, geom = "point", color = "red", radius = 14, trackWidth = 3, grid = TRUE, aes(y = score)) +
  layout_circle(gr, geom = "link", linked.to = "to.gr", radius = 6, trackWidth = 1)
#### - I. Loading data ----
setwd("/home/ntmo83/Documents/Test_circos")
Vaestu12_012<-read.csv("Genome_methyl_pos.csv",sep=',',header = T)

circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02,
                            0))
circos.initialize(factors=c("contig_1","contig_2","contig_3","contig_4"), 
                  xlim=matrix(c(rep(0, 4), c(3013191,123120,17179,1121967)), ncol=2))
# genomes
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr)
})

circos.clear()
col_text <- "grey40"
circos.par("track.height"=0.8, gap.degree=5, cell.padding=c(0, 0, 0, 0))
circos.initialize(factors= c(3013191,123120,17179,1121967), xlim=matrix(c(rep(0, 4), c(3013191,123120,17179,1121967)), ncol=2))

# genomes
circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim), mean(ylim), chr, cex=0.5, col=col_text, 
              facing="bending.inside", niceFacing=TRUE)
}, bg.col="grey90", bg.border=F, track.height=0.06)

# genomes x axis
brk <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)*10^6
circos.track(track.index = get.current.track.index(), panel.fun=function(x, y) {
  circos.axis(h="top", major.at=brk, labels=round(brk/10^6, 1), labels.cex=0.4, 
              col=col_text, labels.col=col_text, lwd=0.7, labels.facing="clockwise")
}, bg.border=F)

# meth content
circos.track(factors=Vaestu12_012$seqid, x=Vaestu12_012$position, y=Vaestu12_012$coverage, panel.fun=function(x, y) {
  circos.lines(x, y, col="grey50", lwd=0.6)
  circos.segments(x0=0, x1=max(c(3013191,123120,17179,1121967)), y0=0.3, y1=0.3, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(c(3013191,123120,17179,1121967)), y0=0.5, y1=0.5, lwd=0.6, lty="11", col="grey90")
  circos.segments(x0=0, x1=max(c(3013191,123120,17179,1121967)), y0=0.7, y1=0.7, lwd=0.6, lty="11", col="grey90")
}, ylim=Vaestu12_012$coverage, track.height=0.08, bg.border=F)
# gc y axis
circos.yaxis(at=c(0.3, 0.5, 0.7), labels.cex=0.25, lwd=0, tick.length=0, labels.col=col_text, col="#FFFFFF")
