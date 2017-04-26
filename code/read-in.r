library(stringdist)
library(seqinr)
library(stringr)
library(ape)
require(xtable)
require(RColorBrewer)
require(sfsmisc)
require(TeachingDemos)
require(cluster)
require(foreach)
require(poppr)
require(doParallel)
require(tidyverse)

# Read in the data files generated in the format script
inf <- read.table("../tmp/seqinfo.txt", header = TRUE, stringsAsFactors = FALSE)
aa <- read.table("../tmp/aminoacids.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
nuc <- read.table("../tmp/nucleotides.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")

attach(inf)

#Create an ordered list of monkey names that will be used 
monknames <- c("T98133" , "A99165", "A99039", "A01198" )

#Create a list by id of which positions are polymorphic for each macaque
polyindslist <- list()
for(id in monknames){
    specinds <- which(monkid == id)
    poly <- c()
    minorallele <- 10
    for(i in 1:ncol(nuc)){
        tab <- table(nuc[specinds,i])
        tab.letters <- tab[names(tab) != "-"]
        poly[i] <- sum(tab.letters >= minorallele) > 1
    }
    polyinds <- which(poly == TRUE)
    polyindslist[[id]] <- polyinds
}

# The weekind variable will keep track of which timepoints are a week apart
weekopts <- c("12|13", "15|16", "20|21", "26|27", "29|38|40|44")
weekind <- rep(0, length(weeks))
for(i in 1:length(weekopts)){
    weekind[ grep(weekopts[i], weeks) ] <- i
}

# Create a consensus sequence - this should have the most popular nucleotide at each position
# Dashes (indicating missing data) should not be counted
cons <- rep("-", ncol(nuc))
for(i in 1:(ncol(nuc)  )){
    tmp <- table(nuc[,i])
    nodash <- (tmp[names(tmp) != "-"])
    if(length(nodash) == 0){ cons[i] <- "X" }else{
        cons[i] <- names(nodash[which.max(nodash)])
    }
}

#Read the precomputed pairwise distance matrix
#This is slow, so it's worth only doing once and then loading it
load("../tmp/distmat")
distMat <- as.matrix(distMat)

