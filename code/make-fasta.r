install.packages('stringdist')
library(stringdist)

nuc <- read.table("../tmp/nucleotides.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")

inf <- read.table("../tmp/seqinfo.txt", header = TRUE, stringsAsFactors = FALSE)
aa <- read.table("../tmp/aminoacids.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
nuc <- read.table("../tmp/nucleotides.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
align <- read.dna("../tmp/dnabin.txt")

#Name the inf data appropriately
names(inf) <- c("samp.loc", "monkid", "weeks", "pID", "f.id")

# Create a consensus sequence - this should have the most popular nucleotide at each position
# Dashes (indicating missing data) should not be counted
cons <- rep("-", ncol(nuc))
for(i in 1:(ncol(nuc)  )){
    tmp <- table(nuc[,i])
    nodash <- (tmp[names(tmp) != "-"])
    cons[i] <- names(nodash[which.max(nodash)])
}

#Determine the distance between the consensus sequence and all sequences
alldists <- stringdist(apply(nuc, 1, paste, collapse = ""), paste(cons, collapse = ""))

#Exclude the sequences that are furthest away
exclude <- which(alldists > 190)
inf <- inf[-exclude,]
aa <-  aa[-exclude, ]
nuc <- nuc[-exclude,]
align <- align[-exclude, ]
attach(inf)

#Zandrea update, 8/29/16
# Week 39 was mislabeled and should be week 29 instead
weeks[weeks == 39] <- 29
inf$weeks <- weeks

infToWrite <- inf[, c(2, 1, 3:5) ]

seqs <- apply(nuc, 1, paste, collapse = "")
seqnames <- apply(infToWrite, 1, paste, collapse = "-")

write(paste(paste(">", seqnames, "\n", seqs, sep = ""), collapse = "\n\n"), "~/Desktop/fastatmp.fa")

setdiff(1:3312, c(RNAinds.n, plasmainds, DNAinds.n))

RNAinds <- c(grep("RNA-[0-9]+", seqnames), grep("PLASMA-[0-9]+", seqnames))
DNAinds <- grep("DNA-", seqnames)


write(paste(paste(">", seqnames[RNAinds], "\n", seqs[RNAinds], sep = ""), collapse = "\n\n"), "~/Desktop/macaqueRNA.fa")
write(paste(paste(">", seqnames[DNAinds], "\n", seqs[DNAinds], sep = ""), collapse = "\n\n"), "~/Desktop/macaqueDNA.fa")



require(foreach)
require(tidyverse)

a <- read.dna("~/Desktop/macaqueRNA.fa", format = "fasta")

nuc.new <- foreach(i = 1:nrow(a), .combine = 'rbind') %do% {
    toupper(paste(a[i,]))
}

nuc.new
rnanuc <- nuc[RNAinds,]

matches <- foreach(i = 1:nrow(nuc.new)) %do% {
    sum(rnanuc[i,] == nuc.new[i,])
}

table(unlist(matches))

#Ok, perfect matches, so the fasta is keeping all info for the RNA
#Let's doublecheck this with DNA

a <- read.dna("~/Desktop/macaqueDNA.fa", format = "fasta")

nuc.new <- foreach(i = 1:nrow(a), .combine = 'rbind') %do% {
    toupper(paste(a[i,]))
}

nuc.new
dnanuc <- nuc[DNAinds,]

matches <- foreach(i = 1:nrow(nuc.new)) %do% {
    sum(dnanuc[i,] == nuc.new[i,])
}

table(unlist(matches))





RNA <- read.dna("~/Desktop/macaqueRNA.fa", format = "fasta")
DNA <- read.dna("~/Desktop/macaqueDNA.fa", format = "fasta")

seqnames <- c(rownames(RNA), rownames(DNA))
rnanuc <- foreach(i = 1:nrow(RNA), .combine = 'rbind') %do% {
    toupper(paste(RNA[i,]))
}
dnanuc <- foreach(i = 1:nrow(DNA), .combine = 'rbind') %do% {
    toupper(paste(DNA[i,]))
}

bothnuc <- rbind(rnanuc, dnanuc)
colnames(bothnuc) <- paste("nuc", 1:900, sep = "")


aas <- matrix(data = NA, ncol = 299, nrow = nrow(bothnuc))
for(i in 1:nrow(bothnuc)){
    aas[i,] <- (translate(bothnuc[i,1+3:(900-1)]))
}
colnames(aas) <- paste("AA", 1:299, sep = "")

aaordered <- aa[c(RNAinds, DNAinds), ]

tmp <- foreach(i = 1:nrow(aaordered)) %do% {
    sum(aaordered[i,] == aas[i,])
}

table(unlist(tmp))

#Cool. 

infnew <- foreach(nameval = seqnames, .combine = 'rbind') %do% {
    strsplit(nameval, "-")[[1]][c(2, 1, 3:5)]
}

infordered <- inf[c(RNAinds, DNAinds), ]

matches <- foreach(i = 1:nrow(infordered)) %do% {
    sum(infordered[i,] == infnew[i,])
}

table(unlist(matches))

colnames(infnew) <- c("samp.loc", "monkid", "weeks", "pID", "f.id")
#Cool

rownames(bothnuc) <- NULL
rownames(aas) <- NULL

write.table(infnew, "../tmp/seqinfo.txt")
write.table(aas, "../tmp/aminoacids.txt")
write.table(bothnuc, "../tmp/nucleotides.txt")

inf <- read.table("../tmp/seqinfo.txt", header = TRUE, stringsAsFactors = FALSE)
aa <- read.table("../tmp/aminoacids.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
nuc <- read.table("../tmp/nucleotides.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")


