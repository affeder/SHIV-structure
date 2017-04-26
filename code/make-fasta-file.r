#DEV ONLY FILE (i.e., not to be uploaded)
require(stringdist)

#Read in the data
inf <- read.table("../dat-intermed/seqinfo.txt", header = TRUE, stringsAsFactors = FALSE)
aa <- read.table("../dat-intermed/aminoacids.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
nuc <- read.table("../dat-intermed/nucleotides.txt", header = TRUE, stringsAsFactors = FALSE, colClasses = "character")

# Create a consensus sequence - this should have the most popular nucleotide at each position
# Dashes (indicating missing data) should not be counted
cons <- rep("-", ncol(nuc))
for(i in 1:(ncol(nuc) )){
    tmp <- table(nuc[,i])
    nodash <- (tmp[names(tmp) != "-"])
    cons[i] <- names(nodash[which.max(nodash)])
}

#Determine the distance between the consensus sequence and all sequences
alldists <- stringdist(apply(nuc, 1, paste, collapse = ""), paste(cons, collapse = ""))

#Name the inf data appropriately


#Exclude the sequences that are furthest away
exclude <- which(alldists > 190)
inf <- inf[-exclude,]
aa <-  aa[-exclude, ]
nuc <- nuc[-exclude,]

names(inf) <- c("samp.loc", "monkid", "weeks", "pID", "f.id")
attach(inf)

#Zandrea update, 8/29/16
# Week 39 was mislabeled and should be week 29 instead
weeks[weeks == 39] <- 29
inf$weeks <- weeks

#Now, this is the data we want to write
infToWrite <- inf[, c(2, 1, 3:5) ]

seqs <- apply(nuc, 1, paste, collapse = "")
seqnames <- apply(infToWrite, 1, paste, collapse = "-")

RNAinds <- c(grep("RNA-[0-9]+", seqnames), grep("PLASMA-[0-9]+", seqnames))
DNAinds <- grep("DNA-", seqnames)

#Write two fastas
write(paste(paste(">", seqnames[RNAinds], "\n", seqs[RNAinds], sep = ""), collapse = "\n\n"), "../dat/RT-SHIV-RNA.fa")
write(paste(paste(">", seqnames[DNAinds], "\n", seqs[DNAinds], sep = ""), collapse = "\n\n"), "../dat/RT-SHIV-DNA.fa")


haps <- apply(nuc[, 135:900], 1, paste, collapse = "")
allDists <- stringdistmatrix(haps)
distMat <- allDists
save(distMat, file = "distmat")
