#This file converts fasta files into the format usable for the analysis
require(ape)
require(foreach)
require(seqinr)
require(stringdist)

RNA <- read.dna("../dat/RT-SHIV-RNA.fa", format = "fasta")
DNA <- read.dna("../dat/RT-SHIV-DNA.fa", format = "fasta")

seqnames <- c(rownames(RNA), rownames(DNA))
rnanuc <- foreach(i = 1:nrow(RNA), .combine = 'rbind') %do% {
    toupper(paste(RNA[i,]))
}
dnanuc <- foreach(i = 1:nrow(DNA), .combine = 'rbind') %do% {
    toupper(paste(DNA[i,]))
}

bothnuc <- rbind(rnanuc, dnanuc)
colnames(bothnuc) <- paste("nuc", 1:900, sep = "")
rownames(bothnuc) <- NULL

aas <- matrix(data = NA, ncol = 299, nrow = nrow(bothnuc))
for(i in 1:nrow(bothnuc)){
    aas[i,] <- (translate(bothnuc[i,1+3:(900-1)]))
}
colnames(aas) <- paste("AA", 1:299, sep = "")
rownames(aas) <- NULL

infnew <- foreach(nameval = seqnames, .combine = 'rbind') %do% {
    strsplit(nameval, "-")[[1]][c(2, 1, 3:5)]
}
colnames(infnew) <- c("samp.loc", "monkid", "weeks", "pID", "f.id")
rownames(bothnuc) <- NULL

write.table(infnew, "../tmp/seqinfo.txt")
write.table(aas, "../tmp/aminoacids.txt")
write.table(bothnuc, "../tmp/nucleotides.txt")

#This is slow (hours)
haps <- apply(bothnuc[, 135:900], 1, paste, collapse = "")
allDists <- stringdistmatrix(haps)
distMat <- allDists
save(distMat, file = "../tmp/distmat")


