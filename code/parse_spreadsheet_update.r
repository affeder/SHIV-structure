#DEV ONLY FILE (i.e., not to be uploaded)
library(Biostrings)
library(seqinr)
library(muscle)
library(stringr)
require(ape)

dat <- readDNAStringSet("../dat/0906.aligned.fasta", format = "fasta")
dat <- dat[1:10]
#Note, muscle is slow (hours)
aln <- muscle::muscle(dat)

seqs <- c()
ids <- c()
for(i in 1:nrow(aln)){
    seqs[i] <- as.character((aln@unmasked[i]))
    ids[i] <- names(aln@unmasked[i])
}


#which index holds our consensus sequence?
conind <- which(regexpr("CON", ids) > 0)
ids[conind]
seqs[conind]

tmp <- strsplit(ids, "\\.")
monkid <- c()
num <- c()

monkid <- rep(0, length(ids))
potentialids <- c("T98133", "A99039", "A01198", "A99165", "A99030")

for(i in 1:length(potentialids)){
    monkid[which(regexpr(potentialids[i], ids, ignore.case = TRUE) > 0)] <- potentialids[i]
}
table(monkid)


#Now, we need to parse the rest of it
samp.loc <- rep(0,  length(ids))

#where is plasma present?
samp.loc[regexpr("plas", ids, ignore.case = T) > 0] <- "PLASMA"
samp.loc[regexpr("vag.rna", ids, ignore.case = T) > 0] <- "VAGRNA"
samp.loc[regexpr("vagrna", ids, ignore.case = T) > 0] <- "VAGRNA"
samp.loc[regexpr("vag.dna", ids, ignore.case = T) > 0] <- "VAGDNA"
samp.loc[regexpr("vagdna", ids, ignore.case = T) > 0] <- "VAGDNA"
samp.loc[regexpr("pbmc.rna", ids, ignore.case = T) > 0] <- "PBMCRNA"
samp.loc[regexpr("pbmc.DNA", ids, ignore.case = T) > 0] <- "PBMCDNA"
samp.loc[regexpr("pbmcrna", ids, ignore.case = T) > 0] <- "PBMCRNA"
samp.loc[regexpr("pbmcdna", ids, ignore.case = T) > 0] <- "PBMCDNA"
samp.loc[regexpr("ln.dna", ids, ignore.case = T) > 0] <- "LNDNA"
samp.loc[regexpr("lndna", ids, ignore.case = T) > 0] <- "LNDNA"
samp.loc[regexpr("ln.rna", ids, ignore.case = T) > 0] <- "LNRNA"
samp.loc[regexpr("lnrna", ids, ignore.case = T) > 0] <- "LNRNA"
samp.loc[regexpr("gut.dna", ids, ignore.case = T) > 0] <- "GUTDNA"
samp.loc[regexpr("gutdna", ids, ignore.case = T) > 0] <- "GUTDNA"
samp.loc[regexpr("gut.rna", ids, ignore.case = T) > 0] <- "GUTRNA"
samp.loc[regexpr("gutrna", ids, ignore.case = T) > 0] <- "GUTRNA"
samp.loc[regexpr("ILIUM", ids, ignore.case = T) > 0] <- "ILIUM"

#strains from A0 having PBMC
table(samp.loc[intersect(which(monkid == "A01198"), which(regexpr("PBMC", samp.loc) > 0))])

ids[which(samp.loc == 0)]

#update the ILIUM samp.loc
samp.loc[samp.loc == "ILIUM"] <- "GUT"
samp.loc[samp.loc == "GUT"] <- "ILIUM"

#There appears to be some sort of id starting with W
#Ah! That will be the week that it's from
weeks <- rep(0, length =  length(ids))

wpID <- unlist(str_extract_all(ids, "[Ww][0-9]+"))
weeknumbs <- names(table(substr(wpID, 2, nchar(wpID))))

weeks.to.check <- weeknumbs[-which(nchar(weeknumbs) > 2)]

for( i in weeks.to.check){
    weeks[regexpr(paste("w",i, sep = ""), ids, ignore.case = T) > 0] <- i 
}
weeks <- as.numeric(weeks)


#ok, apart from the consensus, these are all the ilium samples
ids[which(weeks == 0)]

#They're from monkeys A01198 and A99165
table(weeks[monkid == "A99165"])
#max week = 39
#This should actually be 29
weeks[intersect(which(weeks == 0), which(monkid == "A99165"))] <- 39

table(weeks[monkid == "A01198"])
#max week = 44
weeks[intersect(which(weeks == 0), which(monkid == "A01198"))] <- 44

table(weeks)
#good


#Alright, let's try to get some better resolution on some of this sampling stuff.
toRemove <- c()

table(monkid)
#Let's go monkey by monkey

#T98133 should have 860 sequences
#For some reason, we have 869
table(samp.loc[monkid == "T98133"])

#It looks like the vagina and the gutdna match the numbers that are expected

#Ok, it looks like there's a duplicate in the gutrna
length(unique(ids[intersect(which(monkid == "T98133"), which(samp.loc == "GUTRNA"))]))
length(ids[intersect(which(monkid == "T98133"), which(samp.loc == "GUTRNA"))])

tab <- table(ids[intersect(which(monkid == "T98133"), which(samp.loc == "GUTRNA"))])
repeated <- names(which(tab > 1))
toRemove <- c(toRemove, which(ids == repeated)[1])

#There's also a LN RNA that appears twice- “9.T98133.LN.RNA.W13p7a17”
toRemove <- c(toRemove, grep("W13p7a17", ids)[1])

#Ok, now the gut should be ok.

#It also looks like there's a duplicate in the plasma
tab <- table(ids[intersect(which(monkid == "T98133"), which(samp.loc == "PLASMA"))])
repeated <- names(which(tab > 1))

for(rep in repeated){
    toRemove <- c(toRemove, which(ids == rep)[1])
}

#Ok, so, we've identified our plasma duplicates

#We now need to check up on our PBMC
#We have one extra PBMC RNA - let's find the duplicate
tab <- table(ids[intersect(which(monkid == "T98133"), which(samp.loc == "PBMCDNA"))])
repeated <- names(which(tab > 1))
#Hmm, there are unique IDs

#This is the distribution of the monkeys across weeks. 
table(weeks[intersect(which(monkid == "T98133"), which(samp.loc == "PBMCDNA"))])

#There seems to be one extra in week 16.
sort(ids[intersect(intersect(which(monkid == "T98133"), which(samp.loc == "PBMCDNA")), which(weeks == 16))])

#These ids all check out with the spreadsheet, with the exception being p18e5, which I think is P1835
#We'll keep all of these...

#Alright! Now, let's move on to LN DNA.


#Alright. what are all these 0s?
ids[intersect(which(monkid == "T98133"), which(samp.loc == "0"))]
ids0 <- intersect(which(monkid == "T98133"), which(samp.loc == "0"))
table(samp.loc[monkid == "T98133"])

#There are four PBMC sequences
#t98133.pbmc.w13p32a4
#t98133.pbmc.w13p32a17
#t98133.pbmc.w13p28g6
#t98133.pbmc.w16p32c7

#Alright, so, according to the spreadsheet, all the p33s/p31s are LN DNAs
tmp <- ids[ids0]
tmp

length(tmp[regexpr("T98133.LN.W[0-9]*p3", ids[ids0], ignore.case = T) > 0])
length(tmp)
#The difference between these two is 4, which is the number of pbmc questionable strains
#This means that all of these are LN DNAs
sum(regexpr("T98133.LN.W[0-9]*p3", ids, ignore.case = T) > 0)

samp.loc[regexpr("T98133.LN.W[0-9]*p3", ids, ignore.case = T) > 0] <- "LNDNA"

table(samp.loc[monkid == "T98133"], weeks[monkid == "T98133"])

#Still some suspicious 0s left
ids0 <- intersect(which(monkid == "T98133"), which(samp.loc == "0"))
ids[ids0]

#And finally, the LNDNA/RNA samples marked as week 22 need to be week 20 instead
length(intersect(intersect(which(monkid == "T98133"), which(samp.loc == "LNDNA")), which(weeks == 22)))
#ok, we need to change this one
weeks[intersect(intersect(which(monkid == "T98133"), which(samp.loc == "LNDNA")), which(weeks == 22))] <- 20

intersect(intersect(which(monkid == "T98133"), which(samp.loc == "LNRNA")), which(weeks == 22))
#Ok, nothing to change for this one

#Ok, now this looks good except for the things that need to be removed in "toRemove"


#Let's think about a differnent monkey - A99039

ids[intersect(which(monkid == "A99039"), which(samp.loc == "0"))]
table(samp.loc[which(monkid == "A99039")])
#Gut DNA - one extra
#Gut RNA - missing 13
#LN DNA - good
#LN RNA - good
#PBMCDNA - 115 (should be 170)
#PBMCRNA - missing 1
#Plasma - good
#Vag DNA - one extra
#Vag RNA - good


#Maybe some of the sequences that we are missing have that messed up monkid?
table(monkid)
#Looks promising

ids[monkid == "A99030"]
#So, ideally, we'd have 13 + 55 + 1 = 69 sequences
sum(monkid == "A99030")
#We have 68. Pretty good.
ids0 <- which(monkid == "A99030")
tmp <- ids[ids0]

#Which indivs are gut RNA?
sum(regexpr("a99030.gut.*", ids, ignore.case = T) > 0)
samp.loc[regexpr("a99030.gut.*", ids, ignore.case = T) > 0] <- "GUTRNA"

#Which indivs are PBMC?
#Since the individuals line up here, these must all be the DNA, so the one 0 might be the missing PBMCRNA
sum(regexpr("a99030.pbmc.*", ids, ignore.case = T) > 0)
samp.loc[regexpr("a99030.pbmc.*", ids, ignore.case = T) > 0] <- "PBMCDNA"
samp.loc[which(ids == "9.A99039PBMCW20P14A7")] <- "PBMCRNA"
ids[intersect(which(samp.loc == "0"), which(monkid == "A99039"))]

#Look at the doubles - gut DNA, Vag DNA
tab <- table(ids[intersect(which(monkid == "A99039"), which(samp.loc == "GUTDNA"))] )
repeated <- names(which(tab > 1))
toRemove <- c(toRemove, which(ids == repeated)[1])

tab <- table(ids[intersect(which(monkid == "A99039"), which(samp.loc == "VAGDNA"))] )
repeated <- names(which(tab > 1))
#Alright, so, again, there appears to be no repeated stuff

sort(ids[intersect(intersect(which(monkid == "A99039"), which(samp.loc == "VAGDNA")) , which(weeks == "15"))])

#Ok, we need to switch the A99030 to A99039
monkid[monkid == "A99030"] <- "A99039"

#Ok, I think everything is accounted for in A99039
table(samp.loc[monkid == "A99039"])

#Let's look at 
table(samp.loc[monkid == "A99165"])

#GutDNA 149 (152) -3
#GutRNA 107 (136) -29
#lnDNA 63 (145)
#lnRNA 112 (111) +1
#pbmcDNA 31 (143)
#pbmcRNA 34 (34) Good
#Plasma 132 (124)
#VagDNA 69 (64)
#VagRNA 4 (4) Good

#Ok, so, we are missing 32 gut samples, and we have 32 Ilium samples
ids[intersect(which(samp.loc == "ILIUM"), which(monkid == "A99165"))]
#Which of these are DNA and which are RNA?
ids[intersect(which(samp.loc == "GUTRNA"), which(monkid == "A99165"))]

#ids[intersect(which(weeks == "GUTRNA"), which(monkid == "A99165"))]
table(weeks[which(monkid == "A99165")])
table(samp.loc[intersect(which(monkid == "A99165"), which(weeks == 39))])
#Hm, weird
#It seems like we're missing 4 from GutDNA
#It seems like we have 4 extra in VagDNA?
ids[intersect(intersect(which(monkid == "A99165"), which(weeks == 39)), which(samp.loc == "VAGDNA"))]
#Nope, it looks like 4 GutDNA are here that are maybe repeated?

#What 4 are we missing from GutDNA?
sort(ids[intersect(intersect(which(monkid == "A99165"), which(weeks == 39)), which(samp.loc == "GUTDNA"))])

#We're missing these sequences
#P15I13
#P15I18
#P15I19
#P15I23
#For each of these sequences, we have two sequences, and they're all listed as 
for(i in c("P15I13", "P15I18", "P15I19", "P15I23" )){
    tmpids <- which(regexpr(i, ids) > 0)
    samp.loc[tmpids] <- c("GUTDNA", "GUTDNA")
    toRemove <- c(toRemove, tmpids[1])
}
toRemove

table(samp.loc[monkid == "A99165"], weeks[monkid == "A99165"])
#Gut DNA - week 20, extra sequence?
tab <- table(ids[intersect(intersect(which(monkid == "A99165"), which(samp.loc == "GUTDNA")), which(weeks == 20))] )
repeated <- names(which(tab > 1))
toRemove <- c(toRemove, which(ids == repeated)[1])

table(samp.loc[monkid == "A99165"], weeks[monkid == "A99165"])
#week 20, has 29, two repeats?
#week 26, has 29, one repeats?
#week 39, has 0, should have 32

tab <- table(ids[intersect(intersect(which(monkid == "A99165"), which(samp.loc == "GUTRNA")), which(weeks == 20))] )
repeated <- names(which(tab > 1))
#2 repeats, ok!
toRemove <- c(toRemove, which(ids == repeated[1])[1])
toRemove <- c(toRemove, which(ids == repeated[2])[1])

tab <- table(ids[intersect(intersect(which(monkid == "A99165"), which(samp.loc == "GUTRNA")), which(weeks == 26))] )
repeated <- names(which(tab > 1))
#1 repeat, ok!
toRemove <- c(toRemove, which(ids == repeated)[1])

tab <- table(ids[intersect(intersect(which(monkid == "A99165"), which(samp.loc == "GUTRNA")), which(weeks == 39))] )
#Nothing - Hm, ok. I bet these are our ILIUM samples
#There should be 32
#gee
intersect(which(samp.loc == "ILIUM"), which(monkid == "A99165"))
#32 - awesome!
samp.loc[intersect(which(samp.loc == "ILIUM"), which(monkid == "A99165"))] <- "GUTRNA"

#Alright - Gut DNA/RNA is fine
#LN?
table(samp.loc[monkid == "A99165"], weeks[monkid == "A99165"])
ids[which(regexpr("P21", ids) > 0)]

id.tmp <- intersect(which(samp.loc == "0"), which(monkid == "A99165"))
ids[regexpr("LN", ids[id.tmp]) > 0]
foo <- ids[id.tmp][regexpr("LN", ids[id.tmp]) > 0]
sum(regexpr("p23", foo) > 0) + sum(regexpr("p20", foo) > 0) + sum(regexpr("p24", foo) > 0) + sum(regexpr("p25", foo) > 0)
length(foo)
#Ok, so, all of the p23s, p20s, p24s, p25s are LN DNA -

ids[id.tmp[(regexpr("p23[a-zA-Z]+", ids[id.tmp]) > 0)]]
samp.loc[id.tmp[(regexpr("p23[a-zA-Z]+", ids[id.tmp]) > 0)]] <- "LNDNA"
ids[id.tmp[(regexpr("LN.+p20[a-zA-Z]+", ids[id.tmp]) > 0)]]
samp.loc[id.tmp[(regexpr("LN.+p20[a-zA-Z]+", ids[id.tmp]) > 0)]] <- "LNDNA"
ids[id.tmp[(regexpr("p24[a-zA-Z]+", ids[id.tmp]) > 0)]]
samp.loc[id.tmp[(regexpr("p24[a-zA-Z]+", ids[id.tmp]) > 0)]] <- "LNDNA"
ids[id.tmp[(regexpr("LN.+p25[a-zA-Z]+", ids[id.tmp]) > 0)]]
samp.loc[id.tmp[(regexpr("LN.+p25[a-zA-Z]+", ids[id.tmp]) > 0)]] <- "LNDNA"

table(samp.loc[monkid == "A99165"], weeks[monkid == "A99165"])

#LN RNA seems ok

tab <- table(ids[intersect(intersect(which(monkid == "A99165"), which(weeks == 13)), which(samp.loc== "LNRNA"))])
repeated <- names(which(tab > 1))

#No repeats
sort(names(tab))

table(samp.loc[monkid == "A99165"], weeks[monkid == "A99165"])
#Ok, now we need to deal with PBMC

id.tmp <- intersect(which(samp.loc == "0"), which(monkid == "A99165"))
ids[id.tmp]

ids[regexpr("pbmc", ids[id.tmp]) > 0]
foo <- ids[id.tmp][regexpr("pbmc.+p2.+", ids[id.tmp]) > 0]
foo
sum(regexpr("p23", foo) > 0) + sum(regexpr("p20", foo) > 0) + sum(regexpr("p24", foo) > 0) + sum(regexpr("p25", foo) > 0)
length(foo)
#Ok, so, all of the p23s, p20s, p24s, p25s are LN DNA -

ids[id.tmp[(regexpr("pbmc.+p2.+", ids[id.tmp]) > 0)]]
samp.loc[id.tmp[(regexpr("pbmc.+p2.+", ids[id.tmp]) > 0)]] <- "PBMCDNA"

tab <- table(ids[intersect(intersect(which(monkid == "A99165"), which(weeks == 15)), which(samp.loc== "PLASMA"))])
repeated <- names(which(tab > 1))
repeated
#No repeats...
#There are 28, but there should actually be 20
sort(names(tab))



tab <- table(ids[intersect(intersect(which(monkid == "A99165"), which(weeks == 20)), which(samp.loc== "VAGDNA"))])
repeated <- names(which(tab > 1))
toRemove <- c(toRemove, which(ids == repeated)[1])
#
#  VAGDNA   5 23  0 28  4  0  1
# VagDNA week 20, 28 (one double) -> 27
# VagDNA week 39, missing 4?


tab <- table(ids[intersect(intersect(which(monkid == "A99165"), which(weeks == 39)), which(samp.loc== "VAGDNA"))])
tab



table(samp.loc[monkid == "A01198"], weeks[monkid == "A01198"])
#Ok, the first thing, is that there seem to be MANY doubles. 

#GutDNA?
for(k in c("GUTDNA", "GUTRNA", "PBMCDNA", "PBMCRNA", "LNDNA", "LNRNA", "PLASMA", "VAGDNA", "VAGRNA")){

    tab <- table(ids[intersect(which(monkid == "A01198"), which(samp.loc== k))])
    repeated <- names(which(tab > 1))

    for(i in 1:length(repeated)){
        inds <- which(ids == repeated[i])
        toRemove <- c(toRemove, inds[2:length(inds)])
    }
}

tab <- table(ids[intersect(which(monkid == "A01198"), which(samp.loc== "ILIUM"))])
repeated <- names(which(tab > 1))
for(i in 1:length(repeated)){
    inds <- which(ids == repeated[i])
    toRemove <- c(toRemove, inds[2:length(inds)])
}

table(samp.loc[monkid == "A01198"], weeks[monkid == "A01198"])

ids[intersect(which(monkid == "A01198"), which(samp.loc== "0"))]
#These are all LN, and the numbers match perfectly with the LN DNA, and all LN RNA is accounted for

samp.loc[intersect(which(monkid == "A01198"), which(samp.loc== "0"))] <- "LNDNA"

#Let's change the ILIUM
samp.loc[samp.loc == "ILIUM"] <- "GUTRNA"

#GutDNA
sort(ids[intersect(intersect(which(monkid == "A01198"), which(samp.loc== "GUTDNA")), which(weeks == 12))])

#GutCDNA = GutRNA
samp.loc[regexpr("GUTCDNA", ids) > 0] <- "GUTRNA"

#PBMC DNA -> LN DNA
sort(ids[intersect(intersect(which(monkid == "A01198"), which(samp.loc== "PBMCDNA")), which(weeks == 44))])
#All of the P7s are LN DNA.
tmp.ids <- intersect(intersect(which(monkid == "A01198"), which(samp.loc== "PBMCDNA")), which(weeks == 44))
samp.loc[tmp.ids[regexpr("P7", ids[tmp.ids]) > 0]] <- "LNDNA"

#VagDNA
tmp.ids <- intersect(intersect(which(monkid == "A01198"), which(samp.loc== "VAGDNA")), which(weeks == 12))
samp.loc[tmp.ids[regexpr("VAGcDNA", ids[tmp.ids]) > 0]] <- "VAGRNA"

#Ok, we still have an extra for GutDNA?
table(ids[intersect(intersect(which(monkid == "A01198"), which(samp.loc== "PBMCRNA")), which(weeks == 26))])

#this should be fixed
toRem <- toRemove[!is.na(toRemove)]

samp.loc[which(ids == "t98133.pbmc.w13p32a4")] <- "PBMCDNA"
samp.loc[which(ids == "t98133.pbmc.w13p32a17")] <- "PBMCDNA"
samp.loc[which(ids == "t98133.pbmc.w13p28g6")] <- "PBMCDNA"
samp.loc[which(ids == "t98133.pbmc.w16p32c7")] <- "PBMCDNA"

loc <- samp.loc[-toRem]
sampweek <- weeks[-toRem]
p.id <- wpID[-toRem]
f.id <- ids[-toRem]
monk.id <- monkid[-toRem]




table(monkid)
table(loc[monk.id == "A01198"], sampweek[monk.id == "A01198"])
table(loc[monk.id == "A99039"], sampweek[monk.id == "A99039"])
table(loc[monk.id == "A99165"], sampweek[monk.id == "A99165"])
table(loc[monk.id == "T98133"], sampweek[monk.id == "T98133"])


table(monk.id, loc)

#Ok, and then there is some sort of sequencing ID thing
otherID <- unique(toupper(unlist(str_extract_all(ids, "[Pp][0-9]+.*$"))))
pID <- rep(0, length = length(ids))
for( i in otherID){
    pID[regexpr(i, ids, ignore.case = T) > 0] <- i 
}

ids[which(pID == 0)]

testrelatedids <- unique(toupper(unlist(str_extract_all(toupper(ids), "TEST.*$"))))
for( i in testrelatedids){
    pID[regexpr(i, ids, ignore.case = T) > 0] <- i 
}
ids[which(pID == 0)]

seqs <- seqs[-toRem]
ids <- ids[-toRem]
seqs
length(seqs)
length(ids)

alignL <- (strsplit(seqs, split = ""))
align <- matrix(data = NA, nrow = length(alignL), ncol = length(alignL[[1]]))
for(i in 1:length(alignL)){
    align[i,] <- alignL[[i]]
}
conind <- grep("CON", ids)
consens <- as.character(align[conind,2:ncol(align)])
con <- consens[-which(consens == "-")]

#begin-paste
off1 <- translate(con[1:length(con)])
off2 <- translate(con[2:length(con)])
off3 <- translate(con[3:length(con)])


sto <- matrix(data = NA, ncol = 9, nrow = 300)
for(i in 1:300){
    sto[i,] <- (off3[c(1:9) + i])
}

sto[sto[,1] == "L" & sto[,2] == "K",]
#Alright. 
which(sto[,1] == "L" & sto[,2] == "K")
#AA 253 is our AA 100


#Hm, interesting
#Let's just make a one time map mapping each position from con (no -s) to consens (with -s)

reminds <- which(as.character(align[conind,]) == "-")

#763 == 300
width <- 900
nucs <- matrix(data = NA, ncol = width, nrow = nrow(align))
tmp <- align[, -reminds]
nucs <- tmp[,463:(462+width)]
dim(nucs)
#for(i in 1:nrow(align)){
#    if(i %% 50 == 0){ print(i)}
#    nucs[,i] <- as.character(tmp)[,463+i ]
#}
#translate(nucs[1,3:width])[c(100, 101, 103, 106)]
#Ok, cool.

length(translate(nucs[i,1+0:(width-2)]))


aas <- matrix(data = NA, ncol = 299, nrow = nrow(align))
for(i in 1:nrow(align)){
    aas[i,] <- (translate(nucs[i,1+3:(width-1)]))
}
aas[10,c(100, 101, 103, 106, 108)]
#translate(nucs[10,100*3 +1 + 0:20])
#translate(nucs[10,0*3 +1 + 0:500])
#Alright, good.

colnames(nucs) <- paste("nuc", 1:(ncol(nucs)), sep = "")
colnames(aas) <- paste("AA", 1:ncol(aas), sep = "")
toPrint <- cbind(loc, monk.id, sampweek, p.id, f.id)

dim(toPrint)
dim(nucs)
dim(aas)

                                        #end-paste


write.table(cbind(toPrint, nucs, aas), "../tmp/alldat.txt", row.names = FALSE, quote = FALSE)
write.table(toPrint, "../tmp/seqinfo.txt", row.names = FALSE, quote = FALSE)
write.table(aas, "../tmp/aminoacids.txt", row.names = FALSE, quote = TRUE)
write.table(nucs, "../tmp/nucleotides.txt", row.names = FALSE, quote = TRUE)
write.dna(align, "../tmp/dnabin.txt")



haps <- apply(nuc[, 135:900], 1, paste, collapse = "")
allDists <- stringdistmatrix(haps)
distMat <- allDists
save(distMat, file = "distmat")




