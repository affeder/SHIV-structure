dirsToCheck <- c("../tmp/weekly_trees/", "../tmp/hyphyout/", "../tmp/weekly_trees/auto/", "../tmp/hyphyout/auto")
for(dir in dirsToCheck){
    if(!dir.exists(dir)){ dir.create(dir) }
}

#In paper
#numReps <- 100

#For testing
numReps <- 10
set.seed(2039120)
compartments = c("PLASMA", "GUT", "LN", "VAG", "PBMC");
compartments.withRNA = c("PLASMA", "GUTRNA", "LNRNA", "VAGRNA", "PBMCRNA");
weeknames <- c("12|13", "15|16", "20|21", "26|27", "29|38|39|40|44")


for(j in monknames){
     for(i in 1:length(weeknames)){
         for(k in 1:(length(compartments))){
             for(l in (k):(length(compartments))){
                 print(paste(j, i, compartments[k], compartments[l], sep = "-"))
                 relDir <- paste("../tmp/weekly_trees/auto/",
                                 compartments[k],"-",compartments[l],"/", sep = "")
                 outDir <- paste("../tmp/hyphyout/auto/",
                                 compartments[k],"-",compartments[l],"/", sep = "")
                 if(!dir.exists(relDir)){ dir.create(relDir) }
                 if(!dir.exists(outDir)){ dir.create(outDir) }
                 monkWeek <- intersect(which(monkid == j),
                                       grep(weeknames[i], weeks))
                 comp1 <- intersect(monkWeek,
                                    grep(paste(compartments.withRNA[k]), samp.loc))
                 comp2 <- intersect(monkWeek,
                                    grep(paste(compartments.withRNA[l]), samp.loc))
                 samps1 <- apply(nuc[comp1,135:900], 1, paste, collapse = "")
                 s1 <- as.DNAbin(lapply(unique(samps1), function(x){
                     strsplit(x, split = "")[[1]]}))
                 samps2 <- apply(nuc[comp2,135:900], 1, paste, collapse = "")
                 s2 <- as.DNAbin(lapply(unique(samps2), function(x){
                     strsplit(x, split = "")[[1]]}))
                 for(downSampNum in c(3, 7, 10, 15, 20)){
                     newRelDir <- paste(relDir, "n", downSampNum, "/",sep = "")
                     if(!dir.exists(newRelDir)){ dir.create(newRelDir) }
                     if(length(s1) >= downSampNum & length(s2) >= downSampNum){
                         for(m in 1:numReps){
                             inds1 <- sample(1:length(s1), downSampNum, replace = FALSE)
                             inds2 <- sample(1:length(s2), downSampNum, replace = FALSE)
                             labsToAdd <- c(rep(compartments[k], downSampNum),
                                            rep(compartments[l], downSampNum))
                             if(k == l){
                                 #If the names are the same, we want to give a fake other label 
                                 # so as not to confuse hyphy
                                 labsToAdd[1:downSampNum] <- rep("OTHER", downSampNum)
                             }
                             tmp <- c(s1[inds1], s2[inds2])
                             sto <- hclust(dist.dna(tmp))
                             sto <- as.phylo(sto)
                             labs <- labsToAdd
                             sto$tip.label <- labs
#                             ggtree(sto) + geom_tiplab()
                             write.tree(as.phylo(sto), paste(newRelDir, j, ".", i,".",m, sep = ""))
                         }
                     }
                 }
             }
        }
    }
}

