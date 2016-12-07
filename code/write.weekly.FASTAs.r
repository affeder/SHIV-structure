#Let's also generate fasta files for combinations that we've observed?
labs <- paste(samp.loc, weekind, monkid, sep = "_")

#Function writes a fasta file with the indices provided and writes it to file.out
write.my.fasta <- function(inds, file.out){
    outfile <- file(description = file.out, open = "w")
    for(i in 1:length(inds)){
        writeLines(paste(">",labs[inds[i]], "_", i, sep = ""), outfile)
        writeLines(paste(nuc[inds[i], ], collapse = ""), outfile)
        writeLines("\n", outfile)
    }
    close(outfile)
}

monks <- c("T98133", "A99039", "A99165", "A01198")
RNAcomps <- c("PLASMA", "GUTRNA", "LNRNA", "VAGRNA", "PBMCRNA")
for(j in monks){
    for(i in 1:5){
        for(k in 1:(length(RNAcomps) - 1)){
            lab1 <- paste(RNAcomps[k], i, j, sep = "_")
            matchlab1 <- which(labs == lab1)
            for(l in (k+1):(length(RNAcomps))){
                lab2 <- paste(RNAcomps[l], i, j, sep = "_")
                relDir <- paste("../tmp/paired_fastas/", RNAcomps[k],"-",
                                RNAcomps[l],"/", sep = "")
                if(!dir.exists(relDir)){ dir.create(relDir) }
                matchlab2 <- which(labs == lab2)
                if(length(matchlab1) > 2 & length(matchlab2) > 2){
                    inds <- c(matchlab1, matchlab2)
                    write.my.fasta(inds, paste(relDir, j,".", i, ".fa",sep = ""))
                }          
            }
        }
    }
}


#table(as.character(align[,605]))
#table(as.character(align[,1675]))

compartments = c("PLASMA", "GUT", "LN", "VAG", "PBMC");
compartments.withRNA = c("PLASMA", "GUTRNA", "LNRNA", "VAGRNA", "PBMCRNA");
monkeys = c("T98133", "A99165", "A99039", "A01198");
weeknames <- c("12|13", "15|16", "20|21", "26|27", "29|38|39|40|44")
for(j in monknames){
    mymonkname <- j
     for(i in 1:length(weeknames)){
         for(k in 1:(length(compartments) - 1)){
             for(l in (k+1):(length(compartments))){
                 print(paste(j, i, compartments[k], compartments[l], sep = "-"))
                 relDir <- paste("../tmp/weekly_trees/auto/",
                                 compartments[k],"-",compartments[l],"/", sep = "")
                 outDir <- paste("../tmp/hyphyout/auto/",
                                 compartments[k],"-",compartments[l],"/", sep = "")
                 if(!dir.exists(relDir)){ dir.create(relDir) }
                 if(!dir.exists(outDir)){ dir.create(outDir) }
                 monkWeek <- intersect(which(monkid == mymonkname),
                                       grep(weeknames[i], weeks))
                 comp1 <- intersect(monkWeek,
                                    grep(paste(compartments.withRNA[k]), samp.loc))
                 comp2 <- intersect(monkWeek,
                                    grep(paste(compartments.withRNA[l]), samp.loc))
                 if(length(comp1) > 0 & length(comp2) > 0){
                     tmpinds <- c(comp1, comp2)
                     tmp <- align[tmpinds,605:1675]
                     sto <- hclust(dist.dna(tmp))
                     sto <- as.phylo(sto)
                     labs <- paste(samp.loc[tmpinds])
                     sto$tip.label <- labs
                     write.tree(as.phylo(sto), paste(relDir, j, ".", i, sep = ""))
                 }
             }
        }
    }
}

printDat <- matrix(data = NA, nrow = 200, ncol = 6)
ind <- 0
for(j in monknames){
    for(k in 1:5){
        for(i in 1:(length(compartments) - 1)){
            for(l in (i+1):length(compartments)){
                ind <- ind + 1
                rightMonkTime <- intersect(which(monkid == j),
                                           grep(weeknames[k], weeks))
                rightComp1 <- length(intersect(
                    rightMonkTime, which(samp.loc == compartments.withRNA[i])))
                rightComp2 <- length(intersect(
                    rightMonkTime, which(samp.loc == compartments.withRNA[l])))
                printDat[ind, ] <-
                    c(j, k, compartments[i], compartments[l], rightComp1, rightComp2)
            }
        }
    }
}

write.table(printDat, "../tmp/pairwise.details.txt",
            col.names = F, row.names = F, quote = FALSE, sep = "\t")


