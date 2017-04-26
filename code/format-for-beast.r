#Create fastas for BEAST


source("muller-plot-functions.r")
#source("shared-functions.r")
#source("read.in.data.r")

require(hash)

if(!dir.exists("../tmp/beast/")){ dir.create("../tmp/beast/") }
if(!dir.exists("../tmp/beast/fasta/")){ dir.create("../tmp/beast/fasta") }
if(!dir.exists("../tmp/beast/mcc/")){ dir.create("../tmp/beast/mcc") }
if(!dir.exists("../tmp/beast/xml/")){ dir.create("../tmp/beast/xml") }
if(!dir.exists("../tmp/beast/trees/")){ dir.create("../tmp/beast/trees") }

#Our trees are going to have annotation based on the color things appear in Muller
# plots, and the mutations that are present (which are of course related). The first 
# step in printing these fasta files for BEAST will be tto give each sequence
# a tip name that will allow the reconstruction of both of these things

#First, for each monkey, to determine which mutations are at frequency at least 1% in each macaque

returnColsAndMuts <- function(monk){

   #First, choose ALL variants that are going to appear on the plot
    tmpnuc <- nuc[which(monkid == monk),]
    #Any location that has a minor allele at frequency > 1% should be included
    maf <- nrow(tmpnuc)/100
    hapinds <- c()
    for(j in 1:ncol(tmpnuc)){
        tmp <- table(tmpnuc[,j])
        dashind <- which(names(tmp) == "-")
        if(length(dashind) > 0){
            tmp <-tmp[-dashind]
        }
        tmp <- sort(tmp, decreasing = TRUE)
        if(length(tmp) > 1){
            if(tmp[2] > maf){
                hapinds[length(hapinds)+1] <- j
            }
        }
    }
    #hapinds now has all of the variable haplotypes
    hapinds <<- hapinds
    haps <<- nuc[,hapinds]

    #haps are all of the haplotypes that are possible
    #Note, each monkey will use its own indices (that are relevent)
    indsToProcess <- which(monkid == monk)
    #Call the main function to determine the hierchy
    hierchinf <- setup.hierchy(indsToProcess)
    #Store the returned information 
    dat <- hierchinf$dat
    haps.over.time <- hierchinf$haps.over.time
    mut.relat <- hierchinf$mut.relat
    all.comp.inds <- hierchinf$all.comp.inds
    #Compute names (i.e., the list of mutations) for each of the sequences
    allNames <- c()
    for(i in 1:nrow(mut.relat)){
        if(sum(which(mut.relat[i,] == 1)) == 0){
            allNames[i] <- "WT"
        }else{
            allNames[i] <- labelmaker(names(which(mut.relat[i,] == 1)), all.comp.inds)
        }
    }
    #Assign colors for each sequence
    cols <- rep("white", nrow(haps.over.time))
    #These are our base colors palettes 
    colset <-  c(brewer.pal(8, "Set2")[1:6], brewer.pal(12, "Paired")[c(2, 6, 9, 10, 11)])
    #How many haplotypes are going to have "base colors" versus variations?
    maxuse <- 10
    #First, set some sort of base colors for the most frequent haplotypes
    popnames <- allNames[order(apply(haps.over.time, 1, sum), decreasing = TRUE)][1:maxuse]
    popinds <- c()
    for(i in popnames){
        cols[ which(allNames == i) ] <- colset[which(i == popnames)]
        popinds <- c(popinds, which(allNames == i))
    }
    #Now, go through the descendents and choose colors for each descendent
    for(i in 1:nrow(haps.over.time)){
        if(sum(i == popinds) == 0){
            inq <- matrix(rep(mut.relat[i,], length(popinds)),
                          nrow = length(popinds), byrow = TRUE)
            distvals <- (apply(mut.relat[popinds, ] - inq, 1, function(x){sum(abs(x))}))
            bestAnc <- which(distvals == min(distvals))
            if(length(bestAnc) > 1){
                mutlens <- apply(mut.relat[popinds,], 1, sum)
                bestAnc <- bestAnc[which(mutlens[bestAnc] == max(mutlens[bestAnc]))[1]]
            }
            cols[i] <- paste("anc", bestAnc, sep = "-")
        }
    }
    for(i in 1:maxuse){
        colsReassign <- which(cols == paste("anc", i, sep = "-"))
        ind <- 1
        compCols <- setdiff(1:maxuse, i)
        numPerComp <- ceiling(length(colsReassign)/7)
        for(j in compCols){
            colmake <- colorRampPalette(c(colset[i], colset[j]))
            tmpCols <- colmake( 2*(numPerComp + 1)+ 1)
            for(k in 1:numPerComp){
                cols[colsReassign[ind]] <- tmpCols[1 + k]
                ind <- ind + 1
            }
        }
    }

    keys <- rownames(haps.over.time)

    allNames <- gsub(",", "-", allNames)

    nameHash <- hash(keys = keys, values = allNames)
    colHash <-  hash(keys = keys, values = cols)

    idsToRet <- apply(nuc[indsToProcess,hapinds], 1, paste, collapse = "")

    returnNames <- c()
    returnCols <- c()
    for(i in 1:length(idsToRet)){
        returnNames[i] <- values(nameHash[idsToRet[i]])
        returnCols[i] <- values(colHash[idsToRet[i]])
    }

    return(list(cols = returnCols, muts = returnNames))
}

monknames <- c("T98133", "A99165", "A99039","A01198")
seqsToPrint <- apply(nuc[, 135:900], 1, paste, collapse = "")


for(monk in monknames){

    colAndNameInf <- returnColsAndMuts(monk)

    fastaLabels <- apply(cbind(
        inf[which(monkid == monk), c('samp.loc', 'monkid', 'weeks')],
        paste(which(monkid == monk)),
        colAndNameInf$muts, 
        colAndNameInf$cols), 1, paste, collapse = "_")

    rnaInds <- intersect(which(inf$monkid == monk), grep("RNA|PLASMA", inf$samp.loc))
    dnaInds <- intersect(which(inf$monkid == monk), grep("DNA", inf$samp.loc))

    table(unlist(lapply(strsplit(fastaLabels, " "), length)))
    printRNA <- paste(paste(">", fastaLabels[grep('RNA|PLASMA', fastaLabels)] ,
                             sep = ""), seqsToPrint[rnaInds], sep = "\n")
    write(printRNA, file = paste("../tmp/beast/fasta/", monk, "_RNA.fasta", sep = ""))

    printDNA <- paste(paste(">", fastaLabels[grep('DNA', fastaLabels)] ,
                             sep = ""), seqsToPrint[dnaInds], sep = "\n")
    write(printDNA, file = paste("../tmp/beast/fasta/", monk, "_DNA.fasta", sep = ""))
    
}







