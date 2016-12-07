source("kst-funs.r")

#We're going to be generating random numbers, so we'll set a random seed
set.seed(7183318)

rnacomps <- c("VAGRNA", "PBMCRNA", "PLASMA", "LNRNA", "GUTRNA")

# First, for each macaque, and RNA from each pair of compartments,
# we'll compute KST - then, we can compare those values to the randomized ones
# kstlist will be a nested list organized as follows
# kstlist[[ name of macaque ]][[ timepoint in question ]][[ RNA comps in question ]]
kstlist <- list()
# stats.nulls will keep track of the same information for randomized trials
kstlist.rand <- list()
#nullNum species the num. of randomized computations from each compartment comparison
nullNum <- 10000
#K = 2 specifies we will be doing pairwise comparisons
K <- 2
for(monk in monknames){
    print(monk)
    #allInfForMonk will be a temporary list to store all the timepoint info
    allInfForMonk <- list()
    allInfForMonk.rand <- list()
    for(t in 1:5){
        #Determine all inds that are relevant for that macaque and that time point
        weekinds <- which(weekind == t)
        monkinds <- which(monkid == monk)
        #store them in weekmonk
        weekmonk <- intersect(weekinds, monkinds)
        #allCompsAtT is a temporary list to start all pairwise comparison info
        allCompsAtT <- list()
        allCompsAtT.rand <- list()
        #Which comparisons need to be checked?
        AllCombs <- combn(rnacomps,K)
        for(i in 1:ncol(AllCombs)){
            #Create a list with indices from one compartment in one entry
            #and indices from the other compartment in the other entry
            listOfInds <- list()
            for(l in 1:K){
                listOfInds[[l]] <- intersect(weekmonk, which(samp.loc == AllCombs[l,i]))
            }
            #If both compartments have at least 3 sequences
            if( sum(unlist(lapply(listOfInds, function(x){length(x) >= 3}))) == K){
                comparisonName <- paste(AllCombs[,i], collapse = "-")
                #compute kst
                toSto <- kstList(listOfInds)
                allCompsAtT[[comparisonName]] <- toSto
                #Compute randomized Kst 
                #First, make a vector to store the values
                toStoRand <- c()
                for(j in 1:nullNum ){
#                    if(j %% 50 == 0){ print(j) }
                    #shuffle the assignments in the list keeping track of
                    # which indices belong to which compartment
                    listOfInds <- randIt(listOfInds)
                    #Recompute the statistic
                    toStoRand[j] <- kstList(listOfInds)
                }
                allCompsAtT.rand[[comparisonName]] <- toStoRand
            }
        }
        allInfForMonk[[paste(t)]] <- allCompsAtT
        allInfForMonk.rand[[paste(t)]] <- allCompsAtT.rand
    }
    kstlist[[monk]] <- allInfForMonk        
    kstlist.rand[[monk]] <- allInfForMonk.rand
}


#Benjamini-Hochberg procedure to determine the pvalue cutoff
#First, we'll make a list of all the pvalues, and then use that list to compute 
# our p-value cutoff
comparisons <- c("PLASMA-GUTRNA", "PLASMA-LNRNA", 'VAGRNA-PLASMA', 'PBMCRNA-PLASMA',
                 "PLASMA-GUTRNA", "LNRNA-GUTRNA", "PBMCRNA-GUTRNA", "VAGRNA-GUTRNA")
allvals <- c()
for(monk in monknames){
    for(t in paste(1:5)){
        for(comp in comparisons){
            statval <-  kstlist[[monk]][[t]][[comp]]
            nullvals <- kstlist.rand[[monk]][[t]][[comp]]
            if(length(statval) > 0){
                allvals <- c(allvals, sum(statval < nullvals)/length(nullvals))
            }
        }
    }
}
lenval <- length(allvals)
cm <- sum(1/(1:lenval))
correct <- .05/(lenval*cm)*(1:lenval)
rejInd <- max(which(correct > sort(allvals)))
pvalcut <- correct[rejInd]
print(pvalcut)


