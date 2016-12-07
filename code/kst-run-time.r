source("kst-funs.r")

numreps <- 10000
set.seed(12093580)

#Create a list to keep track of the true KST value between time points
listToPlotRNA <- list()
#Create a list to keep track of the KST values from randomized assignments to time points 
listToPlotRNA.rand <- list()
#For both of the above lists, there will be an organization as follows:
#listToPlotRNA[[macaque name]][[compartment]][[timepoint comparison]]
rnacomps <- c( "GUTRNA", "VAGRNA", "PBMCRNA", "LNRNA", "PLASMA")
for(monk in monknames){
    monkinds <- which(monkid == monk)
    #create temporary lists to accept values (compartment level for each macaque)
    monksToPlot <- list()
    monksToPlot.rand <- list()
    for(comp in rnacomps){
        compinds <- which(samp.loc == comp)
        #create further temporary lists (timepoint level for each compartment)
        compsToPlot <- list()
        compsToPlot.rand <- list()
        for(t in 1:4){
            #We are going to try to match compartments between adjacent timepoints
            t1s <- which(weekind == t)
            inds1 <- intersect(t1s, intersect(compinds, monkinds))
            #If there are at least three sequences at the first time point, look for a second
            # time point against which it can be compared
            if(length(inds1) >= 3){
                #Start by trying to compare it to the next timepoint
                toAdd <- 1
                t2s <- which(weekind == t+toAdd)
                inds2 <- intersect(t2s, intersect(compinds, monkinds))
                #If the next time point doesn't have at least 3 sequences, 
                #look at the time point after that one. 
                while(length(inds2) < 3 & (t + toAdd + 1) <= 5){
                    toAdd <- toAdd + 1
                    t2s <- which(weekind == t+toAdd)
                    inds2 <- intersect(t2s, intersect(compinds, monkinds))
                }
                #If two time points with sufficient sequences were found
                if(length(inds1) >= 3 & length(inds2) >= 3){
                    #record the comparison name
                    compName <- paste(t, t + toAdd, sep = "-")
                    #and then compute and store Kst between them
                    datVals <- kstList(list(inds1, inds2))
                    compsToPlot[[compName]] <- datVals

                    #and between randomized versions of those comparisons
                    randVals <- c()
                    for(i in 1:numreps){
                        randVals[i] <- kstList(randIt(list(inds1, inds2)))
                    }
                    compsToPlot.rand[[compName]] <- randVals
                }
            }
        }
        #Percolate down the storage of the nested lists
        monksToPlot[[comp]] <- compsToPlot
        monksToPlot.rand[[comp]] <- compsToPlot.rand
    }
    listToPlotRNA[[monk]] <- monksToPlot
    listToPlotRNA.rand[[monk]] <- monksToPlot.rand
}

#Now, repeat this with vDNA
#Create a list to keep track of the true KST value between time points
listToPlotDNA <- list()
#Create a list to keep track of the KST values from randomized assignments to time points 
listToPlotDNA.rand <- list()
#For both of the above lists, there will be an organization as follows:
#listToPlotDNA[[macaque name]][[compartment]][[timepoint comparison]]
dnacomps <- c("GUTDNA", "VAGDNA", "PBMCDNA", "LNDNA")
for(monk in monknames){
    monkinds <- which(monkid == monk)
    #create temporary lists to accept values (compartment level for each macaque)
    monksToPlot <- list()
    monksToPlot.rand <- list()
    for(comp in dnacomps){
        compinds <- which(samp.loc == comp)
        #create further temporary lists (timepoint level for each compartment)
        compsToPlot <- list()
        compsToPlot.rand <- list()
        for(t in 1:4){
            #We are going to try to match compartments between adjacent timepoints
            t1s <- which(weekind == t)
            inds1 <- intersect(t1s, intersect(compinds, monkinds))
            #If there are at least three sequences at the first time point, look for a second
            # time point against which it can be compared
            if(length(inds1) >= 3){
                #Start by trying to compare it to the next timepoint
                toAdd <- 1
                t2s <- which(weekind == t+toAdd)
                inds2 <- intersect(t2s, intersect(compinds, monkinds))
                #If the next time point doesn't have at least 3 sequences, 
                #look at the time point after that one. 
                while(length(inds2) < 3 & (t + toAdd + 1) <= 5){
                    toAdd <- toAdd + 1
                    t2s <- which(weekind == t+toAdd)
                    inds2 <- intersect(t2s, intersect(compinds, monkinds))
                }
                #If two time points with sufficient sequences were found
                if(length(inds1) >= 3 & length(inds2) >= 3){
                    #record the comparison name
                    compName <- paste(t, t + toAdd, sep = "-")
                    #and then compute and store Kst between them
                    datVals <- kstList(list(inds1, inds2))
                    compsToPlot[[compName]] <- datVals

                    #and between randomized versions of those comparisons
                    randVals <- c()
                    for(i in 1:numreps){
                        randVals[i] <- kstList(randIt(list(inds1, inds2)))
                    }
                    compsToPlot.rand[[compName]] <- randVals
                }
            }
        }
        #Percolate down the storage of the nested lists
        monksToPlot[[comp]] <- compsToPlot
        monksToPlot.rand[[comp]] <- compsToPlot.rand
    }
    listToPlotDNA[[monk]] <- monksToPlot
    listToPlotDNA.rand[[monk]] <- monksToPlot.rand
}


#At this point, compute all p-values, so we can assess what our multiple hypothesis correction
# should be
allvals <- c()
for(monk in monknames){
    for(comp in names(listToPlotDNA.rand[[monk]])){    
        for(t in names(listToPlotDNA.rand[[monk]][[comp]])){
            statval <-  listToPlotDNA[[monk]][[comp]][[paste(t)]]
            nullvals <- listToPlotDNA.rand[[monk]][[comp]][[paste(t)]]
            if(length(statval) > 0){
                allvals <- c(allvals, sum(statval < nullvals)/length(nullvals))
            }
        }
    }
    for(comp in names(listToPlotRNA.rand[[monk]])){    
        for(t in names(listToPlotRNA.rand[[monk]][[comp]])){
            statval <-  listToPlotRNA[[monk]][[comp]][[paste(t)]]
            nullvals <- listToPlotRNA.rand[[monk]][[comp]][[paste(t)]]
            if(length(statval) > 0){
                allvals <- c(allvals, sum(statval < nullvals)/length(nullvals))
            }
        }
    }
}
allvals <- allvals[!is.na(allvals)]
#Benjamini-Hochberg procedure
cm <- sum(1/(1:length(allvals)))
correct <- (.05/(length(allvals) * cm) )*(1:length(allvals))
pvalcut.time <- correct[rejInd]
print(pvalcut.time)

