#kstList takes a list of vectors and computes Kst assuming that
# each list entry has a vector of indices belonging to a different compartments
kstList <- function(listOfInds){
    #sample size across the entries
    fullSize <- length(unlist(listOfInds))
    #Compute the average pairwise diversity within a compartment
    intra <- 0
    for(k in 1:length(listOfInds)){
        tmpinds <- listOfInds[[k]]
        #distMat stores the number of differences between different indices
        intra <- intra + mean(distMat[tmpinds, tmpinds]) * length(tmpinds)
    }
    intra <- intra/fullSize
    #Computer the average pairwise diversity overall (ignoring compartments)
    inter <- mean((distMat[unlist(listOfInds), unlist(listOfInds)]))
    return( (inter - intra)/inter )
}

#randIt takes a list of index vectors and scrambles the compartment assignments
randIt <- function(listOfInds){
    #combine all of the entries to shuffle
    unlisted <- unlist(listOfInds)
    listToRet <- list()
    #Shuffle the order of the list
    toSplit <- sample(unlisted, length(unlisted), replace = FALSE)
    splits <- c(0, cumsum(unlist(lapply(listOfInds, length))))
    for(i in 1:length(listOfInds)){
        #Assign entries based on the new order
        listToRet[[i]] <- toSplit[(splits[i]+1):splits[i+1]]
    }
    return(listToRet)
}

