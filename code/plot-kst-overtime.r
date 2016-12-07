
#meanMe computes the average Kst for each macaque and timepoint
# It accepts a nested list structure as follows: 
# Indlist[[macaque name]][[time period (1-5)]][[RNA compartment]]
# If a given macaque at a given time point has compartments A, B, C and D
# meanMe will compute the average of Kst between 
#   A and (B, C and D)
#   B and (A, C and D)
#   C and (A, B and D)
#   D and (A, B and C)
# and returns the average of the four comparisons
meanMe <- function(IndList){
    return(lapply(IndList, meanAcrossTs))
}

#This is a helper function to compute the mean Kst across timepoints (within a macaque)
meanAcrossTs <- function(IndListForMonk){
    return(unlist(lapply(IndListForMonk, meanpairwiseForThisT)))
}

#This is a helper function to compute the mean Kst across compartments (for a given time)
meanpairwiseForThisT <- function(timeListAcrossComps){
    return(mean(unlist(lapply(timeListAcrossComps,  pairwiseForThisComp, unlist(timeListAcrossComps))), na.rm = TRUE))
}

#This is a helper function to compute the mean Kst within a compartment 
# (as compared to other compartments at that time point)
pairwiseForThisComp <- function(focalInds, allIndsForThisT){
    #How many sequences are necessary for a compartment to be included?
    lowerlim <- 3
    nonfocalInds <- setdiff(allIndsForThisT, focalInds)
    if(length(focalInds) >= lowerlim & length(nonfocalInds) >= lowerlim){
        return(kstList(list(focalInds, nonfocalInds)))
    }else{
        return(NA)
    }
}


#This is a function that takes a nested index list (as described above)
# and randomizes the list such that the compartment and time structure is maintained
# and then returns a list of the same structure (on which meanMe can be called)
# For example, a monkey may have a list like this:
# IndList[['fakemonkey']][['1']]
# $GUTRNA
# 1 2 3 4 5
# $PLASMA
# 6 7 8 9 10
# $VAGRNA
# 11 12 13
# IndList[['fakemonkey']][['2']]
# $GUTRNA
# 14 15 16 17 18 19 20
# $PLASMA
# 21 22 23
#
# And this function will return a list that looks like this:
# IndList[['fakemonkey']][['1']]
# $GUTRNA
# 4 9 12 3 6
# $PLASMA
# 2 11 13 1 8
# $VAGRNA
# 10 5 7
# IndList[['fakemonkey']][['2']]
# $GUTRNA
# 19 16 22 17 21 14 23
# $PLASMA
# 20 18 15
# 
# Each macaque will have the same time points, and the same 
# compartments of the same size, but indices will be shuffled within timepoints
IndListRand <- function(IndList){
    
    RandListForMonk <- function(IndListForMonk){
        return(lapply(IndListForMonk, RandListForT))
    }
    
    RandListForT <- function(IndListForT){
        compSizes <- unlist(lapply(IndListForT, length))
        allInds <- unlist(IndListForT)
        #Randomize order
        sampToAdd <- sample(allInds, size = length(allInds))
        ListToRet <- list()
        cumS <- c(1, cumsum(compSizes))
        for(i in 1:length(IndListForT) ){
            if(cumS[i] != cumS[i+1]){
                ListToRet[[i]] <- sampToAdd[cumS[(i)]:cumS[(i+1)]]
            }else{
                ListToRet[[i]] <- integer(0)
            }
        }
        names(ListToRet) <- names(IndListForT)
        return(ListToRet)
    }
    return(lapply(IndList, RandListForMonk))
}


rnaTypes <- c("PLASMA", "PBMCRNA", "GUTRNA", "LNRNA", "VAGRNA")

#First, we're just going to compute a list of indices by macaque,
# timepoint and RNA compartment
IndList <- list()
for(monk in monknames){
    mInds <- which(monkid == monk)
    timeList <- list()
    for(t in 1:5){
        tInds <- which(weekind == t)
        typeList <- list()
        for(type in rnaTypes){
            cInds <- which(samp.loc == type)
            typeList[[type]] <- intersect(intersect(mInds, cInds), tInds)
        }
        timeList[[paste(t)]] <- typeList
    }
    IndList[[monk]] <- timeList
}
#IndList should be organized as follows
# Indlist[[macaque name]][[time period (1-5)]][[RNA compartment]]


set.seed(102093510)
#How many randomizations should happen?
numRands <- 10000
#Set up a list and pre-population it was arrays of the desired size
randLists <- list()
for(i in monknames){
    randLists[[i]] <- matrix(data = NA, nrow = numRands, ncol = 5)
}
#For numRands number of times
for(i in 1:numRands){
    #Compute the means (using meanMe()) on the randomized indexed lists (using IndlistRand())
    # of the true index lists
    stoResults <- meanMe(IndListRand(IndList))
    #store the results for each macaque
    for(monk in monknames){
        randLists[[monk]][i,] <- stoResults[[monk]]
    }
}
#Compute the central most randomizations as a measure of confidence
ranges99 <- lapply(randLists, function(x){apply(x,2, quantile, c(0.005, .995), na.rm= TRUE)})
ranges95 <- lapply(randLists, function(x){apply(x,2, quantile, c(0.025, .975), na.rm= TRUE)})

#Plot helper function
setUpPlot <- function(monk, ylims, xlims, xmax){
    plot(0, type = "n", xlim = c(12, xmax), ylim = ylims, axes = FALSE, ylab = "", xlab = "")
    polygon(x = c(12, 20, 20, 12), y= c(-1, -1, 1, 1), col ="lightgrey", border = NA)
    if(monk == "A01198"){     polygon(x = c(12, 20, 20, 12), y= c(-1, -1, 1, 1), col ="grey95", border = NA) } 
    if(monk!= "T98133"){ polygon(x = c(26, xlims[5], xlims[5], 26), y= c(-1, -1, 1, 1), col ="lightgrey", border = NA)}
    box()
    mtext(paste(monk), 3)
    abline(v = xlims[5], lty = "dotted")
    if(monk == "A99165" | monk == "A01198"){ axis(1, xlims ) }
    if(monk == "T98133" | monk == "A01198"){ axis(2, seq(0, max(ylims), .01)) }
    if(monk == "A99165" | monk == "A99039"){ axis(2, seq(0, max(ylims), .05)) }
    if(monk == "A99165"){
        legend("topright", c("Observed data", "Randomized\ndata"), col = c(redcol,tcol(redcol, 80)), pch = c(16, NA), lwd = c(2, 10), title = expression("Average K"["ST"]),  cex = .85, bty = "n")
    }
}

pdf("../out/graphs/migvsel.pdf", width = 6, height = 4.5)
layout(matrix(1:4, nrow = 2))
par(mar = c(1.5, 2.5, 0, 0))
par(oma = c(2, 2, 1, 1))
approxx <- rbind(c(13, 15, 20, 26, 26),
                 c(13, 15, 20, 26, 29),
                 c(13, 15, 20, 26, 38),
                 c(13, 15, 20, 26, 44))
maxx <- c(29, 29, 44, 44)
redcol <- brewer.pal(3, "Set1")[1]
for(monk in monknames){
    xax <- approxx[which(monknames == monk),]
    #use meanMe() to determine the true values of mean Kst at a timepoint
    trueVals <- meanMe(IndList)[[monk]]
    maxy <- max(c(ranges95[[monk]], trueVals), na.rm = TRUE)
    setUpPlot(monk = monk, ylims = c(-maxy/10,  maxy),
              xlims = approxx[which(monk == monknames), ],
              xmax =  maxx[which(monk == monknames)])
    inc <- which(!is.na(ranges99[[monk]][1,]))
    #determine the values for the ranges 
    ys95 <- c(ranges95[[monk]][1,inc], rev(ranges95[[monk]][2,inc]))
    xs95 <- c(xax[inc], rev(xax[inc]))
    #and plot
    polygon(xs95, ys95, col = tcol(redcol, 80), border= NA)
    #then, plot the true values as points and a line
    points(xax[inc], trueVals[inc], col = redcol, pch = 16)
    lines(xax[inc],trueVals[inc], col = redcol, lwd = 1.5)
    abline(h = 0)
    plotBottom(monk, c(-1, 0), -maxy/15)
}
mtext("Average Kst between compartments", side = 2, outer = TRUE)
mtext("Week Post-Infection", side = 1, outer = TRUE, line = 1)
dev.off()

