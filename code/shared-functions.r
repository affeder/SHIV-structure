#function to plot the treatment label on the bottom of the plot
plotBottom <- function(monk, boxyvals, laby, trunc = FALSE){
    print(paste("plotting ",monk))
    polygon(x = c(12, 20, 20, 12),
            y = c(boxyvals[2], boxyvals[2], boxyvals[1], boxyvals[1]),
            col = "black")
    if(monk == "A99165"){
        text(16, laby, "EFV", col = "white", cex = 1)    
        polygon(x = c(26, 29, 29, 26),
                y = c(boxyvals[2], boxyvals[2], boxyvals[1], boxyvals[1]),
                col = "black")
        if(trunc == FALSE){
            text(27.5, laby, "Rx1", col = "white", cex = 1) }
    }
    if(monk == "T98133"){
        text(16, laby, "FTC", col = "white", cex = 1)    
    }
    if(monk == "A99039"){
        polygon(x = c(26, 38, 38, 26),
                y = c(boxyvals[2], boxyvals[2], boxyvals[1], boxyvals[1]),
                col = "black")
        if(trunc == FALSE){
            text(32, laby, "Rx2", col = "white", cex = 1)
        }
        text(16, laby, "FTC", col = "white", cex = 1)    
    }
    if(monk == "A01198"){
        if(trunc == FALSE){
            polygon(x = c(26, 44, 44, 26),
                    y = c(boxyvals[2], boxyvals[2], boxyvals[1], boxyvals[1]),
                    col = "black")
            text(35, laby, "Rx1 ", col = "white", cex = 1)
        }
        text(16, laby, "EFV", col = "white", cex = 1)    
    }
}

#Function to determine what percentage of the RNA is drug resistant at a certain time
percentageDR <- function(monk, weekvalind, mutIdent, mutPos){
    #First get the relevant indices. They should be all inds at which
    relInds <- intersect(
        #The monkey is correct
        which(monkid == monk),
        intersect(
            #The sample location is either plasma or vRNA
            grep("RNA|PLASMA", samp.loc),
            #The week ind is correct
            which(weekind == weekvalind )
        )
    )
    #Return the proportion of the compartment that has the correct aa (relative to all samples at that week)
    return(length(grep(mutIdent, aa[relInds,mutPos]))/length(relInds))
}


#Function to determine what percentage of the RNA is drug resistant at a certain time
percentageDR.wkind <- function(monk, weekvalind, mutIdent, mutPos){
    relInds <- intersect( which(monkid == monk), intersect(
                                    grep("RNA|PLASMA", samp.loc),
                                    which(weekind == weekvalind )
    )
                         )
    return(length(grep(mutIdent, aa[relInds,mutPos]))/length(relInds))
}

tcol <- function(color, percent) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  return(t.col)
}

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
