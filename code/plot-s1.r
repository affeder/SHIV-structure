#Pairwise distances between each pair of sequences (as a proportion of total variable sites)
distMatFrac <- distMat/length(135:900)

#The pi function will use distMatFrac (the distance matrix) to quickly compute pi
pi <- function(inds){
    if(length(inds) == 0){ return(NULL) }
    return(mean(apply((distMatFrac[inds, inds]), 1, mean)))
}


#All compartments over which diversity should be examined
allcomps <- c("VAGRNA", "PBMCRNA", "PLASMA", "LNRNA", "GUTRNA",
              "VAGDNA", "PBMCDNA", "LNDNA", "GUTDNA")
stoList <- list()
for(monk in monknames){
    for(weekval in 1:5){
        for(comp in allcomps){
            #For a particular week, macaque and compartment
            relInds <- intersect(which(monkid == monk & weekind == weekval),
                                 which(samp.loc == comp))
            #If that set of sequences has at least three sequences
            if(length(relInds) >= 3){
                #Compute pi
                stoList[[paste(monk, weekval, comp, sep = "-")]] <- pi(relInds)
            }
        }
        #Repeat the same process separating out samples by RNA only
        relInds <- intersect(grep("PLASMA|RNA", samp.loc),
                             which(monkid == monk & weekind == weekval))
        if(length(relInds) >= 3){
            stoList[[paste(monk, weekval, "allrna", sep = "-")]] <- pi(relInds)
        }
        #Repeat the same process separating out samples by DNA only
        relInds <- intersect(grep("DNA", samp.loc),
                             which(monkid == monk & weekind == weekval))
        if(length(relInds) >= 3){
            stoList[[paste(monk, weekval, "alldna", sep = "-")]] <- pi(relInds)
        }
    }
}
#Now, stoList should have a list of pi over time among 
# all sequences, and among just RNA and just DNA

#plotfeatures just puts a lot of the plotting functionality together
plotfeatures <- function(monk, xvals, dnaorrna){
    polygon(c(12, 20, 20, 12), c(-1, -1, 1, 1), col = "lightgrey", border = FALSE)
    if(monk != "T98133"){ polygon(c(26, xvals[5], xvals[5], 26),
           c(-1, -1, 1, 1), col = "lightgrey", border = FALSE)}
    text(paste(dnaorrna), x = xvals[5]*1.025, y =  ylimmax*.95, pos = 2)
    if(dnaorrna == "vDNA"){
        axis(1, xvals)
    }
    if(monk == "T98133"){
        axis(2, las = 2)
    }
    if(dnaorrna == "vRNA"){
        mtext(monk, 3)
    }
}



pdf("../out/graphs/S1.pdf", height = 4.75, width = 7.5)
layout(matrix(1:8, nrow = 2))
par(mar = c(1, 1, 0, 0))
par(oma = c(4, 4, 2, 2))
#Order: "Plasma", "PBMC",  "LN", "Gut", "Vagina"
pal <- c("#228833", "#EE6677", "#AA3377", "#4477AA", "#66CCEE")
#rnacomps <- c("VAGRNA", "PBMCRNA", "PLASMA", "LNRNA", "GUTRNA", "allrna")
rnacomps <- c("PLASMA", "PBMCRNA", "LNRNA", "GUTRNA", "VAGRNA", "allrna")
#dnacomps <- c("VAGDNA", "PBMCDNA", "LNDNA", "GUTDNA", "alldna")
dnacomps <- c("PBMCDNA", "LNDNA", "GUTDNA", "VAGDNA", "alldna")
rnacols <- c(pal, "black")
dnacols <- c(pal[2:5], "black")
ylimmax <- .0135
ylimmin <- -.0005
approxx <- rbind(c(13, 15, 20, 26, 26),
                 c(13, 15, 20, 26, 29),
                 c(13, 15, 20, 26, 38),
                 c(13, 15, 20, 26, 44))

#We are going to use these two lists to keep track of data points for some regression later on
RNAcrossmonks <- list()
DNAcrossmonks <- list()

for(monk in monknames[1:4]){

    toPlotX <- approxx[which(monknames == monk), ]

    #Set up the RNA plot for the macaque
    plot(0, type = "n", xlim = c(12, toPlotX[5]), ylim = c(ylimmin,ylimmax), axes = FALSE)
    plotfeatures(monk, toPlotX, "vRNA")
    if(monk == "T98133"){
        #Plot a legend in the T98133 plot
        legend("topleft", c("Plasma",  "PBMC", "LN", "Gut", "Vagina"), col = c(rnacols),
               pch = c(16), lty = c("solid", "solid", "solid", "solid"),
               pt.cex = .75, cex = .85, bg = "white") }
    #We are going to use allRNAs to do some linear regression comparison later in this script
    allRNAs <- matrix(data= NA, nrow = 5, ncol = length(toPlotX))
    for(comp in rnacomps){
        #For each compartment, toPlotY will hold the pi values at the different times
        toPlotY <- rep(NA, 5)
        for(weekval in 1:5){
            if(!is.null(stoList[[paste(monk, weekval, comp, sep = "-")]])){
                toPlotY[weekval] <- stoList[[paste(monk, weekval, comp, sep = "-")]]
            }
        }

        rnaindsto <- which(comp == rnacomps)
        #This is to determine that we're actually looking at RNA
        if(rnaindsto <= nrow(allRNAs)){
            #if sto, store the pi values in a vector that will eventually be
            # put in a list and re-examined via regression
            allRNAs[rnaindsto, ] <- toPlotY
        }
        
        #all non-NA points should be plotted, with lines
        inc <- which(!is.na(toPlotY))
        lines(toPlotX[inc], toPlotY[inc], col = rnacols[which(comp == rnacomps)],
              lty = "solid")
        points(toPlotX[inc], toPlotY[inc], col = rnacols[which(comp == rnacomps)], pch = 16)

        #There's one RNA point that we want to plot an asterisk 
        # near to indicate viral suppression
        if(monk == "A99039"){
            text(38, .0032, "*", cex = 2.5)
        }
    }
    box()

    #Set up the DNA plot for the macaque
    plot(0, type = "n", xlim = c(12,  toPlotX[5]), ylim = c(ylimmin,ylimmax), axes = FALSE)
    plotfeatures(monk, toPlotX, "vDNA")

    # Similar to allRNAs, we'll use allDNAs for some later regression
    allDNAs <- matrix(data= NA, nrow = 4, ncol = length(toPlotX)) 
    for(comp in dnacomps){

        #toPlotY holds the pi values that will be plotted for a given compartment
        toPlotY <- rep(NA, 5)
        for(weekval in 1:5){
            if(!is.null(stoList[[paste(monk, weekval, comp, sep = "-")]])){
                toPlotY[weekval] <- stoList[[paste(monk, weekval, comp, sep = "-")]]
            }
        }
        dnaindsto <- which(comp == dnacomps)
        if(dnaindsto <= nrow(allDNAs)){
            allDNAs[dnaindsto, ] <- toPlotY
        }

        #plot all non-NA points/lines
        inc <- which(!is.na(toPlotY))
        lines(toPlotX[inc], toPlotY[inc], col = dnacols[which(comp == dnacomps)], lty = "solid")
        points(toPlotX[inc], toPlotY[inc], col = dnacols[which(comp == dnacomps)], pch = 16)
    }
    box()
    plotBottom(monk, c(0, -1), -.0005)
    abline(h = 0)

    #store alRNAs and allDNAs for later use
    colnames(allRNAs) <- toPlotX
    RNAcrossmonks[[monk]] <- allRNAs
    colnames(allDNAs) <- toPlotX
    DNAcrossmonks[[monk]] <- allDNAs
}
mtext( "Week Post-Infection", side = 1, outer = TRUE, line = 2)
mtext( "Average Pairwise Nucleotide Diversity", side = 2, outer = TRUE, line = 2.5)
dev.off()


#This is somewhat unintuitively written piece of code that flattens RNAcrossmonks
# into a list of matrices of timepoints and pi values
ListOfMatrixOfPisRNA <- lapply(RNAcrossmonks, function(monk.x){
    return( matrix(unlist(apply( monk.x , 1, function(x){
        toRet <- cbind(as.numeric(names(x)), x)[!is.na(x),]
        ordering <- c()
        if(length(toRet) > 2){
            for(i in 1:nrow(toRet)){ ordering <- c(ordering, toRet[i, 1], toRet[i, 2]) }
            return(ordering)
        }else if(length(toRet) == 2){
            return(toRet)
        }else{
            return(c(NA, NA))
        }
    } )), ncol = 2, byrow = TRUE)) } )

ListOfMatrixOfPisDNA <- lapply(DNAcrossmonks, function(monk.x){
    return( matrix(unlist(apply( monk.x , 1, function(x){
        toRet <- cbind(as.numeric(names(x)), x)[!is.na(x),]
        if(length(toRet) > 2){
            ordering <- c()
            for(i in 1:nrow(toRet)){ ordering <- c(ordering, toRet[i, 1], toRet[i, 2]) }
            return(ordering)
        }else if(length(toRet) == 2){
            return(toRet)
        }else{
            return(c(NA, NA))
        }
    } )), ncol = 2, byrow = TRUE)) } )


#Set up the table to print the regression analysis
tabToPrint <- matrix(data = NA, nrow = 4, ncol = 9)
for(monk in monknames){
    #Set up the matrices that will be used to fit the models for both RNA and DNA
    # The third column will have a boolean value for isDNA?
    rnadat <- as.data.frame(cbind(ListOfMatrixOfPisRNA[[monk]],
                                  rep(0, nrow(ListOfMatrixOfPisRNA[[monk]]))))
    dnadat <- as.data.frame(cbind(ListOfMatrixOfPisDNA[[monk]],
                                  rep(1, nrow(ListOfMatrixOfPisDNA[[monk]]))))
    names(rnadat) <- c( "week", "pi", "dna")
    names(dnadat) <- c( "week", "pi", "dna")
    #straightforward regression of diversity on week
    rna.m <- lm(pi ~ week, data =rnadat)
    dna.m <- lm(pi ~ week, data= dnadat)
    rna.model <- rna.m$coefficients
    dna.model <- dna.m$coefficients

    #ANOVA of nested models to determine whether there was a difference between RNA and DNA
    tmpdat <- as.data.frame(rbind(rnadat, dnadat))
    #fm = full model, nm = null model
    fm <- lm(pi ~ week*dna, data =tmpdat)
    nm <- lm(pi ~ week, data= tmpdat)

    #rowvals will be what will be printed - it stores the coefficients and pvalues
    # for both the intercepts and the coefficients for RNA and DNA model fits
    # and also the p-value for the ANOVA comparing the two
    rowvals <- c((summary(rna.m)$coefficients[1,c(1,4)]),
                 (summary(rna.m)$coefficients[2,c(1,4)]),
                 (summary(dna.m)$coefficients[1,c(1,4)]),
                 (summary(dna.m)$coefficients[2,c(1,4)]),
                 (anova(fm, nm)["Pr(>F)"])[2,1])
    tabToPrint[which(monk == monknames), ] <- rowvals
}

#semi-informative names
colnames(tabToPrint) <- c("int-r", "p-int-r", "weeks-r", "p-weeks-r",
                          "int-d", "p-int-d", "weeks-d", "p-weeks-d",
                          "anova")

#Add macaque names
tabToPrint.withrownames <- data.frame( row = monknames, tabToPrint)
#and print
print.xtable(xtable(tabToPrint.withrownames),  type = "html",
             paste("../out/tables/piovertime.html", sep = ""), include.rownames=FALSE)









#I am not sure if I use this part below:


#comparison (fit as part of the same model)
layout(matrix(1:4, nrow = 2))
for(monk in monknames){
    rnadat <- cbind(ListOfMatrixOfPisRNA[[monk]], rep(0, nrow(ListOfMatrixOfPisRNA[[monk]])))
    dnadat <- cbind(ListOfMatrixOfPisDNA[[monk]], rep(1, nrow(ListOfMatrixOfPisDNA[[monk]])))
    tmpdat <- as.data.frame(rbind(rnadat, dnadat))
    names(tmpdat) <- c( "week", "pi", "dna")
    fm <- lm(pi ~ week*dna, data =tmpdat)
    nm <- lm(pi ~ week, data= tmpdat)
    fullmodel <- fm$coefficients
    nullmodel <- nm$coefficients
    plot(0, type = "n", xlim = c(0, 44), ylim = c(0, .01), main = monk)
    abline(fullmodel['(Intercept)'], fullmodel['week'], col = "blue")
    abline(fullmodel['(Intercept)'] + fullmodel['dna'], fullmodel['week'] + fullmodel['week:dna'] , col = "red")
    print(monk)
    print(fm)
    print(anova(fm, nm))
}




