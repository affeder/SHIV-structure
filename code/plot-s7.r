
#Helper plotting function
setUpPlot <- function(ids, ytext, topOrBottom, xvals){

    #Plot a base plot on which elements will be added
    plot(0, type = "n", xlim = c(12, xvals[5]), ylim = c(0 - offset, 1 + offset), ylab = "", xlab = "", main = "", axes = FALSE)

    #Add background graphical elements
    if(ids != "T98133"){polygon(c(26, 26, xvals[5], xvals[5]), c(-.9, 1.5, 1.5, -.9), col = "grey90", border = NA) }
    polygon(c(12, 12, 20, 20), c(-.9, 1.5, 1.5, -.9), col = "grey90", border = NA)
    abline(h = seq(0, 1, by = .25), col = "grey")
    box()

    #If this is the first plot, the y-axis and legend should also be plotted
    if(ids == "T98133"){
        axis(2, at = seq(0, 1, by = .25), labels = paste(seq(0, 1, by = .25)*100, "%", sep = ""), las = 2)
        mtext(ytext, side = 2, line = 3.5)
        if(topOrBottom == "top"){
            legend("topleft", c("PBMC/Plasma", "PBMC", "LN", "Gut", "Vagina", "Mean"),
                   col = c(cols,"black"),
                   pch = c(16),
                   lty = c("solid", "solid", "solid", "solid", "solid"),
                   pt.cex = .75, cex = .85, bg = "white")
        }
    }

    if(topOrBottom == "top"){
        mtext(paste(ids), side = 3)
    }
    if(topOrBottom == "bottom"){
        abline(h = 0)
        plotBottom(ids, c(-1, 0), -.05)
        axis(1, at = xvals)
    }
}



#First, create a list indexed by macaque name 
# of all positions with a >5 count minor allele frequency
# Exclude "-"s
minorallele <- 5
polyindslist <- list()
for(id in monknames){
    specinds <- which(monkid == id)
    poly <- c()
    #Determine which nucleotide position have a minor allele of frequency at least 5
    for(i in 1:ncol(nuc)){
        tab <- table(nuc[specinds,i])
        tab.letters <- tab[names(tab) != "-"]
        poly[i] <- sum(tab.letters >= minorallele) > 1
    }
    polyinds <- which(poly == TRUE)
    polyindslist[[id]] <- polyinds
}


pdf("../out/graphs/S7.pdf", height = 5.25, width = 7.5)

#Set up the layout of the plot
layout(matrix(1:8, nrow = 2))
par(mar = c(0, 1, 0, 0))
par(oma = c(4, 4, 2, 2))

#Set graphical parameters
lwdv <- 1.5
cexv <- 1.5
offset <- .08

#This plot will be by monkey
RNA.inds.all.monks <- c(grep("PLASMA", samp.loc), grep("RNA", samp.loc))
DNA.inds.all.monks <- grep("DNA", samp.loc)

#Set x-axis values based on macaque
approxx <- rbind(c(13, 15, 20, 26, 26),
                 c(13, 15, 20, 26, 29),
                 c(13, 15, 20, 26, 38),
                 c(13, 15, 20, 26, 44))

compsToCheck <- c("PBMC_PLASMA", "PBMC","LN", "GUT", "VAG" )
cols <- c("#228833", "#EE6677", "#AA3377", "#4477AA", "#66CCEE")


for(ids in monknames){

    #Which weeks are relevant for the x-axis
    xvals <- approxx[which(ids == monknames), ]

    potentialinds <- which(monkid == ids)
    RNA.monk.inds <- intersect(potentialinds, RNA.inds.all.monks)
    DNA.monk.inds <- intersect(potentialinds, DNA.inds.all.monks)
    polyinds <- polyindslist[[ids]]

    #Set up matrices to record the percentage of RNA without DNA and vice versa
    toPlotRNA <- matrix(data = NA, nrow = 5, ncol = 5)
    toPlotDNA <- matrix(data = NA, nrow = 5, ncol = 5)

    for(compartment in compsToCheck){

        #Determine the relevant inds for RNA and DNA

        RNAinds <- c()
        DNAinds <- c()
        if(compartment == "LN"){
            RNAinds <- which(samp.loc == "LNRNA")
            DNAinds <- which(samp.loc == "LNDNA")
        }
        if(compartment == "PBMC"){
            RNAinds <- which(samp.loc == "PBMCRNA")
            DNAinds <- which(samp.loc == "PBMCDNA")
        }
        if(compartment == "VAG"){
            RNAinds <- which(samp.loc == "VAGRNA")
            DNAinds <- which(samp.loc == "VAGDNA")
        }
        if(compartment == "GUT"){
            RNAinds <- which(samp.loc == "GUTRNA")
            DNAinds <- which(samp.loc == "GUTDNA")
        }
        if(compartment == "PBMC_PLASMA"){
            RNAinds <- which(samp.loc == "PLASMA")
            DNAinds <- which(samp.loc == "PBMCDNA")
        }
        RNAinds <- intersect(RNAinds, potentialinds)
        DNAinds <- intersect(DNAinds, potentialinds)

        #Set up the plot for that compartment
        colInd <- which(compartment == compsToCheck)

        #For each of these time periods, we're going to break down
        # which vRNA and vDNA is present in a compartment 
        for(timep in 1:5){

            correctweeks <- which(weekind == timep)
            
            #RNA and DNA from the right macaque, location and week
            RNA.comp.week <- intersect(RNAinds, correctweeks)
            DNA.comp.week <- intersect(DNAinds, correctweeks)

            #Turn these sequences into strings so they can be compared more easily
            DNA.seq <- apply(nuc[DNA.comp.week, polyinds], 1, paste, collapse = "")
            RNA.seq <- apply(nuc[RNA.comp.week, polyinds], 1, paste, collapse = "")

            #If there's enough RNA and DNA, we compare them:
            if(length(DNA.seq) >= 3  & length(RNA.seq) >= 3){
                DNAs.with.no.RNAs <- setdiff(DNA.seq, RNA.seq)
                RNAs.with.no.DNAs <- setdiff(RNA.seq, DNA.seq)
                DNAs.with.RNAs <- intersect(DNA.seq, RNA.seq)

                #How many DNAs match no RNAs 
                D.no.R.num <- sum(grepl( paste(DNAs.with.no.RNAs, collapse = "|"), DNA.seq))
                #How many RNAs match with DNAs
                D.with.R.num <- sum(grepl( paste(DNAs.with.RNAs, collapse = "|"), DNA.seq))

                #How many RNAs match no DNAs 
                R.no.D.num <- sum(grepl( paste(RNAs.with.no.DNAs, collapse = "|"), RNA.seq))
                #How many RNAs match with DNAs
                R.with.D.num <- sum(grepl( paste(DNAs.with.RNAs, collapse = "|"), RNA.seq))

                #Compute the proportions of RNA without DNA and DNA without RNA and store it
                R.no.D.prop <- R.no.D.num/(R.with.D.num + R.no.D.num)
                D.no.R.prop <- D.no.R.num/(D.with.R.num + D.no.R.num)

                toPlotRNA[colInd, timep] <- R.no.D.prop
                toPlotDNA[colInd, timep] <- D.no.R.prop
            }
        }
    }

    setUpPlot(ids, "% vRNA without vDNA", "top", xvals)

    #For each compartment
    for(i in 1:5){
        #Determine which timepoints have enough data to plot
        nonNAs <- which(!is.na(toPlotRNA[i,]))
        #And then plot them, if possible
        if(length(nonNAs) > 0){
            lines(xvals[nonNAs], toPlotRNA[i,][nonNAs],
                  col = cols[i], lty = "solid", lwd = lwdv)
            points(xvals[nonNAs], toPlotRNA[i,][nonNAs],
                   pch = 16, col = cols[i], cex = cexv)
        }
     }

    #Also plot the means
    rnameans <- apply(toPlotRNA, 2, mean, na.rm = TRUE)
    nonNAs <- which(!is.nan(rnameans))
    if(length(nonNAs) > 0){
        lines(xvals[nonNAs], rnameans[nonNAs], col = "black", lty = "solid", lwd = lwdv)
        points(xvals[nonNAs], rnameans[nonNAs], pch = 16, col = "black", cex = cexv)
    }

    setUpPlot(ids, "% vDNA without vRNA", "bottom", xvals )
    
    #For each compartment
    for(i in 1:5){
        #Determine which time points have enough data to plot
        nonNAs <- which(!is.na(toPlotDNA[i,]))
        if(length(nonNAs) > 0){
            lines(xvals[nonNAs], toPlotDNA[i,][nonNAs],
                  col = cols[i], lty = "dashed", lwd =lwdv)
            points(xvals[nonNAs], toPlotDNA[i,][nonNAs],
                   pch = 17, col = cols[i], cex = cexv)
        }
    }
    #And plot the means
    dnameans <- apply(toPlotDNA, 2, mean, na.rm = TRUE)
    nonNAs <- which(!is.nan(dnameans))
    if(length(nonNAs) > 0){
        lines(xvals[nonNAs], dnameans[nonNAs], col = "black", lty = "dashed", lwd = lwdv)
        points(xvals[nonNAs], dnameans[nonNAs], pch = 17, col = "black", cex = cexv)
    }
}
mtext("Week Post-Infection", side = 1, outer = TRUE, line = 2.5)
dev.off()


