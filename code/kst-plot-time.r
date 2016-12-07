#This script produces a plot of Kst computed between adjacent timepoints

pdf("../out/graphs/betterovertime.pdf", height = 2.5, width = 8)

#This is a helper plotting function to style the arrows 
plotLines <- function(xvals, plotMeVals, colVals, constval){
    cols <- c(NA, brewer.pal(3, "Set1")[1], "darkgrey")
    for(i in 1:(length(xvals) - 1)){
        arrows(rep(xvals[i]+ constval, length(plotMeVals)), plotMeVals, rep(xvals[i+1] - constval, length(plotMeVals)), plotMeVals, length = .025, code = 3, lwd = 3, col = cols[colVals][i])
    }
}


#A few plotting values
nonsigpoint <- 0
lowoff <- .4

#Set graphical parameters
layout(matrix(1:(length(monknames)), nrow = 1))
par(mar = c(0, 1, 0, 0))
par(oma = c(4, 7, 2, 1))

#The different macaques need different x-axis values
approxx <- rbind(c(13, 15, 20, 26, 26), c(13, 15, 20, 26, 29), c(13, 15, 20, 26, 38), c(12, 15, 20, 26, 44))

comparisons <- c("PLASMA", "LNRNA", 'PBMCRNA', 'VAGRNA', "GUTRNA", "LNDNA", "PBMCDNA", "VAGDNA", "GUTDNA")
colvals <- brewer.pal(6, "Set2")

#These will be used to plot the horizontal lines
plotMeValsLow <- nonsigpoint + (0:4)/15

allSumm <- c("NA", "NA", 0, 0)
for(monk in monknames){


    #Determine what the max index of approxx should be used and enter appropriate x axis vals
    upperxlim <- c(4,5, 5, 5)[which(monk == monknames)]
    xvals <- approxx[which(monknames == monk), ]

    for(i in 1:length(comparisons)){

        #For the first comparison plotted, the plot needs to be made
        if(i == 1 ){

            plot(0, type = "n",
                 xlim = c(12, xvals[upperxlim]+.25), ylim = c(-.3875, .3), axes = FALSE)
            mtext(side = 3, paste(monk))

            #If it's the first macaque, it also needs a y-axis
            if(monk == "T98133"){

                #Inner y-axis
                axis(2, at = c(plotMeValsLow[2:5] - lowoff, plotMeValsLow), c(rev(c( "Gut", "Vagina", "PBMC", "Lymph node")), rev(c("Gut", "Vagina", "PBMC", "Lymph node", "Plasma"))), las = 2)
                #Outer y-axis (including the lines below the names)
                mtext(c("vDNA", "vRNA"),2,
                      las = 0, at = c(.2, .75), line = 5.2, outer = TRUE, xpd = NA)
                xvala <- 3.8
                arrows(xvala, -.02, xvala, .3, length = 0, col = "black", xpd = NA)
                arrows(xvala, -.36, xvala, -.1, length = 0, col = "black", xpd = NA)
            }

            #Plot backgound elements of the plots
            if(monk != "T98133"){ polygon(c(26, 26, xvals[5], xvals[5]), c(-1, 2, 2, -1), col = "grey90", border = NA) }
            polygon(c(12, 12, 20, 20), c(-1, 2, 2, -1), col = "grey90", border = NA)
            abline(h = plotMeValsLow, col = "lightgrey")
            abline(h = plotMeValsLow[2:5] - lowoff, col = "lightgrey")
            box()

            #Plot graphical elements related to the x axis
            axis(1, xvals)
            abline(h = -.36)
            plotBottom(monk, c(-1, -.36), c(-.385))
            
        }

        #For each compartment, iterate through all the timepoint comparisons

        plotline <- c()
        t.ind <- 1
        #If the compartment is RNA, it should draw values from the vRNA data
        if(i <= 5){
            tvals <- names(listToPlotRNA[[monk]][[comparisons[i]]])
            for(t in tvals){
                statval <-  listToPlotRNA[[monk]][[comparisons[i]]][[t]]
                nullvals <- listToPlotRNA.rand[[monk]][[comparisons[i]]][[t]]
                pval <- 1-sum(statval >= nullvals)/length(nullvals)
                if(is.na(pval)){
                    plotline[t.ind] <- 1
                }else{
                    if(pval < pvalcut.time){
                        plotline[t.ind] <-  2 
                    }else{
                        plotline[t.ind] <- 3
                    }
                }
                t.ind <- t.ind + 1
            }
        }
        ##If the compartment is vDNA, it should access the vDNA data
        if(i > 5){
            tvals <- names(listToPlotDNA[[monk]][[comparisons[i]]])
            for(t in tvals){
                statval <-  listToPlotDNA[[monk]][[comparisons[i]]][[t]]
                nullvals <- listToPlotDNA.rand[[monk]][[comparisons[i]]][[t]]
                pval <- 1-sum(statval >= nullvals)/length(nullvals)
                if(is.na(pval)){
                    plotline[t.ind] <- 1
                }else{
                    if(pval < pvalcut.time){
                        plotline[t.ind] <-  2 
                    }else{
                        plotline[t.ind] <- 3
                    }
                }
                t.ind <- t.ind + 1
            }
        }

        #Determine at which height to plot these values 
        yvalToPlot <- plotMeValsLow[(i %% 6) + (i > 5)*2 ] - as.numeric(i > 5)*lowoff

        #Determine which weeks are relevant for the particular comparison
        samptimes <- (weeks[intersect(which(samp.loc == comparisons[i]), which(monkid == monk))])
        samptimes <- names(table(samptimes)[table(samptimes) >= 3])

        if(sum(grepl("15", samptimes)) & sum(grepl("16", samptimes)) ){ samptimes[which(samptimes == "15")] <- "16" } #If both 15 and 16 are present, use 16
        if(sum(grepl("20", samptimes)) & sum(grepl("21", samptimes)) ){ samptimes[which(samptimes == "20")] <- 21 } #If both 20 and 21 are present, use 21

#Update xvals so that they reflect the correct sampling weeks for that particular compartment in that particular macaque
        for(k in samptimes){
            if(k == 12){xvals[which(xvals == 13)] <- 12}
            if(k == 16){xvals[which(xvals == 15)] <- 16}
            if(k == 21){xvals[which(xvals == 20)] <- 21}
            if(k == 27){xvals[which(xvals == 26)] <- 27}
        }

        print(paste(c(monk, comparisons[i],  samptimes)))

        #choose the offset length for the different plots (plots with larger spacing in the x-axis need different offsets).
        constval <- c(.45, .65, .65, .8)[which(monk == monknames)]


        if(length(tvals) > 0){
            #indsForXvals stores which points should be plotted adjacently
            indsForXvals <- as.numeric(sort(unique(unlist(strsplit(tvals, "-")))))

            #Plot the points (black dots)
            points(xvals[indsForXvals], rep(yvalToPlot, length(indsForXvals)), pch = 16, col = "black", cex = .65)
            #PLot arrows connecting them 
            plotLines(xvals[indsForXvals], yvalToPlot, plotline, constval)

            allSumm <- rbind(allSumm, c(monk, i, sum(plotline == 2, na.rm = TRUE), sum(plotline == 3, na.rm = TRUE)))
            
        }
    }
}
mtext("Week Post-Infection", side = 1, outer = TRUE, line = 2.5)

dev.off()


#Do some analysis on the number of compartmentalized comparisons (These numbers
# are presented in the results section of the text)
allSumm <- allSumm[-1,]

#RNA compartmentalized
sum(as.numeric(allSumm[allSumm[,2] <= 5,3]))
#RNA not compartmentalized
sum(as.numeric(allSumm[allSumm[,2] <= 5,4]))
#Total RNA
sum(as.numeric(allSumm[allSumm[,2] <= 5,3])) + sum(as.numeric(allSumm[allSumm[,2] <= 5,4]))

#DNA compartmentalized
sum(as.numeric(allSumm[allSumm[,2] > 5,3]))
#RNA not compartmentalized
sum(as.numeric(allSumm[allSumm[,2] > 5,4]))
#Total DNA comparisons
sum(as.numeric(allSumm[allSumm[,2] > 5,3])) + sum(as.numeric(allSumm[allSumm[,2] > 5,4]))

#GutRNA versus plasma RNA
plasmaind <- which(comparisons == "PLASMA")
gutrnaind <- which(comparisons == "GUTRNA")

#plasma compartmentalized
sum(as.numeric(allSumm[allSumm[,2] == plasmaind,3]))
#plasma not compartmentalized
sum(as.numeric(allSumm[allSumm[,2] == plasmaind,4]))
#Total plasma comparisons
sum(as.numeric(allSumm[allSumm[,2] == plasmaind,3])) + 
    sum(as.numeric(allSumm[allSumm[,2] == plasmaind,4]))

#gut compartmentalized
sum(as.numeric(allSumm[allSumm[,2] == gutrnaind,3]))
#gut not compartmentalized
sum(as.numeric(allSumm[allSumm[,2] == gutrnaind,4]))
#Total gut comparisons
sum(as.numeric(allSumm[allSumm[,2] == gutrnaind,3])) + 
    sum(as.numeric(allSumm[allSumm[,2] == gutrnaind,4]))

