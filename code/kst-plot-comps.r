

#The first half of the script plots boolean significance values of Kst by compartment, macaque and timepoint
#The second of the script plots percentage of values significant between pairs of compartments collapsed by time and macaque

#Choose graphical values
lwdval <- 2
ltys <- "solid"
pchs <- 16
darkval <- .6
#These largey affect how close lines appear to each other and can be tweaked 
# for different visual effect
offsetv <- .35
offsetv2 <- .08
#tval and bval represent the "top" offset for significant points and "bottom" offset
# for non-significant points
tval <- .5
bval <- -.02

#because the comparisons are listed in a specific order, it makes sense to just
#enter all of the comparisons manually
comparisons <- c("PBMCRNA-GUTRNA","VAGRNA-GUTRNA","PLASMA-GUTRNA","LNRNA-GUTRNA",
                 "VAGRNA-LNRNA",'VAGRNA-PLASMA',"VAGRNA-PBMCRNA","VAGRNA-GUTRNA",
                 "PLASMA-LNRNA",'PBMCRNA-PLASMA','VAGRNA-PLASMA',"PLASMA-GUTRNA",
                 "PBMCRNA-LNRNA","PLASMA-LNRNA","VAGRNA-LNRNA","LNRNA-GUTRNA",
                 "PBMCRNA-LNRNA",'PBMCRNA-PLASMA',"VAGRNA-PBMCRNA","PBMCRNA-GUTRNA")

#Set up indices that will allow us to assign the 
# correct colors to the correct comparisons easily
gutind <- 1
lnind <- 2
pbmcind <- 3
plasmaind <- 4
vagind <- 6

#set the colors to match the comparisons
colvals <- rep(brewer.pal(6, "Set2"), 2)
cols <- colvals[c(pbmcind, vagind, plasmaind, lnind,
                  lnind, plasmaind, pbmcind, gutind,
                  lnind, pbmcind, vagind, gutind,
                  pbmcind, plasmaind, vagind, gutind,
                  lnind, plasmaind, vagind, gutind)]

#Set up plot the print
pdf("../out/graphs/bettercomp-new.pdf", height = 6, width = 8)

#Write layout values
par(mar = c(0, .5, 0, 0))
par(oma = c(4, 4, 2, 11.5))
layout(matrix(1:(length(monknames)*5), nrow = 5))

#Each macaque will need slightly different x-axis values - this controls that
approxx <- rbind(c(13, 15, 20, 26, 26), c(13, 15, 20, 26, 29), c(13, 15, 20, 26, 38), c(13, 15, 20, 26, 44))

#allSumm keeps track of some information about the comparisons that will be used for summary
allSumm <- c("NA", 0, "NA", 0)

for(monk in monknames){

    #We don't want to plot the 'fifth' time point for macaque T98133, because it doesn't have one. However, for convenience reasons, it's easier to just exclude that final point.
    upperxlim = 5
    if(monk == "T98133"){ upperxlim = 4 }

    #set the xvals based on which macaque we are plotting
    xvals <- approxx[which(monk == monknames), ]

    #iterate over all the comparisons
    for(i in 1:length(comparisons)){

        #Every four comparisons, we want to make a new plot
#        if(sum(i == c(1, 5, 9, 13, 17)) > 0){
        if(i %% 4 == 1){

            #Create a blank plot on which we will plot points and lines
            plot(0, type = "n", xlim = c(12, xvals[upperxlim]+.25), ylim = c(-.35, .8), axes = FALSE)

            #Add axes depending on where in the plot layout it appears 
            if(monk == "T98133"){
                axis(2, at = c(bval, tval)+.1, c("N.S.", "Sig"), las = 2)
            }
            if(i == 1){
                mtext(side = 3, paste(monk))
            }
            if(i == 17){
                axis(1, at = xvals[1:upperxlim])
            }

            #Plot background boxes, including dark boxes for treatment timing
            if(monk != "T98133"){polygon(c(26, 26, xvals[5], xvals[5]), c(-1, 2, 2, -1), col = "grey90", border = NA) }
            polygon(c(12, 12, 20, 20), c(-1, 2, 2, -1), col = "grey90", border = NA)
            polygon(c(0, 0, 50, 50), c(tval - offsetv2, tval + offsetv, tval + offsetv, tval- offsetv2), col = rgb(0,0,0,darkval))
            polygon(c(0, 0, 50, 50), c(0- offsetv2, 0 + offsetv, 0 + offsetv, 0- offsetv2), col = rgb(0,0,0,darkval))

            box()

            #If the rightmost column is being plotted, add text describing what the focal 
            # focal compartment is
            if(monk == "A01198"){
                texttop <- ""
                coltp <- ""
                if(i == 1){
                    texttop <- "Gut vRNA"
                    coltp <- colvals[gutind]
                }
                if(i == 5){
                    texttop <- "Vagina vRNA"
                    coltp <- colvals[vagind]
                }
                if(i == 9){
                    texttop <- "Plasma vRNA"
                    coltp <- colvals[plasmaind]
                }
                if(i == 13){
                    texttop <- "LN vRNA"
                    coltp <- colvals[lnind]
                }
                if(i == 17){
                    texttop <- "PBMC vRNA"
                    coltp <- colvals[pbmcind]
                }
                if(coltp != ""){ mtext(texttop, 4, las = 3, line = 1, col = coltp)}
            }
        }

        plotline <- c()
        for(t in 1:upperxlim){
            #determine, for that macaque, compartment and timepoint the values of 
            #the actual data
            statval <- kstlist[[monk]][[paste(t)]][[comparisons[i]]]
            #the randomized ksts
            nullvals <- kstlist.rand[[monk]][[paste(t)]][[comparisons[i]]]
            #and how they compare to each other
            pval <- 1-sum(statval >= nullvals)/length(nullvals)

            #if there's no data, don't plot a point
            if(is.nan(pval)){
                plotline[t] <- NA
            }else{
                #if there is data, the y-axis should be determined by whether or not the
                #p-value is below the cutoff
                if(pval < pvalcut){
                    plotline[t] <- tval 
                    allSumm <- rbind(allSumm, c(monk, t, comparisons[i], 1))
                }else{
                    plotline[t] <- bval 
                    allSumm <- rbind(allSumm, c(monk, t, comparisons[i], 0))
                }
            }
        }

        #reorganize p vals (we don't want points to overlap with each other, so we offset them
        # slightly
        plotline <- plotline + (i %% 5)/15

        #Plot all points (NAs handled as non-plots)
        points(xvals[1:upperxlim], plotline, col = cols[i], pch = pchs, cex = 1.5)
        #Plot a line connecting the points
        
        #If there are missing values, remove them and plot a non-NA line
        if(sum(is.na(plotline) > 0)){
            lines((xvals[1:upperxlim])[!is.na(plotline)], plotline[!is.na(plotline)], col = cols[i], lty = ltys, lwd =lwdval)
            #otherwise, just plot the line as normal
        }else{
            lines(xvals[1:upperxlim], plotline, col = cols[i], lty = ltys, lwd = lwdval)
        }
    }
    
    #plot extra features
    abline(h = -.2)
    #plotBottom is a general function that makes the treatment plotting at the bottom
    plotBottom(monk, c(-1, -.2), -.3)
}
mtext("Week Post-Infection", side = 1, outer = TRUE, line = 2.5)
#Plot the legend
legend(55, 3.85 ,c("Gut vRNA", "Vagina vRNA", "Plasma vRNA", "LN vRNA",  "PBMC vRNA" ), col = c(colvals[c(gutind, vagind, plasmaind,lnind, pbmcind  )]), pch = 16, lty = "solid", xpd = NA, lwd = lwdval, pt.cex = 1.5, box.col = "white", title = "Compared to...")
dev.off()





#For this plot, we are going to flatten the data across all time points and macaques
#and only summarize by the compartments being compared

#ocomp = ordered comparisons (so that we access the right entries in our list)
# i.e., GUTNA-VAGRNA will return something, but VAGRNA-GUTRNA will not
ocomps <- c("VAGRNA", "PBMCRNA", "PLASMA","LNRNA", "GUTRNA" )

compList <- list()
#iterate over all pairs of comparisons
for(rnaind1 in 1:(length(ocomps) - 1)){
    for(rnaind2 in (rnaind1+1):(length(ocomps))){

        #For a particular combination, we're going to create a new list
        compName <- paste(ocomps[rnaind1], ocomps[rnaind2], sep = "-")

        #and keep track of the total number of significant differences ('sig') and 
        # the total number of comparisons ('tot')

        compList[[compName]] <- list(sig = 0, tot = 0)

        #for each macaque and each timepoint
        for(monk in monknames){
            for(t in 1:5){

                #Store the truth and the random distribution
                statval <- kstlist[[monk]][[paste(t)]][[compName]]
                nullval <- kstlist.rand[[monk]][[paste(t)]][[compName]]

                if(!is.null(statval)){
                    #If the comparison can be made, increment the "tot" count
                    compList[[compName]]$tot <- compList[[compName]]$tot + 1
                    if(sum(statval > nullval)/length(nullval) > (1 - pvalcut)){
                        #If truth is below the p-value threshold, increment the "sig" count
                        compList[[compName]]$sig <- compList[[compName]]$sig + 1
                    }  
                }
            }
        }
    }
}

#These are helper functions for plotting the squares used in the figure
ttpx <- function(x){
    return(c(x + .5, x-.5, x -.5, x + .5))
}
ttpy <- function(x){
    return(c(x + .5, x+.5, x -.5, x - .5))
}
#Set up the plot
pdf("../out/graphs/compoverall.pdf", width = 5, height = 4.5)
par(mar = c(0,0,0,0))
par(oma = c(0, 0,0,0))
#Set up a vector of color intensities 
cols <- tcol(brewer.pal(3, "Set1")[1], seq(100, 0, by = -1))
plot(0, type = "n", xlim = c(-0.5, 5.5), ylim = c(.5, 6), axes = FALSE)
#Let's list the compartments in a different order to be plotted in
ocompsnew <-  c( "GUTRNA", "VAGRNA", "PBMCRNA", "PLASMA","LNRNA" )
for(rnaind1 in 1:(length(ocompsnew))){
    for(rnaind2 in (rnaind1):(length(ocompsnew))){

        #Access the values computed above 
        compName <- paste(ocompsnew[rnaind1], ocompsnew[rnaind2], sep = "-")
        if(rnaind1 != rnaind2 & is.null(compList[[compName]] )){
            compName <- paste(ocompsnew[rnaind2], ocompsnew[rnaind1], sep = "-")
        }

        vals <- compList[[compName]]
        sigval <- vals$sig
        totval <- vals$tot

        #If there's no data available (i.e., the diagonals), it should be plotted grey
        if(is.null(totval)){
            colval <- "grey"
        }else{
            #Otherwise, its intensity should be proportional to the percentage significant
            colval <- cols[round(sigval/totval*100)]
        }
        #Plot a square using the color intensity 
        polygon(ttpx(rnaind1), ttpy(rnaind2), col = colval)

        #If there's no data available, let's also plot an X across the box
        if(is.null(totval)){
            rx <- range(ttpx(rnaind1))
            ry <- range(ttpy(rnaind1))
            arrows(rx[1], ry[1], rx[2], ry[2], length = 0)
            arrows(rx[2], ry[1], rx[1], ry[2], length = 0)
        }else{

            #If data is available, plot the percentage of comparisons significant
            text(rnaind1, rnaind2, paste(round(sigval/totval*100), "%", sep = ""))
        }
    }
}
#Plot labels
ocomplabs <- c("Gut", "Vagina", "PBMC", "Plasma", "Lymph\nNode")
text(rep(1, 5) - .65, 1:5, ocomplabs, adj = 1)
text((1:5) , rep(5.5, 5), ocomplabs, pos = 3)
dev.off()

#We'll now use allSumm to compute some summary stats
#delete the first row (which was tmp data)
allSumm <- allSumm[-1,]
colnames(allSumm) <- c("monk", "time", "comp", "sig")
#Delete duplicate entries (every entry is added twice because of how the figure is 
#plotted
allSumm <- (allSumm[duplicated(allSumm),])
#What percentage of tests were significantly compartmentalized?
table(allSumm[,4])/length(allSumm[,4])

