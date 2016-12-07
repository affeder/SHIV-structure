#Read in data of the viral loads at each week
dat <- read.table("../dat/Vls2.csv", fill = TRUE, sep = ",", stringsAsFactors = FALSE, header = TRUE)
#Convert to numeric
for(i in 2:ncol(dat)){ 
    dat[,i] <- as.numeric(unlist(lapply(strsplit(dat[,i], ","),
                                        paste, collapse = ""))) }


#Determine the relevant week and compartmental information
allweeks <- sort(unique(weeks))[-1]
allcomps <- sort(unique(samp.loc))[-1]

#Set up lists to keep track of the percentage drug resistant
drm103list <- sapply(monknames,function(x) NULL)
drm184list <- sapply(monknames,function(x) NULL)
toAdd <- rep(NA, length(allweeks))
names(toAdd) <- paste(allweeks)
for(i in  monknames){
    drm103list[[i]] <-  toAdd
    drm184list[[i]] <-  toAdd
}
for(ids in monknames){
    for(week.ind in 1:length(allweeks)){
        #for the plasma froma  certain macaque at a certain week
        relset <- which((monkid == ids & weeks == allweeks[week.ind]) & samp.loc == "PLASMA")
        if(length(relset) > 0){
            #determine the percentage of DRMs present
            drmdist <- newDRM.breakdown(relset)
            drm103list[[ids]][week.ind] <- drmdist[1]
            drm184list[[ids]][week.ind] <- drmdist[2]
        }
    }
}


pdf("../out/graphs/VL.pdf", width = 8, height = 6)
par(mar = c(1, 1,1,0))
par(oma = c(3, 3, 1, 5))
layout(matrix(1:4, nrow = 2))
#max ylim value
maxvl <- max(dat, na.rm = TRUE)
#max xlim value

#This is indexed by the data matrix column headers (this is why the first item is NA- 
#The first column is the week from which the data was taken 
#This shows what the final plotted x value should be. Note: this is not what the final 
# timepoint for each macaque actually is
xlimmaxopts <- c(NA, 44, 29, 44, 29)
#We'll plot a line for the final time point as well
finaltimepoint <- c(NA, 44, 29, 38, 26)
cols <- brewer.pal(8, "Set2")
#This is the highlight color for the opposite axis
altcol <- cols[2]

for(monk in monknames){
    #Take the correct index for the data matrix
    i <- which(monk == names(dat))
    #Use this index to derive graphical parameters
    xlimmax <- xlimmaxopts[i]
    finaltp <- finaltimepoint[i]

    #Set up plot and background:
    plot(0, type = 'n', xlim = c(0, xlimmax), ylim = c(5, maxvl),
         log = "y", axes = FALSE, xlab = "", cex.lab = 1.5, ylab = "")
    mtext(paste(monk), 3)
    polygon(c(12, 20, 20, 12), c(.1, .1, 10^7, 10^7), col = "lightgrey", border = FALSE)
    abline(v = finaltp, col = "black", lty = "dotted")
    
    #If the macaque is T98133, plot the legend
    if(monk != "T98133"){
        polygon(c(26, finaltp, finaltp, 26), c(.1, .1, 10^7, 10^7),
                col = "lightgrey", border = FALSE)
    }else{
        legend("topleft", c("Viremia", "K103N", "M184I/V"),
               col = c("black", altcol, altcol), pch = c(16, 16, 17),
               lty = c("solid", "solid", "dashed"), box.col = "white", cex = .75)
    }

    #If the plot is on the bottom, set up the appropriate x axis
    if(monk == "A99165"){
        axis(1, cex.axis = 1.25, at = c(0, 12,20, 26, 39))
    }
    if(monk == "A01198"){
        axis(1, cex.axis = 1.25, at = c(0, 12,20, 26, 38, 44))
    }
    #If the plot is on the left, set up the y-axis
    if(monk == "T98133"| monk == "A99165"){
        eaxis(2, cex.axis = 1.25, at = 10^c(2:7))
    }

    #xToPlot is always the week in which the sample was taken
    xToPlot <- dat[,1]
    #yToPlot is macaque specific
    yToPlot <- dat[,i]

    #Plot all-non NA data as points and a line
    lines(xToPlot[!is.na(yToPlot)], yToPlot[!is.na(yToPlot)], lwd = 1.5, col = "black" )
    points(xToPlot[!is.na(yToPlot)], yToPlot[!is.na(yToPlot)], pch = 16,
           cex = 1, col = "black")

    #Now, plot the second set of axes
    par(new = T)
    plot(0, type = "n", axes = FALSE, ylim = c(-.085, 1), xlim = c(0, xlimmax),
         xlab = "", ylab = "")
    
    #First, we'll plot K103N using one set up graphical parameters
    toPlot <- drm103list[[monk]]
    toExc <- which(is.na(toPlot))
    weeksToPlot <- as.numeric(names(toPlot)[-toExc])
    pointsToPlot <- toPlot[-toExc]
    points(weeksToPlot, pointsToPlot, pch = 16, col = altcol, lty = "solid")
    lines(weeksToPlot, pointsToPlot, col = altcol, lty = "solid")

    #Second, we'll plot M184V/I using a different set of graphical parameters
    toPlot <- drm184list[[monk]]
    toExc <- which(is.na(toPlot))
    weeksToPlot <- as.numeric(names(toPlot)[-toExc])
    pointsToPlot <- toPlot[-toExc]
    points(weeksToPlot, pointsToPlot, pch = 17, col = altcol, lty = "dashed")
    lines(weeksToPlot, pointsToPlot, col = altcol, lty = "dashed")

    #If the macaque is on the right, plot the second set of y-axes
    if(monk != "T98133" & monk != "A99165"){
        axis(4, at = seq(0, 1, .25 ),
             labels = paste(seq(0, 100, by = 25), "%", sep = ""),
             las = 2, col.axis = altcol, col = altcol)
    }

    #Plot a few final graphical elements
    box()
    abline(h = -.05)
    plotBottom(monk, c(-5, -.05), -.085)
}
mtext("Percentage of Plasma vRNA with Resistance Mutation", 4,
      line = 3.5, outer = TRUE, col = altcol)
mtext("Plasma RT-SHIV RNA (copies/ml)", 2, line = 2, outer = TRUE, col = "black")
mtext("Week Post-Infection", 1, line = 2, outer = TRUE)
dev.off()
