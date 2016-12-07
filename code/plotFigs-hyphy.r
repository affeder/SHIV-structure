require(RColorBrewer)

dat <- read.table("../out/SM.allresults.txt", stringsAsFactors = FALSE, fill = TRUE)

allvals <- dat[, 8]
bigenough <- apply(dat[,c(5,6)] > 3, 1, function(x){ sum(x == TRUE) == 2 })
allvals <- allvals[bigenough]


lenval <- length(allvals)
cm <- sum(1/(1:lenval))
correct <- .05/(lenval*cm)*(1:lenval)
rejInd <- max(which(correct > sort(allvals)))
pvalcut <- correct[rejInd]


vagind <- 6
pbmcind <- 3
plasmaind <- 4
lnind <- 2
gutind <- 1
lwdval <- 1
pdf("../out/graphs/bettercomp-new-SM.pdf", height = 6, width = 8)
inputvals <- monknames
par(mar = c(0, .5, 0, 0))
par(oma = c(4, 4, 2, 11.5))
layout(matrix(1:(length(inputvals)*5), nrow = 5))
approxx <- rbind(c(13, 15, 20, 26, 26),
                 c(13, 15, 20, 26, 29),
                 c(13, 15, 20, 26, 38),
                 c(13, 15, 20, 26, 44))
offsetv <- .35
offsetv2 <- .08
tval <- .5
bval <- -.02
for(monk in inputvals){
    upperxlim = 4
    if(monk != "T98133"){ upperxlim = 5 }
    comparisons <- c("PBMC-GUT", "VAG-GUT", "PLASMA-GUT", "LN-GUT",
                     "VAG-LN", 'VAG-PLASMA', "VAG-PBMC", "VAG-GUT",
                     "PLASMA-LN", 'PBMC-PLASMA', 'VAG-PLASMA',  "PLASMA-GUT",
                     "PBMC-LN", "PLASMA-LN", "VAG-LN", "LN-GUT",
                     "PBMC-LN", 'PBMC-PLASMA', "VAG-PBMC", "PBMC-GUT")
    ltys <- "solid"
    pchs <- 16
    colvals <- rep(brewer.pal(6, "Set2"), 2)
    cols <- colvals[c(pbmcind, vagind, plasmaind, lnind,
                      lnind, plasmaind, pbmcind, gutind,
                      lnind, pbmcind, vagind, gutind,
                      pbmcind, plasmaind, vagind, gutind,
                      lnind, plasmaind, vagind, gutind)]
    xvals <- approxx[which(monk == monknames), ]
    for(i in 1:length(comparisons)){
        if(sum(i == c(1, 5, 9, 13, 17)) > 0){
            plot(0, type = "n", xlim = c(12, xvals[5]+.25),
                 ylim = c(-.35, .8), axes = FALSE)
            if(monk == "T98133"){
                axis(2, at = c(bval, tval)+.1, c("N.S.", "Sig"), las = 2)
            }
            if(i == 1){
                mtext(side = 3, paste(monk))
            }
            if(i == 17){
                axis(1, at = xvals[1:upperxlim])
            }
            if(monk != "T98133"){polygon(c(26, 26, xvals[5], xvals[5]), c(-1, 2, 2, -1), col = "grey90", border = NA) }
            polygon(c(12, 12, 20, 20), c(-1, 2, 2, -1), col = "grey90", border = NA)
            darkval <- .6
            polygon(c(0, 0, 50, 50), c(tval - offsetv2, tval + offsetv, tval + offsetv, tval- offsetv2), col = rgb(0,0,0,darkval))
            polygon(c(0, 0, 50, 50), c(0- offsetv2, 0 + offsetv, 0 + offsetv, 0- offsetv2), col = rgb(0,0,0,darkval))
            box()
        }
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
        plotline <- c()
        for(t in 1:upperxlim){
            possibleMatches <- dat[intersect(which(dat[,1] == monk), which(dat[,2] == t)),3:4]
            opts <- strsplit(comparisons[i], "-")[[1]]
            relDatInd <- as.numeric(names(which(apply(possibleMatches, 1, function(x){ (sum(grepl(opts[1], x)) > 0) & (sum(grepl(opts[2], x)) > 0) }) == TRUE)))
            pval <- dat[relDatInd ,8]
            if(sum(dat[relDatInd,c(5,6)] >= 3) == 2){
                if(pval == "NULL"){
                    plotline[t] <- NA
                }else{
                    if(pval < pvalcut){
                        plotline[t] <- tval #1
                    }else{
                        plotline[t] <- bval #-.02
                    }
                }
            }else{
                plotline[t] <- NA
            }
        }
#reorganize p vals
        plotline <- plotline + (i %% 5)/15
        lwdval <- 2
        lines(xvals[1:upperxlim], plotline, col = cols[i], lty = ltys, lwd = lwdval)
        if(sum(is.na(plotline) > 0)){
            lines((xvals[1:upperxlim])[!is.na(plotline)],
                  plotline[!is.na(plotline)], col = cols[i], lty = ltys, lwd =lwdval)
        }
        points(xvals[1:upperxlim], plotline, col = cols[i], pch = pchs, cex = 1.5)
    }
    abline(h = -.2)
    plotBottom(monk, c(-1, -.2), -.3)
}
#mtext("p-value", side = 2, outer = TRUE, line = 2.5)
mtext("Week Post-Infection", side = 1, outer = TRUE, line = 2.5)
legend(55, 3.85 ,c("Gut vRNA", "Vagina vRNA", "Plasma vRNA", "LN vRNA",  "PBMC vRNA" ), col = c(colvals[c(gutind, vagind, plasmaind,lnind, pbmcind  )]), pch = 16, lty = "solid", xpd = NA, lwd = lwdval, pt.cex = 1.5, box.col = "white", title = "Compared to...")
dev.off()



ttpx <- function(x){
    return(c(x + .5, x-.5, x -.5, x + .5))
}
ttpy <- function(x){
    return(c(x + .5, x+.5, x -.5, x - .5))
}





colnames(dat) <- c("monkid", "week", "rna1", "rna2", "n1",
                   "n2", "estmig", "p", "q25", "q50", "q75")


#For this plot, we are going to flatten the data across all time points and macaques
#and only summarize by the compartments being compared

#ocomp = ordered comparisons (so that we access the right entries in our list)
# i.e., GUTNA-VAGRNA will return something, but VAGRNA-GUTRNA will not
ocomps <- c("GUTRNA" , "VAGRNA", "PBMCRNA", "PLASMA","LNRNA")
shortened.comps <- c("GUT", "VAG", "PBMC", "PLASMA", "LN")

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
                timemonkinds <- intersect(which(dat$week == t),
                                          which(dat$monk == monk))
                #this handles the rna compartments being able to occur in either order
                compinds <- union(
                    intersect(which(dat$rna1 == shortened.comps[rnaind1]),
                                which(dat$rna2 == shortened.comps[rnaind2])), 
                      intersect(which(dat$rna1 == shortened.comps[rnaind2]),
                                which(dat$rna2 == shortened.comps[rnaind1]))
                )
                relind <- intersect(timemonkinds, compinds)
                #Store the pvalue
                pval <- dat[relind, 'p']
                #If the compartment is big enough
                if(sum(dat[relind, c('n1', 'n2')] >=3 ) == 2){
                    #and the p-value is not null
                    if(pval != "NULL"){
                        pval <- as.numeric(pval)
                                        #If the pvalue is not null, increment the total
                        compList[[compName]]$tot <- compList[[compName]]$tot + 1
                        if(pval < pvalcut ){
                                        #If truth is below the p-value threshold, increment the "sig" count
                            compList[[compName]]$sig <- compList[[compName]]$sig + 1
                        }  
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

pdf("../out/graphs/compoverall-SM.pdf", width = 5, height = 4.5)
par(mar = c(0,0,0,0))
par(oma = c(0, 0,0,0))
#Set up a vector of color intensities 
cols <- tcol(brewer.pal(3, "Set1")[1], seq(100, 0, by = -1))
plot(0, type = "n", xlim = c(-0.5, 5.5), ylim = c(.5, 6), axes = FALSE)
for(rnaind1 in 1:(length(ocomps))){
    for(rnaind2 in (rnaind1):(length(ocomps))){
        #Access the values computed above 
        compName <- paste(ocomps[rnaind1], ocomps[rnaind2], sep = "-")
        vals <- compList[[paste(ocomps[rnaind1], ocomps[rnaind2], sep = "-")]]
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


























#This plot is for the supplemental figure showing different values of Kst
storemat <- c(NA, NA, NA, NA, NA, NA)
ocomps <- c("VAGRNA", "PBMCRNA", "PLASMA","LNRNA", "GUTRNA" )
compList <- list()
#iterate over all pairs of comparisons
for(rnaind1 in 1:(length(ocomps) - 1)){
    for(rnaind2 in (rnaind1+1):(length(ocomps))){
        compName <- paste(ocomps[rnaind1], ocomps[rnaind2], sep = "-")
        for(monk in monknames){
            for(t in 1:5){
                monktimes <- which(monkid == monk & weekind == t)
                minsize <- min(c(
                    length(intersect(monktimes, which(samp.loc == ocomps[rnaind1]))),
                    length(intersect(monktimes, which(samp.loc == ocomps[rnaind2]))))
                               )
                statval <- kstlist[[monk]][[paste(t)]][[compName]]
                nullval <- kstlist.rand[[monk]][[paste(t)]][[compName]]
                pvalval <- NA
                if(!is.null(statval)){
                    pvalval <- 1 - sum(statval > nullval)/length(nullval)
                }
                storemat <- rbind(storemat, c(monk, t, ocomps[rnaind1], ocomps[rnaind2], minsize, pvalval))
            }
        }
    }
}


opac <- 0
pdf("../out/graphs/sm-v-kst.pdf", height = 4.5, width = 8)
smcols <- rep(tcol("black", opac), nrow(dat))
smcols[which(dat[,3] == "PLASMA" & dat[,4] == "LN")] <-tcol("red",opac)
smcols[which(dat[,3] == "LN" & dat[,4] == "VAG")] <-tcol("blue",opac)
kstcols <- rep(tcol("black", opac), nrow(storemat))
kstcols[which(storemat[,3] == "PLASMA" & storemat[,4] == "LNRNA")] <- tcol("red", opac)
kstcols[which(storemat[,3] == "VAGRNA" & storemat[,4] == "LNRNA")] <- tcol("blue", opac)
xjitval <- .1
#shared functionality
layout(matrix(1:2, nrow = 1, byrow = TRUE))
par(mar = c(4, 4, 2, .5))
par(oma = c(0, 0,0, 0))
minsize <- apply(dat[,5:6], 1, min)
pvals <- as.numeric(dat[, 8])
pvals[pvals == 0] <- .00001
#Plot 1: Slatkin-Maddison, highlight LN versus plasma comparison
enough <- which(minsize >= 3)
xjit <- rnorm(length(enough), 0, xjitval)
yjit <- rnorm(length(enough), 0, .000001)
plot(minsize[enough] + xjit, pvals[enough] + yjit, log = "y", axes = FALSE, xlab = "Minimum sample size", ylab = "p-value", main = expression("Slatkin-Maddison"), pch =16, col = smcols[enough], cex = .75)
#Replort the red and blue points
points((minsize[enough] + xjit)[which(dat[enough,3] == "PLASMA" & dat[enough,4] == "LN")], (pvals[enough] + yjit)[which(dat[enough,3] == "PLASMA" & dat[enough,4] == "LN")], pch = 16, col = tcol("red", opac), cex = .75)
points((minsize[enough] + xjit)[which(dat[enough,3] == "LN" & dat[enough,4] == "VAG")], (pvals[enough] + yjit)[which(dat[enough,3] == "LN" & dat[enough,4] == "VAG")], pch = 16, col = tcol("blue", opac), cex = .75)
eaxis(1)
eaxis(2)
abline(h = pvalcut, col = "grey", lty = "dashed")
legend("topright", c("LN-Plasma comparison", "LN-Vagina comparison", "All other comparisons"), col = c(tcol("red", opac), tcol("blue", opac), tcol("black", opac)), pch = 16, pt.cex = .75, cex = .75, box.col = "white")
#Plot 2:  Kst, highlighting LN versus plasma comparison
enough <- which(as.numeric(storemat[,5]) >= 3)
xjit <- rnorm(length(enough), 0, xjitval)
yjit <- rnorm(length(enough), 0, .000001)
storemat[storemat[,6] == "0",6] <- .00001
plot(as.numeric(storemat[,5][enough])+xjit, as.numeric(storemat[,6][enough]) + yjit, log = "y", axes = FALSE, xlab = "Minimum sample size", ylab = "p-value", main = expression('K'[ST]), pch =16, cex = .75, col = kstcols[enough])
#Then, plot the red and blue points on top so they won't be missed
points((as.numeric(storemat[,5][enough])+xjit)[which(storemat[enough,3] == "PLASMA" & storemat[enough,4] == "LNRNA")] ,(as.numeric(storemat[,6][enough])+yjit)[which(storemat[enough,3] == "PLASMA" & storemat[enough,4] == "LNRNA")] , col = tcol("red",opac), pch = 16, cex = .75)
points((as.numeric(storemat[,5][enough])+xjit)[which(storemat[enough,3] == "VAGRNA" & storemat[enough,4] == "LNRNA")] ,(as.numeric(storemat[,6][enough])+yjit)[which(storemat[enough,3] == "VAGRNA" & storemat[enough,4] == "LNRNA")] , col = tcol("blue",opac), pch = 16, cex = .75)
eaxis(1)
eaxis(2)
abline(h = pvalcut, col = "grey", lty = "dashed")
dev.off()



