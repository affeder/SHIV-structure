
#Helper function that returns the frequency of 103N and 184V/I
newDRM.breakdown <- function(inds){
    freq103N <- sum(aa[inds, 103] == "N")/length(inds)
    freq184V <- length(grep("V|I", aa[inds, 184]))/length(inds)
    c(freq103N, freq184V)
}

#Make a list of all weeks and compartments to check
allweeks <- sort(unique(weeks))[-1]
allcomps <- sort(unique(samp.loc))[-1]

#Create a list by macaque of compartments versus time
oneDRMlist <- sapply(monknames,function(x) NULL)
#Helper method to add matrices to a list
toAdd <- matrix(data = NA, ncol = length(allcomps), nrow = length(allweeks))
rownames(toAdd) <- paste(allweeks)
colnames(toAdd) <- paste(allcomps)
for(i in  monknames){
    oneDRMlist[[i]] <-  toAdd
}


for(ids in monknames){
   for(week.ind in 1:length(allweeks)){
       for(comp.ind in 1:length(allcomps)){
           #Iterate over each compartment, week and macaque
            relset <- which((monkid == ids & weeks == allweeks[week.ind]) & samp.loc == allcomps[comp.ind])
            if(length(relset) > 0){
                #Determine the percentage of the compartment/time point that is drug resistant
                props <- newDRM.breakdown(relset)
                #Save the relevant number (dependent by macaque)
                if(ids == "A99039" | ids == "T98133"){
                    drmdist <- props[2]
                }else if(ids == "A99165" | ids == "A01198"){
                    drmdist <- props[1]
                }
                oneDRMlist[[ids]][week.ind, comp.ind] <- drmdist
            }
        }
    }
}



pdf("../out/graphs/percent_drugresistant.pdf", width = 5, height =4)
#Graphical parameters
layout(matrix(1:12, nrow = 4))
par(oma = c(4, 5, 2, 7))
par(mar = c(0, 1, 0 ,0))
lwd1 <- 2
lwd2 <- 2
pchr <- 15
pchd <- 16
pal <- brewer.pal(8, "Set2")

#These are compartment indexed vectors
cols <- c(pal[1], pal[1], pal[2], pal[2], pal[3], pal[3], pal[4], pal[6], pal[6])
pchvals <- c(pchd, pchr, pchd, pchr, pchd, pchr, pchr, pchd, pchr)
ltys <- c("solid", "dashed", "solid", "dashed", "solid", "dashed", "dashed", "solid", "dashed")

#This is a macaque indexed vector
maxes <- c(26,26,27,26)

plotcomps <- rbind(c("PLASMA", "PBMCDNA"), c("LNRNA", "LNDNA"), c("GUTRNA", "GUTDNA"), c("VAGRNA", "VAGDNA"))
#For each plotted macaque
for(k in c(1, 3, 2)){
    maxweeks <- maxes[k]
    minweeks <- 12
    
    for(comp.ind in 1:nrow(plotcomps)){
        #Set up the plot for each macaque/compartment and background elements
        plot(0, type = "n", xlim = c(minweeks, maxweeks), ylim = c(-.15, 1.05), xlab = "", ylab = "", axes = FALSE)
        polygon(x = c(12, 20, 20, 12), y = c(-10, -10, 10, 10),
                col = "lightgrey", border = NA)
        abline(h = c(0, 1), col = "darkgrey")
        box()

        #If this is the first plot, set up the y-axis
        if(k == 1){   axis(2, seq(0, 1, by = .25), labels = paste(seq(0, 100, by = 25), "%", sep = ""), cex.axis = 1, las = 2)  }
        #For each plot, setup the x-axis (but only once)
        if(comp.ind == 4){
            axis(1, at = c(13, 15, 20, 26),labels = rep("", 4))
            mtext(text = c(13, 15, 20, 26),at = c(12.5, 15.5, 20, 26), side = 1,cex = .75, line = 1)
        }
        #If this is a plot at the top of a column, plot the macaque name
        if(comp.ind == 1){ mtext(paste(monknames[k]), side = 3) }
        #If this is a plot at the bottom of a column, plot the treatment details
        if(comp.ind == 4){ plotBottom(monknames[k], c(-5, 0), -.1, trunc = TRUE)}

        #For both the rna and dna listed in each column         
        for(rnaOrDnaCol in c(1, 2)){

            #rnaind will hold the index that allows access to the correct column in oneDRMlist
            #and also the correct vector entry in various plotting functions
            rnaind <- which(plotcomps[comp.ind, rnaOrDnaCol] ==  colnames(oneDRMlist[[k]]))
plotcomps[comp.ind, rnaOrDnaCol]

            #Use that column index to determine which weeks we have information from
            oneDRM <- c()
            weekstoplot <- c()
            for(week.ind in 1:length(allweeks)){
                sto <- oneDRMlist[[k]][week.ind, rnaind]
                if(!is.na(sto)){
                    weekstoplot <- c(weekstoplot, rownames(oneDRMlist[[k]])[week.ind])
                    oneDRM <- c(oneDRM, oneDRMlist[[k]][week.ind, rnaind] )
                }
            }
            
            if(length(weekstoplot) > 0){
                #Plot both points and lines
                weekstoplot <- as.numeric(weekstoplot)
                lines(weekstoplot, oneDRM , col = cols[rnaind],
                      lty = ltys[rnaind], lwd = lwd1)
                points(weekstoplot, oneDRM , col = cols[rnaind],
                       pch = pchvals[rnaind])
            }
        }
    }
}
mtext("Week Post-Infection", side = 1, outer = TRUE, line = 2.5)
mtext("Frequency of drug-resistant vRNA or vDNA", side = 2, outer = TRUE, line = 3)
legend(27, 3.5, c("Gut", "LN", "PBMC", "Plasma", "Vagina", "", "DNA", "RNA"), col = c(pal[1:4], pal[6],  "white", "black", "black"), pch = c(16, 16, 16, 16, 16,  NA, NA, NA ), lty = c(NA, NA, NA, NA, NA, NA, "solid", "dashed"), xpd = NA ,lwd = lwd1, pt.cex = 1.5, box.col = "white")
dev.off()
