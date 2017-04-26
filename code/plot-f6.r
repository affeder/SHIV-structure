
plotBackground <- function(monk.name, comp1){
    #Background parameter values
    pchs <- 16
    darkval <- .6

    if(monk != "T98133"){
        xlimval <- maxX[which(monk.name == monknames)]
        polygon(c(26, 26, xlimval, xlimval),  c(bbox.upper, 1 + y.padding, 1 + y.padding, bbox.upper), col = "grey90", border = NA) }
    polygon(c(12, 12, 20, 20), c(bbox.upper, 1 + y.padding, 1 + y.padding, bbox.upper), col = "grey90", border = NA)

    polygon(x = c(12 - x.padding, 12 - x.padding, xlimval + x.padding, xlimval + x.padding), y = c(bbox.upper, 1 + y.padding, 1 + y.padding, bbox.upper))

    #If it's in the top row, plot the monkey name
    if(comp1 == ordlocrna[1]){
        mtext(monk.name, line = -.18)
    }

    #If it's in the last row, plot the x axis and treatment detailsn
    if(comp1 == ordlocrna[5]){
        axis(1, at = c(12, 15, 20, 26, xlimval))
        plotBottom(monk.name, c(bbox.lower, bbox.upper), (bbox.lower - bbox.upper)+ .05)
        polygon(x = c(12 - x.padding, 12 - x.padding, xlimval + x.padding, xlimval + x.padding), y = c(bbox.upper, bbox.lower+.055, bbox.lower+.055, bbox.upper))
    }

    #If it's in the first column, plot the y axis
    if(monk.name == monknames[1]){
        axis(2, at = c(0, .5, 1))
    }

    #If it's in the last column, plot the treatment compartment name
    if(monk.name == monknames[4]){

        text(x = 49, y = .5,
             namedordlocrna[which(ordlocrna == comp1)],
              col = colvals[which(ordlocrna == comp1)], cex = 1.4, srt = -90, xpd = NA)
    }
}



#Shared graphical parameters
colvals <- brewer.pal(6, "Set2")[c(4, 3, 2, 1, 6)]
colvals <- c("#228833", "#EE6677", "#AA3377", "#4477AA", "#66CCEE")
#colvals is ordered to be the right colors to match up with ordlocrna
lwdval <- 1.5
cexval <- 1.5
maxX <- c(26, 29, 39, 44)

bbox.lower <- -.35
bbox.upper <- -.1
y.padding <- .1
x.padding <- 1.5
x.padding.plot <- .358
ylims <- c(bbox.lower+.1, 1 + y.padding)
xlims <- c(13 - x.padding - x.padding.plot)
#I introduce random jitter, so I set a seed
set.seed(9284084)

for(test in c("KST", "SM", "AMOVA")){
    dataToPlot <-  doubleddat %>% 
        filter(stattype == test)  %>%
        filter(ds == 10, !is.na(p)) %>%
        group_by(monk, ordloc1, ordloc2, week.t) %>%
        summarize(pval = mean(p < .05, na.rm = TRUE))

    pdf(paste("../out/graphs/F6-",test,".pdf", sep = ""), width = 7.5, height = 5)
    par(oma = c(3.5, 3.5, 1, 9))
    par(mar = c(0, 0, 0, 0))
    layout(matrix(1:(length(unique(dataToPlot$ordloc1)) * length(unique(dataToPlot$monk))),
                  nrow = 5))
    for(monk.name in monknames){
        for(comp1 in ordlocrna){
            plot(0, type = "n",
                 xlim = c(xlims, maxX[which(monk.name == monknames)] + x.padding +
                              x.padding.plot),
                 ylim = ylims, axes = FALSE)
            #Plot background polygons and features
            plotBackground(monk.name, comp1)
            for(comp2 in rnacomps){
                lToPlot <- dataToPlot %>%
                    filter(monk == monk.name, ordloc1 == comp1, ordloc2 == comp2)
                lToPlot <-  lToPlot %>% mutate(week.t = week.t + rnorm(nrow(lToPlot), 0, .25))
                points(lToPlot$week.t, lToPlot$pval, col = colvals[which(ordlocrna == comp2)],
                       pch = 16, cex = cexval)
                lines(lToPlot$week.t, lToPlot$pval, col = colvals[which(ordlocrna == comp2)],
                      lwd = lwdval)
            }
        }
    }
    mtext("Week Post-Infection", side = 1, outer = TRUE, line = 2.25)
    mtext("Probability of significant test", side = 2, outer = TRUE, line = 2.25)
    #Plot the legend
    legend(51.5, 3.85 ,paste(namedordlocrna), col = c(colvals), pch = 16, lty = "solid", xpd = NA, lwd = lwdval, pt.cex = 1.5, box.col = "white", title = "Compared to...")
    dev.off()
}



