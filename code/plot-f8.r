kstout.time$loc <- factor(kstout.time$loc, levels = c("PLASMA", "LNRNA", "PBMCRNA",
                                               "GUTRNA", "VAGRNA", "LNDNA", "PBMCDNA",
                                               "GUTDNA", "VAGDNA"))

plotComps <- c("PLASMA", "LNRNA", "PBMCRNA", "GUTRNA", "VAGRNA", " ",
               "LNDNA", "PBMCDNA", "GUTDNA", "VAGDNA")
readPlotComps <- c("Plasma", "LN", "PBMC", "Gut", "Vagina", " ",
               "LN", "PBMC", "Gut", "Vagina")



datToPlot <-kstout.time %>%
    filter(!is.na(p)) %>%
        group_by(monk, week1, week2, loc) %>%
            summarize(pval = mean(p < .05)) 

#graphical parameters
colorPal <- colorRampPalette(c("white", brewer.pal(3, "Set1")[1]))
cols <- colorPal(101)


plotBackgroundTime <- function(monk.name){
    pchs <- 16
    darkval <- .6

    if(monk != "T98133"){
        xlimval <- maxX[which(monk.name == monknames)]
        polygon(c(26, 26, xlimval, xlimval),  c(bbox.upper, 11, 11, bbox.upper), col = "grey90", border = NA) }
    polygon(c(12, 12, 20, 20), c(bbox.upper, 11, 11, bbox.upper), col = "grey90", border = NA)

    #plot the monkey name
    mtext(monk.name, line = 0)
    plotBottom(monk.name, c(bbox.lower, bbox.upper), 0)
#    polygon(x = c(12, 12 , xlimval , xlimval), y = c(bbox.upper, bbox.lower, bbox.lower, bbox.upper))
    polygon(c(0, 0, 50, 50), c(bbox.upper, bbox.lower, bbox.lower, bbox.upper))
    axis(1, at = c(12, 15, 20, 26, xlimval))

    #If it's in the first column, plot the y axis
    if(monk.name == monknames[1]){
        axis(2, at = c(1:5, 7:10), labels = readPlotComps[-6], las = 1)
    }

    box()
    
}


pdf("../out/graphs/F8.pdf", height = 2.5, width = 7.5)
bbox.upper <- .5
bbox.lower <- -.5
maxX <- c(26, 29, 38, 44)
lengthval <- .03
lwdval <- 3
par(mar = c(.5,0, .5, .5))
par(oma = c(3, 6, 1.5, 3))
offsets <- c(.5, .58, 1.05, 1.25)
layout(matrix(1:4, nrow = 1))
for(monk.name in monknames){
    xlims = c(12, maxX[which(monk.name == monknames)]  )
    xvals <- c(13, 15, 20, 26, xlims[2])
    offset <- offsets[which(monk.name == monknames)]
    #Plot values
    plot(0, type = "n", xlim = xlims, ylim = c(bbox.lower+.5, 10), axes = FALSE)
    plotBackgroundTime(monk.name)
    abline(h = c(1:5, 7:10), lwd = .25, col = "grey")
    for(comp in plotComps){
        tmpxvals <- xvals
        lToPlot <- datToPlot %>% filter(monk == monk.name, loc == comp )
        if(nrow(lToPlot) > 0){
            relNames <- table(inf[inf$samp.loc == comp & monkid == monk.name, 'weeks'])
            if(sum(grepl('12', names(relNames)))){ tmpxvals[which(tmpxvals == 13)] <- 12 }
            if(sum(grepl('16', names(relNames)))){ tmpxvals[which(tmpxvals == 15)] <- 16 }
            if(sum(grepl('21', names(relNames)))){ tmpxvals[which(tmpxvals == 20)] <- 21 }
            if(sum(grepl('27', names(relNames)))){ tmpxvals[which(tmpxvals == 26)] <- 27 }
            if(sum(grepl('40', names(relNames)))){ tmpxvals[which(tmpxvals == 44)] <- 40 }
            x0vals <- tmpxvals[lToPlot$week1] + offset
            x1vals <- tmpxvals[lToPlot$week2] - offset
            plotCols <- cols[floor(lToPlot$pval*100) + 1]
            arrows(x0 = x0vals, x1 = x1vals, y0 = which(plotComps == comp), code = 3, col = "black", lwd = lwdval + 1, length = lengthval)
            arrows(x0 = x0vals, x1 = x1vals, y0 = which(plotComps == comp), code = 3, col = plotCols, lwd = lwdval, length = lengthval)
            pointsToPlot <- unique(c(tmpxvals[lToPlot$week1], tmpxvals[lToPlot$week2]))
            points(pointsToPlot,
                   rep(which(plotComps == comp), length(pointsToPlot)), pch = 16)
        }
    }
}
mtext("Week Post-Infection", side = 1, outer = TRUE, line = 2)
mtext(c("vRNA", "vDNA"),2,
                      las = 0, at = c(.35, .8), line = 4.5, outer = TRUE, xpd = NA)
xvala <- -109.5
arrows(xvala, .5 , xvala, 5.75, length = 0, col = "black", xpd = NA)
arrows(xvala, 6.75, xvala, 10.25, length = 0, col = "black", xpd = NA)
#plot legend
scale.x <- 50
scale.y.upper <- 10
scale.y.lower <- .5
len <- 100
scale.len <- 1
tinyoffset <- .075
yscale <- seq(scale.y.lower, scale.y.upper, length = len)
segments(rep(scale.x, len) -  scale.len, yscale  , rep(scale.x, len) + scale.len, yscale , col = cols, lwd =2, xpd = NA)
polygon(x = c(scale.x-scale.len-tinyoffset, scale.x-scale.len-tinyoffset,
            scale.x+scale.len+tinyoffset, scale.x+scale.len+tinyoffset),
        y = c(scale.y.upper+tinyoffset, scale.y.lower-tinyoffset,
            scale.y.lower-tinyoffset, scale.y.upper+tinyoffset), xpd =NA)
text(x = scale.x - scale.len, y = c(scale.y.lower, scale.y.upper), c(0, 1), pos = 2, xpd = NA)
segments(x0 = scale.x - scale.len - tinyoffset, y0 = c(scale.y.lower, scale.y.upper), x1 = scale.x - scale.len - 1 - tinyoffset, y1 = c(scale.y.lower, scale.y.upper), xpd = NA)
text(scale.x + 2.5, scale.y.upper - (scale.y.upper - scale.y.lower)/2 - .2,
     "P(significantly compartmentalized)", pos = 3, srt = -90, xpd = NA)
dev.off()





                           
toTest <- kstout.time %>%
    filter(!is.na(p)) %>%
        group_by(monk, week1, week2, loc) %>%
            summarize(pval = mean(p < .05)) %>%
                mutate(DNA = ifelse(grepl("DNA", loc), 1, 0)) %>%
                    ungroup()

#DNA versus RNA test
toTest %>% group_by(DNA) %>%
    summarize(mean(pval))

DNAvals <- toTest %>% filter(DNA == 1)  %>% select(pval)
RNAvals <- toTest %>% filter(DNA == 0)  %>% select(pval)

DNAvals <- c(DNAvals)$pval
RNAvals <- c(RNAvals)$pval
t.test(DNAvals, RNAvals)

#Gut vRNA versus plasma test
gutvplasmaComp <- toTest %>% mutate(plasma = ifelse(grepl("PLASMA", loc), 1, 0), gutrna = ifelse(grepl("GUTRNA", loc), 1, 0))

plasmavals <- gutvplasmaComp %>% filter(plasma == 1)  %>% select(pval)
gutrnavals <- gutvplasmaComp %>% filter(gutrna == 1)  %>% select(pval)

plasmavals <- c(plasmavals)$pval
gutrnavals <- c(gutrnavals)$pval
t.test(plasmavals, gutrnavals)



diffsto <- foreach(monk.name = monknames, .combine = 'rbind') %do% {
    xlims = c(12, maxX[which(monk.name == monknames)]  )
    xvals <- c(13, 15, 20, 26, xlims[2])
    foreach(comp = plotComps, .combine = 'rbind')%do%{
        tmpxvals <- xvals
        lToPlot <- datToPlot %>% filter(monk == monk.name, loc == comp )
        if(nrow(lToPlot) > 0){
            relNames <- table(inf[inf$samp.loc == comp & monkid == monk.name, 'weeks'])
            if(sum(grepl('12', names(relNames)))){ tmpxvals[which(tmpxvals == 13)] <- 12 }
            if(sum(grepl('16', names(relNames)))){ tmpxvals[which(tmpxvals == 15)] <- 16 }
            if(sum(grepl('21', names(relNames)))){ tmpxvals[which(tmpxvals == 20)] <- 21 }
            if(sum(grepl('27', names(relNames)))){ tmpxvals[which(tmpxvals == 26)] <- 27 }
            if(sum(grepl('40', names(relNames)))){ tmpxvals[which(tmpxvals == 44)] <- 40 }
        }
        x0vals <- tmpxvals[lToPlot$week1] 
        x1vals <- tmpxvals[lToPlot$week2]
        diff <- x1vals - x0vals
        cbind(rep(monk.name, length(diff)), rep(comp, length(diff)), diff)
    }
}

diffs <- tbl_df(diffsto)
colnames(diffs) <- c("monk", "comp", "diff")
diffs <- diffs %>% mutate(diff = as.numeric(diff), dna = ifelse(grepl("DNA",comp),1, 0))

diffs %>% group_by(dna) %>%
    summarize(mean(diff))

rnadiffs <- diffs %>% filter(dna == 0)  %>% select(diff)
dnadiffs <- diffs %>% filter(dna == 1)  %>% select(diff)

rnadiffs <- c(rnadiffs)$diff
dnadiffs <- c(dnadiffs)$diff

t.test(rnadiffs,dnadiffs)
