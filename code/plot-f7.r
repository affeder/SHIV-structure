
halfCircle <- function(x,y,r,start=0,end=pi,nsteps=30, colval, xpdval, dashed = NULL){
    if(!is.null(colval)){
        rs <- seq(start,end,len=nsteps)
        xc <- x+r*cos(rs)
        yc <- y+r*sin(rs)
        lines(xc, yc)
        xpdval <- NA
        polygon(c(xc, rep(x, nsteps)), c(yc, rep(y, nsteps)), col = colval, xpd = xpdval)
        if(!is.null(dashed)){
            polygon(c(xc, rep(x, nsteps)), c(yc, rep(y, nsteps)), xpd = xpdval, density = 25, col = "black")
        }
    }
}

dataToPlot <-  dat %>% 
    filter(ds == 10, !is.na(p)) %>%
    group_by(monk, ordloc1, ordloc2, week.t, stattype) %>%
    summarize(pval = mean(p < .05, na.rm = TRUE))
pdf("../out/graphs/F7.pdf", height = 3.35, width = 3.75)
orderedTests <- c("SM", "KST", "AMOVA")
par(oma = c(0, 0,0,0))
par(mar = c(0,0,0,0))
rval <- .45
colorPal <- colorRampPalette(c("white", brewer.pal(3, "Set1")[1]))
cols <- colorPal(101)
namesNoRNA <- unlist(lapply(strsplit(namedordlocrna, " "), function(x){x[1]}))
plot(0, type = "n", xlim = c(-.5, 5.5), ylim = c(.5, 5.75), axes = FALSE, ylab = "", xlab = "")
polygon(y = c(0, 0, 2.5, 2.5), x = c(2.5, 6, 6, 2.5))
for(rnaind1 in 1:(length(ordlocrna))){
    for(rnaind2 in (rnaind1):(length(ordlocrna))){
        summaries <- dataToPlot %>%
            filter(ordloc1 == ordlocrna[rnaind1],
                   ordloc2 == ordlocrna[rnaind2]) %>%
                       ungroup() %>% group_by(stattype) %>%
                           summarize(avP = floor(mean(pval)* 100 + 1), numObs = n())
        rotateInd <- 0
        for(stat in orderedTests){
            pval <- summaries %>% filter(stattype == stat) %>% select(avP)
            #We want the areas of the circles to be proportional to volume. 
            #Therefore, 
            widthval <- summaries %>% filter(stattype == stat) %>%  mutate(sqnum = sqrt(numObs/base::pi)) %>% select(sqnum)* rval/2
            colval <- cols[ as.numeric(pval) ]
            halfCircle(rnaind1, rnaind2, r = as.numeric(widthval), start = 0 + rotateInd,
                       end = (2*base::pi)/3 + rotateInd, col = colval)
            rotateInd <-  rotateInd  + (2*base::pi)/3
        }
    }
}
text(1:5, 5.65, namesNoRNA)
text(.45, 1:5 , namesNoRNA, pos = 2)
rval <- .25
offset <- list("SM" = c(rval + .075, rval + .075),
               "KST" =  c(-rval - .4, 0),
               "AMOVA" = c(rval + .4, -rval - .05))
leg.x <- 3.5
leg.y <- 1.75
for(stat in orderedTests){
    halfCircle(leg.x, leg.y, r = rval, start = 0 + rotateInd,
               end = (2*base::pi)/3 + rotateInd, col = "lightgrey")
    rotateInd <-  rotateInd  + (2*base::pi)/3
    text(leg.x + offset[stat][[1]][1], leg.y +offset[stat][[1]][2], stat, cex = .85)
}
leg.y <- .85
offset <- list("SM" = c(rval + .05, rval + .0025, sqrt(2/base::pi)*.45/2, 2),
               "KST" =  c(-rval - .3, 0, sqrt(8/base::pi)*.45/2, 8),
               "AMOVA" = c(rval + .35, -rval - .05, sqrt(15/base::pi)*.45/2, 15))
for(stat in orderedTests){
    halfCircle(leg.x, leg.y, r = offset[stat][[1]][3], start = 0 + rotateInd,
               end = (2*base::pi)/3 + rotateInd, col = "lightgrey")
    rotateInd <-  rotateInd  + (2*base::pi)/3
    text(leg.x + offset[stat][[1]][1], leg.y +offset[stat][[1]][2], paste(offset[stat][[1]][4]), cex = .85)
}
scale.x <- 5
scale.y.upper <- 2.25
scale.y.lower <- .4
len <- 100
scale.len <- .15
yscale <- seq(scale.y.lower, scale.y.upper, length = len)
segments(rep(scale.x, len) -  scale.len, yscale  , rep(scale.x, len) + scale.len, yscale , col = cols, lwd = 2)
polygon(x = c(scale.x-scale.len, scale.x-scale.len, scale.x+scale.len, scale.x+scale.len),
        y = c(scale.y.upper, scale.y.lower, scale.y.lower, scale.y.upper))
text(x = scale.x - scale.len, y = c(scale.y.lower, scale.y.upper), c(0, 1), pos = 2)
segments(x0 = scale.x - scale.len, y0 = c(scale.y.lower, scale.y.upper), x1 = scale.x - scale.len - .1, y1 = c(scale.y.lower, scale.y.upper))
text(scale.x + .25, scale.y.upper - (scale.y.upper - scale.y.lower)/2 - .2,
     "P(significant test)", pos = 3, srt = -90, cex = .85)
text(x = 2.5, y = 2.25, "Test:", pos = 4, font = 2)
text(x = 2.5, y = 1.25, "n:", pos = 4, font = 2)
box()
axis(1)
axis(2)
dev.off()



#Determine the probability of significance for each test type
dataToPlot <-  dat %>% 
    filter(ds == 10, !is.na(p)) %>%
    group_by(monk, ordloc1, ordloc2, week.t, stattype) %>%
    summarize(pval = mean(p < .05, na.rm = TRUE))

dataToPlot %>% filter(ordloc1 != ordloc2) %>% ungroup() %>% group_by(stattype) %>%
    summarize(mean(pval), n())
