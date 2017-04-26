possibleGtoAs <- length(which(cons == "G"))
rnainds <- grep("PLASMA|RNA", inf$samp.loc)
dnainds <- grep("DNA", inf$samp.loc)
gtoaprop <- apply(nuc[,which(cons == "G")], 1, function(x){length(which(x == "A"))})/possibleGtoAs
layout(matrix(1:2, nrow = 1))
hist(gtoaprop, xlab = "Proportion of G->As (observed out of possible)", main = "")
hist(gtoaprop[rnainds],  add = TRUE, col = rgb(0, 0, 1, .2))
hist(gtoaprop[dnainds], add = TRUE, col = rgb(1, 0, 0, .2))
legend("topright", c("All", "RNA", "DNA"),fill = c("white", rgb(0, 0, 1, .2), rgb(1, 0, 0, .2)), pch = 20)
cdf <- quantile(gtoaprop, seq(0, 1, by = .005))
cdf.r <- quantile(gtoaprop[rnainds], seq(0, 1, by = .005))
cdf.d <- quantile(gtoaprop[dnainds], seq(0, 1, by = .005))
xvals <- as.numeric(gsub("%", "", (names(cdf))))
plot(xvals, cdf, type = "l", xlab = "quantile")
lines(xvals, cdf.r, type = "l", xlab = "quantile", col = rgb(0, 0, 1, .7))
lines(xvals, cdf.d, type = "l", xlab = "quantile", col = rgb(1, 0, 0, .7))
legend("topleft", c("All", "RNA", "DNA"),col = c("black", rgb(0, 0, 1, .7), rgb(1, 0, 0, .7)), lty = "solid")

prop <- .02
fisher.test(cbind( c(sum(gtoaprop[dnainds] >= prop), length(dnainds)),
      c(sum(gtoaprop[rnainds] >= prop), length(rnainds))))


prop <- .02
fisher.test(cbind( c(sum(gtoaprop[dnainds] >= prop), length(dnainds)),
      c(sum(gtoaprop[rnainds] >= prop), length(rnainds))))



numStops <- apply(aa, 1, function(x){length(which(x == "*"))})
fisher.test(cbind( c(sum(numStops[dnainds] > 0), length(dnainds)),
      c(sum(numStops[rnainds] > 0), length(rnainds))))



possibleGtoAs <- length(which(cons == "G"))
rnainds <- grep("PLASMA|RNA", inf$samp.loc)
dnainds <- grep("DNA", inf$samp.loc)
gtoaprop <- apply(nuc[,which(cons == "G")], 1, function(x){length(which(x == "A"))})/possibleGtoAs
rvd <- rep("RNA", length(gtoaprop))
rvd[dnainds] <- "DNA"

plotGtoAs <- tbl_df(cbind(gtoaprop, rvd, inf$samp.loc)) %>% mutate(gtoaprop <- as.numeric(gtoaprop)) 
p1 <- plotGtoAs %>%
    ggplot() + geom_histogram(mapping = aes(x = gtoaprop, fill = rvd))

plotGtoAs %>% ggplot() + stat_ecdf(mapping = aes(x = gtoaprop, col = rvd))

layout(matrix(1:2, nrow = 1))
hist(gtoaprop, xlab = "Proportion of G->As (observed out of possible)", main = "")
hist(gtoaprop[rnainds],  add = TRUE, col = rgb(0, 0, 1, .2))
hist(gtoaprop[dnainds], add = TRUE, col = rgb(1, 0, 0, .2))
legend("topright", c("All", "RNA", "DNA"),fill = c("white", rgb(0, 0, 1, .2), rgb(1, 0, 0, .2)), pch = 20)
cdf <- quantile(gtoaprop, seq(0, 1, by = .005))
cdf.r <- quantile(gtoaprop[rnainds], seq(0, 1, by = .005))
cdf.d <- quantile(gtoaprop[dnainds], seq(0, 1, by = .005))
xvals <- as.numeric(gsub("%", "", (names(cdf))))
plot(xvals, cdf, type = "l", xlab = "quantile")
lines(xvals, cdf.r, type = "l", xlab = "quantile", col = rgb(0, 0, 1, .7))
lines(xvals, cdf.d, type = "l", xlab = "quantile", col = rgb(1, 0, 0, .7))
legend("topleft", c("All", "RNA", "DNA"),col = c("black", rgb(0, 0, 1, .7), rgb(1, 0, 0, .7)), lty = "solid")


hasStop <- (apply(aa, 1, function(x){length(which(x == "*"))}))
table(hasStop, rvd)
table(hasStop, inf$samp.loc)
sum(aa[which(hasStop == 1),71] == "*")

aa[which(hasStop == 1),229]

aa[which(hasStop == 1),212]


hasStop <- apply(aa, 1, function(x){which(x == "*")})
table(unlist(hasStop))

nuc[1:10, 214:216]

table(inf$samp.loc[which(hasStop >= 1)], inf$week[which(hasStop >= 1)])


table(monkid[which(aa[,71] == "*")])

table(apply(nuc[, 214:216], 1, paste, collapse = ""))



#Function to determine what percentage of the RNA is drug resistant at a certain time
percentageDRbyComp <- function(monk, weekvalind, comp, mutIdent, mutPos){
    #First get the relevant indices. They should be all inds at which
    relInds <- intersect(
        #The monkey is correct
        which(monkid == monk),
        intersect(
            #The sample location is either plasma or vRNA
            grep(comp, samp.loc),
            #The week ind is correct
            which(weekind == weekvalind )
        )
    )
    #Return the proportion of the compartment that has the correct aa (relative to all samples at that week)
    return(length(grep(mutIdent, aa[relInds,mutPos]))/length(relInds))
}




percentDR <- foreach(monk = monknames, .combine= "rbind") %do% {
    pos <- 184
    type <- "V|I"
    if(monk == "A99615"){
        pos <- 103
        type <- "N"
    }
    foreach(weekindex = 1:5, .combine= "rbind") %do% {
        var.LN <- percentageDRbyComp(monk, weekindex, "LNRNA", type, pos)
        foreach(comps = c("PLASMA", "GUTRNA", "PBMCRNA", "VAGRNA"),
             .combine= "rbind") %do% {
            var.other <- percentageDRbyComp(monk, weekindex, comps, type, pos)
            c(monk, comps, weekindex, var.other, var.LN)
        }
    }
}


percDR <- tbl_df(percentDR)
colnames(percDR) <- c("monk", "comp", "weekind", "drperc", "lnperc")
percDR <- percDR %>% mutate(drperc = as.numeric(drperc), lnperc = as.numeric(lnperc), weekind = as.numeric(weekind))


DRtest <- percDR %>% filter(!is.nan(drperc),  !(drperc == 1 & lnperc == 1), !(drperc == 0 & lnperc == 0) ) %>%
    mutate(lnless = as.numeric(lnperc < drperc),
           lneq = as.numeric(lnperc == drperc),
           lnmore = as.numeric(lnperc > drperc)) %>%
        group_by(comp)%>%
            summarize(less = sum(lnless,  na.rm = TRUE),
                      eq = sum(lneq,  na.rm = TRUE), 
                      more = sum(lnmore,  na.rm = TRUE))


DRtest


percentDR.dna <- foreach(monk = monknames, .combine= "rbind") %do% {
    pos <- 184
    type <- "V|I"
    if(monk == "A99615"){
        pos <- 103
        type <- "N"
    }
    foreach(weekindex = 1:5, .combine= "rbind") %do% {
        var.gut <- percentageDRbyComp(monk, weekindex, "GUTDNA", type, pos)
        foreach(comps = c("PBMCDNA", "VAGDNA", "LNDNA"),
             .combine= "rbind") %do% {
            var.other <- percentageDRbyComp(monk, weekindex, comps, type, pos)
            c(monk, comps, weekindex, var.other, var.LN)
        }
    }
}




percDRrVd <- tbl_df(percentDR.dna)
colnames(percDRrVd) <- c("monk", "comp", "weekind", "dnaperc", "rnaperc")
percDRrVd <- percDRrVd %>% mutate(rnaperc = as.numeric(rnaperc), dnaperc = as.numeric(dnaperc), weekind = as.numeric(weekind))

percDRrVd  %>% filter(!is.nan(dnaperc),  !(dnaperc == 1 & rnaperc == 1), !(dnaperc == 0 & rnaperc == 0) ) %>%
    mutate(lnless = as.numeric(dnaperc < rnaperc),
           lneq = as.numeric(dnaperc == rnaperc),
           lnmore = as.numeric(dnaperc > rnaperc)) %>%
        group_by(comp)%>%
            summarize(less = sum(lnless,  na.rm = TRUE),
                      eq = sum(lneq,  na.rm = TRUE), 
                      more = sum(lnmore,  na.rm = TRUE))



binom.test(0, 6)
