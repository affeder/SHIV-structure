#We're going to be generating random numbers, so we'll set a random seed
set.seed(4019284)

rnacomps <- c("VAGRNA", "PBMCRNA", "PLASMA", "LNRNA", "GUTRNA")

#In paper
#numreps <- 100
#nullNum <- 1000

#For testing
numreps <- 10
nullNum <- 100

K <- 2
AllCombs <- combn(rnacomps,K)
AllCombs <- cbind(AllCombs,rbind(rnacomps, rnacomps) )
amovaresults.power <- foreach(monk = monknames, .combine = "rbind", .errorhandling="remove") %do% {
    foreach(t = 1:5, .combine = "rbind", .errorhandling="remove") %do% {
        foreach(coldind = 1:ncol(AllCombs), .combine = "rbind", .errorhandling="remove") %do% {
            foreach(downSampNum = c(3, 7, 10, 15, 20) , .combine = "rbind", .errorhandling="remove") %do% {
                foreach(rep = 1:numreps, .combine = "rbind", .errorhandling="remove") %do% {
                    samp1 <- intersect(which(monkid == monk & samp.loc == AllCombs[1,coldind]),
                                       which(weekind == t))
                    samp2 <- intersect(which(monkid == monk & samp.loc == AllCombs[2,coldind]),
                                       which(weekind == t))
                    pvalToReport <- NA
                    if(length(samp1) >= downSampNum & length(samp2) >= downSampNum){
                        downsamp1 <- sample(samp1, downSampNum, replace = FALSE)
                        downsamp2 <- sample(samp2, downSampNum, replace = FALSE)
                        haps1 <- apply(nuc[downsamp1, 135:900], 1, paste, collapse = "")
                        haps2 <- apply(nuc[downsamp2, 135:900], 1, paste, collapse = "")
                        hap1Table <- table(haps1)
                        haps1Freqs <- as.data.frame(hap1Table)
                        rownames(haps1Freqs) <- names(hap1Table)
                        hap2Table <- table(haps2)
                        haps2Freqs <- as.data.frame(cbind(hap = names(hap2Table),
                                                      abund = as.numeric(hap2Table)))
                        hap1Table <- table(haps1)
                        haps1Freqs <- as.data.frame(cbind(hap = names(hap1Table),
                                                          abund = as.numeric(hap1Table)))
                        abundances <- merge(haps1Freqs, haps2Freqs, by = "hap", all = TRUE)
                        abundances$abund.y <- as.numeric.factor(abundances$abund.y)
                        abundances$abund.x <- as.numeric.factor(abundances$abund.x)
                        abundances[is.na(abundances)] <- 0
                        distIndsToSto <- c()
                        for(i in 1:nrow(abundances)){
                            distIndsToSto[i] <- union(downsamp1[which(haps1 == abundances[i,1])],
                                                  downsamp2[which(haps2 == abundances[i,1])])[1]
                        }
                        distsForFunc <- as.dist(distMat[distIndsToSto, distIndsToSto])
                        abundsForFunc <- abundances[,2:3]
                        amovaOut <- amova(abundsForFunc, quasieuclid(distsForFunc))
                        randTestOut <- randtest(amovaOut, nrepet = nullNum)
                        pvalToReport <- randTestOut$p
                    }
                    c(monk, t, t(AllCombs[,coldind]), length(samp1), length(samp2), downSampNum, rep, pvalToReport)
                }
            }
        }
    }
}



rownames(amovaresults.power) <- NULL
colnames(amovaresults.power) <- c("monk", "week", "loc1", "loc2", "size1", "size2", "ds", "trial", "p")
amova.power <- tbl_df(amovaresults.power)

amova.power$p <- as.numeric(amova.power$p)
amova.power$ds <- as.numeric(amova.power$ds)

amova.power %>% filter(!is.na(p)) %>%
    group_by(monk, week, loc1, loc2, ds) %>%
        summarize(prob = mean(p <= .05)) %>%
            ggplot() + geom_point(mapping = aes(x = ds, y = prob, col = week)) +
                facet_grid(loc1 ~ loc2)
    
write.table(amova.power, "../out/amova/amovaout-power.txt", row.names = FALSE, col.names = TRUE)





