#We're going to be generating random numbers, so we'll set a random seed
set.seed(4019284)

rnacomps <- c("VAGRNA", "PBMCRNA", "PLASMA", "LNRNA", "GUTRNA")

#In paper:
#numreps <- 1000
#nullNum <- 1000

#For testing
numreps <- 10
nullNum <- 100

K <- 2
AllCombs <- combn(rnacomps,K)
AllCombs <- cbind(AllCombs,rbind(rnacomps, rnacomps) )
amovaresults <- foreach(monk = monknames, .combine = "rbind", .errorhandling="remove") %do% {
    foreach(t = 1:5, .combine = "rbind", .errorhandling="remove") %do% {
        foreach(coldind = 1:ncol(AllCombs), .combine = "rbind", .errorhandling="remove") %do% {
            downSampNum = 10
            foreach(rep = 1:numreps, .combine = "rbind", .errorhandling="remove") %do% {
                samp1 <- intersect(
                    which(monkid == monk & samp.loc == AllCombs[1,coldind]),
                    which(weekind == t))
                samp2 <- intersect(
                    which(monkid == monk & samp.loc == AllCombs[2,coldind]),
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
                        distIndsToSto[i] <- union(
                            downsamp1[which(haps1 == abundances[i,1])],
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


rownames(amovaresults) <- NULL
colnames(amovaresults) <- c("monk", "week", "loc1", "loc2", "size1", "size2", "ds", "trial", "p")
amovaout <- tbl_df(amovaresults)

amovaout$p <- as.numeric(amovaout$p)
amovaout$ds <- as.numeric(amovaout$ds)

write.table(amovaout, "../out/amova/amovaout.txt", row.names = FALSE, col.names = TRUE)






