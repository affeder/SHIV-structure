set.seed(2741926)
source("kst-funs.r")

rnacomps <- c("VAGRNA", "PBMCRNA", "PLASMA", "LNRNA", "GUTRNA")

#For paper
#numreps <- 100
#nullNum <- 1000

#For testing
numreps <- 10
nullNum <- 10

K <- 2
AllCombs <- combn(rnacomps,K)
AllCombs <- cbind(AllCombs,rbind(rnacomps, rnacomps) )
kstresults.power <- foreach(monk = monknames, .combine = "rbind", .errorhandling="remove") %do% {
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
                        statval <- kstList(list(downsamp1, downsamp2))
                        randvals <- c()
                        for(j in 1:nullNum ){
                                        #shuffle the assignments in the list keeping track of
                                        # which indices belong to which compartment
                            listOfInds <- randIt(list(downsamp1, downsamp2))
                                        #Recompute the statistic
                            randvals[j] <- kstList(listOfInds)
                        }
                        pvalToReport <- sum(statval > randvals)/nullNum
                    }
                    c(monk, t, t(AllCombs[,coldind]), length(samp1), length(samp2), downSampNum, rep, pvalToReport)
                }
            }
        }
    }
}


rownames(kstresults.power) <- NULL
colnames(kstresults.power) <- c("monk", "week", "loc1", "loc2", "size1", "size2", "ds", "trial", "p")
kst.power <- tbl_df(kstresults.power)

kst.power$ds <-  as.numeric(kst.power$ds)
kst.power$p <-  as.numeric(kst.power$p)
kst.power <- mutate(kst.power, p = 1-p)

write.table(kst.power, "../out/kst/kst.power.txt", row.names = FALSE, col.names = TRUE)

