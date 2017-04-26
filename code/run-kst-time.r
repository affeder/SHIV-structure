set.seed(9418247)
source("kst-funs.r")

rnacomps <- c("VAGRNA", "PBMCRNA", "PLASMA", "LNRNA", "GUTRNA")

#For paper
#numreps <- 1000
#nullNum <- 1000

#For testing
numreps <- 10
nullNum <- 10

#Consecutive timepoints
allcomps <- unique(inf$samp.loc)
#allcomps <- allcomps[-length(allcomps)]
K <- 2

kstresults.time <- foreach(monk = monknames, .combine = "rbind", .errorhandling="remove") %do% {
    foreach(comp = allcomps, .combine = "rbind", .errorhandling = "remove") %do% {
        foreach(downSampNum = c(10) , .combine = "rbind", .errorhandling="remove") %do% {
            foreach(t = 1:4, .combine = "rbind", .errorhandling="remove") %do% {
                foreach(rep = 1:numreps, .combine = "rbind", .errorhandling="remove") %do% {

                    pvalToReport <- NA
                    keepChecking <- TRUE
                    nextTime <- t+1
                    
                    samp1 <- intersect(which(monkid == monk & samp.loc == comp),
                                       which(weekind == t))
                    samp2 <- intersect(which(monkid == monk & samp.loc == comp),
                                       which(weekind == nextTime))

                    while(keepChecking){
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
                            keepChecking <- FALSE
                        }else{
                            nextTime <- nextTime + 1
                            samp2 <- intersect(which(monkid == monk & samp.loc == comp),
                                               which(weekind == nextTime ))
                        }
                        if(nextTime > 5){
                            keepChecking <- FALSE
                        }
                    }
                    if(!is.na(pvalToReport)){
                        c(monk, t, nextTime, comp, downSampNum, rep, pvalToReport)
                    }
                }
            }
        }
    }
}



rownames(kstresults.time) <- NULL
colnames(kstresults.time) <- c("monk", "week1", "week2", "loc",  "ds", "trial", "p")
kstout.time <- tbl_df(kstresults.time)

kstout.time$ds <- as.numeric(kstout.time$ds)
kstout.time$p <- as.numeric(kstout.time$p)
kstout.time <- kstout.time %>% mutate(p = 1-p)

kstout.time$week1 <- as.numeric(kstout.time$week1)
kstout.time$week2 <- as.numeric(kstout.time$week2)

write.table(kstout.time, "../out/kst/kstout.time.txt", row.names = FALSE, col.names = TRUE)

