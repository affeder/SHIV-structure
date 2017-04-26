#format and plot all data

ordloc <- c("PLASMA", "PBMC", "LN", "GUT", "VAG")
ordlocrna <- c("PLASMA", "PBMCRNA", "LNRNA", "GUTRNA", "VAGRNA")
namedordlocrna <- c("Plasma", "PBMC vRNA", "LN vRNA", "Gut vRNA", "Vagina vRNA")


#Read in SM data
#smres_in <- read.table("../out/sm/SM.allresults-ds-t2.txt", header = FALSE)
smres_in <- read.table("../out/sm/SM.out.txt", header = FALSE)
smres <- tbl_df(smres_in)
colnames(smres) <- c("monk", "week", "loc1", "loc2","ds", "trial", "p")

smres <- smres %>% mutate(bothComps = paste(loc1,loc2, sep = "|")) 

ordloc1 <- rep(NA, nrow(smres))
ordloc2 <- rep(NA, nrow(smres))
for(i in 1:length(ordloc)){
    doubleMatch <- which(smres$bothComps == paste(ordloc[i], ordloc[i], sep = "|"))
    ordloc1[doubleMatch] <- ordloc[i]
    ordloc2[doubleMatch] <- ordloc[i]
    matchesComp <- grep(ordloc[i], smres$bothComps)
    ordloc1[intersect(matchesComp, which(is.na(ordloc1)))] <- ordloc[i]
    ordloc2[intersect(intersect(matchesComp, which(is.na(ordloc2))),
                      which(ordloc1 != ordloc[i]))] <- ordloc[i]
}

#for SM, I also want to replace "VAG" with  "VAGRNA", "LN" with "LNRNA", etc. 
#It's surprisingly MUCH easier to do this in base than in a dplyr framework
ordloc1[ordloc1 == "VAG"] <- "VAGRNA"
ordloc2[ordloc2 == "VAG"] <- "VAGRNA"
ordloc1[ordloc1 == "LN"] <- "LNRNA"
ordloc2[ordloc2 == "LN"] <- "LNRNA"
ordloc1[ordloc1 == "PBMC"] <- "PBMCRNA"
ordloc2[ordloc2 == "PBMC"] <- "PBMCRNA"
ordloc1[ordloc1 == "GUT"] <- "GUTRNA"
ordloc2[ordloc2 == "GUT"] <- "GUTRNA"
ol2 <- ordloc2
ol1 <- ordloc1

#Assign the new ordered functions
smres <- smres %>%
    mutate(ordloc1 = ol1, ordloc2 = ol2)


smres$p <- as.numeric(smres$p)
smres$ds <- as.numeric(smres$ds)

smres$ordloc1 <- factor(smres$ordloc1, levels = ordlocrna)
smres$ordloc2 <- factor(smres$ordloc2, levels = ordlocrna)


#KST results

#Read in KST data
kstout <- tbl_df( read.table("../out/kst/kstout.txt", header = TRUE) )

#This orders the factors correctly
kstout <- kstout %>% mutate(bothComps = paste(loc1,loc2, sep = "|"))
ordloc1 <- rep(NA, nrow(kstout))
ordloc2 <- rep(NA, nrow(kstout))
for(i in 1:length(ordlocrna)){
    doubleMatch <- which(kstout$bothComps == paste(ordlocrna[i], ordlocrna[i], sep = "|"))
    ordloc1[doubleMatch] <- ordlocrna[i]
    ordloc2[doubleMatch] <- ordlocrna[i]
    matchesComp <- grep(ordloc[i], kstout$bothComps)
    ordloc1[intersect(matchesComp, which(is.na(ordloc1)))] <- ordlocrna[i]
    ordloc2[intersect(intersect(matchesComp, which(is.na(ordloc2))),
                      which(ordloc1 != ordlocrna[i]))] <- ordlocrna[i]
}
ol1 <- ordloc1
ol2 <- ordloc2
kstout <- kstout %>%
    mutate(ordloc1 = ol1, ordloc2 = ol2)
kstout$ordloc1 <- factor(kstout$ordloc1, levels = ordlocrna)
kstout$ordloc2 <- factor(kstout$ordloc2, levels = ordlocrna)



#AMOVA results

#Read in AMOVA data
amovaout <- tbl_df( read.table("../out/amova/amovaout.txt", header = TRUE) )

#This orders the factors correctly
amovaout <- amovaout %>% mutate(bothComps = paste(loc1,loc2, sep = "|"))
ordloc1 <- rep(NA, nrow(amovaout))
ordloc2 <- rep(NA, nrow(amovaout))
for(i in 1:length(ordlocrna)){
    doubleMatch <- which(amovaout$bothComps == paste(ordlocrna[i], ordlocrna[i], sep = "|"))
    ordloc1[doubleMatch] <- ordlocrna[i]
    ordloc2[doubleMatch] <- ordlocrna[i]
    matchesComp <- grep(ordloc[i], amovaout$bothComps)
    ordloc1[intersect(matchesComp, which(is.na(ordloc1)))] <- ordlocrna[i]
    ordloc2[intersect(intersect(matchesComp, which(is.na(ordloc2))),
                      which(ordloc1 != ordlocrna[i]))] <- ordlocrna[i]
}
ol1 <- ordloc1
ol2 <- ordloc2
amovaout <- amovaout %>%
    mutate(ordloc1 = ol1, ordloc2 = ol2)
amovaout$ordloc1 <- factor(amovaout$ordloc1, levels = ordlocrna)
amovaout$ordloc2 <- factor(amovaout$ordloc2, levels = ordlocrna)



dat <- smres %>%
    select(monk, week, ordloc1, ordloc2, ds, trial, p) %>%
    mutate(stattype = "SM")
dat <- amovaout %>% 
    select(monk, week, ordloc1, ordloc2, ds, trial, p) %>%
    mutate(stattype = "AMOVA") %>%
    do(rbind(dat, .))
dat <- kstout %>%
    select(monk, week, ordloc1, ordloc2, ds, trial, p) %>%
    mutate(stattype = "KST") %>%
    do(rbind(dat, .))

weekval <- rep(NA, length = nrow(dat))
weekval[dat$week == 1] <- 13
weekval[dat$week == 2] <- 15
weekval[dat$week == 3] <- 20
weekval[dat$week == 4] <- 26
weekval[dat$week == 5 & dat$monk == "A99039"] <- 38
weekval[dat$week == 5 & dat$monk == "A99165"] <- 29
weekval[dat$week == 5 & dat$monk == "A01198"] <- 44

dat <- mutate(dat, week.t = weekval)

#All time points, all macaques
doubleddat <- dat %>%
    mutate(tmp = ordloc1, ordloc1 = ordloc2, ordloc2 = tmp) %>%
        select(-tmp) %>%
            do(rbind(dat,. ))



