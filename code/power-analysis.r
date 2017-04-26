#Power analysis 

#Read in and format data
kstpower <- tbl_df(read.table("../out/kst/kst.power.txt", header = TRUE))
amovapower <- tbl_df(read.table("../out/amova/amovaout-power.txt", header = TRUE))
smpower <- tbl_df(read.table("../out/sm/SM.allresults-ds-t2.txt", header = FALSE))
colnames(smpower) <- c("monk", "week", "loc1", "loc2","ds", "trial", "p")

ordlocrna <- c("PLASMA", "PBMCRNA", "LNRNA", "GUTRNA", "VAGRNA")

#First, we convert the SM names (abbreviated) into the names used in the
# other two tests
ordloc1 <- as.character(smpower$loc1)
ordloc2 <- as.character(smpower$loc2)
ordloc1[ordloc1 == "VAG"] <- "VAGRNA"
ordloc2[ordloc2 == "VAG"] <- "VAGRNA"
ordloc1[ordloc1 == "LN"] <- "LNRNA"
ordloc2[ordloc2 == "LN"] <- "LNRNA"
ordloc1[ordloc1 == "PBMC"] <- "PBMCRNA"
ordloc2[ordloc2 == "PBMC"] <- "PBMCRNA"
ordloc1[ordloc1 == "GUT"] <- "GUTRNA"
ordloc2[ordloc2 == "GUT"] <- "GUTRNA"
smpower <- smpower %>% mutate(loc1 = ordloc1, loc2 = ordloc2)

kToMerge <- kstpower %>% select(-starts_with("size")) %>% mutate(type = "KST")
aToMerge <- amovapower %>% select(-starts_with("size")) %>% mutate(type = "AMOVA")
sToMerge <- smpower %>% mutate(type = "SM")

fullPower <- kToMerge %>% do(rbind(., aToMerge)) %>% do(rbind(., sToMerge))

fullPower <- fullPower %>% mutate(bothComps = paste(loc1,loc2, sep = "|"))

ol1 <- rep(NA, nrow(fullPower))
ol2 <- rep(NA, nrow(fullPower))
for(i in 1:length(ordlocrna)){
    doubleMatch <- which(fullPower$bothComps == paste(ordlocrna[i], ordlocrna[i], sep = "|"))
    ol1[doubleMatch] <- ordlocrna[i]
    ol2[doubleMatch] <- ordlocrna[i]
    matchesComp <- grep(ordloc[i], fullPower$bothComps)
    ol1[intersect(matchesComp, which(is.na(ol1)))] <- ordlocrna[i]
    ol2[intersect(intersect(matchesComp, which(is.na(ol2))),
                      which(ol1 != ordlocrna[i]))] <- ordlocrna[i]
}
fullPower <- fullPower %>%
    mutate(ordloc1 = ol1, ordloc2 = ol2)
fullPower$ordloc1 <- factor(fullPower$ordloc1, levels = ordlocrna)
fullPower$ordloc2 <- factor(fullPower$ordloc2, levels = ordlocrna)


fullPower <- fullPower %>% select(-starts_with('loc')) %>%
    select(-matches('bothComps'))


listOfComps <- list(
    'PLASMA'="Plasma",
    'PBMCRNA'="PBMC",
    'GUTRNA'="Gut",
    'LNRNA'="LN",  
    'VAGRNA'="Vagina"
)

comp_labeller <- function(variable,value){
  return(listOfComps[value])
}

for(testtype in c("KST", "AMOVA", "SM")){
    pdf(paste("../out/graphs/S8-",testtype,".pdf", sep = ""), width = 7.5, height = 5)
    g <- fullPower %>% filter(!is.na(p)) %>%
        filter(type == testtype) %>%
        group_by(monk, week, ds, ordloc1, ordloc2) %>%
        summarize(prob = mean(p <= .05)) %>%
        ggplot(mapping = aes(x = ds, y = prob, col = monk, pch = as.factor(week), lty = as.factor(week))) +
            geom_point() + geom_line() + 
            facet_grid(ordloc1 ~ ordloc2, labeller = comp_labeller) +
            labs(y = "Probability of a significant test", x = "Downsampled compartment size", lty = "Sampling\nTime point", col = "Macaque", pch = "Sampling\nTime point") +
                theme_bw()
    print(g)
    dev.off()  
}

