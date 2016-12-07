
#For a given monkey, weekval (i.e., 1 = weeks 12/13, 2 = weeks 15/16, etc), mut identity and position, determine the percentage of the vRNA (including plasma) that is drug resistant

#This was to calculate the week at which most of the RNA was resistant, but I ended up hacking it a bit (see maxInds)
maxDRs <- list()
for(monk in monknames){
    toAdd <- c()
    for(i in sort(unique(weekind[monkid == monk]))){
        if(monk == "T98133" | monk == "A99039"){
            toAdd[i] <- percentageDR.wkind(monk, i, "V|I", c(184))
        }else if(monk == "A99165" | monk == "A01198"){
            toAdd[i] <- percentageDR.wkind(monk, i, "N", c(103))
        }
    }
    maxDRs[[monk]] <- toAdd
}


#maxInds returns the first point at which the RNA is more than 90% DR OR the maximally DR if the RNA is never 90+% DR.
maxInds <- unlist(lapply(maxDRs, function(x){
    #if the vRNA is never observed to be at frequency over 90%,
    if(sum(x > .9) == 0){
        # just return the ind at which it's maximized
        return(which.max(x))
    }else{
        #otherwise, return the point at which it's first greater than .9
        return(min(which(x > .9)))
    }}
))



#For the ordered monkeys (T98, 165, 039, 198), what are the relevant DRMs
mutPos <- c(184, 103, 184, 103)
mutIdent <- c("I|V", "N", "I|V", "N")

#thetas will store the values of theta that will eventually be returned
thetas <- matrix(data = NA, nrow = 4, ncol = 4)
colnames(thetas) <- c("Theta_nuc (DR vRNA)", "Theta_hap (DR vRNA)", "Theta_nuc (DR all)", "Theta_hap (DR all)")

#for the vRNA from each macaque
for(monk in monknames){
    #get all the indices for that monk
    monkInd <- which(monk == monknames)
    Nes <- c()
    #For both RNA and RNA+DNA
    for(rnaOrAll in c("RNA|PLASMA", "RNA|DNA|PLASMA")){
        monkRNAinds <- intersect(which(monkid == monk), grep(rnaOrAll, samp.loc))
        relTime <- maxInds[monkInd]
        #Determine ALL the relevant indices for the timepoint we found above
        relIndsAll <- intersect(monkRNAinds, which(weekind == relTime))

        #But now, we only want the drug resistant variants.
        #relIndsDR are our drug resistant, rna or rna/dna inds from a monkey at the timepoint
        relIndsDR <- relIndsAll[grep(mutIdent[monkInd], aa[relIndsAll,mutPos[monkInd]])]

#Now, onto computing H
#Method 1: different identities (i.e., I versus V) count as different variants

        #store the relevant position (i.e., 184 or 103) in pos
        pos <- mutPos[which(monk == monknames)]
        #Look at all the way the nucleotides encode that position 
        AAencodes <- apply(nuc[relIndsDR, (3*(pos+1) -2):(3*(pos+1))], 1, paste, collapse = "")
        #Compute heterozygosity by determining the count of each variant (table) divided by the length to get frequencies, which are then squared and summed. We need to subtract from 1 to go from homozygosity to heterozygosity
        H <- 1-sum((table(AAencodes)/length(relIndsDR))^2)
        theta.est.vars <- H/(1-H)

        #Method 2: different linked variants
        #This is pretty ugly, but I'm determining which positions (135-875, because the others are sort of messy) have a minor allele at frequency greater than 10% of the sample size
        polymorphicInds <- (135:875)[which(unlist(lapply(apply(nuc[relIndsDR, 135:875], 2, table), function(x){
            #sort all the different alleles at a position among the relevant inds and store the second most frequent one
            tmp <- sort(x, decreasing = TRUE)[2]
            #If that second most frequent allele is at frequency higher than 10% (and also exists), keep track of that index
                (!is.na(tmp) & tmp > length(relIndsDR)/10)
#            (!is.na(tmp) & tmp >= 2)
        } )) == TRUE)]

        #Compute H among THESE different frequencies, as before
        H <- 1-sum((table(apply(nuc[relIndsDR,polymorphicInds], 1, paste, collapse = ""))/length(relIndsDR))^2)
        theta.est.haps <- H/(1-H)
        #Store
        Nes <- c(Nes, c(theta.est.vars, theta.est.haps))
    }
    thetas[which(monk == monknames),] <- Nes
}



#Andy Leigh Brown approach
ALB.ests <- c()
for(monk in monknames){
    monkRNAinds <- intersect(which(monkid == monk), grep("PLASMA", samp.loc))
    relTime <- maxInds[which(monk == monknames)]
    relIndsAll <- intersect(monkRNAinds, which(weekind == relTime))
    ALB.ests[which(monk == monknames)] <- mean(distMat[relIndsAll,relIndsAll][lower.tri(distMat[relIndsAll,relIndsAll])])/length(135:900)
}


fulltab <- cbind(ALB.ests, thetas)
rownames(fulltab) <- monknames
colnames(fulltab)[1] <- "Theta_d (Plasma)"

#Convert from thetas to Nes 
mu.transition <- 1.4*10^(-5)
mu.transversion <- 2*10^(-6)
fulltab[1,] <- fulltab[1,]/(2*mu.transition)
fulltab[2,] <- fulltab[2,]/(2*mu.transversion)
fulltab[3,] <- fulltab[3,]/(2*mu.transversion)
fulltab[4,] <- fulltab[4,]/(2*mu.transition)

#print the results
print.xtable(xtable(fulltab, digits = -2),  type = "html", paste("../out/tables/theta.html", sep = ""))



 
