

#Megaplot.gen is a function that makes the multi-paneled figures for each macaque
#It works in a few stages, using primarily the helper function setup.hierchy
#setup.hierchy determines and then establishes an order of mutations (i.e., a sequence with muatation "A" should be a parent to sequences with mutations "A" and "B" (more examples are worked below).
#setup hierchy is called on all of the data first - the establishes a consistent color scheme which is used throughout
#Then, setup hierchy is called on each compartment in question. 
megaplot.gen <- function(monk, vectorofcomps, ylabs, maintitles){

    #First, choose ALL variants that are going to appear on the plot
    tmpnuc <- nuc[which(monkid == monk),]
    #Any location that has a minor allele at frequency > 1% should be included
    maf <- nrow(tmpnuc)/100
    hapinds <- c()
    for(j in 1:ncol(tmpnuc)){
        tmp <- table(tmpnuc[,j])
        dashind <- which(names(tmp) == "-")
        if(length(dashind) > 0){
            tmp <-tmp[-dashind]
        }
        tmp <- sort(tmp, decreasing = TRUE)
        if(length(tmp) > 1){
            if(tmp[2] > maf){
                hapinds[length(hapinds)+1] <- j
            }
        }
    }

    #hapinds now has all of the variable haplotypes
    hapinds <<- hapinds
    haps <<- nuc[,hapinds]
    #haps are all of the haplotypes that are possible
    #Note, each monkey will use its own indices (that are relevent)
    
    indsToProcess <- which(monkid == monk)
    #Call the main function to determine the hierchy
    hierchinf <- setup.hierchy(indsToProcess)

    #Store the returned information 
    dat <- hierchinf$dat
    haps.over.time <- hierchinf$haps.over.time
    mut.relat <- hierchinf$mut.relat
    all.comp.inds <- hierchinf$all.comp.inds

    #Compute names (i.e., the list of mutations) for each of the sequences
    allNames <- c()
    for(i in 1:nrow(mut.relat)){
        if(sum(which(mut.relat[i,] == 1)) == 0){
            allNames[i] <- "WT"
        }else{
            allNames[i] <- labelmaker(names(which(mut.relat[i,] == 1)), all.comp.inds)
        }
    }

    #Assign colors for each sequence
    cols <- rep("white", nrow(haps.over.time))
    #These are our base colors palettes 
    colset <-  c(brewer.pal(8, "Set2")[1:6], brewer.pal(12, "Paired")[c(2, 6, 9, 10, 11)])
    #How many haplotypes are going to have "base colors" versus variations?
    maxuse <- 10
    
    #First, set some sort of base colors for the most frequent haplotypes
    popnames <- allNames[order(apply(haps.over.time, 1, sum), decreasing = TRUE)][1:maxuse]
    popinds <- c()
    for(i in popnames){
        cols[ which(allNames == i) ] <- colset[which(i == popnames)]
        popinds <- c(popinds, which(allNames == i))
    }

    #Now, go through the descendents and choose colors for each descendent
    for(i in 1:nrow(haps.over.time)){
        if(sum(i == popinds) == 0){
            inq <- matrix(rep(mut.relat[i,], length(popinds)),
                          nrow = length(popinds), byrow = TRUE)
            distvals <- (apply(mut.relat[popinds, ] - inq, 1, function(x){sum(abs(x))}))
            bestAnc <- which(distvals == min(distvals))
            if(length(bestAnc) > 1){
                mutlens <- apply(mut.relat[popinds,], 1, sum)
                bestAnc <- bestAnc[which(mutlens[bestAnc] == max(mutlens[bestAnc]))[1]]
            }
            cols[i] <- paste("anc", bestAnc, sep = "-")
        }
    }

    for(i in 1:maxuse){
        colsReassign <- which(cols == paste("anc", i, sep = "-"))
        ind <- 1
        compCols <- setdiff(1:maxuse, i)
        numPerComp <- ceiling(length(colsReassign)/7)

        for(j in compCols){
            colmake <- colorRampPalette(c(colset[i], colset[j]))
            tmpCols <- colmake( 2*(numPerComp + 1)+ 1)
            for(k in 1:numPerComp){
                cols[colsReassign[ind]] <- tmpCols[1 + k]
                ind <- ind + 1
            }
        }
    }

    allNames
    colMap <- c("")
    nameMap <- c("")
    initshape <- dat$mypoly
    kidsToPlot <- dat$children
    access.kidsToPlot <- 1:length(kidsToPlot)
    currep <- 1
    while(length(access.kidsToPlot) > 0){
        currInd <- kidsToPlot[[ access.kidsToPlot[1] ]]$index
        mutname <- labelmaker(names(which(mut.relat[currInd,] == 1)), all.comp.inds)
        nameMap[currep] <- mutname
        colMap[currep] <- cols[which(allNames == mutname)]
        kidsList <- kidsToPlot[[access.kidsToPlot[1]]]$children
        access.kidsToPlot <- access.kidsToPlot[-1]
        if(length(kidsList) > 0){
            for(i in 1:length(kidsList)){
                kidsToPlot[[length(kidsToPlot)+1]] <- kidsList[[i]]
                access.kidsToPlot <- c(length(kidsToPlot), access.kidsToPlot)
            }
        }
        currep <- currep + 1
    }

    modifier <- 1.15
    didntprint <- c()
    png(paste("../out/graphs/",monk,".png", sep = ""),
        height = 2250*modifier, width = 2625*.6, res = 300) #*length(ylabs))
    layout(matrix(1:length(vectorofcomps), ncol = 2, byrow = TRUE))
    par(oma = c(3,6, 2.65, 0))
    countind <- 1
    for(samplocname in vectorofcomps){
        if(length(retInd.id.comp(monk, samplocname)) > 15){
            plotComp.new(monk, samplocname, colMap, nameMap)
            if(countind %% 2 == 1){
                axis(2, at = seq(0, 30, by = 7.5),
                     labels = paste(seq(0, 1, by = .25)*100, "", sep = ""),
                     cex.axis = 1, las = 2)
                #ghee
                linenum <- 4
                linenum2 <- 2.5
                ylabval <- "% of Compartment"
                mtext(ylabval, side = 2, line = linenum2, cex = 1)
                mtext(ylabs[ceiling(countind/2)],  side = 2, line = linenum,  cex = 1.5)
            }
            if(countind == 1){
                mtext(maintitles[1], at = .25, side = 3, line = 1, outer = TRUE, cex = 1.5)
            }
            if(countind == 2){
                mtext(maintitles[2], at = .75, side = 3, line = 1, outer = TRUE, cex = 1.5)
                
            }
            if(countind == length(vectorofcomps)){
                    mtext("Week Post-Infection", side =1, line =3, cex = 1.5)
            }
            if(countind >= (length(vectorofcomps) -1)){
                xaxlabs.tmp <- as.numeric(colnames(haps.over.time))
                xaxlabs <- c()
                for(i in 1:length(xaxlabs.tmp)){
                    if(i == 1){ xaxlabs[i] <- xaxlabs.tmp[i] }
                    else{
                        if(xaxlabs.tmp[i] - xaxlabs.tmp[i -1] > 1){
                            xaxlabs <- c(xaxlabs, xaxlabs.tmp[i])
                        }
                    }
                }
                axis(1, at = xaxlabs, cex.axis = 1)
            }
        }else{
            didntprint <- c(didntprint, samplocname)
        }
        countind <- countind + 1
    }
    print(paste("Did not print:",didntprint))
    dev.off()
}

nucToAA <- function(pos){
    c(floor((pos-1)/3), ((pos-1)%%3) + 1)
}

#takes a nucleotide position and returns the three nuc positions making up the AA
Nucs <- function(pos){
    relaind <- nucToAA(pos)[2]
    if(relaind == 3){ toAdd <- (-2):0 }
    if(relaind == 2){ toAdd <- (-1):1 }
    if(relaind == 1){ toAdd <- 0:2 }
    return(pos + toAdd)
}

labelmaker <- function(listOfMuts, all.comp.inds){

    labList <- c()
    for(i in 1:length(listOfMuts)){

        mut.dec <- strsplit(listOfMuts[i], "-")[[1]]

        #This is the position in allmuts that was modified
        mut.dec.num <- as.numeric(mut.dec[1])

        #This is the position that was modified
        nucInd <- hapinds[all.comp.inds[mut.dec.num]]

        AAnum <- nucToAA(nucInd)[1]

        #This is what it was modified to
        newIdent <- mut.dec[2]

        #Again, we need a function that maps things into amino acids
        ancestral.nucs <- (cons[Nucs(nucInd)])
        derived.nucs <- ancestral.nucs
        derived.nucs[nucToAA(nucInd)[2]] <- newIdent

        ancestralAA <- translate(ancestral.nucs)
        derivedAA <- translate(derived.nucs)

        tmplab <- paste(ancestralAA, AAnum, derivedAA, sep = "")
        if(tmplab == "K103N"){
            tmplab <- paste("K103N-", toupper(newIdent), sep = "")
        }
        labList[i] <- tmplab
        
    }
    couldReturn <- paste(sort(labList), collapse = ",")
    return(couldReturn)
}

isDRM <- function(label){

    hasK103N <- sum(regexpr("K103N", label) > 0)
    hasK103S <- sum(regexpr("K103S", label) > 0)
    hasM184V <- sum(regexpr("M184V", label) > 0)
    hasM184I <- sum(regexpr("M184I", label) > 0)

    if((hasK103N + hasK103S + hasM184V + hasM184I) > 0){ return(1) }else{
        return(0) }
}

#This is the true workhorse function of these plots
setup.hierchy <- function(indsToProcess){

    #This function is going to take some indices, and report a few things
    #haps.over.time (How frequencies of different haplotypes change over time
    #whoismyparent
    #whoismychild
    #directchildren

   #Which weeks do we need to look at
    relweeks <- sort(as.numeric(names(table(inf[indsToProcess,'weeks']))))

    all.comp.inds <- sort(which(apply(haps[indsToProcess,], 2, function(x){c(length(table(x)) > 1)}) == TRUE))

    #if there is only one polymorphic ind?
    if(length(all.comp.inds) == 1){
        all.comp.haps <- unique(haps[indsToProcess,all.comp.inds])
    }else{
        all.comp.haps <-  names(table(apply(haps[indsToProcess, all.comp.inds], 1, paste, collapse = "")))
    }

    haps.over.time <- matrix(data = NA, nrow = length(all.comp.haps), ncol = length(relweeks))
    colnames(haps.over.time)<- relweeks
    rownames(haps.over.time) <- all.comp.haps

    for(i in 1:ncol(haps.over.time)){
        relIndsForWeek <- indsToProcess[inf[indsToProcess,]$weeks == as.numeric(colnames(haps.over.time)[i])]
        if(length(all.comp.inds) == 1){
            hapsThatWeek <- haps[relIndsForWeek, all.comp.inds]
        }else{
            hapsThatWeek <- apply(haps[relIndsForWeek, all.comp.inds], 1, paste, collapse = "")
        }
        for(j in 1:nrow(haps.over.time)){
            haps.over.time[j,i] <- sum(all.comp.haps[j] == hapsThatWeek)
        }
    }

    #Let's also have a matrix that keeps track of all of the expanded forms of the haplotypes (i.e., str split)
    all.comp.haps.split <- matrix(data = NA, nrow = length(all.comp.haps), ncol = nchar(all.comp.haps))
    for(i in 1:length(all.comp.haps)){
        all.comp.haps.split[i,] <- strsplit(all.comp.haps[i], split = "")[[1]]
    }

    
    if(length(all.comp.inds) == 1){
        reference <- names(sort(table(haps[,all.comp.inds]), decreasing = TRUE))[1]
    }else{
        reference <- paste(apply(haps[,all.comp.inds], 2, function(x){ c(names(sort(table(x), decreasing = TRUE))[1]) } ), collapse = "")
    }
    
    ref.split <- strsplit(reference, "")[[1]]

    #Ok, before we figure out how mutations were accumulated, let's figure out WHAT mutations were accumulated by what

    muts.per.hap <- apply(all.comp.haps.split , 1, function(x){ paste( which(x != ref.split),  x[which(x != ref.split)], sep = "-") })

    allmuts <- unique(unlist(muts.per.hap))
    mut.relat <- matrix(data = 0, ncol = length(allmuts), nrow = length(muts.per.hap))
    colnames(mut.relat) <- allmuts

    for(i in 1:length(muts.per.hap)){
        mut.relat[i, muts.per.hap[[i]]] <- 1
    }

    whoismyparent <- rep(NA, length(muts.per.hap))
    ref.ind <- which(unlist(lapply(muts.per.hap, function(x){length(x) == 0})))

    #all singletons are direct descendents of null. easy
    singles <- unlist(lapply(muts.per.hap, function(x){ c(length(x) == 1 )}))
    whoismyparent[ref.ind] <- 0
    whoismyparent[singles] <- ref.ind

    i <- 2
    while( sum(is.na(whoismyparent )) > 0){
        ntons <- which(unlist(lapply(muts.per.hap, function(x){ c(length(x) == i )})))

        #check if length(ntons) > 1
        for(j in 1:length(ntons)){

            #these are the current mutations
            mutsToMatch <- which(mut.relat[ntons[j],] == 1)


            #We want to know 1) how many mutations are shared with all other current children
            #seq 1 has A, B, C, E
            #seq 2 has B, C, D, E
                                        #seq 3 has A
            ##A B

            #seq 1 should be seq 3's child.

            #Let's list our options for parents
            potInds <- which(!is.na(whoismyparent))
            mutlist.pot <- apply(mut.relat[potInds,], 1, function(x){which(x == 1)})

            sharemuts <- unlist(lapply(mutlist.pot, function(x){
                needed.muts <- setdiff(mutsToMatch, x)
                unneeded.muts <- setdiff(x, mutsToMatch)
                                        #if no unneeded muts
                if(length(unneeded.muts) == 0){
                    return(length(needed.muts))
                }else{
                    return(100)
                }
            }))


            #sharemuts has the number of mutations necessary to convert parent to child
            #sharemuts == 0 means that the parent has a mutation that the parent doesn't have
            #this could be modified to allow for backmutation, but not currently

            #of the things that already have parents, do any of them share mutations?

            bestmatches <- which(sharemuts == min(sharemuts))


            #if sharemuts only has one thing, great, take that.
            if(length(bestmatches) == 1){
                whoismyparent[ntons[j]] <- potInds[bestmatches]
            }else if(length(bestmatches) > 1){
                 #otherwise, we will need to tiebreak.
                 
                #which haplotype had the highest frequency the earliest?
                #If there is only one timepoint, we need to do this differently...
                if(ncol(haps.over.time) == 1){
                    earliest.week <- 1
                }else{
                    earliest.week <- apply(haps.over.time[potInds[bestmatches],], 1, function(x){min(which(x != 0))})
                }
                matches.earliest <- bestmatches[which(earliest.week  == min(earliest.week))]
                if(length(matches.earliest) > 1){

                    #If multiple things match the earliest week, choose the one with the highest initial frequency
                    init.freq <- haps.over.time[potInds[matches.earliest], min(earliest.week)]
                    earliest.and.highest.freq <- matches.earliest[which(init.freq == max(init.freq))]

                    if(length(earliest.and.highest.freq) > 1){
                        #just choose one
                        #let's make it the first one. Why not?
                        whoismyparent[ntons[j]] <- potInds[earliest.and.highest.freq[1]]
                    }else{
                        whoismyparent[ntons[j]] <- potInds[earliest.and.highest.freq]
                    }
                }else{
                    whoismyparent[ntons[j]] <- potInds[matches.earliest]
                }
                
            }
        }
        i <- i + 1
    }

    #Let's also reverse this. Whoismychild
    whoismychild <- matrix(data = 0, nrow = nrow(mut.relat), ncol = nrow(mut.relat))
    #0 = col is not row's child
    #1 = col is row's child
    #there will certianly be a better way to do this, but let's just do it the easy way for now

    nochild <- setdiff(1:nrow(mut.relat), whoismyparent)
    whoismychild[nochild,] <- 0
    currchildren.thisround <- nochild

    while(length(currchildren.thisround) > 0){
        currchildren.nextround <- c()
        for(i in currchildren.thisround){
            newinf <- rep(0, nrow(mut.relat))
            newinf[i] <- 1
            whoismychild[whoismyparent[i],] <- whoismychild[whoismyparent[i],] +  whoismychild[i,] + newinf
            currchildren.nextround <- c(currchildren.nextround, whoismyparent[i])
        }
        currchildren.thisround <- unique(currchildren.nextround)
        anyzeros <- which(currchildren.thisround == 0)
        if(length(anyzeros) > 0){ currchildren.thisround <- currchildren.thisround[-anyzeros] } 
        #If we reach a point in which we have only one parent in question and it's the reference, we can stop
        if(length(currchildren.thisround) == 0){
            currchildren.thisround <- c()
        }
    }

    #I also want a direct children list, not an all children
    directchildren <- list()
    for(i in 1:length(whoismyparent)){
        directchildren[[i]] <- 0
    }
    for(i in 1:length(whoismyparent)){
        if(whoismyparent[i] != 0){
            directchildren[[whoismyparent[i]]] <- c(directchildren[[whoismyparent[i]]], i)      
        }
    }

    compute.polygon <- function(ind){

        poly <- list( y = c(rep(0, length(haps.over.time[ind,])), c(haps.over.time[ind,])), x = as.numeric(c( rev(colnames(haps.over.time)), (colnames(haps.over.time) ))))
        
        childlist <- list()
        if(length(directchildren[[ind]]) > 1){
                                        #Report the polygons of each of your children
            dirchild <- directchildren[[ind]]
            for(j in 2:length(dirchild) ){
                print(paste("calling", dirchild[j]))
                observedpoly <- compute.polygon(dirchild[j])
                childlist[[length(childlist) + 1]] <- observedpoly
            }
        }else{
            childlist <- NULL
        }
        c(list(mypoly = poly, children = childlist, index = ind))

    }

    compute.polygon.onecol <- function(ind){

        y = c(rep(0, length(haps.over.time[ind,])+1), rep(haps.over.time[ind,], 2))
        xval <- as.numeric(colnames(haps.over.time))
        x = c(xval, xval + 1, xval + 1, xval )
        poly <- list(y =y , x =x)

        childlist <- list()
        if(length(directchildren[[ind]]) > 1){
            #Report the polygons of each of your children
            dirchild <- directchildren[[ind]]
            for(j in 2:length(dirchild) ){
                print(paste("calling", dirchild[j]))
                observedpoly <- compute.polygon.onecol(dirchild[j])
                childlist[[length(childlist) + 1]] <- observedpoly
            }
        }else{
            childlist <- NULL
        }
        c(list(mypoly = poly, children = childlist, index = ind))

    }


    ind <- which(whoismyparent == 0)
    if(ncol(haps.over.time) == 1){
        dat <- compute.polygon.onecol(ind)
    }else{
        dat <- compute.polygon(ind)
    }

    return(list(dat = dat, haps.over.time = haps.over.time,
                mut.relat = mut.relat, all.comp.inds = all.comp.inds))

}


plotComp.new <- function(monk, subcomp, colors, labels){

    print(c(monk, subcomp))
    relInds <- retInd.id.comp(monk, subcomp) 

    #Import info from hierchy function
    hier <- setup.hierchy(relInds)
    haps.over.time <- hier$haps.over.time
    dat <- hier$dat
    mut.relat <- hier$mut.relat
    all.comp.inds <- hier$all.comp.inds

    colMap <- colors
    nameMap <- labels
    
    weeksizes <- apply(haps.over.time, 2, sum)
    mheight <- max(weeksizes)
    mheight <- 33

    xlims.with.room <- range(as.numeric(colnames(haps.over.time)))
    xlims.with.room[2] <- xlims.with.room[2]

    if(monk == "A99039" & subcomp == "LNRNA"){
        xlims.with.room[2] <- 38
    }
    
    uppery <- 30
    par(mar = c(1, 0, 0, 0))
    par(cex.main = 3)
    plot(0, type = "n", xlim = xlims.with.room, ylim = c(-1, uppery), xlab  = "", ylab = "", axes = FALSE, main = "", cex.lab = 2)

    comboend <- as.numeric(max(colnames(haps.over.time)))
    polygon(x = c(12, 20, 20, 12), y = c(-10, -10, 0, 0), col = "black")
    if(monk != "T98133"){ polygon(x = c(26, comboend, comboend, 26), y = c(-10, -10, 0, 0), col = "black")} 
    
    initshape <- dat$mypoly
    halflen <- length(dat$mypoly$y)/2

    initshape$y <- initshape$y/c(rev(weeksizes), weeksizes)*uppery
    #Plot reference
    polygon(initshape, col = "grey")

    yWT <- initshape$y[initshape$y > 0][1]
    xWT <- as.numeric(names(yWT)[1])

    shapeToAdd <- initshape$y[(halflen+1):(halflen * 2)]
    kidsToPlot <- dat$children
    #There are insufficient list functions to do what I want, so we'll have to keep a vector of list indices

    access.kidsToPlot <- 1:length(kidsToPlot)
    offset <- 0.25

    labels.to.eventually.plot <- c(xWT, yWT/2, xWT + offset, yWT/2,"WT")

    currep <- 1

    while(length(access.kidsToPlot) > 0){

        print(access.kidsToPlot)
        currInd <- kidsToPlot[[ access.kidsToPlot[1] ]]$index

        #Ok, this currInd index is with respect to its internal mut.relat
        kidsToPlot[[access.kidsToPlot[1]]]$index

        tmpyplot <- (kidsToPlot[[ access.kidsToPlot[1] ]]$mypoly$y)

        yToPlot <- tmpyplot/c(rev(weeksizes), weeksizes)*uppery + c(rev(shapeToAdd), shapeToAdd)
        xToPlot <- kidsToPlot[[ access.kidsToPlot[1] ]]$mypoly$x
        currlab <- labelmaker(names(which(mut.relat[currInd,] == 1)), all.comp.inds)
       
        polygon(xToPlot, yToPlot, col = colMap[nameMap == currlab], lwd = 1) #lwd = 2
        if(isDRM(currlab) == 1){
            polygon(xToPlot, yToPlot, col = "black", density = 4, border = NULL, lwd = 1) #lwd = 2
        }
        
       #first, we have to figure out where it appears first
        nomatchfirst <- which.min((rev(yToPlot[1:(halflen)]) ==   (yToPlot[(halflen + 1):(2*halflen)])))
        week.first.appeared <- as.numeric(colnames(haps.over.time))[nomatchfirst]
        y.big <- yToPlot[halflen + nomatchfirst]
        y.lit <- yToPlot[halflen - nomatchfirst + 1]
        yval <- (y.big +  y.lit)/2

        labels.to.eventually.plot <- rbind(labels.to.eventually.plot, c(week.first.appeared, yval, week.first.appeared + offset, yval ,currlab))
        
        shapeToAdd <- yToPlot[(halflen+1):(halflen * 2)]

         #Delete myself from the list (actually, the list index)
         #add any of my children to the BACK of the list, but the FRONT of the access order
        kidsList <- kidsToPlot[[access.kidsToPlot[1]]]$children
        access.kidsToPlot <- access.kidsToPlot[-1]
        if(length(kidsList) > 0){
            for(i in 1:length(kidsList)){
                #add a child to the end of kidsToPlot
                kidsToPlot[[length(kidsToPlot)+1]] <- kidsList[[i]]

                #add its index to the front of access.kidsToPlot
                access.kidsToPlot <- c(length(kidsToPlot), access.kidsToPlot)
            }
        }
        currep <- currep + 1

    }

    abline(v = as.numeric(colnames(haps.over.time)), col = "black", lwd = 3)

    for(i in 1:nrow(labels.to.eventually.plot)){
        
        ct <- labels.to.eventually.plot[i,5]
        print(ct)

        #Only print these labels
        goodLabs <- c("M184V",
                      "M184I",
                      "M184V,N255N",
                      "M184V,N255N,V179I",
                      "M184V,V179I",
                      "D177N,M184V,N255N",
                      "L187L",
                      "K103N-C,L187L",
                      "K103N-T",
                      "K103N-T,L187L",
                      "K103N-C",
                      "L74V",
                      "K103N-C,K249K,L187L,L205L,Q174R",
                      "K249K,L187L,L205L,Q174R",
                      "K249K,L187L,Q174R",
                      "K223K",
                      "L187L,M184V",
                      "K223K,L109L,M184V,V75L",
                      "L187L,M184V,M245T",
                      "L205L,M184V,R277K",
                      "K249K,L187L,L205L",
                      "I270I,L214F",
                      "L214F,M184V",
                      "L214F",
                      "L74V,T58T",
                      "K223K,L187L,M184I",
                      "L205L,M184V,R277K",
                      "G262G,L205L",
                      "K223K,L187L",
                      "L187L,M245T",
                      "M245T",
                      "L205L",
                      "K66K,L187L",
                      "K66R",
                      "WT")
        if(sum(ct == goodLabs) > 0){

            #The plotting functionality is extremely convoluted because I wanted many 
            # small changes to improve label visibility 
            dashlwd <- 2
            ca <- as.numeric(labels.to.eventually.plot[i,1:4 ])
            if(!(length(grep("DNA", subcomp)) > 0 & ct == "K103N-C,K249K,L187L,L205L,Q174R")){

                #This portion handles the text to the left of each line
                if(ca[1] == 26 & monk == "T98133"){
                    arrows(ca[1], ca[2], ca[1] - .15, ca[4], length = 0, lwd = dashlwd)
                    colorLabel(ca[1] - .5 - nchar(ct)*.28, ca[2], ct)
                }else if((ca[1] <= 29 & ca[1] >= 26) & monk == "A99165"){
                    if(ct == "L214F,M184V" ){
                        negoffset <- 4
                        arrows(ca[1], ca[2]-negoffset, ca[1] - .3, ca[4]-negoffset,
                               length = 0, lwd = dashlwd)
                        colorLabel(ca[1] - 1.1 - nchar(ct)*.3, ca[2]-negoffset, ct)
                    }else if(subcomp == "PLASMA" & ct == "K249K,L187L,L205L"){
                        negoffset <- 1
                        arrows(ca[1], ca[2]-negoffset, ca[1] - .3, ca[4]-negoffset,
                               length = 0, lwd = dashlwd)
                        colorLabel(ca[1] - 1.1 - nchar(ct)*.3, ca[2]-negoffset, ct)
                    }else if(subcomp == "GUTRNA" & ct == "K249K,L187L,L205L"){
                        negoffset <- -1
                        arrows(ca[1], ca[2]-negoffset, ca[1] - .3, ca[4]-negoffset,
                               length = 0, lwd = dashlwd)
                        colorLabel(ca[1] - 1.1 - nchar(ct)*.3, ca[2]-negoffset, ct)
                    }else if(subcomp == "GUTRNA" & ct == "K223K"){
                        
                    }else{
                        arrows(ca[1], ca[2], ca[1] - .3, ca[4], length = 0, lwd = dashlwd)
                        colorLabel(ca[1] - .9 - nchar(ct)*.3, ca[2], ct)
                    }
                }else if(ca[1] == 38 & monk == "A99039"){
                    arrows(ca[1], ca[2], ca[1] - .5, ca[4], length = 0, lwd = dashlwd)
                    colorLabel(ca[1] - 2.5 - nchar(ct)*.4, ca[2], ct)
                #This portion handles the text to the right of each line
                }else{
                    if(monk == "A99165" & ct == "K223K"){
                        negoffset <- -2
                        arrows(ca[1], ca[2]-negoffset, ca[3] , ca[4]-negoffset,
                               length = 0, lwd = dashlwd)
                        colorLabel(ca[3] + .15, ca[2]-negoffset, ct)
                    }else if((monk == "A99165" & ct == "K103N-C,K249K,L187L,L205L,Q174R") &
                             subcomp == "PLASMA"){
                        negoffset <- -1
                        arrows(ca[1], ca[2]-negoffset, ca[3] , ca[4]-negoffset,
                               length = 0, lwd = dashlwd)
                        colorLabel(ca[3] + .15, ca[2]-negoffset, ct)
                    }else if((monk == "A99165" & ct == "K249K,L187L,L205L,Q174R") &
                             subcomp == "PLASMA"){
                        negoffset <- .5
                        arrows(ca[1], ca[2]-negoffset, ca[3] , ca[4]-negoffset,
                               length = 0, lwd = dashlwd)
                        colorLabel(ca[3] + .15, ca[2]-negoffset, ct)
                    }else if((monk == "A99165" & ct == "K103N-C,L187L") &
                             subcomp == "PBMCDNA"){
                        negoffset <- 1.5
                        arrows(ca[1], ca[2]-negoffset, ca[3] , ca[4]-negoffset,
                               length = 0, lwd = dashlwd)
                        colorLabel(ca[3] + .15, ca[2]-negoffset, ct)
                    }else if((monk == "A99039" & ct == "K223K,L109L,M184V,V75L") &
                             subcomp == "LNDNA"){
                        arrows(ca[1], ca[2], ca[3] - 1 , ca[4],
                               length = 0, lwd = dashlwd)
                        colorLabel(ca[3] - 13, ca[2], ct)
                    }else{
                        mainoffset <- 0
                        if(monk == "A99039"){ mainoffset <- .3 }
                        arrows(ca[1], ca[2], ca[3] + mainoffset, ca[4], length = 0, lwd = dashlwd)
                        colorLabel(ca[3] + mainoffset + .15, ca[2], ct)
                    }

                }
            }
                ## print(c(monk, ct))
                ## #I want some extra custom plotting adjustments 
                ## if( monk == "A99165" & ct == "L214F,M184V"){
                ##     negoffset <- 10
                ##     arrows(ca[1], ca[2] - negoffset, ca[3], ca[4] - negoffset, length = 0, lwd = dashlwd)
                ##     colorLabel(ca[3] + .15, ca[2] - negoffset, ct)
                ##     print("asdfasdfasdfasdf")
                ## }else{

                ##     arrows(ca[1], ca[2], ca[3], ca[4], length = 0, lwd = dashlwd)
                ##     colorLabel(ca[3] + .15, ca[2], ct)

                ## }

        }
    }
    plotTops <- apply(haps.over.time, 2, sum)
    mtext(paste(plotTops), side = 3, line = 0, at = as.numeric(names(plotTops)), cex = .5)

    treatsize <- 1.25
    if(monk == "A99165"){
        text(16, -1, "EFV", col = "white", cex = treatsize)    
        text(27.5, -1, "Rx1", col = "white", cex = treatsize)
    }
    if(monk == "T98133"){
        text(16, -1, "FTC", col = "white", cex = treatsize)    
    }
    if(monk == "A99039"){
        text(16, -1, "FTC", col = "white", cex = treatsize)    
        text(32, -1, "Rx2 ", col = "white", cex = treatsize)
    }
}

#Return indices for a given week, id and comparture
returnInds <- function(id, week, compartment){
    c(intersect(intersect(which(monkid == id), which(samp.loc == compartment)), which(weeks == week)))
}

#Return indices for a given id and comparture
retInd.id.comp <- function(id, compartment){
    c(intersect(which(monkid == id), which(samp.loc == compartment)))
}

#Is it a synonymous mutation?
is.syn <- function(lab){
    letters <- strsplit(lab, split = "")[[1]]
    if(letters[1] == letters[nchar(lab)]){
        return(TRUE)
    }else{
        return(FALSE)
    }
}

#Plot the labels
colorLabel <- function(xplot, yplot, lab){
    sizeval <- .75
    rval <- .1
    allmuts <- strsplit(lab, split = ",")[[1]]
    for(i in length(allmuts):1){
        currmut <- allmuts[i]
        if(isDRM(currmut) == 1){
            currcol <- "pink"
        }else if(is.syn(currmut)){
            currcol <- "white"
        }else{
            currcol <- "white"
        }
        splitaroundmut <- strsplit(lab, split = currmut)[[1]]
        textToPlot <- paste0(splitaroundmut, currmut)[1]
        shadowtext(xplot, yplot, textToPlot, col = currcol, adj = c(0,0.5), cex = sizeval, r = rval)
    }
}


