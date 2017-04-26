require(ggtree)

if(!dir.exists("../out/beast_out/")){ dir.create("../out/beast_out/") }

#Warning: there are a lot of (acceptable) warnings here

allMonksDNAorRNA <- c("A99165_RNA", "A99165_DNA",
                      "A99039_RNA", "A99039_DNA",
                      "T98133_RNA", "T98133_DNA",
                      "A01198_RNA", "A01198_DNA")


for(plotName in allMonksDNAorRNA){


    sto <- read.beast(paste("../tmp/beast/mcc/",plotName,".mcc", sep = ""))

    cols <- gsub(".*_#", "#", sto@phylo$tip.label)
    muts <- gsub("_#.*|.*_[0-9]+_[0-9]+_", "", sto@phylo$tip.label)
    groupInfo <- gsub("_.*", "", sto@phylo$tip.label)

    classList <- list()
    for(i in unique(muts)){
        classList[[ i ]] <- sto@phylo$tip.label[which(muts == i)]
    }

    classList[[ '0' ]] <- NA

    colKey <-  unique(cols)
    names(colKey) <- unique(muts)
    colKey[which(names(colKey) == "WT")] <- "lightgrey"

    colKey <- c(colKey, '0' = "black")

    positions <- sort(unique(as.numeric(gsub('[A-Z]|\\*', "",
                                             unique(unlist(strsplit(muts, split = "-")))))))

    seqDat <- matrix(data = "WT", nrow = length(sto@phylo$tip.label), ncol = length(positions))
    for(i in 1:length(muts)){
        seqmuts <- strsplit(muts[i], "-")[[1]]
        if(seqmuts[1] != "WT"){
            for(k in 1:length(seqmuts)){
                pos <- as.numeric(gsub("[A-z]|\\*","", seqmuts[k]))
                if(!is.na(pos)){
                    toAdd <- which(positions == pos)
                    if(length(toAdd) > 0){
                        lastval <- substr(seqmuts[k], nchar(seqmuts[k]) , nchar(seqmuts[k]))
                        firstval <- substr(seqmuts[k], 1 , 1)
                        if(firstval == "K" & lastval == "N"){
                            seqDat[i, toAdd] <- paste("K103N-", seqmuts[k + 1], sep = "")
                        }else if((firstval == "M" & lastval == "V") & pos == 184){
                            seqDat[i, toAdd] <- "M184V"
                        }else if(firstval == "M" & lastval == "I"){
                            seqDat[i, toAdd] <- "M184I"
                        }else if(lastval == "X" ){
                            seqDat[i, toAdd] <- "Nonsense/Ambiguous"
                        }else if(lastval == firstval){
                            seqDat[i, toAdd] <- "Synonymous"
                        }else if(lastval == "*"){
                            seqDat[i, toAdd] <- "Stop"
                        }else{
                            seqDat[i, toAdd] <- "Nonsynonymous"
                        }
                    }
                }
            }
        }
    }
    seqDat <- as.data.frame(seqDat)
    colnames(seqDat) <- positions
    rownames(seqDat) <- sto@phylo$tip.label

    seqDat <- cbind(groupInfo, " ", seqDat)
    names(seqDat)[1] <- "Loc"
    names(seqDat)[2] <- " "


    rnacols <- c("#228833", "#EE6677", "#AA3377", "#4477AA", "#66CCEE")
    m184cols <- tcol(brewer.pal(3, "Set1")[1], seq(100, 0, by = -20))
    k103cols <- tcol(brewer.pal(3, "Set1")[2], seq(100, 0, by = -20))
    stopcol <- brewer.pal(6, "Set1")[6]
    othercols <- brewer.pal(6, "Set1")[6]

    vals <- unique(c(as.matrix(seqDat)))
    #Order: "Plasma", "PBMC",  "LN", "Gut", "Vagina"
    if(length(grep("PLASMA", vals)) > 0){ names(vals)[grep("PLASMA", vals)] <- rnacols[1]  }
    if(length(grep("PBMC", vals)) > 0){ names(vals)[grep("PBMC", vals)] <- rnacols[2]  }
    if(length(grep("LN", vals)) > 0){ names(vals)[grep("LN", vals)] <- rnacols[3]  }
    if(length(grep("GUT", vals)) > 0){ names(vals)[grep("GUT", vals)] <- rnacols[4]  }
    if(length(grep("VAG", vals)) > 0){ names(vals)[grep("VAG", vals)] <- rnacols[5]  }

    if(length(grep("Nonsynonymous", vals)) > 0){
        names(vals)[grep("Nonsynonymous", vals)] <- "black"  }
    if(length(grep("Synonymous", vals)) > 0){
        names(vals)[grep("Synonymous", vals)] <- "grey50"  }
    if(length(grep("WT", vals)) > 0){
        names(vals)[grep("WT", vals)] <- "grey85"  }
    if(length(grep("Nonsense", vals)) > 0){
        names(vals)[grep("Nonsense", vals)] <- "beige"  }
    if(length(grep("Stop", vals)) > 0){
        names(vals)[grep("Stop", vals)] <- stopcol  }
    if(length(grep("M184V", vals)) > 0){
        names(vals)[grep("M184V", vals)] <- m184cols[6]  }
    if(length(grep("M184I", vals)) > 0){
        names(vals)[grep("M184I", vals)] <- m184cols[3]  }
    if(length(grep("K103N-C", vals)) > 0){
        names(vals)[grep("K103N-C", vals)] <- k103cols[6]  }
    if(length(grep("K103N-T", vals)) > 0){
        names(vals)[grep("K103N-T", vals)] <- k103cols[3]  }
    names(vals)[which(vals == " ")] <- "white"

    allcols <- names(vals)
    names(allcols) <- vals
    sto <-groupOTU(sto, focus = classList)

    modifier <- 1.15
    png(paste("../out/beast_out/", plotName,".png", sep = ""), 2250*modifier, width = 2625*.2, res = 300)
    #We want our tree branches encoded by color based on the group list defined by classList
    p <-  ggtree(sto, aes(color = group ), lwd = .4) + 
        #We want the colors to be drawn from colKey
        scale_color_manual(values = c(colKey, "black", "black") )+ theme_tree2() + 
            geom_tiplab(align=TRUE, linesize = .025, col = "black", size = 0, linetype = "solid") 
    #Our heat map will have the nucleotides and the sample location
    pp<- gheatmap(p, seqDat, offset = 0.25, width=0.1, font.size=1.5, colnames_angle=-90,
                  hjust=-0.1, colnames = TRUE) + scale_fill_manual(values= allcols, breaks = vals)
    #And the legend throws off the sizing, so we print without the legend
    print(pp + theme(legend.position = "none"))
    dev.off()

    modifier <- 1.15
    png(paste("../out/beast_out/", plotName,"_long.png", sep = ""), 2250*modifier, width = 2625*.4, res = 300)
    #We want our tree branches encoded by color based on the group list defined by classList
    p <-  ggtree(sto, aes(color = group ), lwd = .4) + 
        #We want the colors to be drawn from colKey
        scale_color_manual(values = c(colKey, "black", "black") )+ theme_tree2() + 
            geom_tiplab(align=TRUE, linesize = .025, col = "black", size = 0, linetype = "solid") 
    #Our heat map will have the nucleotides and the sample location
    pp<- gheatmap(p, seqDat, offset = 0.25, width=0.1, font.size=1.5, colnames_angle=-90,
                  hjust=-0.1, colnames = TRUE) + scale_fill_manual(values= allcols, breaks = vals)
    #And the legend throws off the sizing, so we print without the legend
    print(pp + theme(legend.position = "none"))
    dev.off()

    
    #Since we still want access to the legend in compiling composite figures, 
    # we also plot the same thing again with the legend and a much wider number of pixels
    png(paste("../out/beast_out/", plotName,"-legends-1.png", sep = ""),
        2250*modifier, width = 2625, res = 300)
    p <-  ggtree(sto, aes(color = group ), lwd = .4) + 
        scale_color_manual(values = c(colKey, "black") )+ theme_tree2() + 
            geom_tiplab(align=TRUE, linesize = .025, col = "black", size = 0, linetype = "solid") #+ 
    pp<- gheatmap(p, seqDat, offset = 0.25, width=0.5, font.size=2, colnames_angle=-90,
                  hjust=-0.1, colnames = TRUE) + scale_fill_manual(values= allcols, breaks = vals)
    print(pp + theme(legend.position = "none") )
    dev.off()

    png(paste("../out/beast_out/", plotName,"-legends-2.png", sep = ""),
        2250*modifier, width = 2625*.3, res = 300)
    p <-  ggtree(sto, aes(color = group ), lwd = .4) + 
        scale_color_manual(values = c(colKey, "black") )+ theme_tree2() + 
            geom_tiplab(align=TRUE, linesize = .025, col = "black", size = 0, linetype = "solid") #+ 
    pp<- gheatmap(p, seqDat, offset = 0.25, width=0.5, font.size=2, colnames_angle=-90,
                  hjust=-0.1, colnames = TRUE) + scale_fill_manual(values= allcols, breaks = vals)
    print(pp + theme(legend.position = "none") )
    dev.off()

}






