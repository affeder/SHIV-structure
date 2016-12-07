
#This code produces a table with all of the information about the number of samples stratified by 
# macaque, timepoint and sample
fulltab <- rbind(table(weeks[which(monkid == "T98133")], samp.loc[which(monkid == "T98133")]),
                 table(weeks[which(monkid == "A99039")], samp.loc[which(monkid == "A99039")]),
                 table(weeks[which(monkid == "A99165")], samp.loc[which(monkid == "A99165")]),
                 table(weeks[which(monkid == "A01198")], samp.loc[which(monkid == "A01198")]))

printtab <- data.frame( row = rownames(fulltab), fulltab)
print.xtable(xtable(printtab),  type = "html", paste("../out/tables/full.html", sep = ""), include.rownames=FALSE)

