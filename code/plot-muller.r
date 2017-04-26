source("muller-plot-functions.r")

megaplot.gen("T98133",
             c("PLASMA","PBMCDNA", "LNRNA", "LNDNA", "GUTRNA", "GUTDNA", "VAGRNA", "VAGDNA"),
             c("Plasma/PBMC", "LN", "Gut", "Vagina"), c("A. vRNA", "B. vDNA"))

megaplot.gen("A99165",
             c("PLASMA","PBMCDNA", "LNRNA", "LNDNA", "GUTRNA", "GUTDNA"),
             c("Plasma/PBMC", "LN", "Gut"), c("A. vRNA", "B. vDNA"))

megaplot.gen("A99039",
             c("PLASMA","PBMCDNA", "LNRNA", "LNDNA", "GUTRNA", "GUTDNA"),
             c("Plasma/PBMC", "LN", "Gut"), c("A. vRNA", "B. vDNA"))





