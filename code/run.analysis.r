#Run all analysis 

#Read in all data and shared functions and loads all required packages
source("read.in.data.r")
source("shared-functions.r")

#produce the table describing the frequencies of samples by time, location and macaque
#results in the file: out/tables/full.html
source("produce-sample-table.r")

#estimate the selection coefficients of the mutations in each compartment and print
#results in the file: out/selection.coefs.txt
source("compute-s.r")

#estimate theta in each macaque (by the three different means described in the text) and produce a table
#results in the file: out/tables/theta.html
source("compute-theta.r")

#Perform the compartment v compartment phi_st analysis and plot the results
#together, these functions result in the files: out/graphs/bettercomp-new.pdf
#                                               out/graphs/compoverall.pdf
source("kst-run-comps.r")
source("kst-plot-comps.r")

#Perform the compartment over time phi_st analysis and plot the results
#together, these functions result in the files: out/graphs/betterovertime.pdf
source("kst-run-time.r")
source("kst-plot-time.r")

#Plot the vDNA and vRNA overlap plots
#results in the file: out/graphs/RNADNAoverlap.pdf
source("rna.dna.overlap.r")

#Plot the percentage drug resistant plot (by compartment)
#results in the file: out/graphs/percent_drugresistant.pdf
source("percentDR.r")

#Plot the drug resistance in the plasma (and the VL)
#results in the file: out/graphs/VL.pdf
source("DR.and.VL.fig.r")

#Look at diversity over time
#results in the file: out/graphs/diversity.pdf
#                     out/tables/piovertime.html
source("pi.over.time.r")

#Look at changes in Kst over time
#results in the file: out/graphs/migvsel.pdf
source("plot-kst-overtime.r")

#Make the muller plots
#results in the files: out/graphs/T98133.pdf
#                      out/graphs/A99165.pdf
#                      out/graphs/A99039.pdf
source("mullerplots.r")



#Below this point, we will be setting up and running the Slatkin-Maddison 
#test in HyPhy. Requires HyPhy to be installed (reachable via the command HYPHY)
#and perl.

source("write.weekly.FASTAs.r")
system("perl SM.pairwise.runner.pl")
#note: plotFigs-hyphy.r produces a supplemental figure that requires it be run after
# kst-run-comps.r
source("plotFigs-hyphy.r")
