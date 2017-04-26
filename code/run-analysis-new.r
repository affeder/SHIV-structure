#Run all analysis 

#NOTE:
#These lines must be uncommented and run once:
#It is expecting the files ../dat/RT-SHIV-RNA.fa and ../dat/RT-SHIV-DNA.fa
source("set-up-directory-structure.r")
source("one-time-processing.r")

#Read in all data and shared functions and loads all required packages
source("read-in.r")
source("shared-functions.r")

################################################################################
#This section produces tables unconnected to other analyses
################################################################################

#produce the table describing the frequencies of samples by time, location and macaque
#results in the file: out/tables/full.html
source("produce-sample-table.r")

#estimate the selection coefficients of the mutations in each compartment and print
#results in the file: out/selection.coefs.txt
source("compute-s.r")

#estimate theta in each macaque (by the three different means described in the text) and produce a table
#results in the file: out/tables/theta.html
source("compute-theta.r") #

################################################################################
#We perform three tests for compartmentalization: KST, AMOVA and Slatkin-Maddison (SM)
################################################################################

#1) KST analysis
#This produces a written file (kstout.txt) that has 1000 iterations
# of subsampling to compartment size 10
source("run-kst-ds10.r") 
#This produces a written file (kstout.power.txt) that has fewer (100) interations
# of subsampling to compartment sizes 3, 7, 10, 15, 20 for the power analysis
source("run-kst-power.r") 
#This produces a written file (kstout.time.txt) that has 1000 iterations of 
# subsampling to compartment size 10 from consecutive timepoints
source("run-kst-time.r") 

#2) AMOVA analysis
#This produces a written file (amovaout.txt) that has 1000 iterations of 
# subsampling to a compartment size of 10
source("run-amova-ds10.r") 
#This produces a written file (amovaout.power.txt) that has fewer (100) interations
# of subsampling to compartment sizes 3, 7, 10, 15, 20 for the power analysis
source("run-amova-power.r") 

#3) Slatkin-Maddison analysis

#In this next section, we will be setting up and running the Slatkin-Maddison 
#test in HyPhy. Requires HyPhy to be installed (reachable via the command HYPHYMP)
#and perl. 
#There is a place in run-SM.pl that points to the batch file SlatkinMaddison.bf. 
#You can either change the path to point to the path file in your HyPhy install,
# or you can move a batchfile of SlatkinMaddison.bf into your dat file

source("write-trees-for-SM.r")
system("perl run-SM.pl")

################################################################################
#This next section plots compartmentalization results:
################################################################################

#This file reads in and reformats the previously run analysis, and must be run before
# any of the plot sections
source("load-comp-results.r")

#Produce Figure 6 for all three tests:
#results in the file: out/graphs/F6-KST.pdf
#                     out/graphs/F6-amova.pdf
#                     out/graphs/F6-SM.pdf
source("plot-f6.r")

#Produce Figure 7:
#results in the file: out/graphs/F7.pdf
source("plot-f7.r") 

#Produce Figure 8:
#results in the file: out/graphs/F8.pdf
source("plot-f8.r")

#Produce power analysis plots:
source("plot-S9.S9.S10.r")

################################################################################
#This next section plots non-compartmentalization results:
################################################################################


#Plot the drug resistance in the plasma (and the VL)
#results in the file: out/graphs/f1.pdf
source("plot-f1.r")

#Plot the percentage drug resistant plot (by compartment)
#results in the file: out/graphs/f2.pdf
source("plot-f2.r")

#Look at diversity over time
#results in the file: out/graphs/s1.pdf
#                     out/tables/piovertime.html
source("plot-s1.r") 

#Plot the vDNA and vRNA overlap plots
#results in the file: out/graphs/S7.pdf
source("plot-s7.r")


#Make the muller plots
#results in the files: out/graphs/T98133.pdf
#                      out/graphs/A99165.pdf
#                      out/graphs/A99039.pdf
source("plot-muller.r")

################################################################################
#This section sets up files to be run in BEAST and plots them.
################################################################################


#In this section, we will be setting up FASTA files to be run through BEAST
# and then plotting the resulting files. 
source("format-for-beast.r")

#Here, we run BEAST with parameters as detailed in the paper
# separately by macaque and by RNA and DNA. mccs are compiled by TreeAnnotator
#This following plotting code will look for mccs entitled monkey-[DNA|RNA].mcc
#In the following folder: ../tmp/beast/mcc/
#source("plot-phylos.r")
