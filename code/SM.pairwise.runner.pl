#!/usr/bin/perl

#These are the compartments, times and macaques to be checked
@compartments = ("PLASMA", "GUT", "LN", "VAG", "PBMC");
@times = (1, 2, 3, 4, 5);
@monkeys = ("T98133", "A99165", "A99039", "A01198");


#This is an info file with numbers of sequences
open(INF, "<","../tmp/pairwise.details.txt") or die;
my @info = <INF>;
close(INF);

open(OUTFI, ">","../out/SM.allresults.txt") or die;

my $ind = 0;
foreach $monkey (@monkeys){
    foreach $time (@times){
	for (my $c1=0; $c1 <  ( @compartments) - 1; $c1++){
	    for (my $c2=$c1 + 1; $c2 < @compartments; $c2++){

		chomp(@info[$ind]);
		@toParse = split("\t", @info[$ind]);

		$comp1 = @toParse[2];
		$comp2 = @toParse[3];
		$compsize1 = @toParse[4];
		$compsize2 = @toParse[5];

#What stage is it? 1 = 12/13, 2 = 15/16, 3 = 20/21, 4 = 26/27, 5 = end
#what monkey is it?
#What is compartment 1 name
#What is compartment 2 name

#Let's try to read in the file, so we can determine the number of sequences

#How many sequences in comparmtent 1
#How many sequences in compartment 2
#If the number of sequences in compartments 1 and 2 are both >= 1

		#Let's initialize all these variables
		$p = NULL;
		$IQR1 = NULL;
		$IQR2 = NULL;
		$median = NULL;
		$infe = NULL;

		if($compsize1 > 0  && $compsize2 > 0){

		    #construct the paths to the input and output files

		    $infile = "../tmp/weekly_trees/auto/$comp1\-$comp2\/$monkey\.$time";
		    $endfile = "../tmp/hyphy/auto/$comp1\-$comp2\/$comp1\-$comp2\/$monkey\.$time";
		    $regexp = $comp1;
		    $numTrials = 10000;
		      
		    #Run HYPHY
		    #We provide an input tree (of both compartments)
		    #And an output file
		    #A regexp to match some but not all tree entries (i.e., "PLASMA"
		    #And a number of trials
		    system(" HYPHYMP ../dat/SlatkinMaddison_af.bf <<< \$\'$infile\n$endfile\n$regexp\n$numTrials' > tmpfile.txt");

		    #Open that file and parse it
		    open(TMP, "<","tmpfile.txt") or die;
		    @randTrials = ();
		    while(<TMP>){
#			print $_;
			$tmpstring = $_;

			#if a line ends with '}', do something extra
			if($tmpstring =~ /\"([0-9]+)\":([0-9]+)/  ){
			    $toRep = $1;
			    $repNum = $2;
			    push(@randTrials, ($toRep) x $repNum );
			}
			#if p=, start differently
			if($tmpstring =~ '^p.'){
			    @toSplit = split("=", $tmpstring);
			    $p = @toSplit[1];
			    chomp($p);
			}
			#if inferred=, start differently
			if($tmpstring =~ '^Inferred.'){
			    @toSplit = split("=", $tmpstring);
			    $infe = @toSplit[1];
			    chomp($infe);
			}
		    }
		    @sortedTrials = sort {$a <=> $b} @randTrials;
		    $IQR1 = @sortedTrials[$numTrials*.025];
		    $IQR2 = @sortedTrials[$numTrials*.975];
		    $median = @sortedTrials[$numTrials*.5];

		    close(TMP);

		}
		#Print all the results 
		$finalstring = join("\t", @toParse)."\t$infe\t$p\t$IQR1\t$median\t$IQR2\n";
		print OUTFI $finalstring;
		$ind = $ind + 1;
	    }
	}
    }
}

close(OUTFI);
