#!/usr/bin/perl

#These are the compartments, times and macaques to be checked
@compartments = ("PLASMA", "GUT", "LN", "VAG", "PBMC");
@times = (1, 2, 3, 4, 5);
@monkeys = ("T98133", "A99165", "A99039", "A01198");

open(OUTFI, ">>","../out/sm/SM.out.txt") or die;

foreach $monkey (@monkeys){
    foreach $time (@times){
	for (my $c1=0; $c1 <  ( @compartments) ; $c1++){
	    for (my $c2=$c1 ; $c2 < @compartments; $c2++){
		$comp1 = @compartments[$c1];
		$comp2 = @compartments[$c2];

		print("$monkey $time $comp1 $comp2 $compsize1 $compsize2 \n");

		foreach $downSampInd ((3, 7, 10, 15, 20)){

		    print("$downSampInd\n");
		    #in paper, trialnum = 100
		    for (my $trialnum = 1; $trialnum <= 10; $trialnum++){
			
			#Let's initialize all these variables
			$p = NULL;

		    #construct the paths to the input and output files

			$infile = "../tmp/weekly_trees/auto/$comp1\-$comp2\/n$downSampInd\/$monkey\.$time\.$trialnum";
			$endfile = "../tmp/hyphy/auto/$comp1\-$comp2\/n$downSampInd\/$comp1\-$comp2\/$monkey\.$time\.$trialnum";
			$regexp = $comp1;

			#In paper, numTrials = 1000
			$numTrials = 10;

			if (-f $infile){

			    print("running hyphy");
		    #Run HYPHY
		    #We provide an input tree (of both compartments)
		    #And an output file
		    #A regexp to match some but not all tree entries (i.e., "PLASMA"
		    #And a number of trials
#		      system("HYPHYMP ../dat/SlatkinMaddison_af.bf <<< \$\'$infile\n$endfile\n$regexp\n$numTrials' > tmpfile.txt");
			    system("(echo $infile; echo 2;  echo $regexp; echo y; echo $regexp; echo Other; echo tmpfile.txt; echo 2; echo $numTrials;) | hyphymp ../dat/SlatkinMaddison.bf > tmpout.txt");

			    open(TMP, "<","tmpout.txt") or die;
			    @tmp = <TMP>;
			    @split = split(" = ", "@tmp[-2]");
			    $p = @split[1];
			    chomp($p);

			    #Open that file and parse it
			    #Print all the results 

			    if($p eq ''){ $p = "NULL";  }
			    $finalstring = "$monkey\t$time\t$comp1\t$comp2\t$downSampInd\t$trialnum\t$p\n";
			    print OUTFI $finalstring;

			}
		    }
		}
	    }
	}
    }
}


close(OUTFI);
