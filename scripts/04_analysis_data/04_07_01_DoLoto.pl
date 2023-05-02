#!/usr/bin/perl
$| = 1;
use strict;
#~ use warnings;

my $f = '';#DREAM5_NetworkInference_GoldStandard_Network3.tsv
my $ff = '';#DREAM5_NetworkInference_inferred_Network3.tsv
my $threshold = 0;
my $case = 2;
my $folder = '.';
my $outfile = 'outfile.txt';
my $fg = '';
my $ffg = '';
my $id = '';


my %net = ();#tf->gene
%{$net{"true"}} = ();#tf->gene 1, 0 else
%{$net{"pred"}} = ();#p(tf->gene)>threshold 1, 0 else
#~ my %inDegree = ();#indegree
#~ %{$inDegree{"true"}} = ();
#~ %{$inDegree{"pred"}} = ();
my %outDegree = ();#outdegree
%{$outDegree{"true"}} = ();
%{$outDegree{"pred"}} = ();
my @tf = ();#true TFs
my @tf0 = ();#pred TFs
my %tfall= ();#all TFs
my %nodeall= ();#all TFs
my $TNG = 0;
my %done = ();#all nodes 
%{$done{"true"}} = ();
%{$done{"pred"}} = ();
%{$done{"all"}} = ();
my @keysdone = ();
my %dk = ();
my $node = 0;
my %res = ();#motives
%{$res{"true"}} = ();
%{$res{"pred"}} = ();
#~ %{$resmot} = ();
my %resmot = ();
#~ %{$resmotNode} = ();
my %resmotNodeTrue = ();
my %resmotNodePred = ();
@{$res{"gm"}} = (0,0,0,0,0);#global motives
@{$res{"global"}} = (0,0,0,0);#single edge res
my @countgl = (0,0);
my %countt = ();#motive counts
my %countp = ();
my %countlist = ();#to store lists of graphlets TP, FP and FN
my %outnet = ();
my %outnetnodes = ();

my $colorTFTP = '#ff6f00';
my $colorTFFN = '#bf008f';
my $colorTFFP = '#ffb300';

my $colorEFTP = '#0276cf';
my $colorEFFN = '#3b3176';
my $colorEFFP = '#ffffff';


my $coloredgeTP = '#000000';
my $coloredgeFN = '#bf008f';
my $coloredgeFP = '#ffb300';

my $outtxt = '';
my $outtxttable = '';

for my $i(0..12){
	$countt{$i} = 0;
	$countp{$i} = 0;
	%{$resmot{$i}} = ("TP",0,"FP",0,"TN",0,"FN",0);
}




#~ $f = "/home/alberto/Dropbox/Work/Network_predictions/DREAM5_NetworkInference_GoldStandard_Network3.tsv";		
#~ $ff = "/home/alberto/Dropbox/Work/Network_predictions/DREAM5_NetworkInference_GoldStandard_Network3.tsv";
#~ $ff = "/home/alberto/Dropbox/Work/Network_predictions/Community_integration/DREAM5_NetworkInference_Community_Network3.txt";
#~ $f = "/home/alberto/Dropbox/Work/Network_predictions/test_gold01.txt";
#~ $ff = "/home/alberto/Dropbox/Work/Network_predictions/test_gold01.txt";

sub getopt{
	my $usage = "
	runit options
	-g file\t: reference file (tsv TF\\tGENE\\t1)
	-i file\t: input file (tsv TF\\tGENE\\t1)
	-t num\t: threshold to consider an edge in input file as existing (default 0)
	-c num\t: case 1 (ref and input contain all TP and TN, only those interactions in ref are considered)
			 case 2 (ref and input contain only all TP, negatives are assumed to happen among all nodes without a true between them, default)
	-f folder\t: folder where the output file will be placed
	-o name\t: name of output file
	
	-h||--help : print help

";
	while(@ARGV){
		if($ARGV[0] eq '-g'){
			$f = $ARGV[1];
			shift @ARGV;
			shift @ARGV;
		}
		if($ARGV[0] eq '-i'){
			$ff = $ARGV[1];
			shift @ARGV;
			shift @ARGV;
		}
		if($ARGV[0] eq '-t'){
			$threshold = $ARGV[1];
			shift @ARGV;
			shift @ARGV;
		}
		if($ARGV[0] eq '-c'){
			$case = $ARGV[1];
			shift @ARGV;
			shift @ARGV;
		}
		if($ARGV[0] eq '-f'){
			$folder = $ARGV[1];
			shift @ARGV;
			shift @ARGV;
		}
		if($ARGV[0] eq '-o'){
			$outfile = $ARGV[1];
			shift @ARGV;
			shift @ARGV;
		}
		if($ARGV[0] eq '-h' || $ARGV[0] eq '--help'){
			die $usage;
		}
	}	
	
	#~ if($f eq '' || $ff eq ''){
	if($f eq ''){
		die $usage;
	}
	elsif($case != 1 && $case != 2){
		die $usage;
	}
	else{
		my @tfg = split('/',$f);
		$fg = $tfg[$#tfg];
		@tfg = split('/',$ff);
		$ffg = $tfg[$#tfg];
		#~ print "\nGOLD       \t:\t$f\nINPUT      \t:\t$ff\nThreshold\t:\t$threshold\nCase       \t:\t$case\n\n";
		print "\nREFERENCE       \t:\t$fg\nINPUT            \t:\t$ffg\nThreshold        \t:\t$threshold\nCase             \t:\t$case\n\n";
		$outtxttable .= "\nREFERENCE       \t:\t$fg\nINPUT            \t:\t$ffg\nThreshold        \t:\t$threshold\nCase            \t:\t$case\n\n";
		if($ff ne ''){
			$id = "$fg\_$ffg";
		}
		else{
			$id = "$fg";
		}
	} 
}

sub compute_stats{
	my @tmp = ($_[0],$_[1],$_[2],$_[3]);
	for (my $ii = 4; $ii < 23; $ii++){
		$tmp[$ii] = 0;
	}
	if($tmp[0] + $tmp[3]){
		$tmp[4] = ($tmp[0] / ($tmp[0] + $tmp[3]));#tp/(tp+fn),R
	}
	if($tmp[2] + $tmp[1]){
		$tmp[5] = ($tmp[2] / ($tmp[2] + $tmp[1]));#tn/(tn+fp),TNR
	}
	if($tmp[0] + $tmp[1]){
		$tmp[6] = ($tmp[0] / ($tmp[0] + $tmp[1]));#tp/(tp+fp),P
	}
	if($tmp[2] + $tmp[3]){
		$tmp[7] = ($tmp[2] / ($tmp[2] + $tmp[3]));#tn/(tn+fn),NPV
	}
	if($tmp[2] + $tmp[1]){
		$tmp[8] = ($tmp[1] / ($tmp[2] + $tmp[1]));#fp/(fp+tn),FPR
	}
	if($tmp[0] + $tmp[1]){
		$tmp[9] = ($tmp[1] / ($tmp[0] + $tmp[1]));#fp/(fp+tp),FDR
	}
	if($tmp[0] + $tmp[3]){
		$tmp[10] = ($tmp[3] / ($tmp[0] + $tmp[3]));#fn/(fn+tp),FNR
	}
	if($tmp[0] + $tmp[1] + $tmp[2] + $tmp[3]){
		$tmp[11] = (($tmp[0] + $tmp[2]) / ($tmp[0] + $tmp[1] + $tmp[2] + $tmp[3]));#(tp+tn)/(tp+fp+tn+fn),ACC
	}
	if($tmp[4] + $tmp[6]){
		$tmp[12] = (2* ($tmp[4] * $tmp[6]) / ($tmp[4] + $tmp[6]));#(2PR)/(P+R),F1
	}
	if(0 != (($tmp[0] + $tmp[1])  * ($tmp[0] + $tmp[3]) * ($tmp[2] + $tmp[1])  * ($tmp[2] + $tmp[3]))){
		$tmp[13] = (($tmp[0] * $tmp[2] - $tmp[1] * $tmp[3]) / sqrt(($tmp[0] + $tmp[1])  * ($tmp[0] + $tmp[3]) * ($tmp[2] + $tmp[1])  * ($tmp[2] + $tmp[3])));#MCC
	}
	$tmp[14] =  ($tmp[4] + $tmp[5] - 1);#R+TNR-1, infor
	$tmp[15] =  ($tmp[6] + $tmp[7] - 1);#P+NPV-1,marked
	if($tmp[1] + $tmp[3] + $tmp[0]){
		$tmp[16] = ($tmp[0]/($tmp[1] + $tmp[3] + $tmp[0]));#TP/(TP+FN+FP),J
	}
	if($tmp[0] + $tmp[0] + $tmp[1] + $tmp[3]){
		$tmp[17] = ($tmp[0]  / ($tmp[0] + $tmp[0] + $tmp[1] + $tmp[3]));#(tp)/(2tp+tn+fn),SD
	}
	if($tmp[1] + $tmp[3]){
		$tmp[18] = ($tmp[0]  / ($tmp[1] + $tmp[3]));#(tp)/(fp+fn),K1
	}
	$tmp[19] =  ($tmp[6] + $tmp[7] )/2;#(R+P)/2,K2
	my $temp = ($tmp[0] + $tmp[1]) * ($tmp[0] + $tmp[3]);
	if($temp > 0){
		$tmp[20] = ($tmp[0]  / sqrt($temp));#(tp)/sqrt((tp+fn)(tp+fp)),O
	}
	if($temp > 0){
		$tmp[21] = (($tmp[0] * $tmp[0])  / $temp);#(tp*tp)/(tp+fn)(tp+fp),CR
	}
	$tmp[22] = ($tmp[1] + $tmp[3]);#tp+fn,Hammin
	return @tmp;
}

sub dotriplet_noij{#if($net{"true"}{$k[$i]}{$k[$j]}) then $net{"true"}{$k[$i]}{$k[$l]} && $net0{"true"}{$k[$j]}{$k[$l]} && $net0{"true"}{$k[$l]}{$k[$j]} && $net0{"true"}{$k[$j]}{$k[$i]} && $net0{"true"}{$k[$l]}{$k[$i]}
	if($_[6] ne $_[7] && $_[5] ne $_[6] && $_[5] ne $_[7]){
		if($_[0] == 1 && $_[1] == 0 && $_[2] == 0 && $_[3] == 0 && $_[4] == 0){#type 1
			return (0,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 0 && $_[1] == 0 && $_[2] == 0 && $_[3] == 0 && $_[4] == 1){#type 2
			return (1,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 1 && $_[1] == 0 && $_[2] == 0 && $_[3] == 1 && $_[4] == 0){#type 3
			return (2,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 0 && $_[1] == 0 && $_[2] == 1 && $_[3] == 0 && $_[4] == 0){#type 4
			return (3,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 1 && $_[1] == 1 && $_[2] == 0 && $_[3] == 0 && $_[4] == 0){#type 5
			return (4,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 1 && $_[1] == 1 && $_[2] == 0 && $_[3] == 1 && $_[4] == 0){#type 6
			return (5,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 0 && $_[1] == 0 && $_[2] == 0 && $_[3] == 1 && $_[4] == 1){#type 7
			return (6,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 1 && $_[1] == 0 && $_[2] == 0 && $_[3] == 1 && $_[4] == 1){#type 8
			return (7,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 0 && $_[1] == 1 && $_[2] == 0 && $_[3] == 0 && $_[4] == 1){#type 9
			return (8,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 1 && $_[1] == 0 && $_[2] == 1 && $_[3] == 1 && $_[4] == 0){#type 10
			return (9,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 1 && $_[1] == 1 && $_[2] == 1 && $_[3] == 0 && $_[4] == 0){#type 11
			return (10,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 1 && $_[1] == 0 && $_[2] == 1 && $_[3] == 1 && $_[4] == 1){#type 12
			return (11,$_[5],$_[6],$_[7]);
		}
		elsif($_[0] == 1 && $_[1] == 1 && $_[2] == 1 && $_[3] == 1 && $_[4] == 1){#type 13
			return (12,$_[5],$_[6],$_[7]);
		}
	}
}

sub loadGold{
	my $di = '';
	open (F, "<$f") or die "can't open $f\n$!\n";
	{ local $/=undef;  $di=<F>; }
	my @d=split /[\r\n]+/, $di;
	close F;
	my $cc = 0;
	my $cc0 = 0;
	my @nodes = ();
	my %tfs = ();
	if($case == 1){
		while(@d){
			if($d[0] =~ /\t/){
				my @t =  split("\t",shift @d);
				$tfs{$t[0]} = 1;
				$tfall{$t[0]} = 1;
				$nodeall{$t[0]} = 1;
				$nodeall{$t[1]} = 1;
				if($t[2] > $threshold){
					$net{"true"}{$t[0]}{$t[1]} = 1;
					$outnet{$t[0]}{$t[1]}[0] = 'P';
					$outnet{$t[0]}{$t[1]}[1] = 'A';
					if(!defined $outnetnodes{$t[0]}){
						$outnetnodes{$t[0]}[0] = 'TF';
					}
					elsif($outnetnodes{$t[0]}[0] ne 'TF'){
						$outnetnodes{$t[0]}[0] = 'TF';
					}
					if(!defined $outnetnodes{$t[1]}){
						$outnetnodes{$t[1]}[0] = 'nTF';
					}
					$outnetnodes{$t[0]}[1] = 'P';
					$outnetnodes{$t[1]}[1] = 'P';
					$resmotNodeTrue{$t[0]} = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,"rec",0,"tot",0};
					$resmotNodeTrue{$t[1]} = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,"rec",0,"tot",0};
					$cc++;
					#~ if(!defined $inDegree{"true"}{$t[1]}){
						#~ @{$inDegree{"true"}{$t[1]}} = ();
					#~ }
					if(!defined $outDegree{"true"}{$t[0]}){
						@{$outDegree{"true"}{$t[0]}} = ();
						#~ push (@tf,$t[0]);
					}
					push (@{$outDegree{"true"}{$t[0]}},$t[1]);
					#~ push (@{$inDegree{"true"}{$t[1]}},$t[0]);
				}
				else{
					$net{"true"}{$t[0]}{$t[1]} = 0;
					$cc0++;
				}
				if(!defined $done{"true"}{$t[0]}){
					$node++;
					$done{"true"}{$t[0]} = $t[0] ;
					$done{"all"}{$t[0]} = '';
				}
				if(!defined $done{"true"}{$t[1]}){
					$node++;
					$done{"true"}{$t[1]} = $t[1] ;
					$done{"all"}{$t[1]} = '';
				}
			}			
			else{
				shift @d;
			}
		}
			#~ for my $i(keys %{$done{"true"}}){
				#~ for my $j (keys %{$done{"true"}}){
					#~ if(!defined $net{"true"}{$i}{$j}){
						#~ $net{"true"}{$i}{$j} = 0;
						#~ $cc0++;
					#~ }
				#~ }
			#~ }
	}
	elsif($case == 2){
		while(@d){
			if($d[0] =~ /\t/){
				my @t =  split("\t",shift @d);
				$tfs{$t[0]} = 1;
				$tfall{$t[0]} = 1;
				$nodeall{$t[0]} = 1;
				$nodeall{$t[1]} = 1;
				if($t[2] > $threshold){
					$net{"true"}{$t[0]}{$t[1]} = 1;
					$outnet{$t[0]}{$t[1]}[0] = 'P';
					$outnet{$t[0]}{$t[1]}[1] = 'A';
					if(!defined $outnetnodes{$t[0]}){
						$outnetnodes{$t[0]}[0] = 'TF';
					}
					elsif($outnetnodes{$t[0]}[0] ne 'TF'){
						$outnetnodes{$t[0]}[0] = 'TF';
					}
					if(!defined $outnetnodes{$t[1]}){
						$outnetnodes{$t[1]}[0] = 'nTF';
					}
					$outnetnodes{$t[0]}[1] = 'P';
					$outnetnodes{$t[1]}[1] = 'P';
					$resmotNodeTrue{$t[0]} = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,"rec",0,"tot",0};
					$resmotNodeTrue{$t[1]} = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,"rec",0,"tot",0};
					$cc++;
					#~ if(!defined $inDegree{"true"}{$t[1]}){
						#~ @{$inDegree{"true"}{$t[1]}} = ();
					#~ }
					if(!defined $outDegree{"true"}{$t[0]}){
						@{$outDegree{"true"}{$t[0]}} = ();
						#~ push (@tf,$t[0]);
					}
					push (@{$outDegree{"true"}{$t[0]}},$t[1]);
					#~ push (@{$inDegree{"true"}{$t[1]}},$t[0]);
				}
				if(!defined $done{"true"}{$t[0]}){
					$node++;
					$done{"true"}{$t[0]} = $t[0] ;
					$done{"all"}{$t[0]} = '';
				}
				if(!defined $done{"true"}{$t[1]}){
					$node++;
					$done{"true"}{$t[1]} = $t[1] ;
					$done{"all"}{$t[1]} = '';
				}
			}
			else{
				shift @d;
			}
		}
		#~ for my $i(keys %{$done{"true"}}){
			#~ for my $j (keys %{$done{"true"}}){
				#~ if(!defined $net{"true"}{$i}{$j}){
					#~ # $net{"true"}{$i}{$j} = 0;
					#~ $cc0++;
				#~ }
			#~ }
		#~ }
	}
	@keysdone = keys %{$done{"all"}};
	$node = @keysdone;
	$cc0 = ($node*$node) - $cc;
	#~ @tf = keys %{$net{"true"}};
	@tf = sort keys %tfs;
	my $tmp = @tf;
	$TNG  = $cc0;
	print "           \t#TFs\tn\t#P\t#N\n";
	print "REFERENCE :\t$tmp\t$node\t$cc\t$cc0\n";
	$outtxttable .= "      \t#TFs\tn\t#P\t#N\n";
	$outtxttable .= "REFERENCE :\t$tmp\t$node\t$cc\t$cc0\n";
	$countgl[0] = $cc;
}



sub loadPred{
	my $di = '';
	open (F, "<$ff") or die "can't open $ff\n$!\n";
	{ local $/=undef;  $di=<F>; }
	my @d=split /[\r\n]+/, $di;
	close F;
	my $cc = 0;
	my $cc0 = 0;
	my @nodes = ();
	my %tfs = ();
	$node = 0;
	if($case == 1){
		while(@d){
			if($d[0] =~ /\t/){
				my @t =  split("\t",shift @d);
				$tfs{$t[0]} = 1;
				$tfall{$t[0]} = 1;
				$nodeall{$t[0]} = 1;
				$nodeall{$t[1]} = 1;
				if(exists $net{"true"}{$t[0]}{$t[1]}){
					if($t[2] > $threshold){
						$net{"pred"}{$t[0]}{$t[1]} = 1;
						$outnet{$t[0]}{$t[1]}[1] = 'P';
						if(!defined $outnetnodes{$t[0]}){
							$outnetnodes{$t[0]}[0] = 'TF';
						}
						elsif($outnetnodes{$t[0]}[0] ne 'TF'){
							$outnetnodes{$t[0]}[0] = 'TF';
						}
						if(!defined $outnetnodes{$t[1]}){
							$outnetnodes{$t[1]}[0] = 'nTF';
						}
						$outnetnodes{$t[0]}[2] = 'P';
						$outnetnodes{$t[1]}[2] = 'P';
						if($outnetnodes{$t[0]}[1] ne 'P'){
							$outnetnodes{$t[0]}[1] = 'A';
						}
						if($outnetnodes{$t[1]}[1] ne 'P'){
							$outnetnodes{$t[1]}[1] = 'A';
						}
						$resmotNodePred{$t[0]} = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,"rec",0,"tot",0};
						$resmotNodePred{$t[1]} = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,"rec",0,"tot",0};
						$cc++;
						#~ if(!defined $inDegree{"pred"}{$t[1]}){
							#~ @{$inDegree{"pred"}{$t[1]}} = ();
						#~ }
						if(!defined $outDegree{"pred"}{$t[0]}){
							@{$outDegree{"pred"}{$t[0]}} = ();
							#~ push (@tf0,$t[0]);
						}
						push (@{$outDegree{"pred"}{$t[0]}},$t[1]);
						#~ push (@{$inDegree{"pred"}{$t[1]}},$t[0]);
					}
					else{
						$net{"pred"}{$t[0]}{$t[1]} = 0;
						$cc0++;
					}				
					if(!defined $done{"pred"}{$t[0]}){
						$node++;
						$done{"pred"}{$t[0]} = $t[0] ;
						$done{"all"}{$t[0]} = '';
					}
					if(!defined $done{"pred"}{$t[1]}){
						$node++;
						$done{"pred"}{$t[1]} = $t[1] ;
						$done{"all"}{$t[1]} = '';
					}
				}
			}
			else{
				shift @d;
			}
		}
		for my $i(keys %{$net{"true"}}){
			for my $j (keys %{$done{"true"}}){
				if(!defined $net{"pred"}{$i}{$j} && exists $net{"true"}{$i}{$j}){
					$net{"pred"}{$i}{$j} = 0;
					$cc0++;
				}
				if($outnetnodes{$j}[1] eq 'P' && $outnetnodes{$j}[2] ne 'P'){
					$outnetnodes{$j}[2] = 'A';
				}
			}
			if($outnetnodes{$i}[1] eq 'P' && $outnetnodes{$i}[2] ne 'P'){
				$outnetnodes{$i}[2] = 'A';
			}
		}
	}
	elsif($case == 2){
		while(@d){
			if($d[0] =~ /\t/){
				my @t =  split("\t",shift @d);
				$tfs{$t[0]} = 1;
				$tfall{$t[0]} = 1;
				$nodeall{$t[0]} = 1;
				$nodeall{$t[1]} = 1;
				if($t[2] > $threshold){
					$net{"pred"}{$t[0]}{$t[1]} = 1;
					$outnet{$t[0]}{$t[1]}[1] = 'P';
					if(!defined $outnetnodes{$t[0]}){
						$outnetnodes{$t[0]}[0] = 'TF';
					}
					elsif($outnetnodes{$t[0]}[0] ne 'TF'){
						$outnetnodes{$t[0]}[0] = 'TF';
					}
					if(!defined $outnetnodes{$t[1]}){
						$outnetnodes{$t[1]}[0] = 'nTF';
					}
					$outnetnodes{$t[0]}[2] = 'P';
					$outnetnodes{$t[1]}[2] = 'P';
					if($outnetnodes{$t[0]}[1] ne 'P'){
						$outnetnodes{$t[0]}[1] = 'A';
					}
					if($outnetnodes{$t[1]}[1] ne 'P'){
						$outnetnodes{$t[1]}[1] = 'A';
					}
					$resmotNodePred{$t[0]} = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,"rec",0,"tot",0};
					$resmotNodePred{$t[1]} = {0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,"rec",0,"tot",0};
					$cc++;
					#~ if(!defined $inDegree{"pred"}{$t[1]}){
						#~ @{$inDegree{"pred"}{$t[1]}} = ();
					#~ }
					if(!defined $outDegree{"pred"}{$t[0]}){
						@{$outDegree{"pred"}{$t[0]}} = ();
						#~ push (@tf0,$t[0]);
					}
					push (@{$outDegree{"pred"}{$t[0]}},$t[1]);
					#~ push (@{$inDegree{"pred"}{$t[1]}},$t[0]);
				}
				if(!defined $done{"pred"}{$t[0]}){
					$node++;
					$done{"pred"}{$t[0]} = $t[0] ;
					$done{"all"}{$t[0]} = '';
				}
				if(!defined $done{"pred"}{$t[1]}){
					$node++;
					$done{"pred"}{$t[1]} = $t[1] ;
					$done{"all"}{$t[1]} = '';
				}
			}
			else{
				shift @d;
			}
		}
		#~ for my $i(keys %{$done{"pred"}}){
			#~ for my $j (keys %{$done{"pred"}}){
				#~ if(!defined $net{"true"}{$i}{$j}){
					#~ # $net{"true"}{$i}{$j} = 0;
				#~ }
				#~ if($net{"pred"}{$i}{$j} == 1 && $outnet{$i}{$j}[0] ne 'P'){
					#~ $outnet{$i}{$j}[0] = 'A';
				#~ }
				#~ if(!defined $net{"pred"}{$i}{$j}){
					#~ $net{"pred"}{$i}{$j} = 0;
					#~ $cc0++;
				#~ }
			#~ }
		#~ }
		#~ for my $i(keys %{$net{"true"}}){
			#~ for my $j (keys %{$done{"true"}}){
				#~ if(!defined $net{"pred"}{$i}{$j} && exists $net{"true"}{$i}{$j}){
					#~ $net{"pred"}{$i}{$j} = 0;
					#~ $cc0++;
				#~ }
				#~ if($outnetnodes{$j}[1] eq 'P' && $outnetnodes{$j}[2] ne 'P'){
					#~ $outnetnodes{$j}[2] = 'A';
				#~ }
			#~ }
			#~ if($outnetnodes{$i}[1] eq 'P' && $outnetnodes{$i}[2] ne 'P'){
				#~ $outnetnodes{$i}[2] = 'A';
			#~ }
		#~ }
		#~ for my $i(keys %{$done{"true"}}){
			#~ if(!defined $done{"pred"}{$i}){
				#~ $node++;
			#~ }
		#~ }
	}
	@keysdone = keys %{$done{"all"}};
	$node = @keysdone;
	for (my $i = 0; $i < @keysdone; $i++){
		$dk{$keysdone[$i]} = $i;
	}
	#~ @tf0 = keys %{$net{"pred"}};
	$cc0 = ($node*$node) - $cc;
	@tf0 = sort keys %tfs;
	my $tmp = @tf0;
	print "INPUT :   \t$tmp\t$node\t$cc\t$cc0\n";
	$outtxttable .= "INPUT :   \t$tmp\t$node\t$cc\t$cc0\n";
	#~ print "@tf0\n";
	$countgl[1] = $cc;
}


sub printNets{
	my $nodesx = '';
	my $nodestable = '';
	my $nodestableD3 = '';
	my $edgesx = '';
	my $edgestable = '';
	my $edgestableD3 = '';
	my $numnodes = (keys %outnetnodes);
	my $xdim = $numnodes;
	my $ydim = 2*$numnodes/3;
	my $xcent = $xdim/2;
	my $ycent = $ydim/2;
	my $headnet ="<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>
<graph id=\"9908\" label=\"$id\" directed=\"1\" cy:documentVersion=\"3.0\" xmlns:dc=\"http://purl.org/dc/elements/1.1/\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" xmlns:cy=\"http://www.cytoscape.org\" xmlns=\"http://www.cs.rpi.edu/XGMML\">
  <att name=\"networkMetadata\">
    <rdf:RDF>
      <rdf:Description rdf:about=\"http://www.cytoscape.org/\">
        <dc:type>Regulatory Interaction</dc:type>
        <dc:description>LoTo comparison of directed networks</dc:description>
        <dc:identifier>N/A</dc:identifier>
        <dc:date>2015</dc:date>
        <dc:title>$id</dc:title>
        <dc:source>http://www.cytoscape.org/</dc:source>
        <dc:format>Cytoscape-XGMML</dc:format>
      </rdf:Description>
    </rdf:RDF>
  </att>
  <att name=\"shared name\" value=\"$id\" type=\"string\"/>
  <att name=\"selected\" value=\"1\" type=\"boolean\"/>
  <att name=\"name\" value=\"$id\" type=\"string\"/>
  <att name=\"__Annotations\" type=\"list\">
  </att>
  <graphics>
    <att name=\"NETWORK_TITLE\" value=\"$id\" type=\"string\"/>
    <att name=\"NETWORK_CENTER_Y_LOCATION\" value=\"-$ycent\" type=\"string\"/>
    <att name=\"NETWORK_HEIGHT\" value=\"$ydim\" type=\"string\"/>
    <att name=\"NETWORK_DEPTH\" value=\"0.0\" type=\"string\"/>
    <att name=\"NETWORK_CENTER_Z_LOCATION\" value=\"0.0\" type=\"string\"/>
    <att name=\"NETWORK_BACKGROUND_PAINT\" value=\"#ffffff\" type=\"string\"/>
    <att name=\"NETWORK_SCALE_FACTOR\" value=\"2.5\" type=\"string\"/>
    <att name=\"NETWORK_NODE_SELECTION\" value=\"true\" type=\"string\"/>
    <att name=\"NETWORK_CENTER_X_LOCATION\" value=\"-$xcent\" type=\"string\"/>
    <att name=\"NETWORK_EDGE_SELECTION\" value=\"true\" type=\"string\"/>
    <att name=\"NETWORK_WIDTH\" value=\"$xdim\" type=\"string\"/>
  </graphics>";


my $nodenet = "
  <node id=\"XXXXX\" label=\"YYYYY\">
    <att name=\"selected\" value=\"0\" type=\"boolean\"/>
    <att name=\"name\" value=\"YYYYY\" type=\"string\"/>
    <att name=\"type\" value=\"PPPPP\" type=\"string\"/>
    <att name=\"graphlet degree in reference\" value=\"MMMMM\" type=\"string\"/>
    <att name=\"RGD\" value=\"OOOOO\" type=\"string\"/>
    <att name=\"F1\" value=\"QQQQQ\" type=\"string\"/>
    <att name=\"graphlet degree in compared\" value=\"NNNNN\" type=\"string\"/>
    <att name=\"class\" value=\"LLLLL\" type=\"string\"/>
    <att name=\"presence in reference\" value=\"JJJJJ\" type=\"string\"/>
    <att name=\"presence in compared\" value=\"KKKKK\" type=\"string\"/>
    <att name=\"color\" value=\"ZZZZZ\" type=\"string\"/>
    <graphics y=\"-$ycent\" fill=\"ZZZZZ\" x=\"-$xcent\" outline=\"#333333\" h=\"35.0\" z=\"0.0\" width=\"3.0\" w=\"35.0\" type=\"ROUND_RECTANGLE\">
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_1\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_4\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_5\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_SELECTED_PAINT\" value=\"#ffff00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_2\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_4\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_9\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_DEPTH\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_BORDER_STROKE\" value=\"SOLID\" type=\"string\"/>
      <att name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\" value=\"true\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_3\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_9\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_7\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_TRANSPARENCY\" value=\"255\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_3\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_9\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_LABEL_FONT_SIZE\" value=\"12\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_4\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_7\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_LABEL_COLOR\" value=\"#000000\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_6\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_6\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_2\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_2\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_7\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_6\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_LABEL_WIDTH\" value=\"200.0\" type=\"string\"/>
      <att name=\"NODE_VISIBLE\" value=\"true\" type=\"string\"/>
      <att name=\"NODE_SELECTED\" value=\"false\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_8\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_POSITION_5\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_1\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_LABEL_POSITION\" value=\"C,C,c,0.00,0.00\" type=\"string\"/>
      <att name=\"NODE_TOOLTIP\" value=\"\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_8\" value=\"org.cytoscape.ding.customgraphics.NullCustomGraphics,0,[ Remove Graphics ],\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_5\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_LABEL\" value=\"YYYYY\" type=\"string\"/>
      <att name=\"NODE_LABEL_TRANSPARENCY\" value=\"255\" type=\"string\"/>
      <att name=\"NODE_LABEL_FONT_FACE\" value=\"Dialog.plain,plain,12\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_1\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_8\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_CUSTOMGRAPHICS_SIZE_3\" value=\"0.0\" type=\"string\"/>
      <att name=\"NODE_BORDER_TRANSPARENCY\" value=\"255\" type=\"string\"/>
    </graphics>
  </node>";
  
my $edgenet = "
  <edge id=\"XXXXX\" label=\"YYYYY\" source=\"JJJJJ\" target=\"KKKKK\" cy:directed=\"1\">
    <att name=\"interaction\" value=\"1\" type=\"string\"/>
    <att name=\"selected\" value=\"0\" type=\"boolean\"/>
    <att name=\"name\" value=\"YYYYY\" type=\"string\"/>
    <att name=\"source\" value=\"PPPPP\" type=\"string\"/>
    <att name=\"target\" value=\"QQQQQ\" type=\"string\"/>
    <att name=\"reference presence\" value=\"RRRRR\" type=\"string\"/>
    <att name=\"compared presence\" value=\"SSSSS\" type=\"string\"/>
    <att name=\"source presence\" value=\"MMMMM\" type=\"string\"/>
    <att name=\"target presence\" value=\"NNNNN\" type=\"string\"/>
    <att name=\"class\" value=\"OOOOO\" type=\"string\"/>
    <att name=\"color\" value=\"ZZZZZ\" type=\"string\"/>
    <graphics fill=\"ZZZZZ\" width=\"2.0\">
      <att name=\"EDGE_SOURCE_ARROW_SHAPE\" value=\"NONE\" type=\"string\"/>
      <att name=\"EDGE_TRANSPARENCY\" value=\"255\" type=\"string\"/>
      <att name=\"EDGE_LABEL_TRANSPARENCY\" value=\"255\" type=\"string\"/>
      <att name=\"EDGE_BEND\" value=\"\" type=\"string\"/>
      <att name=\"EDGE_TARGET_ARROW_SELECTED_PAINT\" value=\"#ffff00\" type=\"string\"/>
      <att name=\"EDGE_LINE_TYPE\" value=\"SOLID\" type=\"string\"/>
      <att name=\"EDGE_TOOLTIP\" value=\"\" type=\"string\"/>
      <att name=\"EDGE_SOURCE_ARROW_SELECTED_PAINT\" value=\"#ffff00\" type=\"string\"/>
      <att name=\"EDGE_LABEL_FONT_FACE\" value=\"Dialog.plain,plain,10\" type=\"string\"/>
      <att name=\"EDGE_STROKE_SELECTED_PAINT\" value=\"#0033ff\" type=\"string\"/>
      <att name=\"EDGE_LABEL_COLOR\" value=\"#000000\" type=\"string\"/>
      <att name=\"EDGE_CURVED\" value=\"true\" type=\"string\"/>
      <att name=\"EDGE_TARGET_ARROW_UNSELECTED_PAINT\" value=\"ZZZZZ\" type=\"string\"/>
      <att name=\"EDGE_LABEL\" value=\"\" type=\"string\"/>
      <att name=\"EDGE_SELECTED\" value=\"false\" type=\"string\"/>
      <att name=\"EDGE_LABEL_FONT_SIZE\" value=\"10\" type=\"string\"/>
      <att name=\"EDGE_SOURCE_ARROW_UNSELECTED_PAINT\" value=\"ZZZZZ\" type=\"string\"/>
      <att name=\"EDGE_TARGET_ARROW_SHAPE\" value=\"DELTA\" type=\"string\"/>
      <att name=\"EDGE_VISIBLE\" value=\"true\" type=\"string\"/>
    </graphics>
  </edge>";
  
	if($ff eq ''){
		$nodestable = "name\tcolor\tgraphlet_degree\ttype\n";
		$edgestable = "name\tsource\ttarget\tinteraction\tcolor\n";
		$nodenet =~ s/
    <att name=\"RGD\" value=\"OOOOO\" type=\"string\"\/>
    <att name=\"F1\" value=\"QQQQQ\" type=\"string\"\/>
    <att name=\"graphlet degree in compared\" value=\"NNNNN\" type=\"string\"\/>
    <att name=\"class\" value=\"LLLLL\" type=\"string\"\/>
    <att name=\"presence in reference\" value=\"JJJJJ\" type=\"string\"\/>
    <att name=\"presence in compared\" value=\"KKKKK\" type=\"string\"\/>//;
		$edgenet =~ s/
    <att name=\"reference presence\" value=\"RRRRR\" type=\"string\"\/>
    <att name=\"compared presence\" value=\"SSSSS\" type=\"string\"\/>    
    <att name=\"source presence\" value=\"MMMMM\" type=\"string\"\/>
    <att name=\"target presence\" value=\"NNNNN\" type=\"string\"\/>
    <att name=\"class\" value=\"OOOOO\" type=\"string\"\/>//;
		my $counter = 1;
		my %numnodes = ();
		for my $i (keys %outnetnodes){
			$numnodes{$i} = $counter;
			my $color = $colorEFTP;
			if($outnetnodes{$i}[0] eq 'TF'){
				$color = $colorTFTP;
			}
			$nodestable .= "$i\t$color\t$outnetnodes{$i}[3]\t$outnetnodes{$i}[0]\n";
			$nodestableD3 .= "\t\t{id: $counter, label: \'$i\', color:\'$color\'},\n";
			my $tmpnode = $nodenet;
			$tmpnode =~ s/XXXXX/$counter/g;
			$tmpnode =~ s/YYYYY/$i/g;
			$tmpnode =~ s/PPPPP/$outnetnodes{$i}[0]/g;
			$tmpnode =~ s/MMMMM/$outnetnodes{$i}[3]/g;
			$tmpnode =~ s/ZZZZZ/$color/g;
			$nodesx .= $tmpnode;
			$counter++;
		}
		for my $i (keys %outnet){
			for my $j (keys %{$outnet{$i}}){
				my $color = $coloredgeTP;
				$edgestable .= "$i(1)$j\t$i\t$j\t1\t$color\n";
				$edgestableD3 .= "\t\t{from: $numnodes{$i}, to: $numnodes{$j}, color:{color:\'$color\'}},\n";
				my $tmpedge = $edgenet;
				$tmpedge =~ s/XXXXX/$counter/g;
				$tmpedge =~ s/YYYYY/$i(1)$j/g;
				$tmpedge =~ s/JJJJJ/$numnodes{$i}/g;
				$tmpedge =~ s/KKKKK/$numnodes{$j}/g;
				$tmpedge =~ s/PPPPP/$i/g;
				$tmpedge =~ s/QQQQQ/$j/g;
				$tmpedge =~ s/ZZZZZ/$color/g;
				$edgesx .= $tmpedge;
				$counter++;
			}
		}
	}
	else{
		$nodestable = "name\tcolor\ttype\tpresence reference\tpresence compared\tclass\tgraphlet_degree reference\tgraphlet_degree compared\tREC\n";
		$edgestable = "name\tsource\ttarget\tinteraction\tcolor\tpresence reference\tpresence compared\tsource presence\ttarget presence\tclass\n";
		my $counter = 1;
		my %numnodes = ();
		for my $i (keys %outnetnodes){
			$numnodes{$i} = $counter;
			my $color = $colorEFTP;
			my $class = "TP";
			if($outnetnodes{$i}[0] eq 'TF'){
				if($outnetnodes{$i}[1] eq 'A' && $outnetnodes{$i}[2] eq 'P'){
					$color = $colorTFFP;
					$class = "FP";
				}
				elsif($outnetnodes{$i}[1] eq 'P' && $outnetnodes{$i}[2] eq 'A'){
					$color = $colorTFFN;
					$class = "FN";
				}
				else{
					$color = $colorTFTP;
				}
			}
			else{
				if($outnetnodes{$i}[1] eq 'A' && $outnetnodes{$i}[2] eq 'P'){
					$color = $colorEFFP;
					$class = "FP";
				}
				elsif($outnetnodes{$i}[1] eq 'P' && $outnetnodes{$i}[2] eq 'A'){
					$color = $colorEFFN;
					$class = "FN";
				}
			}
			if($outnetnodes{$i}[1] eq 'A'){
				$outnetnodes{$i}[3] = 0;
				$outnetnodes{$i}[4] = 0;
			}
			if($outnetnodes{$i}[2] eq 'A'){
				$outnetnodes{$i}[5] = 0;
			}
			$nodestable .= "$i\t$color\t$outnetnodes{$i}[0]\t$outnetnodes{$i}[1]\t$outnetnodes{$i}[2]\t$class\t$outnetnodes{$i}[3]\t$outnetnodes{$i}[5]\t$outnetnodes{$i}[4]\n";
			$nodestableD3 .= "\t\t{id: $counter, label: \'$i\', color:\'$color\'},\n";
			my $tmpnode = $nodenet;
			$tmpnode =~ s/XXXXX/$counter/g;
			$tmpnode =~ s/YYYYY/$i/g;
			$tmpnode =~ s/LLLLL/$class/g;
			$tmpnode =~ s/PPPPP/$outnetnodes{$i}[0]/g;
			$tmpnode =~ s/JJJJJ/$outnetnodes{$i}[1]/g;
			$tmpnode =~ s/KKKKK/$outnetnodes{$i}[2]/g;
			$tmpnode =~ s/MMMMM/$outnetnodes{$i}[3]/g;
			$tmpnode =~ s/OOOOO/$outnetnodes{$i}[4]/g;
			$tmpnode =~ s/QQQQQ/$outnetnodes{$i}[6]/g;
			$tmpnode =~ s/NNNNN/$outnetnodes{$i}[5]/g;
			$tmpnode =~ s/ZZZZZ/$color/g;
			$nodesx .= $tmpnode;
			$counter++;
		}
		for my $i (keys %outnet){
			for my $j (keys %{$outnet{$i}}){
				my $color = $coloredgeTP;
				my $class = "TP";
				if($outnet{$i}{$j}[0] eq 'A'&& $outnet{$i}{$j}[1] eq 'P'){
					$color = $coloredgeFP;
					$class = "FP";
				}
				elsif($outnet{$i}{$j}[0] eq 'P'&& $outnet{$i}{$j}[1] eq 'A'){
					$color = $coloredgeFN;
					$class = "FN";
				}
				$edgestable .= "$i(1)$j\t$i\t$j\t1\t$color\t$outnet{$i}{$j}[0]\t$outnet{$i}{$j}[1]\t$outnetnodes{$i}[1]$outnetnodes{$i}[2]\t$outnetnodes{$j}[1]$outnetnodes{$j}[2]\t$class\n";
				$edgestableD3 .= "\t\t{from: $numnodes{$i}, to: $numnodes{$j}, color:{color:\'$color\'}},\n";
				my $tmpedge = $edgenet;
				$tmpedge =~ s/XXXXX/$counter/g;
				$tmpedge =~ s/YYYYY/$i(1)$j/g;
				$tmpedge =~ s/JJJJJ/$numnodes{$i}/g;
				$tmpedge =~ s/KKKKK/$numnodes{$j}/g;
				$tmpedge =~ s/PPPPP/$i/g;
				$tmpedge =~ s/QQQQQ/$j/g;
				$tmpedge =~ s/OOOOO/$class/g;
				$tmpedge =~ s/RRRRR/$outnet{$i}{$j}[0]/g;
				$tmpedge =~ s/SSSSS/$outnet{$i}{$j}[1]/g;
				$tmpedge =~ s/MMMMM/$outnetnodes{$i}[1]$outnetnodes{$i}[2]/g;
				$tmpedge =~ s/NNNNN/$outnetnodes{$j}[1]$outnetnodes{$j}[2]/g;
				$tmpedge =~ s/ZZZZZ/$color/g;
				$edgesx .= $tmpedge;
				$counter++;
			}
		}
	}
	open (F, ">$folder/$outfile.xgmml") or die "can't open $outfile.xgmml\n";
	print F $headnet;
	print F $nodesx;
	print F $edgesx;
	print F "</graph>";
	close F;
	open (F, ">$folder/$outfile\_nodes.tsv") or die "can't open $outfile\_nodes.tsv\n";
	print F $nodestable;
	close F;
	open (F, ">$folder/$outfile\_edges.tsv") or die "can't open $outfile\_edges.tsv\n";
	print F $edgestable;
	close F;
	open (F, ">$folder/$outfile\_D3.tsv") or die "can't open $outfile\_nodesD3.tsv\n";
	print F "\tvar nodes = new vis.DataSet([\n";
	print F $nodestableD3;
	print F "\t]);\n\tvar edges = new vis.DataSet([\n";
	print F $edgestableD3;
	print F "\t]);";
	close F;
	#~ open (F, ">$folder/$outfile\_nodesD3.tsv") or die "can't open $outfile\_nodesD3.tsv\n";
	#~ print F $nodestableD3;
	#~ close F;
	#~ open (F, ">$folder/$outfile\_edgesD3.tsv") or die "can't open $outfile\_edgesD3.tsv\n";
	#~ print F $edgestableD3;
	#~ close F;
}

sub findMotives_gold{
	my %donegraph = ();
	for my $i(@tf){#tfs in gold
		for my $l (keys %{$done{"true"}}){
			for my $j (@{$outDegree{"true"}{$i}}){
				if(defined$net{"true"}{$i}{$l} && defined$net{"true"}{$j}{$l} && defined$net{"true"}{$l}{$j} && defined$net{"true"}{$j}{$i} && defined$net{"true"}{$l}{$i}){
					my (@tmp) = dotriplet_noij($net{"true"}{$i}{$l}, $net{"true"}{$j}{$l}, $net{"true"}{$l}{$j}, $net{"true"}{$j}{$i}, $net{"true"}{$l}{$i}, $i, $j, $l);
					
					if($tmp[1] && !exists $donegraph{"$tmp[1] $tmp[2] $tmp[3]"} &&  !exists $donegraph{"$tmp[1] $tmp[3] $tmp[2]"} &&  
					!exists $donegraph{"$tmp[2] $tmp[1] $tmp[3]"} &&  !exists $donegraph{"$tmp[2] $tmp[3] $tmp[1]"} &&  
					!exists $donegraph{"$tmp[3] $tmp[1] $tmp[2]"} &&  !exists $donegraph{"$tmp[3] $tmp[2] $tmp[1]"}){
						$donegraph{"$tmp[1] $tmp[2] $tmp[3]"} = $tmp[0];
						$resmotNodeTrue{$tmp[1]}{$tmp[0]}++;
						$resmotNodeTrue{$tmp[2]}{$tmp[0]}++;
						$resmotNodeTrue{$tmp[3]}{$tmp[0]}++;
						$resmotNodeTrue{$tmp[1]}{"tot"}++;
						$resmotNodeTrue{$tmp[2]}{"tot"}++;
						$resmotNodeTrue{$tmp[3]}{"tot"}++;
						$res{"true"}{$tmp[1]}{$tmp[0]}{"$tmp[2] $tmp[3]"} = 1;
						$countt{$tmp[0]}++;
					}
				}
			}
		}
	}
	 %donegraph = ();
}


sub findMotives_edgeres{	
	my $count = 0;
	my %donegraphpred = ();
	print "count: $count\n" . `date` . "\n";
	for my $i(@tf0){#tfs in ture
		for my $l (keys %{$done{"pred"}}){
			for my $j (@{$outDegree{"pred"}{$i}}){
				my @t = (0,0,0,0,0);
				if(defined $net{"pred"}{$i}{$l}){
					$t[0] = $net{"pred"}{$i}{$l};
				}
				if(defined $net{"pred"}{$j}{$l}){
					$t[1] = $net{"pred"}{$j}{$l};
				}
				if(defined $net{"pred"}{$l}{$j}){
					$t[2] = $net{"pred"}{$l}{$j};
				}
				if(defined $net{"pred"}{$j}{$i}){
					$t[3] = $net{"pred"}{$j}{$i};
				}
				if(defined $net{"pred"}{$l}{$i}){
					$t[4] = $net{"pred"}{$l}{$i};
				}
				my (@tmp) = dotriplet_noij($t[0], $t[1], $t[2], $t[3], $t[4], $i, $j, $l);
				
				if($tmp[1] && !exists $donegraphpred{"$dk{$tmp[1]} $dk{$tmp[2]} $dk{$tmp[3]}"} &&  !exists $donegraphpred{"$dk{$tmp[1]} $dk{$tmp[3]} $dk{$tmp[2]}"} &&  
				!exists $donegraphpred{"$dk{$tmp[2]} $dk{$tmp[1]} $dk{$tmp[3]}"} &&  !exists $donegraphpred{"$dk{$tmp[2]} $dk{$tmp[3]} $dk{$tmp[1]}"} &&  
				!exists $donegraphpred{"$dk{$tmp[3]} $dk{$tmp[1]} $dk{$tmp[2]}"} &&  !exists $donegraphpred{"$dk{$tmp[3]} $dk{$tmp[2]} $dk{$tmp[1]}"}){
					$donegraphpred{"$dk{$tmp[1]} $dk{$tmp[2]} $dk{$tmp[3]}"} = $tmp[0];
					$resmotNodePred{$tmp[1]}{$tmp[0]}++;
					$resmotNodePred{$tmp[2]}{$tmp[0]}++;
					$resmotNodePred{$tmp[3]}{$tmp[0]}++;
					$resmotNodePred{$tmp[1]}{"tot"}++;
					$resmotNodePred{$tmp[2]}{"tot"}++;
					$resmotNodePred{$tmp[3]}{"tot"}++;
					$countp{$tmp[0]}++;
					#~ $res{"pred"}{$tmp[1]}{$tmp[0]}{"$tmp[2] $tmp[3]"} = 1;
				}
				#~ $count++;
			}
			if(defined $net{"true"}{$i}{$l}){
				if(defined $net{"pred"}{$i}{$l}){
					$res{"global"}[0]++;#tp
				}
				else{
					$res{"global"}[3]++;#fn
				}
			}
			elsif(defined $net{"pred"}{$i}{$l}){
				$res{"global"}[1]++;#fp
			}
			else{
				$res{"global"}[2]++;#tn
			}
		}
		$count++;
		print " $count";
	}
	print "\ncount: $count\n" . `date` . "\n";
	$count = 0;
	#~ %donegraphpred = ();
	my %donegraph = ();
	for my $i(@tf){#tfs in gold
		for my $l (keys %{$done{"true"}}){
			for my $j (@{$outDegree{"true"}{$i}}){
				
				my @t = (0,0,0,0,0);
				if(defined $net{"true"}{$i}{$l}){
					$t[0] = $net{"true"}{$i}{$l};
				}
				if(defined $net{"true"}{$j}{$l}){
					$t[1] = $net{"true"}{$j}{$l};
				}
				if(defined $net{"true"}{$l}{$j}){
					$t[2] = $net{"true"}{$l}{$j};
				}
				if(defined $net{"true"}{$j}{$i}){
					$t[3] = $net{"true"}{$j}{$i};
				}
				if(defined $net{"true"}{$l}{$i}){
					$t[4] = $net{"true"}{$l}{$i};
				}
				my (@tmp) = dotriplet_noij($t[0], $t[1], $t[2], $t[3], $t[4], $i, $j, $l);
				
				if($tmp[1] && !exists $donegraph{"$dk{$tmp[1]} $dk{$tmp[2]} $dk{$tmp[3]}"} &&  !exists $donegraph{"$dk{$tmp[1]} $dk{$tmp[3]} $dk{$tmp[2]}"} &&  
				!exists $donegraph{"$dk{$tmp[2]} $dk{$tmp[1]} $dk{$tmp[3]}"} &&  !exists $donegraph{"$dk{$tmp[2]} $dk{$tmp[3]} $dk{$tmp[1]}"} &&  
				!exists $donegraph{"$dk{$tmp[3]} $dk{$tmp[1]} $dk{$tmp[2]}"} &&  !exists $donegraph{"$dk{$tmp[3]} $dk{$tmp[2]} $dk{$tmp[1]}"}){
					$donegraph{"$dk{$tmp[1]} $dk{$tmp[2]} $dk{$tmp[3]}"} = $tmp[0];
					$resmotNodeTrue{$tmp[1]}{$tmp[0]}++;
					$resmotNodeTrue{$tmp[2]}{$tmp[0]}++;
					$resmotNodeTrue{$tmp[3]}{$tmp[0]}++;
					$resmotNodeTrue{$tmp[1]}{"tot"}++;
					$resmotNodeTrue{$tmp[2]}{"tot"}++;
					$resmotNodeTrue{$tmp[3]}{"tot"}++;
					$countt{$tmp[0]}++;
					#~ $res{"true"}{$tmp[1]}{$tmp[0]}{"$tmp[2] $tmp[3]"} = 1;
					if((defined $donegraphpred{"$dk{$tmp[1]} $dk{$tmp[2]} $dk{$tmp[3]}"} && $donegraphpred{"$dk{$tmp[1]} $dk{$tmp[2]} $dk{$tmp[3]}"} == $tmp[0])||
					(defined $donegraphpred{"$dk{$tmp[1]} $dk{$tmp[3]} $dk{$tmp[2]}"} && $donegraphpred{"$dk{$tmp[1]} $dk{$tmp[3]} $dk{$tmp[2]}"} == $tmp[0])||
					(defined $donegraphpred{"$dk{$tmp[2]} $dk{$tmp[3]} $dk{$tmp[1]}"} && $donegraphpred{"$dk{$tmp[2]} $dk{$tmp[3]} $dk{$tmp[1]}"} == $tmp[0])||
					(defined $donegraphpred{"$dk{$tmp[2]} $dk{$tmp[1]} $dk{$tmp[3]}"} && $donegraphpred{"$dk{$tmp[2]} $dk{$tmp[1]} $dk{$tmp[3]}"} == $tmp[0])||
					(defined $donegraphpred{"$dk{$tmp[3]} $dk{$tmp[2]} $dk{$tmp[1]}"} && $donegraphpred{"$dk{$tmp[3]} $dk{$tmp[2]} $dk{$tmp[1]}"} == $tmp[0])||
					(defined $donegraphpred{"$dk{$tmp[3]} $dk{$tmp[1]} $dk{$tmp[2]}"} && $donegraphpred{"$dk{$tmp[3]} $dk{$tmp[1]} $dk{$tmp[2]}"} == $tmp[0])){
						$resmot{$tmp[0]}{"TP"}++;
						$resmot{$tmp[0]}{"TN"}--;
						$resmot{$tmp[0]}{"rec"} += 6;
						$resmotNodeTrue{$tmp[1]}{"rec"}+=6;
						$resmotNodeTrue{$tmp[2]}{"rec"}+=6;
						$resmotNodeTrue{$tmp[3]}{"rec"}+=6;
					}
					else{
						$resmot{$tmp[0]}{"FN"}++;
						$resmot{$tmp[0]}{"TN"}--;
						
						if((defined $net{"pred"}{$i}{$j} && defined $net{"true"}{$i}{$j}) || (!defined $net{"pred"}{$i}{$j} && !defined $net{"true"}{$i}{$j})){
							$resmot{$tmp[0]}{"rec"}++;
							$resmotNodeTrue{$tmp[1]}{"rec"}++;
							$resmotNodeTrue{$tmp[2]}{"rec"}++;
							$resmotNodeTrue{$tmp[3]}{"rec"}++;
						}
						if((defined $net{"pred"}{$i}{$l} && defined $net{"true"}{$i}{$l}) || (!defined $net{"pred"}{$i}{$l} && !defined $net{"true"}{$i}{$l})){
							$resmot{$tmp[0]}{"rec"}++;
							$resmotNodeTrue{$tmp[1]}{"rec"}++;
							$resmotNodeTrue{$tmp[2]}{"rec"}++;
							$resmotNodeTrue{$tmp[3]}{"rec"}++;
						}
						if((defined $net{"pred"}{$j}{$l} && defined $net{"true"}{$j}{$l}) || (!defined $net{"pred"}{$j}{$l} && !defined $net{"true"}{$j}{$l})){
							$resmot{$tmp[0]}{"rec"}++;
							$resmotNodeTrue{$tmp[1]}{"rec"}++;
							$resmotNodeTrue{$tmp[2]}{"rec"}++;
							$resmotNodeTrue{$tmp[3]}{"rec"}++;
						}
						if((defined $net{"pred"}{$l}{$j} && defined $net{"true"}{$l}{$j}) || (!defined $net{"pred"}{$l}{$j} && !defined $net{"true"}{$l}{$j})){
							$resmot{$tmp[0]}{"rec"}++;
							$resmotNodeTrue{$tmp[1]}{"rec"}++;
							$resmotNodeTrue{$tmp[2]}{"rec"}++;
							$resmotNodeTrue{$tmp[3]}{"rec"}++;
						}
						if((defined $net{"pred"}{$j}{$i} && defined $net{"true"}{$j}{$i}) || (!defined $net{"pred"}{$j}{$i} && !defined $net{"true"}{$j}{$i})){
							$resmot{$tmp[0]}{"rec"}++;
							$resmotNodeTrue{$tmp[1]}{"rec"}++;
							$resmotNodeTrue{$tmp[2]}{"rec"}++;
							$resmotNodeTrue{$tmp[3]}{"rec"}++;
						}
						if((defined $net{"pred"}{$l}{$i} && defined $net{"true"}{$l}{$i}) || (!defined $net{"pred"}{$l}{$i} && !defined $net{"true"}{$l}{$i})){
							$resmot{$tmp[0]}{"rec"}++;
							$resmotNodeTrue{$tmp[1]}{"rec"}++;
							$resmotNodeTrue{$tmp[2]}{"rec"}++;
							$resmotNodeTrue{$tmp[3]}{"rec"}++;
						}						
					}
					#~ delete $res{"pred"}{$tmp[1]}{$tmp[0]}{"$tmp[2] $tmp[3]"};
					#~ delete $res{"pred"}{$tmp[1]}{$tmp[0]}{"$tmp[3] $tmp[2]"};
					#~ delete $res{"pred"}{$tmp[2]}{$tmp[0]}{"$tmp[1] $tmp[3]"};
					#~ delete $res{"pred"}{$tmp[2]}{$tmp[0]}{"$tmp[3] $tmp[1]"};
					#~ delete $res{"pred"}{$tmp[3]}{$tmp[0]}{"$tmp[1] $tmp[2]"};
					#~ delete $res{"pred"}{$tmp[3]}{$tmp[0]}{"$tmp[2] $tmp[1]"};
					delete $donegraphpred{"$dk{$tmp[1]} $dk{$tmp[2]} $dk{$tmp[3]}"};
					delete $donegraphpred{"$dk{$tmp[1]} $dk{$tmp[3]} $dk{$tmp[2]}"};
					delete $donegraphpred{"$dk{$tmp[2]} $dk{$tmp[1]} $dk{$tmp[3]}"};
					delete $donegraphpred{"$dk{$tmp[2]} $dk{$tmp[3]} $dk{$tmp[1]}"};
					delete $donegraphpred{"$dk{$tmp[3]} $dk{$tmp[1]} $dk{$tmp[2]}"};
					delete $donegraphpred{"$dk{$tmp[3]} $dk{$tmp[2]} $dk{$tmp[1]}"};
				}
				#~ $count++;
			}
		}
		$count++;
		print " $count";
	}
	
	#~ for(my $j = 0; $j < 13; $j++){
		#~ for my $i (@tf0){
			#~ for my $l (keys %{$res{"pred"}{$i}{$j}}){
				#~ $resmot{$j}{"FP"}++;
				#~ delete $res{"pred"}{$i}{$j}{$l};
			#~ }
		#~ }
	#~ }
	for my $i (keys %donegraphpred){
		#~ for my $l (keys %{$res{"pred"}{$i}{$j}}){
		$resmot{$donegraphpred{$i}}{"FP"}++;
		delete $donegraphpred{$i};
		#~ }
	}
	print "\ncount: $count\n" . `date` . "\n";
	#~ $outtxt .= "count: $count\n" . `date` . "\n";
}

sub printres_gold{
	
	my $tmpoutall = '';
	my $tmpout = '';
	my $tmpNOG = "TFs not in graphlets:\n";
	my $tmpoutn = '';
	my $tmpNOGn = "Genes (no TF) not in graphlets:\n";
	
	print "\nGRAPHLETS:\n";
	print "Type:";
	$tmpoutall  .= "\nGRAPHLETS:\n";
	$tmpoutall  .= "Type:";
	for my $j (0..12){
		print "\t";
		print $j+1;
		$tmpoutall  .= "\t";
		$tmpoutall  .= $j+1;
	}
	print "\ttotal\ncount";
	$tmpoutall  .= "\ttotal\ncount";
	my $i = 0;
	for my $j (0..12){
		print "\t";
		print $countt{$j};
		$tmpoutall  .= "\t";
		$tmpoutall  .= $countt{$j};
		$i += $countt{$j};
	}
	print "\t$i\n\n";
	$tmpoutall  .= "\t$i";
	
	$tmpout  .= "TFs in graphlets:\nTF";
	for my $j (0..12){
		$tmpout  .= "\t";
		$tmpout  .= $j+1;
	}
	$tmpout  .= "\ttotal\n";
	for my $i (@tf){
		
		if(0 != $resmotNodeTrue{$i}{"tot"}){
			$tmpout  .= "$i"; 
			for my $j (0..12){
				$tmpout  .= "\t$resmotNodeTrue{$i}{$j}";
			}
			$tmpout  .= "\t$resmotNodeTrue{$i}{tot}\n";
			$outnetnodes{$i}[3] =  $resmotNodeTrue{$i}{"tot"};
		}
		else{
			$tmpNOG .= "$i\n";
			$outnetnodes{$i}[3] = 0;
		}
		delete $resmotNodeTrue{$i};
	}
	
	$tmpoutn.= "Genes (no TF) in graphlets:\ngene";
	for my $j (0..5){
		$tmpoutn.= "\t";
		$tmpoutn.= $j+1;
	}
	$tmpoutn.= "\ttotal\n";
	for my $i (sort keys %resmotNodeTrue){
		if(0 != $resmotNodeTrue{$i}{"tot"}){
			$tmpoutn.= "$i"; 
			for my $j (0..5){
				$tmpoutn.= "\t$resmotNodeTrue{$i}{$j}";
			}
			$tmpoutn.= "\t$resmotNodeTrue{$i}{tot}\n";
			$outnetnodes{$i}[3] = $resmotNodeTrue{$i}{"tot"};
		}
		else{
			$tmpNOGn .= "$i\n";
			$outnetnodes{$i}[3] = 0;
		}
	}
	
	
	open(F, ">$folder/$outfile") or die "can't open $folder/$outfile";#$countlist{$tmp[0]}{"TP"} .= "$tmp[1] $tmp[2] $tmp[3]\n";
	print F "$tmpoutall\n\n";
	print F "$tmpNOG\n\n";
	print F "$tmpNOGn\n\n";
	print F "$tmpout\n\n";
	print F "$tmpoutn\n";
	
	print F "GRAPHLETS:\n\n";
	for my $i (0..12){
		print F "Type ";
		print F $i+1;
		print F ":\n";
		for my $j (@tf){
			for my $k (keys %{$res{"true"}{$j}{$i}}){
				print  F "$j $k\n";
			}
		}
		print F "\n";
	}
	close F;
	
}

sub printres{

	my $tfn = (keys %tfall);
	$node =  (keys %nodeall);

		$resmot{"0"}{"TN"} = - ($resmot{"0"}{"TP"}+$resmot{"0"}{"FP"}+$resmot{"0"}{"FN"}) + ($tfn * ($node - 1) * ($node - 2) * 0.5);#total max of posssible motives - #TP = #TN, with self loops, ABC not the same as ACB
		$resmot{"1"}{"TN"} = - ($resmot{"1"}{"TP"}+$resmot{"1"}{"FP"}+$resmot{"1"}{"FN"}) +  ($tfn * ($tfn - 1) * ($node - 2));
		$resmot{"2"}{"TN"} = - ($resmot{"2"}{"TP"}+$resmot{"2"}{"FP"}+$resmot{"2"}{"FN"}) +  ($tfn * ($tfn - 1) * ($node - 2));
		$resmot{"3"}{"TN"} = - ($resmot{"3"}{"TP"}+$resmot{"3"}{"FP"}+$resmot{"3"}{"FN"}) +  ($tfn * ($tfn - 1) * ($node - 2) * 0.5);
		$resmot{"4"}{"TN"} = - ($resmot{"4"}{"TP"}+$resmot{"4"}{"FP"}+$resmot{"4"}{"FN"}) +  ($tfn * ($tfn - 1) * ($node - 2));
		$resmot{"5"}{"TN"} = - ($resmot{"5"}{"TP"}+$resmot{"5"}{"FP"}+$resmot{"5"}{"FN"}) +  ($tfn * ($tfn - 1) * ($node - 2) * 0.5);
		$resmot{"6"}{"TN"} = - ($resmot{"6"}{"TP"}+$resmot{"6"}{"FP"}+$resmot{"6"}{"FN"}) +  ($tfn * ($tfn - 1) * ($tfn - 2));
		$resmot{"7"}{"TN"} = - ($resmot{"7"}{"TP"}+$resmot{"7"}{"FP"}+$resmot{"7"}{"FN"}) +  ($tfn * ($tfn-1) * ($tfn - 2) * 0.5);
		$resmot{"8"}{"TN"} = - ($resmot{"8"}{"TP"}+$resmot{"8"}{"FP"}+$resmot{"8"}{"FN"}) +  ($tfn * ($tfn-1) * ($tfn - 2) * (1/3));
		$resmot{"9"}{"TN"} = - ($resmot{"9"}{"TP"}+$resmot{"9"}{"FP"}+$resmot{"9"}{"FN"}) +  ($tfn * ($tfn-1) * ($tfn - 2));
		$resmot{"10"}{"TN"} = - ($resmot{"10"}{"TP"}+$resmot{"10"}{"FP"}+$resmot{"10"}{"FN"}) +  ($tfn * ($tfn-1) * ($tfn - 2) * 0.5);
		$resmot{"11"}{"TN"} = - ($resmot{"11"}{"TP"}+$resmot{"11"}{"FP"}+$resmot{"11"}{"FN"}) +  ($tfn * ($tfn-1) * ($tfn - 2));
		$resmot{"12"}{"TN"} = - ($resmot{"12"}{"TP"}+$resmot{"12"}{"FP"}+$resmot{"12"}{"FN"}) +  ($tfn * ($tfn-1) * ($tfn - 2) * (1/6));
	my $ctt = 0;
	my $ctp = 0;
	my $mt = 0;
	for (my $j = 0; $j < 13; $j++){
		$ctt += $countt{$j};
		$ctp += $countp{$j};
		if($countt{$j}||$countp{$j}){$mt++;}
	}
	
	print "\n\nG\t#T\t#P\tTP\tFP\tTN\tFN\tR\tP\tFPR\tACC\tF1\tMCC\tREC\n";
	$outtxttable .= "\n\nG\t#T\t#P\tTP\tFP\tTN\tFN\tR\tSPC\tP\tNPV\tFPR\tFDR\tFNR\tACC\tF1\tMCC\tinfor\tmarked\tJ\tSD\tK1\tK2\tO\tCR\tH\tREC\n";
	for (my $j = 0; $j < 13; $j++){
		$res{"gm"}[0] += $resmot{$j}{"TP"};
		$res{"gm"}[1] += $resmot{$j}{"FP"};
		$res{"gm"}[2] += $resmot{$j}{"TN"};
		$res{"gm"}[3] += $resmot{$j}{"FN"};
		print $j+1 . "\t$countt{$j}\t$countp{$j}\t";
		$outtxttable .= $j+1 . "\t$countt{$j}\t$countp{$j}\t";
		my @tmp = compute_stats($resmot{$j}{"TP"},$resmot{$j}{"FP"},$resmot{$j}{"TN"},$resmot{$j}{"FN"});
		for(my $i = 0; $i < 23; $i++){
			if($i > 3 && $tmp[$i] != 0 && abs($tmp[$i]) != 1 ){
				$outtxttable .= sprintf("%.4f\t",$tmp[$i]);
			}
			else{
				$outtxttable .= "$tmp[$i]\t";
			}
		}
		for my $i(0,1,2,3,4,6,8,11,12,13){
			if($i > 3 && $tmp[$i] != 0 && abs($tmp[$i]) != 1 ){
				printf("%.4f\t",$tmp[$i]);
			}
			else{
				print "$tmp[$i]\t";
			}
		}
		if(defined $countt{$j} && $countt{$j} > 0){
			$res{"gm"}[4] += (($resmot{$j}{"rec"}/6)/$countt{$j});
			printf("%.4f",(($resmot{$j}{"rec"}/6)/$countt{$j}));
			$outtxttable .= sprintf("%.4f",(($resmot{$j}{"rec"}/6)/$countt{$j}));
		}
		else{
			print "0";
			$outtxttable .=  "0";
		}
		print "\n";
		$outtxttable .=  "\n";
	}
	my @tmp = compute_stats($res{"gm"}[0],$res{"gm"}[1],$res{"gm"}[2],$res{"gm"}[3]);
	print "all\t$ctt\t$ctp\t";
	$outtxttable .= "all\t$ctt\t$ctp\t";
	for(my $i = 0; $i < 23; $i++){
		if($i > 3 && $tmp[$i] != 0 && abs($tmp[$i]) != 1 ){
			#~ printf("%.4f\t",$tmp[$i]);
			$outtxttable .= sprintf("%.4f\t",$tmp[$i]);
		}
		else{
			#~ print "$tmp[$i]\t";
			$outtxttable .= "$tmp[$i]\t";
		}
	} 
	#~ for(my $i = 0; $i < 23; $i++){
	for my $i(0,1,2,3,4,6,8,11,12,13){
		if($i > 3 && $tmp[$i] != 0 && abs($tmp[$i]) != 1 ){
			printf("%.4f\t",$tmp[$i]);
			#~ $outtxt .= sprintf("%.4f\t",$tmp[$i]);
		}
		else{
			print "$tmp[$i]\t";
			#~ $outtxt .= "$tmp[$i]\t";
		}
	} 
	#~ print "\t";
	#~ $outtxt .= "\t";
	if($res{"gm"}[4]){
		printf("%.4f",($res{"gm"}[4]/$mt));
		$outtxttable .= sprintf("%.4f",($res{"gm"}[4]/$mt));
	}
	else{
		print "0";
		$outtxttable .= "0";
	}
	print "\n";
	$outtxttable .= "\n";
	print "gl\t$countgl[0]\t$countgl[1]\t" ;
	$outtxttable .= "gl\t$countgl[0]\t$countgl[1]\t" ;
	#~ @tmp = compute_stats($res{"global"}[0],$res{"global"}[1],$res{"global"}[2],$res{"global"}[3]);
	@tmp = compute_stats($res{"global"}[0],$res{"global"}[1],$TNG,$res{"global"}[3]);
	for(my $i = 0; $i < 23; $i++){
		if($i > 3 && $tmp[$i] != 0 && abs($tmp[$i]) != 1 ){
			#~ printf("%.4f\t",$tmp[$i]);
			$outtxttable .= sprintf("%.4f\t",$tmp[$i]);
		}
		else{
			#~ print "$tmp[$i]\t";
			$outtxttable .= "$tmp[$i]\t";
		}
	}
	#~ for(my $i = 0; $i < 23; $i++){
	for my $i(0,1,2,3,4,6,8,11,12,13){
		if($i > 3 && $tmp[$i] != 0 && abs($tmp[$i]) != 1 ){
			printf("%.4f\t",$tmp[$i]);
			#~ $outtxt .= sprintf("%.4f\t",$tmp[$i]);
		}
		else{
			print "$tmp[$i]\t";
			#~ $outtxt .= "$tmp[$i]\t";
		}
	}
	print "nan\n\n";
	$outtxttable .= "nan\n\n\n";
	#~ print "$ctt $ctp $mt\n";
	my %gdd = ();
	my %gbin = ();
	#~ for my $i (@tf){
	for my $i (@tf){
		@{$gdd{"gold"}{"TF"}{$i}} = (0,0,0,0,0,0,0,0,0,0,0,0,0);
		@{$gbin{"gold"}{"TF"}{$i}} = (0,0,0,0);
	}
	for my $i (@tf0){
	#~ for my $i (@tf){
		@{$gdd{"inf"}{"TF"}{$i}} = (0,0,0,0,0,0,0,0,0,0,0,0,0);
		@{$gbin{"inf"}{"TF"}{$i}} = (0,0,0,0);
	}
	for my $type ("TP"){
		$outtxt .= "\n$type (graphlets in both networks):\n";
		for my $j (0..12){
			$outtxt .=  $j+1 . "\n$countlist{$j}{$type}\n";
			my @tmp = split("\n",$countlist{$j}{$type});
			while(@tmp){
				my @t = split(" ", shift @tmp);
				#~ print "@t||";
				for my $k (0,1,2){
					#~ print "$t[$k]!!";
					if(exists $gdd{"gold"}{"TF"}{$t[$k]}){
						$gdd{"gold"}{"TF"}{$t[$k]}[$j]++;
						$gbin{"gold"}{"TF"}{$t[$k]}[0]++;
						#~ print "$t[$k]!kk!";
					}
					else{
						if(!exists $gdd{"gold"}{"noTF"}{$t[$k]}){
							@{$gdd{"gold"}{"noTF"}{$t[$k]}} = (0,0,0,0,0,0,0,0,0,0,0,0,0);
							@{$gbin{"gold"}{"noTF"}{$t[$k]}} = (0,0,0,0);
						}
						$gdd{"gold"}{"noTF"}{$t[$k]}[$j]++;
						$gbin{"gold"}{"noTF"}{$t[$k]}[0]++;
						#~ print "$t[$k]!\\!\n";
					}
					if(exists $gdd{"inf"}{"TF"}{$t[$k]}){
						$gdd{"inf"}{"TF"}{$t[$k]}[$j]++;
						$gbin{"inf"}{"TF"}{$t[$k]}[0]++;
					}
					else{
						if(!exists $gdd{"inf"}{"noTF"}{$t[$k]}){
							@{$gdd{"inf"}{"noTF"}{$t[$k]}} = (0,0,0,0,0,0,0,0,0,0,0,0,0);
							@{$gbin{"inf"}{"noTF"}{$t[$k]}} = (0,0,0,0);
						}
						$gdd{"inf"}{"noTF"}{$t[$k]}[$j]++;
						$gbin{"inf"}{"noTF"}{$t[$k]}[0]++;
						#~ print "$t[$k]!!\n";
					}
				}
				#~ print "\n";
			}
			#~ print "\n";
		}
	}
	for my $type ("FN"){
		$outtxt .= "\n$type (graphlets only in REFERENCE):\n";
		for my $j (0..12){
			$outtxt .=  $j+1 . "\n$countlist{$j}{$type}\n";
			my @tmp = split("\n",$countlist{$j}{$type});
			while(@tmp){
				my @t = split(" ", shift @tmp);
				#~ print "@t||";
				for my $k (0,1,2){
					#~ print "$t[$k]!!";
					if(exists $gdd{"gold"}{"TF"}{$t[$k]}){
						$gdd{"gold"}{"TF"}{$t[$k]}[$j]++;
						$gbin{"gold"}{"TF"}{$t[$k]}[3]++;
					}
					else{
						if(!exists $gdd{"gold"}{"noTF"}{$t[$k]}){
							@{$gdd{"gold"}{"noTF"}{$t[$k]}} = (0,0,0,0,0,0,0,0,0,0,0,0,0);
							@{$gbin{"gold"}{"noTF"}{$t[$k]}} = (0,0,0,0);
						}
						$gdd{"gold"}{"noTF"}{$t[$k]}[$j]++;
						$gbin{"gold"}{"noTF"}{$t[$k]}[3]++;
					}
					if(exists $gdd{"inf"}{"TF"}{$t[$k]}){
						#~ $gdd{"inf"}{"TF"}{$t[$k]}[$j]++;
						$gbin{"inf"}{"TF"}{$t[$k]}[3]++;
					}
					else{
						if(!exists $gbin{"inf"}{"noTF"}{$t[$k]}){
							#~ @{$gdd{"inf"}{"noTF"}{$t[$k]}} = (0,0,0,0,0,0,0,0,0,0,0,0,0);
							@{$gbin{"inf"}{"noTF"}{$t[$k]}} = (0,0,0,0);
						}
						#~ $gdd{"inf"}{"noTF"}{$t[$k]}[$j]++;
						$gbin{"inf"}{"noTF"}{$t[$k]}[3]++;
						#~ print "$t[$k]!!\n";
					}
				}
			}
		}
	}
	for my $type ("FP"){
		$outtxt .= "\n$type (graphlets only in INPUT):\n";
		for my $j (0..12){
			$outtxt .=  $j+1 . "\n$countlist{$j}{$type}\n";
			my @tmp = split("\n",$countlist{$j}{$type});
			while(@tmp){
				my @t = split(" ", shift @tmp);
				#~ print "@t||";
				for my $k (0,1,2){
					#~ print "$t[$k]!!";
					if(exists $gdd{"inf"}{"TF"}{$t[$k]}){
						$gdd{"inf"}{"TF"}{$t[$k]}[$j]++;
						$gbin{"inf"}{"TF"}{$t[$k]}[1]++;
					}
					else{
						if(!exists $gdd{"inf"}{"noTF"}{$t[$k]}){
							@{$gdd{"inf"}{"noTF"}{$t[$k]}} = (0,0,0,0,0,0,0,0,0,0,0,0,0);
							@{$gbin{"inf"}{"noTF"}{$t[$k]}} = (0,0,0,0);
						}
						$gdd{"inf"}{"noTF"}{$t[$k]}[$j]++;
						$gbin{"inf"}{"noTF"}{$t[$k]}[1]++;
					}
					if(exists $gbin{"gold"}{"TF"}{$t[$k]}){
						#~ $gdd{"gold"}{"TF"}{$t[$k]}[$j]++;
						$gbin{"gold"}{"TF"}{$t[$k]}[1]++;
						#~ print "$t[$k]!kk!";
					}
					else{
						if(!exists $gbin{"gold"}{"noTF"}{$t[$k]}){
							#~ @{$gdd{"gold"}{"noTF"}{$t[$k]}} = (0,0,0,0,0,0,0,0,0,0,0,0,0);
							@{$gbin{"gold"}{"noTF"}{$t[$k]}} = (0,0,0,0);
						}
						#~ $gdd{"gold"}{"noTF"}{$t[$k]}[$j]++;
						$gbin{"gold"}{"noTF"}{$t[$k]}[1]++;
						#~ print "$t[$k]!\\!\n";
					}
				}
			}
		}	
	}
	
	
	my $tmpout = '';
	my $tmpNOG = "$fg TFs not in graphlets:\n";
	my $tmpoutn = '';
	my $tmpNOGn = "$fg Genes (no TF) not in graphlets:\n";
		
	$tmpout  .= "$fg TFs in graphlets:\nTF";
	for my $j (0..12){
		$tmpout  .= "\t";
		$tmpout  .= $j+1;
	}
	$tmpout  .= "\ttotal\tRGD\tF1\n";
	for my $i (@tf){
		
		if(0 != $resmotNodeTrue{$i}{"tot"}){
			$tmpout  .= "$i"; 
			for my $j (0..12){
				$tmpout  .= "\t$resmotNodeTrue{$i}{$j}";
			}
			$tmpout  .= "\t$resmotNodeTrue{$i}{tot}\t";
			$outnetnodes{$i}[3] = $resmotNodeTrue{$i}{"tot"};
			$outnetnodes{$i}[4] = sprintf("%.4f",($resmotNodeTrue{$i}{"rec"}/($resmotNodeTrue{$i}{"tot"}*6)));
			$tmpout  .= sprintf("%.4f\t",($resmotNodeTrue{$i}{"rec"}/($resmotNodeTrue{$i}{"tot"}*6)));
			if ($gbin{"gold"}{"TF"}{$i}[0] != 0){
				$tmpout  .= sprintf("%.4f",((2*($gbin{"gold"}{"TF"}{$i}[0]/($gbin{"gold"}{"TF"}{$i}[0]+$gbin{"gold"}{"TF"}{$i}[1]))*($gbin{"gold"}{"TF"}{$i}[0]/
				($gbin{"gold"}{"TF"}{$i}[0]+$gbin{"gold"}{"TF"}{$i}[3])))/
				(($gbin{"gold"}{"TF"}{$i}[0]/($gbin{"gold"}{"TF"}{$i}[0]+$gbin{"gold"}{"TF"}{$i}[1]))+($gbin{"gold"}{"TF"}{$i}[0]/($gbin{"gold"}{"TF"}{$i}[0]+$gbin{"gold"}{"TF"}{$i}[3])))));
				$outnetnodes{$i}[6] = sprintf("%.4f",((2*($gbin{"gold"}{"TF"}{$i}[0]/($gbin{"gold"}{"TF"}{$i}[0]+$gbin{"gold"}{"TF"}{$i}[1]))*($gbin{"gold"}{"TF"}{$i}[0]/
				($gbin{"gold"}{"TF"}{$i}[0]+$gbin{"gold"}{"TF"}{$i}[3])))/
				(($gbin{"gold"}{"TF"}{$i}[0]/($gbin{"gold"}{"TF"}{$i}[0]+$gbin{"gold"}{"TF"}{$i}[1]))+($gbin{"gold"}{"TF"}{$i}[0]/($gbin{"gold"}{"TF"}{$i}[0]+$gbin{"gold"}{"TF"}{$i}[3])))));
			}
			else{
				$tmpout  .= 0;
				$outnetnodes{$i}[6] = 0;
			}
			$tmpout .= "\n";
		}
		else{
			$tmpNOG .= "$i\n";
			$outnetnodes{$i}[3] = 0;
			$outnetnodes{$i}[4] = 0;
			$outnetnodes{$i}[6] = 0;
		}
		delete $resmotNodeTrue{$i};
	}
	
	$tmpoutn.= "$fg Genes (no TF) in graphlets:\ngene";
	for my $j (0..5){
		$tmpoutn.= "\t";
		$tmpoutn.= $j+1;
	}
	$tmpoutn.= "\ttotal\tRGD\tF1\n";
	for my $i (sort keys %resmotNodeTrue){
		if(exists $resmotNodeTrue{$i}){
			if(0 != $resmotNodeTrue{$i}{"tot"}){
				$tmpoutn.= "$i"; 
				for my $j (0..5){
					$tmpoutn.= "\t$resmotNodeTrue{$i}{$j}";
				}
				$tmpoutn.= "\t$resmotNodeTrue{$i}{tot}\t";
				$outnetnodes{$i}[3] = $resmotNodeTrue{$i}{"tot"};
				$outnetnodes{$i}[4] = sprintf("%.4f",($resmotNodeTrue{$i}{"rec"}/($resmotNodeTrue{$i}{"tot"}*6)));
				$tmpoutn  .= sprintf("%.4f\t",($resmotNodeTrue{$i}{"rec"}/($resmotNodeTrue{$i}{"tot"}*6)));
				if($gbin{"gold"}{"noTF"}{$i}[0] != 0){
					$tmpoutn  .= sprintf("%.4f",((2*($gbin{"gold"}{"noTF"}{$i}[0]/($gbin{"gold"}{"noTF"}{$i}[0]+$gbin{"gold"}{"noTF"}{$i}[1]))*($gbin{"gold"}{"noTF"}{$i}[0]/
				($gbin{"gold"}{"noTF"}{$i}[0]+$gbin{"gold"}{"noTF"}{$i}[3])))/
				(($gbin{"gold"}{"noTF"}{$i}[0]/($gbin{"gold"}{"noTF"}{$i}[0]+$gbin{"gold"}{"noTF"}{$i}[1]))+($gbin{"gold"}{"noTF"}{$i}[0]/($gbin{"gold"}{"noTF"}{$i}[0]+$gbin{"gold"}{"noTF"}{$i}[3])))));
					$outnetnodes{$i}[6] = sprintf("%.4f",((2*($gbin{"gold"}{"noTF"}{$i}[0]/($gbin{"gold"}{"noTF"}{$i}[0]+$gbin{"gold"}{"noTF"}{$i}[1]))*($gbin{"gold"}{"noTF"}{$i}[0]/
				($gbin{"gold"}{"noTF"}{$i}[0]+$gbin{"gold"}{"noTF"}{$i}[3])))/
				(($gbin{"gold"}{"noTF"}{$i}[0]/($gbin{"gold"}{"noTF"}{$i}[0]+$gbin{"gold"}{"noTF"}{$i}[1]))+($gbin{"gold"}{"noTF"}{$i}[0]/($gbin{"gold"}{"noTF"}{$i}[0]+$gbin{"gold"}{"noTF"}{$i}[3])))));
				}
				else{
					$tmpoutn  .= 0;
					$outnetnodes{$i}[6] = 0;
				}
				$tmpoutn .= "\n";
			}
			else{
				$tmpNOGn .= "$i\n";
				$outnetnodes{$i}[3] = 0;
				$outnetnodes{$i}[4] = 0;
				$outnetnodes{$i}[6] = 0;
			}
		}
	}
	
	my $tmpouti = '';
	my $tmpNOGi = "$ffg TFs not in graphlets:\n";
	my $tmpoutni = '';
	my $tmpNOGni = "$ffg Genes (no TF) not in graphlets:\n";
	
	
	$tmpouti  .= "$ffg TFs in graphlets:\nTF";
	for my $j (0..12){
		$tmpouti  .= "\t";
		$tmpouti  .= $j+1;
	}
	$tmpouti  .= "\ttotal\n";
	for my $i (@tf0){
		
		if(0 != $resmotNodePred{$i}{"tot"}){
			$tmpouti  .= "$i"; 
			for my $j (0..12){
				$tmpouti  .= "\t$resmotNodePred{$i}{$j}";
			}
			$tmpouti .= "\t$resmotNodePred{$i}{tot}\n";
			$outnetnodes{$i}[5] = $resmotNodePred{$i}{"tot"};
		}
		else{
			$tmpNOGi .= "$i\n";
			$outnetnodes{$i}[5] = 0;
		}
		delete $resmotNodePred{$i};
	}
	
	$tmpoutni .= "$ffg Genes (no TF) in graphlets:\ngene";
	for my $j (0..5){
		$tmpoutni .= "\t";
		$tmpoutni .= $j+1;
	}
	$tmpoutni .= "\ttotal\n";
	for my $i (sort keys %resmotNodePred){
		if(exists $resmotNodePred{$i}){
			if(0 != $resmotNodePred{$i}{"tot"}){
				$tmpoutni .= "$i"; 
				for my $j (0..5){
					$tmpoutni .= "\t$resmotNodePred{$i}{$j}";
				}
				$tmpoutni .= "\t$resmotNodePred{$i}{tot}\n";
				$outnetnodes{$i}[5] = $resmotNodePred{$i}{"tot"};
			}
			else{
				$tmpNOGni .= "$i\n";
				$outnetnodes{$i}[5] = 0;
			}
		}
	}
	
	
	open(F, ">$folder/$outfile") or die "can't open $folder/$outfile";#$countlist{$tmp[0]}{"TP"} .= "$tmp[1] $tmp[2] $tmp[3]\n";
	print F $outtxttable;
	print F "$tmpNOG\n\n";
	print F "$tmpNOGn\n\n";
	print F "$tmpout\n\n";
	print F "$tmpoutn\n\n";
	print F "$tmpNOGi\n\n";
	print F "$tmpNOGni\n\n";
	print F "$tmpouti\n\n";
	print F "$tmpoutni\n\n";
	print F $outtxt;
	close F;
}

#MAIN

getopt();
loadGold();

if($ff ne ''){
	loadPred();
	findMotives_edgeres();
	printres();
}
else{
	findMotives_gold();
	printres_gold();
}

#~ printNets();
