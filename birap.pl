#!/usr/bin/perl
use strict;
use warnings;
#use diagnostics;
use Data::Dumper qw (Dumper); 
use Scalar::Util qw(looks_like_number);
#open(STDOUT, '>', 'results.log') or die "Can't open log";
#open(STDERR, '>', 'results_error.log') or die "Can't open log";
my $pileup = [];
my $intlen = 70;
my $output;
my $analysis=0;
my $fasta;
my $anno;
my $term;
my $pro;
my $pgm;
my $help = 0;
my $cur_variable = undef;
foreach my $arg (@ARGV)
{
	if($arg =~ m/--output/i)
	{
		$cur_variable = \$output;
	}elsif($arg =~ m/--analysis/i)
	{
		$cur_variable = \$analysis;
	}elsif($arg =~ m/--fasta/i)
	{
		$cur_variable = \$fasta;
	}elsif($arg =~ m/--anno/i)
	{
		$cur_variable = \$anno;
	}elsif($arg =~ m/--pro/i)
	{
		$cur_variable = \$pro;
	}elsif($arg =~ m/--term/i)
	{
		$cur_variable = \$term;
	}elsif($arg =~ m/--pgm/i)
	{
		$cur_variable = \$pgm;
	}elsif($arg =~ m/--help/i)
	{
		$help = 1;
	}elsif($arg =~ m/--intlen/i)
	{
		$cur_variable = \$intlen;
	}elsif($arg =~ m/--pileup/i)
	{
		$cur_variable = $pileup;
	}elsif($cur_variable && ref($cur_variable) eq "ARRAY")
	{
		push(@{$cur_variable}, $arg);
	}elsif($cur_variable && ref($cur_variable) eq "SCALAR")
	{
		$$cur_variable = $arg;
	}
}
if (!(scalar @ARGV) || $help) 
{
	print "\n---Bacterial Intergenic Region Analysis Pipeline Help---\n\n";
	print "Basic usage:-\n";
	print "\$ perl birap.pl --analysis [1 or 2 or 3] --fasta [genome.fasta] --pileup [file1 file2 ... or *.pileup] --pgm [proteogenomic.gff] --anno [annotation.gff OR annotation.ptt] --output [output filename base)]  --pro [promoters.coords] --term [terminators.coords] --intlen [value (default:70)]\n\n";
	print "*********************************************************\n";
	print "Description of parameters:\n\n";
	print "Mandatory parameters: --analysis, --fasta, --pileup and/or --pgm based on \"analysis\", --anno, --output\n\n";
	print "--analysis\tChoose the type of analysis you wish to perform from among the following options:\n";
	print "\t\t 1 for analysis of intergenic regions using RNA-seq data\n";
	print "\t\t 2 for analysis of intergenic regions using Proteomics data\n";
	print "\t\t 3 for analysis of intergenic regions using integrated omics approach\n";
	print "--fasta\tgenome file\n";
	print "--pileup\tone or more bowtie alignment output files in pileup format (obtained using samtools 'mpileup') \n\t(please provide space-separated absolute path for each file with respect to current working directory. 'pathtodir/*.pileup' can also be used)\n";
	print "--pgm\tGFF file from Proteogenomic mapping tool or similar tools providing the locus of expressed peptides in GFF format\n";
	print "--anno\tannotation file in \.gff or \.ptt format\n\n";
	print "--output\tbase filename to be used for all output files\n";	
	print "Optional arguments:\n\n";
	print "--pro\tfile containing loci of promoters in the genome in \.coords format (see below)**\n";
	print "--term\tfile containing loci of terminators in the genome in \.coords format (see bewlow)**\n";
	print "--intlen\tminimum length of an intergenic region to be considered as a region of interest for downstream RNA-seq data analysis (default cutoff is 70bp)\n";
	print "\n**.coords FORMAT:\n";
	print "<START>\\t<STOP>\\t<STRAND>\\t<DESCRIPTION>\n";
	print "Sample entry promoter locus in \.coords file:\n";
	print "14024	14052	+	sigmaA_15bp\n\n";	
	print "Sample entry for terminator locus in \.coords file:\n";
	print "3046632	3046596	-	TERM_001\n";
	print "\nIf files containing locus of promoters and terminators in the genome are provided, the program will identify putative non-coding RNA based on the existing annotation and the expression profile generated from RNA-seq\n\n";
	print "*********************************************************\n";	
	exit;
}
else
{
	if($analysis == 1)
	{
		print "\n******Starting RNA-seq data analyis******\n";
		if (scalar @{$pileup}==0)
		{
			die "******\nYou didnt provide any pileup files!\n******\n";
		}
	}
	elsif($analysis == 2)
	{	
		print "\n******Starting Proteogenomic data analyis******\n";
		if(!$pgm)
		{
			die "******\nYou didnt provide PGM files!\n******\n";	
		}	
	}
	elsif($analysis == 3)
	{
		print "\n******Starting integrated -omics data analyis******\n";
		if ((scalar @{$pileup}==0)||(!$pgm))
		{
			die "******\nYou didnt provide pileup or PGM files!\n******\n";
		}
		print "\nPerforming validation of annotation using integrated omics data!\n";
	}
	else
	{
		die "******\nYou didnt choose method of data analysis (1,2,3)!\n******\n";
	}
	if (!$output)
	{
		die "******\nYou didn't provide a base file name for output!\n******\n";
	}
	if (!$analysis)
	{
		die "******\nYou didn't specify an analysis option: 1, 2, or 3!\n******\n";
	}
	if (!$fasta)
	{
		die "******\nYou didn't provide a genome sequence!\n******\n";
	}if (!$anno)
	{
		die "******\nYou didn't provide an annotation!\n******\n";
	}	
}
my $chr="";		#For keeping track of chromosome/genome name
my @pileup_filled=();	#For containing normalized expression profile of all RNA-seq samples
my @genex=();		#For containing gene expression profile of annotated regions
my @intergenic=();	#For expressed intergenic regions (RNA-seq)
my @protex=();		#For containing protein expression profile of annotated regions
my @intergenic_pep=();	#For expressed intergenic peptides (Proteomics)
my @integrated=();	#For validating annotation using both RNA-seq and Proteomics data
my @pgmex=();		#For validation of proteogenomic mapping data using RNA-seq
my @putative_sRNA=();	#For novel sRNA
my @peptides=();	#For containing proteogenomic mapping data (input) 
my @annotation=();	#For containing annotation (input)

print "\nProcessing genome file....\n";
open GNM, "$fasta" or die "Cannot find the genome sequence file!";
my @genome=<GNM>; 
close(GNM); 
chomp(@genome);
my $genome_header=shift(@genome);
my $genome_seq=join('',@genome);
$genome_seq=~s/\r//g;
my $genome_len=length($genome_seq);
print "Size of genome :".$genome_len."bp\n";
print "\nProcessing annotation file....\n";
#Converting annotation file (.ptt or .gff) to this format : <START>\t<STOP>\t<STRAND>\t<DESCRIPTION>
if($anno=~/ptt/)	#Check if annotation file is a .ptt file
{
	open ANN, "$anno" or die "Cannot find PTT file!";
	while(<ANN>)
	{
		if($_=~/^\#/)
		{
			next;
		}
		chomp($_);	
		$_=~s/\r//g;
		if(($_ =~ /^\d/) && ($_=~ /\.\./))
		{
			my @prot=split(/\t/,$_);
			if(scalar @prot!=9) {print "Unaccepted PTT file format, please verify!";exit;}		
			my @coord=split(/\.\./,$prot[0]); if(($coord[0]!~/\d/)||($coord[1]!~/\d/)){print "Unaccepted PTT file format, please verify!\n";exit;}
			my $description=join(";",$prot[3],$prot[8]);
			my $locus=join("\t",$coord[0],$coord[1],$prot[1],$description);
			push(@annotation,$locus);
		}
	}close(ANN);	
	print "Number of annotated regions:".scalar(@annotation)."\n";
}
elsif($anno=~/gff/)	#Check if annotation file is a .gff file
{
	open ANN, "$anno" or die "Cannot find GFF file!\n";
	while(<ANN>)
	{
		if($_=~/^\#/)
		{
			next;
		}
		chomp($_);
		$_=~s/\r//g;
		my @gene=split(/\t/,$_); if((scalar @gene!=9)||($gene[3]!~/\d/)||($gene[4]!~/\d/)) {print "Unaccepted GFF file format, please verify!\n";exit;}		
		my $locus=join("\t",$gene[3],$gene[4],$gene[6],$gene[8]);
		push(@annotation,$locus);
	}close(ANN);	
	print "Number of annotated regions:".scalar(@annotation)."\n";
}
else
{
	print "Please enter a valid annotation file (Filename should countain \.ptt or \.gff to differentiate annotation file type!!\n";
	exit;
}
if($analysis == 1)	#CHOICE 1, RNA-seq data analysis
{
	processPileup();
	rna_analysis();
	print "Analysis completed!!\n";	
}
elsif($analysis == 2)	#CHOICE 2, Proteogenomic data analysis
{
	processPeptides();
	pgm_analysis();
	print "Analysis completed!!\n";
}
else #CHOICE 3, Intergrated -omics data analysis
{
	processPileup();
	rna_analysis();
	processPeptides();
	pgm_analysis();
	integrated_analysis();
	print "Analysis completed!!\n";
}
print "Time taken: ";
print time-$^T; print "s\n\n";
#Processing PGM file
sub processPeptides
{
	print "\nProcessing PGM file....\n\n";
	open FILE, "$pgm" or die "Cannot find PGM file!!";
	while(<FILE>)
	{
		if ($_ =~ /^\#/)
		{	
			next;
		}
		chomp($_);
		$_=~s/\r//g;		
		my @region=split(/\t/,$_);
		my $type="N/A";
		if($region[1]=~/RTP/)
		{
			$type="RTP";
		}
		elsif($region[1]=~ /ePST/)
		{
			$type="ePST";
		}
		my $description= join(";",$region[8],$type);
		my $loc="";	
		if($region[3]<$region[4])
		{	
			$loc=join("\t",$region[3],$region[4],$region[6],$description);	
		}
		else
		{
			$loc=join("\t",$region[4],$region[3],$region[6],$description);	
		}
		push(@peptides,$loc);	
	}close(FILE); 
	#print "Number of peptides and ePSTs:".scalar(@peptides)."\n";
}
#Processing pileup file/s
sub processPileup
{
	print "\nProcessing pileup file/s....\n\n";
	if((scalar @{$pileup}) == 1)		#If only one pileup file is provided
	{	
		open FILE, "<@$pileup[0]" or die $!;
		my $firstline=<FILE>;
		close(FILE);
		my @line=split(/\t/,$firstline); 
		my $chr=$line[0];
		my $samples=(scalar(@line)-3)/3;#Counting number of replicates in the given pileup file
		if($samples==1)			
		{
			print "Working on file: @$pileup[0]\n";
			print "Pileup file contains $samples replicate\n";
			open FILE, "@$pileup[0]", or die $!;
			my $start=1; my $end=$genome_len;
			my @expressed=();	
			my @temp_pileup_filled=();
			while(<FILE>)
			{
				chomp($_);
				my @data=split("\t",$_);
				if($data[0] ne $chr)
				{
					print "Multiple genome/chromosome identifiers detected in pileup file, exiting now!!";
					exit;
				}				
				while($start<$data[1])
				{
					push(@temp_pileup_filled,0);
					$start++;
				}
				if(looks_like_number($data[3]))
				{				
					if($data[3] !=0){push(@expressed,$data[3]);}
					push(@temp_pileup_filled,$data[3]);$start++;
				}
				else
				{
					push(@temp_pileup_filled,0);
					$start++;
				}
			}close(FILE);
			print "Completed reading file: @$pileup[0]\n";
		#	print "Number of bases expressed in pileup file: ".scalar(@expressed).",genome size: $genome_len\n";
			while($start<=$end)
			{
				push(@temp_pileup_filled,0);
				$start++;
			}							
			my @sorted= sort {$a <=> $b} @expressed;
			my $percentile=int(($#sorted*10)/100);
			my $cutoff=$sorted[$percentile];
		#	print "Basal expression cutoff is: ".$cutoff."\n";
		#	print "Lowest RPB in the sample : ".$sorted[0]."\n"; 
		#	print "Highest RPB in the sample : ".$sorted[$#sorted]."\n"; 
		#	print "Performing binary transformation...\n\n";	
			foreach my $base(@temp_pileup_filled)	#Transform expression profile into a binary notation, 1:Expressed, 0:Not Expressed
			{
				if($base>=$cutoff)
				{
					push (@pileup_filled,1);
				}
				else
				{
					push(@pileup_filled,0);
				}
			}			
		}
		elsif($samples>1)		#Multiple replicates merged into one pileup file
		{
			my @replicates=();
			my @temp_pileup_filled=();
			print "Working on file: @$pileup[0]\n";	
			print "Pileup file contains $samples replicates\n";
			print "Reading from file, please be patient, this may take a while depending on the number of replicates in your pileup file!!\n";
			open FILE, "@$pileup[0]", or die $!;
			while(<FILE>)
			{
				chomp($_);
				my @data=split(/\t/,$_);
				if($data[0] ne $chr)
				{
					print "Multiple genome/chromosome identifiers detected in pileup file, exiting now!!\n";
					exit;
				}
				my $i=0; #replicate count				
				for(my $j=3;$j<$#data;)	#iterating through columns in pileup file for each replicate and identifying RPB 
				{
					if(looks_like_number($data[$j]))
					{
						if($data[$j] !=0){$replicates[$i][$data[1]]=$data[$j];}	#actual RPB
						$temp_pileup_filled[$i][$data[1]]=0; #actual RPB
					}$j+=3;$i++;
				}
			}close(FILE);
			print "Completed reading from file: @$pileup[0]\n";
			print "\n\nAnalyzing expression across replicates \n";
			for(my $i=0;$i<$samples;$i++)
			{
				print "\nProcessing expression profile of replicate ".($i+1)."\n";
				my @expressed=();
				for(my $j=1;$j<=$genome_len;$j++)				
				{
					if($replicates[$i][$j])
					{
						push(@expressed,$replicates[$i][$j]);
					}
					else
					{
						$temp_pileup_filled[$i][$j]=0;
					}
				}
			#	print "Number of bases expressed in pileup file: ".(scalar @expressed)."\n";					
				my @sorted= sort {$a <=> $b} @expressed;
				my $percentile=int(($#sorted*10)/100);
				my $cutoff=$sorted[$percentile];
			#	print "Basal expression cutoff for sample ".($i+1).": ".$cutoff."\n";
			#	print "Lowest RPB in sample ".($i+1).": ".$sorted[0]."\n"; 
			#	print "Highest RPB in the sample ".($i+1).": ".$sorted[$#sorted]."\n";
			#	print "Performing binary transformation...".($i+1)."\n\n";				
				for(my $j=1;$j<=$genome_len;$j++)				
				{
					if($replicates[$i][$j])
					{	
						if($replicates[$i][$j]>=$cutoff)
						{
							$temp_pileup_filled[$i][$j]=1;
						}
						else
						{
							$temp_pileup_filled[$i][$j]=0;
						}
					}
				}
			}
						
			for(my $j=1;$j<=$genome_len;$j++)			#Transforming expression into binary representation
			{
				my $expr=0;
				my %counts=();
				for(my $i=0;$i<$samples;$i++)
				{
					$counts{$temp_pileup_filled[$i][$j]}++;	#Calculate MODE of binary representation of expression across samples 
				}
				if(!$counts{0})
				{
					$expr=1;
				}elsif(!$counts{1})
				{
					$expr=0;
				}else
				{
					if($counts{1}>=$counts{0})
					{
						$expr=1;
					}
					else
					{
						$expr=0;
					}
				}
				push (@pileup_filled,$expr);
			}
		}
	}
	elsif(scalar @{$pileup}>=2)	#If more than one pileup file is provided (one per replicate)
	{
		my $chr="";	
		open FILE, "<@$pileup[0]" or die $!;
		my $firstline=<FILE>;
		my @line=split(/\t/,$firstline); 
		$chr=$line[0];
		close(FILE);
		my $samples=scalar @{$pileup};
		my @temp_pileup_filled=();	
		for(my $i=0;$i<$samples;$i++) 
		{
			print "Working on file: @$pileup[$i]\n";
			my @expressed=();	
			my $start=1; my $end=$genome_len;
			my @temp=();			
			open FILE, "@$pileup[$i]" or die $!;
			while(<FILE>)
			{
				chomp($_);
				my @data=split(/\t/,$_);
				if($data[0] ne $chr)
				{
					print "Multiple genome/chromosome identifiers detected in pileup file, exiting now!!";
					exit;
				}				
				while($data[1]>$start)
				{
					push(@temp,0);
					$start++;
				}	
				if(looks_like_number($data[3]))
				{				
					if($data[3] !=0){push(@expressed,$data[3]);}
					push(@temp,$data[3]);$start++;
				}
				else
				{
					push(@temp,0);
					$start++;
				}
			}close(FILE);
			while($start<=$end)
			{
				push(@temp,0);
				$start++;
			}							
			my @sorted= sort {$a <=> $b} @expressed;
			my $percentile=int(($#sorted*10)/100);
			my $cutoff=$sorted[$percentile];
		#	print "Basal expression cutoff is: ".$cutoff."\n";
		#	print "Lowest RPB in the sample : ".$sorted[0]."\n"; 
		#	print "Highest RPB in the sample : ".$sorted[$#sorted]."\n";
			my $express=0;
			for(my $j=1;$j<=$genome_len;$j++)	#Transform expression profile into a binary notation, 1:Expressed, 0:Not Expressed
			{
				if($temp[$j-1]>=$cutoff)
				{
					$temp_pileup_filled[$i][$j]=1;$express++;
				}
				else
				{
					$temp_pileup_filled[$i][$j]=0;
				}
			}
		#	print "Total number of expressed bases in sample : $express\n";	
		#	print "Performing binary transformation...\n\n";			
		}
		for(my $j=1;$j<=$genome_len;$j++)	#Pooling expression across files
		{
			my %counts=();
			my $expr=0;
			for(my $i=0;$i<$samples;$i++)
			{
				$counts{$temp_pileup_filled[$i][$j]}++;	#Calculate MODE of binary representation of expression across samples 
			}
			if(!$counts{0})
			{
				$expr=1;
			}elsif(!$counts{1})
			{
				$expr=0;
			}else
			{
				if($counts{1}>=$counts{0})
				{
					$expr=1;
				}
				else
				{
					$expr=0;
				}
			}
			push (@pileup_filled,$expr);
		}
	}
}
sub rna_analysis
{
	print "Computing genome expression profile from RNA-seq data!\n";
	my @temp=@pileup_filled; 
	my @eir=();
	my $filename=join("_",$output,"GenEx.tab");	
	open OUT, ">$filename" or die $!;
	print OUT "\#start\tstop\tstrand\tdescription\texpression\n";
	my $expressed_count=0;
	foreach my $reg(@annotation)		#Checking if annotated region is expressed
	{
		my @gene=split(/\t/,$reg);
		my $expr_bases=0;
		my $gene_len=($gene[1]-$gene[0]+1);
		my $status="";
		for(my $i=($gene[0]-1);$i<$gene[1];$i++)
		{
			if($pileup_filled[$i]==1)
			{
				$expr_bases++;
				$temp[$i]=2;	#Marking location as annotated with a "2"
			}
		}
		if(($expr_bases/$gene_len) >= 0.6)	#Checking if at least 60% gene is expressed
		{
			$status="Expressed";
			$expressed_count++;
		}	
		else
		{
			$status="Not_Expressed";
		}
		print OUT $reg."\t".$status."\n";
		my $gen_ex=join("\t",$reg,$status);
		push(@genex,$gen_ex);
	}close(OUT); print "Done!!\n";
	print "Expressed annotated regions: $expressed_count\n";	
	my $start=0; 
	my $stop=0;
	my $count=0; 
	my $eir_count=0;
	print "Identifying expressed intergenic regions (EIRs)...\n";
	for(my $i=0;$i<$#temp;$i++)		
	{	
		if($temp[$i]==1)
		{
			$start=$i;
			$count=0;	
			while($temp[$i]==1)
			{
				$i++;$count++;
			}$i--;
			if($count>=$intlen)
			{
				$stop=$i; 
				$eir_count++;
				my $eir_name=join("","EIR",$eir_count);
				my $region=join("\t",$eir_name,$start,$stop);
				push(@eir,$region);
			}
		}
	}
	print "Checking association of EIR with existing annotation...\n";
	foreach (@eir)
	{	
		my @reg=split(/\t/,$_);
		my $rs=$reg[1];
		my $re=$reg[2];
		my $flag=0;
		my $note="";
		foreach my $gen(@annotation)
		{
			my @line=split(/\t/,$gen);
			my $gs=$line[0];
			my $ge=$line[1];
			if(($rs<$gs) && (($re+10)>=$gs))
			{
				$note="5' UTR of $line[3]";	
				$flag=1;
				last;
			}
			elsif(($ge<$re) && (($ge+10)>=$rs))
			{
				$note="3' UTR of $line[3]";	
				$flag=1;
				last;
			}
		}
		if($flag==0)
		{
			$note="Intergenic";
		}
		my $eir_desc=join("\t",$_,$note);
		push(@intergenic,$eir_desc);	
	}
	if(scalar(@intergenic)>0)
	{
		my $filename1=join("_",$output,"EIR.tab");	
		open OUT1, ">$filename1" or die $!;
		print OUT1 "\#EIR_description\tEIR_start\tEIR_stop\tEIR_association\n";
		my $filename2=join("_",$output,"EIR.fasta");
		open OUT2, ">$filename2" or die $!;
		foreach (@intergenic)
		{
			print OUT1 "$_\n";
			my @line=split(/\t/,$_);
			print OUT2 "\>$line[0]\;$line[1]\.\.$line[2]\n";
			my $len=$line[2]-$line[1];
			my $seq=substr($genome_seq, $line[1],$len);
			for(my $j=0;$j<length($seq);)
			{	
				my $sub_seq= substr $seq,$j,60;
				print OUT2 "$sub_seq\n";	
				$j+=60;
			}			
		}close(OUT1);close(OUT2);
	}
	my @promoters=();
	my @terminators=();	
	if($pro)	
	{
		open FILE, "$pro" or die "Couldn't find Promoter file!!";
		@promoters=<FILE>;
		close(FILE);chomp(@promoters);
	}
	if($term)
	{
		open FILE, "$term" or die "Couldn't find Terminator file!!";
		@terminators=<FILE>;
		close(FILE);chomp(@terminators);
	}
	if((scalar(@promoters)>0) || (scalar(@terminators)>0)) 
	{
		print "Identifying sRNA , if any...!!\n";
		foreach my $re(@intergenic)
		{
			if($re !~ /Intergenic/i)
			{
				next;
			}		
			my $signal=0;
			my @regn=split(/\t/,$re);
			my $r_start=$regn[1];
			my $r_end=$regn[2];
			my $sRNA=join("\t",$regn[0],$regn[1],$regn[2]);
			if(scalar(@promoters)>0)
			{
				foreach my $pr(@promoters)
				{
					if ($pr !~ /^\d/)
					{
						next;
					}
					$pr=~s/\r//g;
					my @prom=split(/\t/,$pr);
					my $start=$prom[0];
					my $end=$prom[1];
					my $strand=$prom[2];
					if($strand eq "+")	
					{
						if(($start+35)>=$r_start && ($end-20)<=$r_start)
						{
							$signal=1; 
							$sRNA.="\t$pr";last;
						}		
					}	
					elsif($strand eq "-")
					{
						if(($start-35)>=$r_start && ($end+20)<=$r_end)
						{
							$signal=1; 
							$sRNA.="\t$pr";last;
						}		
					}
				}
			}
			if(scalar(@terminators)>0)
			{
				foreach my $te(@terminators)
				{
					if ($te !~ /^\d/)
					{
						next;
					}
					$te=~s/\r//g;
					my @term=split(/\t/,$te);
					my $start=$term[0];
					my $end=$term[1];	
					my $strand=$term[2];
					if($strand eq "+")
					{
						if(($start-20)>=$r_end && ($end+20)<=$r_end)
						{
							$signal=1;
							$sRNA.="\t$te";last;
						}
					}
					elsif($strand eq "-")
					{
						if($start+20>=$r_start && ($end-20)<=$r_start)
						{
							$signal=1;
							$sRNA.="\t$te";last;
						}
					}
				}	
			}	
			if($signal == 1)	###RESULTS:sRNA
			{
				push(@putative_sRNA,$sRNA);		
			}
		}
		if(scalar(@putative_sRNA>0))
		{
			my $filename=join("_",$output,"putative_sRNA.tab");	
			open OUT, ">$filename" or die $!;
			print OUT "\#EIR_description\tEIR_start\tEIR_stop\tsignal_start\tsignal_stop\tsignal_strand\tsignal_description\tsignal_start\tsignal_stop\tsignal_strand\tsignal_description\n";
			foreach (@putative_sRNA){ print OUT $_."\n";} close(OUT);
			$filename=join("_",$output,"putative_sRNA.fasta");	
			open OUT, ">$filename" or die $!;	
			foreach (@putative_sRNA)
			{	
				my @line=split("\t",$_);
				my $id=$line[0];
				my $len=$line[2]-$line[1]+1;
				my $seq=substr($genome_seq, $line[1],$len);
				print OUT "\>$line[0]\;$line[1]\.\.$line[2]\n";
				for(my $j=0;$j<length($seq);)
				{
					my $sub_seq= substr $seq,$j,60;
					print OUT "$sub_seq\n";
					$j+=60;
				}
			}close(OUT);	
		}
	}	
	return;
}
sub pgm_analysis
{
	print "Computing genome expression profile from proteomics data!\n";
	my %pep_expr=();
	my @eip_id=();
	my $expressed_count=0;
	foreach my $gene(@annotation)
	{
		my @line=split(/\t/,$gene);	
		my $gs=$line[0];
		my $ge=$line[1];
		my $count=0;
		my $g_strand=$line[2];
		my $status="";
	#	my $direction="";
		foreach my $peptide(@peptides)
		{
			if($peptide !~ /RTP/)		#Excluding ePSTs
			{
				next;
			}
			my @pep=split(/\t/, $peptide);
			my $ps=$pep[0];
			my $pe=$pep[1];
			my $p_strand=$pep[2];
	#		if($p_strand eq $g_strand)
	#		{
	#			$direction="sense";
	#		}
	#		else
	#		{
	#			$direction="antisense";
	#		}
			if($pe<$gs || $ps>$ge)	#0
			{
				if(exists ($pep_expr{$peptide}))
				{
					next;
				}
				else
				{
					$pep_expr{$peptide}=0;
					next;
				}
			}
			elsif($ps>=$gs && $pe<=$ge)	#1
			{
				$count++;$pep_expr{$peptide}=1;
			}	
			elsif($ps<=$gs && $pe>=$gs && $pe<=$ge)	#2
			{
				$count++;$pep_expr{$peptide}=1;
			}
			elsif($ps>=$gs && $ps<=$ge && $pe>$ge)	#3
			{
				$count++;$pep_expr{$peptide}=1;
			}
			elsif($ps<$gs && $pe>$ge)	#4
			{
				$count++;$pep_expr{$peptide}=1;
			}
			elsif($ps<=$ge && $pe>$ge)	#5
			{
				$count++;$pep_expr{$peptide}=1;
			}
			else
			{
				print "shouldn't be here!!";
			}
		}
		if($count==0)
		{
			$status="Not_Expressed";
		}
		else	
		{
			$expressed_count++;
			$status="Expressed";		
		}
		my $result=join("\t",$gene,$count,$status);		####RESULTS Protein annotation validation
		push(@protex,$result);
	}
	print "Expressed annotated regions: $expressed_count\n";
	my $filename=join("_",$output,"ProtEx.tab");
	open OUT, ">$filename" or die $!;
	print OUT "\#start\tstop\tstrand\tdescription\tpeptide_count\texpression\n";
	foreach (@protex){ print OUT $_."\n";} close(OUT);
	print "Identifying expressed intergenic peptides (EIPs)...\n";
	print "Checking EIP association with existing annotation...\n";
	for my $pep (keys %pep_expr)
	{
		if($pep_expr{$pep} == 0)			#Intergenic peptide
		{
			my @reg=split(/\t/,$pep);
			push(@eip_id,$reg[3]);
			my $rs=$reg[0];
			my $re=$reg[1];
			my $flag=0;
			my $note="";
			foreach my $gen(@annotation)
			{
				my @line=split(/\t/,$gen);
				my $gs=$line[0];
				my $ge=$line[1];
				if($reg[2] eq $line[2])		#Check to see if peptide is in the same strand as the gene 
				{
					if(($rs<$gs) && (($re+30)>=$gs))
					{
						$note="Expressed at 5' end $line[3]";	
						$flag=1;
						last;
					}
					elsif(($ge<$re) && (($ge+30)>=$rs))
					{
						$note="Expressed at 3' end $line[3]";	
						$flag=1;
						last;
					}
				}
			}
			if($flag==0)
			{
				$note="Intergenic";
			}
			my $eip_desc=join("\t",$pep,$note);
			push(@intergenic_pep,$eip_desc);			####RESULTS Intergenic Peptides 		
		}
	}
	if(scalar(@intergenic_pep)>0)
	{
		my $filename=join("_",$output,"EIP.tab");	
		open OUT, ">$filename" or die $!;
		print OUT "\#EIP_start\tEIP_stop\tEIP_strand\tEIP_description\tEIP_association\n";
		foreach (@intergenic_pep){ print OUT $_."\n";} close(OUT);	
	}
	if (scalar(@intergenic_pep)>0)		#Generating fasta sequences intergenic ePST
	{
		my $filename=join("_",$output,"Intergenic_ePST.fasta");	
		open OUT, ">$filename" or die $!;	
		foreach my $inpep(@intergenic_pep)
		{
			my @data=split(/\t/,$inpep);
			my @in_id=split(";",$data[3]);
			foreach my $pep(@peptides)
			{
				if(( $pep =~ $in_id[0]) && ($pep =~ /ePST/))
				{
					my @line=split(/\t/,$pep);
					my $len=$line[1]-$line[0]+1;
					my $seq=substr($genome_seq, $line[0],$len);
					print OUT "\>$line[3]\;$line[0]\.\.$line[1]\n";
					for(my $j=0;$j<length($seq);)
					{
						my $sub_seq= substr $seq,$j,60;
						print OUT "$sub_seq\n";
						$j+=60;
					}
					
				}
			}
		}close(OUT);		
	}	
	return;
}
##
sub integrated_analysis
{
	if((scalar(@annotation) == scalar(@genex)) && (scalar(@annotation) == scalar(@protex)))		#Comparing expression pattern of annotated regions
	{
		print "Intergrating expression profiles from -omics data!\n";
		for(my $i=0;$i<=$#annotation;$i++)
		{	
			my @gex=split(/\t/,$genex[$i]);
			my @pex=split(/\t/,$protex[$i]);
			my $result=join("\t",$annotation[$i],$gex[$#gex],$pex[$#pex]);
			push(@integrated,$result);
		}
		my $filename=join("_",$output,"IntegratedEx.tab");	
		open OUT, ">$filename" or die $!;
		print OUT "\#start\tstop\tstrand\tdescription\trna-seq_expression\tproteomics_expression\n";
		foreach (@integrated){ print OUT $_."\n";} close(OUT);	
	}	
	else
	{
		print "Something terribly wrong!! Size of arrays (annotation, genex and protex) don't match!!";
	}
	foreach my $pep(@peptides)		#Validation of expression of peptides and ePSTs at the transcript level
	{
		my @data=split(/\t/,$pep);
		my $expr_bases=0;
		my $status="";
		for(my $i=($data[0]-1);$i<$data[1];$i++)
		{
			if($pileup_filled[$i]==1)
			{
				$expr_bases++;
			}
		}
		my $gene_len=$data[1]-$data[0]+1;
		if(($expr_bases/$gene_len)>=0.6)	#Checking if at least 60% is expressed
		{
			$status="Expressed";
		}	
		else
		{
			$status="Not_Expressed";
		}
		my $pgm_ex=join("\t",$pep,$status);
		push(@pgmex,$pgm_ex);
	}
	my $filename=join("_",$output,"PGMEx.tab");	#####RESULTS PGM EXPRESSION
	open OUT, ">$filename" or die $!;
	print OUT "\#start\tstop\tstrand\tdescription\texpression\n";
	foreach (@pgmex){ print OUT $_."\n";} close(OUT);	
	return;
}

