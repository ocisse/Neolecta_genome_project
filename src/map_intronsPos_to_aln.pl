#!/opt/perl/5.16.3/bin/perl -w



use Data::Dumper; 
use Carp; 
use feature 'say'; 
use strict; 
use IO::All;

# data


# /map_intronsPos_to_aln.pl test.fa > test.map

# read multiple alignemt 
my  $aln = io("$ARGV[0]");
$aln->autoclose(0);
while (my $inLine = $aln->getline || $aln->getline){
chomp $inLine;
	if ($inLine =~ m/^>/){
	   my @line = split /\|/, $inLine;

	   my ($species,$prot) = @line[0..1];
		     $species  =~s/>//;
		if (-e "$species\_$prot.gff"){
			my ($GeneIdStart,$IntStart) = getIntronsNameAndPos("$species\_$prot.gff");
			my %GeneIdStart = %$GeneIdStart;
			my %IntStart = %$IntStart;
	
			# do gene start pos - intro start pos => to corr to CDS pos and genome pos
			# require for Malin
		 	my $adjustedIntronPos = adjust($GeneIdStart{"$species\_$prot.gff"},%IntStart);
			say ">$species\_$prot\t/organism=$species\t {i $adjustedIntronPos i}";

		}
		else {
			say ">$species\_$prot\t/organism=$species\t {i  i}";
		}
	}
	else {
		say $inLine;	
	}
}

# ------------------------------------------# 
# sub
# ------------------------------------------#

sub adjust {
	my ($lenToremove,%int) = @_; 

	my @Len = (); 

	# loop over intron start and remove the gene start pos
	my $x = ''; 
	foreach $x (keys %int){
		my $adjustedLen = ($int{$x} - $lenToremove);
		push(@Len,$adjustedLen);
	}

my @sort = sort {$a <=> $b} @Len;

my $join = join(',', @sort); 
return $join;
}

sub getIntronsNameAndPos {
	my ($gff) = @_;

	my (%geneData,%intronsData) = ();	
	
	open GFF,$gff || croak "cannot open $gff from here:$!\n";
	while(<GFF>){
	chomp; 
		my @line = split /\t/;
		
		if(/gene/){
			my $gStart = $line[3];	
			$geneData{$gff} = $gStart;
		}
		elsif ( $line[2] ~~ 'intron'){
			 my ($iStart,$iId) = ($line[3],$line[8]);
		         $intronsData{$iId} = $iStart;
		}
	}
return(\%geneData,\%intronsData);
close GFF;
}
