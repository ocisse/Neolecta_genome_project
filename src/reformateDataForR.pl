#!/opt/perl/5.16.3/bin/perl -w

use strict; 
use Data::Dumper; 
use feature 'say'; 
use IO::All; 

my $i = 0; 

say "ID\tAVELEN\tNUM";
my $in = io("$ARGV[0]"); 
   $in->autoclose(0);
   while (my $line = $in->getline || $in->getline){
   chomp $line;
	
	my @dat = split /\t/, $line; 
	my $id = $dat[0];
	
	# compute average intron length for a given gene
	shift @dat;


	if ($dat[0] != 0){
		my ($averageLen,$numOfIntrons) = computeAverage(@dat);
	 	say "$id\t$averageLen\t$numOfIntrons";  
		}
	else {
		# no intron => no need
		#say "$id\t0\t0";
	}
}
	

# subs
sub computeAverage {
	my (@line) = @_; 
	
	my $i = 0;
	my $count = 0; 
	my $t = 0;


	foreach $i (@line){
		$count += $i;
		$t++;
	}
	my $ave = ($count/$t);
	
	return ($ave,$t);
}
