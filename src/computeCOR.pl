#!/opt/perl/5.16.3/bin/perl  -w

use strict; 
use Data::Dumper; 
use feature 'say'; 
use IO::All; 

my $length = 0;

my $in = io("$ARGV[0]"); 
   $in->autoclose(0); 
   while (my $line = $in->getline || $in->getline){
   chomp $line; 
		my @dat = split /\t/, $line;
		my ($start,$end) = ($dat[3],$dat[4]);
		
		if ( $end > $start){
			my $size = ($end - $start);
	        	$length += $size;
		
		} else {
			my $size = ($start - $end);
			$length += $size;
		}
	}
my $lengthMb = ($length / 1000000);
my $rouned = sprintf("%.2f",$lengthMb); 

my $lengthKb = ($length / 1000);
my $rounedlengthKb = sprintf("%.2f",$lengthKb); 

say "COR\t$ARGV[0]\t$rouned (Mb)\t$rounedlengthKb (Kb)\t$length (bp)";

