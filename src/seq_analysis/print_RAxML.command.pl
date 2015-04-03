#!/opt/perl/5.16.3/bin/perl -w

use strict; 
use Carp; 
use feature 'say'; 
use Data::Dumper;



my ($protest, $phy) = @ARGV;

my $model = extractModel($protest); 
my $name = $phy;
   $name =~s/.phy//;

say "raxml -f a  -s $phy -n $name -x 1234 -m PROTGAMMA$model -# 100 -T 8 -p 1023";

# sub
 sub extractModel {
	my ($file) = @_; 

	my $tmp = ''; 
	open FILE, $file || croak "cannot open $file:$!\n";
	while(<FILE>){
	chomp; 
		if (/^Best/){
			my ($mod) = $_ =~/Best Model\s+\:\s+(\w+)/; 
			$tmp = $mod;
		}
	}

	return($tmp);
	close FILE;
}
