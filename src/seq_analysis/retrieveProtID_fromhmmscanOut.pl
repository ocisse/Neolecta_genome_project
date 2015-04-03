#!/opt/perl/5.16.3/bin/perl -w

use Data::Dumper;
use Carp; 
use feature 'say';
use strict; 


my $PfamSearch = $ARGV[0]; 

open P,$PfamSearch || croak "cannot open $PfamSearch:$!\n";
while(<P>){
chomp;
	if(/^#/){
		next;
	}
	else{
	my @line = split /\s+/;
	say $line[2];
	}
}
close P;

