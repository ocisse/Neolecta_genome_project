#!/opt/perl/5.16.3/bin/perl -w

use strict; 
use Data::Dumper; 
use Carp; 
use feature 'say';

my $in = shift; 

open I,$in || croak "$!\n";
while(<I>){
chomp; 

	if (/^>/){
		my @line = split /\|/,$_;
		say ">$line[1]";
	}
	else {
	say;
	}
}
close I; 
