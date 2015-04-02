#!/opt/perl/5.16.3/bin/perl -w

use strict;
use feature 'say'; 
use Data::Dumper; 
use Carp; 

my $fasta = shift; 

open F,$fasta || croak "$!\n";
while(<F>){
chomp; 

	if (/^>/){
		my @line = split /\|/;
		say $line[0];
	}
	else {
	say;
	}
}
close F;

