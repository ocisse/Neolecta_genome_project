#!/opt/perl/5.16.3/bin/perl -w

use feature 'say'; 
use Carp; 
use Data::Dumper; 
use strict; 

my $tre = shift; 

my $i = ''; 

open T,$tre || croak "cannot read $tre:$!\n";
while(<T>){
chomp; 
	my @line = split /\:/;
	
	foreach $i (@line){

		if ( $i =~ m/\|/){
			my @tmp = split /\|/, $i;
			print "$tmp[0]:";
		}
		else{
			print "$i:";
		}
	}	
	


}
close T;
 print "\n";
