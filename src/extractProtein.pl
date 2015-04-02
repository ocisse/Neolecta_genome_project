#!/opt/perl/5.16.3/bin/perl -w


use Data::Dumper; 
use Carp; 
use strict; 
use feature 'say'; 


my ($parsed,$originalOrthoMCL) = @ARGV; 


# read and store og for which protein should be retrieved

my @groups = (); 

open FILE, $parsed || croak "cannot open $parsed:$!\n";
while(<FILE>){
chomp; 
	push(@groups, $_);
}
close FILE;

# extract Aspergillus and Neurospora seqs

my $i = ''; 

foreach $i (@groups){
	extract($i,$originalOrthoMCL);
}

# sub

sub extract{
	my ($item,$file) = @_; 

	my $x = ''; 
	open DTA,$file || croak "cannot open$file:$!\n";
	while(<DTA>){
	chomp; 
		if (/^$item:/){
			my @tmp = split /\s+/,$_;
			foreach $x (@tmp){
				my @tmp2 = split /\|/, $x;
				my $species = $tmp2[0];
  					if ( $species ~~ /NIRR/){
						#say "$item\t$x";
						
						# remove isoforms
						if ( $tmp2[1] ~~ /\./){
								#say "isoform detected:\t$tmp2[1]";
							}
						else {
						say "NIRR|$tmp2[1]";
					 	}
					}
				}
			}
	}
}
