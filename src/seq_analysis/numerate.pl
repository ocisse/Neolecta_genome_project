#!/usr/bin/perl -w


use strict;
use Carp;
use feature 'say';

my $dir = $ARGV[0];
my $count = 0;

opendir(DIR,$dir) || croak "cannot open $dir:$!\n";
while( my $file = readdir (DIR)){
	
	if ( $file ~~ '.'){
		next;	
	} 
	elsif ( $file ~~ '..') {
		next;
	}
	else {
	my $cmd1="cp $file file-$count";
	execute($cmd1);	
	$count++;
	
	my $cmd2="rm $file";
 	execute($cmd2);	
	}
}
closedir(DIR);


# sub 

sub execute {
	my ($commd) = @_;
	warn"$commd\n";
	system($commd)==0 || croak "cannot execute $commd:$!\n";
}
