#!/opt/perl/5.16.3/bin/perl -w



use IO::All;
use feature 'say';
use strict;


my $ls = io("$ARGV[0]");
$ls->autoclose(0);
while (my $line = $ls->getline || $ls->getline){
chomp $line; 
	my @dat = split /\s+/, $line;
	my ($size,$name) = ($dat[4],$dat[8]);
	
 	unless ($size > '400'){
		say "$size\t$name";
	}

}
