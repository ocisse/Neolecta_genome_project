#!/opt/perl/5.16.3/bin/perl -w


use IO::All;
use feature 'say'; 
use Data::Dumper; 
use Carp; 


my $in = io("$ARGV[0]"); 
$in->autoclose(0); 
while (my $line = $in->getline || $in->getline){
chomp $line; 
	push(@tmp,$line); 

}
my %hash = (); 

my $i = ""; 
foreach $i (@tmp){
	$hash{$i}++;
}
#say Dumper \%hash;
foreach my $name ( sort {$hash{$b} <=> $hash{$a}} keys %hash){
	printf "%-8s %s\n", $name, $hash{$name};
}
