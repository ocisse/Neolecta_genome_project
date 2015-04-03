#!/usr/bin/perl -w

=head1 NAME

renamePhylipFile.pl

=head1 USAGE

renamePhylipFile.pl < MSA in phylip format>

=head1 DESCRIPTION

The script takes an MSA in phylip format as input and outputs a file where the name of each sequence
is  replaced by 0001, 00002 (SPECIES| number).

The script also automatically outputs a file called "input.mapid.txt", which contain the original
name (just in case).

=head1 PURPOSE

RAxML and other software complain that two sequences have the same name in the MSA. This name truncation
often occurs in the MSA, because of the phylip format requirement.

=head1 AUTHOR

OC, ousmanecis@gmail.com

=cut
use Carp; 
use Data::Dumper; 
use feature 'say'; 
use Getopt::Long; 


my $debug = 0; 

GetOptions(
	'v|verbose' 	=> \$debug,
	'h|help' 	=> sub {exec('perldoc',$0),
				exit(0);}


);

my $phy = shift; 
my %hash = (); 

# counter
my $i = 0;

open FILE,$phy || croak "cannot open $phy:$!\n"; 
while(<FILE>){
chomp; 
	
	if(/\|/){
		my @line = split /\|/,$_;
		   my ($spe,$id) = ($line[0],$line[1]); 
		
		my @line1 = split /\s/, $id;
		my ($tochange) = ($line1[0]);

		my @tmp = ($line1[1],$line1[2],$line1[3],$line1[4],$line1[5]);
		my $other = join (" ", @tmp);	
		 
		if ( $i < 10){

			# store in hash
			$hash{$tochange} = "0000$i"; 
			print "$spe".'|0000'.$i.' '.$other."\n";

			# increment
			$i++;
		 }
		elsif ( $i < 100){

			$hash{$tochange} = "000$i";
			print "$spe".'|000'.$i.' '.$other."\n";
			# increment
			$i++;
			}
		elsif ( $i < 1000){
			$hash{$tochange} = "00$i";
			 print "$spe".'|00'.$i.' '.$other."\n";
			 $i++;
			}	
		}

	else {
		say;
		}
	}

close FILE;

open (DTA,">$phy.mapid.txt") || croak "cannot write  on $phy.mapid.txt:$!\n";
say DTA Dumper (\%hash);

close DTA;
