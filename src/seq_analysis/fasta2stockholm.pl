#!/opt/perl/5.16.3/bin/perl


use Carp; 
use Data::Dumper; 
use feature 'say'; 


my $usage = "$0 < gapped FASTA alignment File(s)>\n";

my @argv; 

while (@ARGV){
	my $arg = shift; 
	if ($arg ~~ /^-/){
		if ($arg ~~ "-h"){ print $usage; exit;}
		else { die $usage}
	}
	else{
		push @argv,$arg;
	}
}
push @argv, "-" unless @argv;

# loop throught FASTA files

foreach my $fasta (@argv){
	# read fasta file
	my (%seq,@name,$name); 
	
	open FASTA,"<$fasta" || croak "cannot open $fasta:$!\n"; 
	while(<FASTA>){
		if  (/^\>.*\|(\S+)/){
			$name = $1;
			croak "Duplicate name:$name" if $seq{$name};
			push @name, $name;
		}else{
			if (/\S/ && !defined $name){
				carp "Ignoring: $_";
			}else{
				s/\s//g;
				$seq{$name} .= $_;

			}
		}

	}
	close FASTA; 

# check all seqs have the same length

	my $length;
	my $lname; 
	foreach my $name (@name){
		my $l = length $seq{$name};
		if (defined $length){
			croak "sequences not all the same length ($lname is length, $name is $l)" unless $length ~~ $l;
		}else{
			$length = length $seq{$name};
			$lname = $name; 
		}
	}
# print Stockholm output	
	say "# STOCKHOLM 1.0";
	foreach my $name (@name){
		print $name, " ", $seq{$name},"\n";
	}
	say "//";
}
