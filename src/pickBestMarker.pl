#!/opt/perl/5.16.3/bin/perl -w

=head1 USAGE

pickBestMarker.pl -i <hmmscan output> -r (optional)

=head1 DESCRIPTION

this simple script expects the result of the following cmd

example: hmmscan -E 0.01 --tblout <outFile> <hmm database> <protein sequence>

=head1 OPTIONS

-r => output best markers e-values

=head1 AUTHORS

Ousmane Cisse, ousmane.cisse@ucr.edu

=cut

use Data::Dumper; 
use feature 'say'; 
use Carp; 
use strict; 
#use Method::Signatures;
use Getopt::Long; 

# data
my ($report,$hmmscan); 

# Options
GetOptions(
    'h|help' => 	sub {exec('perldoc',$0);
	exit(0);},
    'i|input' =>    \$hmmscan,
    'r|report' => 	\$report,
);

#
unless( $ARGV[0] ) {
    die("need an input file using -i or --input on the cmdline. See usage with -h\n");
} else {
    $hmmscan = $ARGV[0];
}

# parse the hmmscan ouput and store data 
my %AllMarkers = (); 

open H,$hmmscan || croak "cannot open $hmmscan:$!\n";
while(<H>){
    chomp; 
    if(/^#/){
	next; 
    }
    else{
	my @line = split /-/; 
	my ($maker,$prot) = @line[0..1];
	my ($value) = $_ =~/-\s+.+\s+\-\s+(.+)\s+/;
	my $eval = extractEvalue($value); 

	# just making a fake key : $maker.'-'.$prot
	$AllMarkers{"$maker-$prot"} = $eval;
    }
}
close H;

# select best markers based on e-values
my %BestMarkers = (); 
my %EvaluesBestMarkers = (); 

my ($i,$t) = ''; 
foreach $i (sort keys %AllMarkers){
    my $markerRec1 = retrieveOri($i);
    my ($t,$eval) = pickTheLowesEval($markerRec1,%AllMarkers); 

    # clean the marker name
    $markerRec1 =~s/\s+$//;
    $BestMarkers{$markerRec1} = $t;
    $EvaluesBestMarkers{$markerRec1}{$t} = $eval;	
}


# print report
report(%BestMarkers); 

# report e-value if requested
if ($report){
    say Dumper \%EvaluesBestMarkers;
}
else {
    say"\n#\tno e-value report requested - if you want to see e-values, please add -r tag";
    print "\n";
}

# sub
sub report {
    my (%bm) = @_; 
    my $b = ''; 

    foreach $b (keys %bm){
	if ($bm{$b}){
	    say "$b\t$bm{$b}";
	}
	else{
	    say "$b\tno marker found";
	}
    }
}

sub pickTheLowesEval {
    my ($item, %data) = @_; 

    my %h1 = (); 

    my $x = ''; 
    foreach $x ( keys %data){
	my $short = retrieveOri($x);

	if ( $item ~~ $short){
	    $h1{$x} = $data{$x};	
	}
    }

    # ordering by values
    my $key = ''; 
    my $lowest = "0.01"; 

    foreach $key (keys %h1){
	if ( $h1{$key} < $lowest){
	    $lowest = $h1{$key};
	}	
    }
    # retrieving the key which has the lowest e-value
    my ($key2,$cleanKey2) = ''; 

    foreach $key2 ( keys %h1){
	if ( $h1{$key2} ~~ $lowest){

	    # clean to retain only the sequence header
	    $cleanKey2 = clean($key2);	
	    return($cleanKey2,$lowest);
	    last;
	}
    }
}
sub clean {
    my ($toclean) = @_; 
    my $clean = ''; 
    my @tmp = split /-/, $toclean;

    $clean = $tmp[1];
    $clean =~s/^\s+//;
    $clean =~s/\s+$//;
    return($clean);
}

sub retrieveOri{
    my ($in) = @_; 
    my @name = split /-/, $in;
    return($name[0]);
}
sub extractEvalue {
    my ($data) = @_;
    my @data1 = split /\s+/,$data; 
    return($data1[0]);
}

