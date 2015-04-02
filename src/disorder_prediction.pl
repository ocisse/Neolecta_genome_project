#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $file = shift;

my $in = Bio::SeqIO->new(-format => 'fasta', -file => $file);


while( my $seq = $in->next_seq ) {
    my $str = $seq->seq;
    print join("\t", $seq->display_id, &fold_index($str)), "\n";
}



sub fold_index {		#as in Prilusky et al 2005
    my ($seq) = @_;
    $seq =~ s/X//g;
    return (2.785*avg_hydropathy($seq) - abs(avg_charge($seq)) - 1.151);
}

sub avg_charge {		#as in Uversky et al. 2000
    my ($seq) = @_;    
    #my %ch = ( "R"=>"1", "K"=>"1", "D"=>"-1", "E"=>"-1" ); #as in Uversky et al. 2000
    my $c=0;
    $c++ while ($seq =~ m/[RK]/g);
    $c-- while ($seq =~ m/[DE]/g);
    return ($c/length($seq));
}

sub avg_hydropathy { #as in Uversky et al. 2000 - scale between 0 and 1 by dividing by 9.0
    my ($seq) = @_;
    my $h=0;
    my %HP = (			#Kyte & Doolittle 1982
				"R"=>"0.0","K"=>"0.6","D"=>"1.0","B"=>"1.0","N"=>"1.0","S"=>"3.6","E"=>"1.0",
				"H"=>"1.3","Z"=>"1.0","Q"=>"1.0","T"=>"3.8","G"=>"4.1","X"=>"4.1","A"=>"6.3",
				"P"=>"2.9","V"=>"8.7","Y"=>"3.2","C"=>"7.0","M"=>"6.4","I"=>"9.0","L"=>"8.2",
				"W"=>"3.6","F"=>"7.2");

    $h+=($HP{$1}/9.0) while ($seq =~ m/([A-Z])/g);
    return ($h/length($seq));
}
