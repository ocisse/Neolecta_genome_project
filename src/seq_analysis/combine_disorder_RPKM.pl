#!/usr/bin/perl -w
use strict;
use Getopt::Long;
my ($disorder,$rpkm) = (
    'disorder/PASA_Neolecta1.gene_structures_post_PASA_updates.13923.disorder',
    'rnaseq/tophat/cufflinks/isoforms.fpkm_tracking');

GetOptions(
    'd|disorder:s' => \$disorder,
    'e|f|fpkm:s'   => \$rpkm,
    );

my %genes;
open(my $fh => $disorder) || die $!;
while(<$fh>) {
    my ($gene,$d) = split;
    $genes{$gene}->{disorder} = $d;
}

open($fh => $rpkm) || die $!;
my $header = <$fh>;
while(<$fh>) {
    my ($isoform,@rest) = split;
    if( $rest[-1] ne 'OK' ) {
	warn("Skipping $isoform, FPKM is not ok\n");
	next;
    }
    my $rpkm = $rest[-4];
    $genes{$isoform}->{rpkm} = $rpkm;
}

print join("\t", qw(GENE FPKM DISORDER)), "\n";
for my $gene ( sort { ($genes{$b}->{rpkm} || 0) <=> ($genes{$a}->{rpkm} || 0) } keys %genes ) {
    print join("\t", $gene, $genes{$gene}->{rpkm} || 0 ,$genes{$gene}->{disorder} || 'NA'),"\n";
}
