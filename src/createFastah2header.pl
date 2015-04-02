#!/opt/perl/5.16.3/bin/perl -w

use strict; 
use IO::All; 
use Carp; 
use strict;
use feature 'say'; 

my %map = (
        'ANID' => 'Aspergillus nidulans',
        'ATHA' => 'Arabidopsis thaliana',
        'CALB' => 'Candida albicans',
        'CCIN' => 'Coprinopsis cinerea',
        'CNEO' => 'Cryptococcus neoformans',
        'HSAP' => 'Homo sapiens',
        'DDIS' => 'Dictyostelium discoideum',
        'MBRE' => 'Monosiga brevicollis',
        'NCRA' => 'Neurospora crassa',
        'NIRR' => 'Neolecta irregularis',
        'PCON' => 'Pyronema confluens',
        'PGRA' => 'Puccinia graminis',
        'PJIR' => 'Pneumocystis jirovecii',
        'SCER' => 'Saccharomyces cerevisiae',
        'SCOMP' => 'Saitoella complicata',
        'SPOM' => 'Schizosaccharomyces pombe',
        'SROS' => 'Sporobolomyces roseus',
        'TDEF' => 'Taphrina deformans',
        'TMEL' => 'Tuber melanosporum',
        'UMAY' => 'Ustilago maydis',
        'YLIP' => 'Yarrowia lipolytica'
);

my %tags = (); 

my $msa = io("$ARGV[0]"); 
$msa->autoclose(0); 
while(my $line = $msa->getline || $msa->getline){
chomp $line; 
	if ( $line =~ m/^>/){
		my @tmp = split /\_/, $line;
		my ($species,$protid) = ($tmp[0],$tmp[1]);
		$species =~s/>//;
		$line =~s/>//;
		say "$line\:\"$map{$species}\"";			
	}

}
