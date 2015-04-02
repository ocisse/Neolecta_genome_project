#!/opt/perl/5.16.3/bin/perl -w 


use IO::All;
use feature 'say'; 
use Data::Dumper; 
use Carp; 

my $nirrInt = 0;
my $spomInt = 0;
my $saicoInt = 0;
my $pjirInt = 0;  
my $ncraInt = 0; 
my $ccinInt = 0; 
my $scerInt = 0; 
my $cneoInt = 0;

my %intPos = ();

my $genepainters = io("$ARGV[0]");
$genepainters->autoclose(0); 
while (my $line = $genepainters->getline || $genepainters->getline){
chomp $line; 
	
	# find all Neolecta genes that contains introns '|'
#	if ($line =~ m/NIRR_/){
 	if ($line =~ m/SPOM_/){
		$nirrInt = () = $line =~ /\-\|/g; 

		# record the position of intron
		my @neol = split /\s+/, $line; 
		
		my $i = ""; 
		my $count = 1; 
		foreach $i (@neol){
			if ($i =~ m/\-\|/){
				$intPos{"$i-$count"} = $count;
				$count++; 						
			} else {
				$count++;
			}
		}
	}
	
	# same for spom
	if ($line =~ m/SPOM_/){ $spomInt = () = $line =~ /\-\|/g;}
	if ($line =~ m/SCOMP_/){ $saicoInt = () = $line =~ /\-\|/g;}
	if ($line =~ m/PJIR_/){ $pjirInt = () = $line =~ /\-\|/g;}
	if ($line =~ m/NCRA_/){ $ncraInt = () = $line =~ /\-\|/g;}
	if ($line =~ m/CNEO_/) { $cneoInt = () = $line =~ /\-\|/g;}
	if ($line =~ m/CCIN_/){ $ccinInt = () = $line =~ /\-\|/g;}	
	if ($line =~ m/SCER_/) { $scerInt  = () = $line =~ /\-\|/g;}
}

say "ORGN\tNUM OF INTRONS";
say '-' x 80;
say "XXX\tNIRR\t$nirrInt";
say "XXX\tSPOM\t$spomInt";
say "XXX\tSAICO\t$saicoInt";
say "XXX\tPJIR\t$pjirInt";
say "XXX\tNCRA\t$ncraInt";
say "XXX\tCCIN\t$ccinInt";
say "XXX\tCNEO\t$cneoInt";
say "XXX\tSCER\t$scerInt";
 
say '-' x 80;

#say Dumper \%intPos;

# Now I want know how the introns are classiy

say "*****\tNEOLECTA INTRONS\t*****";
say '-' x 80;

say "Alignment\tIntron Pos\tOrigin";
say '-' x 80;

my $introns = ""; 
foreach $introns (keys %intPos){
	my $level = extractClass("$intPos{$introns}","$ARGV[0]"); 
	my $clean = $introns; 
           $clean =~s/\-\|\-//;
	say "#\t$genepainters\t$clean\t$level";
}

# sub 
sub extractClass {
	my ($pos,$file) = @_; 

	my $genepainterStat = io($file);
	$genepainterStat->autoclose(0);
	while ($line1 = $genepainterStat->getline || $genepainterStat->getline){
	chomp $line1;
		next if $line1 =~ m/^>/;
		next if $line1 =~ m/^Intron/; 
		
		my @dat = split /\t/, $line1;
		if ($pos == $dat[0]){
			my $lca = $dat[2];
			return $lca;
		}
	}
}












