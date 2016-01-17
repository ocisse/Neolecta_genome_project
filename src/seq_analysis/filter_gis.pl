#!/opt/linux/centos/7.x/x86_64/pkgs/perl/5.20.2/bin/perl  -w 

use strict;
use feature 'say';
use Data::Dumper;
use IO::All;

analyze('gis/Aspergilus_gis.txt','ANID');
analyze('gis/Chytridiomycetes_gis.txt','Batde5');
analyze('gis/Saccharomycetes_gis.txt','CALB');
analyze('gis/Agaricomycetes_gis.txt','CCIN');
analyze('gis/Agaricomycetes_gis.txt','LACBI');
analyze('gis/Sordariomycetes_gis.txt','MAGO7');
analyze('gis/Choanoflagellida_gis.txt','MBRE');
analyze('gis/Sordariomycetes_gis.txt','NCRA');
analyze('gis/null_gis.txt','NIRR');
analyze('gis/Pyronemataceae_gis.txt','PCON');
analyze('gis/Pucciniomycotina_gis.txt','PGRA');
analyze('gis/Sordariomycetes_gis.txt','PODAN');
analyze('gis/Mucoromycotina_gis.txt','Rhior3');
analyze('gis/Saccharomycetes_gis.txt','SCER');
analyze('gis/null_gis.txt','SCOMP');
analyze('gis/Sordariomycetes_gis.txt','SORMK');
analyze('gis/Schizosaccharomycetes_gis.txt','SPOM');
analyze('gis/Pucciniomycotina_gis.txt','SROS');
analyze('gis/Taphrinomycetes_gis.txt','TAPDE');
analyze('gis/Tuberaceae_gis.txt','TMEL');
analyze('gis/Ustilaginomycotina_gis.txt','UMAY');

# sub
sub analyze {
	my ($gis,$name) = @_;
	my %gis = loadgis($gis);

 	my %unclustered = extract_uncl('not_clustered.txt',$name);
	say "$name\tUNCL\t:".scalar (keys %unclustered);

	#say Dumper \%unclustered;

 	my %rawUnclwithBlasthits = extract_blasthits(\%unclustered,$name); # key ANID , value = @ ' hits gis"
	say "$name\tRAW BLAST\t".scalar (keys %rawUnclwithBlasthits);

	my %filteredUncle = filter(\%rawUnclwithBlasthits,\%gis); # I remove those with gis from the same taxonimic linages
	say "$name\tFILTERED BLAST\t".scalar (keys %filteredUncle);

	# those have real blast hits
 	my %orphans = findOrphans(\%unclustered,\%filteredUncle);
	say "$name\tREAL ORPHANS\t".scalar (keys %orphans);
}

sub findOrphans {
	my ($uncl,$filter) = @_;
	my %uncl = %$uncl;
	my %filter = %$filter;

	my %orph = (); 
	my $u = "";
	foreach $u  (keys %uncl ) {
		unless ($filter{$u}){
			$orph{$u} = $u;
		}
	}
	return (%orph);
}
sub filter {
	my ($rawblast,$gis) = @_;
	
	my %rawblast = %$rawblast;
	my %gis = %$gis;

	my $r = "";
	foreach $r ( keys %rawblast ) {
		my @giNum = @{$rawblast{$r}};
		
		my $g = "";
		foreach $g ( @giNum ) {

			if (defined($gis{$g})){
				delete $rawblast{$r};
				last;
			}
		}	
	}
	return(%rawblast);

}
sub extract_blasthits {
	my ($hash,$nam) = @_; 

	say $nam;	
	my %hash = %$hash;

	my %blasthits = ();
	my $b = io('not_clustered.fa_vs_nr.tab');
	   $b->autoclose(0);
	   while ( my $bl = $b->getline || $b->getline ) {
	   chomp $bl; 
		next if $bl =~ m/^#/;
		if ( $bl =~ m/^$nam/) {
			my @bdat = split /\t/, $bl;
			if ( $hash{$bdat[0]}){
				#say "TEST\t$bdat[0]\t$bl";
				# OK it is an unclustered gene
				my ($gi) = $bdat[1] =~/gi\|(\d+)\|/;
				
				my @liste = ();
				if ($blasthits{$bdat[0]}) { # I have an entry already
					@liste = @{$blasthits{$bdat[0]}};
					push(@liste,$gi);
					@{$blasthits{$bdat[0]}} = @liste;
				} else {
					@liste = ();
					push(@liste,$gi);
					@{$blasthits{$bdat[0]}} = @liste;
				}
			}
		}		
	}
	return(%blasthits);
}

sub extract_uncl {
	my ($unclTtal,$nam) = @_;

	my %h1 = ();	
	my $f1 = io($unclTtal);
	   $f1->autoclose(0);
	   while ( my $fl1 = $f1->getline || $f1->getline ) {
	   chomp $fl1;
		if ( $fl1 =~ m/$nam/){
			my @dat = split /\t/, $fl1;
			$h1{$dat[1]} = 'x';
		}
	}
	return(%h1);
}
sub loadgis {
	my %h = (); 
	my $f = io(@_);
	   $f->autoclose(0);
	   while ( my $l = $f->getline || $f->getline ) {
	   chomp $l;
	    $h{$l} = $l;
	}
	return(%h);
}

