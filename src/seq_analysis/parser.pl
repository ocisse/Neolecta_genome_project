#!/opt/linux/centos/7.x/x86_64/pkgs/perl/5.20.2/bin/perl -w

use strict;
use Data::Dumper;
use feature 'say';
use IO::All;

my @gp = ('ANID','Batde5','CALB','CCIN','LACBI','MBRE','NCRA','NIRR','PGRA','Rhior3','SCER','SCOMP','SPOM','SROS','TAPDE','TMEL','UMAY','MAGO7','PODAN','SORMK');

my $r = "";
print "gp";
foreach $r ( @gp ) {
	print ",$r";
}
print "\n";

my $o = io('all_orthomcl.out');
   $o->autoclose(0);
   while ( my $l = $o->getline || $o->getline ) {
   chomp $l;
	my @data = split /\s+/,$l;
#	say "$data[0]\t$data[3]";
	my $group = $data[0];
	my ($groupClean) = $group =~/(ORTHOMCL\d+)\(/;	

	shift @data;
	
	my $i = ""; 
	my @buff = ();

	foreach $i ( @data ) {
		my @tmp = split /\(/,$i;
		my ($spe,$id) = $tmp[0] =~/(\w+)\|(\.*)/;
		push(@buff,$spe);
	}
	
	my %compte = (); 
	my $el = ""; 
	foreach $el ( @buff ) {
		if (defined($el)) {
			$compte{$el}++;
		}
	}
	
#	say Dumper \%compte;

	print "$groupClean";
	
	my $x = "";
	foreach $x ( @gp ){
		if ($compte{$x}){
			print",$compte{$x}";
		} else {
			print",0";
		}
	}
	print "\n";	
}
