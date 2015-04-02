#!/opt/perl/5.16.3/bin/perl -w 


use Data::Dumper; 
use Carp; 
use strict; 
use feature 'say'; 

my $in = shift; 

open I,$in || croak "cannot opne $in:$!\n"; 
while(<I>){
chomp; 

	my @line = split /\t/;
	
	if (/^FAMILYDESC/){
	say;
	my ($FAMILYDESC,$FAMILY,$ANID,$ATHA,$CALB,$CCIN,$CNEO,$DDIS,$HSAP,$MBRE,$NCRA,$NIRR,$PCON,$PGRA,$PJIR,$SCER,$SCOMP,$SCRY,$SJAP,$SOCT,$SPOM,$SROS,$TDEF,$TMEL,$UMAY,$YLIP) = @line[0..25];
	}
	else {
	my ($FAMILYDESC,$FAMILY,$ANID,$ATHA,$CALB,$CCIN,$CNEO,$DDIS,$HSAP,$MBRE,$NCRA,$NIRR,$PCON,$PGRA,$PJIR,$SCER,$SCOMP,$SCRY,$SJAP,$SOCT,$SPOM,$SROS,$TDEF,$TMEL,$UMAY,$YLIP) = @line[0..25];
	
		if ($ANID !~ 0 && $ATHA !~ 0 && $CALB !~ 0 && $CCIN !~ 0 && $CNEO !~ 0 && $DDIS !~ 0 && $HSAP !~ 0 && $MBRE !~ 0 && $NCRA !~ 0 && $NIRR !~ 0 && $PCON !~ 0 && $PGRA !~ 0 && $PJIR !~ 0 && $SCER !~ 0 && $SCOMP !~ 0 && $SPOM !~ 0 &&  $SROS !~ 0 && $TDEF !~ 0 && $TMEL !~ 0 && $UMAY !~ 0 && $YLIP !~ 0)
		{
			say;
		}

	}
}

close I;	
