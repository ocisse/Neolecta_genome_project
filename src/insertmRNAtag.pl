#!/opt/perl/5.16.3/bin/perl -w

use IO::All; 
use feature 'say'; 
use Data::Dumper; 
use List::Util qw(max min);

my @start = ();
my @end = (); 
my %names = (); 


my $file = io("$ARGV[0]"); 
   $file->autoclose(0);  
   while(my $line = $file->getline || $file->getline){
   chomp $line; 
	next if $line =~ m/^#/; 
		
		my ($start,$end) = $line =~/CDS\t(\d+)\t(\d+)\t/;
		my @tmp = split /\t/,$line;
		   my $id = $tmp[8];
	              $id =~s/ID=4;Query=//;
		push(@start, $start);
	        push(@end,$end);
                $names{name} = $id;		
}

my $mRNAstart = min(@start);
my $mRNAend = max(@end);

# report
say "##gff-version	3";
say "##$names{name}\tx";


my @buffer = ();
my $exp = ();

my $file1 = io("$ARGV[0]");
   $file1->autoclose(0);
   while(my $line1 = $file1->getline || $file1->getline){
   chomp $line1;

	#say $line1 if $line1 =~ m/^#/;
	
	if ($line1 !~ m/^#/){
		my @data = split /\t/, $line1;
		my ($id,$scipio,$match, $startCDS,$endCDS,$tmp,$strand,$tmp2,$ann) = @data[0..8];
	#	say "$id\t$scipio\tmRNA\t$mRNAstart\t$mRNAend\t$tmp\t$strand\t$tmp2\t$ann";
		   $exp{"$id\t$scipio\tmRNA\t$mRNAstart\t$mRNAend\t$tmp\t$strand\t.\tID=4;Query=x"} = 'x';
		push(@buffer,$line1);
	}
}

my @test = keys %exp;
say @test;

my $in = ""; 
foreach $in (@buffer){
	say $in;
}
