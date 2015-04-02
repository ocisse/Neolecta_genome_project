#!/opt/perl/5.16.3/bin/perl -w

use strict; 
use feature 'say'; 
use IO::All; 
use Carp; 


my $file = $ARGV[0]; 


my $fileSize = -s ($file);
#say "SIZE\t$fileSize";
if ($fileSize < '273'){
	say $file;
}

