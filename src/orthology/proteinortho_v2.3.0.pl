#!/usr/bin/perl
# based on proteinortho.pl version 2.2.0-beta
# tested with perl version v5.8.8 built for i386-linux-thread-multi and i686-linux-thread-multi
# author: Sonja J. Prohaska
# date: 26.11.2007

###
# VERAENDERUNGEN V 2.3.0 - marcus:
#
# Bugfix:
# ---------------------------------------------------
#  Format-DB legt keine psi und psq Files an, daher aus do_formatdb rausgenommen
#  mfa2pairs schiesst filehandles jetzt nach Verarbeitung (Sonst u.U. too many open filehandles)
#  $blastoutdir Bestimmung ueberarbeitet
#  per default working-directory
#  "/" am Ende von $blastoutdir wird ggf. entfernt
#
# Sonstiges:
# ---------------------------------------------------
#  -remove eingefuerht, die Blastoutputs loescht sobald nicht mehr benoetigt
#  Statt fastas in Befehlsaufruf zu uebergeben, kann auch Liste uebergeben werden
#  Prozentuale Fortschrittsanzeige mit Zeitvorhersage fuer fast alle Schritte
#  Ausgaben herausgekuerzt
#
# Threading:
# ---------------------------------------------------
#  -cpu=x legt nun Anzahl der Threads fest, blastall Parameter ist immer a=1
#  Programmablauf umstrukturiert, nach Aufgaben unterteielt
#  erst formatdb, dann $seq-Belegung, blasts, Header ausschreiben parallel zu mfa2pairs
#  formatdb und blasten gethreaded
#  
# Lastverteilung auf mehrere Rechner:
# ---------------------------------------------------
#  -blastonly zur Verteilung eingefuert, erledigt nur blasts
#  ermoeglicht blasten auf mehrere Rechner zu verteilen
#  dafuer dort einfach Script mit selben Parametern und Pfaden starten
#  Syncronisation der Verteielten Prozesse ueber Datei (netzwerkfaehig)
#  Prozesse koennen mit kill sauber beendet werden
#  Anlegen einer Datei stop im blastout-Verzeichniss beendet alle Prozesse sauber
#  Prozesse brechen sich bei Problemen (File nicht erreichbar) selbst ab
###


###
# VERAENDERUNGEN V 2.1.1 - Sven:
# ---------------
# * -dir Option um Output directory der *.bla files
# * neue file2path sub
# * avoid indexing:
#   `formatdb -i $faa -p T -o T`; -> `formatdb -i $faa -p T -o F`;
#   `formatdb -i $faa -p F -o T`; -> `formatdb -i $faa -p F -o F`;
###

#updates: version 2.0.1 --> option -a
#updates: version 2.1.0 --> option -r, -m, -f, -ff
#                           massive restructuring, adaptions all over the place
#                           new global variables:
#                             $blastoutpath, $force, $fforce, %seqs
#                           new subs:
#                             blastall
#                             do_formatdb
#                             formatdb
#                             mfa2pairs (replaced pair_po and NP)
#                             reciprocal_pairs
#                             _delete_pairs
#                             complete_pairs
#                             blastout2pairs
#                             _write2pairs
#                             print_pairs
#                             file2path_name
#                           formats databases only if necessary or forced
#                           writes blast output to files
#                           allows reading of blast output from files
#                           allows unidirectional (non-reciprocal) links (-r=0)
#                           allows more then one link per gene (-m)
#updates: version 2.1.1 --> bug: program crashes if -plog and/or -clog
#                                are not selected -- fixed

use strict;
use File::Basename;
 use threads;
 use threads::shared;
use vars qw/@query $cpu $blastoutpath $force $fforce $keep $blastonly %seqs/;

# catch SIGINT and SIGTERM for clean stop:
$SIG{'INT'}  = 'sighandler';
$SIG{'TERM'} = 'sighandler';

my $usage= << "JUS";
  usage: proteinortho_v2.1.pl [OPTIONS] MULTIFASTA1 MULTIFASTA2 [MULTIFASTA...]
  or:    proteinortho_v2.1.pl [OPTIONS] LIST_OF_MULTIFASTAS
  needs: blastall, formatdb
options: -e=        e-value [default: 1e-10]
         -a=        number of processors to use [default: 1]
         -p=        blast program {blastn|blastp} [default: blastp]
         -r=0|1     only *r*eciprocal blast hits [default: 1 (on)]
         -m=        maximum number of best blast hits [default: 1]
         -plog      write *p*airwise *b*last hits to logfile
         -plog=     change logfile-name [default: PB.LOG]
         -clog      write *c*onnected *c*omponents to logfile
         -clog=     change logfile-name [default: CC.LOG]
         -ulog      write ultimate log: postprocess plog and clog to logfile
         -ulog=     change logfile-name [default: ultimate.LOG]
         -f         force blastall (even if blast output is found)
         -ff        force formatdb (even if databases are found)
         -dir=      blast-output directory (default is the working directory)
         -remove    removes blastoutputs from blast-output directory after use
         -blastonly only do the blasts: capable to run on multiple machines
                    syncronized via file in blast-output directory, clean break
                    with 'touch stop' or SIGTERM (kill pid of the script)

         -r=1 -m=1 for pairwise reciprocal best blast
         -r=1 -m=x for fuzzy reciprocal blast (where x is degree of fuzzyness)

			 only one entry per in line in LIST_OF_MULTIFASTAS
                   DIFFERENT LOGFILES ARE RECOMMANDED!

-----------------------------------------------------------------------------
IMPORTANT! ensure that your protein IDs are unique ACROSS species! IMPORTANT!
IMPORTANT!   (eventually append species IDs to your protein IDs)   IMPORTANT!
-----------------------------------------------------------------------------
JUS

my @options = grep {substr($_, 0, 1) eq '-'} @ARGV;
my @others  = grep {substr($_, 0, 1) =~ /^[^-]/} @ARGV;

##default settings:
my $evalue = "1e-10";
my $blast = "blastp";
my $clog; #logfile
my $clogref;
my $plog; #logfile
my $plogref;
my $ulog; #logfile
my $ulogref;
my $maxl = 1;
my $reciprocal = 1;
my $bhold = 1;
our @joinable = (); # 4threading
share(@joinable);   # 4threading
$cpu = 1;           # 4threading
$keep = 1;
$force=0;
$fforce=0;
$blastonly = 0;
our $gotlock = 0;	# Stellt sicher das lock nicht von anderen weggenommen wird
##----------

##setting parameters:
foreach my $option (@options){
  if ($option =~ m/^-remove$/)     { $keep = 0;   }
  elsif ($option =~ m/^-blastonly$/)  { $blastonly = 1;   }
  elsif ($option =~ m/^-f$/)          { $force = 1;  }
  elsif ($option =~ m/^-ff$/)         { $fforce = 1; }
  elsif ($option =~ m/^-r=([01])$/)   { $reciprocal = $1; }
  elsif ($option =~ m/^-p=(.*)$/)     { $blast = $1; }
  elsif ($option =~ m/^-e=(.*)$/)     { $evalue= $1; }
  elsif ($option =~ m/^-a=(\d*)$/)    { $cpu   = $1; }
  elsif ($option =~ m/^-m=(\d*)$/)    { $maxl  = $1; }
  elsif ($option =~ m/^-plog$/)       { $plog  = "PB.LOG"; }
  elsif ($option =~ m/^-plog=(.*)$/)  { $plog  = $1; }
  elsif ($option =~ m/^-clog$/)       { $clog  = "CC.LOG"; }
  elsif ($option =~ m/^-clog=(.*)$/)  { $clog  = $1; }
  elsif ($option =~ m/^-ulog$/)       { $ulog  = "ultimate.LOG"; }
  elsif ($option =~ m/^-ulog=(.*)$/)  { $ulog  = $1; }
  elsif ($option =~ m/^-dir=(.*)$/)   { $blastoutpath  = $1; }
  else  { die "Invalid command line option: \'$option\'!\n$usage"; }
}


# A list can also be given
if (scalar(@others)==1) {
	open(LIST,"<$others[0]") || die ("At least two sequence files expected: MULTIFASTA1 MULTIFASTA2!\n$usage");
	@others = <LIST>;
	chomp(@others);
	close(LIST);
	foreach my $line (@others) {
		unless ($line =~ /.\../) {
			die ("$line \n does not seem to be a file!\n$usage");
		}
	}
}

die "At least two sequence files expected: MULTIFASTA1 MULTIFASTA2!\n$usage"
    unless scalar(@others)>=2;

if ($blastonly == 1 && $keep == 0) {
	die "Combination of remove and blastonly does not make much sense\n$usage";
}

print STDERR "Checking accessibility of files...";
for (my $i=0; $i<scalar(@others); $i++) {
  die "Can\'t read file \'$others[$i]\'!\n"
      unless open(FH, $others[$i]);
  close FH;
}
print STDERR "done\n";

#initializing filehandles
if ($plog) {
  die "Can\'t open/write to file \'$plog\'!\n"
      unless open (PLOG, ">$plog");
  $plogref = \*PLOG;
}
if ($clog) {
  die "Can\'t open/write to file \'$clog\'!\n"
      unless open (CLOG, ">$clog");
  $clogref = \*CLOG;
}
if ($ulog) {
  unless ($plog ne '' & $clog ne '') {
    die "you have to set \'-plog -clog\' in combination with \'-ulog\'!\n";
  }
  die "Can\'t open/write to file \'$ulog\'!\n"
      unless open (ULOG, ">$ulog");
  $ulogref = \*ULOG;
}

##begin collect pairwise blast hits for all species pairs:
#------------------------------------------------------------------------------
my @all;                    #square matrix (#species x #species),
                            #holds lists of reciprocal best blast hits
                            #for all species pairs

my $stat_runs = scalar(@others)*(scalar(@others)+1)/2;
my $stat_counter = 0;
my $glob_counter = 0;
my $stime;

#if (0) { #DEBUG
print STDERR "Formating databases...\n";

#format database (Marcus)
for (my $i=0; $i<scalar(@others); $i++) {
  if ($force == 1 | do_formatdb($others[$i],$blast) == 1) {
	# Wait for free Thread-Slots
	my @tliste = threads->list();
	while(defined($tliste[0]) && scalar@tliste >= $cpu*2) {
		sleep(1);
		if (scalar(@joinable) > 0) {
			lock(@joinable);

			# Fetch finished threads
			while(scalar@joinable>0) {
				my $tid = shift(@joinable);
				my $thr = threads->object($tid);
				$thr->join();
				@tliste = threads->list();
			}
		}
	}

	my $dbthread = threads->create('formatdb', ($others[$i],$blast));
  }
  else {
	print STDERR $others[$i]."\n-> formated database exists.\n";
  }
}

# Gather unfinished threads
foreach my $dbthread (threads->list()) {
	$dbthread->join();
}
@joinable = ();

print STDERR "Doing the blasts...\n";

# Blastoutpath belegen
unless (defined($blastoutpath)) {
	my $prepath=qx('pwd');
	chomp($prepath);
	$blastoutpath=$prepath;
}
$blastoutpath =~ s/\/$//; # Slash am Ende entfernen

# Syncronisation for blastonly (Marcus)
my $s_i = 0; # only needed for sync in blastonly case
my $s_j = 0;
if ($blastonly == 1) {
	unless (-e $blastoutpath."/sync") {
		print STDERR "Syncfile not found. Adding new one: ".$blastoutpath."/sync\n";
		&fileLock;
		open(SYNC,">$blastoutpath"."/sync") || die "sync: $!";
		print SYNC "0 0";
		close(SYNC);
		&fileUnlock;
		print STDERR "Sync creation done\n";
	}
	else {
		print STDERR "Syncfile found. Sycronizing...\n";
	}
	($s_i, $s_j) = &sync();
	$glob_counter = ($s_i-1*$s_i)/2+$s_j-$s_i;
}

# start the blastruns
$stime = time;
for (my $i=$s_i; $i<scalar(@others)-1; $i++) {
  for (my $j=$i+1; $j<scalar(@others); $j++) {
	$glob_counter++;

	if ($stat_counter%20 == 0) {
  		my $takentime = time - $stime;
		my $ttg = ($takentime/$glob_counter*$stat_runs-$takentime)/60;
		my $days = int($ttg/(60*24));
		my $hours = int(($ttg-60*24*$days)/60);
		my $minutes = int($ttg-60*24*$days-60*$hours);
		printf STDERR ("%.2f%% blasted overall -> ", ($glob_counter/$stat_runs*100));
		print STDERR "I think $days days $hours hours $minutes minutes to go\n";
	}


	if ($blastonly == 1) { # Syncronisation for blastonly
	  	if ($i < $s_i || ($i == $s_i && $j <= $s_j)) {next;}	 	 # skip what others should have

	    	($s_i, $s_j) = &sync($i,$j); 		 	 # Tell all I do a new job
	    	if ($s_i == -1) {$i = scalar(@others);last;} # Order received: stop
	    	unless ($i > $s_i || ($i == $s_i && $j > $s_j)) {next;}	 # skip what others have
	    	$stat_counter++;
		if ($stat_counter%20 == 0) {
			printf STDERR ("%.2f%% by this process \n", ($stat_counter/$stat_runs*100));
		}
	}

    ## Precheck if files already exsist, saves waiting for threadslots in this case
    if (&blast_needed($others[$i], $others[$j])) {

    # Wait for free Thread-Slots
    my @tliste = threads->list();
	while(defined($tliste[0]) && scalar@tliste >= $cpu) {
		sleep(1);
		if (scalar(@joinable) > 0) {
			lock(@joinable);
			# Fetch finished threads
			while(scalar@joinable>0) {
				my $tid = shift(@joinable);
				my $thr = threads->object($tid);
				$thr->join();
				@tliste = threads->list();
			}
		}
	}
	my ($blastthread) = threads->create('blast', ($others[$i], $others[$j], $blast, $evalue, $bhold, $maxl, $plogref, $reciprocal));
    } ## if (&blast_needed)
    else {
	    print STDERR "$others[$i] vs. \n$others[$j] \n-> blastall output exists.\n\n";
	    print STDERR "$others[$j] vs. \n$others[$i] \n-> blastall output exists.\n\n";
    }
  }

}

 # Gather unfinished threads
foreach my $thr (threads->list()) {
	$thr->join();
}
@joinable = ();

  print STDERR "Blastoutput generated\n";
  if ($s_i == -1) { print STDERR "\nCalculation terminated due to existance of stop-file\n";}

if ($blastonly == 1) {
  exit;
}
#} # Debug if 0

# Seqs belegen
for (my $j=0; $j<scalar(@others); $j++) {
	$seqs{$others[$j]}=$j;
}


##print output header
print STDERR "Printing output header...\n";
my $h_thread = threads->create('print_output_header', (\@others, *STDOUT));
# $h_thread->detach();
#	print_output_header(\@others, *STDOUT);

 
print STDERR "Loading pairwise blast data...\n";

### Gather information
#for all species pairs
# threading is contraproductie here
$stat_counter = 0;
$stime = time;
for (my $i=0; $i<scalar(@others)-1; $i++) {
  for (my $j=$i+1; $j<scalar(@others); $j++) {
     
    #find reciprocal best blast hits with pair_po
    #hash: protegin in species i (key) protein in species j (value) (href1)
    #hash: protein in species j (key) protein in species i (value) (href2)
    my ($t_i, $t_j, $pairs1, $pairs2) = mfa2pairs($i,$j,$others[$i], $others[$j], $blast, $evalue, $bhold, $maxl, $plogref, $reciprocal);
    @all->[$t_i][$t_j]=$pairs1;
    @all->[$t_j][$t_i]=$pairs2;
    $stat_counter++;
    if ($stat_counter%20 == 0) {
	my $takentime = time - $stime;
	my $ttg = ($takentime/$stat_counter*$stat_runs-$takentime)/60;
	my $days = int($ttg/(60*24));
	my $hours = int(($ttg-60*24*$days)/60);
	my $minutes = int($ttg-60*24*$days-60*$hours);
	printf STDERR ("%.2f%% loaded -> ", ($stat_counter/$stat_runs*100));
	print STDERR "I think $days days $hours hours $minutes minutes to go\n";
     }
  }
}

#------------------------------------------------------------------------------
##end collect pairwise blast hits for all species pairs;

#begin print debug
#------------------------------------------------------------------------------
#
#  print_all(\@all, *STDERR);
#
#------------------------------------------------------------------------------
#end print debug

close $plogref if $plogref;

# Join Headerthread (makes shure it has finished now)
$h_thread->join();

#==============================================================================
##main
#==============================================================================
##find valid sets of orthologs among connected components in graph G:
##nodes: proteins/genes/sequences (terms used synonymous)
##edges: [reciprocal] pairwise [best] blast hits
#------------------------------------------------------------------------------

if ($clog) {
  print STDERR "Writing connected components to \'$clog\'...";
}
else {
  print STDERR "Working hard... \n";
}

$stat_counter = 0;
$stime = time;
# go through upper triangle matrix (all pairwise species comparisons)
for (my $i=0; $i<scalar(@others)-1; $i++) {
  for (my $j=$i+1; $j<scalar(@others); $j++) {
    foreach my $gen_i (keys(%{@all->[$i][$j]})) {
	if (exists(@all->[$i][$j]->{$gen_i})) {
	#get the connected component with gene i as an element
	my $ccref=connected_component($gen_i, $i, \@all);
	if ($clog) {
	  print_connected_component($ccref, $clogref);
	}
	#the connected component is a valid set of orthologs
	#if it has not more then one gene per species
	#orthoset exchanges keys and values
	#(key=species, value=gene)
	my $orthoset=orthoset($ccref);
	if ($orthoset) {
	  for (my $k=0; $k<scalar(@others); $k++) {
	    if ($orthoset->{$k} eq "") {
	      print "*\t";	#FOUTPUT
	    }
	    else {
		print $orthoset->{$k},"\t"; #FOUTPUT
	    }
	  }
	  print "\n";		#FOUTPUT
	}
	#delete connected component from @all
	delete_cc_from_all($ccref,\@all);
      }
    }
		$stat_counter++;
		unless ($stat_counter%20) {
	  		my $takentime = time - $stime;
			my $ttg = ($takentime/$stat_counter*$stat_runs-$takentime)/60;
			my $days = int($ttg/(60*24));
			my $hours = int(($ttg-60*24*$days)/60);
			my $minutes = int($ttg-60*24*$days-60*$hours);
			printf STDERR ("%.2f%% processed -> ", ($stat_counter/$stat_runs*100));
			print STDERR "I think finally $days days $hours hours $minutes minutes to go\n";
		}
   }
}
close $clogref if $clogref;
#------------------------------------------------------------------------------

print STDERR "Writing sets of mono-orthologous genes to STDOUT... done.\n";


#begin compute ultimate log
#------------------------------------------------------------------------------
if ($ulog) {
  print STDERR "Writing to \'$ulog\'...";
  open CLOG, "<$clog" or die "Can\'t read/open \'$clog\'!\n";
  my $k;
  my $k_wc;
  my $out;
  while (<CLOG>) {
    if ($_ =~ m/^\#\s*(\d+)/) {
      $k_wc=$k=$1;
      next;
    }
    elsif ($_ =~ m/^\s*$/) {
      next;
    }
    else {
      $k-=1;
      my @line=split;
      $out .= `grep "$line[0]" $plog`;
      if ($k == 0) {
	open TMP, ">tmp" or die  "Can\'t open/write to \'tmp\'!\n";
	print TMP $out;
	close TMP;
	my $sort=`sort tmp | uniq`;
	my $l +=`sort tmp | uniq | wc -l`;
	print $ulogref "# $l $k_wc\n$sort";
	$out='';
      }
    }
  }
  close $ulogref;
  print STDERR "done.\n";
}
#------------------------------------------------------------------------------
#end compute ultimate log

#------------------------------------------------------------------------------

# Sync what to fetch on with blastonly-switch, handles locking
sub sync {
    my $i = shift;	# Was thread machen will
    my $j = shift;

	if (-e $blastoutpath."/stop") {
		print STDERR "Recived order to terminate procedure...\n Finishing running threads\n";
		return(-1,-1);
	}

	&fileLock;
	open(SYNC,"+<","$blastoutpath"."/sync") || die "sync: $!";
	my @content = <SYNC>;
	chomp(@content);
	(my $s_i, my $s_j) = split(/ /,$content[0],2);
	if (defined($i) && defined($j)) {
		if ($i > $s_i || ($i == $s_i && $j > $s_j)) {
			seek(SYNC,0,0);
			truncate(SYNC, 0);
			print SYNC "$i $j";
			print STDERR "Taking Job No $i x $j\n";
		}
	}
	close(SYNC);
	&fileUnlock;
	return ($s_i, $s_j);
}

#  takes: list of filenames (species) and filehandle
#returns: nothing
#writes to filehandle
sub print_output_header {	#FOUTPUT return;
  my ($filesref,$out)=@_;
  # (1) print species names (column heads for table)
  print $out "# " . join("\t", @$filesref) . "\n";
  # (2) get and print number of sequences in the species file
  print $out "# ";
  foreach my $i (@$filesref) {
    my $wc = `grep ">" $i | wc`;
    my @line=split(" ",$wc);
    print $out "$line[0]\t";
  }
  print $out "\n";
}

sub print_all {
  my ($allref, $log)=@_;
  my $x=scalar(@$allref);
  for (my $i=0; $i<$x; $i++) {
    for (my $j=0; $j<$x; $j++) {
      if ($i==$j) {next};
      my $ref= ref($allref->[$i][$j]);
#      print $log "# i:$i j:$j " . $ref . " ";
      if ($ref eq "HASH") {
	print $log "# " . scalar(keys(%{$allref->[$i][$j]})) ." $i $j\n";
	foreach my $k (keys(%{$allref->[$i][$j]})) {
	  foreach my $v (keys(%{$allref->[$i][$j]->{$k}})) {
	    print $log "$k $i $v $j\n";
	  }
	}
      }
      else { print $log "\n" }
    }
  }
}

#    takes: gene, gene's species and matrix (@all)
#  returns: hash of first order neighboring genes (with species as value)
#called by: connected_component
sub find_neighbors {
  my ($gene, $sp, $allref) = @_;
  my %neighbors;
  # gene can only have *more than one* neighbor in each species
  for (my $i=0; $i<scalar(@$allref); $i++) {
    if (exists($allref->[$sp][$i]->{$gene})) {
      foreach my $n (keys(%{$allref->[$sp][$i]->{$gene}})) {
	$neighbors{$n}=$i;
      }
    }
  }
  return \%neighbors;
}

#    takes: gene, gene's species and matrix (@all)
#  returns: hash of connected genes (with species as value)
sub connected_component {
  my ($gene, $sp, $allref) = @_;

  my %cc;    #holds genes (key) and species (value)
             #belonging to the contected component with 'gene'
             #already searched for their neighbors
  my %stack; #temporarely holds genes (key) [and species (value)]
             #for which neighbors need to be examined

  #load stack with first element
  $stack{$gene}=$sp;

  #start recursion:
  do {
    #take *a* gene off the stack
    #(this is neither depth first nor breadth first search, and doesn't matter)
    my @tmp=keys(%stack);
    my $g=shift(@tmp);               #gene
    my $s=$stack{$g};                #species
    #find the gene's neighbors
    my $nref=find_neighbors($g, $s, $allref);
    #add the neigboring genes not yet search for their neigbors to the stack
    foreach my $k (keys(%$nref)) {
      #elements of cc were already searched for neighbors
      unless (exists($cc{$k})) {
	#if the gene is on the stack
	if (exists($stack{$k})) {
	  #check uniqueness of gene IDs
	  unless ($stack{$k} eq $nref->{$k}) { # eq --> == for these are ints now
	    die "Error in gene ID: gene ID \'$k\' in more then one species!\n";
	  }
	}
	#if the gene is not on the stack yet, add it
	else {
	  $stack{$k}=$nref->{$k};
	}
      }
    }
    #move current gene to cc and remove it from the stack
    $cc{$g}=$s;
    delete($stack{$g});
  }
  until scalar(keys(%stack)) <= 0;
  return \%cc;
}

#  takes:connected component (%cc) and ref to filehandle
#returns: nothing
#writes to filehandle
sub print_connected_component {
  my ($ccref, $log) = @_;
  print $log "# " . scalar(keys(%$ccref)) . "\n";
  foreach my $k (sort(keys(%$ccref))) {
    print $log "$k " . $ccref->{$k} . "\n";
  }
}

#CAUTION! sub modifies existing datastructure @all
#  takes: connected component (%cc) and @all
#returns: @all: cc removed
sub delete_cc_from_all {
  my ($ccref, $allref) = @_;
  foreach my $gene (keys(%$ccref)) {
    my $sp=$ccref->{$gene};
    for (my $i=0; $i<scalar(@$allref); $i++) {
      delete($allref->[$sp][$i]->{$gene});
      delete($allref->[$i][$sp]->{$gene});
    }
  }
  return $allref;
}

#  takes: connected component (%cc)
#returns: (1) connected component with keys and values exchanged
#             if cc is a valid orthoset (no species has more then one gene)
#             or
#         (2) 0
#             otherwise
sub orthoset {
  my $ccref=shift;
  my %fercc;
  #exchange keys and values
  foreach my $k (keys(%$ccref)) {
    $fercc{$ccref->{$k}}=$k;
  }
  #if connected component is a valid orthoset
  #no species is represented with more the one gene =>
  if (scalar(keys(%$ccref)) == scalar(keys(%fercc))) {
    return \%fercc;
  }
  #if connected component is NOT a valid orthoset
  else {
    return 0;
  }
}

#writes to harddrive
sub blastall {
  my ($query, $db, $blast, $evalue, $out)=@_;
#  print STDERR "blastall -d $db -i $query -p $blast -a $cpu -m8 -e $evalue >  $out\n";
  `blastall -d $db -i $query -p $blast -a 1 -m8 -e $evalue >  $out`;
}

sub do_formatdb {
  my($faa,$blast)=@_;
#	print STDERR "File: $faa\n";
  my $val=1;
#  if ($blast eq "blastn" | $blast eq "tblastn") {
  if ($blast eq "blastn") {
    if (-e "$faa.nhr" &
	-e "$faa.nin" &
#	-e "$faa.nsd" &
#	-e "$faa.nsi" &
	-e "$faa.nsq") {
      $val=0;
    }
  }
#  elsif ($blast eq "blastp" | $blast eq "blastx") {
  elsif ($blast eq "blastp") {
    if (-e "$faa.phr" &
	-e "$faa.pin" &
#	-e "$faa.psd" &
#	-e "$faa.psi" &
	-e "$faa.psq") {
      $val=0;
    }
  }
  else {
    die "$blast is not supported yet! Feel free to adapt the code to your needs. The programmer.\n"
  }
  return $val;
}

sub formatdb {
  my($faa,$blast)=@_;
##(1) format multi-fasta files with formatdb
#------------------------------------------------------------------------------
  if ($blast =~ m/^blastp$/i) {
    print STDERR "Formating database $faa...\n";
    `formatdb -i $faa -p T -o F`;
  #  print STDERR "done\n";
  }
  elsif ($blast =~ m/^blastn$/i) {
    print STDERR "Formating database $faa...\n";
    `formatdb -i $faa -p F -o F`;
   # print STDERR "done.\n";
  }
  else {
    die "$blast is not supported yet! Feel free to adapt the code to your needs. The programmer.\n"
  }
#------------------------------------------------------------------------------

  # Thread-Ende melden (Marcus)
  lock(@joinable);
  push(@joinable,threads->tid());
  return;
}

# Extracted Path-Check from sub blast
sub blast_needed {
  my ($faa1, $faa2)=@_;
  my ($path1,$name1)=file2path_name($faa1);
  my ($path2,$name2)=file2path_name($faa2);

  my $blastout1=$blastoutpath."/".$name1.".".$name2.".bla";
  my $blastout2=$blastoutpath."/".$name2.".".$name1.".bla";

  if ($force == 0 && -e "$blastout1" && -e "$blastout2") {
	return 0;
  }
  return 1;
}

sub blast {
  my ($faa1, $faa2, $blast, $evalue, $bhold, $maxl, $plogref, $reciprocal)=@_;
  my ($path1,$name1)=file2path_name($faa1);
  my ($path2,$name2)=file2path_name($faa2);
 
  #do blastall: faa1 as query, faa2 as database
  my $blastout1=$blastoutpath."/".$name1.".".$name2.".bla";
  print STDERR "$faa1 vs.\n$faa2:\ngenerating blastall output...\n";
  blastall($faa1,$faa2,$blast,$evalue,$blastout1);
#  print STDERR "done\n";

  #do blastall: faa2 as query, faa1 as database
  my $blastout2=$blastoutpath."/".$name2.".".$name1.".bla";
  print STDERR "$faa2 vs.\n$faa1:\ngenerating blastall output...\n";
  blastall($faa2,$faa1,$blast,$evalue,$blastout2);
#  print STDERR "done\n";

  # Thread-Ende melden (Marcus)
  lock(@joinable);
  push(@joinable,threads->tid());

  return;
}

sub mfa2pairs {
  my ($t_i, $t_j, $faa1, $faa2, $blast, $evalue, $bhold, $maxl, $plogref, $reciprocal)=@_;
  my ($path1,$name1)=file2path_name($faa1);
  my ($path2,$name2)=file2path_name($faa2);

  # Get blastout location
  my $blastout1=$blastoutpath."/".$name1.".".$name2.".bla";
  my $blastout2=$blastoutpath."/".$name2.".".$name1.".bla";

  my $pairs1=blastout2pairs($blastout1, $evalue, $bhold, $maxl);

  #process blast output 1 & 2
  my $pairs2=blastout2pairs($blastout2, $evalue, $bhold, $maxl);

  # remove Blast-Output
  if ($keep == 0) {
	unlink($blastout1,$blastout2);
  }

  #delete directed edges (non-reciprocal hits) if this is what the user wants
  if ($reciprocal) {
#	print STDERR "Rec_pairs\n";    ##
	reciprocal_pairs($pairs1,$pairs2);
  }

  #print pairs
  #print STDERR ">>$plogref<<\n";

  if ($plogref) {
    print_pairs($pairs1,$faa1,$faa2,$plogref);
    print_pairs($pairs2,$faa2,$faa1,$plogref);
  }

  #convert directed to undirected links (a programming requirement)
  unless ($reciprocal) {
#    print STDERR complete_pairs"Comp Pairs\n";    ##
    complete_pairs($pairs1,$pairs2);
  }

	### Remove uneeded values
#	foreach my $key (keys %{$pairs1}) {
#		foreach my $subkey (keys %{$pairs1->{$key}}) {
#			# $pairs->{$line->[0]}->{$line->[1]}=[$line->[10],$line->[11],$line->[2]];
#			$pairs1->{$key}->{$subkey} = undef;
#		}
#	}

	### Remove uneeded values
#	foreach my $key (keys %{$pairs2}) {
#		foreach my $subkey (keys %{$pairs2->{$key}}) {
#			# $pairs->{$line->[0]}->{$line->[1]}=[$line->[10],$line->[11],$line->[2]];
#			$pairs2->{$key}->{$subkey} = undef;
#		}
#	}

  # Thread-Ende melden, nicht noetig (Marcus)
   return($t_i,$t_j,$pairs1,$pairs2);
}

#CAUTION! this subroutine modifies %pairs!
#called by mfa2pairs
sub reciprocal_pairs {
  my ($pairs1,$pairs2)=@_;

  foreach my $i1 (keys(%$pairs1)) {
    foreach my $j1 (keys(%{$pairs1->{$i1}})) {
      #if link is unidirectional
      unless (exists($pairs2->{$j1}) && exists($pairs2->{$j1}->{$i1})) {
	#delete link and nodes loosing all outgoing links
	_delete_pairs($pairs1,$i1,$j1);
      }
    }
  }
  foreach my $i2 (keys(%$pairs2)) {
    foreach my $j2 (keys(%{$pairs2->{$i2}})) {
      unless (exists($pairs1->{$j2}) && exists($pairs1->{$j2}->{$i2})) {
	#delete link and nodes loosing all outgoing links
	_delete_pairs($pairs2,$i2,$j2);
      }
    }
  }
}

#CAUTION! this subroutine modifies %pairs!
#CAUTION! this subroutine is privat (first char: '_')!
#called by reciprocal_pairs
sub _delete_pairs {
  my ($pairs,$ikey,$jkey)=@_;
  #delete link
#  print STDERR "first: deleted pair $ikey $jkey\n";
  delete($pairs->{$ikey}->{$jkey});
  #if degree of the node is now zero
  if (scalar(keys(%{$pairs->{$ikey}}))<1) {
    #delete node
#    print STDERR "deleted node $ikey\n";
    delete($pairs->{$ikey});
  }
}

#CAUTION! this subroutine modifies %pairs!
#called by mfa2pairs
sub complete_pairs {
  my ($pairs1,$pairs2)=@_;
  foreach my $i1 (keys(%$pairs1)) {
    foreach my $j1 (keys(%{$pairs1->{$i1}})) {
      unless (exists($pairs2->{$j1}->{$i1})) {
	$pairs2->{$j1}->{$i1}=[$pairs1->{$i1}->{$j1}->[0],
			       $pairs1->{$i1}->{$j1}->[1],
			       $pairs1->{$i1}->{$j1}->[2]];
      }
    }
  }
  foreach my $i2 (keys(%$pairs2)) {
    foreach my $j2 (keys(%{$pairs2->{$i2}})) {
      unless (exists($pairs1->{$j2}->{$i2})) {
	$pairs1->{$j2}->{$i2}=[$pairs2->{$i2}->{$j2}->[0],
			       $pairs2->{$i2}->{$j2}->[1],
			       $pairs2->{$i2}->{$j2}->[2]];
      }
    }
  }
}

#  takes: blast output location (either file or reference to filehandle)
#         thresholds for e-value (ehold) and bit score (bhold)
#         maximum number of links (maxl)
#returns: reference to hash:
#         key = gene in query
#         value = hash: key(s) = gene(s) in target
#                       value(s) = [0]=e-val, [1]=bit score, [2]=%ID
sub blastout2pairs {
  my ($blastout, $ehold, $bhold, $maxl)=@_;
  my %pairs;

  #$blastout is either a reference to a filehandle or a filename
  #BEGIN: prepair to read from whatever it is:
  my $fhref;
  if (ref($blastout) eq "GLOB") { #there is a blast output in a filehandle
    $fhref=$blastout;
  }
  else {                          #there is a blast output in a file
    open($fhref,"<$blastout") or die "Can\'t open/read \'$blastout\'!\n";
  }
  #END: prepair

  #what to write to %pairs (return value)
	
  while (<$fhref>) {
	my @line=split;			#[0]=seq1_query(-i) [1]=seq2_db(-d)
                    			#[2]=%ID [10]=e-value [11]=bitscore

    #hash structure only allows one one occurence
    #per pair of gene_i(query) and gene_d(database)
    #therefore, how should multiple occurences be treated?
    if (exists(${$pairs{$line[0]}}{$line[1]})) {
      #chose best scoring entry:
      #best scoring in terms of bitscore
      #(because I am not sure that perl reads the scientific notation
      #of e-values right)
      if (%pairs->{$line[0]}->{$line[1]}->[2] >= $line[11]) {
	#better scoring line is already in %pairs
      }
      else {
	#write better scoring line to %pairs
	_write2pairs(\%pairs,\@line,$ehold,$bhold,$maxl);
      }
    }
    else {
      #write line to %pairs
      _write2pairs(\%pairs,\@line,$ehold,$bhold,$maxl)
    }
  }

   unless (ref($blastout) eq "GLOB") { #there is not a blast output in a filehandle
	close($fhref); # neu von Marcus um to many open files zu vermeiden
   }

  return \%pairs;
}

#CAUTION! this subroutine modifies %pairs!
#called by blastout2pairs
sub _write2pairs {
  my ($pairs, $line, $ehold, $bhold, $maxl) = @_;
# scientific notation of $ehold and $line->[10]: '>=' comparison does not work
# if ($line->[10]>=$ehold & $line->[11]>=$bhold &
  if ($line->[11]>=$bhold &
      (scalar(keys(%{$pairs->{$line->[0]}}))<$maxl | #check space for new entry
       exists($pairs->{$line->[0]}->{$line->[1]}))) {#overwrite existing entry
    $pairs->{$line->[0]}->{$line->[1]}=[$line->[10],$line->[11],$line->[2]];
    return "1";
  }
  else {
    return "0";
  }
}

sub print_pairs {
  my ($pairs,$i,$j,$to)=@_;
  my $iID=$seqs{$i};
  my $jID=$seqs{$j};
  my $nr_pairs;
  my $hold;
  foreach my $ikey (keys(%$pairs)) {
    foreach my $jkey (keys(%{$pairs->{$ikey}})) {
      $nr_pairs++;
      $hold .= sprintf "$ikey $iID $jkey $jID " . ${$pairs->{$ikey}->{$jkey}}[0] . " " .
      ${$pairs->{$ikey}->{$jkey}}[1] . "\n";
    }
  }
  print STDERR ">$nr_pairs< >$iID< >$jID<\n";
  print $to "# $nr_pairs $iID $jID\n";
  print $to $hold;
}

sub file2path_name {
  my $file =shift;
  my ($name, $path) = fileparse($file);
  return $path, $name;
}

sub sighandler {
    my $signal = shift;    # signal-nummer
    $SIG{'INT'}  = 'sighandler';   # reinstall sighandler
    $SIG{'TERM'} = 'sighandler';   # reinstall sighandler

    print STDERR "Signal: SIG$signal caught! Prepering for termination...\n Fetching active threads...";

    $cpu = 0; # Dont start new threads while waiting
    if ($gotlock == 1) {
    	&fileUnlock;
    }

    if (defined(threads->list())) {
	    foreach my $thread (threads->list()) {
			$thread->join();
	    }
    }

    print STDERR " finished.\nTerminating...\n\n";
    exit;
}

sub fileLock   {
#	print STDERR "Requesting lock on sync...\n";
	my $counter = 0;
	while (!mkdir($blastoutpath."/lock",0755)) {   # if there already is a lock
		sleep(1);     # sleep for 1 sec and try again
		$counter++;
		if ($counter > 30) {
			die("Termination due to locking problem. Directory $_[0]");
		}
     }
    $gotlock = 1;
	if ($counter > 0) {
	   	print STDERR "Got lock on sync after $counter trys\n";
	}
    return $counter;
}

sub fileUnlock {
	if ($gotlock == 0) {return;}	 # Only Remove if it's mine
	rmdir($blastoutpath."/lock");     # remove file lock dir
    $gotlock = 0;
#	print STDERR "Released lock on sync\n";
}
