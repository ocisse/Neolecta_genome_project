#!/usr/bin/perl
# get_fasta_stats - Get statistics of contigs in Fasta format.
#
# Written by: James D. White, University of Oklahoma, Advanced Center
#   for Genome Technology
#
#Date Written: Jul 27, 1999
#
# 20001025 JDW - Added GC% and non-ACGT base calculations.
# 20011016 JDW - Set quality score for non-ACGT bases to zero.
# 20020318 JDW - Add longest run column to output.
# 20050927 JDW - Add -d, -n, and -x options to compute extended dna
#                statistics.  Use Getopt::Std to process options.
# 20060321 JDW - Expand search rules for finding a ".qual" file with
#                -q.  Modify -q -t headers (Non-ACGT -> NonACGT) to
#                make headers line up better with data.
# 20060405 JDW - Fix divide by zero when all bases are non-ACGT and
#                reduce output clutter.
# 20070904 JDW - Add -g and -T options for histogram and totals only.
#		 Add average sequence length to totals.
# 20071010 JDW - Add -4 option to print number of 454 Newbler
#		 assembled reads for each contig.
# 20080222 JDW - Add -a option to print contig assembly stats.
# 20080226 JDW - Fix numreads recognition for -4 option and add number
#		 of 454 reads to -a option stats.
# 20090312 JDW - Add -e flag to print minimal stats on empty input
#		 file, instead of ending with an error.  Add -r flag
#		 to call input sequences "reads" instead of "contigs"
#		 in the output messages.
# 20090325 JDW - Fix help to add info about -e flag.  Correct check
#		 of option conflicts with -a.  Get rid of divide by
#		 zero when -4 used with no 454 reads.
# 20090331 JDW - Print 'mode=None' when there is no mode.
# 20090624 JDW - Add Longest run of Ns and non-ACGT stats.
# 20090810 JDW - Fix bug in -a option stats.
#
$Date_Last_Modified = "August 10, 2009";


###################### Customization Parameters ######################

$NUMREADS_454 = 0;		# -4, print 454 numReads= field for
				#     newbler assembled contigs?
$ASSEMBLY_STATS = 0;		# -a, print assembly stats
$DOUBLE_STRANDED = 0;		# -d, computed extended statistics on
				#     both strands.  Implies $EXTENDED.
$EMPTY_OK = 0;			# -e  Produce minimal stats on empty
				#     input file instead of error?
$GRAPH_HISTOGRAM = 0;		# -g, print histogram of sequence
				#     sizes?
$NUM_BINS = 50;			# -G, number of bins for histogram
$PRINT_EXTENDED_N = 0;		# -n, use non-ACGT bases as Ns in
				#     extended statistics.
$PHRED_SCORE = 20;		# -p #, phred/phrap score considered
				#     "good"?  (-p is used only if -q
				#     is specified)
$QUAL = 0;			# -q, use fasta qual file for quality
				#     info?
$CONTIG = 'contig';		# -r, Input sequences are to be called
				#     'contigs' or 'reads'
$SHORTEN_CONTIG_NAMES = 0;	# -s, shorten contig names?
$TOTALS = 0;			# -t, (=1) output totals and headers?
				# -T, (=2) output totals and headers
				#     without individual sequence
				#     stats?
$EXTENDED = 0;			# -x, computed extended statistics

$COUNT_WIDTH = 60;		# maximum number of count symbols per
				#     line of histogram output

########## Operating System specific ##########

$directory_separator = '/';	# For Unix

################## End of Customization Parameters ##################


$full_command_name = $0;
if ($full_command_name =~ m"(^|${directory_separator})([^${directory_separator}]*)$")
  {
  $my_name = $2;
  }
else
  {
  $my_name = 'get_fasta_stats';
  }

# get command line options
use Getopt::Std;
($::opt_a, $::opt_d, $::opt_e, $::opt_g, $::opt_G, $::opt_h,
 $::opt_n, $::opt_p, $::opt_q, $::opt_r, $::opt_s, $::opt_t,
 $::opt_T, $::opt_x, $::opt_4) = (undef) x 15;
display_help('') unless( getopts('adegG:hnp:qrstTx4') );
display_more_help() if ($::opt_h);
$NUMREADS_454 =  1		if (defined $::opt_4);
$ASSEMBLY_STATS = 1		if (defined $::opt_a);
$EXTENDED = $DOUBLE_STRANDED = 1 if (defined $::opt_d);
$EMPTY_OK = 1			if (defined $::opt_e);
$GRAPH_HISTOGRAM = 1		if (defined $::opt_g);
if (defined $::opt_G)
  {
  $GRAPH_HISTOGRAM = 1;
  $NUM_BINS = $::opt_G;
  }
$PRINT_EXTENDED_N = 1 		if (defined $::opt_n);
$PHRED_SCORE = $::opt_p		if (defined $::opt_p);
$QUAL = 1			if (defined $::opt_q);
$CONTIG = 'read'		if (defined $::opt_r);
$SHORTEN_CONTIG_NAMES = 1	if (defined $::opt_s);
$TOTALS = 1			if (defined $::opt_t);
$TOTALS = 2			if (defined $::opt_T);
$EXTENDED = 1			if (defined $::opt_x);

# check options for sanity
display_help('Cannot use -q or -4 with -d or -x')
  if (($QUAL || $NUMREADS_454) && $EXTENDED);
display_help("Invalid Phred_score: -p $PHRED_SCORE")
  if ($QUAL && $PHRED_SCORE !~ /^\d\d?$/);
display_help("Invalid value for -G option='$NUM_BINS'\n$USAGE")
  unless $NUM_BINS =~ /\A\d+\Z/ && $NUM_BINS > 2;
display_help('Cannot use -a with -g, -G, or -r')
  if ($ASSEMBLY_STATS && ($GRAPH_HISTOGRAM || $CONTIG eq 'read'));
$GRAPH_HISTOGRAM = 2 if ($ASSEMBLY_STATS);

$fasta_input_file = shift @ARGV;
$fasta_input_file = '-' if (!$fasta_input_file);
open(FASTAIN, $fasta_input_file) || die("Can't open fasta_input_file: '$fasta_input_file'\n");

$fasta_qual_file = "$fasta_input_file.qual";
if ($QUAL)
  {
  if ($fasta_input_file ne '-')
    {
    # assume fasta_input_file.qual, e.g., xxx + xxx.qual,
    #   xxx.fa + xxx.fa.qual, xxx.fna + xxx.fna.qual, or
    #   xxx.fasta + xxx.fasta.qual
    unless (-f $fasta_qual_file)
      {
      $fasta_qual_file = $fasta_input_file;
      # Just adding .qual didn't work, so now try removing the last
      #   qualifier, e.g., xxx.fa + xxx.qual, xxx.fna + xxx.qual, or
      #   xxx.fasta + xxx.qual
      unless ($fasta_qual_file =~ s/\.f(ast|n)?a$/.qual/ &&
              -f $fasta_qual_file)
        {
        print STDERR "Can't find a fasta quality file for '$fasta_input_file'.\n  Quality reporting turned off\n";
        $QUAL = 0;
        } 
      }
    if ($QUAL && !open(QUALIN, $fasta_qual_file))
      {
      print STDERR "Can't open fasta quality input file: '$fasta_qual_file'.\n  Quality reporting turned off\n";
      $QUAL = 0;
      }
    }
  else
    {
    print STDERR "Can't open fasta quality input file for standard input.\n  Quality reporting turned off\n";
    $QUAL = 0;
    }
  }


$Num_Contigs = 0;
$Total_Bases = 0;
$Total_ACGT = 0;
$Total_GC = 0;
$Max_Bases = -1;
$Min_Bases = -1;
$Total_Reads = 0;
$Max_Ns = 0;
$Total_N_Ends = 0;
$Max_Nons = 0;
$Total_Non_ACGT_Ends = 0;
if ($QUAL)
  {
  $Min_Qual = -1;
  $Max_Qual = -1;
  $Phred_Bases = 0;
  $Error_Sum = 0;
  $Error_Bases = 0;
  $ErrorsP10K = 0;
  $Edited_Bases = 0;
  $Longest_Run = 0;
  @QualBins = ();
  for ($q = 99; $q >= 0; $q--)
    {
    $QualBins[$q] = 0;
    }
  }
if ($EXTENDED)
  {
  # initialize extended dna counters
  @mono = qw(a c g t);
  @di = map { ( "${_}a", "${_}c", "${_}g", "${_}t" ) } @mono;
  @tri = ( (map { "a$_" } @di), (map { "c$_" } @di),
           (map { "g$_" } @di), (map { "t$_" } @di) ); 
  @x_list = (@mono, 'n', 'total1', @di, 'nn', 'total2',
             @tri, 'nnn', 'total3');
  @x_print = (
    (grep { $_ le revcomp($_) } @mono), 'n',   'total1',
    (grep { $_ le revcomp($_) } @di  ), 'nn',  'total2',
    (grep { $_ le revcomp($_) } @tri ), 'nnn', 'total3' );
  # reset counters
  %x_totals = %x_counters = map { ($_, 0) } @x_list;
  }
my(@CONTIG_LENGTHS, @NUM_READS,			# for -g, -G, and -a
   %READS_BY_LENGTH);

$found_data = -1;
$header = '';
$contig = '';
$sequence = '';
$line_num = 0;
if ($QUAL)
  {
  $Qline = <QUALIN>;
  $Qline_num = 1;
  }
while ($line = <FASTAIN>)
  {
  chomp $line;
  $line_num++;
  if ($line =~ /^>/)
    {
    if ($found_data < 0)	# first input line?
      {
      print_headers(0) if $TOTALS;
      }
    else			# else process previously read contig
      {
      if ($QUAL && length($sequence) != @Quality)
        {
        print STDERR "Fasta sequence and quality file lengths do not match on\n  $CONTIG='$contig'\n";
        exit 2;
        }
      process_contig($contig, $sequence, $num_reads) if ($contig);
      }
    $found_data = 0;		# Now set up new contig
    $Num_Contigs++;
    $header = $line;
    if ($header =~ m/^>(\S+)/)
      {
      $contig = $1;
      $num_reads = ($header =~ m/\snum[Rr]eads=(\d+)/) ? $1 : 0;
      }
    else
      {
      print STDERR "Error: Invalid fasta_input_file format: '$fasta_input_file'\n  Fasta input line number=$line_num\n";
      exit 2;
      }
    $sequence = '';
    if ($QUAL)
      {
      chomp($Qline);
      $Qheader = $Qline;
      if ($Qheader =~ m/^>(\S+)/)
        {
        $Qcontig = $1;
        if ($contig ne $Qcontig)
          {
          print STDERR "Fasta sequence and quality files do not match on $CONTIG header number $Num_Contigs\n  Sequence $CONTIG='$contig', Quality $CONTIG='$Qcontig'\n";
          exit 2;
          }
        }
      else
        {
        print STDERR "Error: Invalid fasta_qual_input_file format: '${fasta_input_file}.qual'\nQline_num=$Qline_num, Qline='$Qline'\n";
        exit 2;
        }
      $Qline = <QUALIN>;
      $Qline_num++;
      $Quality = '';
      while (length($Qline) && $Qline !~ /^>/)
        {
        chomp($Qline);
        $Quality .= ' ' . $Qline;
        $Qline = <QUALIN>;
        $Qline_num++;
        }
      $Quality =~ s/^\s+//;
      $Quality =~ s/\s+$//;
      @Quality = split(' ', $Quality);
      }
    }
  else # ($line !~ /^>/)
    {
    if ($found_data < 0)
      {
      print STDERR "Error: Invalid fasta_input_file format: '$fasta_input_file'\n";
      exit 2;
      }
    $line =~ s/\s//g;
    $sequence .= $line;
    $found_data = 1;
    }
  } # end while
close(FASTAIN);
close(QUALIN) if $QUAL;

if ($contig)
  {
  process_contig($contig, $sequence, $num_reads);
  }
elsif ($EMPTY_OK)
  {
  print_headers(1);
  }
else
  {
  print STDERR "Error: Empty fasta_input_file: '$fasta_input_file'\n";
  exit 2;
  }

print_totals() if $TOTALS;
print_histogram() if ($GRAPH_HISTOGRAM && $Num_Contigs);
print STDOUT "\n";

exit 0;


######################################################################
# process_contig - compute contig statistics and accumulate overall
#   stats.
######################################################################

sub process_contig
  {
  my($contig, $sequence, $num_reads) = @_;
  my($len, $ACGTbases, $ATbases, $GCbases, $nonACGTbases);

  # Remove Contig name prefix?
  $contig =~ s/^.*([Cc]ontig)/$1/ if $SHORTEN_CONTIG_NAMES;
  $len = length($sequence);
  push @CONTIG_LENGTHS, $len;
  push @NUM_READS, $num_reads;
  $Total_Bases += $len;
  $Max_Bases = $len if $Max_Bases < $len;
  $Min_Bases = $len if $Min_Bases > $len || $Min_Bases < 0;

  $ATbases = ($sequence =~ tr/aAtT/aAtT/);
  $GCbases = ($sequence =~ tr/cCgG/cCgG/);
  $ACGTbases = $ATbases + $GCbases;
  $nonACGTbases = $len - $ACGTbases;
  if ($ACGTbases)
    {
    $GC_per_cent = sprintf "%.1f%%", 100 * $GCbases / $ACGTbases;
    }
  else
    {
    $GC_per_cent = '-';
    }
  $Total_GC += $GCbases;
  $Total_ACGT += $ACGTbases;
  $Total_Reads += $num_reads;
  if ($nonACGTbases)
    {
    my $more_Max_Nons = $Max_Nons + 1;
    my @Nons = ($sequence =~ /[^acgtACGT]{$more_Max_Nons,}/g);
    foreach (@Nons)
      {
      my $l = length $_;
      $Max_Nons = $l if ($Max_Nons < $l);
      }
    $Total_Non_ACGT_Ends += length $1 if ($sequence =~ /^([^acgtACGT]+)/);
    if (substr($sequence, -1) =~ /[^acgtACGT]+$/)
      {
      my $rs = reverse $sequence;
      $Total_Non_ACGT_Ends += length $1 if ($rs =~ /^([^acgtACGT]+)[acgtACGT]/);
      }
    my $more_Max_Ns = $Max_Ns + 1;
    my @Ns = ($sequence =~ /[nN]{$more_Max_Ns,}/g);
    foreach (@Ns)
      {
      my $l = length $_;
      $Max_Ns = $l if ($Max_Ns < $l);
      }
    $Total_N_Ends += length $1 if ($sequence =~ /^([nN]+)/);
    if (uc substr($sequence, -1) eq 'N' && uc substr($sequence, 0, 1) ne 'N')
      {
      my $rs = uc reverse $sequence;
      $Total_N_Ends += length $1 if ($rs =~ /^(N+)/);
      }
    }

  if ($QUAL)
    {
    if ($nonACGTbases)	# set qualities for non-ACGT bases to zero
      {
      $sequence =~ s/[^acgtACGT]/N/g;	# set non-ACGT to 'N' for index
      my $found = 0;
      while (($found = index($sequence, 'N', $found)) >= 0)
        {
        $Quality[$found++] = 0;
        }
      }
    my($min_qual, $max_qual, $phred_bases, $q, $error_sum,
       $error_bases, $errorsp10k, $edited_bases, $current_run,
       $longest_run);
    $longest_run = $current_run = 0;
    my @qualbins = ();
    for ($q = 99; $q >= 0; $q--)
      {
      $qualbins[$q] = 0;
      }
    foreach $q (@Quality)
      {
      $qualbins[$q]++;
      if ($q >= $PHRED_SCORE)
        {
        $current_run++;
        $longest_run = $current_run if $longest_run < $current_run;
        }
      else
        {
        $current_run = 0;
        }
      }
    $min_qual = -1;
    $max_qual = -1;
    $phred_bases = 0;
    $error_sum = 0;
    $error_bases = 0;
    for ($q = 97; $q >= 0; $q--)
      {
      $qb = $qualbins[$q];
      if ($qb)
        {
        $min_qual = $q;
        $max_qual = $q if $max_qual < 0;
        $phred_bases += $qb if ($q >= $PHRED_SCORE);
        $error_sum += $qb * (10.0 ** ($q / -10.0));
        $error_bases += $qb;
        $QualBins[$q] += $qb;
        }
      }
    if ($error_bases > 0)
      {
      $errorsp10k = ($error_sum / $error_bases) * 10000;
      }
    else
      {
      $errorsp10k = 0;
      }
    if ($errorsp10k >= 1000)
      {
      $errorsp10k = int($errorsp10k + 0.5);
      }
    else
      {
      $errorsp10k = int(100.0 * $errorsp10k + 0.5) / 100.0;
      }
    $edited_bases = $qualbins[98] + $qualbins[99];
    $contig .= ' ' if length($contig) == 7;
    print STDOUT "$contig\t$len\t$errorsp10k\t$edited_bases\t" .
                 "$phred_bases\t$min_qual\t$max_qual\t$GC_per_cent\t" .
		 "$nonACGTbases\t$longest_run" .
    		 ($NUMREADS_454 ? "\t$num_reads" : '') . "\n"
      unless ($TOTALS > 1);
    $Min_Qual = $min_qual if $Min_Qual > $min_qual || $Min_Qual < 0;
    $Max_Qual = $max_qual if $Max_Qual < $max_qual;
    $Phred_Bases += $phred_bases;
    $Longest_Run = $longest_run if $Longest_Run < $longest_run;
    $QualBins[98] += $qualbins[98];
    $QualBins[99] += $qualbins[99];
    }
  elsif ($EXTENDED)
    {
    $sequence =~ s/[^ACGTacgtn]/n/g;
    $sequence =~ tr/ACGT/acgt/;
    process_extended_contig($contig, $sequence, 1 - $DOUBLE_STRANDED);
    process_extended_contig($contig, revcomp($sequence), 1)
      if ($DOUBLE_STRANDED);	# process reverse strand?
    }
  else
    {
    print STDOUT "$contig\t$len\t$GC_per_cent\t$nonACGTbases" .
    		 ($NUMREADS_454 ? "\t$num_reads" : '') . "\n"
      unless ($TOTALS > 1);
    }
  } # end process_contig


######################################################################
# process_extended_contig - compute extended contig statistics and
#   accumulate overall extended stats.
######################################################################

sub process_extended_contig
  {
  my($contig, $sequence, $print) = @_;
  my $len = length($sequence);

  # count mono-nucleotides and total
  $x_counters{'total1'} += $len;
  my $newa = ($sequence =~ tr/a/a/);
  $x_counters{'a'} += $newa;
  my $newc = ($sequence =~ tr/c/c/);
  $x_counters{'c'} += $newc;
  my $newg = ($sequence =~ tr/g/g/);
  $x_counters{'g'} += $newg;
  my $newt = ($sequence =~ tr/t/t/);
  $x_counters{'t'} += $newt;
  $x_counters{'n'} += $len - ($newa + $newc + $newg + $newt);

  # count di- and tri-nucleotides
  for ($i = 0; $i < $len - 2; $i++)
    {
    print STDERR "$contig ... $i\n"
      if (($i % 500000) == 0) && ($i != 0);
    $x_counters{'total2'}++;
    $x_counters{'total3'}++;
    my $tri = substr($sequence, $i, 3);
    my $di = substr($tri, 0, 2);
    if ($di !~ /n/)
      {
      $x_counters{$di}++;
      }
    else
      {
      $x_counters{'nn'}++;
      }
    if ($tri !~ /n/)
      {
      $x_counters{$tri}++;
      }
    else
      {
      $x_counters{'nnn'}++;
      }
    } # end for ($i = 0; $i < $len - 2; $i++)
  if ($len > 1)		# take care of last di-nucleotide
    {
    $x_counters{'total2'}++;
    my $last2 = substr($sequence, -2, 2);
    if ($last2 !~ /n/)
      {
      $x_counters{$last2}++;
      }
    else
      {
      $x_counters{'nn'}++;
      }
    }

  # remove non-ACGT from counts unless requested
  unless ($PRINT_EXTENDED_N)
    {
    $x_counters{'total1'} -= $x_counters{'n'};
#    $x_counters{'n'} = 0;
    $x_counters{'total2'} -= $x_counters{'nn'};
    $x_counters{'nn'} = 0;
    $x_counters{'total3'} -= $x_counters{'nnn'};
    $x_counters{'nnn'} = 0;
    }

  if ($print)
    {
    print STDOUT "\n$CONTIG=$contig\n\n" unless ($TOTALS > 1);
    print_extended(%x_counters) unless ($TOTALS > 1);
    # update counters
    $x_totals{$_} += $x_counters{$_} foreach keys %x_counters;
    %x_counters = map { ($_, 0) } keys %x_counters;
    }
  } # end process_extended_contig


######################################################################
# print_extended - print extended dna counters.
######################################################################

sub print_extended
  {
  my (%counters) = @_;
  my $total_total = $counters{'total1'} + $counters{'n'};
  unless ($total_total)
    {
    printf "\u$CONTIG length = 0\n";
    return;
    }
  printf "Non-ACGT bases = %ld = %3.2f%% of total\n", $counters{'n'},
    100.0 * $counters{'n'} / $total_total unless ($PRINT_EXTENDED_N);
  print STDOUT "\nSeq         Count Percent     Seq         Count Percent            Sum    Sum%\n\n"
    if $TOTALS;
  unless ($counters{'total1'})
    {
    print "ACGT bases     = 0 = 0% of total\n";
    return;
    }
  if ($PRINT_EXTENDED_N && ($counters{'total1'} == $counters{'n'}))
    {
    printf STDOUT "%-6s %10ld %6.2f%%\n", 'n', $counters{'n'}, 100.0;
    return;
    }
  foreach my $seq (@x_print)
    {
    $seqn = $counters{$seq};
    if ($seq =~ /total/)
      {
      printf STDOUT "\n%6s %10ld\n\n",
        $seq, $seqn;
      next;
      }
    my $nlen = length $seq;
    my $total = $counters{'total' . $nlen};
    if ($seq =~ /n/)
      {
      if ($PRINT_EXTENDED_N)
        {
        if ($total)
          {
          printf STDOUT "\n%-6s %10ld %6.2f%%\n",
            $seq, $seqn, 100.0 * $seqn / $total;
          }
        else
          {
          printf STDOUT "\n%-6s %10ld    -\n", $seq, $seqn;
          }
        }
      next;
      }
    my $rev = revcomp($seq);
    if ($seq eq $rev)
      {
      if ($total)
        {
        printf STDOUT "%-6s %10ld %6.2f%%\n",
          $seq, $seqn, 100.0 * $seqn / $total;
        }
      else
        {
        printf STDOUT "%-6s %10ld    -\n", $seq, $seqn;
        }
      }
    else 
      {
      $revn = $counters{$rev};
      if ($total)
        {
        printf STDOUT "%-6s %10ld %6.2f%%     %-6s %10ld %6.2f%%     %10ld %6.2f%%\n",
          $seq, $seqn, 100.0 * $seqn / $total,
          $rev, $revn, 100.0 * $revn / $total,
          $seqn + $revn, 100.0 * ($seqn + $revn) / $total;
        }
      else
        {
        printf STDOUT "%-6s %10ld    -        %-6s %10ld    -        %10ld    -\n",
          $seq, $seqn, $rev, $revn, $seqn + $revn;
        }
      }
    } # end foreach my $seq (@x_print)
  } # end print_extended


######################################################################
# print_headers($no_contigs) - print output headers.
######################################################################

sub print_headers
  {
  my($no_contigs) = @_;
  print STDOUT "\n$my_name - Last Modified: $Date_Last_Modified\n\n";
  print STDOUT "\u$CONTIG statistics for fasta file:'$fasta_input_file'\n";
  return if ($TOTALS > 1 || $no_contigs);

  if ($QUAL)
    {
    print STDOUT "      and for fasta quality file:'$fasta_qual_file'\n";
    print STDOUT "\n\u$CONTIG     \t\u$CONTIG\tConsed\tEdited\tPhred$PHRED_SCORE\tMinimum\tMaximum\t\tNonACGT\tLongest" .
    		 ($NUMREADS_454 ? "\tNumber" : '') . "\n";
    print STDOUT "Name    \tLength\tErr/10K\tBases\tBases\tQuality\tQuality\tGC%\tbases\tPhred$PHRED_SCORE" .
    		 ($NUMREADS_454 ? "\tReads" : '') . "\n";
    }
  elsif (! $EXTENDED)
    {
    print STDOUT "\n\u$CONTIG    \t\u$CONTIG\t\tNonACGT" .
    		 ($NUMREADS_454 ? "\tNumber" : '') . "\n";
    print STDOUT "Name    \tLength\tGC%\tbases" .
    		 ($NUMREADS_454 ? "\tReads" : '') . "\n";
    }
  } # end print_headers


######################################################################
# print_totals - print total statistics.
######################################################################

sub print_totals
  {
  my($ave_len);
  unless ($Num_Contigs)
    {
    print STDOUT "\nNumber of \u${CONTIG}s=0, Total bp=0, Shortest=0, Longest=0,\n";
    print STDOUT "Average length=0, Average GC%=NA, Non-ACGT bases=0\n";
    return;
    }

  print STDOUT "\nNumber of \u${CONTIG}s=$Num_Contigs, Total bp=$Total_Bases, Shortest=$Min_Bases, Longest=$Max_Bases,\n";
  $ave_len = sprintf "%.1f", $Total_Bases / $Num_Contigs;
  if ($Total_ACGT)
    {
    $GC_per_cent = sprintf "%.1f%%", 100 * $Total_GC / $Total_ACGT;
    }
  else
    {
    $GC_per_cent = '-';
    }
  $Total_nonACGT = $Total_Bases - $Total_ACGT;
  print STDOUT "Average length=$ave_len, Average GC%=$GC_per_cent," .
	       " Non-ACGT bases=$Total_nonACGT";
  print STDOUT ",\nLongest Run of non-ACGT Bases=$Max_Nons, Total non-ACGT bases on ${CONTIG} ends=$Total_Non_ACGT_Ends,\n" .
	       "Longest Run of Ns=$Max_Ns, Total Ns on ${CONTIG} ends=$Total_N_Ends"
    if ($Total_nonACGT);

  if ($QUAL)
    {
    my($q, $qb);
    for ($q = 97; $q >= 0; $q--)
      {
      $qb = $QualBins[$q];
      if ($qb)
        {
        $Error_Sum += $qb * (10.0 ** ($q / -10.0));
        $Error_Bases += $qb;
        }
      }
    if ($Error_Bases > 0)
      {
      $ErrorsP10K = ($Error_Sum / $Error_Bases) * 10000;
      }
    else
      {
      $ErrorsP10K = 0;
      }
    if ($ErrorsP10K >= 1000)
      {
      $ErrorsP10K = int($ErrorsP10K + 0.5);
      }
    else
      {
      $ErrorsP10K = int(100.0 * $ErrorsP10K + 0.5) / 100.0;
      }
    $Edited_Bases = $QualBins[98] + $QualBins[99];
    print STDOUT ",\nErrors/10K=$ErrorsP10K, Edited bp=$Edited_Bases, Total Phred/Phrap$PHRED_SCORE bp=$Phred_Bases,\n";
    print STDOUT "Minimum \u$CONTIG Quality=$Min_Qual, Maximum \u$CONTIG Quality=$Max_Qual,\n";
    print STDOUT "Longest Run of Phred/Phrap$PHRED_SCORE bp=$Longest_Run";
    }
  elsif ($EXTENDED)
    {
    print STDOUT "\n\nTotals for file='$fasta_input_file'\n\n";
    print_extended(%x_totals);
    }
  print STDOUT ",\nTotal Newbler Assembled 454 Reads=$Total_Reads"
    if ($NUMREADS_454 && ! $EXTENDED);
  print STDOUT "\n";
  } # end print_totals


######################################################################
# print_histogram - print a histogram of contig lengths
######################################################################
sub print_histogram
  {
  my($range, $median, $mode, $mode_max, $sum, $len, $real_divisor,
     $divisor, $min_bases, $max_bases, $factor, $bin,
     $count, $max_count, $count_divisor, $real_count_divisor,
     $real_width, $width, $reads, $i, 
     @modes, @sorted_lengths, @unique_lengths,
     %length_counts, %scale_counts, %bins);
  $range = $Max_Bases - $Min_Bases;
  return unless $range;	# don't waste time on trivial cases
  return unless $Num_Contigs > 2;

# accumulate counts and 454 reads by contig length
  for ($i = 0; $i < scalar @CONTIG_LENGTHS; $i++)
    {
    $len = $CONTIG_LENGTHS[$i];
    $reads = $NUM_READS[$i];
    $length_counts{$len}++;
#    $READS_BY_LENGTH{$len} = 0 unless defined $READS_BY_LENGTH{$len};
    $READS_BY_LENGTH{$len} += $reads;
    }

  # compute the median
  @sorted_lengths = sort { $a <=> $b } @CONTIG_LENGTHS;
  if ($Num_Contigs & 1)
    {
    $median = $sorted_lengths[$Num_Contigs >> 1];
    }
  else
    {
    $sum = $sorted_lengths[$Num_Contigs >> 1] +
      $sorted_lengths[($Num_Contigs >> 1) - 1];
    $median = ($sum & 1) ? sprintf("%.1f", $sum / 2) : $sum >> 1;
    }

  # compute the mode(s)
  $mode_max = -1;
  @unique_lengths = sort { $a <=> $b } keys %length_counts;
  foreach $len (@unique_lengths)
    {
    if ($length_counts{$len} > $mode_max)
      {
      @modes = ($len);
      $mode_max = $length_counts{$len};
      }
    elsif ($length_counts{$len} == $mode_max)
      {
      push @modes, $len;
      }
    }
  $mode = (@modes > 1) ? sprintf("(%s)", join(", ", @modes)) : $modes[0];

  printf STDOUT "Median=%d, Mode=%s\n\n", $median,
    ($mode_max <= 1 && @modes > 1) ? 'None' :
      "$mode occurs $mode_max times";

  # go elsewhere to print assembly stats
  if ($GRAPH_HISTOGRAM > 1)
    {
    print_assembly_stats(%length_counts);
    return;
    }

  print STDOUT "Histogram of sequence lengths\n\n";

  # rescale lengths by K or M if necessary
  $factor = '';
  $scale = 1;
  if ($Max_Bases > 1000 * $NUM_BINS)	# scale down lengths?
    {
    if ($Max_Bases > 1000000 * $NUM_BINS) # scale by 1M
      {
      $factor =  'M ';
      $scale = 1000000;
      }
    else # ($Max_Bases > 1000 * $NUM_BINS) # scale by 1K
      {
      $factor =  'K ';
      $scale = 1000;
      }
    }
  %scale_counts = ();	# collapse countss into rescaled bins
  $scale_counts{int($_ / $scale)} += $length_counts{$_}
    foreach keys %length_counts;
#printf STDOUT "min_bases=%d, max_bases=%d, divisor=%d, %d scale_counts=('%s')\n\n", $min_bases, $max_bases, $divisor, scalar keys %scale_counts, join("', '", map { "$_=$scale_counts{$_}" } (sort { $a <=> $b } keys %scale_counts));

  # put length_counts into bins
  if (scalar @unique_lengths > $NUM_BINS)
    {
    $real_divisor = $range / $scale / $NUM_BINS;
    $divisor = int($real_divisor);
    $divisor++ if (($real_divisor - $divisor) > 0);
    $min_bases = $divisor * int(($Min_Bases - 0.000001) / $divisor);
    $max_bases = $divisor * int(($Max_Bases + $divisor - 1) / $divisor);
    foreach ($bin = $min_bases; $bin <= $max_bases; $bin += $divisor)
      {
      $bins{$bin} = 0;
      }
    foreach $length (sort { $a <=> $b } keys %scale_counts)
      {
      $bin = $divisor * int(($length + $divisor - 1) / $divisor);
      $bins{$bin} += $scale_counts{$length};
      }
    }
  else # we already have few enough bins, so just copy
    {
    %bins = %scale_counts;
    }

  # delete zero count bins on the ends
  foreach $bin (sort { $a <=> $b } keys %bins)
    {
    last if ($bins{$bin} > 0);
    delete $bins{$bin};
    }
  foreach $bin (sort { $b <=> $a } keys %bins)
    {
    last if ($bins{$bin} > 0);
    delete $bins{$bin};
    }
#printf STDOUT "min_bases=%d, max_bases=%d, divisor=%d, %d bins=('%s')\n\n", $min_bases, $max_bases, $divisor, scalar keys %bins, join("', '", map { "$_=$bins{$_}" } (sort { $a <=> $b } keys %bins));

  # now compute the width scale
  $max_count = 0;
  foreach $count (values %bins)
    {
    $max_count = $count if $count > $max_count;
    }
  if ($max_count <= $COUNT_WIDTH)
    {
    $count_divisor =  1;
    }
  else
    {
    $real_count_divisor = $max_count / $COUNT_WIDTH;
    $count_divisor = int($real_count_divisor);
    $count_divisor++ if $real_count_divisor > $count_divisor;
    }
#print STDOUT "max_count=$max_count, real_max_count=$real_max_count, count_divisor=$count_divisor\n";
    
  # finally print the histogram
  printf STDOUT "Maximum\t\t (Each '*' represents %d occurrences,\n",
    $count_divisor;
  printf STDOUT "%sBases\tCount\t  '+' represents 1/2 '*')\n",
    $factor;
  foreach $bin (sort { $a <=> $b } keys %bins)
    {
    $real_width = $bins{$bin} / $count_divisor;
    $width = int($real_width);
    printf STDOUT "%d\t%d\t|%s%s\n", $bin, $bins{$bin}, '*' x $width,
      ($real_width - $width > 0.499999) ? '+' : '';
    }
  } # end print_histogram


######################################################################
# print_assembly_stats - print a histogram of contig lengths
######################################################################
sub print_assembly_stats
  {
  my(%length_counts) = @_;

  my($length, $bin, $count, $bases, $reads, $total_contigs,
     $total_bases, $total_reads);
  my @bin_maxes = (999, 1999, 2999, 3999, 4999, 9999, 19999, 29999,
    39999, 49999, 99999, 1e50);
  my @bin_titles = ('< 1 kb  ', '1 - 2 kb', '2 - 3 kb', '3 - 4 kb',
    '4 - 5 kb', '5 - 10 kb', '10 - 20 kb', '20 - 30 kb', '30 - 40 kb',
    '40 - 50 kb', '50 - 100 kb', '>= 100 kb');
  my @bin_counts = (0) x scalar @bin_titles;
  my @bin_bases = (0) x scalar @bin_titles;
  my @bin_reads = (0) x scalar @bin_titles;

  print STDOUT "Assembled contig length summary\n\n";

  # put length_counts into bins
  $total_bases = 0;
  $total_contigs = 0;
  $total_reads = 0;
  $bin = 0;
  foreach $length (sort { $a <=> $b } keys %length_counts)
    {
    $bin++ while ($length > $bin_maxes[$bin]);
    $count = $length_counts{$length};
    $bin_counts[$bin] += $count;
    $total_contigs += $count;
    $bases = $length * $count;
    $bin_bases[$bin] += $bases;
    $total_bases += $bases;
    $reads = $READS_BY_LENGTH{$length};
    $bin_reads[$bin] += $reads;
    $total_reads += $reads;
    }
    
  # finally print the assembly stats
  print STDOUT "        \t  Total\t        \t% of";
  print STDOUT "\t        \t% of" if ($NUMREADS_454);
  print STDOUT "\nContig Size\t\u${CONTIG}s\tTotal Bases\tBases";
  print STDOUT "\tTotal Reads\tReads" if ($NUMREADS_454);
  for ($bin = 0; $bin < scalar @bin_titles; $bin++)
    {
    $count = $bin_counts[$bin];
    $bases = $bin_bases[$bin];
    $reads = $bin_reads[$bin];
    printf STDOUT "\n%s\t%7d\t%9d\t%.1f%%", $bin_titles[$bin], $count,
      $bases, 100.0 * $bases / $total_bases;
    printf STDOUT "\t%9d\t%.1f%%",
      $reads, 100.0 * $reads / $total_reads if ($NUMREADS_454 && $total_reads);
    }
  printf "\n\nCumulative\t%7d\t%9d\t%.1f%%",
    $total_contigs, $total_bases, 100.0;
  printf "\t%9d\t%.1f%%", $total_reads, ($total_reads ? 100.0 : 0)
    if ($NUMREADS_454);
  $count = $total_contigs - $bin_counts[0];
  $bases = $total_bases - $bin_bases[0];
  $reads = $total_reads - $bin_reads[0];
  printf "\nCumulative>1kb\t%7d\t%9d\t%.1f%%",
    $count, $bases, 100.0 * $bases / $total_bases;
  printf "\t%9d\t%.1f%%", $reads, 100.0 * $reads / $total_reads
    if ($NUMREADS_454 && $total_reads);
  $count -= $bin_counts[1];
  $bases -= $bin_bases[1];
  $reads -= $bin_reads[1];
  printf "\nCumulative>2kb\t%7d\t%9d\t%.1f%%",
    $count, $bases, 100.0 * $bases / $total_bases;
  printf "\t%9d\t%.1f%%", $reads, 100.0 * $reads / $total_reads
    if ($NUMREADS_454 && $total_reads);
  print STDOUT "\n\n";
  } # end print_assembly_stats


######################################################################
# revcomp - compute reverse complement of a dna sequence
######################################################################

sub revcomp
  {
  my($sequence) = @_;
  $sequence =~ tr/acgt/tgca/;
  return reverse $sequence;
  } # end revcomp


######################################################################
# display_help
######################################################################

sub display_help
  {
  my($msg) = @_;
  print STDERR "\n$msg\n" if $msg;
  print STDERR <<EOF;

USAGE: $my_name [-a] [-d] [-e] [-g] [-G num_bars] [-n]
          [-q [-p phred_score]] [-r] [-s] [-t] [-T] [-x] [-4]
	  [fasta_input_file]
              or
       $my_name -h           <== For more information


EOF
  exit 2;
  } # end display_help


######################################################################
# display_help
######################################################################

sub display_more_help
  {
  print STDOUT <<EOF;

$my_name - Get statistics for dna sequences in Fasta format.
The contig names optionally may be shortened by removing everything
before the word "Contig" (-s).  Summary statistics (totals) may be
displayed (-t or -T).  A fasta quality file can also be read to give
error and quality statistics (-q).  A minimum Phred/Phrap score
(-p phred_score) can be specified if quality scores are read.  Base
pairs with a Phred/Phrap quality value of 98 or 99 are counted as
edited bases, but are no used to compute other quality statistics.

If quality scores are not used, only the contig name, length, GC%, and
the number of bases that are not A, C, G, or T are listed for each
contig.  The GC% uses only A, C, G, and T bases in the calculation.
If totals are requested (-t), the number of contigs, total length of
all contigs, lengths of the shortest and longest contigs, average
contig length, average GC%, and total number of non-ACGT bases are
also output.  Column headings are displayed for the contigs, if totals
are requested.

If quality scores are used (-q), the Consed Errors/10Kb, the number
of edited bases (those bases with a quality score of 98 or 99), the
number of Phred/Phrap(#) bases, and the minimum and maximum
Phred/Phrap scores are listed for each contig along with the contig
name and length.  Phred/Phrap(#) bases is the count of the number of
bases with a Phred/Phrap base quality of (#) or greater.  This
Phred/Phrap score (#) defaults to ${PHRED_SCORE}, but a different value can be
specified as '-p new#'.

If both quality scores are used (-q) and totals are requested (-t),
the Consed Errors/10Kb, the number of Phred/Phrap(#) bases, and the
minimum and maximum Phred/Phrap scores are output for the entire file
along with the total contig count and size information listed above.

If extended statistics are requested (-d or -x), then mono-, di-, and
tri-nucleotide statistics are computed for each contig, and for the
entire input file if (-t) is also specified.

If a sequence length histogram is requested (-g), then it will follow
the rest of the output.


USAGE: $my_name [-a] [-d] [-e] [-g] [-G num_bars] [-n]
          [-q [-p phred_score]] [-r] [-s] [-t] [-T] [-x] [-4]
	  [fasta_input_file]
               or
       $my_name -h           <== What you are reading

  where 'num_bars' is the maximum number of frequency bars in the
            sequence length histogram.  Fewer bars may be printed
	    in order to produce bin boundaries with integer values.

	'phred_score' is a threshhold phred/Phrap score used to
	    indicate a minimum "good" quality score.

        'fasta_input_file' is the name of the input sequence file in
            Fasta format, which contains the contigs to be processed.
            \u$my_name will also read a Fasta quality file named:
            "'fasta_input_file'.qual".  If 'fasta_input_file' is
	    omitted, standard input will be used, but a quality file
	    may not be read.


OPTIONS:

  -a  Produce contig assembly stats. -a may not be used with -g, -G,
      or -r.  -a computes similar stats to -g and -G, but contig
      assembly stats are printed instead of a bar chart.

  -d  Compute extended dna statistics on both strands.  -d may not be
      used with -q or -4.  If both -d and -x are specified, then -d is
      used.

  -e  Produce minimal stats on an empty input file instead of an error.

  -g  Produce a histogram graph with up to $NUM_BINS bars of contig lengths
      frequencies for the entire file.  A histogram will not be
      produced if there are fewer then 3 contigs or if all contigs are
      of the same length.

  -G num_bars  Produce a histogram graph with up to 'num_bars' bars of
      contig length frequencies for the entire file.  A histogram will
      not be produced if there are fewer then 3 contigs or if all
      contigs are of the same length.

  -n  Include non-ACGT bases in extended statistics if -d or -x is
      specified;  otherwise this option is ignored.

  -p  Specify a threshhold phred/Phrap score.  For example, the
      default value is ${PHRED_SCORE}.  \u$my_name will then display a
      Phred/Phrap${PHRED_SCORE} score, which is a count of the number of bases
      with a Phred/Phrap score of ${PHRED_SCORE} or better.  '-p 30' would specify
      that a Phred/Phrap30 score should be computed instead, as a
      count of the number of bases with a base quality score of 30 or
      better.  -p is ignored if a Fasta quality file is not read.

  -q  Compute and output quality statistics, using a Fasta quality
      file named: "'fasta_input_file'.qual".  If that file does not
      exist, and 'fasta_input_file' ends in ".fa", ".fna", or
      ".fasta", then a second try is made by replacing the final
      ".fa", ".fna", or ".fasta" with ".qual".  If the Fasta quality
      file cannot be opened or the 'fasta_input_file' is read from
      standard input, then -q is ignored.  Neither -d nor -x may be
      used with -q.

  -r  Output messages are to refer to input sequences as "reads",
      instead of "contigs".  May not be used with -a.
      
  -s  The contig names will be shortened by removing any prefix before
      the word "Contig", i.e., "gono.fasta.screen.Contig26" becomes
      "Contig26".

  -t  Output total fasta file statistics, as well as individual contig
      statistics.  Column headings for the contigs are also displayed.

  -T  Output total fasta file statistics, but not individual contig
      statistics.

  -x  Compute extended dna statistics on the forward strands.  -x may
      not be used with -q or -4.  If both -d and -x are specified,
      then -d is used.

  -4  Print number of 454 Newbler assembled reads for each contig.  If
      the contigs do not have the 454 Newbler " numReads=" comment,
      the number printed will be zero.  -4 may not be used with -d or
      -x.


EXAMPLES:

\$ $my_name b121i21.fasta.screen.contigs

will read the fasta sequence file 'b121i21.fasta.screen.contigs' and
display only the full contig names, lengths, and GC percentages, and
number of bases that are not A, C, G, or T.

   b121i21.fasta.screen.Contig1	68	44.1%	0
   b121i21.fasta.screen.Contig2	3217	53.5%	3
   b121i21.fasta.screen.Contig3	12452	46.2%	0
   b121i21.fasta.screen.Contig4	29839	46.5%	0
   b121i21.fasta.screen.Contig5	65793	46.8%	0

\$ $my_name -q -s 454AllContigs.fna

will read the fasta sequence file '454AllContigs.fna' and the fasta
quality file '454AllContigs.fna.qual'.  If the file
'454AllContigs.fna.qual' does not exist, then the program will try
'454AllContigs.qual' instead.  If either of the two quality file names
can be read, then the program computes contig quality statistics.

   contig00000	1461	4.27	0	1458	5	97	35.9%	0	1439
   contig00001	148	51.98	0	144	4	97	25.0%	0	115
   contig00002	285	6.96	0	283	8	97	34.0%	0	278
   ...
   contig05261	642	76.88	0	617	4	97	58.4%	0	204
   contig05262	146	133.81	0	138	4	97	45.9%	0	118
   contig05263	123	193.92	0	114	3	97	67.5%	0	46


\$ $my_name -q -t -s mtgsp_001c04.fasta.screen.contigs

will read the fasta sequence file 'mtgsp_001c04.contigs' and the
fasta quality file 'mtgsp_001c04.contigs.qual', then compute
contig quality statistics, and summary (total) length and quality
statistics.  Displayed contig names will be shortened by removing any
prefix before the word "Contig".


   get_fasta_stats - Last Modified: $Date_Last_Modified

   Contig statistics for fasta file:'mtgsp_001c04.fasta.screen.contigs'
         and for fasta quality file:'mtgsp_001c04.fasta.screen.contigs.qual'

   Contig   	Contig	Consed	Edited	Phred20	Minimum	Maximum		NonACGT	Longest
   Name    	Length	Err/10K	Bases	Bases	Quality	Quality	GC%	bases  	Phred20
   Contig1 	71	0.16	5	66	33	78	46.5%	0	71
   Contig2 	784	3277	0	488	0	78	34.6%	1	230
   Contig3 	436	5711	0	91	0	49	35.6%	0	7
   Contig4 	1527	1196	0	1204	0	90	46.4%	0	687
   Contig5 	48867	126.15	10	48159	0	90	36.1%	0	8411
   Contig6 	63359	31.79	0	62712	0	90	34.0%	0	12483

   Number of Contigs=6, Total bp=115044, Shortest=71, Longest=63359,
   Average length=19174, Average GC%=35.1%, Non-ACGT bases=1,
   Errors/10K=130.95, Edited bp=15, Total Phred/Phrap20 bp=112720,
   Minimum Contig Quality=0, Maximum Contig Quality=90,
   Longest Run of Phred/Phrap20 bp=12483

\$ $my_name -q -t -s -p 30 mtgsp_001c04.fasta.screen.contigs

is the same as the previous example, except a Phred/Phrap30 score is
computed instead.

   get_fasta_stats - Last Modified: $Date_Last_Modified

   Contig statistics for fasta file:'mtgsp_001c04.fasta.screen.contigs'
         and for fasta quality file:'mtgsp_001c04.fasta.screen.contigs.qual'

   Contig   	Contig	Consed	Edited	Phred30	Minimum	Maximum		NonACGT	Longest
   Name    	Length	Err/10K	Bases	Bases	Quality	Quality	GC%	bases	Phred30
   Contig1 	71	0.16	5	66	33	78	46.5%	0	71
   Contig2 	784	3277	0	436	0	78	34.6%	1	218
   Contig3 	436	5711	0	41	0	49	35.6%	0	5
   Contig4 	1527	1196	0	1103	0	90	46.4%	0	686
   Contig5 	48867	126.15	10	47879	0	90	36.1%	0	5842
   Contig6 	63359	31.79	0	62143	0	90	34.0%	0	12395

   Number of Contigs=6, Total bp=115044, Shortest=71, Longest=63359,
   Average length=19174, Average GC%=35.1%, Non-ACGT bases=1,
   Errors/10K=130.95, Edited bp=15, Total Phred/Phrap30 bp=111668,
   Minimum Contig Quality=0, Maximum Contig Quality=90,
   Longest Run of Phred/Phrap30 bp=12395

\$ $my_name -4 -t 454AllContigs.fna

displays newbler assembled 454 contigs without qualities on the file
'454AllContigs.fna'.

   get_fasta_stats - Last Modified: $Date_Last_Modified

   Contig statistics for fasta file:'454AllContigs.fna'

   Contig          Contig          NonACGT Number
   Name            Length  GC%     bases   Reads
   contig00001     1791    40.4%   0       157
   contig00002     224     44.6%   0       8
   contig00003     297     36.0%   0       5
   contig00004     126     43.7%   0       3
   contig00005     848     39.5%   0       62
   contig00006     664     50.2%   0       15
   contig00007     110     40.0%   0       2
   contig00008     114     56.1%   0       2
   contig00009     126     42.9%   0       2
   contig00010     218     45.9%   0       2
   contig00011     208     42.8%   0       3
   contig00012     186     49.5%   0       3
   contig00013     440     40.9%   0       6
   contig00014     250     51.6%   0       3
   contig00015     159     42.8%   0       2
   contig00016     272     41.2%   0       7
   contig00017     268     41.4%   0       6
   contig00018     205     47.8%   0       2
   contig00019     117     44.4%   0       2
   contig00020     226     37.6%   0       36
   contig00021     3294    40.1%   0       350
   contig00022     991     37.4%   0       107
   contig00023     263     51.7%   0       5
   contig00024     221     51.1%   0       4
   contig00025     629     49.6%   0       13
   contig00026     191     50.3%   0       5
   contig00027     972     42.7%   0       40
   contig00028     1057    41.4%   0       65
   contig00029     1224    39.7%   0       210

   Number of Contigs=29, Total bp=15691, Shortest=110, Longest=3294,
   Average length=541.1, Average GC%=42.2%, Non-ACGT bases=0,
   Total Newbler Assembled 454 Reads=1127

\$ $my_name -x -t Contigs3

computes extended DNA mono-, di-, and tri-nucleotide statistics on the
sequence file 'Contigs3'.

   get_fasta_stats - Last Modified: $Date_Last_Modified

   Contig statistics for fasta file:'Contigs3'

   contig=Contig1

   Non-ACGT bases = 0 = 0.00% of total

   Seq         Count Percent     Seq         Count Percent            Sum    Sum%

   a             771  30.49%     t             653  25.82%           1424  56.31%
   c             554  21.91%     g             551  21.79%           1105  43.69%

   total1       2529

   aa            304  12.03%     tt            227   8.98%            531  21.00%
   ac            153   6.05%     gt            132   5.22%            285  11.27%
   ag            102   4.03%     ct             83   3.28%            185   7.32%
   at            211   8.35%
   ca            177   7.00%     tg            157   6.21%            334  13.21%
   cc            120   4.75%     gg            117   4.63%            237   9.38%
   cg            174   6.88%
   ga            144   5.70%     tc            123   4.87%            267  10.56%
   gc            158   6.25%
   ta            146   5.78%

   total2       2528

   aaa           107   4.23%     ttt            63   2.49%            170   6.73%
   aac            71   2.81%     gtt            46   1.82%            117   4.63%
   ...
   taa            50   1.98%     tta            71   2.81%            121   4.79%
   tca            44   1.74%     tga            53   2.10%             97   3.84%

   total3       2527


   contig=Contig2

   Non-ACGT bases = 0 = 0.00% of total

   ...

\$ $my_name -g -T spiro.contigs

get_fasta_stats - Last Modified: $Date_Last_Modified

Contig statistics for fasta file:'spiro.contigs'

Number of Contigs=192, Total bp=1693521, Shortest=108, Longest=166615,
Average length=8820.4, Average GC%=26.2%, Non-ACGT bases=8421
Median=1111, Mode=647 occurs 3 times

Histogram of sequence lengths

Maximum          (Each '*' represents 2 occurrences,
K Bases Count     '+' represents 1/2 '*')
0       67      |*********************************+
9       99      |*************************************************+
18      1       |+
27      5       |**+
36      4       |**
45      3       |*+
54      4       |**
63      4       |**
81      1       |+
99      1       |+
108     0       |
117     1       |+
126     0       |
135     0       |
144     0       |
153     0       |
162     1       |+
171     1       |+



DATE LAST MODIFIED: $Date_Last_Modified

EOF
  exit 0;
  } # end display_more_help
