#!/usr/bin/env perl
#
##
##  Author: Tim Sterne-Weiler, timbitz (Oct 2014)
##  e-mail: tim.sterne.weiler@utoronto.ca
##

use warnings;
use strict;

use Parallel::ForkManager;

use Cwd qw(abs_path);
use Digest::MD5 qw(md5_hex md5_base64);

use FindBin;
use lib "$FindBin::Bin/../lib";

use Getopt::Long;

# from ../lib
use SamBasics qw(:all);
use FuncBasics qw(isInt randomSeedRNG max min openFileHandle);
use SequenceBasics qw(gcContent);

# INITIALIZE
my $path = abs_path($0);
$0 =~ s/^.*\///;
$path =~ s/\/$0$//;

mkdir("/tmp/aligater") unless(-e "/tmp/aligater");
my $tmpPath = "/tmp/aligater";
my $libPath = "$path/../lib";

my %toFilter;  #main memory use of the program

# GLOBAL INIT

my $RUNBLAST = 0;
my $RUNRACTIP = 0;
my $GENOMECOORD;

my $blastDb = "human_genomic,other_genomic,nt";

# defaults are the lack there of
my $bpMonoLimit = "Inf";
my $gcLimit = 1;
my $pattFilter;
my $threads = 1;
my $mapqMaxMin = "Inf(Inf)>0";  #no filter.
my $mapqPrefix = ""; # no filter
my $unalnEdge = "Inf"; # no max unaln regions

##########################
my $seqIndex = 10;       # hardcoded... for .lig format.
my $mapqIndex = 12; #
##########################

my $interCrossLimit = 0;
my $interStemLimit = 0;

# set default as strict;
my $strictOpt = 0;
my $fullOpt = 0;
my $looseOpt = 0;

GetOptions("gc=f" => \$gcLimit, 
           "mono=i" => \$bpMonoLimit,
           "mqstd=s" => \$mapqMaxMin,
           "mqpref=s" => \$mapqPrefix,
           "p=i" => \$threads, 
           "strict" => \$strictOpt,
           "loose" => \$looseOpt,
           "full" => \$fullOpt,
           "ractip" => \$RUNRACTIP,
           "blast" => \$RUNBLAST,
);

#set hard filters
if($strictOpt) {
  $gcLimit = 0.8;
  $bpMonoLimit = 7;
  $interCrossLimit = 1;
  $interStemLimit = 5;
  $pattFilter = "Low|Simple_repeat";
  $mapqMaxMin = "10(0)>5";
  $mapqPrefSet = "[protein-coding_NA=1]";
  $unalnEdge = 12;
} elsif($looseOpt) {
  $gcLimit = 0.85;
  $bpMonoLimit = 9;
  $mapqMaxMin = "50(0)>5";
  $mapqPrefSet = "[protein-coding_NA=1]";
  $unalnEdge = 18; 
}

if($fullOpt) {
  $RUNRACTIP = 1;
  $RUNBLAST = 1;
}

my $blastThreads;
if($RUNBLAST) {
  $blastThreads = $threads;
  $threads = 1;
}

# don't die, explode!
sub explode {
  my $str = shift;
  chomp($str);
  die "[aligater post]: (Input Line $.) ERROR $str\n";
}

sub reverb {
  my $str = shift;
  chomp($str);
  print STDERR "[$0]: ($.) $str\n";
}

sub checkSoft {
  my $prog = shift;
  system("bash", "-c", "which $prog > /dev/null 2> /dev/null") and
              explode "Cannot find $prog which is required!";
}

# make sure another instance doesn't try to write to the same files...
randomSeedRNG(); # srand `time ^ $$ ^ unpack "%L*", `ps axww | gzip`;
my $rand = substr(md5_hex(rand), 0, 6);

## CHECK IF SOFTWARE IS INSTALLED ------------------------#
if($RUNBLAST) {
  # check if blastn is installed and BLASTDB is set in ENV;
  checkSoft("blastn");
  explode "BLASTDB environmental variable must be set!" unless defined($ENV{"BLASTDB"});
 
  open(FORBLAST, ">$tmpPath/tmp_$rand.fa") or die "Can't open tmp/tmp_$rand.fa for writing!\n";
}
if($RUNRACTIP) {
  checkSoft("ractip");
  my(undef, $racVer) = split(/\s|\./, `ractip -V`);
  explode "ractip version must be > 1.0.0 !" unless $racVer >= 1;
}
#---------------------------------------------------------#

# Parse LIGQ Options -------------------------------------#
my($mapqIdent, $mapqDiff) = split(/\>/, $mapqMaxMin);
my $mapqSing = ($mapqIdent =~ /\((\d+)\)/);
$mapqIdent =~ s/\(\d+\)//;
$mapqSing = "Inf" if($mapqSing eq "");
unless(defined($mapqIdent) and defined($mapqDiff) and isInt($mapqIdent) and isInt($mapqDiff)) {
  explode "Improper `--mqstd` format! INT>INT or INT(INT)>INT !\n";
}
my $mapqPrefHash = {};
my(@keySets) = split(/\,/, $mapqPrefix);
foreach my $set (@keySets) {
  my($k,$v) = ($set =~ /\[(\S+)\=(\d+)\]/);
  unless(defined($k) and defined($v) and isInt($v)) {
    explode "Improper `--mqpref` format! \'[biotype-here=int],[optional-bio=int]'\n";
  }
  $mapqPrefHash->{$k} = $v;
}
# ------------------------------------------------------- #

## Set up fork manager;
my $pm = Parallel::ForkManager->new($threads, $tmpPath);

# data structure retrieval and handling
$pm -> run_on_finish ( # called BEFORE the first call to start()
  sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
    # retrieve data structure from child
    if (defined($data_ref)) {  # children are not forced to send anything
      my($key, $mid, $val) = @{$data_ref};  # child passed a string reference
      if($RUNBLAST) {
        $toFilter{$key} = $val;
        print FORBLAST ">$key\:$mid";
      } else {
        print $key;
      }
    } else {  # problems occuring will throw a warning
      reverb "No message received from child process $pid!\n";
    }
  }
);

# main loop, collect relevant entries and store into memory if --blast
while(my $l = <>) {
   
  if($threads > 1) {
    $pm->start() and next;  # fork to child
  }

  chomp($l);
  my(@a) = split(/\t/, $l);

  my $seq = $a[$seqIndex];
  my $gcSeq  = $seq;
  $gcSeq =~ s/\_//g;
  my $gcContent = gcContent($gcSeq);

  # HARD FILTERS ----------------------------------------------------#
  # mononucleotide tract filter
  my $hardFilt = 0;
  $hardFilt++ if($seq =~ /[Aa]{$bpMonoLimit}|[Tt]{$bpMonoLimit}|[Cc]{$bpMonoLimit}|[Gg]{$bpMonoLimit}/);
  $hardFilt++ if length($seq) < 45; # need at least 22bp on either side.
  $hardFilt++ if $gcContent >= $gcLimit; # greater than limit of gc content
  $hardFilt++ if defined($pattFilter) and $l =~ /$pattFilter/;  #option specific filter.

  # filter based on ligq
  my($mapqFore, $mapqLeft, $mapqRight) = split(/\>/, $a[$mapqIndex]);
  my $mapqParen = ($mapqLeft =~ /\((\d+)\)/);
  $mapqLeft =~ s/\(\d+\)//;
  $mapqParen = 0 if($mapqParen eq "");
  $hardFilt++ if $mapqLeft > $mapqIdent;
  $hardFilt++ if $mapqRight < $mapqDiff;
  $hardFilt++ if $mapqParen > $mapqSing;

  $pm->finish if($hardFilt and $threads > 1); # if threads > 1
  next if $hardFilt;  # if threads == 1
  #------------------------------------------------------------------#

  my($seqA, $seqB) = split(/\_/, $seq); 

  my($dG, $strA, $strB, $len, $amt) = ("","","","","");

  if($RUNRACTIP) {
    ($dG, $strA, $strB, $len, $amt) = runRactIP($seqA, $seqB, undef);

    # DEPRECATED: for debugging:
    # my($dG_ua, $strA_ua, $strB_ua, $len_ua) = runRactIP($seqA, $seqB, "$libPath/rna_andronescu2007_ua.par");

    # HARD FILTERS POST RACTIP-----------------------------------------#
    if($a[0] eq "I") {  #perhaps just reset the code to R instead of these hard filters TODO!!
      $hardFilt++ unless($len >= $interStemLimit);
      $hardFilt++ unless($amt >= $interCrossLimit);

      $pm->finish if($hardFilt and $threads > 1); # if threads > 1
      next if $hardFilt;  # if threads == 1
    }
    #------------------------------------------------------------------#
  }
  my($key, $mid, $val);
  if($RUNBLAST) { # if we are using blast to filter we need to store ligs in memory
    # loop through each ligation site
    while($seq =~ /\_/g) {
      my $pos = $-[0];
      my($leftCoor, $rightCoor) = ( max(0, $pos - 20), min($pos + 20, length($seq)) );
      my $ligString = substr( $seq, $leftCoor, ($rightCoor + 1) - $leftCoor );
      $ligString =~ s/\_//g;
      my($leftLig, $rightLig) = ( $pos - $leftCoor, $rightCoor - $pos );
      #save read and print for blast
#      $toFilter{"LIG_$."} = "$l\t$gcContent\t$strA\t$strB\t$dG\t$len\t$amt\n";
      $key = "LIG_$.";
      $mid = "$leftLig\:$rightLig\n$ligString\n";
      $val = "$l\t$gcContent\t$strA\t$strB\t$dG\t$len\t$amt\n";
      #print FORBLAST ">LIG_$.\:$leftLig\:$rightLig\n$ligString\n";
    }
  } else { # no need to waste memory, lets just print as we go. 
    #print "$l\t$gcContent\t$strA\t$strB\t$dG\t$len\t$amt\n"; 
    $key = ($RUNRACTIP) ? "$l\t$gcContent\t$strA\t$strB\t$dG\t$len\t$amt\n" : "$l\n";
  }
  if($threads > 1) {
    $pm->finish(0, [$key, $mid, $val]);
  } else {
    if($RUNBLAST) {
      $toFilter{$key} = $val;
      print FORBLAST ">$key\:$mid";
    } else {
      print $key;
    }
  }
} # end main loop
$pm->wait_all_children if $threads > 1;
close FORBLAST;


if($RUNBLAST) { # lets run blast and remove ligations that aren't unique.
  foreach my $db (split(/\,/, $blastDb)) {
    runBlastn($db, "tmp_$rand", $blastThreads, $tmpPath);
    openBlastOutAndRemoveHits("$tmpPath/tmp_$rand.$db.out", \%toFilter);
  }
  system("rm $tmpPath/tmp_$rand.*");
  
  # now print the remaining results.
  foreach my $key (keys %toFilter) {
    print "$toFilter{$key}";
  }
}
## END MAIN ##

#######################################################
#                                                     #
######            BEGIN SUBROUTINES            ########
#                                                     #
#######################################################

# this function runs the RactIP program for RNA-RNA interaction prediction
# using dynamic programming using the -e parameter and an optional -P param file
# returned are: deltaG, the first structure (bracket notation), second structure,
# followed by the maximum intermolecular interaction stem length ( [[[[[ or ]]]]] )
# this doesn't yet make use of the z-score function or check that the program is
# properly installed or of the correct version
sub runRactIP {
  my($seqA, $seqB, $param) = @_;
  $seqA =~ s/T/U/g if($seqA =~ /T/);
  $seqB =~ s/T/U/g if($seqB =~ /T/);
  my $rand = substr(md5_hex($.), 0, 12);   
  system("echo \">seqA\n$seqA\n>seqB\n$seqB\" > $tmpPath/$rand\_seq.fa"); 
  $param = defined($param) ? "-P $param" : "";
  my(@res) = `ractip $tmpPath/$rand\_seq.fa -e $param`;
  system("rm $tmpPath/$rand\_seq.fa");
  chomp @res;
  my($structA, $structB) = ($res[2], $res[5]); #set structures
  $res[6] =~ /JS\= ([\d\-\.]+)/;
  my $deltaG = $1; # parse energy
  my(@strA) = split(/(?<=\.)(?=[\]\[\(\)])|(?<=[\]\[])(?=[\.\(\)])|(?<=[\(\)])(?=[\.\[\]])/, $structA);
  my(@strB) = split(/(?<=\.)(?=[\]\[\(\)])|(?<=[\]\[])(?=[\.\(\)])|(?<=[\(\)])(?=[\.\[\]])/, $structB);
  my $maxInterLenA = maxLength(\@strA, "[\\[\\]]");
  my $maxInterLenB = maxLength(\@strB, "[\\[\\]]");
  my $maxInterLen = max($maxInterLenA, $maxInterLenB);
  my $amtNum = findUAinStem(\@strA, \@strB, $seqA, $seqB, "[\\[\\]]");
  
  return($deltaG, $structA, $structB, $maxInterLen, $amtNum);
}

# look for diagonal U's in stem given ractip
sub findUAinStem {
  my($structAref, $structBref, $seqA, $seqB, $char) = @_;
  my $intSeqA = "";
  my $revSeqB = "";
  my @stemA;
  my @stemB;
  my $cnt = 0;
  foreach my $seg (@$structAref) {
    ($cnt += length($seg) and next) unless $seg =~ /^$char+$/;
    $intSeqA .= substr($seqA, $cnt, length($seg));
    push(@stemA, "." x length($seg));
    $cnt += length($seg);
  }
  $cnt = 0;
  foreach my $seg (@$structBref) {
    ($cnt += length($seg) and next) unless $seg =~ /^$char+$/;
    $revSeqB .= substr($seqB, $cnt, length($seg));
    push(@stemB, "." x length($seg));
    $cnt += length($seg);
  }
  @stemB = reverse @stemB;
  $revSeqB = scalar reverse $revSeqB;
#  print "$intSeqA\n$revSeqB\n"; #debugging
  my %xlinkPos;
  # now check for diagonal Us,  AB [[ to ]] DC
  $cnt = 0;
  foreach my $seg (@stemA) {
    my $segA = substr($intSeqA, $cnt, length($seg));
    my $segB = substr($revSeqB, $cnt, length($seg));
    countDiagU($segA, $segB, \%xlinkPos, $cnt);
    $cnt += length($seg);
  }
  $cnt = 0;
  foreach my $seg (@stemB) {
    my $segA = substr($intSeqA, $cnt, length($seg));
    my $segB = substr($revSeqB, $cnt, length($seg));
    countDiagU($segA, $segB, \%xlinkPos, $cnt);
    $cnt += length($seg);
  }
  return(scalar(keys %xlinkPos));
}

# check for diagonal Us,  AB [[ to ]] DC
#                         .U    to    .U or U. to U.
sub countDiagU {
  my($stemA, $stemB, $xlinkHash, $offset) = @_;
  for(my $i=0; $i < length($stemA) - 1; $i++) {
    my $a = substr($stemA, $i, 1);
    my $b = substr($stemA, $i+1, 1);
    my $c = substr($stemB, $i, 1);
    my $d = substr($stemB, $i+1, 1);
    $xlinkHash->{$offset + $i} = 1 if($a =~ /u|U/ and $d =~ /u|U/);
    $xlinkHash->{$offset + $i+1} = 1 if($b =~ /u|U/ and $c =~ /u|U/);
  }
  #return void. 
}

# used by the runRactIP program.
sub maxLength {
  my($aRef, $char) = @_;
  my $maxLen = 0;
  foreach my $elem (@$aRef) {
    my $l = length($elem);
    next unless $elem =~ /^$char+$/;
    $maxLen = ($l > $maxLen) ? $l : $maxLen;
  }
  return($maxLen);
}

sub runBlastn {
  my($db, $basename, $threads, $tmpPath) = @_;
  system("blastn -query $tmpPath/$basename.fa -task blastn -db $db -word_size 24 -outfmt '6 sseqid sstart send qseqid sstrand pident length qstart qend qseq sseq evalue' -perc_identity 75 -culling_limit 1 -num_alignments 5 -num_threads $threads > $tmpPath/$basename.$db.out");
}

sub openBlastOutAndRemoveHits {
  my($filename, $filterHash) = @_;
  my $hndl = openFileHandle($filename);
  while(my $l = <$hndl>) {
    my(@a) = split(/\t/, $l);
    my(@b) = split(/\:/, $a[3]); #split name
    if($a[8] - $a[7] < 26) {
      if($a[7] > $b[1] - 6) { next; }
      if($a[8] < $b[2] + 6) { next; }
    }
    delete $filterHash->{$b[0]};
  }
}

__END__
