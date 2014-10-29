#!/usr/bin/env perl
#
##
##  Author: Tim Sterne-Weiler, timbitz (Oct 2014)
##  e-mail: tim.sterne.weiler@utoronto.ca
##

use warnings;
use strict;

use Cwd qw(abs_path);
use Digest::MD5 qw(md5_hex md5_base64);

use FindBin;
use lib "$FindBin::Bin/../lib";

use Getopt::Long;

use SamBasics qw(:all);
use FuncBasics qw(randomSeedRNG max min);
use SequenceBasics qw(gcContent);

# INITIALIZE
my $path = abs_path($0);
$0 =~ s/^.*\///;
$path =~ s/\/$0$//;

my $tmpPath = "$path/../tmp";

# GLOBAL INIT
my $STRICT = 1;
my $RUNBLAST = 1;

my $blastDb = "human_genomic,other_genomic,nt";

my $bpMonoLimit = 7;
my $gcLimit = 0.8;

my $threads = 1;

GetOptions("o=s" => \$outputCore);

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

# make sure another instance doesn't try to write to the same files...
randomSeedRNG(); # srand `time ^ $$ ^ unpack "%L*", `ps axww | gzip`;
my $rand = substr(md5_hex(rand), 0, 6);


if($RUNBLAST) {

  # check if blastn is installed and BLASTDB is set in ENV;
  system("bash", "-c", "which blastn > /dev/null 2> /dev/null") and
              die "[aligater post]: Cannot find blastn which is required!\n";
  die "[aligater post]: BLASTDB environmental variable must be set!\n" unless defined($ENV{"BLASTDB"});
 
  open(FORBLAST, ">$tmpPath/tmp_$rand.fa") or die "Can't open tmp/tmp_$rand.fa for writing!\n";
}


# main loop, collect relevant entries and store into memory if --blast
while(my $l = <>) {
  chomp($l);
  my(@a) = split(/\t/, $l);

  # HARD FILTERS ----------------------#
  # mononucleotide tract filter
  next if($seq =~ /[Aa]{$bpMonoLimit}|[Tt]{$bpMonoLimit}|[Cc]{$bpMonoLimit}|[Gg]{$bpMonoLimit}/);
  next if length($seq) < 44; # need at least 22bp on either side.
  next if gcContent($seq) >= $gcLimit; # greater than limit of gc content
  #------------------------------------#
 
  
}
close FORBLAST;


if($RUNBLAST) { # lets run blast and remove ligations that aren't unique.
  foreach my $db (split(/\,/, $blastDb)) {
    runBlastn($db, "tmp_$rand", $threads);
    openBlastOutAndRemoveHits("$tmpPath/tmp_$rand.$db.out");
  }
  system("rm $tmpPath/tmp_$rand.*");
}

#######################################################
#                                                     #
################# BEGIN SUBROUTINES ###################
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
  my $rand = substr(md5_hex(rand), 0, 4);   
  system("echo \">seqA\n$seqA\" > $tmpPath/$rand.seqA.fa"); 
  system("echo \">seqB\n$seqB\" > $tmpPath/$rand.seqB.fa");
  $param = defined($param) ? "-P $param" : "";
  my(@res) = `ractip $tmpPath/$rand.seqA.fa $tmpPath/$rand.seqB.fa -e $param`;
  my undef = `rm $tmpPath/$rand*`;
  chomp @res;
  my($structA, $structB) = $res[qw(1 3)]; #set structures
  $res[6] =~ /JS\= ([\d\-\.]+)/;
  my $deltaG = $1; # parse energy
  my(@strA) = split(/(?<=\.)(?=[\]\[\(\)])|(?<=[\]\[\(\)])(?=\.)/, $structA);
  my(@strB) = split(/(?<=\.)(?=[\]\[\(\)])|(?<=[\]\[\(\)])(?=\.)/, $structB);
  my $maxInterLenA = maxLength(\@strA, "[\[\]]");
  my $maxInterLenB = maxLength(\@strA, "[\[\]]");
  my $maxInterLen = max($maxInterLenA, $maxInterLenB);
 
  return($deltaG, $structA, $structB, $maxInterLen);
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
  my($db, $basename, $threads) = @_;
  system("blastn -query $path/../tmp/$name.fa -task blastn -db $db -word_size 20 \
         -outfmt '6 sseqid sstart send qseqid sstrand pident length qstart qend qseq sseq evalue' \
         -perc_identity 75 -culling_limit 1 -num_threads $threads > $path/../tmp/$name.$db.out");
}

sub openBlastOutAndRemoveHits {
  my($filename, $filterHash) = @_;
  my $hndl = openFileHandle($filename);
  while(my $l = <$hndl>) {
    my(@a) = split(/\t/, $l);
    my(@b) = split(/\:/, $a[3]); #split name
    if($a[8] - $a[7] < 33) {
      if($a[7] > $b[1] - 8) { next; }
      if($a[8] < $b[2] + 8) { next; }
    }
    delete $toFilter{$b[0]};
  }
}

__END__
