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
use FuncBasics qw(randomSeedRNG isInt shove openFileHandle);
use SequenceBasics qw(gcContent);

# INITIALIZE
my $path = abs_path($0);
$0 =~ s/^.*\///;
$path =~ s/\/$0$//;

# GLOBAL INIT
my $STRICT = 1;
my $RUNBLAST = 1;

my $bpMonoLimit = 7;
my $gcLimit = 0.8;

my $threads = 1;

GetOptions("o=s" => \$outputCore);

# don't die, explode!
sub explode {
  my $str = shift;
  chomp($str);
  die "[aligater detect]: (Input Line $.) ERROR $str\n";
}

sub reverb {
  my $str = shift;
  chomp($str);
  print STDERR "[$0]: ($.) $str\n";
}

randomSeedRNG(); # srand `time ^ $$ ^ unpack "%L*", `ps axww | gzip`;
my $rand = substr(md5_hex(rand), 0, 6);

if($RUNBLAST) {

  # check if blastn is installed and BLASTDB is set in ENV;
  system("bash", "-c", "which blastn > /dev/null 2> /dev/null") and
              die "[aligater filter]: Cannot find blastn which is required!\n";

  open(FORBLAST, ">$path/../tmp/tmp_$rand.fa") or die "Can't open tmp/tmp_$rand.fa for writing!\n";
}


# main loop, collect relevant entries and store into memory if --blast
while(my $l = <>) {
  chomp($l);
  my(@a) = split(/\t/, $l);

  # mononucleotide tract filter
  next if($seq =~ /[Aa]{$bpMonoLimit}|[Tt]{$bpMonoLimit}|[Cc]{$bpMonoLimit}|[Gg]{$bpMonoLimit}/);.
  next if(length($seq) < 44) { next; } # need at least 22bp on either side.
  if($STRICT and gcContent($seq) >= $gcLimit) { next; }

}
close FORBLAST;


if($RUNBLAST) { # lets run blast and remove ligations that aren't unique.

  system("blastn -query $path/../tmp/tmp_$rand.fa -task blastn -db human_genomic -word_size 20 -outfmt '6 sseqid sstart send qseqid sstrand pident length qstart qend qseq sseq evalue' -perc_identity 75 -culling_limit 1 -num_threads $threads > tmp_$rand.human_genomic.out");
  system("blastn -query $path/../tmp/tmp_$rand.fa -task blastn -db other_genomic -word_size 20 -outfmt '6 sseqid sstart send qseqid sstrand pident length qstart qend qseq sseq evalue' -perc_identity 75 -culling_limit 1 -num_threads $threads > tmp_$rand.other_genomic.out");
  system("blastn -query $path/../tmp/tmp_$rand.fa -task blastn -db nt -word_size 20 -outfmt '6 sseqid sstart send qseqid sstrand pident length qstart qend qseq sseq evalue' -perc_identity 75 -culling_limit 1 -num_threads $threads > tmp_$rand.nt.out");
  
  openBlastOutAndRemoveHits("$path/../tmp/tmp_$rand.human_genomic.out");
  openBlastOutAndRemoveHits("$path/../tmp/tmp_$rand.other_genomic.out");
  openBlastOutAndRemoveHits("$path/../tmp/tmp_$rand.nt.out");

  system("rm $path/../tmp/tmp_$rand.*");
}

#######################################################
#                                                     #
################# BEGIN SUBROUTINES ###################
#                                                     #
#######################################################

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
