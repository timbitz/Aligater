#!/usr/bin/env perl
#
##
##  Author: Tim Sterne-Weiler, timbitz (Oct 2014)
##  e-mail: tim.sterne.weiler@utoronto.ca
##

# This program accepts SAM format as stdin, and looks for chimeric reads.
# It is designed to be called by `liger` which pipes to samtools for BAM
# format if necessary, though this is not required with SAM input.

use warnings;
use strict;

use Cwd qw(abs_path);

# INITIALIZE
my $path = abs_path($0);
$0 =~ s/^.*\///;
$path =~ s/\/$0$//;

use FindBin;
use lib "$FindBin::Bin/../lib";

use SamBasics qw(:all);
use FuncBasics qw(intersect);

sub explode {
  my $str = shift;
  chomp($str);
  die "[liger detect]: ERROR $str\n";
}


# this stores the current set of alignments
my $alnHash = {}; #reference to hash

my $curRead = "";
my $curCount = 1;

while(my $l = <>) {
  next if($l =~ /^(#|\@(SQ|HD|RG|PG|CO))/); #header line
  chomp($l);
  my(@a) = split(/\t/, $l); # split sam
  
  explode "Invalid SAM format!\n" unless(defined($a[0]));
  if($a[0] ne $curRead and $a[0] ne "") {
    # before moving to a new read, process current read's alignments
    # from %alnHash

    # process algorithm here...

    # re-initialize for new read
    $curRead = $a[0];
    $curCount = 1;
    $alnHash = {}; # empty hash ref.
  }

  # process current read.
  next unless isMapped($a[1]); # ignore unmapped reads 
  explode "Invalid CIGAR format!" unless isCIGAR($a[5]); 
  my($start, $len) = alignPosInRead($a[5]);
  my $strand = getStrand($a[1]); #discard if - strand and strand specific?
  my $optHash = parseOpFields(\@a); 
  my $aScore = $optHash->{"AS"};
  print "$start\t$len\t$strand\t$aScore\n";
  #my $numMiss = $optHash->{"NM"};   
   

  #$sHash{$start} = [$curCount...]
  #$eHash{$end} = [$curCount...]
  #$alnHash{$curCount} = [score, info..]

  # increment for next read;
  $curCount++;
}

# Recursively process/collapse alignments.
sub processAlignRec {
  
}

