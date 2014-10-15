#!/usr/bin/env perl
#
##
##  Author: Tim Sterne-Weiler, timbitz (Oct 2014)
##  e-mail: tim.sterne.weiler@utoronto.ca
##

# This program accepts SAM format as stdin, and looks for chimeric reads.
# It is designed to be called by `aligater` which pipes to samtools for BAM
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
use FuncBasics qw(isInt shove);

our $STRAND_SPECIFIC = 1; # transcriptome mapping.

sub explode {
  my $str = shift;
  chomp($str);
  die "[aligater detect]: (Input Line $.) ERROR $str\n";
}

# set iteraters
# this stores the current set of alignments
my $alnHash = {}; #reference to hash
my $sHash = {}; #alignment starts (in read)
my $eHash = {}; #alignment ends (in read)
my $seHash = {}; #alignment "start:end"

my $curRead = "";
my $curCount = 1;

# iterate through sam alignments
while(my $l = <>) {
  next if($l =~ /^(#|\@(SQ|HD|RG|PG|CO))/); #header line
  chomp($l);
  my(@a) = split(/\t/, $l); # split sam
  
  explode "Invalid SAM format!\n" unless(defined($a[0]) and isInt($a[1]));
  if($a[0] ne $curRead and $a[0] ne "") {
    # before moving to a new read, process current read's alignments
    # from $alnHash

    # process algorithm here...
    processAlignRec($alnHash, $sHash, $eHash);

    #output best alignment here.

    # re-initialize for new read
    $curRead = $a[0];
    $curCount = 1;
    ($alnHash, $sHash, $eHash, $seHash) = ({}, {}, {}, {}); # empty hash refs..
  }

  # process current read.
  next unless isMapped($a[1]); # ignore unmapped reads 
  explode "Invalid CIGAR format!" unless isCIGAR($a[5]); 
  my($start, $len) = alignPosInRead($a[5]);
  my $end = $start + $len;
  my $strand = getStrand($a[1]); 
  next if($STRAND_SPECIFIC and $strand eq "-");#discard if - strand and strand specific?
  my $optHash = parseOpFields(\@a); 
  my $aScore = $optHash->{"AS"};
  print "$start\t$len\t$strand\t$aScore\n";
  #my $numMiss = $optHash->{"NM"};   
   
  # if same start/end exists.. take best;
  my $curAln = shove([$aScore, $start, $len], @a);
  if(alignCmp($curAln, $alnHash{ $seHash{"$start\:$end"} }) > 0) {
    $seHash{"$start\:$end"} = "$curCount";
  } else next;  # redundant alignment with lower rank

  # add start and end records.
  $sHash->{$start} = shove($sHash->{$start}, $curCount);
  $eHash->{$end} = shove($eHash->{$end}, $curCount);

  # record full alignment
  $alnHash->{"$curCount"} = $curAln; 

  # increment for next read;
  $curCount++;
}

# Recursively process/collapse alignments.
sub processAlignRec {
  my($alnHash, $sHash, $eHash) = @_;
  foreach my $end (sort keys %$eHash) {
    if(defined($sHash{$end}) { # penalty -5
      
    } elsif(defined($sHash{$end - 1}) or # pen -6
            defined($sHash{$end + 1})) {
      
    } 
    
  }
}

# compare two alignments first by alignment score
# then by symbol and biotypes.
sub alignCmp {
  my($alignA, $alignB) = @_;

  # explode if neither is defined.
  explode "Invalid object $alignA or $alignB!" unless defined($alignA or $alignB);
  return -1 unless defined($alignA); # a is defined but b is not
  return 1 unless defined($alignB); # b is defined by a is not

  my $cmp = ($alignA->[0] <=> $alignB->[0]); #alignment score
  return($cmp) if $cmp != 0;
  $cmp = ($alignB->[2] <=> $alignA->[2]); #length of alignment
  return($cmp) if $cmp != 0;
  # bias against hybrids when vs non-hybrids
  if($alignA->[5] eq "Hybrid" and $alignB->[5] ne "Hybrid") { return -1; }
  if($alignA->[5] ne "Hybrid" and $alignB->[5] eq "Hybrid") { return 1; }
  if($alignA->[5] eq "Hybrid" and $alignB->[5] eq "Hybrid") { return 0; }
  # look at symbol and biotype
  my(undef, undef, $symA, $biotypeA) = split(/\_/, $alignA->[5]);
  my(undef, undef, $symB, $biotypeB) = split(/\_/, $alignB->[5]);
  $cmp = (isOverAbundant($symA) <=> isOverAbundant($symB)); #is overabundant?
  return($cmp) if $cmp != 0;
  $cmp = (biorank($biotypeA) <=> biorank($biotypeB)); # biotype ranking
  return($cmp) if $cmp != 0;
  $cmp = (length($symB) <=> length($symA)); # gene symbol length
  return $cmp;
}

# rank biotype by..
sub biorank {
  my $biotype = shift;
  return 6 if ($biotype eq "misc-RNA");
  return 4 if ($biotype =~ /RNA/);
  return 3 if ($biotype =~ /intron/);
  return 2 if ($biotype !~ /protein-coding/);
  return 1;
}

# score over abundant symbols highest
sub isOverAbundant {
  my($name) = @_;;
  my $ret = 0;
  $ret += ($name =~ /ncrna|RNA|sno|5S|SNOR|Y-RNA|RNY|U\d|miR|piR/) ? 1 : 0;
  $ret += ($name =~ /RN7|7SL|7SK|SRP|MRP|RPPH|RNase/) ? 2 : 0;
  return $ret;
}

