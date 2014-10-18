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
our $HYBRID_PENALTY = -6;

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
  next unless isMapped($a[1]); # ignore unmapped reads 
  explode "Invalid SAM format!\n" unless(defined($a[0]) and isInt($a[1]));

  if($a[0] ne $curRead and $curRead ne "") {
    # before moving to a new read, process current read's alignments
    # from $alnHash

    # process algorithm here...
    processAlignRec($alnHash, $sHash, $eHash);

    #output best alignment here.
    my($bestK) = sort { alignCmp($alnHash->{$b}, $alnHash->{$a}) } keys %$alnHash;
    if(defined $bestK) {
      if($bestK =~ /\:/) { # chimeric read
        my $hybrid = getHybridFormat($bestK, $alnHash);
      } else {
        my $output = join("\t", @{$alnHash->{$bestK}});
        print "$output\n"; # best non-chimeric alignments
      }
    }

    # re-initialize for new read
    $curRead = $a[0];
    $curCount = 1;
    ($alnHash, $sHash, $eHash, $seHash) = ({}, {}, {}, {}); # empty hash refs..
  }

  # process current read.
  explode "Invalid CIGAR format!" unless isCIGAR($a[5]); 
  my($start, $len) = alignPosInRead($a[5]);
  my $end = $start + $len;
  my $strand = getStrand($a[1]); 
  next if($STRAND_SPECIFIC and $strand eq "-");#discard if - strand and strand specific?
  my $optHash = parseOpFields(\@a); 
  my $aScore = $optHash->{"AS"};
  #my $numMiss = $optHash->{"NM"};   

  #print "$curRead\t$start->$end\n";  

  # if same start/end exists.. take best;
  my $curAln = shove([$aScore, $start, $len], @a);

#  if(!defined($seHash->{"$start\:$end"}) or 
#     alignCmp($curAln, $alnHash->{ $seHash->{"$start\:$end"} }) > 0) {
#    $seHash->{"$start\:$end"} = $curCount;
#  } else { next; } # redundant alignment with lower rank
 
  # add start and end records.
  $sHash->{$start} = shove($sHash->{$start}, $curCount);
  $eHash->{$end} = shove($eHash->{$end}, $curCount);

  # record full alignment
  $alnHash->{$curCount} = $curAln; 

  # increment for next read;
  $curCount++;
  $curRead = $a[0];
}

# Recursively process/collapse alignments.
sub processAlignRec {
  my($alnHash, $sHash, $eHash) = @_;
  my $found = [];
  foreach my $end (sort keys %$eHash) {
    my $curList = $eHash->{$end};
    my $findList = [];
    if(defined($sHash->{$end})) {
      push(@$findList, @{$sHash->{$end}});
      push(@$found, $end);
    } 
    if(defined($sHash->{$end - 1})) {
      push(@$findList, @{$sHash->{$end - 1}});
      push(@$found, $end - 1);
    } 
    if(defined($sHash->{$end + 1})) {
      push(@$findList, @{$sHash->{$end + 1}});
      push(@$found, $end + 1);
    }
    
    # now go through and make hybrids
    foreach my $a (@$curList) {
      foreach my $b (@$findList) {
        my($score, $start, $len);
        explode "$a and $b are not defined in alnHash!" unless 
		defined($alnHash->{$a}) and defined($alnHash->{$b}); 
        $score = $alnHash->{$a}->[0] + $alnHash->{$b}->[0] + $HYBRID_PENALTY;
        $start = $alnHash->{$a}->[1];
        $len = ($alnHash->{$b}->[1] + $alnHash->{$b}->[2]) - $alnHash->{$a}->[1];
        my $readName = $alnHash->{$a}->[3];
        # set alnHash value for new hybrid
        my $geneIdHyb = "$alnHash->{$a}->[5]\:\-\:$alnHash->{$b}->[5]";
        $alnHash->{"$a\:$b"} = [$score, $start, $len, $readName, "hybrid", $geneIdHyb];
        $eHash->{$start+$len} = shove($eHash->{$start+$len}, "$a\:$b");
      }
    }

  } # end itor

  # now clean up those starts so we don't
  # pick them up on the next recurse
  foreach my $start (@$found) {
    next unless defined($sHash->{$start});
    delete $sHash->{$start};
  }
  # end condition if no more hybrids are found.
  (scalar @$found) ? processAlignRec($alnHash, $sHash, $eHash) : return;
}

sub getHybridFormat {
  my($hyb, $alnHash) = @_;
  explode "getHybridFormat ERROR!\n" unless (defined($hyb) and defined($alnHash));
  my(@a) = split(/\:/, $hyb);
  my $alpha = join("", ("A".."Z"));
  my(%used);
  my @struct = @a;
  my @geneStruc = @a;
  # convert chimeric read structures from gene symbol-> segment character
  foreach my $id (@a) {
    my $char = substr($alpha, scalar keys %used, 1);
    my(undef, undef, $geneSym, undef) = split(/\_/, $alnHash->{$id}->[5]);
    if(defined($used{$geneSym})) {
      $char = $used{$geneSym};
    }
    # substitute hyb number for char or sym;
    foreach my $s (@struct) { $s = ($s eq $id) ? $char : $id };
    foreach my $s (@geneStruc) { $s = ($s eq $id) ? $geneSym : $id };
    $used{$geneSym} = $char;
  }
  print "hybrid\t$struct\t$geneStruc\t$hyb\n";
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
  if($alignA->[5] =~ /\:\-\:/ and $alignB->[5] !~ /\:\-\:/) { return -1; }
  if($alignA->[5] !~ /\:\-\:/ and $alignB->[5] =~ /\:\-\:/) { return 1; }
  # if both are hybrids
  if($alignA->[5] =~ /\:\-\:/ and $alignB->[5] =~ /\:\-\:/) {
    my(@splA) = split(/\:\-\:/, $alignA->[5]);
    my(@splB) = split(/\:\-\:/, $alignB->[5]);
    $cmp = (scalar @splB <=> scalar @splA); # prefer less ligation sites
    return($cmp) if $cmp != 0;
    my(%uniqA, %uniqB);
    foreach my $id (@splA) { $uniqA{$id} = 1; }
    foreach my $id (@splB) { $uniqB{$id} = 1; }
    $cmp = (scalar(keys %uniqA) <=> scalar(keys %uniqB)); # prefer intra-molecular
    return $cmp;
  }
  # look at symbol and biotype
  my(undef, undef, $symA, $biotypeA) = split(/\_/, $alignA->[5]);
  my(undef, undef, $symB, $biotypeB) = split(/\_/, $alignB->[5]);
  $cmp = (isOverAbundant($symA) <=> isOverAbundant($symB)); #is overabundant?
  return($cmp) if $cmp != 0;
  $cmp = (biorank($biotypeA) <=> biorank($biotypeB)); # biotype ranking
  return($cmp) if $cmp != 0;
  $cmp = (length($symB) <=> length($symA)); # gene symbol length
  return($cmp) if $cmp != 0;
  $cmp = ($symA cmp $symB); # finally lexographical gene symbol
  return $cmp;
}

# rank biotype by..
sub biorank {
  my $biotype = shift;
  return 6 if ($biotype eq "misc-RNA");
  return 5 if ($biotype =~ /RNA/);
  return 3 if ($biotype =~ /intron/);
  return 2 if ($biotype !~ /protein-coding/);
  return 4;
}

# score over abundant symbols highest
sub isOverAbundant {
  my($name) = @_;;
  my $ret = 0;
  $ret += ($name =~ /ncrna|RNA|sno|5S|SNOR|Y-RNA|RNY|U\d|miR|piR/) ? 1 : 0;
  $ret += ($name =~ /RN7|7SL|7SK|SRP|MRP|RPPH|RNase/) ? 2 : 0;
  return $ret;
}

