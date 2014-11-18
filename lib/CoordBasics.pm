package CoordBasics; 
 
# Copyright (C) 2013 Tim Sterne-Weiler
# e-mail: tim.sterne.weiler@utoronto.ca

use 5.008008; 
use strict; 
use warnings; 

use FuncBasics qw(max min);
 
require Exporter; 
 
our @ISA = qw(Exporter); 
 
# Items to export into callers namespace by default. Note: do not export 
# names by default without a very good reason. Use EXPORT_OK instead. 
# Do not simply export all your public functions/methods/constants. 
 
# This allows declaration       use GeneRegion ':all'; 
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK 
# will save memory. 
our %EXPORT_TAGS = ( 'all' => [ qw( 
       coorOverlap
       coorWithin
       coorMatch
       coorExpand
       coorLength
       coorMidpoint
       coorToString
       coorSampleUnif
       strandMatch
       parseRegion 
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
       coorOverlap
       coorWithin
       coorMatch
       coorExpand
       coorLength
       coorMidpoint
       coorToString
       coorSampleUnif
       strandMatch
       parseRegion              
); 

our $VERSION = '0.01'; 



########################################################################################## 



########################################################################################### 
#                                                                                         # 
#      Subroutines (EXPORTED)                                                             # 
#                                                                                         # 
########################################################################################### 
 
# INPUT "chrX:start-stop:strand" 
# OR "chrX:start-stop" 
# OR [ chrX, start, stop, strand ] 
# RETURN array of elements 
sub parseRegion { 
  my $reg  = shift; 
  return unless defined $reg;
  if(ref($reg)) { return(@{$reg}); }
  my(@a) = split(/\:/, $reg); 
  my(@b) = split(/\-/, $a[1]); 
  if(scalar(@a) == 3) { 
    return($a[0], $b[0], $b[1], $a[2]); 
  } else { 
    return($a[0], $b[0], $b[1]); 
  } 
} 
 
# INPUT two genomeCoordinate 
# RETURN # of bases overlapping between regions. 
sub coorOverlap { 
  my($coorA, $coorB) = @_; 
  my($chrA, $stA, $enA, undef) = parseRegion($coorA); 
  my($chrB, $stB, $enB, undef) = parseRegion($coorB);
  unless($chrA eq $chrB) { return(0); }
  if($stA > $enB || $stB > $enA) { return(0); }
  return(min($enA, $enB) - max($stA, $stB)); 
}

# INPUT two genomeCoordinate 
# RETURN true if $coorA is within $coorB. 
sub coorWithin {
  my($coorA, $coorB) = @_;
  my($chrA, $stA, $enA, undef) = parseRegion($coorA);
  my($chrB, $stB, $enB, undef) = parseRegion($coorB);
  unless($chrA eq $chrB) { return(0); }
  if($stB < $stA && $enA < $enB) { return(0); }
  return(1);
}

# INPUT genomeCoordinate
# RETURN genomeCoordinate extended by $num bases in both directions.
sub coorExpand {
  my($coor, $num) = @_;
  my($chr, $st, $en, $strand) = parseRegion($coor);
  $st -= $num;
  $en += $num;
  return([$chr, $st, $en, $strand]);
}

# INPUT two genomeCoordinates
# RETURN # of bases in region 
sub coorMatch {
  my($coorA, $coorB) = @_;
  my($chrA, $stA, $enA, undef) = parseRegion($coorA);
  my($chrB, $stB, $enB, undef) = parseRegion($coorB);
  unless($chrA eq $chrB) { return(0); }
  if($stA != $stB || $enA != $enB) { return(0); }
  return($enA - $stA);
}

# INPUT genomeCoordinate
# RETURN length of coordinate
sub coorLength {
  my $coor = shift;
  my(undef, $st, $en, undef) = parseRegion($coor);
  return($en - $st);
}

sub coorMidpoint {
  my $coor = shift;
  my(undef, $st, $en, undef) = parseRegion($coor);
  return($st + (($en - $st) / 2));
}

sub coorSampleUnif {
  my($coor, $size) = @_;
  unless(defined($size) and $size >= 1) { $size = 1; }
  my $halfSize = int($size / 2);
  my($chr, $st, $en, $and) = parseRegion($coor);
  my $width = $en - $st;
  unless($width >= 1) { die "Can't sample region size <= 0\n"; }
  my $randMid = int(rand($width));
  $randMid += $st;
  return($chr, int($randMid - $halfSize), int($randMid + $halfSize), $and);
}

sub coorToString {
  my $coor = shift;
  my($chr, $st, $en, $strand) = parseRegion($coor);
  if(defined($strand)) {
    return("$chr\:$st\-$en\:$strand");
  } else {
    return("$chr\:$st\-$en");
  }
}

# INPUT two genomeCoordinates
# RETURN true or false if antisense overlapping region.
sub strandMatch {
  my($coorA, $coorB) = @_;
  my(undef, undef, undef, $randA) = parseRegion($coorA);
  my(undef, undef, undef, $randB) = parseRegion($coorB);  
  if(defined($randA) && defined($randB)) {
    if($randA ne $randB) { return(1); }
  } else {
    return(0);
  }
}



1; 
__END__ 
