package GeneRegion;

# Copyright (C) 2012 Tim Sterne-Weiler
# e-mail: tim.sterne.weiler@utoronto.ca

use 5.008008;
use strict;
use warnings;

#use IPC::Open3::Utils qw(:all); 
use DBI;

use hgSQLBasics qw(:all);
use FuncBasics qw(:all);
use CoordBasics qw(:all);

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use GeneRegion ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	getGenomicRegion
	getSmallTrackDist
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	getGenomicRegion
	getSmallTrackDist
);

our $VERSION = '0.01';



#our $REGION = undef;  # SET THIS MANUALLY BEFORE CALLING EXPORTED FUNCTIONS

our $INTRONRANGE = 100000;
our $EXONRANGE = 75;


###########################################################################################
#                                                                                         #
#      Subroutines (EXPORTED)                                                             #
#                                                                                         #
###########################################################################################

our %RUNQ;

##
## Get the distance to small annotation track of some sort.  for example microRNA sites.
##
sub getSmallTrackDist {
  my($coord, $table) = @_;
  unless(defined($table)) { return undef; }
  unless(defined($RUNQ{$table})) { $RUNQ{$table} = prepareSmallCoordSQLquery(undef, $table, undef, undef, "ORDER BY kg.chromStart"); }
  # return (name, chrom, strand, start, end) in resultsRef
  my($resultsRef) = executeExactCoordSQLquery($RUNQ{$table}, $coord);
  unless(scalar(@$resultsRef) > 0) { return undef; }
  my $point = coorMidpoint($coord);
  my $minDist = "Inf";
  my $ind = "Inf";
  for(my $i=0; $i < scalar(@$resultsRef); $i++) {  #iterate through results in sorted order ASCENDING
    my($curSt, $curEn) = ($resultsRef->[$i]->[3], $resultsRef->[$i]->[4]);
    my $curMid = $curSt + (($curEn - $curSt) / 2);
    my $dist = $point - $curMid;
    if(abs($dist) > abs($minDist)) { 
      last;
    } else {
      $minDist = $dist;
      $ind = $i;
    }
  }
  if($minDist eq "Inf" or $ind eq "Inf") { return undef; }
  if($resultsRef->[$ind]->[2] eq "-") { $minDist *= -1; }
  return($minDist, $resultsRef->[$ind]);
}

##
##  Function to determine intragenic anchor from gene annotation tracks (genePred format?)
##
sub getGenomicRegion {
  my($coord, $tableA, $tableB) = @_;
  unless(defined($tableA)) { return undef; }
  unless(defined($RUNQ{$tableA})) { $RUNQ{$tableA} = prepareCoordSQLquery(undef, $tableA); }
  unless(defined($RUNQ{$tableB})) { $RUNQ{$tableB} = prepareCoordSQLquery(undef, $tableB); }
  # this is probably a more conservative first query.. like CCDS annotations
  my($resultsRefA) = executeCoordSQLquery($RUNQ{$tableA}, $coord, undef);
  my($resultsRefB);
  # nothing more to do here if...
  if(scalar(@$resultsRefA) == 0 && !defined($tableB)) { return undef; }
  # parse through isoforms and choose best region/exon.
  my($cdsOrUtrA, $locZoneA, $distA, $bestA, $nameA) = overlapWithExons($coord, $resultsRefA);
  my($cdsOrUtrB, $locZoneB, $distB, $bestB, $nameB);
  my($startDistA, $stopDistA);
  my($startLeftA, $startRightA, $stopLeftA, $stopRightA);
  if($locZoneA =~ /EXON/ && defined($nameA)) {          # hotwire to resultsRefB as knownGene or someUTR containing results
    $resultsRefB = executeCoordSQLquery($RUNQ{$tableB}, $coord, undef);
    ($startDistA, $startLeftA, $startRightA) = distanceToStartStopCodon($coord, $resultsRefB, $nameA, "START");
    ($stopDistA, $stopLeftA, $stopRightA)    = distanceToStartStopCodon($coord, $resultsRefB, $nameA, "STOP");
  }
  # if resultsRefA <= 0 or there was no best exon...
  if(!defined($cdsOrUtrA) or $locZoneA eq "DEEP_INTER") {
    $resultsRefB = executeCoordSQLquery($RUNQ{$tableB}, $coord, undef);
    ($cdsOrUtrB, $locZoneB, $distB, $bestB, $nameB) = overlapWithExons($coord, $resultsRefB);
    my($startDistB, $stopDistB);
    my($startLeftB, $startRightB, $stopLeftB, $stopRightB);
    if(!defined($startDistA) && $locZoneB =~ /EXON/) {
      ($startDistB, $startLeftB, $startRightB) = distanceToStartStopCodon($coord, $resultsRefB, $nameB, "START");
      ($stopDistB, $stopLeftB, $stopRightB)    = distanceToStartStopCodon($coord, $resultsRefB, $nameB, "STOP");
    }
    if(defined($distB) && !defined($distA)) {
      return ($cdsOrUtrB, $locZoneB, $distB, $bestB, $nameB, $startDistB, $stopDistB,
	      			  $startLeftB, $startRightB, $stopLeftB, $stopRightB);
    } elsif(!defined($distB)) {
      return ($cdsOrUtrA, $locZoneA, $distA, $bestA, $nameA, $startDistA, $stopDistA,
				  $startLeftA, $startRightA, $stopLeftA, $stopRightA);
    } elsif($distA < $distB) {
      return ($cdsOrUtrA, $locZoneA, $distA, $bestA, $nameA, $startDistA, $stopDistA,
	      			  $startLeftA, $startRightA, $stopLeftA, $stopRightA);      
    } else {
      return ($cdsOrUtrB, $locZoneB, $distB, $bestB, $nameB, $startDistB, $stopDistB,
	      			  $startLeftB, $startRightB, $stopLeftB, $stopRightB);
    }
  }
  return ($cdsOrUtrA, $locZoneA, $distA, $bestA, $nameA, $startDistA, $stopDistA,
			      $startLeftA, $startRightA, $stopLeftA, $stopRightA);  
}


###########################################################################################
#                                                                                         #
#      Subroutines (BACKEND)                                                              #
#                                                                                         #
###########################################################################################

sub distanceToStartStopCodon {
  my($coord, $resultRef, $name, $codon, $recurse) = @_;
  my($chrom, $start, $stop, undef) = parseRegion($coord);
  my($midPoint) = int(($stop - $start) / 2) + $start;
  my($minDist);
  my($txLeft, $txRight); #controls for plotting expected based on uniform distribution
#  my($corrTxLength); # txLength corresponding to transcript giving minDist.
  foreach my $rowRef (@$resultRef) {
#    unless($rowRef->[0] eq $name) { next; }
    my(@exStarts) = split(/\,/, $rowRef->[5]);
    my(@exEnds) = split(/\,/, $rowRef->[6]);
    my($endCond, $initItor, $cdsMark, $inc);
    # set loop values for forward or reverse iteration...`
    if(($codon eq "START" && $rowRef->[2] eq "+") ||
       ($codon eq "STOP" && $rowRef->[2] eq "-")) {
      ($endCond, $initItor, $cdsMark, $inc) = (scalar(@exStarts), 0, $rowRef->[3], 1);
    } else {
      ($endCond, $initItor, $cdsMark, $inc) = (scalar(@exStarts), 0, $rowRef->[4], 1);
#      ($endCond, $initItor, $cdsMark, $inc) = (scalar(@exStarts)+1, -1, $rowRef->[4], -1);
    }
    my($subtract) = 0;
    my($subFLAG) = 0;  # once we are in range to compute dist
    for(my $i=$initItor; abs($i) < $endCond; $i += $inc) {
      if(coorOverlap([ $chrom, $midPoint, $midPoint+1 ], 
		     [ $chrom, $exStarts[$i], $exEnds[$i] ])) {
        $subFLAG++;
      }
      if(coorOverlap([ $chrom, $cdsMark-1, $cdsMark+1 ],
  	             [ $chrom, $exStarts[$i], $exEnds[$i] ])) {
        $subFLAG++;
      } 
      if($subFLAG == 2) {
         # COMPUTE DISTANCE - $subtract
         my($dist) = abs($midPoint - $cdsMark) - $subtract;
         if($midPoint < $cdsMark) { $dist *= -1; }
         if($rowRef->[2] eq "-") { $dist *= -1; }
         # short circuit logic here..
         if(!defined($minDist) || abs($minDist) > abs($dist)) { 
           $minDist = $dist;
           if(!defined($recurse)) {
       #       print STDERR "recursing, $recurse, $minDist, $rowRef->[7], $rowRef->[8]\n";
	     if($rowRef->[2] eq "+") {
               ($txLeft, undef, undef) = distanceToStartStopCodon([$chrom, $rowRef->[7], $rowRef->[7], undef], [$rowRef], $name, $codon, 1);
               ($txRight, undef, undef) = distanceToStartStopCodon([$chrom, $rowRef->[8]-1, $rowRef->[8]-1, undef], [$rowRef], $name, $codon, 1);
	     } else {
               ($txRight, undef, undef) = distanceToStartStopCodon([$chrom, $rowRef->[7], $rowRef->[7], undef], [$rowRef], $name, $codon, 1);
               ($txLeft, undef, undef) = distanceToStartStopCodon([$chrom, $rowRef->[8]-1, $rowRef->[8]-1, undef], [$rowRef], $name, $codon, 1);
             }
#	     if($rowRef->[2] eq "-") {
#               ($txLeft, $txRight) = ($txLeft*-1, $txRight*-1);
#             }
           } #else { print STDERR "$recurse\n"; }
         }
	 last; # go to next isoform.
      } elsif($subFLAG == 1) {
        # subtract next intron unless there isn't one..
        unless(abs($i)+1 == $endCond) {
          $subtract += $exStarts[$i+(1+$initItor)] - $exEnds[$i+$initItor];
        } #else { print STDERR "starts[$i+(1+$initItor)] - ends[$i+$initItor], $endCond\n"; }
      }
    } #end loop of current isoform
  }
#  print STDERR "returning $minDist, $txLeft, $txRight,  recurse=$recurse\n";
  return($minDist, $txLeft, $txRight);
}

# INPUT: coord and reference to result array
# OUTPUT: ("CDS"/"UTR, "DEEP/3'SS/5'SS/MID_INT/CDS", distanceToSS, bestExon, isoformName)
sub overlapWithExons {
  my($coord, $resultRef) = @_;
  my($chrom, $start, $stop, undef) = parseRegion($coord);
  my($midPoint) = int(($stop - $start) / 2) + $start;
  my(%bestExon);
  my(%geneName);
  my(%exonInd);
  my(%regOver);
  my($cdsOrUtr) = "INTER";
  my($locZone) = "DEEP_INTER";
  if(scalar(@$resultRef) > 0) { $cdsOrUtr = "UTR"; } # default..
  foreach my $rowRef (@$resultRef) {
    my(@exStarts) = split(/\,/, $rowRef->[5]);
    my(@exEnds) = split(/\,/, $rowRef->[6]);
    # if stop is >= cdsStart && start <= cdsEnd;
    if($stop >= $rowRef->[3] && $start <= $rowRef->[4]) { $cdsOrUtr = "CDS"; }
    for(my $i=0; $i < scalar(@exStarts); $i++) {
      my($exMidpoint) = int(($exEnds[$i] - $exStarts[$i]) / 2) + $exStarts[$i];
      # Set regions for current exon.
      my($leftIntCoor) = [ $chrom, $exStarts[$i] - $INTRONRANGE, $exStarts[$i] ];
      my($leftCdsCoor) = [ $chrom, $exStarts[$i], min($exStarts[$i] + $EXONRANGE, $exMidpoint) ];
      my($midCdsCoor) = [ $chrom, min($exStarts[$i] + $EXONRANGE, $exMidpoint), 
			  	  max($exEnds[$i] - $EXONRANGE, $exMidpoint) ];
      my($rightCdsCoor) = [ $chrom, max($exEnds[$i] - $EXONRANGE, $exMidpoint), $exEnds[$i] ];
      my($rightIntCoor) = [ $chrom, $exEnds[$i], $exEnds[$i] + $INTRONRANGE ];
      # Compute overlap with target coord..
      my($leftIntOver) = coorOverlap($coord, $leftIntCoor);
      my($leftCdsOver) = coorOverlap($coord, $leftCdsCoor);
      my($midCdsOver) = coorOverlap($coord, $midCdsCoor);
      my($rightCdsOver) = coorOverlap($coord, $rightCdsCoor);
      my($rightIntOver) = coorOverlap($coord, $rightIntCoor);
      if($i == 0) { $leftIntOver = 0; }
      if($i == $#exStarts) { $rightIntOver = 0; }
      # this is a speed checkpoint. short circuit logic without comparisons..
      unless($leftIntOver || $leftCdsOver || $midCdsOver || $rightCdsOver || $rightIntOver) { next; }
      # get index and max value..
      my($ind, $val) = maxOfArray( $leftIntOver, $leftCdsOver, $midCdsOver, $rightCdsOver, $rightIntOver );
      # print "$ind\t$val\n";
      $bestExon{"$chrom\:$exStarts[$i]\-$exEnds[$i]\:$rowRef->[2]"} += calcDistance($ind, $exStarts[$i], $exEnds[$i], $midPoint);
      $regOver{"$chrom\:$exStarts[$i]\-$exEnds[$i]\:$rowRef->[2]"} = [ $ind, $val ];  # set to best overlapping region   
      $geneName{"$chrom\:$exStarts[$i]\-$exEnds[$i]\:$rowRef->[2]"} = $rowRef->[0];
      $exonInd{"$chrom\:$exStarts[$i]\-$exEnds[$i]\:$rowRef->[2]"} = ($i+1)."\/".($#exStarts+1);
    }
  }
  if(scalar(keys %bestExon) == 0) { return($cdsOrUtr, $locZone, undef, undef, undef); }
  my($best) = (sort { abs($bestExon{$a}) <=> abs($bestExon{$b}) } keys %bestExon)[0];
  my(undef, $leftSS, $rightSS, $geneStrand) = parseRegion($best);
#  my($distance) = calcDistance($regOver{$best}->[0], $leftSS, $rightSS, $midPoint);
  $locZone = exonRegToString($regOver{$best}->[0], $geneStrand, split(/\//, $exonInd{$best}));
  if($cdsOrUtr eq "CDS" and $locZone =~ /TSS/) { $locZone =~ s/TSS/STARTCODON/; }
  if($cdsOrUtr eq "CDS" and $locZone =~ /PAS/) { $locZone =~ s/PAS/STOPCODON/; }
  return($cdsOrUtr, $locZone, $bestExon{$best}, $best, $exonInd{$best}, $geneName{$best});  ## HACK added exonInd{$best}.. geneName ret undef 
}

# INPUT: region in index form from overlapExons and strand from gene..
# OUTPUT: text version of region in human readable format.

# KINDOF DIRTY... SORRY...
sub exonRegToString {
  my($reg, $strand, $targExInd, $totalExInd) = @_;
  my($class);
  if($reg == 4 or $reg == 0) { $class = "INTRON"; }
  if($reg == 1 or $reg == 2 or $reg == 3) { $class = "EXON"; }
  if($totalExInd == 1) { $class = "EXON"; }
  if(($reg == 1 && $strand eq "+" && $targExInd == 1) or 
     ($reg == 3 && $strand eq "-" && $targExInd == $totalExInd)) { return "TSS_$class"; }
  if(($reg == 1 && $strand eq "-" && $targExInd == 1) or
     ($reg == 3 && $strand eq "+" && $targExInd == $totalExInd)) { return "PAS_$class"; }
  if(($reg == 0 && $strand eq "+") or ($reg == 4 && $strand eq "-")) { return "3'SS_$class"; }
  if(($reg == 1 && $strand eq "+") or ($reg == 3 && $strand eq "-")) { return "3'SS_$class"; }
  if($reg == 2) { return "MID_$class"; }
  if(($reg == 3 && $strand eq "+") or ($reg == 1 && $strand eq "-")) { return "5'SS_$class"; }
  if(($reg == 4 && $strand eq "+") or ($reg == 0 && $strand eq "-")) { return "5'SS_$class"; }
  else { die "ERROR: reg # is > 4!!!\n"; } # this will never happen.
}

# INPUT: region in index form, leftSS of exon, rightSS of exon, coord midpoint. (strand unspecific)
# OUTPUT: distance.
sub calcDistance {
  my($reg, $leftSS, $rightSS, $mid) = @_;
  if($reg == 0) { return($leftSS - $mid); }
  if($reg == 1) { return($mid - $leftSS); }
  if($reg == 2) { return(min(abs($mid - $leftSS), abs($rightSS - $mid))); }
  if($reg == 3) { return($rightSS - $mid); }
  if($reg == 4) { return($mid - $rightSS); }
  else { die "ERROR: reg # is > 4!!!\n"; } # this will never happen.
}


##########################################################################################

1;
__END__

=head1 NAME

GeneRegion - Perl extension for Genome Browser Queries

=head1 SYNOPSIS

  use GeneRegion;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for GeneRegion, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited....  Yes I was.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Tim Sterne-Weiler, E<lt>timsw@localdomainE<gt>

=head1 COPYRIGHT AND LICENSE

# Copyright (C) 2012 Tim Sterne-Weiler
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software 
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
# OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

=cut
