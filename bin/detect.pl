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
use POSIX qw(strftime);
use Digest::MD5 qw(md5_hex md5_base64);

use FindBin;
use lib "$FindBin::Bin/../lib";

use Getopt::Long;

use GeneAnnot; # oop
use GeneElement;
use SamBasics qw(:all);
use FuncBasics qw(isInt shove openFileHandle);
use CoordBasics qw(coorOverlap parseRegion);
use SequenceBasics qw(maskstr revComp);

our $COORDEXPAND = 100; # this is the buffer range for overlapping genomic loci

our $STRAND_SPECIFIC = 1; # transcriptome mapping.-->SETTING TO 0 IS NOT RECOMMENDED!
our $HYBRID_PENALTY = -24; # this should be optimized.
our $ANTISENSE_PENALTY = -24; #so should this  #DEPRECATED VARIABLE 2/2015

our %GENEFAM;  # load this from anno/Species.gene_families.txt;
our $GENEANNO;  # GeneAnnot object if --gtf=s is used to specify genome coordinates
our $LOCALANNO; # GeneElement object if --rmsk=s is used to specify a repeat masker bed track

# INITIALIZE
my $path = abs_path($0);
$0 =~ s/^.*\///;
$path =~ s/\/$0$//;

# set default output prefix based on timestamp
my $outputCore = strftime 'Output_%F_%H.%M.%S', localtime;

my $base64Flag = 0;
my $suppressAlnFlag = 0;

# this variable suppresses chimeric output in favor of non-chimeric reads
# in the same format.
my $nonChimeraFlag = 0; 

my $gtfFile; #undef by default
my $geneFamFile;
my $rmskFile;

GetOptions("o=s" => \$outputCore, 
           "base" => \$base64Flag,
           "noaln" => \$suppressAlnFlag, 
	   "nochim" => \$nonChimeraFlag,
           "noanti" => \$STRAND_SPECIFIC, #TODO
           "gfam=s" => \$geneFamFile,
           "gtf=s" => \$gtfFile,
           "rmsk=s" => \$rmskFile,
           "pen=i" => \$HYBRID_PENALTY
);

# make sure this is negative.
$HYBRID_PENALTY *= -1 if($HYBRID_PENALTY > 0);

my $nonHybHndl;
my $hybHndl; # alignment output filehandles.

unless($suppressAlnFlag) {
  open($nonHybHndl, "| samtools view -bS - > $outputCore.std.bam") or die "Can't write to $outputCore.bam";
  open($hybHndl,"| samtools view -bS - > $outputCore.lig.bam") or die "Can't write to $outputCore.chimera.bam";
}
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

## Load GTF/GFF file if defined;
if(defined($gtfFile)) {
  $GENEANNO = new GeneAnnot;
  $GENEANNO->load_GFF_or_GTF($gtfFile);
  $GENEANNO->initGeneLookup();
}
# done loading GTF;

## Load BED file for RepeatMasker if defined;
if(defined($rmskFile)) {
  $LOCALANNO = new GeneElement;
  $LOCALANNO->load_BED($rmskFile);
}
# done loading BED;

# load gene family file if possible
if(defined($geneFamFile)) {
  my $geneHndl = openFileHandle($geneFamFile);
  while(my $l = <$geneHndl>) {
    next if $l =~ /^#/;
    chomp $l;
    my(@a) = split(/\t/, $l);
    explode "Improper gene family file format!" unless defined($a[0]) and defined($a[1]);
    $GENEFAM{$a[1]} = $a[0];
    for(my $i=1; $i < 1000; $i++) { #add pseudogene symbols
      $GENEFAM{"$a[0]$i"."P"} = $a[0];
      $GENEFAM{"$a[0]-$i"."P"} = $a[0];
      $GENEFAM{"$a[0]-$i"} = $a[0];
    }
  }
  close $geneHndl;
}
# done loading gene family file

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
  if($l =~ /^(#|\@(SQ|HD|RG|PG|CO))/) { #header line
    unless($suppressAlnFlag) {
      print $hybHndl $l;
      print $nonHybHndl $l;
    }
    next;
  }
  chomp($l);
  my(@a) = split(/\t/, $l); # split sam
  next unless isMapped($a[1]); # ignore unmapped reads 
  explode "Invalid SAM format!\n" unless(defined($a[0]) and isInt($a[1]));

  # substitute current readname for cleaner md5
  $a[0] = ($base64Flag) ? md5_base64($a[0]) : md5_hex($a[0]);
  $a[0] = substr($a[0], 0, 16);  # 10 ** 24 should be enough
   
  if($a[0] ne $curRead and $curRead ne "") {
    # before moving to a new read, process current read's alignments
    # from $alnHash

    # process algorithm here...
    processAlignRec($alnHash, $sHash, $eHash);

    #output best alignment here.
    my(@alnKeys) = sort { alignCmp($alnHash->{$b}, $alnHash->{$a}) } keys %$alnHash;
    my($bestK) = $alnKeys[0];
    if(defined $bestK) {
      if($bestK =~ /\:/) { # chimeric read
        my $hybrid = getHybridFormat(\@alnKeys, $alnHash) unless($nonChimeraFlag);
        my(@alns) = split(/\:/, $bestK);
        foreach my $al (@alns) {
          my $aRef = $alnHash->{$al};
          my $samOutput = join("\t", @$aRef[3 .. $#$aRef]); #TODO if $gtfFile convert sam to genomic
          print $hybHndl "$samOutput\n" unless $suppressAlnFlag;
        }
      } else { #non-chimeric read
        my $UNUSED = getHybridFormat(\@alnKeys, $alnHash) if($nonChimeraFlag);
        my $aRef = $alnHash->{$bestK};
        my $samOutput = join("\t", @$aRef[3 .. $#$aRef]);
        print $nonHybHndl "$samOutput\n" unless $suppressAlnFlag; # best non-chimeric alignments
      }
    }
    # re-initialize for new read
    $curRead = $a[0];
    $curCount = 1;
    ($alnHash, $sHash, $eHash, $seHash) = ({}, {}, {}, {}); # empty hash refs..
  }

  # process current read.
  explode "Invalid CIGAR format!" unless isCIGAR($a[5]); 
  my $strand = getStrand($a[1]); 
  #discard if - strand and strand specific?
  next if($STRAND_SPECIFIC and $strand eq "-");
  my($start, $len) = alignPosInRead($a[5]);
  my $end = $start + $len;
  if($strand eq "-") { #if - strand, reverse position
    $start = length($a[9]) - $end; # TODO TEST
    $end = $start + $len;
  }
  my $optHash = parseOpFields(\@a); 
  my $aScore = $optHash->{"AS"};
  #my $numMiss = $optHash->{"NM"};   

  # if same start/end exists.. take best;
  my $curAln = shove([$aScore, $start, $len], @a);  #???

  # add start and end records.
  $sHash->{$start} = shove($sHash->{$start}, $curCount);
  $eHash->{$end} = shove($eHash->{$end}, $curCount);

  # record full alignment
  $alnHash->{$curCount} = $curAln; 

  # increment for next read;
  $curCount++;
  $curRead = $a[0];
}

unless($suppressAlnFlag) {
  close $nonHybHndl;
  close $hybHndl;
}

#######################################################
#						      #
################# BEGIN SUBROUTINES ###################
#						      #
#######################################################

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
        my $antiPenalty = ($alnHash->{$a}->[4] ne "hybrid" and getStrand($alnHash->{$a}->[4]) eq "-") ? $ANTISENSE_PENALTY : 0;
        $antiPenalty += ($alnHash->{$b}->[4] ne "hybrid" and getStrand($alnHash->{$b}->[4]) eq "-") ? $ANTISENSE_PENALTY : 0;
        $score = $alnHash->{$a}->[0] + $alnHash->{$b}->[0] + $HYBRID_PENALTY + $antiPenalty;
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

# this function outputs the `.lig` file format...
sub getHybridFormat {
  my($alnKeys, $alnHash) = @_;
  my $hyb = $alnKeys->[0]; #sorted by alignment best key TODO analyze uniqueness.
  explode "getHybridFormat ERROR!\n" unless (defined($hyb) and defined($alnHash));
  my(@a) = split(/\:/, $hyb);
  my $alpha = join("", ("A".."Z"));
  my(%used);
  # set variables as alignment struct for substitutions.
  my $charStruct = $hyb;
  my $geneSymStruct = $hyb;
  my $ensTranStruct = $hyb;
  my $ensGeneStruct = $hyb;
  my $biotypeStruct = $hyb;
  my $repNameFamStruct = $hyb;
  my $repClassStruct = $hyb;

  my $refPositions = $hyb;
  $refPositions =~ s/:/,/g;
  my $genPositions = $hyb;
  $genPositions =~ s/:/,/g;
  my $rmskPositions = $hyb;
  $rmskPositions =~ s/:/,/g;
  my $alnLengths = $refPositions;

  my $ligSitePos = "";
  my $readName = $alnHash->{$hyb}->[3];
  my $alnScore = $alnHash->{$hyb}->[0];

  # chimera quality here
  my $ligq = chimeraUniqueness($alnKeys, $alnHash);

  my $readSeq = $alnHash->{1}->[12];  # get raw read sequence in forward orientation
  $readSeq = revComp($readSeq) if (getStrand($alnHash->{1}->[4]) eq "-"); # TEST TODO

  # convert chimeric read structures from gene symbol-> segment character
  for(my $i=0; $i < scalar(@a); $i++) {
    my $id = $a[$i];
    my $char = substr($alpha, scalar keys %used, 1);
    my($ensTran, $ensGene, $geneSym, $biotype) = split(/\_/, $alnHash->{$id}->[5]);
    my $refPos = $alnHash->{$id}->[6]; 
    my $readPos = $alnHash->{$id}->[1];
    my $length = $alnHash->{$id}->[2];
    my $strand = getStrand($alnHash->{$id}->[4]); # TODO TEST
    my $rmskPos = "NA";

    # set genomic position if possible
    my $coord = $GENEANNO->toGenomeCoord($ensTran, $refPos) if defined($GENEANNO);
    my($genomeChr, $genomePos, $genomeRan) = parseRegion($coord);
    my $genomeCoord = "NA";
    my($repName, $repFamily, $repClass) = ("NA", "NA", "NA");
    if(defined($genomePos)) {  # if genome position exists...
      if(defined($LOCALANNO)) {
        # as of 3/2015 we look +25bp from the genomePos to check for repeat overlap.
        my $resRef = $LOCALANNO->bedOverlap([$genomeChr, $genomePos, $genomePos+25, $genomeRan]);
        # get name of local bed entry if exists
        my($resCoord, $local) = split(/\t/, $resRef->[0]) if defined($resRef);
        # NOTE: ordered based on selected fields at UCSC!!
        if(defined($resCoord) and defined($local)) { # if there is overlap results from rmsk
          ($repName, $repClass, $repFamily) = split(/\,/, $local);
          my $resGenomeStrand = ($genomeRan eq "+") ? "-" : "+"; # reverse strand for..
          # check if the repeat name/class is on the opposite strand
          # if so set as inverted repeat `IR`
          if($resCoord =~ /\:$resGenomeStrand$/) {
            $repName = "IR-$repName";
            $repClass = "IR-$repClass";
          }
          # lets now get the ligation position within the local rmsk element.
          my(undef, $rmskStart, $rmskEnd, undef) = parseRegion($resCoord);
          if($genomeRan eq "+") {
            $rmskPos = $genomePos - $rmskStart;  #TODO check off by one calculation here
          } else {
            $rmskPos = $rmskEnd - $genomePos;  # same here.
          } 
        }
      }
      # if - strand alignment, reverse genomeRan
      $genomeRan = (defined($genomeRan) and $genomeRan eq "+") ? "-" : "+" if($strand eq "-"); # quaternary operator? ;-)
      $genomeCoord = "$genomeChr\:$genomePos\:$genomeRan";    
    } # else can't get genomic locus
    
    $repFamily = ($repName eq $repFamily) ? "" : "\_$repFamily"; #

    # set previously used char for same gene symbol
    if(defined($used{$geneSym})) {
      $char = $used{$geneSym};
    }
    $char = "-$char" if $strand eq "-"; # mark the antisense structure; #TODO TEST

    # substitute hyb number for char or sym etc;
    $charStruct    =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$char/g;
    $geneSymStruct =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$geneSym/g;
    $ensTranStruct =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$ensTran/g;
    $ensGeneStruct =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$ensGene/g;
    $biotypeStruct =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$biotype/g;

    $repNameFamStruct =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$repName$repFamily/g;
    $repClassStruct   =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$repClass/g;

    # push alignment start position in reference.
    $refPositions  =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$refPos/;
    $genPositions  =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$genomeCoord/;
    $rmskPositions =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$rmskPos/;
    $alnLengths    =~ s/\b(?<![\-\.])$id(?![\-\.])\b/$length/;

    # make read sequence structure;
    if($i == 0) { # first alignment
      $readSeq = maskstr($readSeq, 0, $readPos) if ($readPos > 1);
    } elsif($i == $#a) { #last alignment
      $readSeq = maskstr($readSeq, $readPos+$length - 1, length($readSeq) - ($readPos+$length) + 1);
    }
    #insert _ at ligation site. 
    substr($readSeq, $readPos - 1, 0, "_") unless($i == 0);

    $used{$geneSym} = $char;
  }
  # get hybrid code and gene family structure.
  my($hybCode, $familyStruct) = getHybridCode($charStruct, $geneSymStruct, $repNameFamStruct, $genPositions); 
  $charStruct =~ tr/B-Z/A/ if $hybCode eq "S"; # if we changed the hybCode, alter charStruct to match.
  ## Main output format...
  print "$hybCode\t$charStruct\t$hyb\t$geneSymStruct\t$ensGeneStruct\t$ensTranStruct\t$biotypeStruct\t$repNameFamStruct";
  print "\t$repClassStruct\t$readName\t$readSeq\t$alnScore\t$ligq\t$refPositions\t$rmskPositions\t$alnLengths\t$genPositions\n";
}

# I = putative inter-molecular, P = paralogous intra-molecular, S = intra-molecular, A = sense-antisense
sub getHybridCode {
  my($chars, $genes, $repNames, $genomePos) = @_;
  my(@c) = split(/\:/, $chars);
  my(@g) = split(/\:/, $genes); 
  my(@r) = split(/\:/, $repNames);
  my $code = "S"; # intra molecular by default.
  my(@prefix) = @g; #copy array
  # try to get gene families..
  for(my $i=0; $i < scalar(@prefix); $i++) {
    if(defined($GENEFAM{$prefix[$i]})) {  #subst with gene family
      $prefix[$i] = $GENEFAM{$prefix[$i]};
    } else {  # gene fam not identified, try symbol prefix
      my(@c) = split(/[\d\-]/, $prefix[$i]);
      $prefix[$i] = $c[0]; 
    }
  }
  # now test if it is gene family intra-molecular
  my $geneFamStruc = join(":", @prefix);  
  my $testStruc = $geneFamStruc;
  for(my $i=0; $i < scalar(@prefix); $i++) {
    my $fam = $prefix[$i];
    $testStruc =~ s/\b$fam\b/$i/g;
    $repNames =~ s/\b$r[$i]\b/$i/g unless $r[$i] eq "NA";
  }
  if($chars =~ /B/) {  #possibly inter-molecular
    $code = ($testStruc =~ /1/) ? "I" : "P";  # set code 
    # set P if both are exact same repeat from --rmsk
    $code = "P" if ($code eq "I" and $repNames =~ /^[0\:]+$/);
  }
#  $code = "A" if ($chars =~ /\-/); # antisense; #deprecated 1/2015
  # check overlap of genomePos
  if($genomePos) {
    my(@pos) = split(/\,/, $genomePos);
    my $overlap = 0;
    my $compNum = 0;
    my $antisense = 0;
    for(my $i=0; $i < scalar(@pos) - 1; $i++) {
      next if $pos[$i] eq "NA";
      my($aChr, $aPos, $aRan) = split(/\:/, $pos[$i]);
      my $aCoord = [$aChr, $aPos-$COORDEXPAND, $aPos+$COORDEXPAND, $aRan];
      my $aAlias = $prefix[$i];
      for(my $j = $i+1; $j < scalar(@pos); $j++) {
        next if $pos[$j] eq "NA";
        $compNum++;
        my($bChr, $bPos, $bRan) = split(/\:/, $pos[$j]);
        my $bCoord = [$bChr, $bPos-$COORDEXPAND, $bPos+$COORDEXPAND, $bRan];
        my $bAlias = $prefix[$j];
        # check if there is an equivalent alias at each locus.
        ($overlap++ and print STDERR "\n$aAlias,$bAlias\n\n") if $GENEANNO->coorAliasLookup($bCoord, $aAlias);
        ($overlap++ and print STDERR "\n$aAlias,$bAlias\n\n") if $GENEANNO->coorAliasLookup($aCoord, $bAlias);
        # check if coordinates overlap
        $overlap++ if coorOverlap($aCoord, $bCoord);
        $antisense++ if $aRan ne $bRan;
      }
    }
    # deprecated 1/2015 and $code ne "A"
    if($overlap == $compNum and $overlap > 0) { # then this is the same locus
      $code = $antisense ? "A" : "S"; #TODO TEST
    }
  }
  return($code, $geneFamStruc);
}

# this function's purpose is to determine how much better the maximal alignment score is
# from the next best alignment score, and how many if any alignments have the equivalnt score
# as the maximum.  #It returns this as two scalars ( indNum, nextDiff )
sub chimeraUniqueness {
  my($sortAlnArr, $alnHash) = @_;
  my $MINDIFF = 20;
  my $indNum = 0;  # keep track of the number of alignments with the same maximal score
  my $bestScore = $alnHash->{ $sortAlnArr->[0] }->[0]; # record best score.
  my $nextDiff = $bestScore; # find the difference between the max score and the next best score.
  my(@numBest); # this is filled with hashes to record the number of best alns for each mapping segment
  my(@segNum) = split(/\:/, $sortAlnArr->[0]);
  for(my $i = 0; $i < scalar(@segNum); $i++) {
    push(@numBest, {}); # push empty hash
  }
  my $singles = 0;
  for(my $i = 0; $i < scalar(@$sortAlnArr); $i++) {
    my $k = $sortAlnArr->[$i];
    #$uniqNums = "$uniqNums\n$k\t$alnHash->{$sortAlnArr->[$i]}->[5] ";  -- For debugging only.
    my $curScore = $alnHash->{$k}->[0];
    if($curScore >= $bestScore - $MINDIFF or !defined($curScore)) {
      $indNum++;
      my(@alnKeys) = split(/\:-\:/, $alnHash->{$sortAlnArr->[$i]}->[5]); # add the keys to the numBest record
      ($singles++ and next) if(scalar(@alnKeys) == 1); # IF non-hybrid alignment!!
      for(my $keyIt = 0; $keyIt < scalar(@alnKeys); $keyIt++) {
        my(@db) = split(/\_/, $alnKeys[$keyIt]);
        $numBest[$keyIt]->{$db[2]} = "";
      }
      next;
    } else {
      $curScore = 0 if(!defined($curScore)); # just in case.
      $nextDiff = $bestScore - $curScore;
      last;
    }
  }
  my $uniqNums = ""; # combine the number of mappings for each segment for printing
  foreach my $hsh (@numBest) {
    $uniqNums = "$uniqNums," if ($uniqNums ne "");
    $uniqNums = "$uniqNums" . scalar(keys %$hsh);
  }
  $indNum = "$indNum\($singles\)" if($singles);
  return("$uniqNums\>$indNum\>$nextDiff");
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

__END__
