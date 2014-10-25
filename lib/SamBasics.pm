package SamBasics; 
 
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
   mappedLength
   alignPosInRead
   totalLength
   parseOpFields
   isCIGAR
   isPCRdup
   isSecondary
   isMapped
   getStrand
   passedQC
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
   mappedLength
   alignPosInRead
   totalLength
   parseOpFields
   isCIGAR
   isPCRdup
   isSecondary
   isMapped
   getStrand
   passedQC
); 

our $VERSION = '0.01'; 



########################################################################################## 



########################################################################################### 
#                                                                                         # 
#      Subroutines (EXPORTED)                                                             # 
#                                                                                         # 
########################################################################################### 
 
# provides length of alignment in reference
# to be paired with reference position in sam entry
sub mappedLength {
  my $cigar = shift;
  my $sum = 0;
  my(@a) = split(/(?<=[MIDNSHP\=X])(?=\d)/, $cigar);
  foreach my $q (@a) {
    if($q =~ /[SIH]$/) { next; }
    my(@b) = split(/(?<=\d)(?=[MIDNSHP\=X])/, $q);
    $sum += $b[0];
  }
  return($sum);
}

# read-based version of mappedLength, 
# get the start and length of alignment in the read
sub alignPosInRead {
  my $cigar = shift;
  my $length = 0;
  my $start = 1;
  my(@a) = split(/(?<=[MIDNSHP\=X])(?=\d)/, $cigar);  
  foreach my $q (@a) {
     if($q =~ /[HDNP]$/) { next; } # reference inc
     my(@b) = split(/(?<=\d)(?=[MIDNSHP\=X])/, $q);    
     if($q =~ /S$/ and $length == 0) {
       #increment start until match
       $start += $b[0];
     } else {
       $length += $b[0] if($q !~ /S$/);
     }
  }
  return($start, $length)
}

sub isCIGAR {
  my($cigar) = shift;
  my $match = ($cigar =~ /^[MIDNSHIPX\=\d]+$/) ? 1 : 0;
  return($match);
}

# parses optional fields and returns reference to hash
# with key=>value
sub parseOpFields {
  my $samArrRef = shift;
  my %optRef;
  for(my $i=10; $i < scalar(@$samArrRef); $i++) {
    if(defined($samArrRef->[$i])) {
      my(@a) = split(/\:/, $samArrRef->[$i]);
      $optRef{$a[0]} = $a[2];
    }
  }
  return(\%optRef);
} 

sub passedQC {
  my($flag) = shift;
  return(($flag & (0x200)) ? 0 : 1);
}

sub isMapped {
  my($flag) = shift;
  return(($flag & (0x4)) ? 0 : 1);
}

sub isSecondary {
  my($flag) = shift;
  return(($flag & (0x100)) ? 1 : 0);
}

sub isPCRdup {
  my($flag) = shift;
  return(($flag & (0x400)) ? 1 : 0);
}

sub getStrand {
  my($flag) = shift;
  return(($flag & (0x10)) ? "-" : "+");
}


1; 
__END__ 

