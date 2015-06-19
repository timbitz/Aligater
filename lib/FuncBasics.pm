package FuncBasics; 
 
# Copyright (C) 2013 Tim Sterne-Weiler
# e-mail: tim.sterne.weiler@utoronto.ca
 
use 5.008008; 
use strict; 
use warnings; 

use Scalar::Util qw(looks_like_number);
 
require Exporter; 
 
our @ISA = qw(Exporter); 
 
# Items to export into callers namespace by default. Note: do not export 
# names by default without a very good reason. Use EXPORT_OK instead. 
# Do not simply export all your public functions/methods/constants. 
 
# This allows declaration       use GeneRegion ':all'; 
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK 
# will save memory. 
our %EXPORT_TAGS = ( 'all' => [ qw( 
       max
       min
       shove
       isInt
       maxOfArray
       minOfArray
       uniqArray
       uniqArrayCnt
       sumOfArray
       sumOfHash
       intersect
       mean
       stdDev
       openFileHandle
       randomSeedRNG
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
       max
       min
       shove
       isInt
       maxOfArray
       minOfArray
       uniqArray
       uniqArrayCnt
       sumOfArray
       sumOfHash	
       intersect
       mean
       stdDev
       openFileHandle
       randomSeedRNG
); 

our $VERSION = '0.01'; 



########################################################################################## 



########################################################################################### 
#                                                                                         # 
#      Subroutines (EXPORTED)                                                             # 
#                                                                                         # 
########################################################################################### 

sub max { 
  my($a, $b) = @_; 
  return($a > $b ? $a : $b ); 
} 
 
sub min { 
  my($a, $b) = @_; 
  return($a < $b ? $a : $b ); 
} 
 
# will push $val to the end of $aRef
# whether or not $aRef is defined
sub shove {
  my $aRef = shift;
  my(@val) = @_;
  if(ref($aRef) eq 'ARRAY') {
    push(@$aRef, @val);
  } else {
    $aRef = [];
    push(@$aRef, @val);
  }
  return($aRef);
}

sub isInt {
  my $n = shift;
  unless(defined($n)) { return 0; }
  unless(looks_like_number($n)) { return 0; }
  my $ret = ($n == int($n)) ? 1 : 0;
  return $ret;
}

# returns the index of the max point
sub maxOfArray { 
  my(@inp) = @_; 
  my(@index) = sort { $inp[$b] <=> $inp[$a] } 0 .. $#inp; 
  return($index[0], $inp[$index[0]]); 
} 

# returns the index of the min point
sub minOfArray { 
  my(@inp) = @_; 
  my(@index) = sort { $inp[$a] <=> $inp[$b] } 0 .. $#inp; 
  return($index[0], $inp[$index[0]]); 
}

sub uniqArray {
  my(@inp) = @_;
  my(%uniq);
  foreach my $i (@inp) { $uniq{$i} = 1; }
  return(keys %uniq)
}

sub uniqArrayCnt {
  my(@inp) = @_;
  my(%uniq);
  foreach my $i (@inp) { $uniq{$i} = 1; }
  return(scalar keys %uniq)
}

sub sumOfHash {
  my $hRef = shift;
  my $sum = 0;
  foreach my $k (keys %$hRef) {
    if(defined($hRef->{$k})) {
      $sum += $hRef->{$k};
    }
  }
  return $sum;
}

sub sumOfArray {
  my $aRef = shift;
  my $sum = 0;
  foreach my $el (@$aRef) {
    if(defined($el)) {
      $sum += $el;
    }
  }
  return $sum;
}

## set functions
sub intersect {
  my($aRefOne, $aRefTwo) = @_;
  # O(nlogn)2
  my(@sOne) = sort { $a cmp $b } @{$aRefOne};
  my(@sTwo) = sort { $a cmp $b } @{$aRefTwo};
  my(@res);
  ## BESIDES INITIAL SORT, ALGORITHM IS O(2n)
  my($o, $t) = (0, 0); # iterators for one and two
  while($o <= $#sOne and $t <= $#sTwo) {
    my $cmpOT = $sOne[$o] cmp $sTwo[$t];
    if($cmpOT == 0) {
      push(@res, $sOne[$o]);
      $o++;
      $t++;
    } elsif($cmpOT < 0) {
      $o++;
    } elsif($cmpOT > 0) {
      $t++;
    }
  }
  return(\@res);
}
#

# sample std deviation calculation
sub stdDev {
  my(@inp) = @_;
  my $ave = mean(@inp);
  my(@err);
  foreach my $i (@inp) {
    push(@err, ($i - $ave)**2);
  }
  my $total = 0;
  foreach my $i (@err) { $total += $i; }
  my($retVal) = sqrt($total / (scalar(@err)-1)); 
  return defined($retVal) ? $retVal : 0;
}

sub mean {
  my(@inp) = @_;
  my($len) = scalar(@inp);
  if($len == 0) { return undef; }
  my($total) = 0;
  foreach my $i (@inp) { $total += $i; }
  return($total / $len);
}


sub openFileHandle {
  my $fName = shift;
  my $handle;
  if($fName =~ /\.gz$/) {
    open($handle, "gunzip -c $fName |") or die "Cannot open fileA=$fName\n";
  } elsif($fName =~ /\.bam/) {
    open($handle, "samtools view $fName |") or die "Cannot open file=$fName\n";
  } else {
    open($handle, $fName) or die "Cannot open file=$fName\n";
  }
  return($handle);
}

sub randomSeedRNG {
  my $seed = shift;
  my $ret;
  if(defined($seed)) {
    srand($seed);  #pretty much useless to supply the seed.  
  } else {
    $ret = time ^ $$ ^ unpack "%L*", `ps axww | gzip`;
    srand($ret);
  }
  return $ret;
}

1; 
__END__ 


