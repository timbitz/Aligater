package StatFreq; 
 
# Copyright (C) 2014 Tim Sterne-Weiler
# e-mail: tim.sterne.weiler@utoronto.ca
 
use 5.008008; 
use strict; 
use warnings; 

use FuncBasics qw(sumOfHash min max stdDev);
 
require Exporter; 
 
our @ISA = qw(Exporter); 
 
# Items to export into callers namespace by default. Note: do not export 
# names by default without a very good reason. Use EXPORT_OK instead. 
# Do not simply export all your public functions/methods/constants. 
 
# This allows declaration       use GeneRegion ':all'; 
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK 
# will save memory. 
our %EXPORT_TAGS = ( 'all' => [ qw( 
	printFreqHash
	printCatHash 
	makeFreqHash
	makeFreqHashFromFile
	initializeFreqHash
	incFreqHash
	incCatHash
	addPseudoToFreqHash
	normHashToDensity
	calcStdDevHash
	subtractHash
	additionHash
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
        printFreqHash
	printCatHash 
        makeFreqHash
	makeFreqHashFromFile
        initializeFreqHash
        incFreqHash
        incCatHash
	addPseudoToFreqHash
	normHashToDensity
	calcStdDevHash
	subtractHash
	additionHash
); 

our $VERSION = '0.01'; 



########################################################################################## 



########################################################################################### 
#                                                                                         # 
#      Subroutines (EXPORTED)                                                             # 
#                                                                                         # 
########################################################################################### 


# normalize freq distribution
# precomputed sum is optional for more complex normalizations
sub normHashToDensity {
  my($freqRef, $sum) = @_;
  $sum = sumOfHash($freqRef) unless defined($sum);
  if($sum <= 0) { return $freqRef; }  # bad.
  foreach my $k (keys %$freqRef) {
    unless($freqRef->{$k} == 0) {
      $freqRef->{$k} /= $sum;
    }
  }
  # no return value, in place norm
}

# subtract one freq from another (normalize first!!)
# this function returns freqA - freqB
# without altering either freqA or freqB directly.
sub subtractHash {
  my($freqA, $freqB) = @_;
  my %retFreq;
  foreach my $k (keys %$freqA, keys %$freqB) {
    $retFreq{$k} = 0;
  }
  foreach my $k (keys %retFreq) {
    if(defined($freqA->{$k})) {
      $retFreq{$k} += $freqA->{$k};
    }
    if(defined($freqB->{$k})) {
      $retFreq{$k} -= $freqB->{$k};
    }
  }
  return(\%retFreq);
}

# add one freq to another
# this function returns freqA + freqB
# without altering either freqA or freqB directly.
sub additionHash {
  my($freqA, $freqB) = @_;
  my %retFreq;
  foreach my $k (keys %$freqA, keys %$freqB) {
    $retFreq{$k} = 0;
  }
  foreach my $k (keys %retFreq) {
    if(defined($freqA->{$k})) {
      $retFreq{$k} += $freqA->{$k};
    }
    if(defined($freqB->{$k})) {
      $retFreq{$k} += $freqB->{$k};
    }
  }
  return(\%retFreq);
}

# This function takes a set of frequency hashes
# takes the union of all the keys, and then for each key
# calculates the std devation of all the values (s/undef/0/)
sub calcStdDevHash {
  my(@freqs) = @_;
  my %retStd;
  # initialize
  foreach my $hRef (@freqs) {
    foreach my $k (keys %$hRef) {
      $retStd{$k} = 0;
    }
  }
  foreach my $k (keys %retStd) {
    my(@vals);
    foreach my $hRef (@freqs) {
      if(defined($hRef->{$k})) {
        push(@vals, $hRef->{$k});
      } else {
        push(@vals, 0);
      }
    }
    $retStd{$k} = stdDev(@vals);
  }
  return(\%retStd);
}

# increment categorical hash
sub incCatHash {
  my($catHash, $val) = @_;
  unless(defined($catHash->{$val})) {
    $catHash->{$val} = 1;
  } else {
    $catHash->{$val}++;
  }
}

# print a cat hash
sub printCatHash {
  my($hRef, $optString, $stdDevRef) = @_;
  my $opt = "";
  if(defined($optString)) { 
    $opt = "\t$optString";
  }
  foreach my $key (sort { $a cmp $b } keys %$hRef) {
    my $curOpt = $opt;
    if(defined($stdDevRef->{$key})) {
      $curOpt = "\t$stdDevRef->{$key}$opt";
    }
    print "$key\t$hRef->{$key}$curOpt\n";
  }
}

# print a hash in numerical order
sub printFreqHash {
  my($hRef, $optString, $stdDevRef, $stdDevPrintFlag) = @_;
  my $opt = "";
  if(defined($optString)) {
    $opt = "\t$optString";
  }
  foreach my $key (sort { $a <=> $b } keys %$hRef) {
    my $curOpt = $opt;
    if($stdDevPrintFlag) {
      if(defined($stdDevRef->{$key})) {
        $curOpt = "\t$stdDevRef->{$key}$opt";
      } else {
        $curOpt = "\t0$opt";
      }
    }
    print "$key\t$hRef->{$key}$curOpt\n";
  }
}

# increment a freq hash.
sub incFreqHash {
  my($freq, $val, $binSize) = @_;
  my $adj = ($val < 0) ? 1 : 0;
  my $bin = (int($val / $binSize) - $adj) * $binSize + ($binSize / 2);
  unless(defined($freq->{$bin})) {
    $freq->{$bin} = 1;
  } else {
    $freq->{$bin}++;
  }
}

# use this function to add pseudo counts (or 0)
# to every bin within the range of the freq dist
# you can also specify the start and end values
# no return.
sub addPseudoToFreqHash {
  my($freq, $binSize, $pseudo, $beginVal, $endVal) = @_;
  unless(defined($pseudo)) { $pseudo = 0; } # default
  unless(defined($beginVal)) { $beginVal = "inf"; } # upper def
  unless(defined($endVal)) { $endVal = "-inf"; } # lower def
  my(@sortKeys) = sort { $a <=> $b } keys %$freq;
  my($start) = min($sortKeys[0], $beginVal);
  my($stop) = max($sortKeys[-1], $endVal);
  for(my $i=$start; $i <= $stop; $i += $binSize) {
    unless(defined($freq->{$i})) {
      $freq->{$i} = $pseudo;
    } else {
      $freq->{$i} += $pseudo;
    }
  }
}

# calculate freq distribution on filehandle/column
sub makeFreqHashFromFile {
  my($fHndl, $col, $binSize, $offset, $maxVal) = @_;
  my $hRef = initializeFreqHash($binSize, $offset, $maxVal);
  while(my $l = <$fHndl>) {
    chomp $l;
    my(@a) = split(/\t/, $l);
    my $bin = int($a[$col] / $binSize) * $binSize;
    $hRef->{$bin}++;    
  }
  return $hRef;
}


# calculate freq distribution on array.
sub makeFreqHash {
  my($arrRef, $binSize, $offset, $maxVal) = @_;
  my $hRef = initializeFreqHash($binSize, $offset, $maxVal);
  for(my $i=0; $i < scalar(@$arrRef); $i++) {
    my $bin = int($arrRef->[$i] / $binSize) * $binSize;
    $hRef->{$bin}++;
  }
  return $hRef;
}

# Use to initialize freq distribution
sub initializeFreqHash {
  my($size, $off, $end, $set) = @_;
  unless(defined($set)) { $set = 0; }
  my %hash;
  for(my $i=$off; $i <= $end; $i += $size) {
    $hash{$i} = $set;
  }
  return(\%hash);
}

1; 
__END__ 

