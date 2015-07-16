package SequenceBasics; 
 
# Copyright (C) 2014 Tim Sterne-Weiler
# e-mail: tim.sterne.weiler@utoronto.ca
 
use 5.008008; 
use strict; 
use warnings; 
 
use DBI;

use FuncBasics qw(:all);
use CoordBasics qw(:all);

require Exporter; 
 
our @ISA = qw(Exporter); 
 

our %EXPORT_TAGS = ( 'all' => [ qw( 
	getSeq
        maskstr
	revComp
        gcContent
	maxRegex
	shuffleSeq
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
	getSeq
        maskstr
	revComp
        gcContent
	maxRegex
	shuffleSeq
); 

our $VERSION = '0.01'; 

#      Subroutines (EXPORTED)                                                             # 

sub revComp {
  my $seq = shift;
  $seq =~ tr/ATGCatgc/TACGtacg/;
  return(scalar reverse($seq));
}

# this function works like substr, except it masks instead.
# so maskstr("APPLE", 0, 2) returns "apPLE";
sub maskstr {
  my($seq, $offset, $length) = @_;
  return($seq) unless($length);
  my $placeHolder = ">".rand().">";
  my $sub = substr($seq, $offset, $length, $placeHolder);
  $sub = lc $sub;
  $seq =~ s/$placeHolder/$sub/g;
  return($seq);
}

sub gcContent {
  my $seq = shift;
  my $count = () = $seq =~ /[GC]/gi;
  return $count / length($seq);
}

sub getSeq {
  my($coord, $dir, $opts) = @_;
  unless(defined($opts)) { $opts = ""; }
  my($chr, $start, $stop, $strand) = parseRegion($coord);
  if($strand eq "-") { $strand = "m"; }
  my(@stdout);
  put_cmd_in("nibFrag $opts $dir/$chr.nib $start $stop $strand stdout",
            \@stdout, undef);
  shift(@stdout);
  my($seq) = "";
  foreach my $line (@stdout) {
    chomp($line);
    $seq = "$seq$line";
  }
  return($seq);
}


sub maxRegex {
  my($str, $regex) = @_;
  my($maxLen) = 0;
  my($maxMatch);
  while($str =~ /($regex)/gs) {
    my($s, $l) = ($1, length($1));
    if($l > $maxLen) {
      $maxLen = $l;
      $maxMatch = $s;
    }
  }
  return($maxMatch, $maxLen);
}

sub shuffleSeq {
  my $seq = shift;
#  randomSeedRNG(undef); 
  my $shuffled = "";
  while(length($seq) > 0) {
    my $r = int(rand(length($seq)));
    $shuffled .= substr($seq, $r, 1, "");
  }
  return($shuffled);
}

1; 
__END__ 

