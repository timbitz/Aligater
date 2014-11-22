package GeneElement; 
 
# Copyright (C) 2014 Tim Sterne-Weiler
# e-mail: tim.sterne.weiler@utoronto.ca
 
use 5.008008; 
use strict; 
use warnings; 
 
use DBI; 

use CoordBasics qw(:all);
use FuncBasics qw(:all);
use hgSQLBasics qw(binFromRangeExtended);
 
require Exporter; 
 
our @ISA = qw(Exporter); 
 
# Items to export into callers namespace by default. Note: do not export 
# names by default without a very good reason. Use EXPORT_OK instead. 
# Do not simply export all your public functions/methods/constants. 
 
# This allows declaration       use GeneRegion ':all'; 
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK 
# will save memory. 
our %EXPORT_TAGS = ( 'all' => [ qw( 
   load_BED
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
  load_BED
); 

our $VERSION = '0.01'; 



########################################################################################## 

########################################################################################### 
#                                                                                         # 
#      Subroutines (EXPORTED)                                                             # 
#                                                                                         # 
########################################################################################### 


###  THE OBJECT CONSTRUCTOR!!
sub new {
  my($class) = shift;
  my($self) = {};  # SEE BELOW FOR USAGE>

  bless ($self, $class);
  return($self);
}

# THIS FUNCTION FILLS IN THE $SELF HASH AS FOLLOWS:
#
# TO ACCESS: $SELF
#
#   $self->{"CHROM"}->{chrom}->{bin}->{coord} = name;

sub load_BED {
  my($self, $fileName) = @_;

  my(%CHROM);
  my(%COORD);

  my $BEDhndl = openFileHandle($fileName);
  while(my $l = <$BEDhndl>) {
    chomp($l);
    next if($l =~ /^#/);
    my(@t) = split(/\t/, $l);
    # ALTER THIS LINE TO PROVIDE SUPPORT FOR CDS..
    die "GeneElement.pm Input $.: invalid BED format!" unless defined $t[0] and defined $t[1] and defined $t[2];
    $CHROM{$t[0]} = {} and $COORD{$t[0]} = {} unless defined($CHROM{$t[0]});
    die "GeneElement.pm Input $.: invalid BED format col 1 and 2 $t[1] and $t[2] must be integers" unless isInt($t[1]) and isInt($t[2]);
    my $bin = binFromRangeExtended($t[1], $t[2]);
    $CHROM{$t[0]}->{$bin} = {} unless defined($CHROM{$t[0]}->{$bin});
    $COORD{$t[0]}->{$bin} = [] unless defined($COORD{$t[0]}->{$bin});
    my $coord = (defined($t[5]) and $t[5] =~ /[+-]/) ? [$t[0],$t[1],$t[2],$t[5]] : [$t[0],$t[1],$t[2]];
    my $name = (defined($t[3])) ? $t[3] : $fileName;
    $CHROM{$t[0]}->{$bin}->{coorToString($coord)} = $name;
    push(@{$COORD{$t[0]}->{$bin}}, $coord);
  }
  close $BEDhndl;

  $self->{"COORD"} = \%COORD;
  $self->{"CHROM"} = \%CHROM;
}

## This function looks for overlap of the internal BED object with a given coordinate.
## as per parseRegion, strand is not used.
sub bedOverlap {
  my($self, $coord) = @_;
  my($chr, $start, $end, $strand) = parseRegion($coord);
  my $bin = binFromRangeExtended($start, $end);
#  print "$bin\n";
  my(@results);
  foreach my $coordRef (@{$self->{"COORD"}->{$chr}->{$bin}}) {
#    print coorToString($coordRef)."\t".coorToString([$chr, $start, $end, $strand])."\n";
    if(coorOverlap($coordRef, [$chr, $start, $end, $strand])) {
      my $str = coorToString($coordRef);
      my $name = $self->{"CHROM"}->{$chr}->{$bin}->{$str};
      push(@results, "$str\t$name");
    }
  }
  return((scalar(@results) > 0) ? \@results : undef);
}

###########################################################################################
###
###
### 
###

sub printBED {
  ### TO DO
}


1; 
__END__ 
