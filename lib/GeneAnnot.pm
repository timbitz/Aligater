package GeneAnnot; 
 
# Copyright (C) 2013 Tim Sterne-Weiler
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
 

our %EXPORT_TAGS = ( 'all' => [ qw( 
   load_GFF_or_GTF
   initGeneLookup
   coorAliasLookup
   printRefFlat
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
  load_GFF_or_GTF
  initGeneLookup
  coorAliasLookup
  printRefFlat
); 

our $VERSION = '0.01'; 


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
#   $self->{"GENE_CHILD"}->{gene_id}->{transcript_id} = strand;
#   $self->{"GENE_ALIAS"}->{gene_id} = { /gene_alias/ => value };
#   $self->{"GENE_COORD"}->{gene_id} = [ chr, start, stop, strand ];
#
#   $self->{"ISO_EXON"}->{transcript_id} = { "chr:start-stop" => frame, .. }
#   $self->{"ISO_ALIAS"}->{transcript_id} = { /transcript_alias/ => value };
#   $self->{"ISO_ATTR"}->{transcript_id} = { /!(id or alias)/ => value };
#   $self->{"ISO_PARENT"}->{transcript_id} = "gene_id"; 

sub load_GFF_or_GTF {
  my($self, $fileName, $ignoreDecimal) = @_;

  $ignoreDecimal = (defined($ignoreDecimal) and $ignoreDecimal) ? 1 : 0;

  my(%GENE_CHILD);
  my(%GENE_ALIAS);
  my(%GENE_COORD);

  my(%ISO_EXON);
  my(%ISO_ALIAS);
  my(%ISO_ATTR);
  my(%ISO_GENE);

  my $GFFhndl = openFileHandle($fileName);
  my($curGene) = "";
  my($curIso) = "";
  while(my $l = <$GFFhndl>) {
    chomp($l);
    if($l =~ /^#/) { next; }
    my(@t) = split(/\t/, $l);
    ### ALTER THIS LINE TO PROVIDE SUPPORT FOR CDS.. Currently I only need "exon"
    if($t[2] eq "gene" or $t[2] eq "mRNA" or $t[2] eq "CDS") { next; }
    unless($t[2] eq "exon") { next; }  #this does the same thing as the previous line.

    #some basic regex skills.
    $t[-1] =~ /gene_id[\s\"\=]+(\S+)[\s\"]+\;/;
    my($thisGene) = $1;
    $t[-1] =~ /transcript_id[\s\"\=]+(\S+)[\s\"]+\;/;
    my($thisTrans) = $1;

    $thisGene =~ s/\.\d+$//g if $ignoreDecimal;
    $thisTrans =~ s/\.\d+$//g if $ignoreDecimal;

    if($thisGene ne $curGene) { # new gene...
      unless(defined($GENE_CHILD{$thisGene})) {
	# $GENE_CHILD{gene}->{transcript} = strand;
        $GENE_CHILD{$thisGene} = { $thisTrans => $t[6] }
      }
      $GENE_ALIAS{$thisGene} = getAttribs($t[8], "gene_(name|alias)");
      $GENE_COORD{$thisGene} = [ $t[0], $t[3], $t[4], $t[6] ];
      # Don't do this again...
      $curGene = $thisGene;
    } 

    ## Add isoform to gene.
    $GENE_CHILD{$thisGene}->{$thisTrans} = $t[6];

    if(defined($ISO_EXON{$thisTrans})) {
      $ISO_EXON{$thisTrans}->{"$t[0]\:$t[3]\-$t[4]"} = $t[7];
    } else {
      $ISO_EXON{$thisTrans} = { "$t[0]\:$t[3]\-$t[4]" => $t[7] };
    }

    unless(defined($ISO_ALIAS{$thisTrans})) {
      $ISO_ALIAS{$thisTrans} = getAttribs($t[8], "transcript_alias");
    }
    unless(defined($ISO_ATTR{$thisTrans})) {
      $ISO_ATTR{$thisTrans} = getAttribs($t[8], "^(id|alias)");
    }
    unless(defined($ISO_GENE{$thisTrans})) {
      $ISO_GENE{$thisTrans} = $thisGene;
    }

    $GENE_COORD{$thisGene}->[1] = min($GENE_COORD{$thisGene}->[1], $t[3]);
    $GENE_COORD{$thisGene}->[2] = max($GENE_COORD{$thisGene}->[2], $t[4]);

  }
  close $GFFhndl;


  $self->{"GENE_CHILD"} = \%GENE_CHILD;
  $self->{"GENE_ALIAS"} = \%GENE_ALIAS;
  $self->{"GENE_COORD"} = \%GENE_COORD;

  $self->{"ISO_EXON"} = \%ISO_EXON;
  $self->{"ISO_ALIAS"} = \%ISO_ALIAS;
  $self->{"ISO_ATTR"} = \%ISO_ATTR;
  $self->{"ISO_PARENT"} = \%ISO_GENE;
}

###########################################################################################
###
###
### 
###

sub printGTF {
  ### TO DO
}

# TO ACCESS: $SELF
# #
# #   $self->{"GENE_CHILD"}->{gene_id}->{transcript_id} = strand;
# #   $self->{"GENE_ALIAS"}->{gene_id} = { /gene_alias/ => value };
# #   $self->{"GENE_COORD"}->{gene_id} = [ chr, start, stop, strand ];
# #
# #   $self->{"ISO_EXON"}->{transcript_id} = { "chr:start-stop" => frame, .. }
# #   $self->{"ISO_ALIAS"}->{transcript_id} = { /transcript_alias/ => value };
# #   $self->{"ISO_ATTR"}->{transcript_id} = { /!(id or alias)/ => value };
# #   $self->{"ISO_PARENT"}->{transcript_id} = "gene_id";
sub printRefFlat {
  my($self) = shift;

  #internal subroutines
  my $addToBlank = sub {
    my($str, $toAdd) = @_;
    if($toAdd eq ".") { $toAdd = -1; }
    if($str eq "") { return($toAdd); }
    else {
      return("$str\,$toAdd");
    }
  };
  my $parseEns = sub {
    my $id = shift;
    unless($id =~ /^ENS/) { return($id); }
    my(@a) = split(/\./, $id);
    return($a[0]);
  };

  foreach my $gene (keys %{$self->{"GENE_CHILD"}}) {
    my($chr, $start, $stop, $strand) = parseRegion($self->{"GENE_COORD"}->{$gene});
    my(@aliasKeys) = sort keys %{$self->{"GENE_ALIAS"}->{$gene}};
    my($alias);
    if(defined($aliasKeys[0])) {
      $alias = $self->{"GENE_ALIAS"}->{$gene}->{$aliasKeys[0]};
    } else {
      $alias = $gene;
    }
    foreach my $trans (keys %{$self->{"GENE_CHILD"}->{$gene}}) {
      my($exStList, $exEnList, $exFrList) = ("", "", "");
      my $exonCount = 0;
      foreach my $exonCoord (sort { $a cmp $b } keys %{$self->{"ISO_EXON"}->{$trans}} ) {
        my($frame) = $self->{"ISO_EXON"}->{$trans}->{$exonCoord};
        my($c, $s, $e) = parseRegion($exonCoord);
        $exStList = $addToBlank->($exStList, $s);
        $exEnList = $addToBlank->($exEnList, $e);
        $exFrList = $addToBlank->($exFrList, $frame);
        $exonCount++;
      }

      $trans = $parseEns->($trans);
      $alias = $parseEns->($alias);
      my $bin = binFromRangeExtended($start, $stop);
      print "$bin\t$trans\t$chr\t$strand\t$start\t$stop\t$start\t$stop\t";
      print "$exonCount\t$exStList\t$exEnList\t0\t$trans\tnone\tnone\t$exFrList\n";
      
    }
  }
}

### toGenomeCoord -- convert transcript coordinate to genome.
### returns array ref in the form of [ chrom, genome_coord, strand ];
sub toGenomeCoord {
  my($self, $tranId, $coord) = @_;
  my $geneId = $self->{"ISO_PARENT"}->{$tranId};
  return unless defined $geneId;
  my $strand = $self->{"GENE_CHILD"}->{$geneId}->{$tranId};
  return unless defined $strand;
  my(@exons) = sort { 
    ($strand eq "+") ? coorMidpoint($a) <=> coorMidpoint($b) : coorMidpoint($b) <=> coorMidpoint($a) 
                    } keys %{$self->{"ISO_EXON"}->{$tranId}};
  return unless scalar(@exons) > 0;
  my $cur = 0; # running total
  my $ret;
  foreach my $exon (@exons) {
    my($c, $s, $e) = parseRegion($exon);
    # check if coord is within this exon
    if(($e - $s) + $cur > $coord) {
      if($strand eq "+") {
        $ret = ($coord - $cur) + $s;
      } else { #minus strand, trace backwards
        $ret = $e - ($coord - $cur);
      }
      return unless defined $ret;  # RETURN
      return [$c, $ret, $strand];  # HERE
    } else {
      $cur += ($e - $s) + 1;
    }
  }
}

#
# this function sets $self->{"CHR_BIN_GENE"}->{chr}->{bin}->{gene_alias}->[coord]
#
sub initGeneLookup {
  my($self) = shift;
  my(%chrBinLookup);
  foreach my $geneId (keys %{$self->{"GENE_ALIAS"}}) {
    my $coord = $self->{"GENE_COORD"}->{$geneId};
    next unless defined $coord; # test if defined
    my($chr, $pos) = parseRegion($coord);
    my $bin = binFromRangeExtended($pos, $pos+2);
    #initialize hash structure;
    unless(defined($chrBinLookup{$chr})) {
      $chrBinLookup{$chr} = {};
    }
    unless(defined($chrBinLookup{$chr}->{$bin})) {
      $chrBinLookup{$chr}->{$bin} = {};
    }
    # now add each alias for lookup
    foreach my $attr (keys %{$self->{"GENE_ALIAS"}->{$geneId}}) {
      my $alias = $self->{"GENE_ALIAS"}->{$geneId}->{$attr};
      next unless defined $alias;
      $chrBinLookup{$chr}->{$bin}->{$alias} = $coord;
    }
  }
  $self->{"CHR_BIN_GENE"} = \%chrBinLookup;
}

# this function assumes initGeneLookup has been run already.
sub coorAliasLookup {
  my($self, $coord, $alias) = @_;
  my($chr, $pos) = parseRegion($coord);
  my $bin = binFromRangeExtended($pos, $pos+2);
  return undef unless (defined $chr and defined $bin);
  return undef unless defined $alias;
  my $lookupCoor = $self->{"CHR_BIN_GENE"}->{$chr}->{$bin}->{$alias};
  return undef unless defined $lookupCoor;
  return coorOverlap($lookupCoor, $coord);
}


#      Subroutines (BACKEND)                                                              # 


#this function collects attributes from the last column of GTF file.
# the second argument is used as a regex filter for collection.
sub getAttribs {
  my($text, $filtRegex) = @_;
  my(@txt) = split(/\;/, $text);
  my(%return);
  foreach my $i (@txt) {
    unless($i =~ /$filtRegex/) { next; }
    $i =~ /(\S+)[\s\"\=]+(\S+)[\s\"\=]+/;
    my($key, $val) = ($1, $2);
    $return{$key} = $val;
  }
  return(\%return);
}

1; 
__END__ 

