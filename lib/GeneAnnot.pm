package GeneAnnot; 
 
# Copyright (C) 2013 Tim Sterne-Weiler
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
 
use 5.008008; 
use strict; 
use warnings; 
 
use IPC::Open3::Utils qw(:all);  
use DBI; 

use CoordBasics qw(:all);
use FuncBasics qw(:all);
 
require Exporter; 
 
our @ISA = qw(Exporter); 
 
# Items to export into callers namespace by default. Note: do not export 
# names by default without a very good reason. Use EXPORT_OK instead. 
# Do not simply export all your public functions/methods/constants. 
 
# This allows declaration       use GeneRegion ':all'; 
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK 
# will save memory. 
our %EXPORT_TAGS = ( 'all' => [ qw( 
   load_GFF_or_GTF
   printRefFlat
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
  load_GFF_or_GTF
  printRefFlat
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
#   $self->{"GENE_ISOS"}->{gene_id}->{transcript_id} = strand;
#   $self->{"GENE_ALIAS"}->{gene_id} = { /gene_alias/ => value };
#   $self->{"GENE_COORD"}->{gene_id} = [ chr, start, stop, strand ];
#
#   $self->{"ISO_EXON"}->{transcript_id} = { "chr:start-stop" => frame, .. }
#   $self->{"ISO_ALIAS"}->{transcript_id} = { /transcript_alias/ => value };
#   $self->{"ISO_ATTR"}->{transcript_id} = { /!(id or alias)/ => value };

sub load_GFF_or_GTF {
  my($self, $fileName) = @_;

  my(%GENE_ISOS);
  my(%GENE_ALIAS);
  my(%GENE_COORD);

  my(%ISO_EXON);
  my(%ISO_ALIAS);
  my(%ISO_ATTR);

  open(GFF, $fileName) or die "cannot open gff3 file $fileName!\n";
  my($curGene) = "";
  my($curIso) = "";
  while(my $l = <GFF>) {
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

    if($thisGene ne $curGene) { # new gene...
      unless(defined($GENE_ISOS{$thisGene})) {
	# $GENE_ISOS{gene}->{transcript} = strand;
        $GENE_ISOS{$thisGene} = { $thisTrans => $t[6] }
      }
      $GENE_ALIAS{$thisGene} = getAttribs($t[8], "gene_alias");
      $GENE_COORD{$thisGene} = [ $t[0], $t[3], $t[4], $t[6] ];
      # Don't do this again...
      $curGene = $thisGene;
    } 

    ## Add isoform to gene.
    $GENE_ISOS{$thisGene}->{$thisTrans} = $t[6];

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

    $GENE_COORD{$thisGene}->[1] = min($GENE_COORD{$thisGene}->[1], $t[3]);
    $GENE_COORD{$thisGene}->[2] = max($GENE_COORD{$thisGene}->[2], $t[4]);

  }
  close GFF;


  $self->{"GENE_ISOS"} = \%GENE_ISOS;
  $self->{"GENE_ALIAS"} = \%GENE_ALIAS;
  $self->{"GENE_COORD"} = \%GENE_COORD;

  $self->{"ISO_EXON"} = \%ISO_EXON;
  $self->{"ISO_ALIAS"} = \%ISO_ALIAS;
  $self->{"ISO_ATTR"} = \%ISO_ATTR;

}

###########################################################################################
###
###
### 
###

sub printGTF {

}

# TO ACCESS: $SELF
# #
# #   $self->{"GENE_ISOS"}->{gene_id}->{transcript_id} = strand;
# #   $self->{"GENE_ALIAS"}->{gene_id} = { /gene_alias/ => value };
# #   $self->{"GENE_COORD"}->{gene_id} = [ chr, start, stop, strand ];
# #
# #   $self->{"ISO_EXON"}->{transcript_id} = { "chr:start-stop" => frame, .. }
# #   $self->{"ISO_ALIAS"}->{transcript_id} = { /transcript_alias/ => value };
# #   $self->{"ISO_ATTR"}->{transcript_id} = { /!(id or alias)/ => value };
sub printRefFlat {
  my($self) = shift;
  
  foreach my $gene (keys %{$self->{"GENE_ISOS"}}) {
    my($chr, $start, $stop, $strand) = parseRegion($self->{"GENE_COORD"}->{$gene});
    my(@aliasKeys) = sort keys %{$self->{"GENE_ALIAS"}->{$gene}};
    my($alias);
    if(defined($aliasKeys[0])) {
      $alias = $self->{"GENE_ALIAS"}->{$gene}->{$aliasKeys[0]};
    } else {
      $alias = $gene;
    }
    foreach my $trans (keys %{$self->{"GENE_ISOS"}->{$gene}}) {
      my($exStList, $exEnList, $exFrList) = ("", "", "");
      my $exonCount = 0;
      foreach my $exonCoord (keys %{$self->{"ISO_EXON"}->{$trans}} ) {
        my($frame) = $self->{"ISO_EXON"}->{$trans}->{$exonCoord};
        my($c, $s, $e) = parseRegion($exonCoord);
        $exStList = addToBlank($exStList, $s);
        $exEnList = addToBlank($exEnList, $e);
        $exFrList = addToBlank($exFrList, $frame);
        $exonCount++;
      }

      $trans = parseEns($trans);
      $alias = parseEns($alias);
      print "1\t$trans\t$chr\t$strand\t$start\t$stop\t$start\t$stop\t";
      print "$exonCount\t$exStList\t$exEnList\t0\t$trans\tnone\tnone\t$exFrList\n";
      
    }
  }
  sub addToBlank {
    my($str, $toAdd) = @_;
    if($toAdd eq ".") { $toAdd = -1; }
    if($str eq "") { return($toAdd); }
    else {
      return("$str\,$toAdd");
    }
  }
  sub parseEns {
    my $id = shift;
    unless($id =~ /^ENS/) { return($id); }
    my(@a) = split(/\./, $id);
    return($a[0]);
  }
}

########################################################################################### 
#                                                                                         # 
#      Subroutines (BACKEND)                                                              # 
#                                                                                         # 
########################################################################################### 


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

# Copyright (C) 2012 Tim Sterne-Weiler & Jeremy Sanford 
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

