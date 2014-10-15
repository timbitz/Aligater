package SequenceBasics; 
 
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

use FuncBasics qw(:all);
use CoordBasics qw(:all);

require Exporter; 
 
our @ISA = qw(Exporter); 
 
# Items to export into callers namespace by default. Note: do not export 
# names by default without a very good reason. Use EXPORT_OK instead. 
# Do not simply export all your public functions/methods/constants. 
 
# This allows declaration       use GeneRegion ':all'; 
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK 
# will save memory. 
our %EXPORT_TAGS = ( 'all' => [ qw( 
	getSeq
	revComp
        gcContent
	maxRegex
	shuffleSeq
) ] ); 
 
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } ); 

our @EXPORT = qw( 
	getSeq
	revComp
        gcContent
	maxRegex
	shuffleSeq
); 

our $VERSION = '0.01'; 



########################################################################################## 



########################################################################################### 
#                                                                                         # 
#      Subroutines (EXPORTED)                                                             # 
#                                                                                         # 
########################################################################################### 

sub revComp {
  my $seq = shift;
  $seq =~ tr/ATGCatgc/TACGtacg/;
  return(scalar reverse($seq));
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

=head1 NAME 

GeneRegion - Perl extension for Genome Browser Queries 

=head1 SYNOPSIS 

  use BasicFunc; 
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

