package hgSQLBasics;

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

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use GeneRegion ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
	useHgLocal
	resetDatabase
	singleSQLquery
	singleCoordSQLquery
	prepareCoordSQLquery
	prepareSmallCoordSQLquery
	executeCoordSQLquery
	executeExactCoordSQLquery
        binFromRangeExtended
	allBinsFromRangeExtended
        loadUserPassFromHgConf
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
        useHgLocal
        resetDatabase
        singleSQLquery
        singleCoordSQLquery
        prepareCoordSQLquery
	prepareSmallCoordSQLquery
        executeCoordSQLquery
        executeExactCoordSQLquery
	binFromRangeExtended
        allBinsFromRangeExtended
        loadUserPassFromHgConf
);

our $VERSION = '0.01';



##########################################################################################

our $GENOME = "hg19";

#our $REGION = undef;  # SET THIS MANUALLY BEFORE CALLING EXPORTED FUNCTIONS

our $MAXBIN;

###########################################################################################
#                                                                                         #
#       Database Variables                                                                #
#                                                                                         #
###########################################################################################

our($db_host)="genome-mysql.cse.ucsc.edu";
our($db_user)="genomep";
our($db_passwd)="password";
our($db)=$GENOME;

our($verbose) = 0;

###  Set up DBI Database....
#our($dbh) = DBI->connect("DBI:mysql:$db;host=$db_host",$db_user, $db_passwd, { RaiseError => 1 } )
#    	         or die ( "Couldn't connect to database: " . DBI->errstr );
our($dbh);

###########################################################################################
#                                                                                         #
#      Subroutines (EXPORTED)                                                             #
#                                                                                         #
###########################################################################################

sub useHgLocal {
  my($user, $pass, $mydb, $verbose) = @_;
  $db_host = "localhost";
  unless(defined($user) and defined($pass)) {
    die "No User or Pass given to useHgLocal!\n";
  } else {
    $db_user = $user;
    $db_passwd = $pass;
  }
  if(defined($mydb)) { $db = $mydb; }
  resetDatabase($db, $verbose);
}

# reset the database to "___"
sub resetDatabase {
  my($ucscGenome, $verbose) = @_;
  if($verbose) { DBI->trace(2); } 
  $dbh = DBI->connect("DBI:mysql:$ucscGenome;host=$db_host",$db_user, $db_passwd, { RaiseError => 1 } )
             or die ( "Couldn't connect to database: " . DBI->errstr );
}

# recursive function to load ~/.hg.conf
sub loadUserPassFromHgConf {
  my($prefix, $filename, $user, $pass) = @_;  # user/pass should be undef, used for recursion
  # end if defined
  if(defined($user) and defined($pass)) {
    return($user, $pass);
  } else {
    # else try to find them 
    my $hndl = openFileHandle("$prefix$filename");
    my @lines = <$hndl>;
    foreach my $l (@lines) {
      next if($l =~ /^#/);
      chomp($l);
      if($l =~ /^include\s(.*)/) {
        # recursively open included file... look for user/pass there
        ($user, $pass) = loadUserPassFromHgConf($prefix, $1, $user, $pass);
      }
      if($l =~ /db\.user\=(.*)/) {
        $user = $1;
      } elsif($l =~ /db\.password=(.*)/) {
        $pass = $1;
      }
      # end criterion
      if(defined($user) and defined($pass)) {
        return($user, $pass);
      }
    }
  }
  # if we are here then user/pass couldn't be found
  die "Error: Cannot find User/Pass in ~/.hg.conf!\n";
}

# RUN SQL query, RETURN RESULTS ARRAYREF;
sub singleSQLquery {
  my($isoQ, $return, $tables, $where) = @_;
  unless(defined($isoQ)) {
    $isoQ = "SELECT $return FROM $tables WHERE $where;";
  }
  if($verbose) { print STDERR "$isoQ\n"; }
  my $runQ = $dbh->prepare($isoQ);
  $runQ->execute;
  return(fetchResults($runQ));
}

# INPUT: coord and mySQL table.
# OUTPUT: reference to array of rows from query.
# IF multiple tables are given, the last one must not have a type, and must be the coordinate queried table..
sub singleCoordSQLquery {
  my($coor, $return, $tables, $where) = @_;
  my($chrom, $start, $stop, undef) = parseRegion($coor);
  my($bins) = "";

  my(@binsArr) = allBinsFromRangeExtended((($stop-$start) >> 1) + $start);
  $bins = join(" or kg.bin = ", @binsArr);
  $bins = "and (kg.bin = $bins\)";

  unless(defined($return)) { $return = "kg.name, kg.chrom, kg.strand, kg.cdsStart, kg.cdsEnd, kg.exonStarts, kg.exonEnds, kg.txStart, kg.txEnd"; }
  unless(defined($tables)) { die("No table given to singleCoordSQLquery!!--- usage: singleCoordSQLquery(\$coord, \$returntypes, \$table, \$optionalWhere)\n"); }
  if(!defined($where)) { $where = ""; } else { $where = "AND $where"; } 
  
  $where = "kg.chrom = \"$chrom\" AND kg.txStart <= $start AND kg.txEnd >= $stop $bins $where";

  return(singleSQLquery(undef, $return, "$tables kg", $where));
}

# INPUT: coord and mySQL table.
# OUTPUT: reference to array of rows from query.
# IF multiple tables are given, the last one must not have a type, and must be the coordinate queried table..
sub executeCoordSQLquery {
  my($runQ, $coor, $optVarRef) = @_;
  my($chrom, $start, $stop, $strand) = parseRegion($coor);
  my($bins) = "";
  unless(defined($runQ)) { die "Must prepare query with prepareCoordSQLquery or run basicCoordSQLquery instead!\n"; }
  my(@binsArr) = allBinsFromRangeExtended((($stop-$start) >> 1) + $start);

  unless(defined($strand)) { 
    $strand = "%";
  }
  if(defined($optVarRef)) {
    $runQ->execute($chrom, $start, $stop, $strand, @binsArr, @$optVarRef);
  } else {
    $runQ->execute($chrom, $start, $stop, $strand, @binsArr);
  }
  return(fetchResults($runQ));    
}

# INPUT: coord and mySQL table.
# OUTPUT: reference to array of rows from query, matching EXACT BIN, not all overlapping bins as executeCoordSQLquery does.
# IF multiple tables are given, the last one must not have a type, and must be the coordinate queried table..
sub executeExactCoordSQLquery {
  my($runQ, $coor, $optVarRef) = @_;
  my($chrom, $start, $stop, $strand) = parseRegion($coor);
  my($bins) = "";
  unless(defined($runQ)) { die "Must prepare query with prepareCoordSQLquery or run basicCoordSQLquery instead!\n"; }
  my($bin) = binFromRangeExtended($start, $stop+5); #hack inc
  unless(defined($strand)) {
    $strand = "%";
  }
  if(defined($optVarRef)) {
    $runQ->execute($chrom, $start, $stop, $strand, $bin, @$optVarRef);
  } else {
    $runQ->execute($chrom, $start, $stop, $strand, $bin);
  }
  return(fetchResults($runQ));
}

sub fetchResults {
  my $runQ = shift;
  my(@results);
  while(my @row = $runQ->fetchrow_array()) {
    unless(scalar(@row) > 0) { return undef; }
    push(@results, \@row);
  }
  return(\@results);
}

# This is used to prepare an sql query to a short track table.. like microRNA target sites.
# The expectation is that the whole track entry is within the bounds of the coordinate.
sub prepareSmallCoordSQLquery {
  my($return, $tables, $where, $bins, $opt) = @_;
  unless(defined($bins)) { $bins = "kg.bin = ?"; }
  unless(defined($return)) { $return = "kg.name, kg.chrom, kg.strand, kg.chromStart, kg.chromEnd"; }
  unless(defined($tables)) { die("No table given to coordSQLquery!!--- usage: prepareSmallCoordSQLquery(\$coord, \$returntypes, \$table, \$optionalWhere)\n"); }
  if(!defined($where)) { $where = ""; } else { $where = "AND $where"; }
  my $isoQ = "SELECT $return FROM $tables kg WHERE kg.chrom = ? AND kg.chromStart >= ? AND kg.chromEnd <= ? AND kg.strand LIKE ? AND $bins $where $opt;";
  my $runQ = $dbh->prepare($isoQ);
  return($runQ);
}

# Use this to prepare a query against a genePred style formatted table (like ensGene or ucscGene or refGene etc)
# The coordinate should be small and within the bounds of the larger gene annotation entry.
sub prepareCoordSQLquery {
  my($return, $tables, $where, $bins) = @_;
  unless(defined($bins)) { $bins = "(kg.bin = ? OR kg.bin = ? OR kg.bin = ? OR kg.bin = ? OR kg.bin = ?)"; }
  unless(defined($return)) { $return = "kg.name, kg.chrom, kg.strand, kg.cdsStart, kg.cdsEnd, kg.exonStarts, kg.exonEnds, kg.txStart, kg.txEnd"; }
  unless(defined($tables)) { die("No table given to coordSQLquery!!--- usage: coordSQLquery(\$coord, \$returntypes, \$table, \$optionalWhere)\n"); }
  if(!defined($where)) { $where = ""; } else { $where = "AND $where"; }
  my $isoQ = "SELECT $return FROM $tables kg WHERE kg.chrom = ? AND kg.txStart <= ? AND kg.txEnd >= ? AND kg.strand LIKE ? AND $bins $where;";
  my $runQ = $dbh->prepare($isoQ);
  return($runQ);
}

# not truely random, biased towards bins with less elements..
sub randomSQLquery {
  my($table) = @_;
  unless(defined($table)) { die "no table supplied to randomGeneQueryMySQL!\n"; }
  unless(defined($MAXBIN)) { 
    my($rQ) = $dbh->prepare("SELECT max(bin) FROM $table;");
    $rQ->execute;
    ($MAXBIN) = $rQ->fetchrow_array();
  }
  my($bin) = int(rand($MAXBIN));
  my($isoQ) = "SELECT kg.name, kg.chrom, kg.strand, kg.cdsStart, kg.cdsEnd, kg.exonStarts, kg.exonEnds, kg.txStart, kg.txEnd FROM $table kg WHERE kg.bin = $bin ORDER BY RAND() LIMIT 1;";
  my($runQ) = $dbh->prepare($isoQ);  #prepare/execute query
  $runQ->execute;
#  my(@results);
#  while(my @row = $runQ->fetchrow_array()) {
#    push(@results, \@row);
#  }
  return($runQ->fetchrow_array());
}


### CODE ADAPTED FROM Kent Source
#our @binOffsetsExtended =(4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0);
our @binOffsetsExtended =(512+64+8+1, 64+8+1, 8+1, 1, 0);
our $binFirstShift = 17;
our $binNextShift = 3;


# int start, int end
sub binFromRangeExtended {
  my($start, $end) = @_;
#/* Given start,end in chromosome coordinates assign it
# * a bin.   There's a bin for each 128k segment, for each
# * 1M segment, for each 8M segment, for each 64M segment,
# * for each 512M segment, and one top level bin for 4Gb.
# *      Note, since start and end are int's, the practical limit
# *      is up to 2Gb-1, and thus, only four result bins on the second
# *      level.
# * A range goes into the smallest bin it will fit in. 
  my $startBin = $start;
  my $endBin = $end-1;
  $startBin >>= $binFirstShift;
  $endBin >>= $binFirstShift;
  for (my $i=0; $i < scalar(@binOffsetsExtended); ++$i)
    {
    if ($startBin == $endBin) {
        return $binOffsetsExtended[$i] + $startBin; }
    $startBin >>= $binNextShift;
    $endBin >>= $binNextShift;
    }
  die sprintf("start %d, end %d out of range in findBin (max is 2Gb)", $start, $end);
  return 0;
}

sub allBinsFromRangeExtended {
  my($position) = @_;

  my $startBin = $position;
  my %output;
  $startBin >>= $binFirstShift;
  for (my $i=0; $i < scalar(@binOffsetsExtended); $i++) {
    $output{$binOffsetsExtended[$i] + int($startBin)} = "";
    $startBin >>= $binNextShift;
  }
  return keys %output;
}

##########################################################################################

1;
__END__
