a<em>ligate</em>r
=================

Software suite for detection of chimeric or circular RNAs from high-throughput sequencing data.
Designed for use on LIGR-seq data (see citation below), but can be applied to other RNA-RNA interaction screens.

```
Sharma, E.*, Sterne-Weiler, T.*, O’Hanlon D., Blencowe, BJ. (2016) “Global mapping of the human RNA:RNA interactome reveals new functions for non-coding RNAs.” Molecular Cell. 62(4):618-26
```

Table of Contents
-----------------

- [Requirements](#requirements)
- [Installation](#installation)
- [Overview](#overview)
- [Usage](#usage)
- [File Formats](#file-formats)

Requirements
------------

_julia v0.4_
 * ArgParse
 * Match
 * Distributions
 * GZip

_perl v5 Packages_
 * Parallel::ForkManager
 * Getopt::Long
 * DBI

_external software_
 * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) - short-read alignment
 * blastn - `ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/`
 * blast databses `ftp://ftp.ncbi.nlm.nih.gov/blast/db/`: nt, human_genomic, other_genomic 

_optional_
 * [RactIP](http://rtips.dna.bio.keio.ac.jp/ractip/) - intermolecular _in silico_ RNA folding
 * [nibFrag](http://hgdownload.soe.ucsc.edu/admin/exe/) - faToNib to make nib genome files from fasta format
 * [faToNib](http://hgdownload.soe.ucsc.edu/admin/exe/) - nibFrag to quickly retrieve sequences from nib formatted chromosomes
 

Installation
------------

- [System](#system)
- [Database Setup](#database-setup)

###System###

Install most recent [julia release here](http://julialang.org/downloads/), which must be = v0.4 (v0.5 is not yet supported).  If you are new to julia, or installing programs via command line, there is a [helpful guide here](https://en.wikibooks.org/wiki/Introducing_Julia/Getting_started).  The packages can then be installed manually by opening the `julia` REPL:
```julia
  pkgs = ("ArgParse", "Match", "Distributions", "GZip")
  map( Pkg.add, pkgs )
  map( Pkg.test, pkgs ) 
```
The `perl` packages can be installed using `cpan -i Parallel::ForkManager`

###Database Setup###

The current database is for hg19, you can download it here: [hg19 transcriptome](http://hgwdev.sdsc.edu/~timsw/GRCh37.v19.bt2.tar.gz)

It is also possible to use a different build or species, but automation of this process is still in development.  In the meantime you should be able to reverse engineer the file formats in the hg19 version. Briefly fasta transcriptome headers need to be in the format of `>ENST00001_ENSG00001_GENESYM_BIOTYPE`, for example `>ENST00000432079.1_ENSG00000116747.8_TROVE2_protein-coding`.  If you are having trouble feel free to open an issue or e-mail me.

Set up the default database like:
```bash
$ git clone https://github.com/timbitz/Aligater.git
$ cd Aligater
$ wget http://hgwdev.sdsc.edu/~timsw/GRCh37.v19.bt2.tar.gz
$ tar xzvf GRCh37.v19.bt2.tar.gz
```

Overview
--------

<div align=center><img src="/lib/overview.png" width="500px" /></div>


Usage
-----

- [Alignment and Detection](#alignment-and-detection)
- [Post Processing](#post-processing)
- [Reclassification](#reclassification)
- [Statistics](#statistics)

The aligater executable is a wrapper for a set of tools that are meant to work
as input/output for one another, often compatible with piping one into the next.
```
$ aligater -h
     aligater [sub-command] [-h]
         - align   : align short-reads to transcriptome
         - detect  : detect chimeric reads by recursive chaining of transcriptome SAM blocks
         - post    : post-process LIG format files with BLAST or RACTIP
         - reclass : create 2D k-mer db and reclassify chimeras using heirarchical type
         - stats   : compare crosslinked to mock-treated samples using multinomial statistics
         - table   : compile interaction results into tabular format
```

###Alignment and Detection###

Using the default alignment parameters is recommended for LIGR-seq and the step
can be run as simply as:
```bash
$ filename="somefile.fastq.gz"
$ nodir=`basename $filename`
$ prefastq=${nodir%%.*}

$ aligater align -x db $filename > sam/$prefastq.sam
```
or `--bam` flag will pipe the output of aligater align into `samtools view -bS`

But feel free to alter other alignment parameters as you see fit for custom purposes:
```bash
$ aligater align -h
```

Example:
```bash
$ align_max=50
$ aligater align -p $cores -k $align_max --bam -x db/GRCh37.v19 $filename > bam/$prefastq.bam
```

The next step is detection, which can be piped from the first: `aligater align | aligater detect > out.lig`
```bash
$ detectparam='--gtf [annoFile.gtf(.gz)] --gfam [gene_fam.txt(.gz)] --rmsk [maskerFile.bed(.gz)]'
$ aligater detect $detectparam < sam/$prefastq.sam > lig/$prefastq.lig
```
for example with the default transcriptome you should have something like:
```bash
$ detectparam='--gtf db/GRCh37.v19.gtf.gz --gfam db/hg.gene_fam.txt.gz --rmsk db/GRCh37.repeatMasker.slim.bed.gz'
```
will print a [lig formatted](#file-formats) file to STDOUT.

or if bam output was specified from the first command.
```bash
$ samtools view -h bam/$prefastq.bam | aligater detect $detectparam > lig/$prefastq.lig
```

###Post Processing###

There are a few parts to the post processing step, a number of filtering flags, a `--blast` (mandatory) step and a folding step using `--ractip` (optional), each requiring the use of external dependencies (installed to your `$path`), `blastn` and `ractip`, see [requirements](#requirements).

```bash
aligater post [filtering_flags] [--blast] [--ractip] 
```

It is recommended that these commands be run separately, to effectively make use system resources. 
For example, `--ractip` takes several hours with substantial CPU on multiple cores (`-p`), but doesn't require much RAM.  
Meanwhile, `--blast` is very CPU heavy and requires significant RAM as well!

The first command `aligater post --blast` should be considered mandatory, it is *not* recommended to skip this step.
This step requires proper installation of `blastn` and the blast databases to the environmental variable `BLASTDB`.

```bash
$ export BLASTDB="/path/to/blast/databases"
$ blastn -version
blastn: 2.2.28+
Package: blast 2.2.28, build Mar 12 2013 16:52:31

$ aligater post [--tmp (def: /tmp/aligater)] [-p] --loose --blast < lig/$prefastq.lig > lig/$prefastq.blast.lig
```
Edit the `tmp` directory to hit and number of threads `-p` to use as necessary.


###Reclassification###

This is a two part step, you need to first create a junction library `.jlz` from the
annotations of all samples and replicates that you plan to compare:
```bash
$ cat lig/*blast.filtered.lig | aligater reclass --uniq --save database.jlz
```

and then the actual reclassification step on each `lig` file:
```bash
$ aligater reclass --load database.jlz -g -b -u < lig/$input.lig > lig/$input.reclass.lig
```

###Statistics###

This step produces the final output of the `aligater` package and takes a number of command line options that require that the input file names be properly formatted! 

The input expects two foreground files, a crosslinking reagent treated (xlink) or mock-treated/control (unxlink) which are identical file names other than some specific identifier that you could specify as a wildcard. For example:

```bash
$ xlinkfile="sample_rep1-xlink.final.lig"
$ unxlinkfile="sample_rep1-unx.final.lig"

$ foregroundfile="sample_rep1-%.final.lig"
```
along with an `--nd xlink,unx` flag which specifies the strings to insert into the `%` wildcard to make the two input files.

Similarly, the two background files containing non-chimeric expression level `.lig` files must be similarly formatted with a `%` wildcard that will be interpolated with the same `--nd` strings.
```bash
$ backgroundfile="sample_rep1-%.expression.lig"

$ nameParam="--fore $foregroundfile --back $backgroundfile --nd xlink,unx"
```

Additionally we have to supply the comma delimited normalization constants (total fastq input read numbers) `--nc` for the foreground files:
```bash
$ totalAMTreads=73217814
$ totalMOCKreads=62648851

$ normParam="--nc $totalAMTreads,$totalMOCKreads"
```

This makes up the core arguments to `aligater stats`:
```bash
$ aligater stats $nameParam $normParam > output.pvl
```

There are a number of other arguments which greatly expand the data compiled by `aligater stats` such as the `--vs` option which allows an arbitrary number of variable columns to summarize for each interaction 'col:type,col:type,etc..'.  Each entry consists of a column number (1-based) followed by a `:` and a data type character `[pcfdn]` where each character stands for:
 * `p` : paired column with colon delimited entries for example BIOTYPEA:BIOTYPEB
 * `s` : string containing column to include
 * `f` : floating point number within the scope of Float64
 * `d` : integer number within the scope of Int64
 * `n` : some member of the Number abstract type, could be complex
 
For default `lig` format you can use: `--vs 18:f,24:p,21:f,22:d,23:d` which should summarize most of the important columns.

NOTE: If you skipped the RactIP folding step in `aligater post --ractip` then you will have a different number of columns than the default expects.  To change command line options for `aligater stats` to be compatible with your input `.lig` files, use `--gi 19` and `--vs 18:p`.

Additionally there is a `--filt` optional flag that can be used to specify a single column and regex pattern for filtering purposes.  For example you may want to filter for only intermolecular interactions using `--filt 1:I`.

###File Formats###

The `aligater detect` step outputs a basic .lig (tab delimited) format file, and the `aligater post` outputs an extended version of the .lig fileformat.  The standard format is:

| Column | Example | Regex | Description | 
| ------ | ------- | ----- | ----------- |
| 1 | I | [IPS] | Single letter code classification: I=Intermolecular, P=Putatively Paralogous, S=Intramolecular |
| 2 | A:B | [A-Z\:]+ | Chimeric Read Structure: A,B,C refer to molecules, : denotes a ligation.  A:A = intramolecular ligation, A:B = intermolecular ligation, A:B:A = intermolecular ligation from A to B and then from B back to A |
| 3 | 1:3 | [\d\:]+ | Local alignments (from SAM) that have been chained to make up the chimera |
| 4 | RNU4-1:RNU6-1251P | \S+ | Gene symbols of transcripts aligned to chimeric read |
| 5 | ENSG00000200795.1:ENSG00000201372.1 | \S+ | Gene ids of chimeric segments |
| 6 | ENST00000363925.1:ENST00000364502.1 | \S+ | Transcript ids of chimeric segments |
| 7 | snRNA:snRNA | \S+ | Biotypes of transcripts in chimeric segments |
| 8 | U4_snRNA:U6_snRNA | \S+ or NA | If alignment overlaps RepeatMasker annotation, RepeatNames or NA are given in A:B order |
| 9 | snRNA:snRNA | \S+ or NA | If alignment overlaps RepeatMasker annotation, then RepeatClass or NA is given in A:B order |
| 10 | 68c3777d9e14bf33 | \S+ or [a-f1-9]+ | Readname or md5 hash of readname |
| 11 | ACTGGCAATTTAAAATTGGAA+ | [ATGCN]+ | Read Sequence, _ indicates ligation site |
| 12 | 124 | \d+ | Alignment Score |
| 13 | 1,30>30>17 | [\d\,\>\(\)]+ | Chimera uniqueness, deprecated |
| 14 | 74,26 | [\d,]+ | Alignment positions in transcripts |
| 15 | 74,27 | [\d,]+ or NA,NA | RepeatMask offset positions |
| 16 | 55,41 | [\d,]+ | Alignment lengths |
| 17 | chr12:120730966:-,chr20:42101679:+ | POS,POS where POS=\S+\:\d+\:[+-] | Genome start position of alignment segment |

The `aligater stats` step outputs a `.pvl` file consistent with `Extended Data Table 1` of Sharma E, Sterne-Weiler T, et al. 2016:

```
Extended Data Table 1
Gene-ids : Comma delimited HUGO or Repeat family name in lexographical order
OE[+amt/-amt] : (+AMT/-AMT) / (Expected(+AMT)/Expected(-AMT))
AMT/Mock : (+AMT/-AMT)
Exp[AMT/Mock] : Expected(+AMT)/Expected(-AMT)
AMTReads : Number of reads in the +AMT sample + 1 pseudo count
OE-amt : Observed / Expected for +AMT sample (see Methods)
MockReads : Number of reads in the -AMT sample + 1 pseudo count
OE-mock : Observed / Expected for -AMT sample (see Methods)
AMT+ pval : Binomial p-value for significance in the +AMT sample
Mock pval : Binomial p-value for significance in the -AMT sample
AMT RPM : The transcript's expression in Reads per Million from -ligase sample
Mock RPM : The transcript's expression in Reads per Million from -ligase sample
Transcript Biotypes : Gene biotypes from GENCODE annotations
... any other summary line from --vs string
```
