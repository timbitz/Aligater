a<em>ligate</em>r
=================

Software suite for detection of chimeric or circular RNAs from high-throughput sequencing data.
Designed for use on LIGR-seq data (citation coming..), but can be applied to regular RNAseq for
detection of circular RNAs as well (experimental).

Table of Contents
-----------------

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [File Formats](#file-formats)

Requirements
------------

_julia 0.4_
 * ArgParse
 * Match
 * Distributions
 * GZip

_perl v5 Packages_
 * Parallel::ForkManager
 * Getopt::Long

_mandatory external software_
 * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) - short-read alignment
 * blastn - `ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/`
 * blast databses `ftp://ftp.ncbi.nlm.nih.gov/blast/db/`: nr, human_genomic, other_genomic 

_optional external software_
 * [RactIP](http://rtips.dna.bio.keio.ac.jp/ractip/) - intermolecular _in silico_ RNA folding
 * [nibFrag](http://hgdownload.soe.ucsc.edu/admin/exe/) - faToNib to make nib genome files from fasta format
 * [faToNib](http://hgdownload.soe.ucsc.edu/admin/exe/) - nibFrag to quickly retrieve sequences from nib formatted chromosomes
 

Installation
------------

- [System](#system)
- [Database Setup](#database-setup)

###System###

If you are new to `julia`, the latest 0.4-dev binaries can be downloaded from the nightly builds [64-bit](https://status.julialang.org/download/linux-x86_64) or [32-bit](https://status.julialang.org/download/linux-i686).  The packages can then be installed manually by opening the `julia` REPL:
```julia
  pkgs = ("ArgParse", "Match", "Distributions", "GZip")
  map( x->Pkg.add(x), pkgs ) 
```
The `perl` packages can be installed using `cpan -i Parallel::ForkManager`

###Database Setup###

The current database is for hg19, you can download it here: [hg19 transcriptome](http://google.com)

It is also possible to use a different build or species, but automation of this process is still in development, so you should e-mail me for a synopsis of the required steps.  


Usage
-----

- [Alignment and Detection](#alignment-and-detection)
- [Post Processing](#post-processing)
- [Reclassification](#reclassification)
- [Statistics](#statistics)

The aligater executable is a wrapper for a set of tools that are meant to work
as input/output for one another, often compatible with piping one into the next.
```bash
$ aligater -h

```

###Alignment and Detection###

Using the default alignment parameters is recommended for LIGR-seq and the step
can be run as simply as:
```bash
$ filename="somefile.fastq.gz"
$ nodir=`basename $filename`
$ prefastq=${nodir%%.*}

$ aligater align -x db [file.fq(.gz)] > sam/$prefastq.sam
or
$ aligater align --bam -x db [file.fq(.gz)] > bam/$prefastq.bam
```

But feel free to alter tham as you see fit for custom purposes:
```bash
$ aligater align -h
```

The next step is detection, which can be piped from the first: `aligater align | aligater detect > out.lig`
```bash
$ detectparam='--gtf [annoFile.gtf(.gz)] --gfam [gene_fam.txt(.gz)] --rmsk [maskerFile.bed(.gz)]'
$ aligater detect $detectparam < alignFile.sam
```
will print a [lig formatted](#file-formats) file to STDOUT.

or if bam output was specified from the first command.
```bash
$ samtools view -h bam/$prefastq.bam | aligater detect $detectparam > lig/$prefastq.lig
```

###Post Processing###

There are a few parts to the post processing step, `--blast` (MANDATORY) and `--ractip` (OPTIONAL), each requiring the use of external dependencies (installed to your `$path`), `blastn` and `ractip`, see [requirements](#requirements).

```bash
aligater post [--blast] [--ractip]
```

It is recommended that these commands be run separately, to effectively make use system resources. For example, `--ractip` takes several hours with multiple cores (`-p`), but doesn't require much RAM.  Meanwhile, `--blast` does require significant RAM as well!

The first command `aligater post --blast` should be considered mandatory, it is *not* recommended to skip this step.
```bash
$ export BLASTDB="/path/to/blastdatabases/"

$ aligater post --loose --blast < lig/$prefastq.lig > lig/$prefastq.blast.lig
```

###Reclassification###

This is a two part step, you need to first create a junction library `.jlz` from the
annotations of all samples and replicates that you plan to compare:
```bash
$ aligater reclass --uniq --save database.jlz < lig/*.lig
```

and then the actual reclassification step on each `lig` file:
```bash
$ aligater reclass --uniq --load database.jlz < lig/$input.lig > lig/$input.reclass.lig
```

###Statistics###

This step produces the final output of the `aligater` package

