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

Installation
------------

If you are new to `julia`, the latest 0.4-dev binaries can be downloaded from the nightly builds [64-bit](https://status.julialang.org/download/linux-x86_64) or [32-bit](https://status.julialang.org/download/linux-i686).  The packages can then be installed manually by opening the `julia` REPL:
```julia
  pkgs = ("ArgParse", "Match", "Distributions", "GZip")
  map( x->Pkg.add(x), pkgs ) 
```
The `perl` packages can be installed using `cpan -i Parallel::ForkManager`

Usage
-----

- [Database Setup](#database-setup)
- [Alignment and Detection](#alignment-and-detection)
- [Post Processing](#post-processing)
- [Statistics](#statistics)
- 
