# `HaploBlocks`

`HaploBlocks` implements an efficient approach for detecting positive selection in population
genomic datasets with hundres of thousands of individuals or more.

Please cite [Kirsch-Gerweck, *et al.* Journal. (Year)](tba.hyper-link.tba) (doi: tba)
when using this program.

The current implementation has been written by Michel T Henrichs, Bastien Cazaux, with minor bug fixes by Benedikt Kirsch-Gerweck.

## Installation

Tested on Ubuntu 20.04 LTS:

```
git clone https://github.com/bekirsch/HaploBlocks.git
cd HaploBlocks
make
```

### Install prerequisites

You will need R (for plotting). To install R on Ubuntu:

```
sudo apt-get install r-base-core
```
We require the `latex2exp` and `stringr` packages. If you don't have these installed in your local R installation (make sure you have one on your system), you should be able to install them from within R via `install.packages(c("latex2exp", "stringr"))`. If you have root access to your machine, you can install packages without requiring any user interaction by
```
sudo R -e 'install.packages(c("latex2exp", "stringr"), repos="https://cran.r-project.org")'
```

## Run HaploBlocks

Before analysing a dataset you need to create a lookup-table for the recent common ancestry significance testing. Therefore run
```
./filter_lookup --N_e <diploid population size> --max_k <number of haploid samples> > <output filename>
```

To perform a simple chromosome-wide scan for selection run
```
./full --vcf_path <path to vcf-file>
       --genetic_map_path <path to genetic map>
       --lookup_path <path to lookup-table>
       --out_folder <path to output>
```
**Note**, that HaploBlocks does not process gzipped VCF-files.

### Command line arguments
`--eff_pop_size <diploid population size>` sets the effective diploid population size (default: 1e4).

`--remove` deletes intermediate files at the end of a run.

`--overwrite` overwrites existing files with the same name.

`--skip_filters` skips significance testing schemes.


`...`
