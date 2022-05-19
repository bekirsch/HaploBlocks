# `HaploBlocks`

`HaploBlocks` implements an efficient approach for detecting positive selection in population
genomic datasets with hundreds of thousands of individuals or more.

Please cite [Kirsch-Gerweck, *et al.* Journal. (Year)](tba.hyper-link.tba) (doi: tba)
when using this program.

The current implementation has been written by Michel T Henrichs, Bastien Cazaux, with minor bug fixes by Benedikt Kirsch-Gerweck.

## Installation

Tested on Ubuntu 20.04 LTS:

```
git clone https://github.com/bekirsch/HaploBlocks.git
cd HaploBlocks/build/
make
```

### Install prerequisites

You will need the command-line tool `convert` by ImageMagick to produce simple plots generated after each run of `HaploBlocks`. The current version of ImageMagick can be installed via:

```
sudo apt-get install imagemagick
```

<!--- More elaborate figures can be produced after running `HaploBlocks` using R. To install R on Ubuntu:

```
sudo apt-get install r-base-core
```
We require the `latex2exp` and `stringr` packages. If you don't have these installed in your local R installation (make sure you have one on your system), you should be able to install them from within R via `install.packages(c("latex2exp", "stringr"))`. If you have root access to your machine, you can install packages without requiring any user interaction by
```
sudo R -e 'install.packages(c("latex2exp", "stringr"), repos="https://cran.r-project.org")'
```
-->

## Run HaploBlock

Before analysing a dataset you need to create a lookup-table for the recent common ancestry significance testing. Therefore run
```
<Path-to-HaploBlocks>/build/filter_lookup --N_e <haploid population size> --max_k <number of haploid samples> > <output filename>
```
Also, make sure you have your genetic map in plink format available. The map file is a headerless, four-column and tab delimited file. The four columns are: 

1. Chromosome code (as specified in your VCF file; usually 'chr1' or '1'),


2. Variant identifier (can be '.'),


3. Position in cM,


4. Base-pair coordinate.

**Note:** Step 3 of the example below converts a map in HapMap format to a map in plink format.


To perform a simple chromosome-wide scan for selection run
```
<Path-to-HaploBlocks>/build/full --vcf_path <path to vcf-file> \
       --genetic_map_path <path to genetic map> \
       --lookup_path <path to lookup-table> \
       --out_folder <path to output> \
```
**Note:** HaploBlocks does not process gzipped VCF-files.

### Example

The `example` folder contains a variant call file `example.vcf.gz`, which contains chromosome 2 of 503 individuals with european ancestry of the `1000 Genomes Project Phase 3` ([The 1000 Genomes Project Consortium (2015)](https://doi.org/10.1038/nature15393)) dataset. The file is subsetted to a 2.5 Mbp region around the LCT locus. A genetic map `CEU_recombination_map_hapmap_format_hg19_chr_2.txt` is located in this folder as well. This was obtained from https://drive.google.com/drive/folders/1Tgt_7GsDO0-o02vcYSfwqHFd3JNF6R06 released as part of [Spence and Song (2019)](https://doi.org/10.1126/sciadv.aaw9206).

1. Navigate to `<Path-to-HaploBlocks>/example/` and run
```
gunzip -k example.vcf.gz
```
to unpack the VCF.

2. Create the lookup-table, via
```
../build/filter_lookup --N_e 2e4 --max_k 1006 > example.lookup
```
where `max_k` is `2*503`.

3. Next, make sure the genetic map has the right formatting:

```
awk 'NR>1 { print "2\t.\t", $4, "\t",$2 }' CEU_recombination_map_hapmap_format_hg19_chr_2.txt > example.map
```

4. To run the scan execute
```
../build/full --vcf_path example.vcf \
       --genetic_map_path example.map \
       --lookup_path example.lookup \
       --out_folder ./ 
```

### Command line arguments
`--eff_pop_size <diploid population size>` sets the effective diploid population size (default: 1e4).

`--remove` deletes intermediate files at the end of a run.

`--overwrite` overwrites existing files with the same name.

`--skip_filters` skips significance testing schemes.


`...`
