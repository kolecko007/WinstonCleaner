# Decross
Decross is a software tool for detecting and removing cross-contaminated 
contigs from assembled transcriptomes. The program uses BLAST to identify 
suspicious contigs and RPKM values to sort these as either correct or 
contamination. 

# Requirements

To run Decross, the following requirements must be satisfied:
* Python 2.7
* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [bbtools](https://jgi.doe.gov/data-and-tools/bbtools/) (`pileup.sh`)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

# Installation

0. Checkout repository

    `git clone git@github.com:kolecko007/decross_final.git`
    
    `cd decross_final`

1. Install pip dependencies:

    `pip2 install --user -r requirements.txt`

2. Initialize settings:

    `cp config/settings.yml.default config/settings.yml`

# Quick Start
1. Prepare the folder with input data and an empty folder for the results
1. Open `config/settings.yml` and specify input and output paths
1. `bin/prepare_data.py`
1. `bin/find_contaminations.py`
1. Inspect the results in the output folder

# Usage
## Input
The input data should be presented as a set of triads of files for each dataset.
For each dataset it is necessary to prepare:
* left reads `.fastq`
* right reads `.fastq`
* assembled transcriptome `.fasta` file

Names of the files must be in the following format:
* `NAME_1.fastq`
* `NAME_2.fastq`
* `NAME.fasta`

For example:
* `brucei_1.fastq`
* `brucei_2.fastq`
* `brucei.fasta`
* `giardia_1.fastq`
* `giardia_2.fastq`
* `giardia.fasta`

For file names only letters, digits and `_` symbols are allowed.

All the files must be placed together in one folder.

## Configuration

All the settings are declared in `config/settings.yml`.

* `decross.paths.input` &mdash; input folder with reads and contigs
* `decross.paths.output` &mdash; output folder with the results
* `decross.paths.tools.pileup_sh` &mdash; (_optional_) bbtools `pileup.sh` execution command
* `decross.paths.tools.bowtie2` &mdash; (_optional_) bowtie2 execution command
* `decross.paths.tools.bowtie2_build` &mdash; (_optional_) bowtie2-build execution command
* `decross.hits_filtering.len_ratio` &mdash; minimal `qcovhsp` for hits filtering
* `decross.hits_filtering.len_minimum` &mdash; minimal hit lenth for hits filtering
* `decross.coverage_ratio.regular` &mdash; coverage ratio for REGULAR dataset pair type 
(lower values make contamination prediction more strict, less contaminations will be found)
* `decross.coverage_ratio.close` &mdash; coverage ratio for CLOSE dataset pair type
* `decross.threads.multithreading` &mdash; enable multithreading (disabling is convenient for debugging purposes)
* `decross.threads.count` &mdash; number of threads if multithreading enabled
* `decross.tools.blast.threads` &mdash; number of threads for BLAST processing
* `decross.tools.bowtie.threads` &mdash; number of threads for bowtie2 processing
* `decross.in_memory_db` &mdash; (`true`|`false`) load coverage database to RAM in the beginning. 
Makes contamination lookup faster, but requires decent amount of memory.

The default configuration can be found in file `config/settings.yml.default`.

```
decross:
  in_memory_db: false

  paths:
    input: /path/to/folder/with/data/
    output: /path/to/output/folder

  hits_filtering:
    len_ratio: 70
    len_minimum: 100

  coverage_ratio:
    REGULAR: 1.1
    CLOSE: 0.04

  threads:
    multithreading:  true
    count:   8

  tools:
    blast:
      threads: 8
    bowtie:
      threads: 8
```

## Data preparation
The first step is to prepare the data for decross processing.

`bin/prepare_data.py`

The result will be stored in the folder, specified in `decross.paths.output` option.

After the preparation the file `types.csv` can be inspected and edited.
It contains all possible combinations of dataset pairs and their types.

The default types are:
* `CLOSE` - taxonomically close organisms
* `REGULAR` - simple pair of organisms

In `types.csv` there can also be specified any amount of custom types.
Their names must be in upper case. 

```
predator,prey,95.0,LEFT_EATS_RIGHT
prey,predator,95.0,RIGHT_EATS_LEFT
``` 

In these case coverage ratio for each custom type must be specified in `decross.coverage_ratio` section of
 `settings.yml` file:
 
```
...
  coverage_ratio:
    REGULAR: 1.1
    CLOSE: 0.04
    LEFT_EATS_RIGHT: 10
    RIGHT_EATS_LEFT: 0.1
...
```


## Contamination cleanup

`bin/find_contaminations.py`

## Output

The results will be saved in the folder, specified in `decross.paths.output` option.

For each datasets there will be the following structure of files.

* **DATASET_NAME_clean.fasta** &mdash; clean contigs
* **DATASET_NAME_deleted.fasta** &mdash; contaminated contigs
* **DATASET_NAME_suspicious_hits.csv** &mdash; all suspicious BLAST hits
* **DATASET_NAME_contamination_sources.csv** &mdash; 
sources of contaminations with a following columns: source contamination dataset name, number of sequences
* **DATASET_NAME_contaminations.csv** &mdash; list of blast hits from which contaminations were detected
* **DATASET_NAME_missing_coverage.csv** &mdash; list of contig ids without a coverage


# TODO
* Moving to python3
* Logging system
* Extended testing
* export to graph format
