Isoforms
========

## Programs ##

### Python

Python-based applications using `isoform.py` for common functions.
'isoform.py' uses itertools for APC (old version).

+ `cmpiso` - compares collections of isoforms in GFF
+ `geniso` - generates isoforms and their probabilities (old version)
+ `modelbuilder` - creates the various model files
+ `optiso` - optimizes model parameters with a genetic algorithm (old version)
+ `run_apc` - builds Makefile to run `isoformer` and `optiso` on apc set

Python-based applications using 'isoform2.py' for common functions.
'isoform2.py' uses a backtracking algorithm for APC (current version).

+ 'geniso2' - generates isoforms and their probabilites
+ 'geniso3' - generates isoforms and their probabilities, with nmd adjustment
+ 'run_geniso' - parellelizes geniso3 for use on the smallgenes dataset
+ 'optiso2' - optimizes model paremeters with a genetic algorithm (current version)
+ 'run_optiso2' - parellelizes optiso2
+ 'nmd-ish.py' - re-scores isoforms that are NMD targets

### C

The genomikon repo contains a couple of faster implementations in the
`isoformer` directory.

+ `isoformer` - this is the same as `geniso` but ~100x faster
+ `isocounter` - as above, but only counting, not calculating probabilities
+ `isorandom` - counting isoforms in random sequences

## Utilities ##

+ `conformity.py` - compares outputs of `geniso` and `isoformer`
+ `optiso-mp` - multi-processing version with some odd bugs
+ `speedo.py` - compares speeds of `geniso` and `isoformer`
+ `summary.py` - creates TSV of the apc set

## Data ##

Data collection is described in `datacore2024/project_splicing`. The 1045 genes
of the smallgenes dataset.

## Models ##

See the `models` directory for standard models and `modelbuilder` for how to
build the models.

## Trivia ##

There are 19.938 billion RNASeq_splice records in WormBase. As a rough estimate
of intron frequency, divide intron counts by 20 billion.

## Files ##
+ 'smallgenes.tar.gz' - contains the smallgenes dataset
+ 'models/worm.splicemodel' - APC models for use with geniso2/3
+ 'results_optiso2.csv' - model weights for each gene
