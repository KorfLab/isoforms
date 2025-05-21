Isoforms
========

## Manifest ##

- APCanalysis - Ismael's sandbox
- archive - old stuff no longer used
- bin - programs
- data - output of geniso runs mostly
- gong - Gong's sandbox
- ian - Ian's sandbox
- lib - libraries (just isoform.py right now)
- models - PWMs and the worm.splicemodel

## Programs ##

- `cmpiso` - compares collections of isoforms in GFF
- `geniso` - generates isoforms and their probabilities (old version)
- `modelbuilder` - creates the various model files
- `optiso` - optimizes model parameters with a genetic algorithm
- `run_geniso` - parellelizes geniso on the smallgenes dataset
- `run_optiso` - parellelizes optiso
- `nmd-ish.py` - re-scores isoforms that are NMD targets

-   `isorandom` - counting isoforms in random sequences


## Data

Data collection is described in `smallgenes`.

- `results_optiso2.csv` - model weights for each gene
- `APCisos.base.tar.gz` - APC results without optiso and without nmd
- `APCisos.optiso.tar.gz` - APC results with optiso only
- `APCisos.nmd.tar.gz` - APC results with nmd only
- `APCisos.optiso.nmd.tar.gz` - APC results with optiso and with nmd

## Models

See the `models` directory for standard models and `modelbuilder` for
how to build the models.

## Trivia

There are 19.938 billion RNASeq_splice records in WormBase. As a rough
estimate of intron frequency, divide intron counts by 20 billion.


## Generating APC isoforms

The APC algorithm generates all possible combinatons of isoforms given
an input gene sequence in fasta format, scores each isoform, and returns
a gff of predicted isoforms.

Other components of the algorithm include `optiso`. `optiso` uses a genetic
algorithm to optimize the weights for each part of the APC model (kmer, pwm,
length).

The smallgenes dataset is used as input. `git clone` this repo in the directory
containing your other github repos, such that:
```
Code/
├──smallgenes/
├──isoforms/
```
```
git clone https://github.com/KorfLab/smallgenes.git
```
Next unpack the smallgenes dataset in `isoforms/`
```
ln -s ../smallgenes/genomes/c.elegans/smallgenes.tar.gz
tar -xvf smallgenes.tar.gz
```
Create the APC model file with `modelbuilder`. This requires the OpenTURNS library
to be installed. Can be installed using conda:
```
conda install -y openturns
```
Next use `modelbuilder` to create `worm.splicemodel`, which contains the models used 
as input for `geniso3`. Individual model files are also created and can be viewed in `models/`
```
./modelbuilder smallgenes/ worm models/ > worm.splicemodel
mv worm.splicemodel models/
```
Next run `optiso` as the weights are a component of `geniso3`. Note that `run_optiso` by 
default uses 1 CPU core. This will take a while, I recommend doing this on
a remote cluster. Same for `run_geniso`. This program creates `results_optiso2.csv`.
```
chmod +x run_optiso2
./run_optiso2 smallgenes/ models/worm.splicemodel --cpu 15
```
`run_geniso` can be used with or without weights (`optiso`), and with or without NMD
detection. By default, one CPU core is used. Make sure the output file names 
describe the chosen parameters.
Without weights and without NMD:
```         
./run_geniso smallgenes/ models/worm.splicemodel --outdir APCisos.base/ --outname APC.base --cpu 15
```
With weights and without NMD:
```         
./run_geniso smallgenes/ models/worm.splicemodel --weights results_optiso2.csv --outdir APCisos.optiso/ --outname APC.optiso --cpu 15
```
Without weights and with NMD:
```         
./run_geniso smallgenes/ models/worm.splicemodel --outdir APCisos.nmd/ --outname APC.nmd --nmd --cpu 15
```
With weights and with NMD:
```         
./run_geniso smallgenes/ models/worm.splicemodel --weights results_optiso2.csv --outdir APCisos.optiso.nmd/ --outname APC.optiso.nmd --nmd --cpu 15
```
The output of `run_geniso` is a directory with gff files containing the APC isoforms for each gene in the smallgenes dataset. 
