*Population Genomics Workshop 2018, University of Sheffield*
#Identification of regions of differentiation using Hidden Markov Models (HMM)
####Victor Soria-Carrasco
The aim of this practical is to estimate differentiation (FST) between a pair of populations and identify contiguous regions of differentiation across the genome using a HMM approach.

## Initial set up
We are going to create a working directory in a dedicated space in the HPC cluster (/data/$USER) and copy the necessary scripts and data files to run this practical.

Change to your data directory:
```bash
cd /data/$USER/
```
Create the working directory for this practical:
```bash
mkdir fst_hmm
```
Copy scripts:
```bash
cp -r /usr/local/extras/Genomics/workshops/January2018/fst_hmm/scripts ./fst_hmm
```
Copy data:
```bash
cp -r /usr/local/extras/Genomics/workshops/January2018/fst_hmm/data ./fst_hmm/
```
Copy results:
```bash
cp -r /usr/local/extras/Genomics/workshops/January2018/fst_hmm/results ./fst_hmm/
```
Change to the working directory we are going to use for this practical:
```bash
cd fst_hmm
```

## Data
Now let's have a look at the data.
You should have the following input files:
```bash
ls data

```
