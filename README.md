*Population Genomics Workshop 2018, University of Sheffield*

# Identification of regions of differentiation using Hidden Markov Models (HMM)
#### Victor Soria-Carrasco
The aim of this practical is to estimate differentiation (FST) between a pair of populations and identify contiguous regions of differentiation across the genome using a HMM approach.

## Initial set up
We are going to create a working directory in a dedicated space in the HPC cluster (/data/$USER) and copy the necessary scripts and data files to run this practical.

Connect to Iceberg HPC cluster (change user by your username, e.g. cs4ab33):
```bash
ssh cs4ab33@iceberg.sheffield.ac.uk
```
Request an interactive session:
```bash
qrsh
```
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
ls -lh data
```

>total 2.5G<br>
>-rw-r--r-- 1 cs4ab33 cs  22K Jan 26 12:31 lg_ord_sca_length.dsv<br>
>-rw-r--r-- 1 cs4ab33 cs 768M Jan 26 12:31 timemaHVA.gl<br>
>-rw-r--r-- 1 cs4ab33 cs 467M Jan 26 12:31 timemaHVA.vcf.gz<br>
>-rw-r--r-- 1 cs4ab33 cs 781M Jan 26 12:31 timemaHVC.gl<br>
>-rw-r--r-- 1 cs4ab33 cs 476M Jan 26 12:31 timemaHVC.vcf.gz<br>

There is a vcf file for each population (HVA and HVC). You can have a look at them:
```bash
gzip -dc timemaHVA.vcf.gz | less -S
# or
samtools view timemaHVA.vcf.gz | less -S
``

They need to be converted to the stripped down format (“gl”, genotype likelihood) required for the program we are going to use for allele frequency estimation. This is done with a script called *bcf2gl.sh*. Let's have a look with the text editor nano:
```bash
nano scripts/bcf2gl.sh
``
