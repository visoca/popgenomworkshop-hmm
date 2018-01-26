*Population Genomics Workshop 2018, University of Sheffield*

# Delimitation of contiguous regions of differentiation using Hidden Markov Models (HMM)

#### Victor Soria-Carrasco

The aim of this practical is to estimate differentiation (i.e. F<sub>ST</sub>) between a pair of populations and identify contiguous regions of differentiation across the genome using a HMM approach.

## 1. Initial set up
We are going to create a working directory in a dedicated space in the HPC cluster (/data/$USER) and copy the necessary scripts and data files to run this practical.

Connect to Iceberg HPC cluster (change user by your username, e.g. cs4ab33):
```bash
ssh cs4ab33@iceberg.sheffield.ac.uk
```
Request an interactive session:
```bash
qrsh
```
#### Important note
***
This tutorial relies on having access to a number of programs. The easiest way is to have your account configured to use the Genomics Software Repository. If that is the case you should see the following message when you get an interactive session with ```qrsh```:
```
  Your account is set up to use the Genomics Software Repository
     More info: http://soria-carrasco.staff.shef.ac.uk/softrepo
```
If you don't get that message, follow the instructions [here](http://soria-carrasco.staff.shef.ac.uk/softrepo/) to set up your account.

In addition, if you want to configure the ```nano``` text editor to have syntax highlighting and line numbering, you can configure it this way:
```bash
cat /usr/local/extras/Genomics/workshops/January2018/.nanorc >> /home/$USER/.nanorc
```
***

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

## 2. Data formatting
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

There is a vcf file containing single nucleotide polymorphisms (SNPs) from whole genome sequences of 20 individuals for each population (HVA and HVC). vcf is a very popular format for genetic variants, you can find more info [here](http://www.internationalgenome.org/wiki/Analysis/vcf4.0/). You can have a look at the file content with the following commands:
```bash
gzip -dc timemaHVA.vcf.gz | less -S
# or with bcftools
bcftools view timemaHVA.vcf.gz | less -S
# excluding long header
bcftools view -H timemaHVA.vcf.gz | less -S
```
They need to be converted to the stripped down format (*gl*, genotype likelihood) required for ```estpEM```, the program we are going to use for allele frequency estimation later. This is done with a custom Perl script called ```bcf2gl.pl```. To process the files in the cluster in parallel in different nodes, we will be using an [SGE array job](http://docs.hpc.shef.ac.uk/en/latest/parallel/JobArray.html) that calls the Perl script for each file. Let's have a look at the ```bcf2gl.sh``` script with the text editor nano:
```bash
nano scripts/bcf2gl.sh
```
```bash
#!/bin/bash
#$ -l h_rt=1:00:00
#$ -l rmem=8g
#$ -t 1-2
#$ -j y
#$ -o bcf2gl.log

# Convert bcf/vcf to gl format
# (required for estpEM)

BCF2GL="/data/$USER/fst_hmm/scripts/bcf2gl.pl"

INDIR="/data/$USER/fst_hmm"

# array of input files
INPUT=('data/timemaHVA.vcf.gz' 'data/timemaHVC.vcf.gz')
# array of output files
OUTPUT=('timemaHVA.gl' 'timemaHVC.gl')
# array of log files
LOG=('timemaHVA.bcf2gl.log' 'timemaHVC.bcf2gl.log')

# Take index from SGE task id
I=$(($SGE_TASK_ID-1))

hostname > ${LOG[$I]}
uname -a >> ${LOG[$I]}
date >> ${LOG[$I]}
echo "----------------------------------------------------------" >> ${LOG[$I]}
echo >> ${LOG[$I]}

$BCF2GL -i ${INPUT[$I]} -o ${OUTPUT[$I]} >> ${LOG[$I]} 2>&1

echo  >> ${LOG[$I]}
echo "----------------------------------------------------------" >> ${LOG[$I]}
date  >> ${LOG[$I]}
```
Submit the job using the command:
```bash
qsub scripts/bcf2gl.sh
```
You can check the status of your jobs with ```Qstat```. Jobs would take ~10 min per file to run, but if you don't want to wait, you can use the files provided in the data folder after deleting the job (replace #jobid with the actual job id):
```bash
qdel #jobid
cp -r data/*.gl ./
```
Output .gl files should look like this:
```bash
less -S timemaHVA.gl
```
>20 4391556<br>
>timemaHVA_8021022 timemaHVA_8021024 timemaHVA_8021026 timemaHVA_8021027 timemaHVA_8021028 timemaHVA_<br>
>1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1<br>
>lg01_ord0000_scaf00353:8612 0 0 0 0 0 0 0 6 75 0 6 7 0 15 85 0 6 41 0 3 34 0 3 4 0 0 0 0 9 46 0 9 6<br>
>lg01_ord0000_scaf00353:8752 0 9 6 0 12 10 0 3 37 0 6 7 29 0 6 0 24 71 0 7 14 0 18 12 0 12 70 0 9 75<br>
>lg01_ord0000_scaf00353:8758 0 12 11 0 18 19 0 3 40 0 6 7 0 6 45 0 24 77 0 39 60 0 21 16 0 9 42 0 12<br>
>lg01_ord0000_scaf00353:8770 0 12 11 0 18 19 0 6 69 0 6 7 0 6 44 0 24 72 0 42 62 0 27 20 0 6 7 0 9 3<br>
>lg01_ord0000_scaf00353:8772 0 12 11 0 18 19 0 6 67 0 6 7 0 6 39 0 24 73 0 13 17 0 27 20 0 6 7 0 9 3<br>
>lg01_ord0000_scaf00353:8797 0 15 14 0 18 43 0 6 72 0 6 7 0 6 45 0 24 48 0 40 26 0 30 22 0 9 40 0 6<br>
>lg01_ord0000_scaf00353:8838 0 12 11 0 12 37 0 3 40 0 6 7 0 0 0 0 18 12 0 18 5 0 24 29 0 9 43 0 9 46<br>
>lg01_ord0000_scaf00353:8845 0 12 11 22 0 0 40 3 0 0 6 7 0 3 37 0 15 12 0 54 82 0 21 23 29 0 0 34 1<br>
>lg01_ord0000_scaf00353:8851 0 12 11 0 12 43 0 3 35 0 9 9 0 3 39 0 15 12 0 51 84 0 18 16 0 9 41 0 6<br>
>lg01_ord0000_scaf00353:8862 0 12 11 0 6 39 0 3 34 0 9 9 0 3 33 0 15 12 0 48 79 0 21 16 0 9 43 0 9 4<br>

The first line shows the number of samples and the number of SNPs, the second line show the samples ids, and the third line specifices what samples will be used for analyses. The lines below represent a SNP per line. The first field is the id (scaffold:position in this case) and it is followed by the phred scaled genotype likelihoods for the three possible genotypes (all variants are biallelic) for each individual. Therefore there are 3 x number of samples fields genotype likelihoods on each line.

## 3. Allele frequency estimation
We will infer allele frequencies from genotype likelihoods by maximum likelihood, using an implementation of the iterative soft expectation-maximization algorithm (EM) described in [Li 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/) and also used in bcftools. This algorithm has been implemented in the program ```estpEM```, withcode kindly provided by [Zach Gompert, Utah State University](https://gompertlab.wordpress.com/)). As before, we are going to use and SGE array job to parallelize jobs. In this case it will run ```estpEM``` to estimate allele frequencies from genotype likelihoods. You can have a look at the script:

```bash
nano scripts/alfreq.sh
```
```bash
#!/bin/bash
#$ -l h_rt=1:00:00
#$ -l rmem=8g
#$ -l mem=16g
#$ -t 1-2
#$ -j y
#$ -o alfreq.log

# Estimate allele frequencies from genotype likelihoods
# by maximum-likelihood; estpEM implements the 
# Li 2012 expectation-maximization (EM) algorithm 
# (also used in bcftools)

ESTPEM='estpEM'

INDIR="/data/$USER/fst_hmm"

# array of input files
INPUT=('timemaHVA.gl' 'timemaHVC.gl')
# array of output files
OUTPUT=('timemaHVA.alfreq.txt' 'timemaHVC.alfreq.txt')
# array of log files
LOG=('timemaHVA.alfreq.log' 'timemaHVC.alfreq.log')
# Take index from SGE task id
I=$(($SGE_TASK_ID-1))

hostname > ${LOG[$I]}
uname -a >> ${LOG[$I]}
date >> ${LOG[$I]}
echo "----------------------------------------------------------" >> ${LOG[$I]}
echo >> ${LOG[$I]}

$ESTPEM -i ${INPUT[$I]} -o ${OUTPUT[$I]} -h 2  >> ${LOG[$I]} 2>&1

echo  >> ${LOG[$I]}
echo "----------------------------------------------------------" >> ${LOG[$I]}
date  >> ${LOG[$I]}
```
Submit the job using the command:
```bash
qsub scripts/alfreq.sh
```
You can check the status of your jobs with ```Qstat```. Each job task should take less than 5 min per file to run. When finished, we can have a look at the log files:
```bash
less -S timemaHVA.alfreq.log
```
It should look like this:
>node211<br>
>Linux node211 2.6.32-696.18.7.el6.x86_64 #1 SMP Wed Jan 3 19:31:16 CST 2018 x86_64 GNU/Linux<br>
>Fri Jan 26 18:36:16 GMT 2018<br>
>----------------------------------------------------------<br>
><br>
>Reading data from timemaHVA.gl<br>
>Number of loci: 4391556<br>
>Number of individuals: 20<br>
>Using EM algorithm to estimate allele frequencies<br>
>Writing results to timemaHVA.alfreq.txt<br>
>Runtime: 0 hr 1 min 11 sec<br>
><br>
>----------------------------------------------------------<br>
>Fri Jan 26 18:37:28 GMT 2018<br>

And to the output text files:
```bash
less -S timemaHVA.alfreq.txt
```
That should look like this:
>lg01_ord0000_scaf00353:8612 0.2074 0.0006<br>
>lg01_ord0000_scaf00353:8752 0.1963 0.0674<br>
>lg01_ord0000_scaf00353:8758 0.1250 0.0001<br>
>lg01_ord0000_scaf00353:8770 0.1199 0.0001<br>
>lg01_ord0000_scaf00353:8772 0.1791 0.0408<br>

The columns are: locus, allele frequency from counts, and allele frequency inferred by maximum likelihood. Now where going to merge the results from both populations into a single file:

```bash
# create header
echo -e "locus\tafHVA\tafHVC" > timemaHVAxHVC.alfreq.txt

# join files by locus (1st field), retain only locus id (1st field) and ML AFs (3rd and 5th fields)
$ join -j 1 timemaHVA.alfreq.txt timemaHVC.alfreq.txt | \
awk '{OFS="\t"; print $1,$3,$5}' >> timemaHVAxHVC.alfreq.txt

# show file
less -S timemaHVAxHVC.alfreq.txt
```
>locus afHVA afHVC<br>
>lg01_ord0000_scaf00353:8612 0.0006 0.0005<br>
>lg01_ord0000_scaf00353:8752 0.0674 0.1185<br>
>lg01_ord0000_scaf00353:8758 0.0001 0.0201<br>
>lg01_ord0000_scaf00353:8770 0.0001 0.0000<br>
>lg01_ord0000_scaf00353:8772 0.0408 0.0263<br>

## 4. F<sub>ST</sub> estimation
We are going to calculate genetic differentiation between two populations and for a large
number of variants using F<sub>ST</sub>. We are going to use the F<sub>ST</sub> Hudson's estimator for every SNP:

where *Hw* is the within-population heterozygosity, *Hb* is the between-population heterozygosity, and *p<sub>1</sub>* and *p<sub>2</sub>* are the allele frequencies in each population.

Notice we are going to use  ```R``` for calculating F<sub>ST</sub> and the rest of downstream analyses. Open an R session simply typing the command ```R```. 

For F<sub>ST</sub> estimating we are going to use the code provided in the script ```fst.R```. We will go here step by step (remember this has be run within ```R```).

First we are going to change the working directory:
```R
setwd("/data/$USER/fst_hmm")
```
and then we are going to load a library that will greatly speed up the loading of bit files:
```R
# library to speed up loading of big tables
library(data.table)
```
This is the function to calculate Hudson's F<sub>ST</sub>:
```R
# function to calculate Hudson's Fst
# --------------------------------------------------
hudsonFst<-function(locus=NA, p1=NA, p2=NA){
    numerator<-p1 * (1 - p1) + p2 * (1 - p2) 
    denominator<-p1 * (1 - p2) + p2 * (1 - p1) 
    fst<-1 - numerator/denominator
    out<-data.frame(locus,numerator,denominator,fst)
    return(out)
}
# --------------------------------------------------
```
Then we load allele frequencies:
```R
# load file with allele frequencies
pops.af<-fread("timemaHVAxHVC.alfreq.txt", header=T, sep="\t")
```
and calculate the F<sub>ST</sub>s:
```R
# Calculate Fst
fst.HVAxHVC<-hudsonFst(locus=pops.af$locus, p1=pops.af$afHVA,p2=pops.af$afHVC)

# Write table to file
write.table(fst.HVAxHVC, file="timemaHVAxHVC.fst.dsv", 
                quote=F, row.names=F, sep="\t")

# plot FST histogram
hist(fst.HVAxHVC$fst, breaks=50)

# Genome-wide FST statistics
# ------------------------------------------------------------------------------
# It is better to calculate mean FST using numerator and denominator
nloci<-length(na.exclude(fst.HVAxHVC$fst))
fst.mean<-1-((sum(fst.HVAxHVC$numerator)/nloci)/(sum(fst.HVAxHVC$denominator)/nloci))
fst.mean
mean(fst.HVAxHVC$fst,na.rm=T)

# min, max, quantiles
min(fst.HVAxHVC$fst,na.rm=T)
max(fst.HVAxHVC$fst,na.rm=T)
fst.quantile<-quantile(fst.HVAxHVC$fst, prob=c(0.5,0.025,0.975), na.rm=T)
fst.quantile

fst.genome<-c("genome",nloci, fst.mean,fst.quantile)
# ------------------------------------------------------------------------------

# FST stats by linkage group
# ------------------------------------------------------------------------------
lgs<-unique(gsub("_ord.*","", fst.HVAxHVC$locus))
fst.stats.lgs<-lapply(lgs, function(lg){
              fst<-fst.HVAxHVC[grep(paste("^",lg,"_",sep=""), fst.HVAxHVC$locus),]
              nloci<-length(na.exclude(fst$fst))
              fst.mean<-1-((sum(fst$numerator)/nloci)/(sum(fst$denominator)/nloci))
              fst.quantile<-quantile(fst$fst, prob=c(0.5,0.025,0.975), names=F, na.rm=T)
              return(c(nloci=nloci,
                       mean=fst.mean,median=fst.quantile[1],
                       lo95=fst.quantile[2],hi95=fst.quantile[3]))
            })
fst.stats.lgs<-data.frame(lg=lgs,do.call("rbind",fst.stats.lgs),stringsAsFactors=F)

# Add genome
fst.stats<-rbind(fst.genome,fst.stats.lgs)
# ------------------------------------------------------------------------------

# Write summary table to file
# ------------------------------------------------------------------------------
write.table(fst.stats,file="timemaHVAxHVC.lgs.fst.dsv",
            quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------
```

## 4. Delimitation of contiguous regions of differentiation using a HMM model
We will use a 3-state discrete homogeneous Hidden Markov Model (HMM) to delimit contiguous regions of genetic differentiation. We will classify the genom in regions of low, medium, and high differentiation (LDI, IDI, and HDI, respectively). This will allow to investigate the number, size, and distribution of regions of differentiation.
