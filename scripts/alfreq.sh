#!/bin/bash
#$ -l h_rt=1:00:00
#$ -l rmem=8g
#$ -t 1-2
#$ -j y
#$ -o alfreq.log

# Estimate allele frequencies from genotype likelihoods
# by maximum-likelihood; estpEM implements the 
# Li 2011 expectation-maximization (EM) algorithm 
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


