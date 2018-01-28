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

