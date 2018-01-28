#!/bin/bash
source /usr/local/extras/Genomics/.bashrc

cd /data/$USER
mkdir fst_hmm
cp -r /usr/local/extras/Genomics/workshops/January2018/fst_hmm/scripts ./fst_hmm/
cp -r /usr/local/extras/Genomics/workshops/January2018/fst_hmm/data ./fst_hmm/

cd fst_hmm
JOBID1=$(qsub scripts/bcf2gl.sh | awk '{print $3}' | perl -pe 's/\..*//g')
JOBID2=$(qsub -hold_jid $JOBID1 scripts/alfreq.sh | awk '{print $3}' | perl -pe 's/\..*//g')

cat > joinalfreq.sh <<'EOF'
echo -e "locus\tafHVA\tafHVC" > timemaHVAxHVC.alfreq.txt
join -j 1 timemaHVA.alfreq.txt timemaHVC.alfreq.txt | \
awk '{OFS="\t"; print $1,$3,$5}' >> timemaHVAxHVC.alfreq.txt
EOF

JOBID3=$(qsub -hold_jid $JOBID2 -j y -o joinalfreq.log joinalfreq.sh | awk '{print $3}' | perl -pe 's/\..*//g')
JOBID4=$(qsub -hold_jid $JOBID3 -b y -j y -o fst.log Rscript scripts/fst.R | awk '{print $3}')
JOBID5=$(qsub -hold_jid $JOBID4 -b y -j y -o fit_hmm.log Rscript scripts/fit_hmm.R | awk '{print $3}')

qsub -hold_jid $JOBID5 -b y -j y -o results.log bash -c 'rm joinalfreq.*; mkdir results; mv *.* results/'

