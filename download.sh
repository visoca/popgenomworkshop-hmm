#!/bin/bash
# Get everything in one go:
# clone repository scripts and download data and results

# get repository
git clone https://github.com/visoca/popgenomworkshop-hmm.git
cd popgenomworkshop-hmm

# get data
wget "https://doc-00-2g-docs.googleusercontent.com/docs/securesc/ha0ro937gcuc7l7deffksulhg5h7mbp1/36uesckkeefiq35qgofh4ub90gn6krob/1517112000000/16281733546219680372/*/1w9g9vHGDmt8LUiqeGeI0DP6aOSo7O_iK?e=download" -O fst_hmm_data.tar.bz2
tar -xvf fst_hmm_data.tar.bz2
rm fst_hmm_data.tar.bz2

# get results
wget "https://doc-0c-2g-docs.googleusercontent.com/docs/securesc/ha0ro937gcuc7l7deffksulhg5h7mbp1/9eflgeai3v4d6plmlud5loaisil3i3q4/1517112000000/16281733546219680372/*/1_FUNj1S3GwerRVPHlhDG5Emp3rAVRo1p?e=download" -O fst_hmm_results.tar.bz2
tar -xvf fst_hmm_results.tar.bz2
rm fst_hmm_results.tar.bz2
