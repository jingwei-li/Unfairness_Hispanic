#!/bin/bash

proj_dir='/data/project/FairAI_Hispanic'
DIR="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

CPUS='1'
RAM='5G'
LOGS_DIR=$proj_dir/scripts/logs/ABCD_read_motion
# create the logs dir if it doesn't exist
[ ! -d "${LOGS_DIR}" ] && mkdir -p "${LOGS_DIR}"

# print the .submit header
printf "# The environment
universe       = vanilla
getenv         = True
request_cpus   = ${CPUS}
request_memory = ${RAM}
# Execution
initial_dir    = $proj_dir/scripts/Unfairness_Hispanic/preparation
executable     = /usr/bin/matlab95
transfer_executable   = False
\n"

# loop through lists
for start in $(seq 1 500 8501); do
    end=$(( start + 499 ))
    curr_ls=$proj_dir/scripts/lists/censor_${start}_${end}.txt
    printf "arguments = -singleCompThread -r ABCD_read_motion('/data/project/FairAI_Hispanic/data/inm7-superds/original/abcd/derivatives/abcd-hcp-pipeline','$curr_ls',1,'/data/project/FairAI_Hispanic/scripts/lists','_${start}_${end}')\n"
    printf "log       = ${LOGS_DIR}/\$(Cluster).\$(Process).${start}.log\n"
    printf "output    = ${LOGS_DIR}/\$(Cluster).\$(Process).${start}.out\n"
    printf "error     = ${LOGS_DIR}/\$(Cluster).\$(Process).${start}.err\n"
    printf "Queue\n\n"
done