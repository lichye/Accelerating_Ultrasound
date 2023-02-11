#!/bin/bash
# Set up input Size

INPUT_SIZE=64

let "USE_LOCAL_MIC_NUMBER=$RANDOM % 7"
echo "Running interactive job on mic$USE_LOCAL_MIC_NUMBER"
micnativeloadex beamform.mic -t 300 -d $USE_LOCAL_MIC_NUMBER -a "$INPUT_SIZE" 
micnativeloadex solution_check.mic -t 30 -d $USE_LOCAL_MIC_NUMBER -a "$INPUT_SIZE" 


 

