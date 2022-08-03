#!/bin/bash

## Script written by Sean Harrington - this script just loops around the fsc-selectbestrun.sh (https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh)
##       (Note that I edited that script to change run* to Run_* because I set this up differently)
##    to get the best likelihoods across replicates for each model - assumes you're in a Reps directory as created by Prep_FSC_reps.sh

source ~/.bashrc

# Loop to run the slightly modified fsc-selectbestrun.sh on each model
for x in */; do
  cd $x
  fsc-selectbestrun.sh
  cd ..
done

# make a directory to put the best bestlhoods file from each model into across all models - but only if it doesn't exist
if [[ ! -e best_L_allMods ]]; then
    mkdir best_L_allMods
elif [[ ! -d best_L_allMods ]]; then
    echo "best_L_allMods already exists but is not a directory" 1>&2
else
	echo "best_L_allMods already exists" 1>&2
fi

# copy the bestlhoods files into there
cp */bestrun/*bestlhoods best_L_allMods


