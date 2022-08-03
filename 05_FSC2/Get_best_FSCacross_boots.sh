#!/bin/bash

## Script written by Sean Harrington - this script just loops around the fsc-selectbestrun.sh (https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh)
##       (Note that I edited that script to change run* to Run_* because I set this up differently)
##    to get the best likelihoods across replicates for each bootstrap replicate - assumes you're in a Reps directory as created by Prep_FSC_reps_of_bootreps.sh

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

# copy the bestlhoods files into there - for the bootstaps, all files are named the same, need to rename them - use a loop
ind=1 # set up an indexing number to add to the file names
for file in */bestrun/*bestlhoods; do # for each file
	base=$(basename $file)  # get the basename of the file
	cp $file best_L_allMods/$base$ind # copy the file over into the correct directory, adding a number onto the end of the file
	ind=$((ind+1)) # increment the index
done











