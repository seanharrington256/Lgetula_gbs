#!/bin/bash -l

### This script generates multiple replicates of fastsimcoal models for replicate runs for estimating parameters of demographic models - see description below for more detail
# takes a command line argument specifying if you are using either multi-dimensional SFS (argument MSFS) or two-dimensional SFS (argument jointMAF)
# usage examples:
	# Prep_FSC_reps.sh MSFS
	# Prep_FSC_reps.sh jointMAF


# Starting with a directory named "Models" that contains the SFS file and multiple directories, one for each model to be tested, 
#  each containing .tpl and .est files for FSC. Run this script while the working directory is Models and it
#  copies these over into a directory called Reps_input in the parent directory of Models. Within Reps, it
#  makes a directory for each model, and within that 50 folders for replicate FSC runs - final line 
#  does some text grabbing to change the names of the .obs SFS file in each folder to match the .est file.

# note that the if/else statement causes this to run ONLY if the directory is named exactly Models - if you want to use other names & run it, remove the if/else, just be careful where you are



if [ -z "$1" ] # first check that there is an argument specifying multi- or two-dimensional SFS and exit if Not
then
	echo "an argument 'MSFS' or 'jointMAF' must be supplied to indicate multidimensional or joint SFS: Prep_FSC_reps.sh MSFS or Prep_FSC_reps.sh jointMAF"
    exit
fi

if [ ${PWD##*/} == "Models" ] # only do the operation if the current directory is Models - don't want to generate weird nonsense we're not supposed to be
then
	mkdir ../Reps     # make a Reps directory in the parent directory of Models
	for dir in */; do mkdir ../Reps/$dir; done   # 
	for dir in */; do for i in ../Reps/$dir/Run_{1..50}; do cp -r $dir "$i" ; done; done
	if [ $1 == "MSFS" ]
	then
		for dir in */; do for i in ../Reps/$dir/Run_{1..50}; do s=$(echo $i/*.est); s=${s##*/}; s=${s%.est}; cp *.obs $i/$s"_MSFS.obs" ; done; done
	fi
	if [ $1 == "jointMAF" ]
	then
		for dir in */; do for i in ../Reps/$dir/Run_{1..50}; do s=$(echo $i/*.est); s=${s##*/}; s=${s%.est}; cp *.obs $i/$s"_jointMAFpop1_0.obs" ; done; done
	fi
else
	echo "Not in a Models directory - make sure you are in the right directory or rename the directory"
fi

