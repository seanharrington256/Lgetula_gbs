#!/bin/bash -l


### This script generates multiple replicates of fastsimcoal models within bootstrap replicates
###     for each bootstrap rep, we want to make make 50 replicate FSC runs
# Run this from inside the directory in which each directory is a bootstrap rep - should also contain the est, tpl, and pv files


for dir in */; do cp *tpl *pv *est $dir; done # first copy the est, tpl, and pv files into each of the bootstrap rep directories 
mkdir ../Reps  # make a directory for the replicates to go
for dir in */; do mkdir ../Reps/$dir; done   # copy over everything into each rep
for dir in */; do for i in ../Reps/$dir/Run_{1..50}; do cp -r $dir "$i" ; done; done
