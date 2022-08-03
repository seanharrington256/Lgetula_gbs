## Complex cycles of divergence and migration shape population structure in the common kingsnake species complex

### Harrington and Burbrink

## IN REVIEW

This repo contains contains code used for analyses in the manuscript "Complex cycles of divergence and migration shape population structure in the common kingsnake species complex"

Raw sequence data are available from the NCBI Sequence Read Archive. 

# OR THEY WILL BE, ONCE I UPLOAD THEM

Note that these data are demultiplexed data for each sample. Our full set of raw data included samples from other snake species for different projects in addition to the samples used here. We have included the scripts and parameter files that we used to demultiplex the data for the sake of completeness and transparency, but note that you will not be able to replicate these initial steps.

# Alignments will be added to Dryad


## ipyrad data processing

Demultiplexing scripts and files are in the `01_ipyrad_step_1_Demux` directory. `.pbs` scripts were used to run demultiplexing and merging of plates on the American Museum of Natural History (AMNH) computing cluster and using the parameter and barcodes files. Note that paths and node architecture are specific to the AMNH cluster and my my account.

Scripts for the next processing steps are in the `02_ipyrad_steps_2_to_5` directory.

Following these steps, `03_ipyrad_branching_steps_6_7` contains scripts to create new ipyrad branches containing only: 1) the full set of L. getula samples, including low quality samples, 2) the set of L. getula samples used in most of analyses after dropping low quality samples, 3) eastern samples only, and 4) western samples only. These branching scripts use the names files to create the the params files, which are then used to run steps 6 and 7 of ipyrad and create the final datasets for each of these sets of individuals. `ipyrad_p123_Lgetula_67_v1.pbs` runs steps 6 & 7 on all samples, then remaining `.pbs` scripts run step 7 on the 3 datasets used in the manuscript.


## Population clustering, IBD, and GDM in R

Directory `04_R_analyses` contains R scripts used to perform analyses.

* `Lgetula_pop_assignment.R ` contains code to run DAPC and sNMF for population clustering.

* `Lgetula_pop_assignment_SUPP_fig.R` is a modification of the previous script to spit out a pdf with sNMF and DAPC results from various values of K instead of just the focal analyses included in the main body of the manuscript.

* `Lgetula_IBD.R` contains code to plot out kernel-density plots of geographic vs. genetic distance. Code to extract environmental data for each sample is also included at the end of this script. Environmental data are then fed into the following script to run GDM analyses.

* `Lgetula_GDM.R` will run generalized dissimilarity analysis on assemblies.




## FastSimCoal2

Directory `05_FSC2` contains scripts used to run FastSimCoal2 (FSC2) analyses.

This is complex to set up and run:


## First step is to make the file that indicates which individual is in which population


This is done using the script: `Lgetula_make_pops_for_SFS.R` This requires that sNMF was already run on the samples to estimate population clustering depends on the way things are set up in the `Lgetula_pop_assignment.R` script and follows the same conventions, including requiring several of the files from that script.

This has to be run somewhat interactively, so that naming of different populations can be done in a way that makes sense and is interpretable. I've already done this, and the resulting `Lgetula_p123_v2_25missFSCpops_k3.txt` file is on Dryad.



## Generate SFS


I then used Isaac Overcast's ultra-convenient `easySFS.py` script from here: [https://github.com/isaacovercast/easySFS](https://github.com/isaacovercast/easySFS). To determine how many alleles to project down (downsample to) to account for missing data
in k=3 assembly run, the first line below, the second is the projection I used:


```
	./easySFS.py -i Lgetula_p123_v2_25miss.vcf -p Lgetula_p123_v2_25missFSCpops_k3.txt --preview
	./easySFS.py -i Lgetula_p123_v2_25miss.vcf -p Lgetula_p123_v2_25missFSCpops_k3.txt --proj=48,14,10
```

The file output from here is on Dryad as `curMK3_MSFS.obs`.

# DRYAD LINK HERE AGAIN!!!!!


## Fixing divergence time

For FSC2 parameter estimation from the folded SFS (as used here), mutation rate or divergence time needs to be fixed. Pyron & Burbrink 2009 (mtDNA) has root at 4.91 Ma, Chen et al. 2017 looks like ~3-4. I used dates from Burbrink & Gehara tree "Tree1\_Beast" & "Tree2\_Beast" from https://datadryad.org/stash/dataset/doi:10.5061/dryad.4qs50. Tree1\_Beast = 4.33 and Tree2\_Beast = 4.49, the mean of these is 4.41. FSC works in generation times, so assuming a generation time of 3 years, so root divergence time, in generations = 1470000.


### Prepping and running FSC models for parameter estimation

* Note that this is the way that I set this up to allow me to use some scripts that I've written to manipulate files--if I was building this from the ground up now, I'd probably do some things differently and more efficiently, but this works.

1. Make the models. Making the models and ensuring that they are specified correctly is the most tedious part of this. Fortunately, I've done this already, so if you want to replicate exactly what I did, just look in the `Models` subdirectory within `05_FSC2` to see the FSC2 `.est` and `.tpl` files that I used. If you want to do this with your own data, you'll unfortunately have to make these files yourself manually--feel free to contact me for help, and I'm happy to give input as I'm able.

2. Get the linux executable for FSC: http://cmpg.unibe.ch/software/fastsimcoal2/ (for running this on a Linux cluster, as I did).

3. Put the `Models` directory (without any name changes or scripts below won't work) onto the cluster. The `Models` directory must contain a directory for each model to be estimated (each containing a `.tpl` and `.est` file) - the current `Models` directory does this - as well as the sfs file. You will need to add the sfs `.obs` file from Dryad if running this.

4. To run multiple replicates of FSC parameter estimation, use the `Prep_FSC_reps.sh` script. Details of how to use this script are in the comments inside it. Briefly, run it from inside the Models folder made in step 6, supplying an argument specifying if you are using either multi-dimensional SFS (argument MSFS) or two-dimensional SFS (argument jointMAF). This will make 50 replicates in a directory "Reps" in the parent directory of Models, each containing the necessary tpl, est, and obs file.

5. Move that "Reps" directory onto the cluster. Submit all of the jobs to the cluster as a big job array on a SLURM cluster using `FSC_Lgetula_p123_v2_25miss_FixRootTime_K3_Mods.slurm`.
	
6. For each of the models find the single run with the best likelihood. The script `Get_best_FSCacross_mods.sh` wraps the `fsc-selectbestrun.sh` script from here: [https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh](https://raw.githubusercontent.com/speciationgenomics/scripts/master/fsc-selectbestrun.sh) (and copied here) across each model. Comments in this script further describe its usage. To get AIC from the bestlhoods files, go inside `best_L_allMods directory`, created by last script, and run the script `Get_AIC_across_mods.R`(!!!NOTE: This AIC calculation will only work if you output 1 and exactly 1 parameter to this file for each estimated parameter otherwise the AIC values will be incorrect -- i.e., if you have a complex parameter, and you output both the estimated parameter and the complex parameter that is a transformation of the estimated parameter, you'll need a different solution here).

7. Convert the parameter estimates into useful units (e.g., years, number of migrants per generation, etc.) -- I still don't have a fully automated solution for this, so has to be done on an ad hoc basis depending on what parameters are included, etc. `Par_conv_FSC_getula.R` will do these conversions here - see comments in the script for more info.  This script will also start some prep for parametric bootstrap estimation of confidence intervals around parameter estimates.

8. Parametric bootstrapping for confidence intervals. This relies on having the `boot_input` directory made by the previous script. Use script `Prep_FSC_reps_of_bootreps.sh` in the directory that contains the directories each of which is a boot rep to generate 50 fsc reps of each bootstrap. `FSC_Lgetula_p123_v2_25miss_K3_boot.slurm` will then run parameter estimation on each bootstrap replicate.

9. Summarize bootstraps - for each bootstrap replicate, we need to get the best of the runs, pull the `bestlhoods` file out, and then summarize these across everything. `Get_best_FSCacross_boots.sh` will get the best run within each bootstrap rep. `Get_pars_across_bootreps.R` will then process these estimates.








