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

* `Lgetula_GDM.R` will run generalized dissimilarity analysis.


