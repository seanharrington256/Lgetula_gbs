## Complex cycles of divergence and migration shape population structure in the common kingsnake species complex

### Harrington and Burbrink

## IN REVIEW

This repo contains contains code used for analyses in the manuscript "Complex cycles of divergence and migration shape population structure in the common kingsnake species complex"

Raw sequence data are available from the NCBI Sequence Read Archive. Note that these data are demultiplexed data for each sample. Our full set of raw data included samples from other snake species for different projects in addition to the samples used here. We have included the scripts and parameter files that we used to demultiplex the data for the sake of completeness and transparency, but note that you will not be able to replicate these initial steps.

### ipyrad data processing

Demultiplexing scripts and files are in the `01_ipyrad_step_1_Demux` directory. `.pbs` scripts were used to run demultiplexing and merging of plates on the American Museum of Natural History (AMNH) computing cluster and using the parameter and barcodes files. Note that paths and node architecture are specific to the AMNH cluster and my my account.

Scripts for the next processing steps are in the `02_ipyrad_steps_2_to_5` directory.

Following these steps, `03_ipyrad_branching_steps_6_7` contains scripts to create new ipyrad branches containing only: 1) the full set of L. getula samples, including low quality samples, 2) the set of L. getula samples used in most of analyses after dropping low quality samples, 3) eastern samples only, and 4) western samples only. These branching scripts use the names files to create the the params files, which are then used to run steps 6 and 7 of ipyrad and create the final datasets for each of these sets of individuals.




