#### Script to make population assignment files to use in making SFS files for FastSimCoal

# requires that Lgetula_pop_assignment.R has already been run to generate snmmf output
#     many of the conventions in this script are copied from that script

## Note that your paths won't match mine and will need to be changed



### load up relevant packages
# library(adegenet)
library(LEA)
library(vcfR)
library(plotrix)
library(mapdata)
library(rworldmap)



main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/pop_assignment"
FSC_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/FSC"
setwd(main_dir)


# Set up assembly
species<-"Lgetula_p123_v2_25miss"
k<-3



obj.at<-load.snmfProject(paste0(species,".u.snmfProject")) # load up sNMF results

## Read in genetic data just to get names
path_vcf<-paste0(main_dir,"/", species,".vcf")
gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
ind_names<-gendata@ind.names ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting



#### For a selected value of k, get the best snmf run - k is specified up in if() statements showing file paths
ce <- cross.entropy(obj.at, K = k)
best.run <- which.min(ce) # find the run with the lowest cross validation error
qmatrix <- Q(obj.at, K = k, run = best.run)   # Get the snmf Q matrix from the best run at the best k
cluster <- apply(qmatrix, 1, which.max) # use the qmatrix to find which population cluster each sample has most membership in

#### Combine individual names with the population membership
fsc_pops<-cbind(ind_names, cluster)

## read in coordinates and plot so we know what to name each
coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't

## sort the coords into the right order
ind_names[which(!ind_names %in% coords[,"number"])]
# match up the coordinates to the order of the individuals from snmf
match_coords<-match(ind_names, coords[,"number"])
snmf_coords<-coords[match_coords,]

col_plot <- fsc_pops # make this into an object to specify colors

if(k==2){
  col_plot[which(col_plot[,2]==1),2]<-"white"
  col_plot[which(col_plot[,2]==2),2]<-"black"
}
if(k==3){
  col_plot[which(col_plot[,2]==1),2]<-"white"
  col_plot[which(col_plot[,2]==2),2]<-"black"
  col_plot[which(col_plot[,2]==3),2]<-"red"
}
if(k==4){
  col_plot[which(col_plot[,2]==1),2]<-"white"
  col_plot[which(col_plot[,2]==2),2]<-"black"
  col_plot[which(col_plot[,2]==3),2]<-"red"
  col_plot[which(col_plot[,2]==4),2]<-"purple"
}
if(k==5){
  col_plot[which(col_plot[,2]==1),2]<-"white"
  col_plot[which(col_plot[,2]==2),2]<-"black"
  col_plot[which(col_plot[,2]==3),2]<-"red"
  col_plot[which(col_plot[,2]==4),2]<-"purple"
  col_plot[which(col_plot[,2]==5),2]<-"blue"
}




map("worldHires", "Mexico", xlim=c(-125,-65), ylim=c(23,53),col="gray90", fill=TRUE)
map("state", xlim=c(-125,-65), ylim=c(23,53), add=TRUE,col="gray90", fill=TRUE)
points(x=snmf_coords[,"lon"], y=snmf_coords[,"lat"], pch=21, bg=col_plot[,"cluster"])


# I did initially play around with other values of K - only used K=3 for the manuscript

if(species=="Lgetula_p123_v2_25miss" && k==2){
  ## here, white = 1 = west - so set this up:
  fsc_pops[which(fsc_pops[,2]==1),2]<-"west"
  fsc_pops[which(fsc_pops[,2]==2),2]<-"east"
}
if(species=="Lgetula_p123_v2_25miss" && k==3){
  ## set up the populatoins based on how it was colored:
  fsc_pops[which(fsc_pops[,2]==1),2]<-"nmex"
  fsc_pops[which(fsc_pops[,2]==2),2]<-"west"
  fsc_pops[which(fsc_pops[,2]==3),2]<-"east"
}
if(species=="Lgetula_p123_v2_25miss" && k==4){
  ## set up the populatoins based on how it was colored:
  fsc_pops[which(fsc_pops[,2]==1),2]<-"east"
  fsc_pops[which(fsc_pops[,2]==2),2]<-"nmex"
  fsc_pops[which(fsc_pops[,2]==3),2]<-"west"
  fsc_pops[which(fsc_pops[,2]==4),2]<-"cent"
}
if(species=="Lgetula_p123_v2_25miss" && k==5){
  ## set up the populatoins based on how it was colored:
  fsc_pops[which(fsc_pops[,2]==1),2]<-"east"
  fsc_pops[which(fsc_pops[,2]==2),2]<-"cent"
  fsc_pops[which(fsc_pops[,2]==3),2]<-"miss"
  fsc_pops[which(fsc_pops[,2]==4),2]<-"nmex"
  fsc_pops[which(fsc_pops[,2]==5),2]<-"west"
}


setwd(FSC_dir)

## write the file - this has to be done in a convoluted way to prevent a newline character on the last line, which FSC reads as a final blank line and then will not run correctly--but will also not throw an informative error-I learned this the hard way
lines <- apply(fsc_pops, 1, function(x) paste(x[[1]], x[[2]], sep="\t"))
lines[1:length(lines)-1] <- paste0(lines[1:length(lines)-1], "\n")
cat(lines, file=paste0(species, "FSCpops_k", k, ".txt"), fill=FALSE, sep="")

