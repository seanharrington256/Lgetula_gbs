#### Plot out the getula sNMF and DAPC results at various values of k
####        into a big pdf for supporting info for the paper


### load up relevant packages
library(adegenet)
library(LEA)
library(plotrix)
library(mapdata)
library(rworldmap)


## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/pop_assignment"
setwd(main_dir)

### Set this up so that major code blocks only need to be written once, with subsitutions based on specific species
######################################################################################################################
######################################################################################################################

## start with an object defining which assembly I'm working with - this is working with just 1 assembly
##    and writing several plots to one multi-page pdf

### NOTE!!!! #### 
## HERE THE ASSEMBLY NAME SHOULD EXACTLY MATCH THE iPYRAD PREFIX FOR OUTPUT FILES!!!!
##     e.g., if the .ugeno file is Lgetula_p123_v2_25miss.ugeno, then species should be defined as Lgetula_p123_v2_25miss


species<-"Lgetula_p123_v2_25miss"


######################################################################################################################
######################################################################################################################

  ###########################################################
  ## Set up paths to input files
  ###########################################################
    path_ugeno<-paste0(main_dir,"/", species,".ugeno")
    path_ustr<-paste0(main_dir,"/", species,".ustr")
    # path_usnps<-paste0(main_dir,"/", species,".usnps")
  
  
  ### Set up some colors for plotting farther down
  colors_2<-c("red", "blue") # colors for plotting 2 populations
  colors_3<-c("red", "blue", "white") # colors for plotting 3 populations
  colors_4<-c("red", "blue", "purple", "green") # colors for plotting 4 populations
  colors_5<-c("red", "blue", "purple", "green", "pink") # colors for plotting 5 populations
  colors_6<-c("red", "blue", "purple", "green", "pink", "yellow") # colors for plotting 6 populations
  colors_7<-c("red", "blue", "purple", "green", "pink", "yellow", "cadetblue2") # colors for plotting 7 populations
  
  ## Read in coordinates for plotting farther down
  setwd(main_dir)
  coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't
  # read in the geno file to get the number of individuals and snps for this assembly
  geno_txt<-readLines(path_ugeno)
  nums_snps<-length(geno_txt)
  num_ind<-length(strsplit(geno_txt[[1]], "")[[1]])
  
  
  
  ## quirk of read.structure function is that it requires the strucure file to have the file extension “.stru” - do some copying to make a new file with this extension
  path_stru<-gsub(".ustr", ".stru", path_ustr)  # Use a regular expression substitution to generate the new file name
  file.copy(path_ustr, path_stru) # make a copy of the file with the new name
  # Now we can read in this file
  DAPC_ustr<-read.structure(path_stru, n.ind=num_ind, n.loc=nums_snps, onerowperind = FALSE, col.lab=1, col.pop=0, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
  DAPC_ustr ## Take a quick look at how the data is structured for Adegenet
  ind_names<-rownames(DAPC_ustr@tab) ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting

  
  ####   this section will run DAPC interactively - I've commented it out to run everything with no user input  
  # ## find clusters & run DAPC
  # grp <- find.clusters(DAPC_ustr, max.n.clust=8, n.pca=num_ind-1)
  # dapc1 <- dapc(DAPC_ustr, grp$grp) # run DAPC
  # 
  # 
  # ## plot the DAPC the ugly way
  # scatter(dapc1, col=colors_4,  bg="white",
  #         legend=FALSE, posi.da = "bottomright",
  #         solid=.5
  # )
  # # another way
  # assignplot(dapc1)
  # # bar chart 
  # compoplot(dapc1, posi="bottomright",
  #           txt.leg=paste("Cluster", 1:length(grp$size)), lab="",
  #           n.col=1, xlab="individuals")
  
  ##########################################################################################################################
  ### Run sNMF  --- and below it DAPC at the same values of k for comparison, in each case, retaining all PCs
  ##########################################################################################################################
  
  # snmf requires the geno file to have the extension .geno - the geno file of unlinked snps has ugeno
  #   as above, copy the geno and make one with the extension .u.geno
  path_geno<-gsub(".ugeno", ".u.geno", path_ugeno)  # Use a regular expression substitution to generate the new file name
  file.copy(path_ugeno, path_geno) # do the copying with the new name
  
  
  # Run sNMF using 1 to 10 ancestral populations and evaluate the fit of different k values to the data using cross entropy criterion
  # before running snmf, check if it's already been run
  if(dir.exists(gsub("geno", "snmf", basename(path_geno)))){
    obj.at<-load.snmfProject(gsub("geno", "snmfProject", basename(path_geno))) # if it has, just load up the results
  }else{ # otherwise, run sNMF
    obj.at <- snmf(input.file = path_geno,  # input file is the .geno format file. We set up the path to this above
                   K = 1:10, # we will test for k=1 through 10
                   ploidy = 2,
                   entropy = T, # use the cross entropy criterion for assessing the best k value
                   repetitions = 10, # Run 10 independent replicate analyses
                   CPU = 2,
                   project = "new", tolerance = 0.00001, iterations = 500)
  }
  
  ## make a directory to put the output plots into if it doesn't already exist
  sNMF_out_dir<-paste0(main_dir, "/Pop_structr_out")  # specify a full path to the directory
  if(!dir.exists(sNMF_out_dir)){ # check if the directory  exists and then only create it if it does not
    dir.create(sNMF_out_dir)
  }
  setwd(sNMF_out_dir)
  
  # make pdf of cross-entropy plot
  # pdf(paste0(species, "_snmf_cross_ent.pdf"), width = 8, height=5)
  pdf(paste0(species, "_snmf_cross_ent.pdf"), width = 8, height=5)
  plot(obj.at, col = "lightblue", cex = 1.2, pch = 19, main="L. getula sNMF cross entropy")
  dev.off()
  
  # look at outstats
  outstats <- summary(obj.at)
  outstats # take a look
  
  # Plot k=2 through k=6 for all
  k_plot<-2:6

  ### For  down below, get the geographic coordinates sorted out
  ## make sure there aren't any individuals that don't have coordinates
  ind_names[which(!ind_names %in% coords[,"number"])]
  # match up the coordinates to the order of the individuals from snmf
  match_coords<-match(ind_names, coords[,"number"])
  snmf_coords<-coords[match_coords,]
  
  pdf(paste0(species, "_sNMF_DAPC_SUPP.pdf"), width = 8, height=9)
  #### use a loop to plot various different k values 
  for(i in k_plot){
    # confirm cross entropy values for K are consist. across runs
    ce <- cross.entropy(obj.at, K = i) 
    ce # pretty similar
    best.run <- which.min(ce) # find the run with the lowest cross validation error
    
    ## Get the snmf Q matrix from the best run at the best k
    qmatrix <- Q(obj.at, K = i, run = best.run)
    admix<-as.data.frame(qmatrix)
    
    ## Run a DAPC, too
    npcs<-num_ind-1  ## use max number of pcs, number of individuals-1
    ndas<-i-1  ## use max number of discriminant axes, which is k-1 (k is contained in i in this loop)
    grp_loop <- find.clusters(DAPC_ustr, n.pca=npcs, n.clust=i)## get groups
    dapc_loop <- dapc(DAPC_ustr, grp_loop$grp, n.pca=npcs, n.da=ndas)# run DAPC
    
    par(mfrow=c(2,1))
    par(mar=c(0.2,0.2,2,2))
    
    ## plot it out
    #### Start with plotting snmf to map
    map("worldHires", "Mexico", xlim=c(-125,-65), ylim=c(23,53),col="gray90", fill=TRUE)
    map("state", xlim=c(-125,-65), ylim=c(23,53), add=TRUE,col="gray90", fill=TRUE)
    # map("worldHires", "usa", xlim=c(-125,-65), ylim=c(23,53), add=TRUE,col="gray90", fill=TRUE) # swapped this out for the above line to get state boundaries
    text(x=-95, y=51, labels=paste0("L. getula SNMF K=", i), cex=0.6)
    ###Plot out the pies, with parameters set up depending on what k is
    if(i==2){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(admix$V1[x],admix$V2[x]), radius=0.6,
                                                   col=colors_2) }
    }
    if(i==3){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(admix$V1[x],admix$V2[x], admix$V3[x]), radius=0.6,
                                                   col=colors_3) }
    }
    if(i==4){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(admix$V1[x],admix$V2[x], admix$V3[x], admix$V4[x]), radius=0.6,
                                                   col=colors_4) }
    }
    if(i==5){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(admix$V1[x],admix$V2[x], admix$V3[x], admix$V4[x], admix$V5[x]), radius=0.6,
                                                   col=colors_5) }
    }
    if(i==6){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(admix$V1[x],admix$V2[x], admix$V3[x], admix$V4[x], admix$V5[x], admix$V6[x]), radius=0.6,
                                                   col=colors_6) }
    }
    if(i==7){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(admix$V1[x],admix$V2[x], admix$V3[x], admix$V4[x], admix$V5[x], admix$V6[x], admix$V7[x]), radius=0.6,
                                                   col=colors_7) }
    }
    ## plot out DAPC to map for same k
    map("worldHires", "Mexico", xlim=c(-125,-65), ylim=c(23,53),col="gray90", fill=TRUE)
    map("state", xlim=c(-125,-65), ylim=c(23,53), add=TRUE,col="gray90", fill=TRUE)
    text(x=-95, y=51, labels=paste0("L. getula DAPC K=", i), cex=0.6)
    if(i==2){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(dapc_loop$posterior[x,1], dapc_loop$posterior[x,2]), radius=0.6,
                                                   col=colors_2)}
    }
    if(i==3){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(dapc_loop$posterior[x,1], dapc_loop$posterior[x,2], dapc_loop$posterior[x,3]), radius=0.6,
                                                   col=colors_3)}
    }
    if(i==4){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(dapc_loop$posterior[x,1], dapc_loop$posterior[x,2], dapc_loop$posterior[x,3], dapc_loop$posterior[x,4]), radius=0.6,
                                                   col=colors_4)}
    }
    if(i==5){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(dapc_loop$posterior[x,1], dapc_loop$posterior[x,2], dapc_loop$posterior[x,3], dapc_loop$posterior[x,4], dapc_loop$posterior[x,5]), radius=0.6,
                                                   col=colors_5)}
    }
    if(i==6){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(dapc_loop$posterior[x,1], dapc_loop$posterior[x,2], dapc_loop$posterior[x,3], dapc_loop$posterior[x,4], dapc_loop$posterior[x,5], dapc_loop$posterior[x,6]), radius=0.6,
                                                   col=colors_6)}
    }
    if(i==7){
      for (x in 1:nrow(snmf_coords)) {floating.pie(snmf_coords$lon[x],snmf_coords$lat[x], # plot the pies
                                                   c(dapc_loop$posterior[x,1], dapc_loop$posterior[x,2], dapc_loop$posterior[x,3], dapc_loop$posterior[x,4], dapc_loop$posterior[x,5], dapc_loop$posterior[x,6], dapc_loop$posterior[x,7]), radius=0.6,
                                                   col=colors_7)}
    }
  }
  dev.off()


