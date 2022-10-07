#### Script to do population assignment for Lampropeltis getula assemblies


### NOTE THAT YOU WILL NEED TO CHANGE FILE PATHS IF USING THIS SCRIPT!!!!!



### load up relevant packages
library(adegenet)
library(LEA)
library(plotrix)
library(mapdata)
library(rworldmap)
library(ggplot2)
library(scatterpie)



## Set up an object to contain the path to the main directory where the data are and then set that as the working directory
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/pop_assignment"
setwd(main_dir)

### Set this up so that major code blocks only need to be written once, with subsitutions based on specific assemblies
######################################################################################################################
######################################################################################################################

## start with an object defining which assemblies we're working on then loop over all of them
## the use of the object name "species" is a remnant of when I was doing this across different species
##  instead of different assemblies within a single species

### NOTE!!!! #### 
## HERE THE ASSEMBLY NAME SHOULD EXACTLY MATCH THE iPYRAD PREFIX FOR OUTPUT FILES!!!!
##     e.g., if the .ugeno file is Lgetula_p123_v2_25miss.ugeno, then species should be defined as Lgetula_p123_v2_25miss




all_assemblies<-c("Lgetula_p123_v2_25miss", "Lgetula_p123_v4_25miss", "Lgetula_p123_v5_WEST_25miss")



#### Some overall setup for mapping and plotting

## Get out the data for states and for Mexico and then combine them together
states <- map_data("state") # US states data
mex <- map_data("worldHires", "Mexico") # Mexico data
mex$group <- mex$group + length(states$group) # have to do this to get rid of weird lines that show up otherwise because of groups in Mexico already being group numbers in states
to_map <- rbind(states, mex) # combine these together
to_map <- dplyr::filter(to_map, lat > 23) # drop off souther coordinates that we don't need, only need northern Mexico

# make a list of colors:
colors_6<-c("V1" = "red", "V2" = "blue", "V3" = "white", "V4" = "purple", "V5" = "pink", "V6" = "yellow")

## Read in coordinates for plotting farther down
setwd(main_dir)
coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't


## make a directory to put the output plots into if it doesn't already exist
sNMF_out_dir<-paste0(main_dir, "/Pop_structr_out")  # specify a full path to the directory
if(!dir.exists(sNMF_out_dir)){ # check if the directory  exists and then only create it if it does not
  dir.create(sNMF_out_dir)
}

  
######################################################################################################################
######################################################################################################################

for(species in all_assemblies){
  ###########################################################
  ## Set up paths to input files
  ###########################################################
  setwd(main_dir)
    path_ugeno<-paste0(main_dir,"/", species,".ugeno")
    path_ustr<-paste0(main_dir,"/", species,".ustr")

  
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
  
  setwd(sNMF_out_dir)
  
  # make pdf of cross-entropy plot
  pdf(paste0(species, "_snmf_cross_ent.pdf"), width = 8, height=5)
  plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)
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
  
  
  #### use a loop to plot various different k values 
  for(i in k_plot){
    # confirm cross entropy values for K are consist. across runs
    ce <- cross.entropy(obj.at, K = i)
    ce # pretty similar
    best.run <- which.min(ce) # find the run with the lowest cross validation error
    
    ## Get the snmf Q matrix from the best run at the best k
    qmatrix <- Q(obj.at, K = i, run = best.run)
    admix<-as.data.frame(qmatrix)
    
    # get the coordinate and admix data into a single dataframe
    for_pies <- cbind(snmf_coords, admix)
    
    # Get the right number of colors
    colors <- colors_6[1:ncol(admix)]
    

    ## Run a DAPC, too
    npcs<-num_ind-1  ## use max number of pcs, number of individuals-1
    ndas<-i-1  ## use max number of discriminant axes, which is k-1 (k is contained in i in this loop)
    grp_loop <- find.clusters(DAPC_ustr, n.pca=npcs, n.clust=i)## get groups
    dapc_loop <- dapc(DAPC_ustr, grp_loop$grp, n.pca=npcs, n.da=ndas)# run DAPC
    
    # combine with coords
    colnames(dapc_loop$posterior) <- paste0("V", colnames(dapc_loop$posterior)) # so that this matches sNMF admix columns
    for_pies_dapc <- cbind(snmf_coords, dapc_loop$posterior)
    
    ## plot it out
    pdf(file=paste0(species,"_SNMF_DAPC_K",i,".pdf"), width=6, height=5)
    #### Start with plotting snmf to map
    snmf_plot <- ggplot(to_map, aes(long, lat, group = group)) + # map out the US & Mexico
      geom_polygon(data = to_map, fill = "grey90", color = "black", size = 0.2) + # make them polygons
      geom_scatterpie(data = for_pies, aes(x=lon, y=lat, group = number, r = 0.6), cols = grep("^V", colnames(for_pies), value = TRUE), size = 0.1) + # plot the pies - use grep to get the column names that start with V, these are the admix proportions
      scale_fill_manual(values = colors) +
      guides(fill="none") + # get rid of the legend for admixture
      theme_minimal() +
      labs(title=paste0(species, " SNMF K=", i), x ="Longitude", y = "Latitude") +
      coord_map("moll") # Mollweide projection
    print(snmf_plot)
    
    
    ## plot out DAPC to map for same k
    dapc_plot <- ggplot(to_map, aes(long, lat, group = group)) + # map out the US & Mexico
      geom_polygon(data = to_map, fill = "grey90", color = "black", size = 0.2) + # make them polygons
      geom_scatterpie(data = for_pies_dapc, aes(x=lon, y=lat, group = number, r = 0.6), cols = grep("^V", colnames(for_pies_dapc), value = TRUE), size = 0.1) + # plot the pies - use grep to get the column names that start with V, these are the admix proportions
      scale_fill_manual(values = colors) +
      guides(fill="none") + # get rid of the legend for admixture
      theme_minimal() +
      labs(title=paste0(species, " DAPC K=", i), x ="Longitude", y = "Latitude") +
      coord_map("moll") # Mollweide projection
    print(dapc_plot)
    dev.off()
  }
}









