# Spatial Bootstrap Resampling Script
# Accompanies Manuscript: "Water Balance Creates a Threshold in Soil pH at the Global Scale"
# Eric Slessarev 11/20/2016

# This script generates a spatially random sample from point data
# the algorithm assumes a 180X360 global sampling grid (1 degree cells)
# sampling nodes are drawn from cell centers, with a distance 
# constraint to avoid edge biases. The algorithm then randomly selects a profile
# from the closest cell containing profiles.

#============================================================================
# INPUT DATA DESCRIPTION

# this function is applied to three data frames:
# "b" soil profile data at 0.5 m depth
# "a" soil profile data at 0.1 m depth
# "nb" NCSS profile data at 0.5 m depth
# It also requires "gm", a global matrix of environmental data.
# Gridded environmental data are stored as columns in these matrices,
# each grid cell has a unique index (cellindex), plus center coordinates.
# Each profile also has a unique index (MID).

# the following fields are required:

# cclon: longitude of cell center
# cclat: latitude of cell center
# cellindex: integer cell index, counting down each column from the west 
# SEA: binary where 1 signifies land, 0 signifies water
# MID: integer, sequential row index for a, b and nb


# example for a, b or nb:

# MID cellindex   DS      PID       cclon    cclat    pH    AWB_PM  ...
#  1    2906     NCSS  NCSS24576    -163.5    64.5    6.7    3.29   ...
#  2    2906     NCSS  NCSS24577    -163.5    64.5    6.6    3.29   ...
#  3    3082     NCSS  NCSS33677    -162.5    68.5    7.1    2.89   ...
# ...

# example for gm:

# cellindex   cclon    cclat    SEA    AWB_PM   ...
#     1      -179.5    89.5      0      NA      ...
#     2      -179.5    89.5      0      NA      ...
#     3      -179.5    89.5      0      NA      ...
#    ...

# These data are available on request
# (contact Eric Slessarev)

#load("C:/Users/.../soildata.RData")

#============================================================================


#============================================================================
# PARAMETERS

#choose a maximum length scale for sampling (km)
length_scale = 100

# choose a sample size
n = 20000

# resolution of the gridded datasets in decimal degrees
resolution = 1

#============================================================================


#============================================================================
# FUNCTIONS
# These functions are stored in separate scripts
# They should be sourced before using the spatial sample function

# geodist()
# uses the Haversine formula to obtain distance between one point and vector of points

# source("C:/Users.../geodist.R")

# cellarea()
# gets area of a grid cell using latitude and resolution

# source("C:/Users.../cellarea.R")
#============================================================================


#============================================================================
# SPATIAL SAMPLING ALGORITHM

# the algorithm is stored in a function that can be applied to a dataset in the format described above
# In outline it works like this:
# (1) loop through coordinates of pedons to generate a geographic search area constrained by length_scale
# (2) randomly sampl n nodes from the search area cell centers
# (3) loop through sampling nodes and find closest cell to each that contains pedons
# (4) randomly draw a pedon from the closest cell

spatial_sample = function(soildata,gm,length_scale,n){
  cat("\nrunning spatial sampling algorithm.\n")

  # all of the cell-center coordinates on landsurface
  cellcoordinates = gm[gm$SEA==1,c("cclon","cclat")]

  # all of the cell-centers that contain pits
  soildata_coords = soildata[!duplicated(soildata$cellindex),c("cellindex","cclon","cclat")]
  
  # logical to select cell-center coordinates within 100km of a cell with a profile
  search_area = logical(nrow(cellcoordinates))

  # loop 1: retrieve search area
  
  cat("\nretrieving search area\n\n")
  for(i in 1:nrow(soildata_coords)){
    
    # distances between cellcenter[i] and every cell center in the search area
    alldists_search =  geodist(soildata_coords$cclon[i],soildata_coords$cclat[i],
                             cellcoordinates$cclon,cellcoordinates$cclat)
    
    # logical defining the area within 100km
    ellipse = alldists_search <= length_scale
    
    # join to search area logical
    search_area = search_area|ellipse
    
    if(!as.logical(i%%500)){
      cat(round(100*i/nrow(soildata_coords)),"% retrieved...\n")
    }
  }
  cat("\n100 % retrieved.\n\n")

  # subset the cell centers that fall within the search area
  cellcoordinates = cellcoordinates[search_area,]
  
  # obtain a vector for the area of these cells
  cellareas = cellarea(cellcoordinates$cclat,resolution)

  # obtain sampling nodes
  # rs = "random sample". Sample cells from the land-surface with replacement (area weighted)
  rs =  cellcoordinates[sample(nrow(cellcoordinates),n,replace=T,prob=cellareas),]
  
  # loop 2: get pedons
  cat("\nsampling pits.\n\n")
  sampled_indices = c()
  for(p in 1:n){
    
    # distances between pedons and sample node (km)
    alldists = geodist(rs$cclon[p],rs$cclat[p],soildata_coords$cclon,soildata_coords$cclat)
    
    # cell index of closest cell containing pit
    closest_cell = soildata_coords$cellindex[which.min(alldists)]

    # pits in cell
    closest_pits = soildata$MID[soildata$cellindex==closest_cell]
      
    # random draw of pit from cell
    sampled_pit = closest_pits[sample(length(closest_pits),1,replace=T)]
    
    # append
    sampled_indices= c(sampled_indices,sampled_pit)
    
    if(!as.logical(p%%1000)){
      cat(round(100*p/n),"% sampled...\n")
    }
    

  }
  
  # get spatial sample
  ss = soildata[sampled_indices,]
  return(ss)
}
#============================================================================


#============================================================================


# a non spatial sample
bs = b[sample(nrow(b),20000,replace=T),] 

# spatial sample of 0.5 m data
ss = spatial_sample(b,gm,length_scale,n)

# spatial sample of 0.1 m data
ssa = spatial_sample(a,gm,length_scale,n)

# spatial sample of 0.5 m NCSS data
ssn = spatial_sample(nb,gm,length_scale,n)


