# This script creates between group metrics (dimension level 2) in dispRity to calculate degree of the functional turnover (or divergence between trait spaces; distance between centroids?) and nestedness 
#(i.e. subset within a subset) to charcterise uniqueness of migrant vs residents in terms of the annual trait space. A seasonal analysis is then also run
# starting with a binary defintion of migratory status present year-round vs present <12 months
# and then testing migratory status as graded scale (from present 1 month to present 12 months) (bootstrapping/rarefying- how much does a metric of functional overlap/nestedness change when we remove each species in each group from the trait space?)



#pre-requisites-----

eco_region_traitspace<- ("") # full or beak-only AVONET PCA scores, status and species (columns) by ecoregion (rows)- can either use PCOA or PCA
annual_traitspace <- ("")


#------------function to get migratory status for each ecoregion based on binary definition------------------


detect.migrants.binary <- function(traitspace, eco) {
  out <- setNames(vector("list", length(eco)), eco)
  
  for (i in eco) {
    d <- traitspace %>%
      filter(ecoregion == i) %>%
      distinct(species, n_year_month, .keep_all = TRUE)
    
    out[[i]] <- list(
      migrant_1  = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_2 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_3 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_4 = d %>% filter(n_year_month == 1 %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               migrant_5  = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               migrant_6 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               migrant_7 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               migrant_8 = d %>% filter(n_year_month == 1) %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               migrant_9  = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               migrant_10 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               migrant_11 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               resident = d %>% filter(n_year_month == 1 %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
                               )
  } 
  out #returns list of PC scores by ecoregion and binary migratory status
}


#call function
mig_groups_binary <- detect.migrants.binary(
  traitspace = eco_region_traitspace,
  eco = unique(eco_region_traitspace$ecoregion)
) # get lists for disprity


#---------------function for identifying migrants based on number of months present in ecoregion------------------------

detect.migrants.monthly <- function(traitspace, eco) {
out <- setNames(vector("list", length(eco)), eco)
  
  for (i in eco) {
    d <- traitspace %>%
      filter(ecoregion == i) %>%
      distinct(species, n_year_month, .keep_all = TRUE)
    
    out[[i]] <- list(
      migrant_1  = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_2 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_3 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_4 = d %>% filter(n_year_month == 1 %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_5  = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_6 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_7 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_8 = d %>% filter(n_year_month == 1) %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_9  = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_10 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      migrant_11 = d %>% filter(n_year_month == 1)  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      resident = d %>% filter(n_year_month == 1 %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
    )
  } 
  out #returns list of PC scores by ecoregion and monthly migratory status
}

#call function
mig_groups_monthly <- detect.migrants.monthly(
  traitspace = eco_region_traitspace,
  eco = unique(eco_region_traitspace$ecoregion)
) # get lists for disprity and assign object



#loop through each ecoregion

##------------Create metric for calcualting distance between resident and migrant centroids, two arguments (matrix 1 = res, matrix 2 = mig)------------
pairwise.distance.between.centroids <- function(matrix1, matrix2) {
    ## calculate the centroid of the matrices of n-dimensional traitspace
    centroid1 <- colMeans(matrix) #centroid of res
    centroid2 <- colMeans(matrix2) #centroid of mig
    ## measuring the distances between them
    return(sum(sqrt(centroid1^2-centroid2^2))) #euclidean distance
}


## Test if disprity understands it
make.metric(pairwise.distance.between.centroids) #in this case should be level 1 (only one output) and between.groups


## use in dispRity, will need to create a for loop to loop round each ecoregion

custom.subsets()

#boot matrix to correct for sample size (species richness)
mig_boots <- boot.matrix(, mig_groups = 100, # needs to be 130 for binary analysis as migrant groups contain up to 122 species in some ecoregions, but for monthly analysis will need to be 75 as species richness as high as 72
                         
                          rarefaction = TRUE) #replace all elements for a sample each time

disparity(mig_boots,  metric = pairwise.distance.between.centroids, between.groups = TRUE) #or SES?

#for species contribution calculations (avg per group?)





#--------------create metric for group uniqueness compared to annual traitspace (measured as average distance of species from group centroid), can also use for the seasonal traitspace---------------

group.distance.to.annual.centroid <- function(matrix1, matrix2) {
  annual_centroid <- colMeans(matrix1) #centroid of annual trait space
  species_coords <- colSums(matrix2) #sum of pc coordinates for all species in group for each dimension
  

  return(sum(avg(sqrt(species_coords^2-annual_centroid^2)))) #average distance of each species from annual trait space centroid in n-dimensional traitspace s(i.e. summed for all dimensions)
}
  
  make.metric(group.distance.to.annual.centroid) #test it performs as expected, should be one output and between groups
  
  #call functions, will need to loop round each ecoregion and groups will be monthly for 2nd analysis
  dispRity(eco_region_traitspace, groups = migrant, resident, metric = group.distance.to.annual.centroid, between.groups = TRUE)
  
  
# metric for calculating species-level contribution to trait space using ellipsoid volume, based on leave one out approach of Cardoso et al. (2026) but using ellipsoid volumes as less sensitive to outliers (Guillerme et al., 2026)
  
  ellipse.contribution <- function(species_row, matrix1, matrix2) {
    species_row <- [.,matrix1]
    annual_centroid <- colMeans(matrix1) #centroid of annual trait space
    species_coords <- colSums(matrix2) #sum of pc coordinates for all species in group for each dimension
    
    
    return(sum(avg(sqrt(species_coords^2-annual_centroid^2)))) #average distance of each species from annual trait space centroid in n-dimensional traitspace s(i.e. summed for all dimensions)
  }
  
  