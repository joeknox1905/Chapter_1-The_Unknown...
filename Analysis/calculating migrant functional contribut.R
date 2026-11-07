# This script creates between group metrics (dimension level 2) in dispRity to calculate degree of the functional turnover (or divergence between trait spaces; distance between centroids?) and nestedness 
#(i.e. subset within a subset) to charcterise uniqueness of migrant vs residents in terms of the annual trait space. A seasonal analysis is then also run
# starting with a binary defintion of migratory status present year-round vs present <12 months
# and then testing migratory status as graded scale (from present 1 month to present 12 months) (bootstrapping/rarefying- how much does a metric of functional overlap/nestedness change when we remove each species in each group from the trait space?)

%need to make simpler and fix lots of bugs/errors%

#pre-requisites-----

eco_region_traitspace<- ("") # full or beak-only AVONET PCA/PCOA scores, status and species (columns) by ecoregion (rows)- can either use PCOA or PCA


#------------function to get migratory status for each ecoregion based on binary definition------------------


detect.migrants.binary <- function(traitspace, eco) {
  out <- setNames(vector("list", length(eco)), eco)  #return list of rows of migrants and residents by binary year-round vs not-year-round defintion for each ecoregion
  
  for (i in eco) {
    d <- traitspace %>%
      filter(ecoregion == i) %>%
      distinct(species, status, .keep_all = TRUE) 
    
    out[[i]] <- list(
      migrant  = d %>% filter(status == "migrant")  %>% select(-ecoregion, -status) %>% column_to_rownames("species"),
      resident = d %>% filter(n_year_month == "resident")  %>% select(-ecoregion, -status) %>% column_to_rownames("species"))
  } 
  out 
}


#call function for each ecoregion
mig_groups_binary <- detect.migrants.binary(
  traitspace = eco_region_traitspace,
  eco = unique(eco_region_traitspace$ecoregion) #get lists for each ecoregion and assign to object for disprity
)


#---------------function for identifying migrants based on number of months present in ecoregion------------------------

detect.migrants.monthly <- function(traitspace, eco) {
out <- setNames(vector("list", length(eco)), eco) #list of rows of migrant by monthly status and residents for each ecoregion
  
  for (i in eco) {
    d <- traitspace %>%
      filter(ecoregion == i) %>%
      distinct(species, n_year_month, .keep_all = TRUE) 
    
    out[[i]] <- list( #list of matrices (one per migrant group)
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
) # #get lists for each ecoregion and assign to object for disprity




##------------Create metric for calcualting distance between resident and migrant centroids, two arguments (matrix 1 = res, matrix 2 = mig)------------

#Do migrants differ on average to residents in local (ecoregional) trait space (i.e. do they occupy different niches?)
#N.b. tells us nothing about the spread of species within their group or the nestedness of migrants and residents


pairwise.distance.between.centroids <- function(matrix1, matrix2) {
    ## calculate the centroid of the matrices of n-dimensional traitspace
    centroid1 <- colMeans(matrix1) #centroid of res
    centroid2 <- colMeans(matrix2) #centroid of mig
    ## measuring the distances between them
    return(sum(sqrt(centroid1^2-centroid2^2))) #euclidean distance
}


## Test if disprity understands it
make.metric(pairwise.distance.between.centroids) #in this case should be level 1 (only one output) and between.groups


#boot matrix to correct for sample size (species richness)
mig_boots_eco1 <- boot.matrix(, mig_groups_binary[eco1] = 100, # call, for each ecoregion, needs to be 130 for binary analysis as migrant groups contain up to 122 species in some ecoregions, but for monthly analysis will need to be 75 as species richness as high as 72. #Annual as reference?
                         
                          rarefaction = TRUE) #replace all elements for a sample each time

disparity(mig_boots_eco1,  metric = pairwise.distance.between.centroids, between.groups = TRUE) #or SES?

# repeat for each ecoregion

mig_boots_eco2 <- boot.matrix(, mig_groups_binary[eco2] = 100, 
                         
                          rarefaction = TRUE)

disparity(mig_boots_eco2,  metric = pairwise.distance.between.centroids, between.groups = TRUE) 

mig_boots_eco3 <- boot.matrix(, mig_groups_binary[eco3] = 100,
                         
                          rarefaction = TRUE) 

disparity(mig_boots_eco3,  metric = pairwise.distance.between.centroids, between.groups = TRUE) 

mig_boots_eco4 <- boot.matrix(, mig_groups_binary[eco4] = 100, 
                          rarefaction = TRUE) 

disparity(mig_boots_eco4,  metric = pairwise.distance.between.centroids, between.groups = TRUE) 

mig_boots_eco5 <- boot.matrix(, mig_groups_binary[eco5] = 100, 
                         
                          rarefaction = TRUE) 

disparity(mig_boots_eco5,  metric = pairwise.distance.between.centroids, between.groups = TRUE) 





#--------------create metric for group uniqueness compared to ecoregional annual traitspace (measured as average distance of species from group centroid), can also use for the seasonal traitspace---------------

#Are migrants on average locally (ecoregionally) unique (distanced from average trait combinations) compared to annual trait space?
#tells us more about how migrants (and residents) differ from the combined average trait combinations of both groups (so would be more indicative of turnover and less nestedness rather than a few extremes skewing the traitspace )
#might also be useful to compare to compare each migrant group to full Britain and Ireland functional trait space because groups may differ locally but not regionally, suggesting that migrants would be more locally important/keystone, butless so regionally 
#i.e. important to look at scale- species endemic to only one ecoregion and functionally unique (see below metric) would then potentially present problems for ecoregions if they are lost under cliamte-change induced shifts. 

group.distance.to.annual.centroid <- function(matrix1, matrix2) {
  annual_centroid <- colMeans(matrix1) #centroid of ecoregional annual trait space
  species_coords <- colSums(matrix2) #sum of pc coordinates for all species in group for each dimension
  

  return(sum(avg(sqrt(species_coords^2-annual_centroid^2)))) #average distance of each species from annual trait space centroid in n-dimensional traitspace s(i.e. summed for all dimensions)
}
  
  make.metric(group.distance.to.annual.centroid) #test it performs as expected, should be one output and between groups
  
  #call functions, will need to loop round each ecoregion and groups will be monthly for 2nd analysis
  dispRity(annual_trait_space, groups = mig_groups_binary, resident, metric = group.distance.to.annual.centroid, between.groups = TRUE)
  
  
# metric for calculating species-level contribution to trait space using ellipsoid volume, based on leave one out approach of Cardoso et al. (2026) but using ellipsoid volumes as less sensitive to outliers (Guillerme et al., 2026)-----------------------------
  
#Are migrants locally (ecoregionally) functionally unique or redundant?
#tells us more about the species-level contribution (redundancy vs complementarity) to trait space and its removal shifts functional richness, less volume loss = less unique


  ellipse.contribution <- function (matrix1) { #requires matrix1 (ecoregional annual trait space), matrix2 (ecoregional annual trait space-1)
  
    For (i in 1:nrow(matrix1)) { #for each species row in ecoregional annual trait space
      species_row <- matrix1[i, ]
      matrix2 <- matrix1[-i, ] #remove species from trait space 
      # Calculate ellipsoid volume for each species
    }
    return(ellipsoid.volume(matrix1)-(ellipsoid.volume(matrix2)) = #ellipsoid volume contribution of each species return ellipsoid volume of species trait space occupation i.e. their functional diversity- larger = more unique/ less redundant contribution as species will become more clustered, suggesting the species is further from other species in trait space

  }

  make.metric(ellipse.contribution) #will be level 1 (one output) and not between groups as species-level
  
  #call function, will be compared to ecoregional annual trait space
  dispRity(metric = ellipse.contribution (ecoregional_annual_trait_space)) #will be eigen as matrix is distance like for both PCOA and PCA
  
  #for ellipsoid.contribution metric will need to do some stats to look at whether migrant vs resident contribution differs in each ecoregional trait space