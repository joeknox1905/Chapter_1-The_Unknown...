# This script creates between group metrics (dimension level 2) in dispRity to calculate degree of the functional turnover (or divergence between trait spaces; distance between centroids?) and nestedness 
#(i.e. subset within a subset) to charcterise uniqueness of migrant vs residents in terms of the annual trait space. A seasonal analysis is then also run
# starting with a binary defintion of migratory status present year-round vs present <12 months
# and then testing migratory status as graded scale (from present 1 month to present 12 months) (bootstrapping/rarefying- how much does a metric of functional overlap/nestedness change when we remove each species in each group from the trait space?)


#pre-requisites-----

BI_traitspace<- ("") # full Britain and Ireland trait space or beak-only AVONET PCA/PCOA scores, status and species (columns) by ecoregion (rows)- can either use PCOA or PCA

#will also need seasonal trait space for seasonal comparisons (e.g. winter/summer/autumn migrants unique compared to residents 
season_eco_region_traitspace<- ("")

#@@@JK to TG: (PC scores are created from the full traitspace for Britain and Ireland) so for each ecoregion, migrant and resident group centroids will be relative to the full Britain and Ireland species traitspace- is this ok if I am focusing on local (ecoregional) functional diversity?



#-----------------------dispRIty metrics-----------------------------------------:

##------------Create metric for calcualting distance between resident and migrant centroids, two arguments (matrix 1 = res, matrix 2 = mig)------------

#Do migrants differ on average to residents in local (ecoregional) trait space (i.e. do they occupy different niches?)
#N.b. tells us nothing about the spread of species within their group or the nestedness of migrants and residents


distance.between.centroids <- function(matrix, matrix2) {
    ## calculate the centroid of the matrices of n-dimensional traitspace
   sqrt(sum((colMeans(matrix) - colMeans(matrix2))^2))
} #sum of the euclidean distance between migrant group and resident centroids across multidimensional traitspace


## Test if disprity understands it
make.metric(pairwise.distance.between.centroids) #in this case should be level 1 (only one output) and between.groups


##--------------create metric for group uniqueness compared to ecoregional annual traitspace (measured as average distance of species from group centroid), can also use for the seasonal traitspace---------------

#Are migrants on average locally (ecoregionally) unique (distanced from average trait combinations) compared to annual trait space?
#tells us more about how migrants (and residents) differ from the combined average trait combinations of both groups (so would be more indicative of turnover and less nestedness rather than a few extremes increasing the distance between group centroids)
#might also be useful to compare to compare each migrant group to full Britain and Ireland functional trait space because groups may differ locally but not regionally, suggesting that migrants would be more locally important/keystone, butless so regionally 
#i.e. important to look at scale- species endemic to only one ecoregion and functionally unique (see below metric) would then potentially present problems for ecoregions if they are lost under cliamte-change induced shifts. 

group.distance.to.annual.centroid <- function(matrix, matrix2) {
  annual_centroid <- colMeans(matrix) #centroid of ecoregional annual trait space
  species_coords <- rowMeans(matrix2) #avg of pc coordinates for all species in group for each dimension
  

  return(sum(avg(sqrt(annual_centroid-species_coords)^2))) #average distance of each species from annual trait space centroid in n-dimensional traitspace s(i.e. averaged for all dimensions)
}
  
  make.metric(group.distance.to.annual.centroid) #test it performs as expected, should be one output and between groups


#------------------------metric for calculating species-level contribution to trait space using ellipsoid volume, based on leave one out approach of Cardoso et al. (2026) but using ellipsoid volumes as less sensitive to outliers (Guillerme et al., 2026)-----------------------------
  
#Are migrants locally (ecoregionally) functionally unique or redundant?
#tells us more about the species-level contribution (redundancy vs complementarity) to trait space and its removal shifts functional richness, less volume loss = less unique


  ellipse.contribution <- function (matrix) { #requires matrix1 (ecoregional annual trait space), matrix2 (ecoregional annual trait space-1)
  


    for (i in 1:nrow(matrix)) { #for each species row in ecoregional annual trait space
      species_row <- matrix[i, ]
      matrix2 <- matrix[-species_row] #remove species from trait space 
      # Calculate ellipsoid volume for each species
    }
    return(ellipsoid.volume(matrix)-(ellipsoid.volume(matrix2)) = #ellipsoid volume contribution of each species return ellipsoid volume of species trait space occupation i.e. their functional diversity- larger = more unique/ less redundant contribution as species will become more clustered, suggesting the species is further from other species in trait space

  }

 
#test that metric works as expected
  make.metric(ellipse.contribution) #will be level 1 (one output) and not between groups as species-level

#for ellipsoid.contribution metric will need to do some stats to look at whether migrant vs resident contribution differs in each ecoregional trait space

#------------function to calcualate disprity metrics by ecoregion and summarise in a single table for all ecoregions------------------


calculate.disprity.by.migrant.group <- function(d, group, metric, trait_cols, between,
                        min_n = 6, bootstraps = 100) {
  d   <- distinct(d, species, .keep_all = TRUE) # get distinct species from 
  mat <- d %>% column_to_rownames("species") %>% select(all_of(trait_cols)) %>% as.matrix() #convert dataframe of pc scores to matrix for disprity
  grp <- split(rownames(mat), d[[group]]) %>% keep(~ length(.x) >= min_n) #detect migrant groups with at least 6 species (number of dimensions)-should be all, and extract as matrix
  if (length(grp) < 2) return(NULL) #for ecoregions with less than 2 migrant groups (should be none), do not perform calculation.
  custom.subsets(mat, group = grp) %>% 
    boot.matrix(bootstraps = bootstraps, rarefaction = min(lengths(grp))) %>% #rarefy to the number of species in the group with the fewest species
    dispRity(metric = metric, between.groups = between) %>%
    summary() #summarise disprity metric
}

calculate.disprity.by.ecoregion <- function(traitspace, site = "ecoregion", ...) {
  split(traitspace, traitspace[[site]]) %>% #list of pc score dataframes split by ecoregion
    map(calculate.disprity.by.migrant.group , ...) %>% #apply disprity by migrant and resident group to each ecoregion dataframe and return list of results
    bind_rows(.id = site) #bind rows of results from list so that disprity is summarised for all ecoregions and pairs in a dataframe
}

#specificy trait columns from trait space
trait_cols <- grep("^PC", names(BI_traitspace), value = TRUE)



#get results e.g. for distance between group centroids:
result_binary  <- calculate.disprity.by.ecoregion(BI_traitspace, group = "status", metric = distance.between.centroids , trait_cols = trait_cols, between = TRUE)
result_monthly <- calculate.disprity.by.ecoregion(BI_traitspace, group = "n_year_month", metric = distance.between.centroids , trait_cols = trait_cols, between = TRUE) # for this analysis, the monthly migrant group pairs would also be compared against eachother which is what we want













run.eco_binary <- function(ecoregion, metric, between, min_n = 6) { #minimum number = 6 as this is the same as the number of pc dimensions in BI traitspace which each species PC scores are derived from

for (i in mig_groups_binary) { #for each ecoregion x species x status pc score dataframe
    df_eco <- mig_groups_binary %>% 
      filter(ecoregion == "ecoregion") #filter rows of traitspace to specific ecoregion 
  
}
  matrix_eco <- (df_eco) %>% #convert dataframe to ecoregional matrix with pc scores 
  column_to_rownames("species") %>% 
  select(-status, -ecoregion)
 
  eco_mig_groups <- keep(mig_groups_binary(df_eco), ~ length(.x) >= min_n) # needs at least as many species as number of dimensions (6 for pca) to run, this shoudln't drop any species
  eco_mig_boot_groups <- custom.subsets(matrix_eco, group = eco_mig_groups) %>%
    boot.matrix(bootstraps = 100, rarefaction = min(lengths(mig_groups_binary))) #rarefy by the group with the lowest number of species #change to mig_groups_monthly for monthly analysis
  dispRity(eco_mig_boot_groups, metric = metric, between.groups = between) %>%
  summary()
}





#------------------putting it together-------------------------------------------------

