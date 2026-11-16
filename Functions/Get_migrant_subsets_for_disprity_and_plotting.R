###################################################################################################################################################################################################
#------------functions to calculate disprity metrics by ecoregion, status and season for binary status and monthly status defintiions and then summarise metric in a single table for all ecoregions------------------
######################################################################################################################################################################################################

calculate.disprity.by.migrant.group <- function(d, group, metric, trait_cols, between,
                        min_n = 6, bootstraps = 1000) {
  d   <- distinct(d, species, .keep_all = TRUE) # get distinct species from 
  mat <- d %>% column_to_rownames("species") %>% select(all_of(trait_cols)) %>% as.matrix() #convert dataframe of pc scores to matrix for disprity
  grp <- split(rownames(mat), d[[group]]) %>% keep(~ length(.x) >= min_n) #detect migrant groups with at least 6 species (number of dimensions)-should be all, and extract as matrix
  if (length(grp) < 2) return(NULL) #for ecoregions with less than 2 migrant groups (should be none), do not perform calculation.
  custom.subsets(mat, group = grp) %>% 
    boot.matrix(bootstraps = bootstraps, rarefaction = min)) %>% #bootstrap with replacement (resample species from ecoregional pool for each migrant group and rarefy to the minimum number of species of a group in the ecoregion
    #rarefaction = min would work as resident vs migrant groups are not too skewed in species richness 
    disparity<- dispRity(metric = metric, groups= groups, between.groups = between) %>% in the group with the fewest species
    disp_result<- summary(disparity) #summarise disprity metric, (50, 95%) quantiles, median boostrap value)
}
#n.b. bootstrapping ensures greater confidence in metric for the group i.e. metric is due to a real biological difference, not just an artefact of sampling variance

#for yearly analysis
calculate.disprity.by.ecoregion <- function(traitspace, site = "ecoregion", ...) {
  split(traitspace, traitspace[[site]]) %>% #list of pc score dataframes split by ecoregion
    map(calculate.disprity.by.migrant.group , ...) %>% #apply disprity by migrant and resident group to each ecoregion dataframe and return list of results
    bind_rows(.id = site) #bind rows of results from list so that disprity is summarised for all ecoregions and pairs in a dataframe
}


#for seasonal analysis with migrants grouped by seasonal presence
calculate.disprity.by.ecoregion.and.season <- function(traitspace, site = "ecoregion", season ="season"...) {
  split(traitspace, traitspace[[site, season]]) %>% #list of pc score dataframes split by ecoregion
    map(calculate.disprity.by.migrant.group , ...) %>% #apply disprity by migrant and resident group to each ecoregion dataframe and return list of results
    bind_rows(.id = site, season) #bind rows of results from list so that disprity is summarised for all ecoregions and pairs in a dataframe
}
