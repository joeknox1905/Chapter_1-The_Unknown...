# This script creates between group metrics (dimension level 2) in dispRity to calculate degree of the functional turnover (or divergence between trait spaces; distance between centroids?) and nestedness 
#(i.e. subset within a subset) to charcterise uniqueness of migrant vs residents in terms of the annual trait space. A seasonal analysis is then also run
# starting with a binary defintion of migratory status present year-round vs present <12 months
# and then testing migratory status as graded scale (from present 1 month to present 12 months) (bootstrapping/rarefying- how much does a metric of functional overlap/nestedness change when we remove each species in each group from the trait space?)

##############################################################################################################################################
#pre-requisites-----
##############################################################################################################################################


BI_traitspace<- read.csv("Data/Processed/eco_region_traitspace__PCA_BI.csv") # full Britain and Ireland trait space or beak-only AVONET PCA/PCOA scores, status and species (columns) by ecoregion (rows)- can either use PCOA or PCA

#specificy trait columns from trait space
trait_cols <- grep("^PC", names(BI_traitspace), value = TRUE)

#will also need seasonal trait space for seasonal comparisons (e.g. winter/summer/autumn migrants unique compared to residents 
season_eco_region_traitspace<- ("")

#@@@JK to TG: (PC scores are created from the full traitspace for Britain and Ireland) so for each ecoregion, migrant and resident group centroids will be relative to the full Britain and Ireland species traitspace- is this ok if I am focusing on local (ecoregional) functional diversity?


##############################################################################################################################################
#-------------------------------------------------------------dispRIty metrics----------------------------------------------------------------
##############################################################################################################################################


##------------Create metric for calcualting distance between resident and migrant centroids, two arguments (matrix 1 = res, matrix 2 = mig)------------

#Do migrants differ on average to residents in local (ecoregional) trait space (i.e. do they occupy different niches?)
#N.b. tells us nothing about the spread of species within their group or the nestedness of migrants and residents, could just be skewed by a few outliers



distance.between.centroids <- function(matrix, matrix2) {
    ## calculate the centroid of the matrices of n-dimensional traitspace
   sqrt(sum((colMeans(matrix) - colMeans(matrix2))^2))
} #sum of the euclidean distance between migrant group and resident centroids across multidimensional traitspace


## Test if disprity understands it
make.metric(pairwise.distance.between.centroids) #in this case should be level 1 (only one output) and between.groups


#-----------------------------metric for calculating displacement to centroid------------------------

#i.e. is migrant groups average trait combinations (centroid) at the center (centroid, 1) or away from centre (<1) of traitspace and therefore does it represent typical/redundant niches (generalist) or more unique niches (specialist)?
#more of a positional measure than a size measure (Guillerme et al., 2020). 
#also may be good to calcualate ecoregional migrant group displacement compared to full BI trait space not just annual trait space as would give an idea of whether 
#migrant group is unique locally but more redundant regionally- if that is the case, migrants are more functionally similar to migrants/residents in other ecoregions, suggesting either local migration from other ecoregions to the ecoregion (could confirm with migratory distance analysis)
#or that residents are actually more functionally unique in this ecoregion and that ecoregion shapes resident communities more than migrant communities.

displacements() in dispRity

  #calculated as:
  
  avg.group.displacements<- function(matrix) {
    reference_centroid<- colMeans(trait_cols)
    group_centroid<-  colMeans(matrix)
  
  sqrt(sum(group_centroid)^2/sqrt(sum(group_centroid-reference_centroid)^2) #@@@JK:TG is this correct?

#need to call between.groups = FALSE
result_binary  <- calculate.disprity.by.ecoregion(BI_traitspace, group = "status", metric = distance.between.centroids , trait_cols = trait_cols, between = FALSE) 



#----------------metric for calculating distribution of species within each group around centroid----------------


#i.e. "functional" arrangement within each group- niche packing vs expansion in each group
#if normal, means species are tend to be within their group (whatever trait combinations represent the centroid)
#if uneven, suggests clumpiness/dispersion from centroid

distribution.of.group.around.group.centroid <- function(matrix) {
group_centroid<- colMeans(matrix) #centroid of group
for (i in 1:nrow(matrix)) { #for each species row in ecoregional annual trait space
species_row <- matrix[i, ]#get pc score for each dimension for each pc

sqrt(sum((species_row - colMeans(matrix2))^2)) #distance from group centroid of each species

##@@@JK- then will need to plot histogram of distribution function and test difference using t test (not permutation as already have done bootstrapping to account for null tests)



#need to call between.groups = FALSE
result_binary  <- calculate.disprity.by.ecoregion(BI_traitspace, group = "status", metric = distance.between.centroids , trait_cols = trait_cols, between = FALSE) 

}



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
  
#Ellpisoid volume = π^(d/2) / Γ(d/2 + 1) × ∏ᵢ₌₁ᵈ √λᵢ

# Where:
 #π^(d/2) / Γ(d/2 + 1) is the volume of unit ball in d dimensions, Γ is the gamma function

#And: 
# ∏ᵢ₌₁ᵈ √λᵢ is the product of the ellipsoids semi-axis lengths, with √λᵢ  as the standard deviation along each pc axis

 
#test that metric works as expected
  make.metric(ellipse.contribution) #will be level 1 (one output) and not between groups as species-level. output should be a 

#call function
result_binary  <- calculate.disprity.by.ecoregion(BI_traitspace, group = "status", metric = distance.between.centroids , trait_cols = trait_cols, between = FALSE) 
#change for ellipsoid contribution (between = FALSE as species level)

# then would do stats e.g. function to get mean ellipse contribution per group from disprity summary and calculate difference between migrants and residents

calculate.mean.ellipse.contrib.sig <- function(disparity_result) {
disparity_result<- disparity_result %>% #need disparity summary from calculate.disprity.by.ecoregion with metric = ellipse.contribtion) 
group_by(group =group) %>%
summarise(mean_ellips_cont = mean(disparity_result$obs))  #get mean ellipse contribution (bootstrapped and rarefied) of each group

t.test((mean_ellips_cont ~ group, data = disparity_result )) # welches t test for non-parametric data

}

#call function
calculate.mean.ellipse.contrib.sig(disparity_result = result_binary) #or result_monthly

##############################################################################################################################################
#------------functions to calculate disprity metrics by ecoregion and summarise metric in a single table for all ecoregions------------------
##############################################################################################################################################

calculate.disprity.by.migrant.group <- function(d, group, metric, trait_cols, between,
                        min_n = 6, bootstraps = 1000) {
  d   <- distinct(d, species, .keep_all = TRUE) # get distinct species from 
  mat <- d %>% column_to_rownames("species") %>% select(all_of(trait_cols)) %>% as.matrix() #convert dataframe of pc scores to matrix for disprity
  grp <- split(rownames(mat), d[[group]]) %>% keep(~ length(.x) >= min_n) #detect migrant groups with at least 6 species (number of dimensions)-should be all, and extract as matrix
  if (length(grp) < 2) return(NULL) #for ecoregions with less than 2 migrant groups (should be none), do not perform calculation.
  custom.subsets(mat, group = grp) %>% 
    boot.matrix(bootstraps = bootstraps, rarefaction = min)) %>% #bootstrap with replacement (resample species from ecoregional pool for each migrant group and rarefy to the minimum number of species of a group in the ecoregion
    #rarefaction = min would work as resident vs migrant groups are not too skewed in species richness 
    disparity<- dispRity(metric = metric, between.groups = between) %>% in the group with the fewest species
    disp_result<- summary(disparity) #summarise disprity metric, (50, 95%) quantiles, median boostrap value)
}
#n.b. bootstrapping ensures greater confidence in metric for the group i.e. metric is due to a real biological difference, not just an artefact of sampling variance

calculate.disprity.by.ecoregion <- function(traitspace, site = "ecoregion", ...) {
  split(traitspace, traitspace[[site]]) %>% #list of pc score dataframes split by ecoregion
    map(calculate.disprity.by.migrant.group , ...) %>% #apply disprity by migrant and resident group to each ecoregion dataframe and return list of results
    bind_rows(.id = site) #bind rows of results from list so that disprity is summarised for all ecoregions and pairs in a dataframe
}


#get results e.g. for distance between group centroids:
result_binary  <- calculate.disprity.by.ecoregion(BI_traitspace, group = "status", metric = distance.between.centroids , trait_cols = trait_cols, between = TRUE) 


#get result monthly
result_monthly <- calculate.disprity.by.ecoregion(BI_traitspace, group = "n_year_month", metric = distance.between.centroids , trait_cols = trait_cols, between = TRUE) 
# for this analysis, the monthly migrant group pairs would also be compared against eachother which is what we want


#@@@ JK will need to update for seasons e.g. group by several variables but can do later

