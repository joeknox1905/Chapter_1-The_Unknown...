




###################################################################################################################################################################################################
#------------calculate disprity metrics by ecoregion, status and season for binary status and monthly status defintiions and then summarise metric in a single table for all ecoregions------------------
###################################################################################################################################################################################################

#load packages
library(dispRity)


#load functions

#making groups function
make.ecoregional.migrant.groups <- function(traitspace_df, by) {
  lapply(split(traitspace_df$species, traitspace_df[by], drop = TRUE), unique) #get list of species per ecoregion, status (and season for seasonal analysis) by splitting dataframe
}

#calcualting distance between centroids function
pairwise.distance.between.centroids <- function(matrix, matrix2) {
    ## calculate the centroid of the matrices of n-dimensional traitspace
   sqrt(sum((colMeans(matrix) - colMeans(matrix2))^2))
}

#calcualting quantile overlap function (proportional overlap, does not identify nestedness only how much of traitspace is common to BOTH groups) so a good start- i.e. can answer, are groups occupying different or common areas of trait space?
quantile.overlap<- function (matrix, matrix2, probs = c(0.25, 0.975)) {
  q1 <- apply(matrix, 2, quantile, probs = probs) 
  q2<- apply(matrix2, 2, quantile, probs = probs)
  intersect <- pmax(0, pmin(q1[2, ], q2[2, ]) - pmax(q1[1, ], q2[1, ])) 
  span <- pmax(q1[2, ], q2[2, ]) - pmin(q1[1, ], q2[1, ])
  intersect / span 
}

# LOad traitspace and convert to matrix

traitspace_df<- read.csv("Data/Processed/eco_region_season_traitspace_PCA_BI.csv") 

pc_cols <- names(traitspace_df)[8:13] #get names of PC columns excluding metadata
species_pcs <- unique(traitspace_df[c("species", pc_cols)]) #get unique pc cols per species # i.e. the full trait space, 1 row per species

#convert to matrix for disprity
traitspace_mat <- as.matrix(species_pcs[, pc_cols]) 
rownames(traitspace_mat) <- species_pcs$species 



#yearly analysis: Do migrants and residents differ in functional space occupation/contribution?-----------------------------

#see ellipse plots for preliminary plotting (mostly small overlaps, close to circle in a circle situation (see nestedness below))

#call function to split groups by ecoregion and migratory status
groups_yearly   <- make.ecoregional.migrant.groups(traitspace_df =traitspace_df, by = c("ecoregion", "status"))


#feed into custom subsets
subsets_yearly  <- custom.subsets(traitspace_mat, group = groups_yearly) #make subsets


#visualise
plot(subsets_yearly, type ="preview") #PC1 V2


plot(subsets_yearly, type="preview", specific.args = list(dimensions = c(3,4))) #PC1 v3

#boot and rarefy subsets to minimum species richness in any group for comparable ecoregional/season/migratory status comparisons
boot_yearly<- boot.matrix(subsets_yearly, bootstraps = 100, rarefaction = "min") #bootstrap with replacement (resample species from ecoregional pool for each migrant group and rarefy to the minimum number of species of a group in the ecoregion


disp_result<- dispRity(boot_yearly, metric = pairwise.distance.between.centroids, between.groups = TRUE) 
# this will compare every ecoregional status group with every other ecoregional status group- which is good because will get effect of ecoregion as well


summary(disp_result)# mean observed disparity per group is reported + bootstrapped median + confidence intervals

pdf(file = "disp_result_dist.between.centroids.pdf",  ) #disrtibution of bootstrapped centroid dist value + median for each pairwise comparison
    width = 4, 
    height = 4)

plot(disp_result, las=2) #boxplot of results
dev.off()


name.subsets(disp_result) #5x5 ecoregion x migratory status groups
size.subsets(disp_result) #correctly rarefied to minimum 60 species per subset


#would probably have to do some kind of test to establish whether the migrant and resident groups differ in each ecoregion- but if we are interested in functional contribution of co-existing species in ecoregions, will probably need to look at each seperately

#but first how do the observed values between centroids of each migrant and resident group within each ecoregion look?

#get observed matrix
get.disparity(disp_result) #will have to extract values from here e.g.


Caledon conifer forests.migrant:Caledon conifer forests.resident = 1.23, bootstrapped = 1.28, rarefied = 1.35
Celtic broadleaf forests.migrant:Celtic broadleaf forests.resident = 0.70, bootstrapped = 0.76, rarefied = 0.88
English Lowlands beech forests.migrant:English Lowlands beech forests.resident = 1.33, bootstrapped = 1.37, rarefied = 1.43
European Atlantic mixed forests.migrant:European Atlantic mixed forests.resident = 2.79, boostrapped =2.81, rarefied = 2.80 # will be similar across all as rarefied to this groups minimum species
North Atlantic moist mixed forests.migrant:North Atlantic moist mixed forests.resident = 0.78, bootstrapped = 0.84, rarefied = 0.94

#in this case looks like observed values are pretty close to boostrapped and rarefied values but will use booted and rared values to be safe
boot_rare_disp<- <- get.disparity(disp_result, observed = FALSE,
rarefaction = TRUE)

# so looks like european atlantic forest especially, and then english lowland and caldeon conifer forests have a large distance between centroids, suggesting divergent average trait combinations between migrants and residents
# in these ecoregions- would have to look into ecosystem types, the diets, phylogeny and kinds of birds to continue to subset first, also look across groups and find a better way of plotting

#or table of result

pdf(file= pairwise_centroid_dist_result_binary_yearly.pdf)
xtabs(get.disparity(disp_result, observed = FALSE,
rarefaction = TRUE) ~ disp_result$subsets)
dev.off()



# (proportional overlap via quantiles)


# nestedness is one group a complete subset of another group- would need directionality 
#i.e. which group is a subset of which group

#seasonal analysis-------------------------



#create migrant and resident group by ecoergion and season

groups_seasonally <- make.ecoregional.migrant.groups(traitspace_df = traitspace_df, by = c("ecoregion", "season", "status")) #seasonal analysis has additional season split

subsets_seaosonally <- custom.subsets(data = traitspace_mat, group = groups_seasonally) 


plot(subsets_yearly, type ="preview") #PC1V2  


plot(subsets_yearly, type="preview", specific.args = list(dimensions = c(3,4)) PC3V4






