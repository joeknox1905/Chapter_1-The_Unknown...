








###################################################################################################################################################################################################
#------------calculate disprity metrics by ecoregion, status and season for binary status and monthly status defintiions and then summarise metric in a single table for all ecoregions------------------
###################################################################################################################################################################################################



traitspace_df<- read.csv("Data/Processed/eco_region_traitspace_PCA_BI.csv") 


traitspace_mat <- as.matrix(traitspace_df[, 8:13]) #convert to matrix and extract PC scores

#yearly analysis-----------------------------

#create migrant and resident groups by ecoregion from row indexes

groups_yearly <- split(
  seq_len(nrow(traitspace_df)), 
  list(traitspace_df$ecoregion, traitspace_df$status), #get rows of all species in each ecoregion and migratory status group
  drop = TRUE
)



#pass to custom.subsets to make groups
subsets_yearly <- custom.subsets(data = traitspace_mat, group = groups_yearly) 


plot(subsets_yearly, type ="preview") #PC1 V2


plot(subsets_yearly, type="preview", specific.args = list(dimensions = c(3,4))) #PC1 v3

#boot and rarefy subsets
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


# so looks like european atlantic forest especially, and then english lowland and caldeon conifer forests have a large distance between centroids, suggesting divergent average trait combinations
# n these ecoregions- would have to llok into ecosystem types, the diets, phylogeny and kinds of birds to continue to subset first, also look across groups and find a better way of plotting








#seasonal analysis-------------------------



#create migrant and resident group by ecoergion and season


groups_seasonally <- split(
  seq_len(nrow(traitspace_df)), 
  list(traitspace_df$ecoregion, traitspace_df$status,  traitspace_df$status), #get rows of all species in each ecoregion and migratory status group
  drop = TRUE
)

subsets_seaosonally <- custom.subsets(data = traitspace_mat, group = groups_seasonally) 


plot(subsets_yearly, type ="preview") #PC1V2  


plot(subsets_yearly, type="preview", specific.args = list(dimensions = c(3,4)) PC3V4






