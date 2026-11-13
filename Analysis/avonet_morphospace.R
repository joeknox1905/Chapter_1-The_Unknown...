# This script first creates a full "morphospace" for the 195 species present in britain and ireland using aoh-refined Birdlife range maps based on 11 morphological traits in AVONET
# HWI, body mass, beak depth, beak width, beak length from the culmen, beak length from the nares, wing length, length to the secondary flight feather, tail length and tarsus length (Keyser et al., 2024)
# beak traits = trophic traits linked to size and type of food consumed, tarsus length, wing traits
# and tail length = locomotory traits linked to foraging mode, dispersal and microhabitat use, body mass is linked to metabolism and foraging behaviour (Stewart et al., 2022).
# HWI may also be useful for calculating dispersal buffers later (Stewart et al., 2022).
# Then, following Williams and Curley (2025), to control for size, it takes the residuals of the averaged and logged PC traits against
# logged body mass to control for the strong correlation of body mass with size and better distinguish between species with similar body masses but different
# traits (Weeks et al., 2026). These size-relative traits are used for all traits except body mass (for obvious reasons).
# It a beak only traitspace and a non-beak only trait space to assess which traits are important in characterising the avian trait space and the migrant vs resident traitspace.
# These comparisons are repeated excluding HWI and kipps distance (as in Keyser et al., 2022) as likely highly correlated and shown to skew PC1 and 2 comapred to non-beak traitspace.
# see chapter 1 script)


#-----load libraries------

library(tidyr)
library(stringr)
library(kableExtra)
library(dispRity)
library(stats)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(ggcorrplot)
library(tibble)
library(ggfortify)

assemb <- readRDS("Data/Processed/seasonal_ecoregion_pools_aoh") # get monthly presence of each species for full annual trait space

# load avonet trait data
avonet <- read.csv("Data/Raw/AVONET1_BirdLife.csv") # trait means in Birdlife V7
avonet_crosswalk <- read.csv("Data/Raw/AVONET_Raw_Data.csv") # raw traits for taxonomic crosswalk

# standardise species name format
avonet <- avonet %>% mutate(species = str_replace_all(Species1, " ", "_"))


# get all species from Britain and Ireland aoh maps
spp_name <- sort(unique(assemb$species))

# species in monthly presence data with no trait data from AVONET
no_avonet <- setdiff(spp_name, avonet$species) %>%
  print()

# view these species in avonet crosswalk to manually correct names
no_avonet <- str_replace_all(no_avonet, "_", " ")
avonet_crosswalk %>%
  filter(if_any(everything(), ~ .x %in% no_avonet)) %>%
  View()

# Old Birdlife names to use to get traits from AVONET database, retaining Birdlife V10 names
# Also get taxonomy match from Rob in future
avonet_lookup <- c(
  Curruca_communis = "Sylvia_communis",
  Curruca_curruca = "Sylvia_curruca",
  Curruca_undata = "Sylvia_undata",
  Luscinia_svecica = "Cyanecula_svecica"
)

# join all ten traits to all species in British Isles,
avonet_traits <- tibble(species = spp_name, ) %>%
  mutate(lookup_name = coalesce(avonet_lookup[species], species)) %>%
  left_join(avonet %>% select(
    species, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
    Tarsus.Length, Wing.Length, Secondary1, Tail.Length, Mass, Kipps.Distance # JK:TG<- kipps distance is also a derived variable (wing length - secondary1); should I remove?
  ), by = c("lookup_name" = "species")) %>%
  dplyr::select(-lookup_name)

# get species with no trait data in avonet (should be 0)
avonet_na_any <- avonet_traits %>%
  filter(if_any(-species, is.na)) %>%
  print()
avonet_na_all <- avonet_traits %>%
  filter(if_all(-species, is.na)) %>%
  print()


# control for size in morphological traits (Williams and CUrley, 2025 residuals method)------------------

# AVONET morphological traits- removing HWI as it is a ratio
morph_traits <- avonet_traits %>%
  select(
    species, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
    Tarsus.Length, Wing.Length, Secondary1, Tail.Length, Mass, Kipps.Distance
  ) %>%
  column_to_rownames("species")

# Log body mass for function call
LogMass <- log(morph_traits$Mass)

# Log and transform morphological traits relative to body mass (residuals)
get.residuals <- function(x) {
  LogTrait <- log(x)
  ResTrait <- lm(LogTrait ~ LogMass)$residuals
}

res_traits <- apply(morph_traits[, -9], 2, get.residuals) # apply to all columns except body mass
colnames(res_traits) <- paste0("Res", colnames(res_traits)) # rename to reflective relative trait values
my_residual_data <- cbind(res_traits, LogMass)

# TG:done Note that the two bits above can probably be functionalised. Like using something like my_residual_data <- apply(my_data, 2, get.residuals)


# Check correlations among logged morphological traits
ggcorrplot(cor(my_residual_data[, 1:10]), hc.order = TRUE)

# body mass not correlated with anything so effectively removes allometry
# main correlations (can use for PCOA groups):
# beak length traits
# beak shape traits (width +depth)
# wing length and pointiness (kipps distance)
# tail length and first flight feather (secondary 1)

#-----------------------Perform PCA on full trait space--------------------------------
pca <- prcomp(my_residual_data, scale. = TRUE) # scale

full_space_axes <- select.axes(pca, threshold = 0.95, inc.threshold = TRUE) # select axes using disprity
# 7 dimensions needed
# TG:done might be good to run this one with your groups. e.g. pca <- custom.subsets(...) so that we can check if 6D is enough for each group.

summary(full_space_axes) # PCs1:7 explain 97.7% of variance

pdf("PCA_axes_variance_explained_full_BI_space", width = 14, height = 10)
plot(full_space_axes)
dev.off()


# create dataframe with tibble for joining metadata
avospace <- as_tibble(pca$x, rownames = "species") %>%
  select(species, starts_with(paste0("PC")))

# check that dimensions are ok for groups in yearly migrant contribution analysis---

# join metadata (migratory status, monthly presence)
eco_pca <- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(avospace, by = "species")

eco_mig_group <- eco_pca %>%
  distinct(ecoregion, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, sep = "_"))
  } # list of species in each migrant group by ecoregion

lengths(eco_mig_group) # 15 to 119 species per group

# select.axes on subsetted ecoregional migrant group matrices
year_eco_mig_axes <- select.axes(pca, group = eco_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

# 7 dimensions ok!

plot(year_eco_mig_axes) # JK:TG:new<-this is cool- can see that PC3 explains most variance in european atlantic forest residents but not migrants and that PC2 explains most variance in caledon confifer forest migrants but not residents!

summary(year_eco_mig_axes) # explains 95.5 to 98.1% of variance across ecoregions and migrant groups

pdf("PCA_axes_variance_explained_by_ecoregion_status_year_group.pdf", width = 14, height = 10)
plot(year_eco_mig_axes)
dev.off()

# check that dimesnsions are ok for groups in seasonal migrant contribution analysis-----

# join metadata (migratory status, monthly presence)

eco_season_pca <- assemb %>%
  distinct(ecoregion, species, season, status, n_month_year, n_month_season) %>%
  inner_join(avospace, by = "species")

eco_season_mig_group <- eco_season_pca %>%
  distinct(ecoregion, season, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, .$season, sep = "_"))
  } # list of species in each migrant group by ecoregion and season

lengths(eco_season_mig_group) # 15 to 102 species per group

# select.axes on subsetted ecoregional migrant group matrices

season_eco_mig_axes <- select.axes(pca, group = eco_season_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

# 7 dimensions ok!

plot(season_eco_mig_axes)

summary(season_eco_mig_axes) # explains 95.5 to 98.2% of variance

pdf("PCA_axes_variance_explained_by_ecoregion_status_season_group.pdf", width = 14, height = 10)
plot(season_eco_mig_axes)
dev.off()

# JK to JK:TG:new<- can see from this crap plot (too many groups) that actually season does not seem to have much of an effect on the variance of trait space across any ecoregion. BUT for European Atlantic forests, residents are much more spread across PC3 than migrants.


# checking monthly migrant defintion groups for the yearly analysis------------------


eco_month_mig_group <- eco_pca %>%
  distinct(ecoregion, n_month_year, species) %>%
  {
    split(.$species, paste(.$ecoregion, .$n_month_year, sep = "_"))
  } # list of species in each group by ecoregion and number of months present (more granular migratory status definition)

lengths(eco_month_mig_group) # JK:TG:new<- only one ecoregion (Caledon conifer forests) has migrant present in only one month of the year
# 1 to 72 species per group- this why the species-level ellpisoid.contribution metric I made will work better than metrics that work between groups for this monthly analysis

# select.axes on subsetted ecoregional migrant group matrices

# 7 dimensions again fine

eco_month_axes <- select.axes(pca, group = eco_month_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

summary(eco_month_axes) # PC7 explains 95.8 to 98.9% of variance in all groups (apart from migrant_month1 where is NA)

pdf("PCA_axes_variance_explained_by_ecoregion_mon_present_group.pdf", width = 14, height = 10)
plot(eco_month_axes)
dev.off() # JK:TG:new<-cannot plot due to lack of species present in only 1 month of year


# may need to check whether number of dimensions is enough for migrant grouping by number of months present in each season if decide to do this--------------

eco_month_season_mig_group <- eco_season_pca %>%
  distinct(ecoregion, season, n_month_season, species) %>%
  {
    split(.$species, paste(.$ecoregion, .$season, .$n_month_season, sep = "_"))
  } # list of species in each group by ecoregion and number of months present (more granular migratory status definition)

lengths(eco_month_season_mig_group) # JK:TG<- 0 to 124 species, notably low species richness in winter, would have to rarefy but to what?
# again species-level ellipsoid contribution + stats e.g. season*season_mon_present ~ ellipsoid.cont probably better metric for these groups


eco_month_season_axes <- select.axes(pca, group = eco_month_season_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

eco_month_season_axes # 7 dimensions needed

summary(eco_month_season_axes) # PC1:7 explain 96.2 to 99.4% of variance in all groups

pdf("PCOA_axes_variance_explained_by_ecoregion_season_mon_present_group.pdf", width = 14, height = 10)
plot(eco_month_season_axes)
dev.off() # can't plot again as NAs


# variance as a kable (whole space, no groups)---------

pca_variance <- summary(pca)$importance %>%
  t() %>%
  as_tibble(rownames = "PC") %>%
  rename(
    SD = `Standard deviation`, Proportion_of_Variance = `Proportion of Variance`,
    Cumulative_Proportion = `Cumulative Proportion`
  ) # PC1-6 explain 95.1% of variance (>95% threshold); PC1 less dominant than when body mass is not controlled for (26.3% of variance) #TG: classic! Nice one though


kable(pca_variance, align = rep("c", ncol(pca_variance))) # save in table

# save 7 axes with metadata----------------------


avospace_7 <- as_tibble(pca$x, rownames = "species") %>%
  select(species, any_of(paste0("PC", 1:7)))


eco_season_pca_7 <- assemb %>%
  distinct(ecoregion, species, season, status, n_month_year, n_month_season) %>%
  inner_join(avospace, by = "species")

# save full traitspace with metadata (ecoregion, season, mig status, annual months present, seasonal month present)
write.csv(eco_season_pca_7, "eco_region_traitspace_PCA_BI.csv")


#----------------------------perform PCA on beakspace only---------------------------

my_residual_beak_data <- my_residual_data[, c(
  "LogMass", "ResBeak.Length_Culmen",
  "ResBeak.Length_Nares", "ResBeak.Width", "ResBeak.Depth"
)] # select size-relative beak traits


# Check correlations among morphological traits
ggcorrplot(cor(my_residual_beak_data), hc.order = TRUE) # beak depth and width and beak length  traits correlated, body mass not correlated with anything as for full space so allometry removed

# Perform PCA--------------
beak_pca <- prcomp(my_residual_beak_data, scale. = TRUE) # pca with scaling

full_beak_axes <- select.axes(beak_pca, threshold = 0.95, inc.threshold = TRUE) # select axes using disprity
# 4 dimensions needed, explain 98.8% of variance

summary(full_beak_axes)

# save full beak trait space variance explained
pdf("Beak_PCA_axes_variance_explained_full_BI_space", width = 14, height = 10)
plot(full_beak_axes)
dev.off()

# beak avospace as dataframe
beak_avospace <- as_tibble(beak_pca$x, rownames = "species") %>%
  select(species, starts_with(paste0("PC")))


# check that 4 dimensions is ok across all ecoregional migrant subsets for yearly analysis-------------------------

eco_pca_beak <- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(beak_avospace, by = "species")

eco_beak_mig_group <- eco_pca_beak %>%
  distinct(ecoregion, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, sep = "_"))
  } # list of species in each migrant group by ecoregion


# select.axes on subsetted ecoregional migrant group matrices
year_eco_beak_mig_axes <- select.axes(beak_pca, group = eco_beak_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

# 4 dimensions ok!
# JK to JK:TG:new<-  this is cool- can see that PC3 explains most variance in european atlantic forest residents but not migrants and that PC2 explains most variance in caledon confifer forest migrants but not residents!

summary(year_eco_beak_mig_axes) # explains 98.4 to 99.1% of variance across ecoregions and migrant groups

pdf("Beak_PCA_axes_variance_explained_by_ecoregion_status_year_group.pdf", width = 14, height = 10)
plot(year_eco_beak_mig_axes)
dev.off() # PC2 explains much more variance in most ecoregions compared to PC1, especially for Caledon conifer forests and European Atlantic Mixed forests
# again there is a discrepancy between migrants and residents, mainly it seems with the variance explained along PC2.

# #check that 4 dimensions is ok across all ecoregional migrant subsets for seasonal analysis---:

# join metadata (migratory status, monthly presence)

eco_beak_season_pca <- assemb %>%
  distinct(ecoregion, species, season, status, n_month_year, n_month_season) %>%
  inner_join(beak_avospace, by = "species")

eco_beak_season_mig_group <- eco_beak_season_pca %>%
  distinct(ecoregion, season, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, .$season, sep = "_"))
  } # list of species in each migrant group by ecoregion and season


# select.axes on subsetted seasonal ecoregional migrant group matrices
beak_season_eco_mig_axes <- select.axes(beak_pca, group = eco_beak_season_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

# 4 dimensions ok!


summary(beak_season_eco_mig_axes) # PCs1:4 explain 97.9 to 99.3% of variance

pdf("Beak_PCA_axes_variance_explained_by_ecoregion_status_season_group.pdf", width = 14, height = 10)
plot(beak_season_eco_mig_axes)
dev.off()

# JK to JK:TG:new<-  again as for full space, very little difference in variance explained across seasons, but PC2 especially explanatory for European Atlantic mixed forest residents, will be interesting to look at what traits these species represent


# checking monthly migrant defintion groups for the yearly analysis------------------


eco_beak_month_mig_group <- eco_pca_beak %>%
  distinct(ecoregion, n_month_year, species) %>%
  {
    split(.$species, paste(.$ecoregion, .$n_month_year, sep = "_"))
  } # list of species in each group by ecoregion and number of months present (more granular migratory status definition)

lengths(eco_beak_month_mig_group) # only one ecoregion (Caledon conifer forests) has migrant present in only one month of the year

# select.axes on subsetted ecoregional migrant group matrices


beak_eco_month_axes <- select.axes(beak_pca, group = eco_beak_month_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

# this time 5 dimensions would be needed as for English lowland forest migrants present for 2 months of the year explains only 94.9% of variance (bit silly but will save 5 dims)

summary(beak_eco_month_axes) # PC7 explains 95.8 to 98.9% of variance in all groups (apart from migrant_month1 where is NA)

pdf("Beak_PCA_axes_variance_explained_by_ecoregion_mon_present_group.pdf", width = 14, height = 10)
plot(beak_eco_month_axes)
dev.off() # cant save due to same issue as before


beak_avospace_5 <- as_tibble(pca$x, rownames = "species") %>%
  select(species, any_of(paste0("PC", 1:5))) # save 5 dimensions


eco_beak_season_pca_5 <- assemb %>%
  distinct(ecoregion, species, season, status, n_month_year, n_month_season) %>%
  inner_join(beak_avospace_5, by = "species")

# save full traitspace with metadata (ecoregion, season, mig status, annual months present, seasonal month present)
write.csv(eco_beak_season_pca_5, "eco_region_traitspace_PCA_beak_BI.csv")


#------------------Plotting interesting patterns that have spotted with variance plots earlier with plot.dispRity----------------------

# create yearly ecoregionxmigratory status subsets
eco_year_subsets <- custom.subsets(pca$x, group = eco_mig_group)

# Can we spot the Caledon forest migrant spread along PC2?


plot(eco_year_subsets, type = "preview") # PC1 V2- can't see an obvious spread of caledon forest residents along PC2, but hard to tell with colours

plot(eco_year_subsets, type = "preview", specific.args = list(dimensions = c(2, 3))) # maybe in PC2 V3? #slightly, seems like the Caledon resdients group could be skewed by an outlier along PC2, so probably not good to use convex hulls!

# what about the European Atlantic mixed forest migrants?

# yep, quite spread along PC3 but not massively obvious

plot(eco_year_subsets, type = "preview", specific.args = list(dimensions = c(3, 4)))


#------------------------------Checking loadings--------------------------------------------

# full traitspace
full_loadings <- as.data.frame(pca$rotation)[1:7] # loadings for PCs 1-6 #for full loadings, body mass has a small influence on trait space as it only correlates with HWI (which is expected), so traitspace is largely independent of allometry


# @@@JK to TG: does this make sense? I think given this, I will use the full trait space, or I could use the sheffield beak space but I think I may get the new dataset from Gavin later
# TG: it's worth checking how you're going to the loading thingy together. Are you just measuring the correlation between your input variables and the output PC axes? I think in general loadings can be useful for exploring the traitspace and having a vague understanding but in effect, to demonstrate properly it requires much more work. TL:DR; PCA ranks the correlations between your variables hierachically, it doesn't try to understand what variable explains what. E.g. PC5 is the major axes of variance-covariance in your dataset once you've removed the 4 first major axes of variance-covariance.
# JK:TG:new<- yes just looking at correlation between variable and PC axes. Also not sure if i would have to plot loadings for each subset or just the whole trait space, I am assuming just the full space as PC scores are relative to the full space?

# plot loadings for full traitspace---------------------

p1_full <- autoplot(pca,
  x = 1, y = 2, data = avospace_7, colour = "species",
  loadings = TRUE, loadings.colour = "red",
  loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black"
) +
  theme(legend.position = "none")

# PC1 is about wing and tail length
# PC2 is essentially beak length vs beak depth and width- so Caledon residents may have unusual/a wide diversity of beaks?


p2_full <- autoplot(pca,
  x = 3, y = 4, data = avospace_7, colour = "species",
  loadings = TRUE, loadings.colour = "red",
  loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black"
) +
  theme(legend.position = "none")

# PC3 is about tarsus and secondary length vs wing pointiness, so European Atlantic forest migrants likely have a diversity of/unusual wing shape/tarsi traits?
# PC4 is about wing length vs all other traits (but not paritcuarly strong loading)

p3_full <- autoplot(pca,
  x = 5, y = 6, data = avospace_7, colour = "species",
  loadings = TRUE, loadings.colour = "red",
  loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black"
) +
  theme(legend.position = "none")

# PC5 is explained entirely by body mass
# PC6 is wing length, beak depth and tarsus length vs tail length and beak length/width

# combine plots
pall_full <- gridExtra::arrangeGrob(p1_full, p2_full, p3_full, nrow = 1) # body mass only loads onto PC5, which is entirely explained by body mass- this does not have a large loading for any ecoregional migrant group, so size not a big factor in differentiating species/groups in traitspace

pdf("PCA_Loadings_plot_PC1:6.pdf", width = 14, height = 10)
plot(pall_full)
dev.off()

# get PC 7 for interest
p4_full <- autoplot(pca,
  x = 7, y = 6, data = avospace_7, colour = "species",
  loadings = TRUE, loadings.colour = "red",
  loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black"
) +
  theme(legend.position = "none")

p4_full

# JK:TG:new<- PC7 is beak width vs beak depth- it is interesting that beak width seperates from beak depth along PC7, despite the fact that the two traits are heavily correlated?
