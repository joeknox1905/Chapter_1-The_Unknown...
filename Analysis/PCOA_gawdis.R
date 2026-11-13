# load packages-----------
library(dplyr)
library(dispRity)
library(gawdis)
library(stats)
library(tidyverse)
library(ggcorrplot)
library(stringr)
library(kableExtra)
library(vegan)
library(ape)


# This script creates an annual AVONET trait matrix for all 11 traits from species present aoh-refined maps by performing a PCOA with gawdis e.g. (Williams and Curley, 2025; de Bello et al., 2021)
# Then, following Williams and Curley (2025), to control for size, it takes the residuals of the averaged and logged PC traits against
# logged body mass to control for the strong correlation of body mass with size and better distinguish between species with similar body masses but different
# traits (Weeks et al., 2026). These size-relative traits are used for all except HWI which is a ratio.
# Then, traits are weighted by group (beak shape, beak length, body mass, tarsus length, wing shape, wing/tail length) based on the correlations of the correlation plot and Williams and Curley (2025) to ensure all traits are
# any one trait is not skewing trait space contributions (Guillerme et al., 2020).
# then, a PCOA is performed using base r. # note NMDS is better for ensuring selected axes capture variation along hidden niche dimensions (Li et al., 2026)
# Then selects.axes is used from dispRity to select the optimal number of axes for explaining the variance of the yearly trait space.


# creating trait matrix based on britain and ireland annual species list (from aoh-refined Birdlife ranges)-------------
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

# join all ten traits to all species in British Isles, excluding HWI which is a ratio of kipps distance/wing length x100
avonet_traits <- tibble(species = spp_name, ) %>%
  mutate(lookup_name = coalesce(avonet_lookup[species], species)) %>%
  left_join(avonet %>% select(
    species, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
    Tarsus.Length, Wing.Length, Secondary1, Tail.Length, Mass, Kipps.Distance
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

# 10 AVONET morphological traits (excluding HWI as this = Kipps.Distance/Wing.Length x100; so would just describe variables already considered)
morph_traits <- avonet_traits %>%
  select(species, Mass, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth, Tarsus.Length, Wing.Length, Kipps.Distance, Secondary1, Tail.Length) %>%
  column_to_rownames("species")


# TG:done In general I would avoid ratios, especially if they are based on some variables already present. E.g. if HWI = BeakWidth/BeakDepth you're effectively counting these two variables twice compared to the other ones not in a ratio. Also, logging is probably not a problem here but it should not be an automatic thing (e.g. why do you need to log the data? If the answer is "because that's just what people do", you probably should not log the data).

# Log body mass for function call
LogMass <- log(morph_traits$Mass)

# Log and transform morphological traits relative to body mass (residuals)
get.residuals <- function(x) {
  LogTrait <- log(x)
  ResTrait <- lm(LogTrait ~ LogMass)$residuals
}

res_traits <- apply(morph_traits[, -1], 2, get.residuals) # apply to all columns except body mass
colnames(res_traits) <- paste0("Res", colnames(res_traits)) # rename to reflective relative trait values
my_residual_data <- cbind(LogMass, res_traits)


# Check correlations among morphological traits
ggcorrplot(cor(my_residual_data[, 1:10]), hc.order = TRUE) # correlated variables ordered together
# get groups for gawdis as allows for grouping of traits (i.e.body mass not correlayed with anything beak length traits correlated, wing length and pointiness correlated, wing width and tail length correlated)


# scale traits JK:TG:new<- not sure if I need to scale if scaling by weights with gawdis below?


# Create distance matrix of species using functional traits (WIlliams and Curley, 2025) based on trait correlations
dist_mat <- gawdis(my_residual_data,
  w.type = "analytic", # weights optimised by mathematical algorithm
  groups = c(
    1, # Body mass
    2, 2, # Beak length traits
    3, 3, # Beak shape (beak depth and width)
    4, # Tarsus length
    5, 5, # Wing length and pointiness (kipps distance, wing length) # hard to differentiate groups here as kipps distance also correlated with wing length but chose based on logic of Williams and Curley 2025 as both measure pointiness of wing. Doesn't change weights much anyway.
    6, 6 # Wing width (first secondary length) and tail length
  ),
  groups.weight = TRUE
) # traits will have different weights within groups as well as between groups
attr(dist_mat, "weights") # get weights of each group

# gawdis scales each trait by dividing by the range of species' dissimilarity values for each trait

# fairly even weights, beak shape given most weight (0.22), followed by body mass, wing length/pointiness, beak length, tarsus length and wing/tail length (0.15)

# join species names to distance matrix
attr(dist_mat, "Labels") <- rownames(my_residual_data)

# Perform PCOA--------------

# pcoa on gawdis distance matrix using base r,
# @@@ JK to TG: hopefully this is the right approach; 23 dimensions seems like quite a lot, 175 is way too many no?
# TG: we can walk through the whole procedure on our meeting on Friday. It's not clear to me why you go through all these steps here (residuals and shit). Not saying it's the wrong way but as my comment above for the logs, it needs to be justified for you and done in consciousness ;). But otherwise, no, 23 dimensions is totally fine. With PCoA, you tend to have $number_of_elements - 2$ dimensions because the distance step tend to spread the variance more evenly.

pcoa1 <- wcmdscale(d = dist_mat, eig = TRUE)

pcoa1$eig
# all negative eigenvalues outside of first 3 dimensions, so adding Calliez correction to ensure distances are euclidean e.g (Li et al., 2026)

pcoa2 <- wcmdscale(d = dist_mat, eig = TRUE, add = TRUE)

pcoa2$eig # eigen values are now all positive


# using ape
pcoa3 <- ape::pcoa(dist_mat) # without correction


select.axes(pcoa1$points, threshold = 0.95, inc.threshold = TRUE) # 23 dimensions
select.axes(pcoa2$points, threshold = 0.95, inc.threshold = TRUE) # 176 dimensions? WHAT??
full_space_axes <- select.axes(pcoa3$vectors, threshold = 0.95, inc.threshold = TRUE)
full_space_axes # 23 dimensions are needed, 'vectors' = extracts all positive eigen values so don't need correction

# Both ape::pcoa and vegan::wmdscale() perform similarly so using ape-------------


summary(full_space_axes) # Axes 1:23 explain 95.3% of variance

pdf("PCOA_axes_variance_explained_full_BI_space", width = 14, height = 10)
plot(full_space_axes)
dev.off() # really seems that axes 1:5 or at least 1:10 explain most variance


# create dataframe with tibble for joining metadata
avospace <- as_tibble(pcoa3$vectors, rownames = "species") %>%
  select(species, starts_with(paste0("Axis")))

# check that dimensions are ok for groups in yearly migrant contribution analysis-------

# join metadata (migratory status, monthly presence)
eco_pcoa <- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(avospace, by = "species")

eco_mig_group <- eco_pcoa %>%
  distinct(ecoregion, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, sep = "_"))
  } # list of species in each migrant group by ecoregion

lengths(eco_mig_group) # 15 to 119 species per group

# select.axes on subsetted ecoregional migrant group matrices
year_eco_mig_axes <- select.axes(pcoa3$vectors, group = eco_mig_group, threshold = 0.95, inc.threshold = TRUE)
year_eco_mig_axes # 24 dimensions needed
# get dimensions necessary for each group


summary(year_eco_mig_axes) # explains 95.1 to 95.9% of variance across ecoregions and migrant groups

pdf("PCOA_axes_variance_explained_by_ecoregion_status_year_group.pdf", width = 14, height = 10)
plot(year_eco_mig_axes)
dev.off() # PC1 is dominant across all ecoregions and migrant groups unlike PCA. Noticeably, PC2 is much less dominant in european atlantic forest residents (and caledon forest residents), whilst PC3 is more dominant in Caledon forest residents than migrants

# check that dimesnsions are ok for groups in seasonal migrant contribution analysis-----

# join metadata (migratory status, monthly presence)

eco_season_pcoa <- assemb %>%
  distinct(ecoregion, species, season, status, n_month_year, n_month_season) %>%
  inner_join(avospace, by = "species")

eco_season_mig_group <- eco_season_pcoa %>%
  distinct(ecoregion, season, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, .$season, sep = "_"))
  } # list of species in each migrant group by ecoregion and season

lengths(eco_season_mig_group) # 15 to 102 species per group, notably species richness is quite similar across seasonal migrant groups

# select.axes on subsetted ecoregional migrant group matrices

season_eco_mig_axes <- select.axes(pcoa3$vectors, group = eco_season_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

season_eco_mig_axes # 24 dimensions needed,

summary(season_eco_mig_axes) # explains 95.1 to 96% of variance

pdf("PCOA_axes_variance_explained_by_ecoregion_status_season_group.pdf", width = 14, height = 10)
plot(season_eco_mig_axes)
dev.off()

# JK to Tg:<- again as with PCA, most variation in variance explained appears to come from differences between migrants and residents rather than seasonal differences i.e. residents have very similar variance patterns across seasons and migrants have very similar variance patterns across seasons but migrants and residents differ.
# I guess, seems that the PCOA + gawdis weighting vs PCA procedure doesn't make too much of a difference.


# checking monthly migrant defintion groups for the yearly analysis------------------


eco_month_mig_group <- eco_pcoa %>%
  distinct(ecoregion, n_month_year, species) %>%
  {
    split(.$species, paste(.$ecoregion, .$n_month_year, sep = "_"))
  } # list of species in each group by ecoregion and number of months present (more granular migratory status definition)

lengths(eco_month_mig_group) # JK:TG:new<- only one ecoregion (Caledon conifer forests) has migrant present in only one month of the year
# 1 to 72 species per group- this why the species-level ellpisoid.contribution metric I made will work better than metrics that work between groups for this monthly analysis

# select.axes on subsetted ecoregional migrant group matrices


eco_month_axes <- select.axes(pcoa3$vectors, group = eco_month_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

eco_month_axes # 36 dimensions needed

summary(eco_month_axes) # PC1:36 explain 96 to 98.9% of variance in all groups (apart from migrant_month1 where is NA)

pdf("PCOA_axes_variance_explained_by_ecoregion_mon_present_group.pdf", width = 14, height = 10)
plot(eco_month_axes)
dev.off() # JK:TG:new<-cannot plot due to lack of species present in only 1 month of year

# may need to check whether number of dimensions is enough for migrant grouping by number of months present in each season if decide to do this

eco_month_season_mig_group <- eco_season_pcoa %>%
  distinct(ecoregion, season, n_month_season, species) %>%
  {
    split(.$species, paste(.$ecoregion, .$season, .$n_month_season, sep = "_"))
  } # list of species in each group by ecoregion and number of months present (more granular migratory status definition)

lengths(eco_month_season_mig_group) # JK:TG<- 0 to 124 species, notably low species richness in winter, would have to rarefy but to what?
# again species-level ellipsoid contribution + stats e.g. season*season_mon_present ~ ellipsoid.cont probably better metric for these groups


eco_month_season_axes <- select.axes(pcoa3$vectors, group = eco_month_season_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

eco_month_season_axes # 35 dimensions needed

summary(eco_month_season_axes) # PC1:35 explain 96.4 to 99% of variance in all groups

pdf("PCOA_axes_variance_explained_by_ecoregion_season_mon_present_group.pdf", width = 14, height = 10)
plot(eco_month_season_axes)
dev.off() # can't plot again


# save trait space with required dimension number------------------

pcoa_axes <- pcoa3$vectors[, 1:36] # save pcoa points for axes with positive eigen values


pcoa_axes_df <- pcoa_axes %>%
  as.data.frame() %>%
  rownames_to_column("species")

# save with ecoregional and migratory status data
eco_pcoa <- assemb %>%
  distinct(species, ecoregion, season, status, n_month_year, n_month_season) %>%
  inner_join(pcoa_axes_df, by = "species")

write_csv(eco_pcoa, "eco_region_traitspace_PCOA_BI.csv")


#-----------------get "loadings"---------------


# getting variance explained by eigen values?
eig_percent <- round(pcoa3$values$Relative_eig * 100, 2) # @@@JK to TG: variance explains adds to more than 100? # TG: hehehe, that's probably not what it should be doing ;) but we can see where it's going wrong. Maybe because of some correction steps above?


# plot with disprity---------------

# create yearly ecoregion x migratory status subsets
eco_year_subsets <- custom.subsets(pcoa3$vectors, group = eco_mig_group)

plot.dispRity(pcoa_axes, type = "preview") # PC1 v 2
plot.dispRity(pcoa_axes, type = "preview", specific.args = list(dimensions = c(2, 4))) # PC2 V 4 etc
