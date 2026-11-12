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

# join all eleven traits to all species in British Isles,
avonet_traits <- tibble(species = spp_name, ) %>%
  mutate(lookup_name = coalesce(avonet_lookup[species], species)) %>%
  left_join(avonet %>% select(
    species, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
    Tarsus.Length, Wing.Length, Secondary1, Tail.Length, Mass, Kipps.Distance, Hand.Wing.Index
  ), by = c("lookup_name" = "species")) %>%
  dplyr::select(-lookup_name)

# get species with no trait data in avonet (should be 0)
avonet_na_any <- avonet_traits %>%
  filter(if_any(-species, is.na)) %>%
  print()
avonet_na_all <- avonet_traits %>%
  filter(if_all(-species, is.na)) %>%
  print()

# save trait matrix yearly
write.csv(avonet_traits, "avonet_traits_BI_avg_raw_yearly.csv")


# control for size in morphological traits (Williams and CUrley, 2025 residuals method)------------------

# 10 AVONET morphological traits (excluding HWI as this = Kipps.Distance/Wing.Length x100; so would just describe variables already considered)
morph_traits <- avonet_traits %>%
  select(
    species, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
    Tarsus.Length, Wing.Length, Secondary1, Tail.Length, Mass, Kipps.Distance
  )


# TG:done In general I would avoid ratios, especially if they are based on some variables already present. E.g. if HWI = BeakWidth/BeakDepth you're effectively counting these two variables twice compared to the other ones not in a ratio. Also, logging is probably not a problem here but it should not be an automatic thing (e.g. why do you need to log the data? If the answer is "because that's just what people do", you probably should not log the data).

# Log body mass for function call
LogMass <- log(morph_traits$Mass)

# Log and transform morphological traits relative to body mass (residuals)
get.residuals <- function(x) {
  LogTrait <- log(x)
  ResTrait<- lm(LogTrait ~ LogMass)$residuals
}

res_traits <- apply(morph_traits[, -9], 2, get.residuals) # apply to all columns except body mass
colnames(res_traits) <- paste0("Res", colnames(res_traits))  #rename to reflective relative trait values
my_residual_data <- cbind(res_traits, LogMass)


# Check correlations among morphological traits
ggcorrplot(cor(my_residual_data[, 1:10]), hc.order = TRUE) # correlated variables ordered together
# get groups for gawdis as allows for grouping of traits (i.e. HWI correlated with body mass, beak length traits correlated, wing traits correlated, wing and tail traits correlated)

my_residual_data <- rel_morph_traits %>%
  select(1:11) %>%
  scale() # scale traits


# Create distance matrix of species using functional traits (WIlliams and Curley, 2025) based on trait correlations
dist_mat <- gawdis(my_residual_data,
  w.type = "analytic", # weights optimised by mathematical algorithm
  groups = c(
    1, # Body mass 
    2, 2, # Beak length
    3, 3, # Beak shape (beak depth and width)
    4, # Tarsus length
    5, 5, # Wing length and pointiness (kipps distance, wing length) # hard to differentiate groups here as kipps distance also correlated with wing length but chose based on logic of Williams and Curley 2025 as both measure pointiness of wing. Doesn't change weights much anyway.
    6, 6 # Wing width (first secondary length) and tail length
  ), 
  groups.weight = TRUE
) # traits will have different weights within groups as well as between groups
attr(dist_mat, "weights") # get weights of each group

# scales each trait by dividing by the range of species' dissimilarity values for each trait

# fairly even, beak shape given most weight (0.22), followed by body mass, beak length, wing length, tarsus length and wing/tail length (0.15)


# Perform PCOA-----------

# pcoa on gawdis distance matrix using base r,
# @@@ JK to TG: hopefully this is the right approach; 23 dimensions seems like quite a lot, 175 is way too many no?
# TG: we can walk through the whole procedure on our meeting on Friday. It's not clear to me why you go through all these steps here (residuals and shit). Not saying it's the wrong way but as my comment above for the logs, it needs to be justified for you and done in consciousness ;). But otherwise, no, 23 dimensions is totally fine. With PCoA, you tend to have $number_of_elements - 2$ dimensions because the distance step tend to spread the variance more evenly.

pcoa1 <- wcmdscale(d = dist_mat, eig = TRUE)

pcoa1$eig
# all negative eigenvalues outside of first 2 dimensions, so adding Calliez correction to ensure distances are euclidean e.g (Li et al., 2026)

pcoa2 <- wcmdscale(d = dist_mat, eig = TRUE, add = TRUE)

pcoa2$eig


# using ape
pcoa3 <- ape::pcoa(dist_mat) # without correction


select.axes(pcoa1$points, threshold = 0.95, inc.threshold = TRUE) # 23 dimensions?-
select.axes(pcoa2$points, threshold = 0.95, inc.threshold = TRUE) # 175 dimensions? WHAT??
select.axes(pcoa3$vectors, threshold = 0.95, inc.threshold = TRUE) # 23 dimensions are needed? #note that this gets the axes with positive eigenvalues


# so both ape::pcoa and vegan::wmdscale() perform similarly

# using ape---------------


# create dataframe with tibble for joining metadata
avospace <- as_tibble(pcoa3$vectors, rownames = "species") %>%
  select(species, all_of(paste0("Axis"))) 

# check that dimensions are ok for groups in yearly migrant contribution analysis:

# join metadata (migratory status, monthly presence)
eco_pcoa <- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(avospace, by = "species")

eco_mig_group <- eco_pcoa %>%
  distinct(ecoregion, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, sep = "_"))
  } # create migrant groups by species, ecoregion and status (1 row)


eco_pcoa_groups <- custom.subsets(pcoa3$vectors, group = eco_mig_group) # split pca by ecoregion and status to keep relative to full BI traitspace
select.axes(eco_pcoa_groups, threshold = 0.95, inc.threshold = TRUE) # get dimensions necessary for each group


# check that dimesnsions are ok for groups in seasonal migrant contribution analysis:

# join metadata (migratory status, monthly presence)
eco_pcoa_season <- assemb %>%avospace <- as_tibble(pcoa3$vectors, rownames = "species") %>%
  select(species, all_of(paste0("Axis"))) 
  distinct(ecoregion, season, species, status, n_month_year, n_month_season) %>%
  inner_join(avospace, by = "species")

eco_season_mig_group <- eco_pcoa_season %>%
  distinct(ecoregion, season, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, sep = "_"))
  } # create migrant groups by species, ecoregion and status (1 row)


eco_season_pca_groups <- custom.subsets(pcoa3$vectors, group = eco_season_mig_group) # split pca by ecoregion and status to keep relative to full BI traitspace
select.axes(eco_season_pca_groups, threshold = 0.95, inc.threshold = TRUE) # get dimensions necessary for each group


pcoa_axes <- pcoa3$vectors[, 1:23] # save pcoa points for axes with positive eigen values



pcoa_axes_df <- pcoa_axes %>%
  as.data.frame() %>%
  rownames_to_column("species")

# save with ecoregional and migratory status data
eco_pcoa <- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(pcoa_axes, by = "species")

write_csv(eco_pcoa, "eco_region_traitspace_PCOA_BI.csv")


# plot PCOA pc1 v 2 and compare to pca-----------------


# getting variance explained by eigen values
eig_percent <- round(pcoa3$values$Relative_eig * 100, 2) # @@@JK to TG: variance explains adds to more than 100? # TG: hehehe, that's probably not what it should be doing ;) but we can see where it's going wrong. Maybe because of some correction steps above?


# TG:done dispRi-tip if you plot them through dispRity you get the dimensions stuff automatically ( plot.dispRity(pcoa_axes, type = "preview") ).



# plot with disprity---------------
 

plot.dispRity(pcoa_axes, type = "preview") # PC1 v 2
plot.dispRity(pcoa_axes, type = "preview", specific.args = list(dimensions = c(2, 4))) # PC2 V 4 etc


#get specific subsets--------------

mig_mat <- eco_pcoa  %>%
distinct(species) %>%
select(paste0("PC"[]) %in% status = "migrant") %>%
column_to_rownames("species") %>%
as.matrix() 

res_mat <- eco_pcoa  %>%
distinct(species) %>%
select(paste0("PC"[]) %in% status = "resident") %>%
column_to_rownames("species") %>%
as.matrix() 

#check that dimensions are enough per group


select.axes(pcoa3$vectors<- custom.subsets(groups = c(mig_res, mig_mat), threshold = 0.95, inc.threshold = TRUE))

#would i need to select axes per ecoregion?

plot.dispRity(pcoa_axes, type = "preview", specific.args = list(dimensions = c(2, 4), matrix= c(mig_mat, res_mat))) #? for subsets