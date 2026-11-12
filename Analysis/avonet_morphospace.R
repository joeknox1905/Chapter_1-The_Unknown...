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

# Perform PCA--------------
pca <- prcomp(my_residual_data, scale. = TRUE) # scale

select.axes(pca, threshold = 0.95, inc.threshold = TRUE) # select axes using disprity
# 7 dimensions needed
# TG:done might be good to run this one with your groups. e.g. pca <- custom.subsets(...) so that we can check if 6D is enough for each group.

# create dataframe with tibble for joining metadata
avospace <- as_tibble(pca$x, rownames = "species") %>%
  select(species, starts_with(paste0("PC")))

# check that dimensions are ok for groups in yearly migrant contribution analysis:

# join metadata (migratory status, monthly presence)
eco_pca <- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(avospace, by = "species")

eco_mig_group <- eco_pca %>%
  distinct(ecoregion, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, sep = "_"))
  } # list of species in each migrant group by ecoregion


# select.axes on subsetted ecoregional migrant group matrices

year_eco_mig_axes <- select.axes(pca, group = eco_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

# 7 dimensions ok!

plot(year_eco_mig_axes) # JK:TG<- this is cool- can see that PC3 explains most variance in european atlantic forest residents but not migrants and that PC2 explains most variance in caledon confifer forest migrants but not residents!

summary(year_eco_mig_axes) # explains 95.5 to 98.1% of variance across ecoregions and migrant groups

pdf("PCA_axes_variance_explained_by_ecoregion_status_year_group.pdf", width = 14, height = 10)
plot(year_eco_mig_axes)
dev.off()

# check that dimesnsions are ok for groups in seasonal migrant contribution analysis:

# join metadata (migratory status, monthly presence)

eco_season_pca <- assemb %>%
  distinct(ecoregion, species, season, status, n_month_year, n_month_season) %>%
  inner_join(avospace, by = "species")

eco_season_mig_group <- eco_season_pca %>%
  distinct(ecoregion, season, species, status) %>%
  {
    split(.$species, paste(.$ecoregion, .$status, .$season, sep = "_"))
  } # list of species in each migrant group by ecoregion and season


# select.axes on subsetted ecoregional migrant group matrices

season_eco_mig_axes <- select.axes(pca, group = eco_season_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

# 7 dimensions ok!

plot(season_eco_mig_axes)

summary(season_eco_mig_axes) # explains 95.5 to 98.2% of variance

pdf("PCA_axes_variance_explained_by_ecoregion_status_season_group.pdf", width = 14, height = 10)
plot(season_eco_mig_axes)
dev.off()

# JK:TG<- can see from this crap plot (too many groups) that actually season does not seem to have much of an effect on the variance of trait space across any ecoregion. BUT for European Atlantic forests, residents are much more spread across PC3 than migrants.


# checking monthly migrant defintion groups for the yearly analysis------------------


eco_month_mig_group <- eco_pca %>%
  distinct(ecoregion, n_month_year, species) %>%
  {
    split(.$species, paste(.$ecoregion, .$n_month_year, sep = "_"))
  } # list of species in each group by ecoregion and number of months present (more granular migratory status definition)

lengths(eco_month_mig_group) # only one ecoregion (Caledon conifer forests) has migrant present in only one month of the year

# select.axes on subsetted ecoregional migrant group matrices

# 7 dimensions again fine

eco_month_axes <- select.axes(pca, group = eco_month_mig_group, threshold = 0.95, inc.threshold = TRUE)
# get dimensions necessary for each group

summary(eco_month_axes) # PC7 explains 95.8 to 98.9% of variance in all groups (apart from migrant_month1 where is NA)

pdf("PCA_axes_variance_explained_by_ecoregion_mon_present_group.pdf", width = 14, height = 10)
plot(eco_month_axes)
dev.off() # TG:JK<- cannot plot due to lack of species present in only 1 month of year


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


# 2. AVONET beakspace only (+size)---------

rel_beak_traits <- rel_morph_traits %>%
  select(
    species, LogBodyMass, RelBeakLengthCulmen,
    RelBeakLengthNares, RelBeakWidth, RelBeakDepth
  ) %>%
  column_to_rownames("species")


# Check correlations among morphological traits
ggcorrplot(cor(rel_beak_traits[, 1:5]), hc.order = TRUE) # beak depth and width vs beak length  traits correlated, can group using gawdis if decide to use beak-only space and pcoa

# Perform PCA--------------
beak_pca <- prcomp(rel_beak_traits, scale. = TRUE)

select.axes(beak_pca, threshold = 0.95, inc.threshold = TRUE) # select axes using disprity
# 4 dimensions needed

# explore variance visually
beak_pca_variance <- summary(beak_pca)$importance %>%
  t() %>%
  as_tibble(rownames = "PC") %>%
  rename(
    SD = `Standard deviation`, Proportion_of_Variance = `Proportion of Variance`,
    Cumulative_Proportion = `Cumulative Proportion`
  ) # PC1-4 explain 98.8% of variance (>95% threshold), again body mass is now less dominant, with PC1 explaining 39.1% of variance

kable(beak_pca_variance, align = rep("c", ncol(beak_pca_variance)))


beak_avospace <- as_tibble(beak_pca$x, rownames = "species") %>%
  select(species, all_of(paste0("PC", 1:4)))

eco_beak_pca
# join ecoregional presence to PC scores
eco_beak_pca <- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(beak_avospace, by = "species")


# save beakspace
write.csv(eco_beak_pca, "eco_region_traitspace_PCA_beak_BI.csv")



# Check trait loadings-----------------

# full traitspace
full_loadings <- as.data.frame(pca$rotation)[1:6] # loadings for PCs 1-6 #for full loadings, body mass has a small influence on trait space as it only correlates with HWI (which is expected), so traitspace is largely independent of allometry

# all traits show positive correlation
# PC1: mainly explained by wing length/pointiness
# PC2: mainly explained by wing length/pointiness and beak length and width
# PC3: mainly explained by beak length and shape and tarsus length
# PC4: mainly explained by beak shape
# PC5: mainly explained by body mass
# PC6: mainly explained by beak depth and tarsus length

# beak traitspace
beak_loadings <- as.data.frame(beak_pca$rotation)[1:4] # loadings for PCs 1-4 # for beak loadings body mass has no influence on PC1, 2 or 4 but complete influence on PC3 as there is no HWI included in the trait space, if want to use beakspace will need to rerun, removing size component,
# but this is not ideal as HWI is known to be an important migratory trait.

# @@@JK to TG: does this make sense? I think given this, I will use the full trait space, or I could use the sheffield beak space but I think I may get the new dataset from Gavin later
# TG: it's worth checking how you're going to the loading thingy together. Are you just measuring the correlation between your input variables and the output PC axes? I think in general loadings can be useful for exploring the traitspace and having a vague understanding but in effect, to demonstrate properly it requires much more work. TL:DR; PCA ranks the correlations between your variables hierachically, it doesn't try to understand what variable explains what. E.g. PC5 is the major axes of variance-covariance in your dataset once you've removed the 4 first major axes of variance-covariance.
# yes just looking at correlation between variable and PC axes

# plot loadings for full traitspace---------------------

p1_full <- autoplot(pca,
  x = 1, y = 2, data = avospace, colour = "species",
  loadings = TRUE, loadings.colour = "red",
  loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black"
) +
  theme(legend.position = "none")

# high values of PC1 = long wings and tails
# high values of PC2 = pointy wings and long beaks


p2_full <- autoplot(pca,
  x = 3, y = 4, data = avospace, colour = "species",
  loadings = TRUE, loadings.colour = "red",
  loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black"
) +
  theme(legend.position = "none")

# high values of PC3= large bodied species with short tarsus and beak length
# high values of PC4= wide and deep beaks

p3_full <- autoplot(pca,
  x = 5, y = 6, data = avospace, colour = "species",
  loadings = TRUE, loadings.colour = "red",
  loadings.label = TRUE, loadings.label.size = 3, loadings.label.colour = "black"
) +
  theme(legend.position = "none")

# high values of PC6 = small bodied species
# high values of PC4 = long tails but short tarsus and shallow beak


pall_full <- gridExtra::arrangeGrob(p1_full, p2_full, p3_full, nrow = 1)

ggsave("Manuscript/Figures/pca_loadings_plots.pdf", pall_full, width = 12, height = 6, dpi = 300, bg = "white")
