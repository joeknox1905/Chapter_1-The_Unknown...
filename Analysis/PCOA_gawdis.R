
#load packages-----------
library(dplyr)
library(dispRity)
library(gawdis)
library(tidyverse)
library(ggcorrplot)
library(stringr)
library(kableExtra)



 #This script creates an annual AVONET trait matrix for all 11 traits from species present aoh-refined maps by performing a PCOA with gawdis e.g. (Williams and Curley, 2025; de Bello et al., 2021)
# Then, following Williams and Curley (2025), to control for size, it takes the residuals of the averaged and logged PC traits against
#logged body mass to control for the strong correlation of body mass with size and better distinguish between species with similar body masses but different 
#traits (Weeks et al., 2026). These size-relative traits are used for all except HWI which is a ratio.
#Then, traits are weighted by group (beak shape, beak length, body mass, tarsus length, wing shape, wing/tail length) based on the correlations of the correlation plot and Williams and Curley (2025) to ensure all traits are 
#any one trait is not skewing trait space contributions (Guillerme et al., 2020).
#then, a PCOA is performed using base r. # note NMDS is better for ensuring selected axes capture variation along hidden niche dimensions (Li et al., 2026)
#Then selects.axes is used from dispRity to select the optimal number of axes for explaining the variance of the yearly trait space.


#creating trait matrix based on britain and ireland annual species list (from aoh-refined Birdlife ranges)-------------
assemb<- readRDS("/Data/Processed/seasonal_ecoregion_pools_aoh") #get monthly presence of each species for full annual trait space

#load avonet trait data
avonet<- read.csv("/Data/Raw/AVONET1_BirdLife.csv") #trait means in Birdlife V7
avonet_crosswalk<- read.csv("/Data/Raw/AVONET_Raw_Data.csv") #raw traits for taxonomic crosswalk

#standardise species name format
avonet <- avonet %>% mutate(species = str_replace_all(Species1, " ", "_"))


#get all species from Britain and Ireland aoh maps
spp_name <- sort(unique(assemb$species))

# species in monthly presence data with no trait data from AVONET
no_avonet <- setdiff(spp_name, avonet$species) %>%
  print()

#view these species in avonet crosswalk to manually correct names
no_avonet<- str_replace_all(no_avonet, "_", " ")
avonet_crosswalk %>%
  filter(if_any(everything(), ~ .x%in% no_avonet)) %>%
  View()

#Old Birdlife names to use to get traits from AVONET database, retaining Birdlife V10 names 
#Also get taxonomy match from Rob in future
avonet_lookup <- c(
  Curruca_communis = "Sylvia_communis",
  Curruca_curruca   = "Sylvia_curruca",
  Curruca_undata    = "Sylvia_undata",
  Luscinia_svecica = "Cyanecula_svecica")

#join all eleven traits to all species in British Isles,
avonet_traits <- tibble(species = spp_name,) %>%
  mutate(lookup_name = coalesce(avonet_lookup[species], species)) %>%
  left_join(avonet %>% select(species, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
                              Tarsus.Length, Wing.Length, Secondary1, Tail.Length, Mass, Kipps.Distance, Hand.Wing.Index), by = c("lookup_name" = "species")) %>%
  dplyr::select(-lookup_name)

#get species with no trait data in avonet (should be 0)
avonet_na_any <- avonet_traits %>% filter(if_any(-species, is.na)) %>%
  print()
avonet_na_all <- avonet_traits %>% filter(if_all(-species, is.na)) %>%
  print()

#save trait matrix yearly
write.csv(avonet_traits, "avonet_traits_BI_avg_raw_yearly.csv")


#control for size in morphological traits (Williams and CUrley, 2025 residuals method)------------------

# AVONET morphological traits
morph_traits<- avonet_traits %>%
select(species, Beak.Length_Culmen, Beak.Length_Nares, Beak.Width, Beak.Depth,
Tarsus.Length, Wing.Length, Secondary1, Tail.Length, Mass, Kipps.Distance, Hand.Wing.Index)


# Log absolute traits (excluding HWI which is a ratio? or excluding completely?)
morph_traits<- morph_traits %>%
  mutate(LogBodyMass = log(Mass),
         LogBeakLengthCulmen = log(Beak.Length_Culmen),
         LogBeakLengthNares = log(Beak.Length_Nares),
         LogBeakWidth = log(Beak.Width),
         LogBeakDepth = log(Beak.Depth),
         LogTarsusLength = log(Tarsus.Length),
         LogKippsDistance = log(Kipps.Distance),
         LogWingLength = log(Wing.Length),
         LogSecondary1 = log(Secondary1),
         LogTailLength = log(Tail.Length))

# Calculate morphological traits relative to body mass, not HWI which is a ratio

rel_morph_traits <- morph_traits%>%
  mutate(RelBeakLengthCulmen = lm(LogBeakLengthCulmen ~ LogBodyMass,
                                  data=morph_traits)$residuals,
         RelBeakLengthNares = lm(LogBeakLengthNares ~ LogBodyMass,
                                 data=morph_traits)$residuals,
         RelBeakWidth = lm(LogBeakWidth ~ LogBodyMass,
                           data=morph_traits)$residuals,
         RelBeakDepth = lm(LogBeakDepth ~ LogBodyMass,
                           data=morph_traits)$residuals,
         RelTarsusLength = lm(LogTarsusLength ~ LogBodyMass,
                              data=morph_traits)$residuals,
         RelKippsDistance = lm(LogKippsDistance ~ LogBodyMass,
                               data=morph_traits)$residuals,
         RelWingLength = lm(LogWingLength ~ LogBodyMass,
                            data=morph_traits)$residuals,
         RelSecondary1 = lm(LogSecondary1 ~ LogBodyMass,
                           data=morph_traits)$residuals,
         RelTailLength = lm(LogTailLength ~ LogBodyMass,
                            data=morph_traits)$residuals) %>%
  select(species, LogBodyMass, RelBeakLengthCulmen,
         RelBeakLengthNares, RelBeakWidth, RelBeakDepth,
         RelTarsusLength, RelKippsDistance, Hand.Wing.Index,
         RelWingLength, RelSecondary1, RelTailLength) 



# Check correlations among morphological traits
ggcorrplot(cor(rel_morph_traits[,2:12]), hc.order = TRUE) # correlated variables ordered together
#get groups for gawdis as allows for grouping of traits (i.e. HWI correlated with body mass, beak length traits correlated, wing traits correlated, wing and tail traits correlated)

rel_morph_traits<- rel_morph_traits %>%
    column_to_rownames("species")
    

# Create distance matrix of species using functional traits 
dist_mat <- gawdis(rel_morph_traits, w.type ="analytic", #weights optimised by mathematical algorithm
                       groups=c(1, # Body mass # also correlates with HWI but makes sense as larger species tend to have a smaller HWI (as well as sedentary species)
                                2,2, # Beak length
                                3,3, # Beak shape (beak depth and width)
                                4, # Tarsus length
                                5,5, # Wing shape (HWI, kipps distance) # hard to differentiate groups here as kipps distance also correlated with wing length but chose based on logic of Williams and Curley 2025 as both measure pointiness of wing. Doesn't change weights much anyway.
                                6,6, 6), # Wing/tail length (wing length, tail length and secondary 1)
                       groups.weight=TRUE) # traits will have different weights within groups as well as between groups 
attr(dist_mat, "weights") # get weights of each group

#fairly even, beak shape given most weight (0.22), followed by body mass, beak length, wing length, tarsus length and wing/tail length (0.15)


#Perform PCOA-----------

#pcoa on gawdis distance matrix using base r, add = TRUE adds a calliez correction (computed constant as 0.45) to account for negative eigen values and ensure distances are euclidean (Li et al., 2026)
pcoa <- cmdscale(dist_mat, eig = TRUE, add = TRUE, k = 12) 
select.axes(pcoa$points, threshold = 0.95, inc.threshold = TRUE) #10 dimensions are needed?

#save with ecoregional data

rownames(pcoa$points) <- labels(dist_mat)

pcoa_scores <- pcoa$points[, 1:10] %>%
  as.data.frame() %>%
  rownames_to_column("species")

#join eco-regional presence and status to pc scores

eco_pcoa<- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(pcoa_scores, by = "species")


write_csv(eco_pcoa, "eco_region_traitspace__PCOA_BI.csv")


#alternatively could do pca on log transformed and scaled traits---------


#scale (Mammola et al., 2021; Keyser et al., 2024) traits and perform PCA

pca <- prcomp(rel_morph_traits, scale. = TRUE)

select.axes(pca, threshold = 0.95, inc.threshold = TRUE) # select axes using disprity
# 6 dimensions needed

rel_morph_traits<-  rel_morph_traits %>%
  as.data.frame() %>%
  rownames_to_column("species")


avospace <- as_tibble(pca$x) %>%
mutate(species = rel_morph_traits$species, .before = 1) %>%
select(species, all_of(paste0("PC", 1:6))) #95.1% of variance explained


#explore variance visually---------

pca_variance <- summary(pca)$importance %>%
  t() %>%
  as_tibble(rownames = "PC") %>%
  rename(SD = `Standard deviation`, Proportion_of_Variance = `Proportion of Variance`,
         Cumulative_Proportion = `Cumulative Proportion`) #PC1-4 explain 97.1% of variance (>95% threshold) 

kable(pca_variance, align = rep("c", ncol(pca_variance)))


eco_pca<- assemb %>%
distinct(ecoregion, species, status, n_month_year) %>%
inner_join(avospace, by = "species")


write.csv(eco_pca, "eco_region_traitspace__PCA_BI.csv")






  







