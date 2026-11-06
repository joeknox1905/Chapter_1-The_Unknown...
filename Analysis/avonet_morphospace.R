#This script first creates a full "morphospace" for the 195 species present in britain and ireland using aoh-refined Birdlife range maps based on 11 morphological traits in AVONET
#HWI, body mass, beak depth, beak width, beak length from the culmen, beak length from the nares, wing length, length to the secondary flight feather, tail length and tarsus length (Keyser et al., 2024)
#beak traits = trophic traits linked to size and type of food consumed, tarsus length, wing traits
#and tail length = locomotory traits linked to foraging mode, dispersal and microhabitat use, body mass is linked to metabolism and foraging behaviour (Stewart et al., 2022). 
#HWI may also be useful for calculating dispersal buffers later (Stewart et al., 2022).
# Then, following Williams and Curley (2025), to control for size, it takes the residuals of the averaged and logged PC traits against
#logged body mass to control for the strong correlation of body mass with size and better distinguish between species with similar body masses but different 
#traits (Weeks et al., 2026). These size-relative traits are used for all except HWI which is a ratio.
#It a beak only traitspace and a non-beak only trait space to assess which traits are important in characterising the avian trait space and the migrant vs resident traitspace.
#These comparisons are repeated excluding HWI and kipps distance (as in Keyser et al., 2022) as likely highly correlated and shown to skew PC1 and 2 comapred to non-beak traitspace.
#see chapter 1 script)



#-----load libraries------

library(tidyr)
library(stringr)
library(kableExtra)
library(dispRity)
library(stats)
library(gridExtra)
library(dplyr)

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
ggcorrplot(cor(rel_morph_traits[,2:12]), hc.order = TRUE) 

#Perform PCA--------------
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

#save full traitspace
write.csv(eco_pca, "eco_region_traitspace_PCA_BI.csv")



#2. AVONET beakspace only (+size)---------

rel_beak_traits<- rel_morph_traits %>%
select(species, LogBodyMass, RelBeakLengthCulmen,
  RelBeakLengthNares, RelBeakWidth, RelBeakDepth) %>%
  column_to_rownames("species")


# Check correlations among morphological traits
ggcorrplot(cor(rel_beak_traits[,2:6]), hc.order = TRUE) #beak depth and width vs beak length  traits correlated, can group using gawdis if decide to use beak-only space and pcoa

#Perform PCA--------------
beak_pca <- prcomp(rel_beak_traits, scale. = TRUE)

select.axes(beak_pca, threshold = 0.95, inc.threshold = TRUE) # select axes using disprity
# 4 dimensions needed

#explore variance visually
beak_pca_variance <- summary(beak_pca)$importance %>%
  t() %>%
  as_tibble(rownames = "PC") %>%
  rename(SD = `Standard deviation`, Proportion_of_Variance = `Proportion of Variance`,
         Cumulative_Proportion = `Cumulative Proportion`) #PC1-4 explain 98.8% of variance (>95% threshold) 

kable(beak_pca_variance, align = rep("c", ncol(beak_pca_variance)))

rel_beak_traits<-  rel_beak_traits %>%
  as.data.frame() %>%
  rownames_to_column("species")

beak_avospace <- as_tibble(beak_pca$x) %>%
  mutate(species = rel_morph_traits$species, .before = 1) %>%
  select(species, all_of(paste0("PC", 1:4))) 


#join ecoregional presence to PC scores
eco_beak_pca<- assemb %>%
  distinct(ecoregion, species, status, n_month_year) %>%
  inner_join(beak_avospace, by = "species")


#save beakspace
write.csv(eco_beak_pca, "eco_region_traitspace_PCA_beak_BI.csv")





  

