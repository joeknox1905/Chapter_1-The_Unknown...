#This script creates a "beakspace" for 264 birds present in the monthly British Isles dataset
#PC scores were extracted from an 8-dimensional hypervolume from (Guillerme et al., 2023) derived from 4 landmarks and 3 curves of semi landmarks from 3D scans of 8748 bird species aligned with Procrustes superimposition (see Cooney et al., 2017; Chira et al., 2020; Hughes et al., 2022). 
#Size axes were removed as may obscure shape components which are more indicative of ecological function (Chira et al., 2018)) and symmetry was controlled for (Guillerme et al., 2023). Dimensions explained 98.7% variance among species. 

#----load libraries---------

library(dplyr)
library(tidyr)
library(stringr)
library(kableExtra)



#read in data------------

#load monthly presence data by species, site and season
mon_presence<- readRDS("/Data/Processed/seasonal_ecoregion_pools_aoh") #shapespace

#load sheffield shapespace (no size; procustes on scaled landmarks ) for 8748 species
get_rda<- load("/Data/Raw/shapespace_allbirds_lvl_superorder_order.rda")

beakspace <- get(get_rda)
beakspace <- beakspace[[1]]
sheffield_species <- beakspace$space$animal

#load sheffield formspace (plus size; procustes to align landmarks but no scaling) for 8748 species
formspace<- readRDS("/Data/Raw//2020_08_07_MMB_MORPHO_FORMSPACE_FULL.rds")
form_scores <- formspace$PCscores #extract pc scores
rownames(form_scores) <- names(formspace$size) #get species names
form_scores <- as_tibble(form_scores, rownames = "species") #create species column


#---correct taxonomic mismatches---------

spp_name <- sort(unique(mon_presence$species))

#get species in monthly presence data not in Sheffield beakspace dataset
no_beak <- setdiff(spp_name, sheffield_species) %>%
  print()

#view rows from avonet crosswalk (BirdTree) to manually correct names
no_beak<- str_replace_all(no_beak, "_", " ")
avonet_crosswalk %>%
  filter(if_any(everything(), ~ .x%in% no_beak)) %>%
  View()

#BirdTree names to use to get traits from beak landmark dataset, retaining Birdlife V7 names from monthly presence data
beak_lookup <- c(
  Acanthis_flammea = "Carduelis_flammea", Alexandrinus_krameri = "Psittacula_krameri",
  Anser_caerulescens = "Chen_caerulescens", Ardea_alba = "Casmerodius_albus",
  Ardenna_grisea = "Puffinus_griseus", Calidris_falcinellus = "Limicola_falcinellus",
  Calidris_pugnax = "Philomachus_pugnax", Calonectris_borealis = "Calonectris_diomedea",
  Chloris_chloris = "Carduelis_chloris", Curruca_communis = "Sylvia_communis",
  Curruca_curruca = "Sylvia_curruca", Curruca_undata = "Sylvia_undata",
  Cyanistes_caeruleus = "Parus_caeruleus", Dryobates_minor = "Dendrocopos_minor",
  Emberiza_calandra = "Miliaria_calandra", Gulosus_aristotelis = "Phalacrocorax_aristotelis",
  Hydrobates_leucorhous = "Oceanodroma_leucorhoa", Hydrocoloeus_minutus = "Larus_minutus",
  Linaria_cannabina = "Carduelis_cannabina", Linaria_flavirostris = "Carduelis_flavirostris",
  Lophophanes_cristatus = "Parus_cristatus", Lyrurus_tetrix = "Tetrao_tetrix",
  Mareca_penelope = "Anas_penelope", Mareca_strepera = "Anas_strepera",
  Periparus_ater = "Parus_ater", Poecile_montanus = "Parus_montanus",
  Poecile_palustris = "Parus_palustris", Pterodroma_deserta = "Pterodroma_feae",
  Spatula_clypeata = "Anas_clypeata", Spatula_querquedula = "Anas_querquedula",
  Spinus_spinus = "Carduelis_spinus", Sternula_albifrons = "Sterna_albifrons",
  Thalasseus_sandvicensis = "Sterna_sandvicensis")

#Pterodroma_cahow left out as no data in beakspace

#------join ordinated trait data to species in monthly presence---------

#join 8 PC scores to species in British Isles to shapespace
shape_traits <- tibble(species = spp_name) %>%
  mutate(lookup_name = coalesce(beak_lookup[species], species)) %>%
  left_join(beakspace$space %>% select(animal, starts_with("PC")), by = c("lookup_name" = "animal")) %>%
 filter(species != "Pterodroma_cahow") %>% #remove Pterodroma cahow from df as no pc scores
  select(-lookup_name)


#join 8 PC scores to species in British Isles to formspace
form_traits <- tibble(species = spp_name) %>%
  mutate(lookup_name = coalesce(beak_lookup[species], species)) %>%
  left_join(form_scores %>% select(species, starts_with("PC")), by = c("lookup_name" = "species")) %>%
  filter(species != "Pterodroma_cahow") %>% #remove Pterodroma cahow from df as no pc scores
  select(-lookup_name)

#get species with missing trait data and no trait data at all
shape_na_any <- shape_traits %>% filter(if_any(-species, is.na)) %>%
  print()
shape_na_all <- shape_traits %>% filter(if_all(-species, is.na)) %>%
  print()

form_na_any <- form_traits %>% filter(if_any(-species, is.na)) %>%
  print()
form_na_all <- form_traits %>% filter(if_all(-species, is.na)) %>%
  print()

#264 species in beakspaces



#-----------explore variance----------

#shapespace

pc_matrix_shape <- shape_traits %>% select(starts_with("PC")) %>% as.matrix()
var_per_axis <- apply(pc_matrix, 2, var, na.rm = TRUE)
pct_var_local <- var_per_axis / sum(var_per_axis) * 100
round(pct_var_local, 1)
round(cumsum(pct_var_local), 1)

#PC1:5 explain 97.5% of variance (>95% threshold; Guillerme et al., 2023). 

#formspace

pc_matrix_form <- form_traits %>% select(starts_with("PC")) %>% as.matrix()
var_per_axis <- apply(pc_matrix, 2, var, na.rm = TRUE)
pct_var_local <- var_per_axis / sum(var_per_axis) * 100
round(pct_var_local, 1)
round(cumsum(pct_var_local), 1)

#PC1:3 explain 98.1% of variance (>95% threshold; Guillerme et al., 2023). 


#repeat for size + shapespace i.e. "formspace"----
#get size pca scores from Thomas

#Add migratory status ------

mig_status_df <- mon_presence %>%
  group_by(species) %>%
  summarise(
    mig_status = case_when(
      all(local_status == "migrant")   ~ "migrant",
      all(local_status == "resident")  ~ "resident",
      any(local_status == "resident") & any(local_status == "migrant") ~ "partial_migrant",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )

table(mig_status_df$mig_status) #110 migrants, 79 partial migrants, 75 residents

sheff_shape   <- left_join(shape_traits,    mig_status_df, by = "species")
sheff_form   <- left_join(form_traits,     mig_status_df, by = "species")


#save trait spaces
trait_spaces <- list(sheff_shape = sheff_shape, #8 PCs
                     sheff_form = sheff_form) #a mammoth 231 PCs

my_path <- "~/Data/Processed/"

save.pca.scores <- function(name, scores, path, suffix = ".csv") {
  write.csv(scores, file = paste0(path, name, "_pca_scores_plus_status", suffix), row.names = TRUE)
}

invisible(Map(save.pca.scores, names(trait_spaces), trait_spaces, path = my_path))


