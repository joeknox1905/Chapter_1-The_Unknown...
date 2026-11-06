#joining niche and phylogenetic data (Birdlife ) to species PC scores to identify phylogenetic clustering and the distribution of migrants
#across dietary and foraging niches

#Load data

#trait data
avonet<- read.csv("/home/bip26jk/Desktop/Chapter 1/Data/Processed/avonet_traits_BI_avg_raw_yearly.csv")
avonet_crosswalk<- read.csv("/home/bip26jk/Desktop/Chapter 1/Data/Raw/AVONET_Raw_Data.csv") #raw traits for taxonomic crosswalk


# ecoregional assemblages
assemb<- readRDS("/home/bip26jk/Desktop/Chapter 1/Data/Processed/seasonal_ecoregion_pools_aoh") 

#Avoniche data (32 foraging niches)
niche<- read.csv("/home/bip26jk/Desktop/Chapter 1/Data/Raw/Avoniche1_Birdlife2020.csv")

niche <- niche %>% mutate(species = str_replace_all(Species1, " ", "_"), Family1 = family, Order1= order) #

#get species not in avoniche from seasonal assemblages
no_avoniche <- setdiff(spp_name, niche$species) %>%
  print()

no_avoniche <- str_replace_all(no_avoniche , "_", " ")
avonet_crosswalk %>%
  filter(if_any(everything(), ~ .x%in% no_avoniche)) %>%
  View()

#change species name to align with Birdlife name in avoniche
avoniche_lookup <- c(
  Curruca_communis = "Sylvia_communis",
  Curruca_curruca   = "Sylvia_curruca",
  Curruca_undata    = "Sylvia_undata",
  Luscinia_svecica = "Cyanecula_svecica")


#sorting taxonomic mismatches and selecting main foraging niche from avoniche (Sayol et al., 2026)
niche <- tibble(species = spp_name,) %>%
  mutate(lookup_name = coalesce(avoniche_lookup[species], species)) %>%
  left_join(niche %>% select(species, order, family,  MainNiche ), by = c("lookup_name" = "species")) %>%
  dplyr::select(-lookup_name)



unique(niche$MainNiche)
#18 unique foraging niches

"Generalist" "IGE"        "APG"        "APP"        "PGR"       
"IGG"        "VAS"        "APD"        "VPE"        "IAE"       
"IVS"        "APA"        "SGR"        "PAS"        "VAE"       
"SEL"        "ISG"        "APL"  

#number of species per foraging niche
#some niches only have one species so will be outliers?
#main diet types better? #n.b. foraging traits + diet type could for migratory analysis but not so good for migrant by diet type analysis as probably will already be baked into variation via PCOA
#and single species niches will be more unique- in this case restrict to only AVONET?


niche %>%
  group_by(MainNiche) %>%
  summarise(n = n_distinct(species))

#most are generalist (no main foraging niche >60%); 92 species. 
#29 are aquatic predator ground (Aquatic Prey Ground—species capturing invertebrates or vertebrates while standing in aquatic habitats (including beaches, estuaries, wetlands and marshes) (e.g., storks, herons, shorebirds). Prey may be captured on the ground or on/under water. This category includes species capturing aquatic prey (e.g., fish) or terrestrial prey in aquatic habitats (e.g., grasshopper))
#16 are aquatic predator dive (species capturing invertebrates or vertebrates under water by diving from the surface (not the air) i.e. lots of freshwater species)
#11 are Invertebrates Glean Ground (species capturing invertebrates on the ground. In contrast to ISG, the search and attack manoeuvres take place on the ground (e.g., thrushes)). This includes species standing on the ground and gleaning insects from vegetation (e.g., tinamous or larks) but excludes species that jump or sally upwards to capture prey from vegetation (ISS) or the air (ISA). The ground is dry and thus excludes aquatic habitats (e.g., beaches, estuaries, wetlands, marshes [APG]).
#all other categories are represented by <=8 species


 write.csv(niche, "species_niches_phylogeny.csv")