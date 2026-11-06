This script uses the gawfis



#=================================#
# Calculate functional space PCoA #
#=================================#

library(tidyverse)
library(mFD)
library(ggcorrplot)
library(gawdis)

# ----- LOAD DATA ----------

# Load winter occurrence data
data_winter <- read.csv("Data/occur_winter.csv") %>%
  select(-X) %>%
  mutate(cell_Year = paste(cell_id,Year,sep="_")) %>%
  column_to_rownames("cell_Year") %>%
  select(-c(cell_id,Year))

# Load breeding occurrence data
data_breeding <- read.csv("Data/occur_breeding.csv") %>%
  select(-X) %>%
  mutate(cell_Year = paste(cell_id,Year,sep="_")) %>%
  column_to_rownames("cell_Year") %>%
  select(-c(cell_id,Year))

# Load species names
Names <- read.csv("Data/Raw/Species_names.csv") %>%
  filter(Species %in% c(colnames(data_winter),
                        colnames(data_breeding)))

# ----- CALCULATE FUNCTIONAL TRAIT PCoA ----------

# Load trait data

# Seasonal diet, SAviTraits
traits_Diet <- read.csv("Data/Raw/SAviTraits_1-0_1.csv") %>%
  mutate(Scientific_name = Species_Scientific_Name) %>%
  filter(Scientific_name %in% Names$Scientific_name) %>%
  mutate(
    # Winter diet as average of Dec/Jan diet (the months of BBS surveys)
    Winter = (Dec + Jan) / 2, 
    # Breeding diet as average of May/Jun diet (the months of CBC surveys)
    Breeding = (May + Jun) / 2,
    # Group all non-fish vertebrates together (Ect, End, Unk),
    #   otherwise use sub categories
    Diet = case_when(Diet_Sub_Cat == "Invertebrate" ~ "Diet.Inv",
                     Diet_Sub_Cat == "Scavenger" ~ "Diet.Scav",
                     Diet_Sub_Cat == "Fish" ~ "Diet.Fish",
                     Diet_Cat == "Vertebrates" ~ "Diet.Vert",
                     Diet_Sub_Cat == "Fruit" ~ "Diet.Fruit",
                     Diet_Sub_Cat == "Nectar" ~ "Diet.Nect",
                     Diet_Sub_Cat == "Seed" ~ "Diet.Seed",
                     Diet_Sub_Cat == "Other" ~ "Diet.Herb")) %>%
  group_by(Scientific_name, Diet) %>%
  mutate(Winter = sum(Winter),
         Breeding = sum(Breeding)) %>%
  slice(1) %>% ungroup() %>%
  select(Scientific_name, Diet, Winter, Breeding) %>%
  pivot_longer(cols=c(Winter,Breeding), names_to="Season") %>%
  pivot_wider(id_cols=c(Scientific_name, Season),
              names_from = Diet) %>%
  select(Scientific_name, Season, Diet.Inv, Diet.Vert, Diet.Fish,
         Diet.Scav, Diet.Fruit, Diet.Nect, Diet.Seed, Diet.Herb)

# Foraging strata, EltonTraits
traits_ForStrat <- read.csv("Data/Raw/BirdFuncDat.csv") %>%
  mutate(Trait_name = Scientific) %>%
  filter(Trait_name %in% Names$Trait_name) %>%
  # Use same foraging statra categorie as EltonTraits
  select(Trait_name, ForStrat.watbelowsurf, ForStrat.wataroundsurf,
         ForStrat.ground, ForStrat.understory, ForStrat.midhigh,
         ForStrat.canopy, ForStrat.aerial)

# Morphological traits, AVONET
# First get absolute traits, then calculate relative traits
traits_Morpho_abs <- read.csv("Data/Raw/AVONET2_eBird.csv") %>%
  mutate(Scientific_name = Species2) %>%
  filter(Scientific_name %in% Names$Scientific_name) %>%
  # Calculate log trait values for all morpho traits except HWI
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
traits_Morpho_rel <- traits_Morpho_abs %>%
  # Calculate morphological traits relative to body mass
  mutate(RelBeakLengthCulmen = lm(LogBeakLengthCulmen ~ LogBodyMass,
                                  data=traits_Morpho_abs)$residuals,
         RelBeakLengthNares = lm(LogBeakLengthNares ~ LogBodyMass,
                                 data=traits_Morpho_abs)$residuals,
         RelBeakWidth = lm(LogBeakWidth ~ LogBodyMass,
                           data=traits_Morpho_abs)$residuals,
         RelBeakDepth = lm(LogBeakDepth ~ LogBodyMass,
                           data=traits_Morpho_abs)$residuals,
         RelTarsusLength = lm(LogTarsusLength ~ LogBodyMass,
                              data=traits_Morpho_abs)$residuals,
         RelKippsDistance = lm(LogKippsDistance ~ LogBodyMass,
                               data=traits_Morpho_abs)$residuals,
         RelWingLength = lm(LogWingLength ~ LogBodyMass,
                            data=traits_Morpho_abs)$residuals,
         RelSecondary1 = lm(LogSecondary1 ~ LogBodyMass,
                           data=traits_Morpho_abs)$residuals,
         RelTailLength = lm(LogTailLength ~ LogBodyMass,
                            data=traits_Morpho_abs)$residuals) %>%
  select(Scientific_name, LogBodyMass, RelBeakLengthCulmen,
         RelBeakLengthNares, RelBeakWidth, RelBeakDepth,
         RelTarsusLength, RelKippsDistance, Hand.Wing.Index,
         RelWingLength, RelSecondary1, RelTailLength)
# Check correlations among morphological traits
ggcorrplot(cor(traits_Morpho_rel[,2:12]), hc.order = TRUE)
# Based on correlations, use 6 groups of mophological traits:
#   Body mass (LogBodyMass)
#   Beak length (RelBeakLengthCulmen, RelBeakLengthNares)
#   Beak shape (RelBeakWidth, RelBeakDepth)
#   Tarsus length (RelTarsusLength)
#   Wing shape (RelKippsDistance, Hand.Wing.Index)
#   Wing/tail length (RelWingLength, RelSecondary1, RelTailLength)

# Combine trait data
spp_traits <- Names %>%
  right_join(traits_Diet, by="Scientific_name") %>%
  left_join(traits_ForStrat, by="Trait_name") %>%
  left_join(traits_Morpho_rel, by = "Scientific_name") %>%
  mutate(Species_Season = paste(Phylo_name, Season, sep="_")) %>%
  column_to_rownames("Species_Season") %>%
  select(-c(Species, Common_name, Scientific_name,
            Phylo_name, Phylo.name, Trait_name, Season))
  
# Create distance matrix of species using gawdis which weights traits whilst accounting for group correlations (de Bello et al., 2021)
dist<- gawdis(spp_traits, w.type ="analytic",
                       groups=c(1,1,1,1,1,1,1,1, # Diet
                                2,2,2,2,2,2,2, # Foraging strata
                                3, # Body mass
                                4,4, # Beak length
                                5,5, # Beak shape
                                6, # Tarsus length
                                7,7, # Wing shape
                                8,8,8), # Wing/tail length
                       groups.weight=TRUE)
attr(sp_dist_bird, "weights") # Weights of each group

# Create PCoA using distance matrix
fspaces_quality_bird <- mFD::quality.fspaces(
  sp_dist             = sp_dist_bird,
  maxdim_pcoa         = 12,
  deviation_weighting = c("absolute", "squared"),
  fdist_scaling       = c(TRUE, FALSE),
  fdendro             = "average")

# Check quality of PCoA functional space
fspaces_quality_bird$quality_fspaces
apply(fspaces_quality_bird$quality_fspaces, 2, which.min)
# Lowest mad scores is 11 axes, but that would be computationally intractable
# Use 4 axes for analysis (maximum allowed for calculating beta diversity)

# Save functional traits PCoA coordinates, 4 axes
# Species have separate winter and breeding functional trait coordinates
func_PCoA <- fspaces_quality_bird$"details_fspaces"$"sp_pc_coord"[,1:4]
saveRDS(func_PCoA, "Data/func_PCoA.RDS")

# ----- INTERPRET FUNCTIONAL TRAIT PCoA (Fig. 3) ----------

# Check correlation between functional axes and traits
# Diet traits
axes_vs_traits_1 <- mFD::traits.faxes.cor(
  sp_tr = spp_traits[,1:8],
  sp_faxes_coord = func_PCoA, plot = TRUE)
View(axes_vs_traits_1$"tr_faxes_stat"[axes_vs_traits_1$"tr_faxes_stat"$"value" > .1,])
axes_vs_traits_1$tr_faxes_plot
# Foraging strata & body mass
axes_vs_traits_2 <- mFD::traits.faxes.cor(
  sp_tr = spp_traits[,9:16],
  sp_faxes_coord = func_PCoA, plot = TRUE)
View(axes_vs_traits_2$"tr_faxes_stat"[axes_vs_traits_2$"tr_faxes_stat"$"value" > .1,])
axes_vs_traits_2$tr_faxes_plot
# Morphological traits
axes_vs_traits_3 <- mFD::traits.faxes.cor(
  sp_tr = spp_traits[,17:26],
  sp_faxes_coord = func_PCoA, plot = TRUE)
View(axes_vs_traits_3$"tr_faxes_stat"[axes_vs_traits_3$"tr_faxes_stat"$"value" > .1,])
axes_vs_traits_3$tr_faxes_plot

# Save plots of functional space including all species
func_df <- data.frame(func_PCoA) %>%
  rownames_to_column("Species_Season") %>%
  separate(Species_Season, c("Genus","Species","Season"), sep="_") %>%
  mutate(Phylo_name = paste(Genus,Species,sep="_")) %>%
  left_join(Names, by="Phylo_name") %>%
  group_by(Common_name) %>%
  mutate(diff_PC1 = abs(PC1[1]-PC1[2]),
         diff_PC2 = abs(PC2[1]-PC2[2]),
         diff_PC3 = abs(PC3[1]-PC3[2]),
         diff_PC4 = abs(PC4[1]-PC4[2]),
         tota_diff = diff_PC1 + diff_PC2 +
           diff_PC3 + diff_PC4)

hull_14 <- func_df[chull(func_df$PC1, func_df$PC4),]
hull_23 <- func_df[chull(func_df$PC2, func_df$PC3),]
min_14 <- min(c(func_df$PC1,func_df$PC4))
max_14 <- max(c(func_df$PC1,func_df$PC4))
min_23 <- min(c(func_df$PC2,func_df$PC3))
max_23 <- max(c(func_df$PC2,func_df$PC3))
# Save plot, axes 1 and 4
svg("Results/func_14.svg", width = 4, height = 4)
ggplot(as.data.frame(func_PCoA), aes(x=PC1, y=PC4)) +
  theme_classic() +
  theme(aspect.ratio=1,
        axis.text = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=10)) +
  geom_polygon(data=hull_14, fill="white", color="gray50") +
  geom_point(color="black", size=.5, alpha=.5) +
  coord_cartesian(xlim=c(min_14, max_14),
                  ylim=c(min_14, max_14))
dev.off()
# Save plot, axes 2 and 3
svg("Results/func_23.svg", width = 4, height = 4)
ggplot(as.data.frame(func_PCoA), aes(x=PC2, y=PC3)) +
  theme_classic() +
  theme(aspect.ratio=1,
        axis.text = element_text(color="black",size=10),
        axis.title = element_text(color="black",size=10)) +
  geom_polygon(data=hull_23, fill="white", color="gray50") +
  geom_point(color="black", size=.5, alpha=.5) +
  coord_cartesian(xlim=c(min_23, max_23),
                  ylim=c(min_23, max_23))
dev.off()
