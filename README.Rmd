---
title: "README"
author: "Joe Knox"
date: "2026-10-26"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Data/

/Raw
- /Users/joeknox/Desktop/Phd Prep/Data/BOTW_2025.gpkg
Birdlife International Range maps for all species in Britain Ireland
- FRC_""_"".tif files 
aoh rasters for year-round resident ranges for resident species (single map) or non-breeding/breeding season maps (up to two maps) for migratory species (Lumbierres et al., 2022)
- gadm.rds files for Britain, Ireland and Crown dependency mapping
- ecoregion .shp/.dbx/.shx files for mapping ecoregions to britain and Ireland extent (Dinerstein et al., 2017)
- /Users/joeknox/Desktop/Phd Prep/Data/AVONET_Raw_Data.csv 
AVONET raw trait data for all bird species (Tobias et al., 2022)
- AVONET1_BirdLife.csv
AVONET 1- Birdlife version 7 taxonomic crosswalk (Tobias et al., 2022)
- shapespace_allbirds_lvl_superorder_order.rda (Guillerme et al., 2023)
pc scores for creating trait space from 3D bird beak landmarks with size component removed (Chira et al., 2018) derived from 4 landmarks and 3 curves of semi landmarks from 3D scans of 8748 bird species aligned with Procrustes superimposition (see Cooney et al., 2017; Chira et al., 2020; Hughes et al., 2022).
- 2020_08_07_MMB_MORPHO_FORMSPACE_FULL.rds (Guillerme et al., 2023) 
pc scores for 8748 bird species derived from 4 landmarks and 3 curves of semi landmarks from 3D scans of 8748 bird species with size component retained
-bird_season_by_month.csv (Bellve et al., 2026)
monthly range types (resident, non-breeding, breeding, passage) for 1528 migratory species for pinpointing presence of seasonal migratory birds

/Processed
- pc_scores_avobeak_aoh.csv  
pc scores for 5-dimensional trait space created from log-transformed, scaled and averaged AVONET traits for all species in Britain and Ireland present from refined area of habitat maps (210 species)
- pc_scores_avofull_aoh.csv 
pc scores for full 9-dimensional trait space, created from log-transformed, scaled and averaged traits for all species in Britain and Ireland present from refined area of habitat maps (210 species)
- monthly_presence_by_gridcell.rds
species presence x grid cell x month
- monthly_presence_by_gridcell.rds
- seasonal_ecoregion_pools.rds
3 tibbles
season: species presence x ecoregion x season (winter, summer, NA (august, march, september) of species using unrefined Birdlife range maps
month: species presence x ecoregion x season (winter, summer, NA (august, march, september) of species using unrefined Birdlife range maps
- seasonal_ecoregion_pools_aoh.rds
3 tibbles
1. season: species presence x ecoregion x season (winter, summer, NA (august, march, september)) using Birdlife ranges refined with Area of Habitat Maps
2. month: species presence x ecoregion x month  using Birdlife ranges refined with Area of Habitat Maps
- aoh_ecoregion_km2.rds
Birdlife Id ("sisid") x ecoregion x range type x total area (weighted sum of 1km grid cells in layer; km2) covered by each range type in ecoregion 

#Results/
any data processed for creating figures/tables

#Functions/
- Chapter_1_functions.R 
make.pools function for making monthly and seasonal assemblages (and probably weekly and daily) from Birdlife range maps and area of habitat maps 
xtable functions for reproducing figures as results/analyses are updated


#Analysis/
- avonet_morphospace.R
Creating a 9 dimensional pca space using prcomp() function on 9 log-transformed and scaled AVONET traits averaged for each species (with kipps distance and HWI removed do to bias towards migrants and likely high correlation with other wing traits). Also selecting axes with dispRity.
- sheffield_beakspace.R
Creating a "beak space" from ordinated 3d landmarks with procrustes alignment and selecting axes with dispRity.
- Chapter 1-Temporal Trait Space Exploration.Rmd
Exploratory analysis of British Isles, comparing migratory species with resident species
-Seasonal_assemblages_by_ecoregion.Rmd
code to create monthly assemblages and seasonal assemblages from Birdlife international seasonal range maps and area of habitat maps, and monthly phenologies for migrant species (Bellve et al., 2026). Migratory status is specified as "migrant" (at least 1 non-resident range map), and "resident" (1 resident map) and "transience" specified as number of months present or the number of months present in that season for seasonal assemblages.
-ebirdST_weekly_assemblages.Rmd
creating weekly assemblages from eBird status and trends (Strimas-Mackey et al., 2025) probability of occurrence records with "migrant" (at least 1 non-resident range map), and "resident" (1 resident map) defined as transience specified as number of weeks present
- Chapter1_Temporal_Traitspaces.Rmd
Code for creating final trait spaces

#Manuscript/

- initial_thesis_plan.tex
Latex tex file for initial thesis plan element of confirmation review
- initial_thesis_plan.pdf
Pdf for initial thesis plan element of confirmation review

/Figures

- ecoregion_seasonal_status_pca_chulls_plot.png
comparing convex hull overlap of residents and migrants across seasons and ecoregions (PC1 V 2; 4 dimensions; AVONET)-  identifying nestedness vs turnover with intial unrefined birdlife range pipeline
- ecoregion_seasonal_status_pca_ellipses_plot.png
comparing ellipsoid volume overlap of residents and migrants across seasons and ecoregions (PC1 V 2; 4 dimensions; AVONET)- nestedness vs turnover
with intial unrefined birdlife range pipeline
- ecoregion_seasonal_status_pca_ellipses_plot_aoh_1.png
comparing ellipsoid volume overlap of residents and migrants across seasons and ecoregions (PC1 V 2; 4 dimensions; AVONET)- nestedness vs turnover
with refined unrefined aoh map pipeline
- ecoregion_seasonal_status_pca_ellipses_plot_aoh_2.png
comparing ellipsoid volume overlap of residents and migrants across seasons and ecoregions (PC1 V 3; 4 dimensions; AVONET)- nestedness vs turnover
with refined unrefined aoh map pipeline
- ecoregion_seasonal_status_pca_ellipses_plot_aoh_3.png
comparing ellipsoid volume overlap of residents and migrants across seasons and ecoregions (PC2 V 3; 4 dimensions; AVONET)- nestedness vs turnover with refined unrefined aoh map pipeline
- box_plot_PERMANOVA_within_group_centroid_dist_all_dims_aoh.png 
comparing spread (among group distance to centroids; AVONET) between migrants and residents in each ecoregion and season) for all 9 dimensions using vegan betadisp() function with Holm correction.
- box_plot_PERMANOVA_within_group_centroid_dist_aoh.png
comparing spread (among group distance to centroids; AVONET) between migrants and residents in each ecoregion and season) for 4 dimensions using vegan betadisp() function with Holm correction.
- pc1v3_beak_vs_full_space_aoh.png
comparing spread (among group distance to centroids; AVONET) between migrants and residents in each ecoregion and season) for 4 dimensions using vegan betadisp() function with Holm correction.
- pc1v2_beak_vs_full_space_aoh.png
pca loadings plots comparing PC1 v 2 between beak space (wing, tail and tarsus trait removed; 5 dimensional) and full space (non-beak traits also retained; 9 dimensional- no kipps distance or HWI)
- pc1v3_beak_vs_full_space_aoh.png
pca loadings plots comparing PC1 v 3 between beak space (wing, tail and tarsus trait removed; 5 dimensional) and full space (non-beak traits also retained; 9 dimensional- no kipps distance or HWI)
- pc2v3_beak_vs_full_space_aoh.png
pca loadings plots comparing PC2 v 3 between beak space (wing, tail and tarsus trait removed; 5 dimensional) and full space (non-beak traits also retained; 9 dimensional- no kipps distance or HWI)
- pc2v3_beak_vs_full_space_aoh.png
pca loadings plots comparing PC2 v 4 between beak space (wing, tail and tarsus trait removed; 5 dimensional) and full space (non-beak traits also retained; 9 dimensional- no kipps distance or HWI)



/Tables

/References
- bibliography_chap1.bib
A latex bib file for references





