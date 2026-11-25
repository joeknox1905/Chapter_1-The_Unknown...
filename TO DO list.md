---
title: "TO DO LIST"
author: "Joe Knox"
date: "2026-10-26"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Chapter 1:


To do:
[] R stats course
[] unit testing and correlation without beak tip to nares?
[] re-run aoh pipeline
[] get raw traits and then plot correlations of each trait- or can do after and get from literature- e.g. HWI - leaning towards inclusing this as opposed to wing length? but then dont have a wing length variable
[] go through comments and enact for PCA (rerun without kipps distance and with HWI), 
[] visualise PCAs in traitspace (ggplots/ordiplot/ellipse plot)- functionalise plots; using custom functions to get subsets for plotting- then can refine more metrics e.g. edge ones
[] Work on report and ask thomas about realised vs fundamental niche (modelling groups relative to full regional space or their local trait space, at the moment it is relative to the full britain and ireland space to get an idea of how species are occupying niches relative to the full niche space availiable (the limits of which cannot be defined by a single arbitrary geogrpahical border))
[] species-level metrics- distribution from other group centroid (unit testing needed) specilisation and uniqueness, scrap ellipsoid contribution?
[]put unit tests in seperate script and simplify fdisprity fucntion scripts to just functions
[] what traits/things explain the difference among migrants/residents in traitspace if there is any- transience (monthly analysis)/ diet/forgaing niche (avoniche)- do literature search to explore possible seasonal effects.- plot migratory status against each trait - could then do linear regression and PGLS
[] Update README
[] migration analysis as ecoregional present year-round vs not status (not taking account of transienece)- disprity (distance between centroids, sum of variance (Guillerme, 2018) or nestedness (controls for group simply being a redundant subset of the residents/annual traitspace)/turnover (beta diversity with betapart e.g. Oliveria et al., 2023, or beta.temp- seasonal overlap; but would want to do with ellipsoid not convex hulls?) or summed kernel contribution of status for each migrant grouped within the annual/seasonal/monthly (maybe omit monthly?) trait space.
[] migration as a monthly status (scale)
[] migratory distance analysis
[] dietary analysis- migratory status by diet type i.e. models/ plot overlap
[] other continents- need to change seasons again!
[] activity period collation


#n.b. use this command to bash each time to update bib file and report at the same time latexmk -pdf main.tex

# later down the line will need to control for pseudoreplication- i.e. spatial autocorrelation if using modelling e.g. species from ecoregions that are more close together (1|ecoregion) for this, may be more similar in their functional uniqueness due to more similar cliamtic, habitat, soil conditions etc...

#general pipeline:

#run aoh maps
#get assemblages
#plot preliminary trait correlations- migratory status vs each trait- PGLS?- ask Thomas
#run PCA with HWI and without, assess differences- but likely leave in HWI as highly predictice of migratory status (i.e. dispersal ability). If not need to leave wing length and body mass as central to variance in many niches.- i.e. wing length verty important for terrestrial birds.
#plot ordi plots and ellipse plots to assess general overlap
#run yearly disprity analysis
#get preliminary results (centroid dists, nestedness vs turnover, plot box plots)
#get disprity as a distribution with bootstrapping (plotting each analysis)- look at which species occupy extremes away from centroid of residents and which are fucntionally unique - what traits/diet types/foraging niches do they tend to represent- gives a more focused idea
#run models of migratory status against functional irreplaceability (specialization x uniqueness)
#repeat for seasonal and monthly analyses
#repeat for other regions


Hey Thomas, 

can we talk about PCA stuff in the meeting also- currently i am debating between including  HWI which is an important predictor of dispersal ability (and thus migration) and ditching wing length which is included in HWI, and including wing length (which is a large predictor of avian diveristy) but ditching HWI (which is derived from wing length). BAsicallly do I keep an important predictor of migration (and the wing shape component) but lose a major predictor of avian phenotypic diversity or retain an important predictor of avian divresity (and the wing size component) and lose an important predictor of migratory ability (and wing shape). Kipps distance the other wing shape variable is derived/uses both variables so I would need to cut it out either way. Also I currently have two beak length variables, one of which is included in the other, do I remove the nested one? This is bad practice from what I have seen in the literature and also may be an important component of modularity and phenotypic variation in aerial birds?