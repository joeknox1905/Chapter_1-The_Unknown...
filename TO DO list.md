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
[] re-run aoh pipeline
[] go through comments and enact for PCA (rerun without kipps distance), 
[] visualise PCAs in traitspace (ggplots/ordiplot/ellipse plot)- functionalise plots; using custom functions to get subsets for plotting- then can refine more metrics e.g. edge ones
[] plot migratory status against each trait - could then do linear regression
[] Work on report and ask thomas about realised vs fundamental niche (modelling groups relative to full regional space or their local trait space, at the moment it is relative to the full britain and ireland space to get an idea of how species are occupying niches relative to the full niche space availiable (the limits of which cannot be defined by a single arbitrary geogrpahical border))
[] species-level metrics- distribution from other group centroid (unit testing needed) specilisation and uniqueness, scrap ellipsoid contribution?
[]put unit tests in seperate script and simplify fdisprity fucntion scripts to just functions
[] what traits/things explain the difference among migrants/residents in traitspace if there is any- transience (monthly analysis)/ diet/forgaing niche (avoniche)- do literature search to explore possible seasonal effects.
[] Update README
[] migration analysis as ecoregional present year-round vs not status (not taking account of transienece)- disprity (distance between centroids, sum of variance (Guillerme, 2018) or nestedness (controls for group simply being a redundant subset of the residents/annual traitspace)/turnover (beta diversity with betapart e.g. Oliveria et al., 2023, or beta.temp- seasonal overlap; but would want to do with ellipsoid not convex hulls?) or summed kernel contribution of status for each migrant grouped within the annual/seasonal/monthly (maybe omit monthly?) trait space.
[] migration as a monthly status (scale)
[] migratory distance analysis
[] dietary analysis- migratory status by diet type i.e. models/ plot overlap
[] other continents- need to change seasons again!
[] activity period collation


#n.b. use this command to bash each time to update bib file and report at the same time latexmk -pdf main.tex

# later down the line will need to control for pseudoreplication- i.e. spatial autocorrelation if using modelling e.g. species from ecoregions that are more close together (will need long/lat) for this, may be more similar in their functional uniqueness due to more similar cliamtic, habitat, soil conditions etc...