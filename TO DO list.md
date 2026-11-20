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
[] re-run aoh pipeline
[] go through comments and enact for PCA (rerun without kipps distance), 
[] visualise PCAs in traitspace (ggplots/ordiplot/ellipse plot)- functionalise plots; using custom functions to get subsets for plotting- then can refine more metrics e.g. edge ones
[] Work on report
[] species-level metrics- distribution from other group centroid, funcitonal specialization and uniquenesss?, scrap ellipse cont?
[]- blombergs k and whether to average dimensions for quantile metrics, discuss with Thomas
[] what traits/things explain the difference among migrants/residents in traitspace if there is any- transience (monthly analysis)/ diet/forgaing niche (avoniche)- do literature search to explore possible seasonal effects.
[] Update README
[] xtable- like this idea, repurposable tables when updating analysis/results
[] migration analysis as ecoregional present year-round vs not status (not taking account of transienece)- disprity (distance between centroids, sum of variance (Guillerme, 2018) or nestedness (controls for group simply being a redundant subset of the residents/annual traitspace)/turnover (beta diversity with betapart e.g. Oliveria et al., 2023, or beta.temp- seasonal overlap; but would want to do with ellipsoid not convex hulls?) or summed kernel contribution of status for each migrant grouped within the annual/seasonal/monthly (maybe omit monthly?) trait space.
[] migration as a monthly status (scale)
[] migratory distance analysis
[] dietary analysis- migratory status by diet type i.e. models/ plot overlap
[] other continents- need to change seasons again!
[] activity period collation


#n.b. use this command to bash each time to update bib file and report at the same time latexmk -pdf main.tex