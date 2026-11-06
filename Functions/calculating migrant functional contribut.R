# This script creates between group metrics (dimension level 2) in dispRity to calculate degree of the functional turnover (or divergence between trait spaces; distance between centroids?) and nestedness 
#(i.e. subset within a subset) to charcterise uniqueness of migrant vs residents in terms of the annual trait space. A seasonal analysis is then also run
# starting with a binary defintion of migratory status present year-round vs present <12 months
# and then testing migratory status as graded scale (from present 1 month to present 12 months) (bootstrapping/rarefying- how much does a metric of functional overlap/nestedness change when we remove each species in each group from the trait space?)



also can make matrix a function to get migrants e.g.

#make mig groups

need to refine to loop for each ecoregion

eco_region_traitspace<- (....)


detect.migrants_binary<- function(eco_region_traitspace) {#traitspace with status

migrant<- traitspace %>%
filter(status == "migrant") %>% 
column_to_rownames("species") %>% 
select (-status)

resident <- traitspace %>% 
filter(status == "resident") %>% 
olumn_to_rownames("species") %>% 
select (-status)


return(migrant, resident)
}

#call function
detect.migrants_binary(eco_region_traitspace) #want it to return 2 matrices, one of migrant, one of resident

need to refine to loop for each ecoregion

detect.migrants_monthly<- function(traitspace) {#traitspace with status

migrant<- traitspace %>%
filter(status == "migrant") %>% 
column_to_rownames("species") %>% 
select (-status)

resident <- traitspace %>% 
filter(status == "resident") %>% 
olumn_to_rownames("species") %>% 
select (-status)

#return 12 matrices for each ecoregion
return(migrant1, migrant2, migrant3, migrant4, migrant5, migrant6, migrant7, migrant8, migrant9, migrant10, migrant11, resident)
}


#call function 



#then can input matrices from detect.migrants ()

#loop through each ecoregion

## Create your metric, two arguments (matrix 1 = res, matrix 2 = mig), later change to all except mig_1 vs mig_1 etc
distance.between.centroids <- function(matrix, matrix2) {
    ## calculate the centroid of the matrices
    centroid1 <- colMeans(matrix) #centroid of res
    centroid2 <- colMeans(matrix2) #centroid of mig
    ## measuring the distances between them
    return(sum(sqrt(centroid1^2-centroid2^2))) #euclidean distance
}


## Test if disprity understands it
make.metric(distance.between.centroids) #in this case should be level 1 (only one output) and between.groups


## use in dispRity
dispRity(..., metric = distance.between.centroids, between.groups = TRUE)

#for species contribution calculations (avg per group?)



boot.matrix(mig_group) = rarefy = TRUE? or distance.between.centroids(boot.matrix = true)
#read manual