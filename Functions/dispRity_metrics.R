# This script creates between group metrics (dimension level 2) in dispRity to calculate degree of the functional turnover (or divergence between trait spaces; distance between centroids?) and nestedness 
#(i.e. subset within a subset) to charcterise uniqueness of migrant vs residents in terms of the annual trait space. A seasonal analysis is then also run
# starting with a binary defintion of migratory status present year-round vs present <12 months
# and then testing migratory status as graded scale (from present 1 month to present 12 months) (bootstrapping/rarefying- how much does a metric of functional overlap/nestedness change when we remove each species in each group from the trait space?)

##############################################################################################################################################
#pre-requisites-----
##############################################################################################################################################

library(dispRity)
library(MCMCglmm)
library(testthat)

traitspace_df<- read.csv("Data/Processed/eco_region_season_traitspace_PCA_BI.csv") # full Britain and Ireland trait space or beak-only AVONET PCA/PCOA scores, status and species (columns) by ecoregion (rows)- can either use PCOA or PCA
#TG: it's better/more common to write things in path_from_the_file (so that "/Data/..." becomes "../Data/...")

pc_cols <- names(traitspace_df)[8:13] #get names of PC columns excluding metadata
species_pcs <- unique(traitspace_df[c("species", pc_cols)]) #get unique pc cols per species # i.e. the full trait space, 1 row per species

#convert to matrix for disprity
traitspace_mat <- as.matrix(species_pcs[, pc_cols]) 
rownames(traitspace_mat) <- species_pcs$species 


##############################################################################################################################################
#-------------------------------------------------------------dispRIty metrics----------------------------------------------------------------
##############################################################################################################################################


##------------Create metric for calcualting distance between resident and migrant centroids, two arguments (matrix 1 = res, matrix 2 = mig)------------

#Do migrants differ on average to residents n local (ecoregional) trait space (i.e. do they occupy different niches?)
#N.b. tells us nothing about the spread of species within their group or the nestedness of migrants and residents, could just be skewed by a few outliers



pairwise.distance.between.centroids <- function(matrix, matrix2) {
    ## calculate the centroid of the matrices of n-dimensional traitspace
   sqrt(sum((colMeans(matrix) - colMeans(matrix2))^2))
} #sum of the euclidean distance between migrant group and resident centroids across multidimensional traitspace


## Test if disprity understands it
make.metric(pairwise.distance.between.centroids) #in this case should be level 1 (only one output) and between.groups

## Unit testing

## Step 1: make sure the function does what you tell it to do
stupid <- matrix(rnorm(25), 5, 5)
testthat::expect_equal(pairwise.distance.between.centroids(stupid, stupid), 0)

matrix <- matrix(0, 5, 2)
matrix2 <- matrix(10, 5, 2)
testthat::expect_equal(pairwise.distance.between.centroids(matrix, matrix2), sqrt(2*(10^2)))

## Step 2: try to understand how the function works in a real multidimensional example/space
set.seed(1)
space_2D <- traitspace_mat[, c("PC1", "PC2")]
random_groups <- replicate(2, sample(1:nrow(traitspace_mat), floor(nrow(traitspace_mat)/2)))
random_test <- custom.subsets(space_2D, group = list("A" = random_groups[,1], "B" = random_groups[,2]))
# plot(random_test)
expect_equal(summary(dispRity(random_test, metric = pairwise.distance.between.centroids, between.groups = TRUE))$obs, 0.088)

#non-overlapping groups
non_random_groups <- list("A" = which(space_2D[,1] <= -0.5 & space_2D[,2] <= -0.5),
"B" = which(space_2D[,1] >= 1 & space_2D[,2] >= 1))
non_random_test <- custom.subsets(space_2D, group = non_random_groups)
# plot(non_random_test)
expect_equal(summary(dispRity(non_random_test, metric = pairwise.distance.between.centroids, between.groups = TRUE))$obs, 5.382)

#completely overlapping groups
complete_overlap <- custom.subsets(space_2D, group = list("A" = random_groups[,1], "B" = random_groups[,1]))
# plot(complete_overlap)
expect_equal(summary(dispRity(complete_overlap, metric = pairwise.distance.between.centroids, between.groups = TRUE))$obs, 0)


#also compare to group.dist (minimum vector norm dist between two groups)
disprity((metric = group.dist, probs = c(colMeans())), groups = "status", between.groups = true)

## Step 3: formalising- but other way is probably better for interpretation and visualisation

## Creating a 2D uniform space
example_space <- space.maker(300, 2, runif)

example_test <- test.metric(example_space, metric = c(pairwise.distance.between.centroids),
 shifts = c("random",  "size", "density", "position")) #position may be better here for this metric i.e. do the groups change position in trait space
 
 #get results
 summary(example_test)
 #runs three replicates of the space reduction and again runs tests sequentially removing 10% of data with specified metrics and linear modelling
# a good metric will have no trend in the size direction, but have a trend for the reduction of interest i.e. if interested in size, should show sensitivity to space decreases?
#for size returns size.inner (contracting space), size.outer (hollowing space), for position returns position.top and position.bottom(removing data from one side of space or the other)


#visualise
plot(example_test)




#-------------metric for calculating proportional overlap of groups-----------------------------------------------

#how much do groups represent eachother in functional space?


quantile.overlap<- function (matrix, matrix2, probs = c(0.025, 0.975)) {
  q1 <- apply(matrix, 2, quantile, probs = probs) # apply quantile calculation to each PC dimension (column) in the group matrix
  q2<- apply(matrix2, 2, quantile, probs = probs) #returns vector of 95% quantiles for for each dimension
  intersect <- pmax(0, pmin(q1[2, ], q2[2, ]) - pmax(q1[1, ], q2[1, ])) #get minimum upper bound of groups and maximum lower bound of groups for each dimension (column)  i.e. get minimum distance between group (nearest edges), this will be the intersection (overlap), if groups do not overlap, set to 0 as distance is negative
  span <- pmax(q1[2, ], q2[2, ]) - pmin(q1[1, ], q2[1, ]) #minus maximum upper bound of groups from minimum lower bound of groups for each dimension to get combined group span across each dimension
  intersect / span # quantile overlap = proportion of overlap compared to full span.
  #no overlap = 0
  #complete overlap = 1 
  #partial overlap = >0 & <1
}

#boot matrix etc

disprity(traitspace_mat, metric = quantile.overlap, between.groups = (c(1,2), c(4,3), c(5,6), c(7,8), c(9,10))) #returns vector of overlap for each dimension or between.groups = (c(1,2), c(2,1)) listing each migrant x resident ecoregional pair


#get.disprity = returns vector of overlap across every dimension
#summary returns median across all dimensions


#unit testing----------

make.metric(quantile.overlap) #outputs matrix of proportional overlap per dimension and pairwise comparsion as expected, dimension level 2 

#test that two equal matrices would return 1 as complete overlap 
set.seed(1)
equal_matrices <- matrix(rnorm(25), 5, 5) # if matrix 1 + 2 = are the same, would expect proportional overlap to be 1 (intersect == span), so overlap = 1
testthat::expect_all_equal(quantile.overlap(equal_matrices, equal_matrices), 1) # correct

#testing on a real multi-dimensional space

#generate two random groups in a hypothetical 2d traitspace
space_2D <- traitspace_mat[, c("PC1", "PC2")]
random_groups <- replicate(2, sample(1:nrow(traitspace_mat), floor(nrow(traitspace_mat)/2)))

#testing complete overlap
test_that("groups with complete overlap return median quantile.overlap of 1", {
complete_overlap <- custom.subsets(space_2D, group = list("A" = random_groups[,1], "B" = random_groups[,1]))
# plot(complete_overlap)
expect_equal(summary(dispRity(complete_overlap, metric = quantile.overlap, between.groups = TRUE))$obs.median, 1) 
}) #yay!

#testing partial overlap

idx<- 1:nrow(space_2D)
half<- floor(length(idx)/2)

set.seed(1)
test_that("quantile.overlap between two partly overlapping groups is between 0,1", {
  partial_overlap <- custom.subsets(space_2D, group = list("A" = idx[1:half], "B" = idx[(half + 1):(2 * half)]))
  # plot(partial_overlap)
  expect_lt(summary(dispRity(partial_overlap, metric = quantile.overlap, between.groups = TRUE))$obs.median, 1)
  expect_gt(summary(dispRity(partial_overlap, metric = quantile.overlap, between.groups = TRUE))$obs.median, 0)
}) #yay!

#testing that nested groups are still counted as not completely overlapping so does not differentiate nestedness from partial overlap
set.seed(1)
test_that("quantile overlap for a nested group = >0 & <1", {
nested_groups<- list("A" = which(space_2D[,1] >-3 & space_2D[,1] <3 & space_2D[,2] >-3 & space_2D[,2] <3),
"B" = which(space_2D[,1] >-2 & space_2D[,1] <2 & space_2D[,2] >-2 & space_2D[,2] <2))
nested_overlap <- custom.subsets(space_2D, group = nested_groups)
 #plot(nested_overlap)
 expect_lt(summary(dispRity(nested_overlap, metric = quantile.overlap, between.groups = TRUE))$obs.median, 1)
  expect_gt(summary(dispRity(nested_overlap, metric = quantile.overlap, between.groups = TRUE))$obs.median, 0)
}) #yay



#non-overlapping groups: should be 0
test_that("quantile overlap = 0 for non-overlapping groups", {
non_random_groups <- list("A" = which(space_2D[,1] <= -0.5 & space_2D[,2] <= -0.5),
"B" = which(space_2D[,1] >= 1 & space_2D[,2] >= 1))
non_random_test <- custom.subsets(space_2D, group = non_random_groups)
# plot(non_random_test)
expect_equal(summary(dispRity(non_random_test, metric = quantile.overlap, between.groups = TRUE))$obs, 0)
}) #good


#-----------------------------------metric for descibing turnover - what proportion of trait space is unique to each group?-------------------------------



quantile.uniqueness<- function (matrix, matrix2, probs = c(0.025, 0.975)) {
  q1 <- apply(matrix, 2, quantile, probs = probs) 
  q2<- apply(matrix2, 2, quantile, probs = probs)
  intersect <- pmax(0, pmin(q1[2, ], q2[2, ]) - pmax(q1[1, ], q2[1, ])) #if no intersect, defualts to 0
  span <- q1[2, ] - q1[1, ]
  turnover<- span-intersect

return(turnover/span)
}
   
#call disprity
dispRity(traitspace_mat, metric=quantile.uniqueness, between.groups = list(c(matrix,matrix2), c(matrix2,matrix))) #would need to list each ecoregion pair

#turnover = span (= complete turnover, group is completely unique)
#turnover < span (= partial turnover (just partial overlap)
#turnover = 0 = no turnover (complete overlap/nestedness)

#should be calculated pairwise for intersection, but return a matrix of values for each dimension for each group- dimesion level 3

unit testing--------------


make.metric(quantile.uniqueness) #correct, outputs as a matrix, performed between groups

set.seed(1)
equal_matrices <- matrix(rnorm(25), 5, 5) # if matrix 1 + 2 = are the same, would expect turnover to be 0 (intersect = span), so turnover = 0/span = 0
testthat::expect_all_equal(quantile.uniqueness(equal_matrices, equal_matrices), 0) # correct, yay

#testing on a real multi-dimensional space

#generate two random groups in a hypothetical 2d traitspace
space_2D <- traitspace_mat[, c("PC1", "PC2")]
random_groups <- replicate(2, sample(1:nrow(traitspace_mat), floor(nrow(traitspace_mat)/2)))

#testing complete overlap, should be 0 turnover as intersect = span
test_that("groups with complete overlap return median quantile.uniqueness of 0", {
complete_overlap <- custom.subsets(space_2D, group = list("A" = random_groups[,1], "B" = random_groups[,1]))
# plot(complete_overlap)
expect_all_equal(summary(dispRity(complete_overlap, metric = quantile.uniqueness, between.groups = list(c("A","B"), c("B","A"))))$obs.median, 0)
}) #yay, all good and returns two matrices for each pairwise comparison 

#testing partial overlap

idx<- 1:nrow(space_2D)
half<- floor(length(idx)/2)

set.seed(1)
test_that("quantile.uniqueness between two partly overlapping groups is between 0,1", {
  partial_overlap <- custom.subsets(space_2D, group = list("A" = idx[1:half], "B" = idx[(half + 1):(2 * half)])) #assign half row values to group, A and half to group B
  # plot(partial_overlap)
  expect_all_true(summary(dispRity(partial_overlap, metric = quantile.uniqueness, between.groups = list(c("A","B"), c("B","A"))))$obs.median <1)
  expect_all_true(summary(dispRity(partial_overlap, metric = quantile.uniqueness, between.groups = list(c("A","B"), c("B","A"))))$obs.median > 0)
}) #yay!

#testing that nested group has no turnover
set.seed(1)
test_that("quantile.uniqueness for a nested group = 0 and the non-nested group = >0", {
nested_groups<- list("A" = which(space_2D[,1] >-3 & space_2D[,1] <3 & space_2D[,2] >-3 & space_2D[,2] <3),
"B" = which(space_2D[,1] >-2 & space_2D[,1] <2 & space_2D[,2] >-2 & space_2D[,2] <2))
nested_overlap <- custom.subsets(space_2D, group = nested_groups)
 expect_all_true(get.disparity(dispRity(nested_overlap, metric = quantile.uniqueness, between.groups = list(c("A","B"), c("B","A"))))$`A:B` >0)
expect_all_equal(get.disparity(dispRity(nested_overlap, metric = quantile.uniqueness, between.groups = list(c("A","B"), c("B","A"))))$`B:A`, 0)
}) #yay- accurately characterises turnover of each group reative to the other group

#non-overlapping groups: should be 1
test_that("quantile uniqueness = 1 for non-overlapping groups", {
non_random_groups <- list("A" = which(space_2D[,1] <= -0.5 & space_2D[,2] <= -0.5),
"B" = which(space_2D[,1] >= 1 & space_2D[,2] >= 1))
no_overlap <- custom.subsets(space_2D, group = non_random_groups)
# plot(non_random_test)
 expect_all_equal(get.disparity(dispRity(no_overlap, metric = quantile.uniqueness, between.groups = list(c("A","B"), c("B","A"))))$`A:B`, 1)
expect_all_equal(get.disparity(dispRity(no_overlap, metric = quantile.uniqueness, between.groups = list(c("A","B"), c("B","A"))))$`B:A`, 1)
}) #good 



#-----------------------------metric for describing nestedness- what proportion of a group's traits are a subset of the other group's?--------------------------------------

quantile.nestedness<- function (matrix, matrix2, probs = c(0.025, 0.975)) {
  q1 <- apply(matrix, 2, quantile, probs = probs) 
  q2<- apply(matrix2, 2, quantile, probs = probs)
  intersect <- pmax(0, pmin(q1[2, ], q2[2, ]) - pmax(q1[1, ], q2[1, ])) 
  span1 <- q1[2, ] - q1[1, ]
  span2 <- q2[2, ] - q2[1, ]
  return(intersect/span1, intersect/span2)
}


#call disprity
disprity(trait_mat, metric=quantile.uniqueness, between.groups = TRUE)

#or this would return 1 value per pairwise comaprison (1 for each group), dimension level 2?

quantile.nestedness<- function (matrix, matrix2, probs = c(0.025, 0.975)) {
  q1 <- apply(matrix, 2, quantile, probs = probs) 
  q2<- apply(matrix2, 2, quantile, probs = probs)
  intersect <- pmax(0, pmin(q1[2, ], q2[2, ]) - pmax(q1[1, ], q2[1, ])) 
  span <- q1[2, ] - q1[1, ]

  return(intersect/span)
}

#call disprity
disprity(trait_mat, metric=quantile.nestedness, between.groups = list(c(1,2), c(2,1))) #should be 20 combinations for 5 ecoregions

#should be calculated pairwise and return 2 matrices, 1 per group of each dimension -dimesion level 3 to see which group is nested in which #summary would return median of all dimensions
#intersect of 0/span = 0 (non-overlapping group)
#intersect =span = 1 (completely nested group; group is a complete subset of the other group)
#intersect < span = partially nested group (just partial overlap)



#unit testing------------------


set.seed(1)
equal_matrices <- matrix(rnorm(25), 5, 5) # if matrix 1 + 2 = are the same, would expect nestedness to be 1 (intersect = span), so nestedness = 1/1
testthat::expect_all_equal(quantile.nestedness(equal_matrices, equal_matrices), 1) #correct

#generate two random groups in a hypothetical 2d traitspace
space_2D <- traitspace_mat[, c("PC1", "PC2")]
random_groups <- replicate(2, sample(1:nrow(traitspace_mat), floor(nrow(traitspace_mat)/2)))

#testing complete overlap, should be 1 as intersect = span so nestedness = 1 
test_that("groups with complete overlap return median quantile.nestedness of 1", {
complete_overlap <- custom.subsets(space_2D, group = list("A" = random_groups[,1], "B" = random_groups[,1]))
# plot(complete_overlap)
expect_all_equal(summary(dispRity(complete_overlap, metric = quantile.nestedness, between.groups = list(c("A","B"), c("B","A"))))$obs.median, 1)
}) #yep, all good

idx<- 1:nrow(space_2D)
half<- floor(length(idx)/2)

#testing partial overlap: should be between 0 and 1 as groups are partially nested
test_that("quantile.nestedness between two partly overlapping groups is between 0,1", {
  partial_overlap <- custom.subsets(space_2D, group = list("A" = idx[1:half], "B" = idx[(half + 1):(2 * half)])) #assign half row values to group, A and half to group B
  # plot(partial_overlap)
  expect_all_true(summary(dispRity(partial_overlap, metric = quantile.nestedness, between.groups = list(c("A","B"), c("B","A"))))$obs.median <1)
  expect_all_true(summary(dispRity(partial_overlap, metric = quantile.nestedness, between.groups = list(c("A","B"), c("B","A"))))$obs.median > 0)
}) #good


#test that non-overlapping groups = 0 
set.seed(1)
test_that("quantile nestedness = 0 for non-overlapping groups", {
non_random_groups <- list("A" = which(space_2D[,1] <= -0.5 & space_2D[,2] <= -0.5),
"B" = which(space_2D[,1] >= 1 & space_2D[,2] >= 1))
no_overlap <- custom.subsets(space_2D, group = non_random_groups)
# plot(non_random_test)
expect_all_equal(get.disparity(dispRity(no_overlap, metric = quantile.nestedness, between.groups = list(c("A","B"), c("B","A"))))$`A:B`, 0)
expect_all_equal(get.disparity(dispRity(no_overlap, metric = quantile.nestedness, between.groups = list(c("A","B"), c("B","A"))))$`B:A`, 0)
}) #good 

set.seed(1)
test_that("quantile.nestedness for a nested group = 0 and the non-nested group = <1", {
nested_groups<- list("A" = which(space_2D[,1] >-3 & space_2D[,1] <3 & space_2D[,2] >-3 & space_2D[,2] <3),
"B" = which(space_2D[,1] >-2 & space_2D[,1] <2 & space_2D[,2] >-2 & space_2D[,2] <2))
nested_overlap <- custom.subsets(space_2D, group = nested_groups)
expect_all_true(get.disparity(dispRity(nested_overlap, metric = quantile.nestedness, between.groups = list(c("A","B"), c("B","A"))))$`A:B` <1)
expect_all_equal(get.disparity(dispRity(nested_overlap, metric = quantile.nestedness, between.groups = list(c("A","B"), c("B","A"))))$`B:A`, 1)
}) #yay- accurately characterises nestedness of each group reative to the other group





#-----------------------------metric for calculating displacement from centroid (deviance pf migrant group from average trait combos)------------------------

#i.e. is migrant groups average trait combinations (centroid) at the center (centroid, 1) or away from centre (<1) of traitspace and therefore does it represent typical/redundant niches (generalist) or more unique niches (specialist)?
#more of a positional measure than a size measure (Guillerme et al., 2020). 
#also may be good to calcualate ecoregional migrant group displacement compared to full BI trait space not just annual trait space as would give an idea of whether 
#migrant group is unique locally but more redundant regionally- if that is the case, migrants are more functionally similar to migrants/residents in other ecoregions, suggesting either local migration from other ecoregions to the ecoregion (could confirm with migratory distance analysis)
#or that residents are actually more functionally unique in this ecoregion and that ecoregion shapes resident communities more than migrant communities.

displacements() in dispRity

  #calculated as:
  
  avg.group.displacements<- function(matrix) {
    reference_centroid<- colMeans(trait_cols)
    group_centroid<-  colMeans(matrix)
  
  sqrt(sum(group_centroid)^2/sqrt(sum(group_centroid-reference_centroid)^2) #@@@JK:TG is this correct?

#need to call between.groups = FALSE
result_binary  <- calculate.disprity.by.ecoregion(BI_traitspace, group = "status", metric = distance.between.centroids , trait_cols = trait_cols, between = FALSE) 

#could also calculate point.dist- distance from centroid of matrix 2 to matrix which may be more what you want
dispRity(point.dist(point=colMeans))

#----------------metric for calculating distribution of species within each group around centroid----------------


#i.e. "functional" arrangement within each group- niche packing vs expansion in each group
#if normal, means species tend to be more clustered  within their group around average trait combinations (whatever trait combinations represent the centroid), niche packing within group- maybe more fine-scale specialisation
#if uneven, suggests dispersion from centroid- niche expansion within group

distribution.of.group.around.group.centroid <- function(matrix) {
group_centroid<- colMeans(matrix) #centroid of group
for (i in 1:nrow(matrix)) { #for each species row in ecoregional annual trait space
species_row <- matrix[i, ]#get pc score for each dimension for each pc

#avg e.g. if smaller avg distance in migrants than migrants more clustered around centroid and more clumped? to compare between groups- i guess kind of similae to 
sqrt(sum((species_row - colMeans(matrix2))^2)) #distance from group centroid of each species

need to call between.groups = FALSE
result_binary  <- calculate.disprity.by.ecoregion(BI_traitspace, group = "status", metric = distribution.of.group.around.group.centroid , trait_cols = trait_cols, between = FALSE) 

# plot histogram of distribution function and test difference using t test (not permutation as already have done bootstrapping to account for null tests)
pdf("distribution_around_centroid_per_"group".pdf", width = 14, height = 10)
ggplot(result_binary, aes(x= $obs)), 
geom_histogram(binwidth = 10) +
facet_grid(group~.) #migrants and residents side by side
dev.off()

} # plot histogram of distribution function and test difference using t test (not permutation as already have done bootstrapping to account for null tests)



#------------------------metric for calculating species-level contribution to trait space using ellipsoid volume, based on leave one out approach of Cardoso et al. (2026) but using ellipsoid volumes as less sensitive to outliers (Guillerme et al., 2026)-----------------------------
  
#Are migrants locally (ecoregionally) functionally unique or redundant?
#tells us more about the species-level contribution (redundancy vs complementarity) to trait space and its removal shifts functional richness, less volume loss = less unique


  ellipse.contribution <- function (matrix) { #requires matrix1 (ecoregional annual trait space), matrix2 (ecoregional annual trait space-1)
  


    for (i in 1:nrow(matrix)) { #for each species row in ecoregional annual trait space
      species_row <- matrix[i, ]
      matrix2 <- matrix[-species_row] #remove species from trait space 
      # Calculate ellipsoid volume for each species
    }
    return(ellipsoid.volume(matrix)-(ellipsoid.volume(matrix2)) = #ellipsoid volume contribution of each species return ellipsoid volume of species trait space occupation i.e. their functional diversity- larger = more unique/ less redundant contribution as species will become more clustered, suggesting the species is further from other species in trait space

  } 
  
#Ellpisoid volume = π^(d/2) / Γ(d/2 + 1) × ∏ᵢ₌₁ᵈ √λᵢ

# Where:
 #π^(d/2) / Γ(d/2 + 1) is the volume of unit ball in d dimensions, Γ is the gamma function

#And: 
# ∏ᵢ₌₁ᵈ √λᵢ is the product of the ellipsoids semi-axis lengths, with √λᵢ  as the standard deviation along each pc axis

 
#test that metric works as expected
  make.metric(ellipse.contribution) #will be level 1 (one output) and not between groups as species-level. output should be a 

#call function
result_binary  <- calculate.disprity.by.ecoregion(BI_traitspace, group = "status", metric = distance.between.centroids , trait_cols = trait_cols, between = FALSE) 
#change for ellipsoid contribution (between = FALSE as species level)

# then would do stats e.g. function to get mean ellipse contribution per group from disprity summary and calculate difference between migrants and residents

calculate.mean.ellipse.contrib.sig <- function(disparity_result) {
disparity_result<- disparity_result %>% #need disparity summary from calculate.disprity.by.ecoregion with metric = ellipse.contribtion) 
group_by(group =group) %>%
summarise(mean_ellips_cont = mean(disparity_result$obs))  #get mean ellipse contribution (bootstrapped and rarefied) of each group

t.test((mean_ellips_cont ~ group, data = disparity_result )) # welches t test accounts for unequal variance of groups, 
#data points are independent so satisfies assumption as different species in each migrant group per ecoseason (partial migrants are residents)

}

#call function
calculate.mean.ellipse.contrib.sig(disparity_result = result_binary) #or result_monthly



#other forms of species-level functional contribution calculations based on (Marino et al., 2024):

#functional specialization

#euclidean distance from centroid of multi-dimesnisonal  traitspace of each species (then model migrant vs resident)


#functional uniqueness
#sum of euclidean distances between species and its 5 nearest neighbours in multi-dimensional traitspace

#fucntional irreplaceability = functional specialisation x functional uniqueness?





