## ------------Distance between two group centroids in a traitspace, two arguments (matrix 1, matrix 2)------------
pairwise.distance.between.centroids <- function(matrix, matrix2) {
    ## calculate the centroid of the matrices of n-dimensional traitspace
    sqrt(sum((colMeans(matrix) - colMeans(matrix2))^2))
} # sum of the euclidean distance between migrant group and resident centroids across multidimensional traitspace


## Test if disprity understands it
make.metric(pairwise.distance.between.centroids) # in this case should be level 1 (only one output) and between.groups

# call argument
groups <- custom.subsets(data, group = c(1, 2))
dispRity(metric = pairwise.distance.between.centroids, groups = groups, between.groups = TRUE)
