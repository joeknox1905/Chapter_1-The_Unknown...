###################################################################################################################################################################################################
#------------function to split a dataframe by multiple column variables to feed into custom subsets------------------
######################################################################################################################################################################################################


# takes a dataframe with species (and other metadata) as rows and PC scores as columns
make.ecoregional.migrant.groups <- function(traitspace_df, by) {
  lapply(split(traitspace_df$species, traitspace_df[by], drop = TRUE), unique) # split traitspace by distinct ecoregion, status (and season for seasonal analysis)
}

# call function
groups_yearly <- make.ecoregional.migrant.groups(traitspace_df = , by = c("", "")) # split by column variables e.g. ecoregion or status

# then would into custom.subsets i.e. custom.subsets(trait matrix, group = groups_yearly)
