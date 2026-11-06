#Functionalising code


#make.pools function-----------
#creating monthly and seasonal pools from Birdlife international range maps and area of habitat maps



#Xtable functions--------
#for repurposing tables when analyses and results are updated (see Thomas' github)


#creating a standardisation function--------
head(iris)
summary(iris) #ranges of different traits are very different so we will
#need to standardise (observed-mean/SD; like SES!)

#standardise petal length:

# Standardize petal length
mean_petal_length <- mean(iris$Petal.Length)
sd_petal_length <- sd(iris$Petal.Length)
standard_petal_length <- (iris$Petal.Length - mean_petal_length)/sd_petal_length

summary(standard_petal_length) #ensure mean is 0 i.e. that observed values are standardised
sd(standard_petal_length) #ensure SD is 1

# Standardize petal width
mean_petal_width <- mean(iris$Petal.Width)
sd_petal_width <- sd(iris$Petal.Width)
standard_petal_width <- (iris$Petal.Width - mean_petal_width)/sd_petal_width
summary(standard_petal_width)


#but we will have to do this for all 7 traits so how can we optimise this?
#we can create a fuction!

#first write comments to outline what steps we need to achieve
standardize <- function(data) {
  # Calculate mean of data
  # Calculate standard deviation of data
  # Calculate standardized scores
  
}

standardize <- function(data) {
  # Calculate mean of data
  mean_data <- mean(data)
  # Calculate standard deviation of data
  # Calculate standardized scores
}

standardize <- function(data) {
  # Calculate mean of data
  mean_data <- mean(data)
  # Calculate standard deviation of data
  sd_data <- sd(data)
  # Calculate standardized scores
  stand_data <- (data - mean_data)/sd_data
  return(stand_data) #what we want function to output
}

#do not refer to iris df as want it to be generalisable to other data too!

#call function to standardise petal length (returns vector for each value(element) in the column
standard_petal_length <- standardize(data = iris$Petal.Length)

#write comments to outline function for users
# Standardize a vector of numerical data
# data: vector of numerical data
# returns: vector of original data standardized to have mean = 0, sd = 1
standardize <- function(data) {
  # Calculate mean of data
  mean_data <- mean(data)
  # Calculate standard deviation of data
  sd_data <- sd(data)
  # Calculate standardized scores
  stand_data <- (data - mean_data)/sd_data
  return(stand_data)
}

#BUT what if we can make this simpler by applying to lots of columns at once?
#use lapply
#
iris_stand <- apply(X = iris[, -5], #x is the element you want to apply the function to
                    #in this case every column in iris df except species column as is not numeric
                    MARGIN = 2, #MARGIN =1 is do calculation for each row,  2 is each column
                    FUN = standardize)

#convert output back to dataframe (rather than vector, to view)
iris_stand <- as.data.frame(iris_stand) 
head(iris_stand)

# Add back species column
iris_stand$Species <- iris$Species

#using source files---------------------

#create new R script which contains the function (could contain all functions for a specific project)
#e.g. 
message("Loading functions from standarization-functions.R")
source(file = "standardization-functions.R")
#will load all functions in r script, then can call in the script when needed

#making functions generalisable--------

#to make function generalizable we can't specify which column numbers are not numeric as will change
#so we use is.numeric instead with sapply
sapply(X = iris, FUN = is.numeric)
#returns logical vector of which columns are/are not numeric

#then clean up function using col_numeric to perform function only on numeric columns within a dataframe

standardize.df <- function(df) {
  # Identify numeric columns
  col_numeric <- sapply(X = df, FUN = is.numeric) #function within a function (is.numeric already exists)
  # Run standardize on all numeric columns
  df[, col_numeric] #here were are not creating a new dataframe, we are applying function to dataframe, without having to convert to dataframe
  <- apply(X = df[, col_numeric], 
           MARGIN = 2,
           FUN = standardize)
  return(df)
}

#if loaded the above frunction from a source file e.g. 
# source(file = "standardization-functions.R"), all you would need is this line in the analysis script:


iris_stand <- standardize.df(df = iris)

#Debugging functions using cat---------
#n.b. can also use testthat package
#e.g. 

standardize <- function(data) {
  # Calculate mean of data
  mean_data <- mean(data)
  cat("Value of mean: ", mean_data, "\n") # cat enables you to write messages about what function should return along the way
  #\n prints the mean value on a new line
  sd_data <- sd(data)
  # Calculate standardized scores
  stand_data <- (data - mean_data)/sd_data
  return(stand_data)
}

petal_length_stand <- standardize(data = iris$Petal.Length)

#better to keep functions small to can debug easily
#eventually will want to remove cat comments (or keep tests in a seperate file to easily user friendly)