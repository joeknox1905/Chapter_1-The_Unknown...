# function for getting the residuals of the log of each integer trait against log body mass

get.residuals <- function(x) {
    LogTrait <- log(x)
    ResTrait <- lm(LogTrait ~ LogMass)$residuals
}
