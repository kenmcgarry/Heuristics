# heuristics.R
# https://cran.r-project.org/web/packages/heuristica/vignettes/
# https://cran.r-project.org/web/packages/heuristica/

library(heuristica)


schools <- data.frame(Name=c("Heathfield", "Breckenbeds", "Harlow Green", "Prior Street", "HomeSchooled"), 
                      Dropout_Rate=c(25.5, 11.8, 28.7, 21.6, 4.5), 
                      Low_Income_Students=c(82.5, 88.8, 63.2, 84.5, 30.3), 
                      Limited_English_Students=c(11.4, 0.1, 0, 28.3, 0.1))
schools


# fit two models to the high school data
# ttbModel, Take The Best, which uses the highest-validity cue that discriminates 
# (more details below).
# regModel, a version of R's "lm" function for linear regression wrapped to fit 
# into heurstica's interface.

criterion_col <- 2
ttb <- ttbModel(schools, criterion_col, c(3:4))
reg <- regModel(schools, criterion_col, c(3:4))


