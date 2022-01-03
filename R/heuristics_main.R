# heuristics.R
# https://machinelearningmastery.com/introduction-to-bayesian-belief-networks/
# https://cran.r-project.org/web/packages/heuristica/vignettes/

library(heuristica)
library(eba)
library(Rgraphviz)
library(bnlearn)
library(gRain)
library(gtools)
library(LaplacesDemon)
library(LearnBayes)

setwd("C:/common_laptop/R-files/heuristics")

source("heuristics_functions.R")
source("heuristics_buildbayes.R")
#source("heuristics_applyheuristics.R")
#source("heuristics_others.R")  # compare/contrast with other methods
#source("heuristics_interest.R")



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rgraphviz")


