# heuristics_buildbayes.R
# https://cran.r-project.org/web/packages/bnclassify/
# https://bookdown.org/robertness/causalml/docs/tutorial-probabilistic-modeling-with-bayesian-networks-and-bnlearn.html
# https://perso.univ-rennes1.fr/valerie.monbet/GM/Sachs.html

## Here we load in data sets and construct Bayesian Networks
## Then extract "rules" from them.

######################### sachs data #########################################
# We can say that Sachs et al used stimulatory and inhibitory interventions applied
# to the BN enabled them to identify causal relationships - which were interesting.
# build BN, discretize data, inference
sachs <- read.table("sachs.data.txt.gz", header = TRUE)
dsachs <-discretize(sachs, method="hartemink", breaks=3, ibreaks=60, idisc ="quantile")

isachs.modelstring <-
  paste("[PKC][PKA|PKC][Raf|PKC:PKA]",
        "[pmek|PKC:PKA:praf][p44.42|pmek:PKA]",
        "[pakts473|p44.42:PKA][P38|PKC:PKA]",
        "[pjnk|PKC:PKA][plcg][PIP3|plcg]",
        "[PIP2|Plcg:PIP3]",sep="")


isachs.modelstring <- paste("[PKC][PKA|PKC][praf|PKC:PKA]",
                "[pmek|PKC:PKA:praf][p44.42|pmek:PKA]",
                "[pakts473|p44.42:PKA][P38|PKC:PKA]",
                "[pjnk|PKC:PKA][plcg][PIP3|plcg]",
                "[PIP2|plcg:PIP3]",sep="")

isachs.modelstring <-
  paste("[PKC][PKA|PKC][Raf|PKC:PKA][Mek|PKC:PKA:Raf]",
        "[Erk|Mek:PKA][Akt|Erk:PKA][P38|PKC:PKA]",
        "[Jnk|PKC:PKA][Plcg][PIP3|Plcg][PIP2|Plcg:PIP3]",sep="")

dag.sachs <- model2network(isachs.modelstring)
isachs <- dsachs[, 1:11]

# Give names to the factor levels
for (i in names(isachs)){
  levels(isachs[, i]) <- c("LOW","AVG","HIGH")
}

sachs.fitted <- bn.fit(dag.sachs, isachs, method="bayes")
graphviz.plot(sachs.fitted)

jtree <- compile(as.grain(sachs.fitted))  # convert BN into gRain type network
jprop <- setEvidence(jtree, nodes ="p44.42", states="LOW")
querygrain(jtree, nodes="pakts473")$pakts473

querygrain(jtree, nodes="PKA")$PKA
querygrain(jprop, nodes="PKA")$PKA

bn.fit.barchart(sachs.fitted$PKA)
bn.fit.barchart(sachs.fitted$Raf)
bn.fit.barchart(sachs.fitted$Plcg)

cpt <- jtree$cptlist  # get cpt's
plist <- compileCPT(cpt)   # Most horsework done here.
sachs_cpt <- attributes(cpt)
sachs_rules <- BN2R(plist)

#--------------------------- inferencing --------------------------------------------
particles <- cpdist(sachs.fitted, nodes = "pakts473",evidence = (p44.42 == "LOW"))
prop.table(table(particles))
particles<- cpdist(fitted, nodes = "PKA",evidence = (p44.42 == "LOW"))
prop.table(table(particles))
cpquery(fitted, event = (pakts473 == "LOW") & (PKA != "HIGH"),evidence = (p44.42 == "LOW") | (praf == "LOW"))


#################### asia data set #################################
# load the data.
load("data_asia.rda")   # variable names changed by me from single letters.
# create and plot the network structure.
dag.asia <- model2network("[asia][smoke][tub|asia][lung|smoke][bronc|smoke][dysp|bronc:either][either|tub:lung][xray|either]")
asia.fitted <- bn.fit(dag.asia, asia, method="bayes")
graphviz.plot(asia.fitted)
jtree <- compile(as.grain(asia.fitted))   # convert BN into gRain type network

cpt <- jtree$cptlist  # get cpt's
plist <- compileCPT(cpt)   # Most horsework done here.
asia_cpt <- attributes(cpt)
asia_rules <- BN2R(plist)


################### Heckerman Fraud data set? ################################
# https://bookdown.org/robertness/causalml/docs/tutorial-probabilistic-modeling-with-bayesian-networks-and-bnlearn.html







