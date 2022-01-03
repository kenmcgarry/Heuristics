# heuristics_buildbayes.R
# https://cran.r-project.org/web/packages/bnclassify/
# https://bookdown.org/robertness/causalml/docs/tutorial-probabilistic-modeling-with-bayesian-networks-and-bnlearn.html
# https://perso.univ-rennes1.fr/valerie.monbet/GM/Sachs.html


# some example dags
dag1 <- model2network("[A][B|A][C|A]")
dag2 <- model2network("[A|B:C][B][C]")
dag3 <- model2network("[A][B|A][C|A][D]")
dag4 <- model2network("[A][B|A][F|A]")
dag5 <- set.arc(dag1, from = "B", to = "C")
dag6 <- drop.arc(dag5, from = "A", to = "C")
par(mfrow = c(3, 2))
graphviz.compare(dag1, dag2, dag5, dag6)

# The style of the graphical comparison can be further modified using the arguments that 
# graphviz.compare() shares with graphviz.plot(), such as shape, layout, main and sub.
par(mfrow = c(2, 2))
graphviz.compare(dag1, dag2, dag5, dag6, shape = "rectangle", layout = "fdp",
                   main = c("ONE", "TWO", "THREE", "FOUR"),
                   sub = paste("SHD =", c("0", shd(dag1, dag2), shd(dag1, dag5), shd(dag1, dag6))))


# Other optional arguments specific to diff = "diff-from-first" can be passed as a list via 
# the argument diff.args. They are:
#
# tp.col, tp.lty, tp.lwd are the colour, line type and line width of true positive arcs;
# fp.col, fp.lty, fp.lwd are the colour, line type and line width of false positive arcs;
# fn.col, fn.lty, fn.lwd are the colour, line type and line width of negative arcs.

par(mfrow = c(2, 2))
graphviz.compare(dag1, dag2, dag5, dag6, shape = "rectangle", layout = "fdp",
                   main = c("ONE", "TWO", "THREE", "FOUR"),
                   sub = paste("SHD =", c("0", shd(dag1, dag2), shd(dag1, dag5), shd(dag1, dag6))),
                   diff.args = list(tp.lwd = 2, tp.col = "green", fn.col = "orange"))

# naive Bayes

data(learning.test)
# this is an in-sample prediction with naive Bayes (parameter learning
# is performed implicitly during the prediction).
bn = naive.bayes(learning.test, "A")
pred = predict(bn, learning.test)
table(pred, learning.test[, "A"])

# another example
data(learning.test)
fitted = bn.fit(hc(learning.test), learning.test)
# the result should be around 0.025.
cpquery(fitted, (B == "b"), (A == "a"))
# programmatically build a conditional probability query...
var = names(learning.test)
obs = 2
str = paste("(", names(learning.test)[-3], " == '",
            sapply(learning.test[obs, -3], as.character), "')",
            sep = "", collapse = " & ")
str
str2 = paste("(", names(learning.test)[3], " == '",
             as.character(learning.test[obs, 3]), "')", sep = "")
str2

cmd = paste("cpquery(fitted, ", str2, ", ", str, ")", sep = "")
eval(parse(text = cmd))
# ... but note that predict works better in this particular case.

attr(predict(fitted, "C", learning.test[obs, -3], prob = TRUE), "prob")
# do the same with likelihood weighting.
cpquery(fitted, event = eval(parse(text = str2)),
        evidence = as.list(learning.test[2, -3]), method = "lw")
attr(predict(fitted, "C", learning.test[obs, -3],
             method = "bayes-lw", prob = TRUE), "prob")
# conditional distribution of A given C == "c".
table(cpdist(fitted, "A", (C == "c")))

dmarks <- discretize(marks, breaks=3, method="interval")

bn <- hc(dmarks)
fitted <- bn.fit(bn, data=dmarks)
plot(bn)

fitted$MECH  # display CPT's
fitted$ALG
fitted$STAT

jtree <- compile(as.grain(fitted))

jprop <- setFinding(jtree, nodes = "MECH")

querygrain(jtree, nodes = "STAT")

######################### sachs data #########################################
# We can say that Sachs et al used stimulatory and inhibitory interventions applied
# to the BN enabled them to identify causal relationships - which were interesting.
# build BN, discretize data, inference
sachs <- read.table("sachs.data.txt.gz", header = TRUE)
dsachs <-discretize(sachs, method="hartemink", breaks=3, ibreaks=60, idisc ="quantile")

#isachs.modelstring <-
#  paste("[PKC][PKA|PKC][Raf|PKC:PKA]",
#        "[pmek|PKC:PKA:praf][p44.42|pmek:PKA]",
#        "[pakts473|p44.42:PKA][P38|PKC:PKA]",
#        "[pjnk|PKC:PKA][plcg][PIP3|plcg]",
#        "[PIP2|Plcg:PIP3]",sep="")


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

jtree <- compile(as.grain(sachs.fitted))
jprop <- setFinding(jtree, nodes ="p44.42", states="LOW")
querygrain(jtree, nodes="pakts473")$pakts473
querygrain(jtree, nodes="PKA")$PKA

bn.fit.barchart(sachs.fitted$PKA)
bn.fit.barchart(sachs.fitted$Raf)
bn.fit.barchart(sachs.fitted$Plcg)

result <- BN2R(sachs.fitted)

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

myrules <- BN2R(plist)

