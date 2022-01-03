# convert asia names from single letters to bettter descriptions
data(asia)

names(asia) <- c("asia","smoke","tub","lung","bronc","either","xray","dysp")
save(asia, file="data_asia.rda")
