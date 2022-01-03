# heuristics_functions.R


##################################### functions #################################
# Convert bnlearn network into gRain network, Convert a gRain CPT into rules. 
# Even a simple net like asia will generate a lot of rules (36 in fact):
# P( asia = no)   0.99
# P( asia = yes)  0.008

# The cpt table is a bit fiddly and to get the rules organised takes a bit conversion.
BN2R <- function(cpt_list){
  numcpts <- length(cpt_list)
  rulelist <- data.frame(rulelist=character(0), stringsAsFactors = FALSE)

  for(i in 1:length(cpt_list)){           # Pass makerule() one cpt at a time
    templist <- makerule(cpt_list[[i]])   # Most work is performed in makerule()
    rulelist <- rbind(rulelist,templist)  # Join rules together
  }
  
  rulelist <- as.data.frame(rulelist,stringsAsFactors = FALSE)
  print(rulelist,right=FALSE)
  return(rulelist)

}

# Receives a CPT (conditional variables and their levels to build a rule), returns 
# a pasted string consisting of one to several rules based on the complexity of the supplied CPT.
makerule <- function(cpt_item){
  prob_list <- list()
  templist <- ""
  rulelist <- "Empty"
  rulecount <- 1
  vcount <- 1
  
  probvars <- attr(cpt_item,"dimnames")  # get all variable names from cpt
  nvars <- length(probvars)   # how many variables in cpt?
  probdata <- as.data.frame(cpt_item)
  problevels <- row.names.data.frame(probdata)
  probdata <- unlist(probdata)
  
  for(i in 1:nvars){prob_list[i] <- list(probvars[[i]]) }
  all <- expand.grid(prob_list, stringsAsFactors = FALSE) 
  names(all) <- names(probvars)   # get the var names and then..
  
  npvars <- length(probvars)
  all <- all[names(probvars[npvars:1])]  # ..change order of vars
 
  nrules <- dim(all)
  for(i in 1:nrules[1]){   # get every rule given the combinations of levels
    for(j in 1:nvars){    # get every variable in each rule
      templist <- paste0(templist,names(all[j]),"=",all[i,j],",")
      if(j==1 && nvars > 1){templist <- stringr::str_replace(templist, ",", ""); templist <- paste0(templist,"|")}
      if(j==1 && nvars ==1){templist <- stringr::str_replace(templist, ",", ""); templist <- paste0(templist," ")}
      if(j==nvars){templist <- paste0(templist,"   ",round(probdata[i],3))}
    } 
    rulelist[rulecount] <- templist 
    templist <- ""
    rulecount <- rulecount+1
  } 
  
  rulelist <- as.data.frame(rulelist,stringsAsFactors =FALSE)
  rulelist <- stringi::stri_c(rulelist[,1], ")")  # add ")" at end of string
  rulelist <- paste0("P(", rulelist)  # add at start of string
  rulelist <- as.data.frame(rulelist,stringsAsFactors = FALSE)  # causes problems with rbind() if not declared as dataframe
  return(rulelist)
}

# some info for extracting attributes from data structures
#rules <- attr(CPT,"dag")
#rules <- rules@edgeData
#rules <- names(rules)


# functions from Avril Cochlan
calcLikelihoodForProportion <- function(successes, total)
{  curve(dbinom(successes,total,x))} # plot the likelihood}

calcPosteriorForProportion <- function(successes, total, a, b)
{
  # Adapted from triplot() in the LearnBayes package
  # Plot the prior, likelihood and posterior:
  likelihood_a = successes + 1; likelihood_b = total - successes + 1
  posterior_a = a + successes;  posterior_b = b + total - successes
  theta = seq(0.005, 0.995, length = 500)
  prior = dbeta(theta, a, b)
  likelihood = dbeta(theta, likelihood_a, likelihood_b)
  posterior  = dbeta(theta, posterior_a, posterior_b)
  m = max(c(prior, likelihood, posterior))
  plot(theta, posterior, type = "l", ylab = "Density", lty = 2, lwd = 3,
       main = paste("beta(", a, ",", b, ") prior, B(", total, ",", successes, ") data,",
                    "beta(", posterior_a, ",", posterior_b, ") posterior"), ylim = c(0, m), col = "red")
  lines(theta, likelihood, lty = 1, lwd = 3, col = "blue")
  lines(theta, prior, lty = 3, lwd = 3, col = "green")
  legend(x=0.1,y=m, c("Prior", "Likelihood", "Posterior"), lty = c(3, 1, 2),
         lwd = c(3, 3, 3), col = c("green", "blue", "red"))
  # Print out summary statistics for the prior, likelihood and posterior:
  calcBetaMode <- function(aa, bb) { BetaMode <- (aa - 1)/(aa + bb - 2); return(BetaMode); }
  calcBetaMean <- function(aa, bb) { BetaMean <- (aa)/(aa + bb); return(BetaMean); }
  calcBetaSd   <- function(aa, bb) { BetaSd <- sqrt((aa * bb)/(((aa + bb)^2) * (aa + bb + 1))); return(BetaSd); }
  prior_mode      <- calcBetaMode(a, b)
  likelihood_mode <- calcBetaMode(likelihood_a, likelihood_b)
  posterior_mode  <- calcBetaMode(posterior_a, posterior_b)
  prior_mean      <- calcBetaMean(a, b)
  likelihood_mean <- calcBetaMean(likelihood_a, likelihood_b)
  posterior_mean  <- calcBetaMean(posterior_a, posterior_b)
  prior_sd        <- calcBetaSd(a, b)
  likelihood_sd   <- calcBetaSd(likelihood_a, likelihood_b)
  posterior_sd    <- calcBetaSd(posterior_a, posterior_b)
  print(paste("mode for prior=",prior_mode,", for likelihood=",likelihood_mode,", for posterior=",posterior_mode))
  print(paste("mean for prior=",prior_mean,", for likelihood=",likelihood_mean,", for posterior=",posterior_mean))
  print(paste("sd for prior=",prior_sd,", for likelihood=",likelihood_sd,", for posterior=",posterior_sd))
}

findBeta <- function(quantile1,quantile2,quantile3)
{
  # find the quantiles specified by quantile1 and quantile2 and quantile3
  quantile1_p <- quantile1[[1]]; quantile1_q <- quantile1[[2]]
  quantile2_p <- quantile2[[1]]; quantile2_q <- quantile2[[2]]
  quantile3_p <- quantile3[[1]]; quantile3_q <- quantile3[[2]]
  
  # find the beta prior using quantile1 and quantile2
  priorA <- beta.select(quantile1,quantile2)
  priorA_a <- priorA[1]; priorA_b <- priorA[2]
  
  # find the beta prior using quantile1 and quantile3
  priorB <- beta.select(quantile1,quantile3)
  priorB_a <- priorB[1]; priorB_b <- priorB[2]
  
  # find the best possible beta prior
  diff_a <- abs(priorA_a - priorB_a); diff_b <- abs(priorB_b - priorB_b)
  step_a <- diff_a / 100; step_b <- diff_b / 100
  if (priorA_a < priorB_a) { start_a <- priorA_a; end_a <- priorB_a }
  else                     { start_a <- priorB_a; end_a <- priorA_a }
  if (priorA_b < priorB_b) { start_b <- priorA_b; end_b <- priorB_b }
  else                     { start_b <- priorB_b; end_b <- priorA_b }
  steps_a <- seq(from=start_a, to=end_a, length.out=1000)
  steps_b <- seq(from=start_b, to=end_b, length.out=1000)
  max_error <- 10000000000000000000
  best_a <- 0; best_b <- 0
  for (a in steps_a)
  {
    for (b in steps_b)
    {
      # priorC is beta(a,b)
      # find the quantile1_q, quantile2_q, quantile3_q quantiles of priorC:
      priorC_q1 <- qbeta(c(quantile1_p), a, b)
      priorC_q2 <- qbeta(c(quantile2_p), a, b)
      priorC_q3 <- qbeta(c(quantile3_p), a, b)
      priorC_error <- abs(priorC_q1-quantile1_q) +
        abs(priorC_q2-quantile2_q) +
        abs(priorC_q3-quantile3_q)
      if (priorC_error < max_error)
      {
        max_error <- priorC_error; best_a <- a; best_b <- b
      }
    }
  }
  print(paste("The best beta prior has a=",best_a,"b=",best_b))
}


