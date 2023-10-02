##############.
## GROUP 19 ##.
##############.


## MEMBERS : Tieas, Rubing, Ishika, Shinichi and Fatma ##


################
##  Project A ##
################.

# Project A : Question 1 #.

library(HDInterval)
# Define variables
n_smoker <- 3435 # The total number of smokers
x_smoker <- 171  # The number of smokers who had a stroke
n_non_smoker <- 4437 # The total number of non-smokers
x_non_smoker <- 117 # The number of non-smokers who had a stroke

# Define prior parameters with unit prior, Beta(1,1)
alpha_prior <- 1
beta_prior <- 1

# Posterior parameters, for θ+ and θ-
alpha_post_smoker <- alpha_prior + x_smoker
beta_post_smoker <- beta_prior + n_smoker - x_smoker
alpha_post_non_smoker <- alpha_prior + x_non_smoker
beta_post_non_smoker <- beta_prior + n_non_smoker - x_non_smoker

# Define density functions
prior_density <- function(x) dbeta(x, alpha_prior, beta_prior)
post_smoker_density <- function(x) dbeta(x, alpha_post_smoker, beta_post_smoker)
post_non_smoker_density <- function(x) dbeta(x, alpha_post_non_smoker, beta_post_non_smoker)

# Define summary measures for the posterior distribution
# for smoker, θ+ 
post_smoker_mean <- alpha_post_smoker / (alpha_post_smoker + beta_post_smoker)
post_smoker_median <- qbeta(0.5, alpha_post_smoker, beta_post_smoker)
post_smoker_mode <- (alpha_post_smoker - 1) / (alpha_post_smoker + beta_post_smoker - 2)
post_smoker_sd <- sqrt((alpha_post_smoker * beta_post_smoker) / ((alpha_post_smoker + beta_post_smoker)^2 * (alpha_post_smoker + beta_post_smoker + 1)))
post_smoker_ci <- qbeta(c(0.025, 0.975), alpha_post_smoker, beta_post_smoker)
post_smoker_hpd <- hdi(rbeta(1000, alpha_post_smoker, beta_post_smoker))

# for non-smoker, θ-
post_non_smoker_mean <- alpha_post_non_smoker / (alpha_post_non_smoker + beta_post_non_smoker)
post_non_smoker_median <- qbeta(0.5, alpha_post_non_smoker, beta_post_non_smoker)
post_non_smoker_mode <- (alpha_post_non_smoker - 1) / (alpha_post_non_smoker + beta_post_non_smoker - 2)
post_non_smoker_sd <- sqrt((alpha_post_non_smoker * beta_post_non_smoker) / ((alpha_post_non_smoker + beta_post_non_smoker)^2 * (alpha_post_non_smoker + beta_post_non_smoker + 1)))
post_non_smoker_ci <- qbeta(c(0.025, 0.975), alpha_post_non_smoker, beta_post_non_smoker)
post_non_smoker_hpd <- hdi(rbeta(1000, alpha_post_non_smoker, beta_post_non_smoker))

# Print summary measures
# for Smoker, θ+
cat("Smoker:\n")
cat("Mean:", post_smoker_mean, "\n")
cat("Median:", post_smoker_median, "\n")
cat("Mode:", post_smoker_mode, "\n")
cat("Standard Deviation:", post_smoker_sd, "\n")
cat("95% Credible Interval:", '[',post_smoker_ci[1], ',', post_smoker_ci[2],']', "\n")
cat("95% HPD:", '[',post_smoker_hpd[1], ',', post_smoker_hpd[2],']\n')
# for non-smoker, θ-
cat("Non-Smoker:\n")
cat("Mean:", post_non_smoker_mean, "\n")
cat("Median:", post_non_smoker_median, "\n")
cat("Mode:", post_non_smoker_mode, "\n")
cat("Standard Deviation:", post_non_smoker_sd, "\n")
cat("95% Credible Interval:", '[',post_non_smoker_ci[1],',',post_non_smoker_ci[2],']', "\n")
cat("95% HPD:", '[',post_non_smoker_hpd[1], ',', post_non_smoker_hpd[2],']\n')

# Set up the plot layout to have 1 row and 2 columns
par(mfrow = c(1, 2))
# Plot for smoker, posterior distribution of θ+ 
plot(0, 0, type = "n", xlim = c(0, 0.1), ylim = c(0, 180), 
     xlab = "Probability of Disease", ylab = "Density", 
     main = "Posterior Distribution of θ+")
curve(post_smoker_density, add = TRUE, col = "blue", lwd = 2)
abline(v = post_smoker_mean, col = "red", lwd = 2, lty = 2)

# Plot for non-smoker, posterior distribution of  θ-
plot(0, 0, type = "n", xlim = c(0, 0.1), ylim = c(0, 180), 
     xlab = "Probability of Disease ", ylab = "Density", 
     main = "Posterior Distribution of θ-")
curve(post_non_smoker_density, add = TRUE, col = "blue", lwd = 2)
abline(v = post_non_smoker_mean, col = "red", lwd = 2, lty = 2)



# Project A : Question 3 #

# Draw samples from posterior distributions
n_samples <- 1000
smoker_samples <- rbeta(n_samples, alpha_post_smoker, beta_post_smoker)
non_smoker_samples <- rbeta(n_samples, alpha_post_non_smoker, beta_post_non_smoker)

# Calculate relative risk
rel_risk_samples <- smoker_samples / non_smoker_samples

# Plot posterior distribution of relative risk
par(mfrow = c(1, 1))
hist(rel_risk_samples, breaks=100, main = "", xlab="Relative Risk",
     col = "lightblue", border = "black")
abline(v=mean(rel_risk_samples), col="red", lwd = 2, lty = 2)  # Line at RR=1

# Calculate and print summary measures
# Mode function to calculate the mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
cat("Mean Relative Risk:", mean(rel_risk_samples), "\n")
cat("Median Relative Risk:", median(rel_risk_samples), "\n")
cat("Mode Relative Risk:", Mode(rel_risk_samples), "\n")  # Mode function to calculate the mode
cat("Standard Deviation of Relative Risk:", sd(rel_risk_samples), "\n")  # Calculate standard deviation
cat("95% Credible Interval for Relative Risk:", '[',quantile(rel_risk_samples, c(0.025, 0.975)),']\n')
cat("95% HPD for Relative Risk:", '[',hdi(rel_risk_samples),']\n')



# Project A : Question 4, 5, and 6 #

library("rjags")

# Data
model.data = list(
  x1 = 171,
  n1 = 3435,
  x2 = 117,
  n2 = 4437
)

model.inits <- list(theta1 = 0.5, theta2 = 0.5)

# Model
modelString = "
model {
  # Priors
  theta1 ~ dbeta(1, 1)
  theta2 ~ dbeta(1, 1)

  # Likelihood
  x1 ~ dbin(theta1, n1)
  x2 ~ dbin(theta2, n2)

  # Relative risk
  theta_RR <- theta1 / theta2
}
"
# specify model, data, number of parallel chains
jags <- jags.model(textConnection(modelString),
                   data = model.data,
                   inits = model.inits,
                   n.chains = 2)

# Generate MCMC samples and save output for specified variables
out <- coda.samples(jags,
                    c("theta1", "theta2", "theta_RR"),
                    n.iter=10000, thin=1)

# Posterior summary statistics 
burnin <- 2000
summary(window(out,start=burnin))
# Iterations = 2000:11000
# Thinning interval = 1 
# Number of chains = 2 
# Sample size per chain = 9001 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean       SD  Naive SE Time-series SE
# theta1   0.05008 0.003714 2.768e-05      3.481e-05
# theta2   0.02658 0.002436 1.816e-05      2.364e-05
# theta_RR 1.90032 0.225936 1.684e-03      2.131e-03
# 
# 2. Quantiles for each variable:
#   
#   2.5%     25%     50%     75%   97.5%
# theta1   0.04309 0.04753 0.04996 0.05249 0.05768
# theta2   0.02198 0.02491 0.02652 0.02817 0.03163
# theta_RR 1.50402 1.74229 1.88404 2.04298 2.38514

#HPD
library(runjags)
out.combined <- combine.mcmc(out)
HPDinterval(out.combined)

# lower      upper
# theta1   0.04306756 0.05772057
# theta2   0.02189134 0.03148368
# theta_RR 1.47313610 2.34421193
# attr(,"Probability")
# [1] 0.95

# Produce general summary of obtained MCMC sampling, 
# History plot & posterior distributions & autocorrelation plot
print(out,digits=3)
plot(out, trace=TRUE, density = TRUE)   
plot(window(out,start=burnin), trace=TRUE, density = TRUE)  

# Convergence tests
out.mcmc <- as.mcmc.list(out)
gelman.diag(out.mcmc)
gelman.plot(out.mcmc,ask=FALSE)

geweke.diag(out.mcmc)
geweke.plot(out.mcmc,ask=FALSE)



# Project A : Question 7 #

# Model
modelString2 = "
model {
  # Priors
  theta1 ~ dbeta(1, 1)
  theta2 ~ dbeta(1, 1)

  # Likelihood
  x1 ~ dbin(theta1, n1)
  x2 ~ dbin(theta2, n2)

  # Relative risk
  theta_AR <- 1 - theta2 / theta1
}
"

# specify model, data, number of parallel chains
jags2 <- jags.model(textConnection(modelString2),
                    data = model.data,
                    inits = model.inits,
                    n.chains = length(model.inits))

# Generate MCMC samples and save output for specified variables
out2 <- coda.samples(jags2,
                     c("theta1", "theta2", "theta_AR"),
                     n.iter=10000, thin=1)

# Posterior summary statistics
burnin2 <- 2000
summary(window(out2,start=burnin2))
# Iterations = 2000:11000
# Thinning interval = 1 
# Number of chains = 2 
# Sample size per chain = 9001 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean       SD  Naive SE Time-series SE
# theta1   0.04997 0.003667 2.733e-05      3.316e-05
# theta2   0.02662 0.002400 1.789e-05      2.288e-05
# theta_AR 0.46438 0.062436 4.653e-04      5.819e-04
# 
# 2. Quantiles for each variable:
#   
#   2.5%     25%     50%     75%  97.5%
# theta1   0.04294 0.04749 0.04989 0.05244 0.0573
# theta2   0.02219 0.02495 0.02653 0.02819 0.0315
# theta_AR 0.33122 0.42544 0.46740 0.50772 0.5763

# History plot & posterior distributions & autocorrelation plot
plot(out2, trace=TRUE, density = TRUE)   
plot(window(out2,start=burnin2), trace=TRUE, density = TRUE) 

#HPD
out2.combined <- combine.mcmc(out2)
HPDinterval(out2.combined)
# lower      upper
# theta1   0.04292285 0.05727797
# theta2   0.02200795 0.03132428
# theta_AR 0.34114014 0.58434425
# attr(,"Probability")
# [1] 0.95


# Convergence tests
out2.mcmc <- as.mcmc.list(out2)
gelman.diag(out2.mcmc)
gelman.plot(out2.mcmc,ask=FALSE)

geweke.diag(out2.mcmc)
geweke.plot(out2.mcmc,ask=FALSE)


################
##  Project B ##
################.

#The packages that we use
library("nimble")
library("coda")
library("readr")

# Project B : Question 1 #

N <- c(5946, 980, 2074, 6101)
z <- c(1319, 134, 241, 193)
model.dataZ <- list('z' = z)
model.dataN <- list('N' = N)

# Set 2 different initials for 2 chains
init1 <- list(p = 0) # 'extreme' choice for prevalence
init2 <- list(p = 1) # 'extreme' choice for prevalence
model.initials <- list(init1, init2)
#other ways:
# model.initials2<- list(list(p = a),
#                        list(p = b)
#                        )

# Model specification
ModelPrevalence <- nimbleCode(
  {
    # likelihood specification : binomial
    for (i in 1:length(N)) {
      z[i] ~ dbinom(p[i], N[i])
    }
    
    # prior information : vague one
    for (i in 1:length(N)) {
      p[i] ~  dunif(0, 1)
    }
    
    # derived quantity : nothing
    # .... can put something here if we have it;
    # like a link function or etc
    
  })

# Run 2 chains using 2 different sets of initial; 5000 iterationns (no burn-in yet)
OutputMCMC_q1 <- nimbleMCMC(
  code=ModelPrevalence,
  data=model.dataZ,
  constants=model.dataN,
  inits=model.initials,
  monitors=c("p"),
  niter=5000,
  nchains=2,
  #nburnin=0,
  setSeed=2023,
  summary=TRUE
)

#Getting the chains
chain1_q1 <- OutputMCMC_q1$samples$chain1
chain2_q1 <- OutputMCMC_q1$samples$chain2

# History plots: p1 (each chain)
traceplot(as.mcmc(chain1_q1[,'p[1]']),
          col = "turquoise",
          main = 'Trace Plot of Prevalence for New York City',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q1[,'p[1]']),
          add = TRUE,
          col = "pink")
#to check if chain1 and chain2 produce the same result: sum(chain1_q1 == chain2_q1)

# History plots: p2 (each chain)
traceplot(as.mcmc(chain1_q1[,'p[2]']),
          col = "turquoise",
          main = 'Trace Plot of Prevalence for Westchester/Rockland Counties',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q1[,'p[2]']),
          add = TRUE,
          col = "pink")

# History plots: p3 (each chain)
traceplot(as.mcmc(chain1_q1[,'p[3]']),
          col = "turquoise",
          main = 'Trace Plot of Prevalence for Long Island',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q1[,'p[3]']),
          add = TRUE,
          col = "pink")

# History plots: p4 (each chain)
traceplot(as.mcmc(chain1_q1[,'p[4]']),
          col = "turquoise",
          main = 'Trace Plot of Prevalence for Rest of NYS',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q1[,'p[4]']),
          add = TRUE,
          col = "pink")

# Autocorrelation plots: p1 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q1[,'p[1]'], chain2_q1[,'p[1]'])),
  main = 'Autocorrelation of Prevalence for New York City')

# Autocorrelation plots: p2 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q1[,'p[2]'], chain2_q1[,'p[2]'])),
  main = 'Autocorrelation of Prevalence for for Westchester/Rockland Counties')

# Autocorrelation plots: p3 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q1[,'p[3]'], chain2_q1[,'p[3]'])),
  main = 'Autocorrelation of Prevalence for Long Island')

# Autocorrelation plots: p4 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q1[,'p[4]'], chain2_q1[,'p[4]'])),
  main = 'Autocorrelation of Prevalence for Rest of NYS')

#Gelman and Rubin's convergence diagnostics updated model
chain.combined_q1 <- mcmc.list(
  as.mcmc(chain1_q1[, c(1,2,3,4)]),
  as.mcmc(chain2_q1[, c(1,2,3,4)])
)
gelman.diag(chain.combined_q1)
gelman.plot(chain.combined_q1)



# Project B : Question 2 #

# Run MCMC and set the burn-in and thinning, in the end, save 5000 iterations
OutputMCMC_q2 <- nimbleMCMC(
  code=ModelPrevalence,
  data=model.dataZ,
  constants=model.dataN,
  inits=model.initials,
  monitors=c("p"),
  nchains=2,
  nburnin=1000,
  niter=16000,
  thin = 3,
  setSeed=2023,
  summary=TRUE
)

#Getting the chains
chain1_q2 <- OutputMCMC_q2$samples$chain1
chain2_q2 <- OutputMCMC_q2$samples$chain2

# History plots: p1 (each chain)
traceplot(as.mcmc(chain1_q2[,'p[1]']),
          col = "turquoise",
          main = 'Trace Plot of Prevalence for New York City',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q2[,'p[1]']),
          add = TRUE,
          col = "pink")
#to check if chain1 and chain2 produce the same result: sum(chain1_q2 == chain2_q2)

# History plots: p2 (each chain)
traceplot(as.mcmc(chain1_q2[,'p[2]']),
          col = "turquoise",
          main = 'Trace Plot of Prevalence for Westchester/Rockland Counties',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q2[,'p[2]']),
          add = TRUE,
          col = "pink")

# History plots: p3 (each chain)
traceplot(as.mcmc(chain1_q2[,'p[3]']),
          col = "turquoise",
          main = 'Trace Plot of Prevalence for Long Island',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q2[,'p[3]']),
          add = TRUE,
          col = "pink")

# History plots: p4 (each chain)
traceplot(as.mcmc(chain1_q2[,'p[4]']),
          col = "turquoise",
          main = 'Trace Plot of Prevalence for Rest of NYS',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q2[,'p[4]']),
          add = TRUE,
          col = "pink")

# Autocorrelation plots: p1 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q2[,'p[1]'], chain2_q2[,'p[1]'])),
  main = 'Autocorrelation of Prevalence for New York City')

# Autocorrelation plots: p2 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q2[,'p[2]'], chain2_q2[,'p[2]'])),
  main = 'Autocorrelation of Prevalence for for Westchester/Rockland Counties')

# Autocorrelation plots: p3 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q2[,'p[3]'], chain2_q2[,'p[3]'])),
  main = 'Autocorrelation of Prevalence for Long Island')

# Autocorrelation plots: p4 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q2[,'p[4]'], chain2_q2[,'p[4]'])),
  main = 'Autocorrelation of Prevalence for Rest of NYS')

#Gelman and Rubin's convergence diagnostics updated model
chain.combined_q2 <- mcmc.list(
  as.mcmc(chain1_q2[, c(1,2,3,4)]),
  as.mcmc(chain2_q2[, c(1,2,3,4)])
)
gelman.diag(chain.combined_q2)
gelman.plot(chain.combined_q2)



# Project B : Question 3 #

#TESTING DIFFERENCE BETWEEN REGIONS
# Get summary statistics for p
(p1.mean <- OutputMCMC_q2$summary$all.chains['p[1]','Mean'])
(p2.mean <- OutputMCMC_q2$summary$all.chains['p[2]','Mean'])
(p3.mean <- OutputMCMC_q2$summary$all.chains['p[3]','Mean'])
(p4.mean <- OutputMCMC_q2$summary$all.chains['p[4]','Mean'])

# Preparation: Get MCMC samples for each p
p1.sample <- c(chain1_q2[,'p[1]'],
               chain2_q2[,'p[1]'])
p2.sample <- c(chain1_q2[,'p[2]'],
               chain2_q2[,'p[2]'])
p3.sample <- c(chain1_q2[,'p[3]'],
               chain2_q2[,'p[3]'])
p4.sample <- c(chain1_q2[,'p[4]'],
               chain2_q2[,'p[4]'])

# Test whether the difference: New York VS Westchester
cat("\n Test whether the difference: New York VS Westchester \n")
d12.sample <- p1.sample - p2.sample 
quantile(d12.sample, 
         probs=c(2.5, 97.5)/100) # Equal tail interval
HPDinterval(mcmc(data=d12.sample),
            alpha = 0.05) # HPD interval

# Test whether the difference: New York VS Long Island
cat("\n Test whether the difference: New York VS Long Island \n")
d13.sample <- p1.sample - p3.sample 
quantile(d13.sample, 
         probs=c(2.5, 97.5)/100) # Equal tail interval
HPDinterval(mcmc(data=d13.sample),
            alpha = 0.05) # HPD interval

# Test whether the difference: New York VS Rest of NYS
cat("\n Test whether the difference:  New York VS Rest of NYS \n")
d14.sample <- p1.sample - p4.sample 
quantile(d14.sample, 
         probs=c(2.5, 97.5)/100) # Equal tail interval
HPDinterval(mcmc(data=d14.sample),
            alpha = 0.05) # HPD interval

# Test whether the difference: Westchester VS Long Island
cat("\n Test whether the difference:  Westchester VS Long Island \n")
d23.sample <- p2.sample - p3.sample 
quantile(d23.sample, 
         probs=c(2.5, 97.5)/100) # Equal tail interval
HPDinterval(mcmc(data=d23.sample),
            alpha = 0.05) # HPD interval

# Test whether the difference: Westchester VS Rest of NYS
cat("\n Test whether the difference:  Westchester VS Rest of NYS \n")
d24.sample <- p2.sample - p4.sample 
quantile(d24.sample, 
         probs=c(2.5, 97.5)/100) # Equal tail interval
HPDinterval(mcmc(data=d24.sample),
            alpha = 0.05) # HPD interval

# Test whether the difference: Long Island VS Rest of NYS
cat("\n Test whether the difference:  Long Island VS Rest of NYS \n")
d34.sample <- p3.sample - p4.sample 
quantile(d34.sample, 
         probs=c(2.5, 97.5)/100) # Equal tail interval
HPDinterval(mcmc(data=d34.sample),
            alpha = 0.05) # HPD interval



# Project B : Question 4 #

# Calculate the MODE (peak) for the Beta distribution
beta.mode<- function(alpha, beta) {
  return((alpha-1) / (alpha+beta-2))}

#sensitivity  
beta.mode(205,29) #result:0.8793103
#specificity =
beta.mode(288,2) #result: 0.9965278

# Create plots of Beta distribution
# Define range
p = seq(0, 1, length=1000)

# Plot for sensitivity:
par(mfrow=c(1,1))
plot(p, dbeta(p, 205 , 29), type='l', ylab="Density", xlab="Sensitivity",
     main="Prior distribution for Sensitivity ~ Beta(205, 29)",
     xlim=c(0.80,1), ylim=c(-0.5,30))
points(0.879, 0.5, col="red", pch=19) #the mode or point estimate
abline(v=0.879, col="purple", lty =2) #the mode or point estimate
text(0.879, 3, "0.879", col="purple")

# Plot for specificity:
par(mfrow=c(1,1))
plot(p, dbeta(p, 288, 2), type='l', ylab="Density", xlab="Specificity",
     main="Prior distribution for Specificity ~ Beta(288, 2)",
     xlim=c(0.80,1), ylim=c(-0.5,120))
points(0.996, 0.5, col="red", pch=19) #the mode or point estimate
abline(v=0.996, col="purple", lty =2) #the mode or point estimate
text(0.996, 6, "0.996", col="purple")



# Project B : Question 5 #

# Set 2 different initials for 2 chains including for se and sp
init1_q5 <- list(p = 0, se=0, sp=0)
init2_q5 <- list(p = 1, se=1, sp=1) # 'extreme' choice for prevalence
model.initials2 <- list(init1_q5, init2_q5)

# Model specification
ModelCumIncidence <- nimbleCode(
  {
    # likelihood specification : binomial
    for (i in 1:length(N)) {
      z[i] ~ dbinom(p[i], N[i])
    }
    
    # derived quantity : pi
    for (i in 1:length(N)) {
      pi[i] <- ( p[i] + sp - 1 ) / (se + sp -1)
    }
    
    # prior information : vague one for the prevalence
    for (i in 1:length(N)) {
      p[i] ~  dunif(0, 1)
    }
    se ~ dbeta(205 , 29)
    sp ~ dbeta(288 , 2)
    
  })

# Run 2 chains using 2 different sets of initial; 5000 iterations in the end
OutputMCMC_q5 <- nimbleMCMC(
  code=ModelCumIncidence,
  data=model.dataZ,
  constants=model.dataN,
  inits=model.initials2,
  monitors=c("p","pi"),
  nchains=2,
  nburnin=1000,
  niter=16000,
  thin = 3,
  setSeed=2023,
  summary=TRUE
)

#Getting the chains
chain1_q5 <- OutputMCMC_q5$samples$chain1
chain2_q5 <- OutputMCMC_q5$samples$chain2

# History plots: pi1 (each chain)
traceplot(as.mcmc(chain1_q5[,'pi[1]']),
          col = "turquoise",
          main = 'Trace Plot of Cummulative Incidence for New York City',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q5[,'pi[1]']),
          add = TRUE,
          col = "pink")

# History plots: pi2 (each chain)
traceplot(as.mcmc(chain1_q5[,'pi[2]']),
          col = "turquoise",
          main = 'Trace Plot of Cummulative Incidence for Westchester/Rockland Counties',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q5[,'pi[2]']),
          add = TRUE,
          col = "pink")

# History plots: pi3 (each chain)
traceplot(as.mcmc(chain1_q5[,'p[3]']),
          col = "turquoise",
          main = 'Trace Plot of Cummulative Incidence for Long Island',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q5[,'p[3]']),
          add = TRUE,
          col = "pink")

# History plots: pi4 (each chain)
traceplot(as.mcmc(chain1_q5[,'p[4]']),
          col = "turquoise",
          main = 'Trace Plot of Cummulative Incidence for Rest of NYS',
          ylab = 'Value')
traceplot(as.mcmc(chain2_q5[,'p[4]']),
          add = TRUE,
          col = "pink")

# Autocorrelation plots: pi1 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q5[,'pi[1]'], chain2_q5[,'pi[1]'])),
  main = 'Autocorrelation of Cummulative Incidence for New York City')

# Autocorrelation plots: pi2 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q5[,'pi[2]'], chain2_q5[,'pi[2]'])),
  main = 'Autocorrelation of Cummulative Incidence for Westchester/Rockland Counties')

# Autocorrelation plots: pi3 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q5[,'pi[3]'], chain2_q5[,'pi[3]'])),
  main = 'Autocorrelation of Prevalence for Long Island')

# Autocorrelation plots: pi4 (each chain)
autocorr.plot(as.mcmc(
  c(chain1_q5[,'pi[4]'], chain2_q5[,'pi[4]'])),
  main = 'Autocorrelation of Cummulative Incidence for Rest of NYS')


#Gelman and Rubin's convergence diagnostics updated model
chain.combined_q5 <- mcmc.list(
  as.mcmc(
    chain1_q5[,c("pi[1]","pi[2]","pi[3]","pi[4]")]),
  as.mcmc(
    chain2_q5[,c("pi[1]","pi[2]","pi[3]","pi[4]")])
)
gelman.diag(chain.combined_q5)
gelman.plot(chain.combined_q5)

# Here comes the conclusion: region-based pi
(pi1.mean <- OutputMCMC_q5$summary$all.chains['pi[1]','Mean'])
(pi2.mean <- OutputMCMC_q5$summary$all.chains['pi[2]','Mean'])
(pi3.mean <- OutputMCMC_q5$summary$all.chains['pi[3]','Mean'])
(pi4.mean <- OutputMCMC_q5$summary$all.chains['pi[4]','Mean'])



