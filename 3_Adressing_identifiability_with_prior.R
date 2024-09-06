## 3. Using an informative prior to address identifiability issues ----

library(nimble) ; library(MCMCvis)
library(tidyverse)
library(cowplot)

### a. Initialize parameters to simulate data ----
N <- 30  ; J <- 36 ; psi <- 0.8
p <- 0.5 ; w_a <- 0.9 ; w_b <- 0.7 

detect_simulate <- function(N, J, psi, p){
  1:N %>%
    map_df(function(i) {
      Z <- rbinom(1, 1, prob=psi)
      Y <- rbinom(J ,1, prob=Z*p)
      return(c(Z=Z,Y=Y))
    })
}

### c.UNIF (0,1) ----
#### i. Building model
occu_unif01 <- nimbleCode({
  
  # Likelihood
  for(i in 1:nsites) {
    # state model
    z[i] ~ dbern(psi)
    # obs model
    p.mix[i] <- z[i] * p 
    
    for(j in 1:nvisits){
      y[i,j] ~ dbern(p.mix[i])
      w.mix[i,j] <- (y[i,j] * wa) + ((1-y[i,j]) * (1-wb))
      w[i,j] ~ dbern(w.mix[i,j])
    }
  }
  
  # priors
  psi ~ dunif(0,1)
  p ~ dunif(0,1)
  wa ~ dunif(0,1)
  wb ~ dunif(0,1)
})

# specify parameters that need to be estimated
# parameters <- c("psi","p", "wa", "wb")
my.constants <- list(nsites = N, nvisits = J)

# pick initial values
initial.values <- function() list(psi = runif(1,0,1),
                                        p = runif(1,0,1),
                                        wa = runif(1,0,1),
                                        wb = runif(1,0,1),
                                        z = as.numeric(apply(Y,1,sum)>0))
# specify MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

#### ii. Run simulations ----
# simulations runs
k = 100

# final_results <- array(dim=c(k,8,4))
unif01 <- matrix(nrow=k, ncol=8)

psi_exp <- vector() ; p_exp <- vector() ; wa_exp <- vector() ; wb_exp <- vector()

set.seed(123) # run before running simulations for each prior
for(i in 1:k){
  Y <- detect_simulate(N,J,psi,p)
  W <- apply(Y[,-1], c(1,2), function(x){
    if(x==0){
      W <- rbinom(n=1, size=1, prob=(1-w_b))
    } else{
      W <- rbinom(n=1, size=1, prob=w_a)
    }
    return(W)
  })
  sumI <- apply(W, 1,sum)
  
  # psi_exp[i] <- sum(Y[,1])/N
  # p_exp[i] <- sum(rowSums(Y[,-1]))/(J*sum(Y[,1]))
  # wa_exp[i] <- sum(W[which(Y[,-1]==1)])/sum(Y[,-1])
  # wb_exp[i] <- length(which(W[which(Y[,-1]==0)]==0)) / length(which(Y[,-1]==0))
  # }
  # #combine the expected parameter according to each dataset (with seeds(123))
  # exp <- cbind(psi_exp, p_exp, wa_exp, wb_exp)
  
  # read in data
  my.data <- list(w = W)
  
  # pick initial values
  # init.values <- initial.values[]
  
  # create model as an R object (uncompiled model)
  res<- nimbleModel(code = occu_unif01, 
                    data = my.data,
                    constants = my.constants,
                    inits = initial.values())
  
  # create a MCMC configuration
  compileNimble(res)
  resConf <- configureMCMC(res)
  resConf$removeSamplers(c('p','wa'))
  resConf$addSampler(target = c('p','wa'),
                     type = 'AF_slice')
  
  # create a MCMC function and compile it
  res_MCMC <- buildMCMC(resConf)
  Cres_MCMC <- compileNimble(res_MCMC, project = res)
  
  # run NIMBLE
  samples <- runMCMC(mcmc = Cres_MCMC,
                     niter = n.iter,
                     nburnin = n.burnin,
                     nchain = n.chains)
  
  smry <- MCMCsummary(object = samples, round = 2, Rhat=F, n.eff=F)
  # ith simu, 8 parameter values (mean and median), jth prior
  unif01[i,1:4] <- smry$mean   # final_results[i,1:4,j]
  unif01[i,5:8] <- smry$`50%`  # final_results[i,5:8,j]
  # }
  
  # follow the advancement
  print(i)
}

#combine the expected parameter according to each dataset (with seeds(123))
exp <- cbind(psi_exp, p_exp, wa_exp, wb_exp)

# save the results of the simulations for each priors
write.csv(unif01, file="Using Informative priors/bayes_NIMBLE_50-100_unif01.csv")

#### iii. Check chains convergence (for the last simulation)
MCMCtrace(samples, pdf=F, ind=T, Rhat=T, n.eff=T, params=c("psi","p"))

par(mfrow=c(2,1))
hist(c(samples$chain1[,"wa"], samples$chain2[,"wa"]), main="wA", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=wa_exp[k], col="blue", lwd=3)
hist(c(samples$chain1[,"p"], samples$chain2[,"p"]), main="p", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=p_exp[k], col="blue", lwd=3)

ggplot()+
  geom_point(aes(x=samples$chain1[,"wa"], y=samples$chain1[,"p"]), col="red")+
  geom_point(aes(x=samples$chain2[,"wa"], y=samples$chain2[,"p"]), col="blue")+
  ylab(expression(hat(p))) + xlab(expression(hat(w[A])))

MCMCsummary(object = samples, round = 2, Rhat=T, n.eff=T)

ggplot(as_tibble(samples$chain1[,"psi"])) +
  geom_histogram(aes(x = value), color = "white") +
  labs(x = "occupancy")


### d. UNIF (0.5,1) ----
#### i. Building model
occu_unif51 <- nimbleCode({
  
  # Likelihood
  for(i in 1:nsites) {
    # state model
    z[i] ~ dbern(psi)
    # obs model
    p.mix[i] <- z[i] * p 
    
    for(j in 1:nvisits){
      y[i,j] ~ dbern(p.mix[i])
      w.mix[i,j] <- (y[i,j] * wa) + ((1-y[i,j]) * (1-wb))
      w[i,j] ~ dbern(w.mix[i,j])
    }
  }
  
  # priors
  psi ~ dunif(0,1)
  p ~ dunif(0,1)
  wa ~ dunif(0.5,1)
  wb ~ dunif(0,1)
})

# specify parameters that need to be estimated
# parameters <- c("psi","p", "wa", "wb")
my.constants <- list(nsites = N, nvisits = J)

# pick initial values
initial.values <- function() list(psi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  wa = runif(1,0.5,1),
                                  wb = runif(1,0,1),
                                  z = as.numeric(apply(Y,1,sum)>0))
# # specify MCMC details
# n.iter <- 5000
# n.burnin <- 1000
# n.chains <- 2

#### ii. Run simulations ----
# simulations runs
k = 100

unif51 <- matrix(nrow=k, ncol=8)
# psi_exp <- vector() ; p_exp <- vector() ; wa_exp <- vector() ; wb_exp <- vector()

set.seed(123) # run before running simulations for each prior
for(i in 1:k){
  Y <- detect_simulate(N,J,psi,p)
  W <- apply(Y[,-1], c(1,2), function(x){
    if(x==0){
      W <- rbinom(n=1, size=1, prob=(1-w_b))
    } else{
      W <- rbinom(n=1, size=1, prob=w_a)
    }
    return(W)
  })
  sumI <- apply(W, 1,sum)


  # read in data
  my.data <- list(w = W)
  
  # create model as an R object (uncompiled model)
  res<- nimbleModel(code = occu_unif51, 
                    data = my.data,
                    constants = my.constants,
                    inits = initial.values())
  
  # create a MCMC configuration
  compileNimble(res)
  resConf <- configureMCMC(res)
  resConf$removeSamplers(c('p','wa'))
  resConf$addSampler(target = c('p','wa'),
                     type = 'AF_slice')
  
  # create a MCMC function and compile it
  res_MCMC <- buildMCMC(resConf)
  Cres_MCMC <- compileNimble(res_MCMC, project = res)
  
  # run NIMBLE
  samples <- runMCMC(mcmc = Cres_MCMC,
                     niter = n.iter,
                     nburnin = n.burnin,
                     nchain = n.chains)
  
  smry <- MCMCsummary(object = samples, round = 2, Rhat=F, n.eff=F)
  # ith simu, 8 parameter values (mean and median), jth prior
  unif51[i,1:4] <- smry$mean   # final_results[i,1:4,j]
  unif51[i,5:8] <- smry$`50%`  # final_results[i,5:8,j]
  # }
  
  # follow the advancement
  print(i)
}
# #combine the expected parameter according to each dataset (with seeds(123))
# exp <- cbind(psi_exp, p_exp, wa_exp, wb_exp)

# save the results of the simulations for each priors
write.csv(unif51, file="Using Informative priors/bayes_NIMBLE_50-100_unif51.csv")

#### iii. Check chains convergence (for the last simulation)
MCMCtrace(samples, pdf=F, ind=T, Rhat=T, n.eff=T, params=c("psi","p"))

par(mfrow=c(2,1))
hist(c(samples$chain1[,"wa"], samples$chain2[,"wa"]), main="wA", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=wa_exp[k], col="blue", lwd=3)
hist(c(samples$chain1[,"p"], samples$chain2[,"p"]), main="p", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=p_exp[k], col="blue", lwd=3)

ggplot()+
  geom_point(aes(x=samples$chain1[,"wa"], y=samples$chain1[,"p"]), col="red")+
  geom_point(aes(x=samples$chain2[,"wa"], y=samples$chain2[,"p"]), col="blue")+
  ylab(expression(hat(p))) + xlab(expression(hat(w[A])))

MCMCsummary(object = samples, round = 2, Rhat=T, n.eff=T)

ggplot(as_tibble(samples$chain1[,"psi"])) +
  geom_histogram(aes(x = value), color = "white") +
  labs(x = "occupancy")


### e. BETA WEAKLY Informative ----
#### i. Elicite prior
# we are looking for beta distribution parameters : alpha and beta
# beta_mode = (alpha - 1) / (alpha+beta-2)
# probability density function : pbeta(interval,alpha, beta)-prob
# we solve the system of 2 equations with 2 unknowns

# weakly informative prior with the mode
w_mode = 0.9 # alpha=9beta-8
w_prob = 0.5 # P(X<0.5) = 0.001

w_f <- function(b){
  answer = pbeta(0.5,((9*b)-8), b) - 0.01 # (9*b-8)
  return(answer)
}

root = uniroot(f=w_f, lower=1, upper=100) #, extendInt="yes")
w_beta=root$root
w_alpha=(9*w_beta)-8 #-8 # 3*w_beta-2

x <- seq(0,1,0.001)
plot(x=x, y=dbeta(x, shape1=w_alpha, shape2=w_beta), xlim=c(0.5,1))
abline(v=w_mode, col="red")

qbeta(0.99, shape1=w_alpha, shape2=w_beta)

#### ii. Building model
occu_beta_weak <- nimbleCode({
  
  # Likelihood
  for(i in 1:nsites) {
    # state model
    z[i] ~ dbern(psi)
    # obs model
    p.mix[i] <- z[i] * p 
    
    for(j in 1:nvisits){
      y[i,j] ~ dbern(p.mix[i])
      w.mix[i,j] <- (y[i,j] * wa) + ((1-y[i,j]) * (1-wb))
      w[i,j] ~ dbern(w.mix[i,j])
    }
  }
  
  # priors
  psi ~ dunif(0,1)
  p ~ dunif(0,1)
  wa ~ dbeta(w_alpha,w_beta)
  wb ~ dunif(0,1)
})

# specify parameters that need to be estimated
# parameters <- c("psi","p", "wa", "wb")
my.constants <- list(nsites = N, nvisits = J, w_alpha=w_alpha, w_beta=w_beta)

# pick initial values
initial.values <- function() list(psi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  wa = rbeta(1,w_alpha, w_beta),
                                  wb = runif(1,0,1),
                                  z = as.numeric(apply(Y,1,sum)>0))
# # specify MCMC details
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2

#### iii. Run simulations ----
# simulations runs
k = 100

# final_results <- array(dim=c(k,8,4))
beta_weak <- matrix(nrow=k, ncol=8)
# psi_exp <- vector() ; p_exp <- vector() ; wa_exp <- vector() ; wb_exp <- vector()

set.seed(123) # run before running simulations for each prior
for(i in 1:k){
  Y <- detect_simulate(N,J,psi,p)
  W <- apply(Y[,-1], c(1,2), function(x){
    if(x==0){
      W <- rbinom(n=1, size=1, prob=(1-w_b))
    } else{
      W <- rbinom(n=1, size=1, prob=w_a)
    }
    return(W)
  })
  sumI <- apply(W, 1,sum)
  
  # read in data
  my.data <- list(w = W)
  
  # create model as an R object (uncompiled model)
  res<- nimbleModel(code = occu_beta_weak, 
                    data = my.data,
                    constants = my.constants,
                    inits = initial.values())
  
  # create a MCMC configuration
  compileNimble(res)
  resConf <- configureMCMC(res)
  resConf$removeSamplers(c('p','wa'))
  resConf$addSampler(target = c('p','wa'),
                     type = 'AF_slice')
  
  # create a MCMC function and compile it
  res_MCMC <- buildMCMC(resConf)
  Cres_MCMC <- compileNimble(res_MCMC, project = res)
  
  # run NIMBLE
  samples <- runMCMC(mcmc = Cres_MCMC,
                     niter = n.iter,
                     nburnin = n.burnin,
                     nchain = n.chains)
  
  smry <- MCMCsummary(object = samples, round = 2, Rhat=F, n.eff=F)
  # ith simu, 8 parameter values (mean and median), jth prior
  beta_weak[i,1:4] <- smry$mean   # final_results[i,1:4,j]
  beta_weak[i,5:8] <- smry$`50%`  # final_results[i,5:8,j]
  # }
  
  # follow the advancement
  print(i)
}

# combine the expected parameter according to each dataset (with seeds(123))
# exp <- cbind(psi_exp, p_exp, wa_exp, wb_exp)

# save the results of the simulations for each priors
write.csv(beta_weak, file="Using Informative priors/bayes_NIMBLE_50-100_betaW.csv")

#### iiii. Check chains convergence (for the last simulation)
MCMCtrace(samples, pdf=F, ind=T, Rhat=T, n.eff=T, params=c("psi","p"))

par(mfrow=c(2,1))
hist(c(samples$chain1[,"wa"], samples$chain2[,"wa"]), main="wA", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=wa_exp[k], col="blue", lwd=3)
hist(c(samples$chain1[,"p"], samples$chain2[,"p"]), main="p", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=p_exp[k], col="blue", lwd=3)

ggplot()+
  geom_point(aes(x=samples$chain1[,"wa"], y=samples$chain1[,"p"]), col="red")+
  geom_point(aes(x=samples$chain2[,"wa"], y=samples$chain2[,"p"]), col="blue")+
  ylab(expression(hat(p))) + xlab(expression(hat(w[A])))

MCMCsummary(object = samples, round = 2, Rhat=T, n.eff=T)

ggplot(as_tibble(samples$chain1[,"psi"])) +
  geom_histogram(aes(x = value), color = "white") +
  labs(x = "occupancy")


### f. BETA HIGHLY Informative ----

#### i. Elicite prior
# highly informative prior with the mean
h_mode=0.9
h_prob = 0.8 # P(X<0.8) = 0.01

h_f <- function(b){
  # knowing the mode alpha=9*b-8
  answer = pbeta(h_prob,((9*b)-8), b)-0.01
  return(answer)
}

root = uniroot(f=h_f, lower=1, upper=100)
h_beta=root$root
h_alpha=(9*h_beta)-8

x <- seq(0,1,0.001)
plot(x=x, y=dbeta(x, shape1=h_alpha, shape2=h_beta))
abline(v=h_mean, col="red")

#### ii. Building model
occu_beta_high <- nimbleCode({
  
  # Likelihood
  for(i in 1:nsites) {
    # state model
    z[i] ~ dbern(psi)
    # obs model
    p.mix[i] <- z[i] * p 
    
    for(j in 1:nvisits){
      y[i,j] ~ dbern(p.mix[i])
      w.mix[i,j] <- (y[i,j] * wa) + ((1-y[i,j]) * (1-wb))
      w[i,j] ~ dbern(w.mix[i,j])
    }
  }
  
  # priors
  psi ~ dunif(0,1)
  p ~ dunif(0,1)
  wa ~ dbeta(h_alpha, h_beta)
  wb ~ dunif(0,1)
})

# specify parameters that need to be estimated
# parameters <- c("psi","p", "wa", "wb")
my.constants <- list(nsites = N, nvisits = J, h_alpha=h_alpha, h_beta=h_beta)

# pick initial values
initial.values <- function() list(psi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  wa = rbeta(1,h_alpha, h_beta),
                                  wb = runif(1,0,1),
                                  z = as.numeric(apply(Y,1,sum)>0))
# # specify MCMC details
# n.iter <- 5000
# n.burnin <- 1000
# n.chains <- 2

#### iii. Run simulations ----
# simulations runs
k = 100

# final_results <- array(dim=c(k,8,4))
beta_high <- matrix(nrow=k, ncol=8)
# psi_exp <- vector() ; p_exp <- vector() ; wa_exp <- vector() ; wb_exp <- vector()

set.seed(123) # run before running simulations for each prior
for(i in 1:k){
  Y <- detect_simulate(N,J,psi,p)
  W <- apply(Y[,-1], c(1,2), function(x){
    if(x==0){
      W <- rbinom(n=1, size=1, prob=(1-w_b))
    } else{
      W <- rbinom(n=1, size=1, prob=w_a)
    }
    return(W)
  })
  sumI <- apply(W, 1,sum)
  
  # read in data
  my.data <- list(w = W)
  
  # create model as an R object (uncompiled model)
  res<- nimbleModel(code = occu_beta_high, 
                    data = my.data,
                    constants = my.constants,
                    inits = initial.values())
  
  # create a MCMC configuration
  compileNimble(res)
  resConf <- configureMCMC(res)
  resConf$removeSamplers(c('p','wa'))
  resConf$addSampler(target = c('p','wa'),
                     type = 'AF_slice')
  
  # create a MCMC function and compile it
  res_MCMC <- buildMCMC(resConf)
  Cres_MCMC <- compileNimble(res_MCMC, project = res)
  
  # run NIMBLE
  samples <- runMCMC(mcmc = Cres_MCMC,
                     niter = n.iter,
                     nburnin = n.burnin,
                     nchain = n.chains)
  
  smry <- MCMCsummary(object = samples, round = 2, Rhat=F, n.eff=F)
  # ith simu, 8 parameter values (mean and median), jth prior
  beta_high[i,1:4] <- smry$mean   # final_results[i,1:4,j]
  beta_high[i,5:8] <- smry$`50%`  # final_results[i,5:8,j]
  # }
  
  # follow the advancement
  print(i)
}

# save the results of the simulations for each priors
write.csv(beta_high, file="Using Informative priors/bayes_NIMBLE_100_betaH.csv")

.#### iiii. Check chains convergence (for the last simulation)
MCMCtrace(samples, pdf=F, ind=T, Rhat=T, n.eff=T, params=c("psi","p"))

par(mfrow=c(2,1))
hist(c(samples$chain1[,"wa"], samples$chain2[,"wa"]), main="wA", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=wa_exp[k], col="blue", lwd=3)
hist(c(samples$chain1[,"p"], samples$chain2[,"p"]), main="p", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=p_exp[k], col="blue", lwd=3)

ggplot()+
  geom_point(aes(x=samples$chain1[,"wa"], y=samples$chain1[,"p"]), col="red")+
  geom_point(aes(x=samples$chain2[,"wa"], y=samples$chain2[,"p"]), col="blue")+
  ylab(expression(hat(p))) + xlab(expression(hat(w[A])))

MCMCsummary(object = samples, round = 2, Rhat=T, n.eff=T)

ggplot(as_tibble(c(samples$chain1[,"wa"], samples$chain2[,"wa"]))) +
  geom_histogram(aes(x = value), color = "white") +
  labs(x = "occupancy")


### g. PLOT ----

#### i. Plot prior distributions ----
x <- seq(0,1,0.001)
ggplot()+
  geom_point(aes(x=x, y=dunif(x, min=0,max=1)), shape=16, color="#264653")+
  geom_point(aes(x=x, y=dunif(x, min=0.5,max=1)), shape=16, color="#2a9d8f")+
  geom_point(aes(x=x, y=dbeta(x, shape1=w_alpha, shape2=w_beta)),shape=16, color="#f4a261")+
  geom_point(aes(x=x, y=dbeta(x, shape1=h_alpha, shape2=h_beta)),shape=16, color="#e76f51")+
  scale_x_continuous(breaks=seq(0,1,0.1), minor_breaks=NULL)+
  ylab("Density") + theme_light() + xlab(NULL)


#### ii. Boxplot of posterior distributions from simulations ----
k=100

# palette 
Palette <- c("#cad2c5", "#84a98c", "#52796f", "#354f52")

##### UNIF (0,1) ----
unif01 <- read.csv(file="Using Informative priors/Review1/bayes_NIMBLE_100_unif01.csv", 
                   col.names=c("X","p_mean", "psi_mean", "w_a_mean", "w_b_mean",  
                               "p_median", "psi_median", "w_a_median", "w_b_median"))

res_u01 <- as.data.frame(unif01[,2:9]) %>% 
  relocate(psi_median, .after=w_b_mean) %>% 
  relocate(psi_mean, .before=p_mean)

bias_u01 <- cbind(as.matrix(res_u01[,1:4]) - as.matrix(exp),
                  as.matrix(res_u01[,5:8]) - as.matrix(exp))

bias_u01_mean <- data.frame(param= c(bias_u01[,1:4]), # stack(bias_u01[,1:4])[,1], 
                            gro=c(rep("psi",k),rep("p",k),rep("wa",k),rep("wb",k)))
bias_u01_mean$gro <- factor(bias_u01_mean$gro, levels=c('wb','wa','p','psi'), ordered=T)

bias_u01_med <- data.frame(param= c(bias_u01[,5:8]), #stack(bias_u01[,5:8])[,1], 
                           gro=c(rep("psi",k),rep("p",k),rep("wa",k),rep("wb",k)))
bias_u01_med$gro <- factor(bias_u01_med$gro, levels=c('wb','wa','p','psi'), ordered=T)

#boxplot
u01 <- ggplot() +
  geom_boxplot(data=bias_u01_med, aes(x=param, y=gro, fill=gro), outlier.alpha = 0.6,
               alpha=0.7) + # fill="#E25E3E", col="#E25E3E", 
  geom_jitter(data=bias_u01_med, aes(x=param, y=gro, col=gro), size=0.4)+ # color="#E25E3E", 
  geom_vline(xintercept=0) + #xlim(-0.1,0.1) +
  theme_light(base_size = 10) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values=Palette) + scale_colour_manual(values=Palette) + 
  scale_y_discrete(labels=c( expression(hat(w[B])),  expression(hat(w[A])),
                             expression(hat(p)), expression(hat(psi)))) + 
  # scale_x_continuous(breaks=seq(-0.4, 0.3,by=0.1)) +
  xlim(-0.75,0.35) +
  theme(axis.text.y = element_text(family="serif", face="bold", size=20, 
                                   angle=360, vjust = 0.5),
        legend.position = "none")


##### UNIF (0.5,1) ----
unif51 <- read.csv(file="Using Informative priors/Review1/bayes_NIMBLE_100_unif51.csv", 
                   col.names=c("X","p_mean", "psi_mean", "w_a_mean", "w_b_mean",  
                               "p_median", "psi_median", "w_a_median", "w_b_median"))

res_u51 <- as.data.frame(unif51[,2:9]) %>% 
  relocate(psi_median, .after=w_b_mean) %>% 
  relocate(psi_mean, .before=p_mean)

bias_u51 <- cbind(as.matrix(res_u51[,1:4]) - as.matrix(exp),
                  as.matrix(res_u51[,5:8]) - as.matrix(exp))

bias_u51_mean <- data.frame(param=c(bias_u51[,1:4]), # stack(bias_u51[,1:4])[,1], 
                            gro=c(rep("psi",k),rep("p",k),rep("wa",k),rep("wb",k)))
bias_u51_mean$gro <- factor(bias_u51_mean$gro, levels=c('wb','wa','p','psi'), ordered=T)

bias_u51_med <- data.frame(param=c(bias_u51[,5:8]), #stack(bias_u51[,5:8])[,1], 
                           gro=c(rep("psi",k),rep("p",k),rep("wa",k),rep("wb",k)))
bias_u51_med$gro <- factor(bias_u51_med$gro, levels=c('wb','wa','p','psi'), ordered=T)

u51 <- ggplot() +
  geom_boxplot(data=bias_u51_med, aes(x=param, y=gro, fill=gro), outlier.alpha = 0.6,
               alpha=0.7) + # fill="#E25E3E", col="#E25E3E", 
  geom_jitter(data=bias_u51_med, aes(x=param, y=gro, col=gro), size=0.4) + #color="#E25E3E",
  geom_vline(xintercept=0) + #xlim(-0.1,0.1) +
  theme_light(base_size = 10) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values=Palette) + scale_colour_manual(values=Palette) + 
  scale_y_discrete(labels=c( expression(hat(w[B])),  expression(hat(w[A])),
                             expression(hat(p)), expression(hat(psi)))) + 
  # scale_x_continuous(breaks=seq(-0.4, 0.3,by=0.1)) + 
  xlim(-0.75,0.35) +
  theme(axis.text.y = element_text(family="serif", face="bold", size=20, 
                                   angle=360, vjust = 0.5),
        legend.position = "none")

##### BETA WEAKLY Informative ----
bWeak <- read.csv(file="Using Informative priors/Review1/bayes_NIMBLE_100_betaW.csv", 
                  col.names=c("X","p_mean", "psi_mean", "w_a_mean", "w_b_mean",  
                              "p_median", "psi_median", "w_a_median", "w_b_median"))

res_bWeak <- as.data.frame(bWeak[,2:9]) %>% 
  relocate(psi_median, .after=w_b_mean) %>% 
  relocate(psi_mean, .before=p_mean)

bias_bWeak <- cbind(as.matrix(res_bWeak[,1:4]) - as.matrix(exp),
                    as.matrix(res_bWeak[,5:8]) - as.matrix(exp))

bias_bWeak_mean <- data.frame(param=c(bias_bWeak[,1:4]), 
                              gro=c(rep("psi",k),rep("p",k),rep("wa",k),rep("wb",k)))
bias_bWeak_mean$gro <- factor(bias_bWeak_mean$gro, levels=c('wb','wa','p','psi'), ordered=T)

bias_bWeak_med <- data.frame(param=c(bias_bWeak[,5:8]), 
                             gro=c(rep("psi",k),rep("p",k),rep("wa",k),rep("wb",k)))
bias_bWeak_med$gro <- factor(bias_bWeak_med$gro, levels=c('wb','wa','p','psi'), ordered=T)

# Boxplot
bW <- ggplot() +
  geom_boxplot(data=bias_bWeak_med, aes(x=param, y=gro, fill=gro), outlier.alpha = 0.6,
               alpha=0.7) + # fill="#E25E3E", col="#E25E3E"
  geom_jitter(data=bias_bWeak_med, aes(x=param, y=gro, col=gro), size=0.4) + # color="#E25E3E",
  geom_vline(xintercept=0) + #xlim(-0.1,0.1) +
  theme_light(base_size = 10) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values=Palette) + scale_colour_manual(values=Palette) + 
  scale_y_discrete(labels=c( expression(hat(w[B])),  expression(hat(w[A])),
                             expression(hat(p)), expression(hat(psi)))) +
  # scale_x_continuous(breaks=seq(-0.4, 0.3,by=0.1)) + 
  xlim(-0.75,0.35) +
  theme(axis.text.y = element_text(family="serif", face="bold", size=20, 
                                   angle=360, vjust = 0.5),
        legend.position = "none")


##### BETA HIGHLY Informative ----
bHigh <- read.csv(file="Using Informative priors/Review1/bayes_NIMBLE_100_betaH.csv", 
                  col.names=c("X","p_mean", "psi_mean", "w_a_mean", "w_b_mean",  
                              "p_median", "psi_median", "w_a_median", "w_b_median"))

res_bHigh <- as.data.frame(bHigh[,2:9]) %>% 
  relocate(psi_median, .after=w_b_mean) %>% 
  relocate(psi_mean, .before=p_mean)

bias_bHigh <- cbind(as.matrix(res_bHigh[,1:4]) - as.matrix(exp),
                    as.matrix(res_bHigh[,5:8]) - as.matrix(exp))

bias_bHigh_mean <- data.frame(param=c(bias_bHigh[,1:4]), 
                              gro=c(rep("psi",k),rep("p",k),rep("wa",k),rep("wb",k)))
bias_bHigh_mean$gro <- factor(bias_bHigh_mean$gro, levels=c('wb','wa','p','psi'), ordered=T)

bias_bHigh_med <- data.frame(param=c(bias_bHigh[,5:8]), 
                             gro=c(rep("psi",k),rep("p",k),rep("wa",k),rep("wb",k)))
bias_bHigh_med$gro <- factor(bias_bHigh_med$gro, levels=c('wb','wa','p','psi'), ordered=T)

# Boxplot
bH <- ggplot() +
  geom_boxplot(data=bias_bHigh_med, aes(x=param, y=gro, fill=gro), outlier.alpha = 0.6,
               alpha=0.7) + # fill="#E25E3E", col="#E25E3E"
  geom_jitter(data=bias_bHigh_med, aes(x=param, y=gro, col=gro), size=0.4) + # color="#E25E3E",
  geom_vline(xintercept=0) + #xlim(-0.1,0.1) +
  theme_light(base_size = 10) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values=Palette) + scale_colour_manual(values=Palette) + 
  scale_y_discrete(labels=c( expression(hat(w[B])),  expression(hat(w[A])),
                             expression(hat(p)), expression(hat(psi)))) +
  # scale_x_continuous(breaks=seq(-0.4, 0.3,by=0.1)) + 
  xlim(-0.75,0.35) +
  theme(axis.text.y = element_text(family="serif", face="bold", size=20, 
                                   angle=360, vjust = 0.5),
        legend.position = "none")


cowplot::plot_grid(u01, u51, bW,bH, nrow=4, labels="AUTO")


### h. Sensitivity analysis ----
# initialize parameters to simulate data
N <- 30  ; J <- 36 ; psi <- 0.8
p <- 0.5 ; w_b <- 0.7 

detect_simulate <- function(N, J, psi, p){
  1:N %>%
    map_df(function(i) {
      Z <- rbinom(1, 1, prob=psi)
      Y <- rbinom(J ,1, prob=Z*p)
      return(c(Z=Z,Y=Y))
    })
}

#### Define w_a ----
w_a <- c(0.2,0.5,0.8) # to change to 0.2, 0.5 and 0.8

i<-3 # 1,2,3

  # simulate data
  set.seed(123)
  
  Y <- detect_simulate(N,J,psi,p)
  W <- apply(Y[,-1], c(1,2), function(x){
    if(x==0){
      W <- rbinom(n=1, size=1, prob=(1-w_b))
    } else{
      W <- rbinom(n=1, size=1, prob=w_a[i])
    }
    return(W)
  })
  sumI <- apply(W, 1,sum)
  
  # psi_exp <- sum(Y[,1])/N
  # p_exp <- sum(rowSums(Y[,-1]))/(J*sum(Y[,1]))
  # wa_exp <- sum(W[which(Y[,-1]==1)])/sum(Y[,-1])
  # wb_exp <- length(which(W[which(Y[,-1]==0)]==0)) / length(which(Y[,-1]==0))
  
  # Elicite highly informative prior
  h_mode = w_a[i]
  h_prob = w_a[i]-0.1
  
  h_f <- function(b){
    # a=(b+3)/4 for m=0.2
    # a=b for m=0.5
    # a=4*b-3for m=0.8
    answer= pbeta(h_prob,((h_mode*b-2*h_mode+1)/(1-h_mode)),b) - 0.01
    return(answer)
  }
  
  root = uniroot(f=h_f, lower=1, upper=100)
  h_beta=root$root
  h_alpha= ((h_mode*h_beta-2*h_mode+1)/(1-h_mode)) # (h_beta+3)/4
  
  x <- seq(0,1,0.001)
  plot(x=x, y=dbeta(x, shape1=h_alpha, shape2=h_beta))
  abline(v=h_mode, col="red")
  
  # Building model
  occu_beta_high <- nimbleCode({
    
    # Likelihood
    for(i in 1:nsites) {
      # state model
      z[i] ~ dbern(psi)
      # obs model
      p.mix[i] <- z[i] * p 
      
      for(j in 1:nvisits){
        y[i,j] ~ dbern(p.mix[i])
        w.mix[i,j] <- (y[i,j] * wa) + ((1-y[i,j]) * (1-wb))
        w[i,j] ~ dbern(w.mix[i,j])
      }
    }
    
    # priors
    psi ~ dunif(0,1)
    p ~ dunif(0,1)
    wa ~ dbeta(h_alpha, h_beta)
    wb ~ dunif(0,1)
  })
  
  # specify parameters that need to be estimated
  my.constants <- list(nsites = N, nvisits = J, h_alpha=h_alpha, h_beta=h_beta)
  
  # pick initial values
  initial.values <- function() list(psi = runif(1,0,1),
                                    p = runif(1,0,1),
                                    wa = rbeta(1,h_alpha, h_beta),
                                    wb = runif(1,0,1),
                                    z = as.numeric(apply(Y,1,sum)>0))
  # specify MCMC details
  n.iter <- 5000
  n.burnin <- 1000
  n.chains <- 2
  
  # read in data
  my.data <- list(w = W)
  
  # create model as an R object (uncompiled model)
  res<- nimbleModel(code = occu_beta_high, 
                    data = my.data,
                    constants = my.constants,
                    inits = initial.values())
  
  # create a MCMC configuration
  compileNimble(res)
  resConf <- configureMCMC(res)
  resConf$removeSamplers(c('p','wa'))
  resConf$addSampler(target = c('p','wa'),
                     type = 'AF_slice')
  
  # create a MCMC function and compile it
  res_MCMC <- buildMCMC(resConf)
  Cres_MCMC <- compileNimble(res_MCMC, project = res)
  
  # run NIMBLE
  samples <- runMCMC(mcmc = Cres_MCMC,
                     niter = n.iter,
                     nburnin = n.burnin,
                     nchain = n.chains)

  # Check chains convergence
  MCMCtrace(samples, pdf=F, ind=T, Rhat=T, n.eff=T, gvals=c(psi,p), params=c("psi","p")) #, type="trace")
  MCMCtrace(samples, pdf=F, ind=T, Rhat=T, n.eff=T, gvals=c(w_a[i], w_b), params=c("wa","wb"))


MCMCsummary(object = samples, round = 2, Rhat=T, n.eff=T)

# density plots for both chains
par(mfrow=c(2,1))
hist(c(samples$chain1[,"wa"], samples$chain2[,"wa"]), main="wA", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=wa_exp, col="blue", lwd=3)
hist(c(samples$chain1[,"p"], samples$chain2[,"p"]), main="p", 
     xlim=c(0.2,1), breaks=25, xlab=NULL)
abline(v=p_exp, col="blue", lwd=3)

ggplot()+
  geom_point(aes(x=samples$chain1[,"wa"], y=samples$chain1[,"p"]), col="red")+
  geom_point(aes(x=samples$chain2[,"wa"], y=samples$chain2[,"p"]), col="blue")+
  ylab(expression(hat(p))) + xlab(expression(hat(w[A])))

ggplot(as_tibble(c(samples$chain1[,"wa"], samples$chain2[,"wa"]))) +
  geom_histogram(aes(x = value), color = "white") +
  labs(x = "occupancy")

