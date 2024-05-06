# Using informative prior to account for identifiability issues 
# in occupancy models with identification errors
#---------------------------------------------------------------#

library(purrr) ; library(dplyr) ; library(ggplot2) ; 
library(cowplot) ; library(ggridges) ; library(R2jags)
set.seed(123)

## 0. Data Simulation ----

# Initial settings to simulate data
N <- 30       # sites
J <- 12 # 36  # visits
psi <- 0.8    # occupancy probability
p <- 0.5      # detection probability
w_a <- 0.9    # correct identification prob (sensitivity)
w_b <- 0.7    # correct non-identification prob
# from https://computo.sfds.asso.fr/published-202204-deeplearning-occupancy-lynx/

# function to simulate detection matrix
detect_simulate <- function(N, J, psi, p){
  1:N %>%
    map_df(function(i) {
      Z <- rbinom(1, 1, prob=psi)
      Y <- rbinom(J ,1, prob=Z*p)
      return(c(Z=Z,Y=Y))
    })
}

# Negative Log Likelihood
NegLik <- function(occ){
  psi <- plogis(occ[1])
  p <- plogis(occ[2])
  w_a <- plogis(occ[3])
  w_b <- plogis(occ[4])
  cp_a <- (((1-p)*(1-w_b)+(p*w_a))^sumI) * (((1-p)*w_b+p*(1-w_a))^(J-sumI))
  cp_b <- (w_b^(J-sumI)) * ((1-w_b)^sumI)
  loglik <- log(cp_a*psi + cp_b*(1-psi))
  return(-sum(loglik))
}



## 1. Classical Estimation from the identification layer ----

### a. Determining initial values for Maximum Likelihood Estimation ----
Y <- detect_simulate(N,J,psi,p)
Z <- Y[,1]
Y <- Y[,-1]

# Identification Matrix from detection matrix
W <- apply(Y, c(1,2), function(x){
  if(x==0){
    W <- rbinom(n=1, size=1, prob=(1-w_b))
  } else{
    W <- rbinom(n=1, size=1, prob=w_a)
  }
  return(W)
})
colnames(W) <- NULL

head(as.data.frame(Y)) ; head(W)

# Formatting data by site
sumI <- apply(W, 1,sum)
sumI

# defining inital values
psi0 <- length(sumI[which(sumI!=0)])/N
p0 <- sqrt(w_a*p)
w_a0 <- sqrt(w_a*p)
w_b0 <- 1- mean(sumI[which(Z==0)])/J

### b. Estimation from the Minimum of the negative log-Likelihood ----
my_optim <- optim(par=c(psi0,p0,w_a0, w_b0), fn=NegLik, hessian=T)

print(paste("PSI=", plogis(my_optim$par[1]),
            "p=", plogis(my_optim$par[2]),
            "w_a=", plogis(my_optim$par[3]),
            "w_b=", plogis(my_optim$par[4])))



### c. Simulations for estimates ----
# Once for J=12 and once for J=36
k = 1000
# res_simu12 <- matrix(nrow=k, ncol=4)
res_simu36 <- matrix(nrow=k, ncol=4)

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

  my_optimS <- optim(par=c(psi0,p0,w_a0, w_b0), fn=NegLik, hessian=T)
  # res_simu12[i,] <- plogis(my_optimS$par)
  res_simu36[i,] <- plogis(my_optimS$par)
}

# colnames(res_simu12) <- c("psi", "p", "w_a", "w_b")
# head(res_simu12)
colnames(res_simu36) <- c("psi", "p", "w_a", "w_b")

### d. Plot ----
res_simu <- res_simu12
#### J=12 ----
denspsi <- density(res_simu12[,1], bw=0.02)
densp <- density(res_simu12[,2], bw=0.02)
denswa <- density(res_simu12[,3], bw=0.02)
denswb <- density(res_simu12[,4], bw=0.02)

psi12 <- ggplot()+
  geom_histogram(aes(res_simu[,1], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray", binwidth = 0.05) + #, binwidth = 0.1 "#B0D9B1"
  geom_density(aes(x = res_simu[,1], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2,linewidth=0.8, bw=0.02)+ #, after_stat(count) "#618264"
  geom_vline(xintercept=psi, col="darkred") + # darkred "#0F2C59"
  geom_point(aes(x=psi, y=mean(denspsi$y[which(round(denspsi$x, digit=2)==psi)]/max(denspsi$y))), 
             col="darkred", size=3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  xlab(NULL) + theme_light(base_size = 10)+ # + xlim(0,1)
  ylab(NULL)
  # ylab(expression(hat(psi))) +
  theme(axis.title.y = element_text(family="serif", face="bold", size=20, 
                                    angle=360, vjust = 0.5))

p12 <- ggplot()+
  geom_histogram(aes(res_simu[,2], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray",binwidth = 0.05) +
  geom_density(aes(x = res_simu[,2], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2, linewidth=0.8, bw=0.02)+ 
  geom_vline(xintercept=p, col="darkred") + #"#eb5e28"
  geom_point(aes(x=p, y=mean(densp$y[which(round(densp$x, digit=2)==p)]/max(densp$y))), 
             col="darkred", size=3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1)) + # n.breaks=10, limits= c(0,1)
  xlab(NULL) + theme_light(base_size = 10) +
  ylab(NULL)
  ylab(expression(hat(p))) +
  theme(axis.title.y = element_text(family="serif", face="bold", size=20, 
                                    angle=360, vjust = 0.5))

wa12 <- ggplot()+
  geom_histogram(aes(res_simu[,3], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray", binwidth = 0.05) +
  geom_density(aes(x = res_simu[,3], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2,linewidth=0.8, bw=0.02)+ 
  geom_vline(xintercept=w_a, col="darkred") +
  geom_point(aes(x=w_a, y=mean(denswa$y[which(round(denswa$x, digit=2)==w_a)]/max(denswa$y))),
            col="darkred", size=3)+
  scale_x_continuous(n.breaks=10, limits=c(0,1)) +
  xlab(NULL) + theme_light(base_size = 10)+
  ylab(NULL)
  ylab(expression(hat(w[A]))) +
  theme(axis.title.y = element_text(family="serif", face="bold", size=20, 
                                    angle=360, vjust = 0.5))
wb12 <- ggplot()+
  geom_histogram(aes(res_simu[,4], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray", binwidth = 0.05) +
  geom_density(aes(x = res_simu[,4], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2,linewidth=0.8, bw=0.02)+
  geom_vline(xintercept=w_b, col="darkred") +
  geom_point(aes(x=w_b, y=mean(denswb$y[which(round(denswb$x, digit=2)==w_b)]/max(denswb$y))), 
             col="darkred", size=3) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1)) +
  ylab(NULL) + theme_light(base_size = 10) +
  xlab(NULL)
  xlab("J=12") +  ylab(expression(hat(w[B]))) +
  theme(axis.title.x = element_text(family="serif", face="bold", size=20),
        axis.title.y = element_text(family="serif", face="bold", size=20, 
                                    angle=360, vjust = 0.5))

J12 <- plot_grid(psi12,p12, wa12, wb12, ncol=1,
          label_x=0, label_y=0, label_size = 8, hjust = -0.5, vjust = -0.5)
J12

#### J=36 ----
denspsi36 <- density(res_simu36[,1], bw=0.02)
densp36 <- density(res_simu36[,2], bw=0.02)
denswa36 <- density(res_simu36[,3], bw=0.02)
denswb36 <- density(res_simu36[,4], bw=0.02)

psi36 <- ggplot()+
  geom_histogram(aes(res_simu36[,1], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray", binwidth = 0.05) + 
  geom_density(aes(x = res_simu36[,1], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2,linewidth=0.8, bw=0.02)+ 
  geom_vline(xintercept=psi, col="darkred") +
  geom_point(aes(x=psi, y=mean(denspsi36$y[which(round(denspsi36$x, digit=2)==psi)]/max(denspsi36$y))), 
             col="darkred", size=3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1)) +
  ylab(NULL) + xlab(NULL) + theme_light(base_size = 10) 

p36 <- ggplot()+
  geom_histogram(aes(res_simu36[,2], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray",binwidth = 0.05) +
  geom_density(aes(x = res_simu36[,2], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2, linewidth=0.8, bw=0.02)+ 
  geom_vline(xintercept=p, col="darkred") +
  geom_point(aes(x=p, y=mean(densp36$y[which(round(densp36$x, digit=2)==p)]/max(densp36$y))), 
             col="darkred", size=3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1)) +
  ylab(NULL) + xlab(NULL) + theme_light(base_size = 10)

wa36 <- ggplot()+
  geom_histogram(aes(res_simu36[,3], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray", binwidth = 0.05) +
  geom_density(aes(x = res_simu36[,3], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2,linewidth=0.8, bw=0.02)+ 
  geom_vline(xintercept=w_a, col="darkred") +
  geom_point(aes(x=w_a, y=mean(denswa36$y[which(round(denswa36$x, digit=2)==w_a)]/max(denswa36$y))),
             col="darkred", size=3)+
  scale_x_continuous(n.breaks=10, limits=c(0,1)) +
  ylab(NULL) + xlab(NULL) + theme_light(base_size = 10)
wb36 <- ggplot()+
  geom_histogram(aes(res_simu36[,4], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray", binwidth = 0.05) +
  geom_density(aes(x = res_simu36[,4], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2,linewidth=0.8, bw=0.02)+
  geom_vline(xintercept=w_b, col="darkred") +
  geom_point(aes(x=w_b, y=mean(denswb36$y[which(round(denswb36$x, digit=2)==w_b)]/max(denswb36$y))), 
             col="darkred", size=3) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1)) +
  ylab(NULL) + theme_light(base_size = 10) +
  xlab(NULL) + draw_plot_label(label="J=12",x=0.5, vjust=-0.5)
  xlab("J=36") +
  theme(axis.title.x = element_text(family="serif", face="bold", size=20))
  # draw_text("J=12", x=0, y=0, vjust=2)

J36 <- plot_grid(psi36,p36, wa36, wb36, ncol=1, 
                 label_x=0, label_y=0, label_size = 8, hjust = -0.5, vjust = -0.5)
J36

J1236 <- plot_grid(J12,J36, byrow=FALSE, ncol=2)
J1236

#### Standard Deviation ----
sqrt(sum((res_simu12[,1]-mean(res_simu12[,1]))^2)/(k-1))
sqrt(sum((res_simu36[,1]-mean(res_simu36[,1]))^2)/(k-1))

#### Biais ----
mean(res_simu36[,2])-0.5 ; mean(res_simu36[,3])-0.9



## 2. Constrained optimisation wa > 1-wb ----

# constraints declaration for constrOprim() function
ui <- c(0,0,1,1)
ci <- 1

k = 1000 # number of simulations
res_val <- matrix(nrow=k, ncol=6) # results with constraint
res_sansval <- matrix(nrow=k, ncol=6) # results without constraint

for(i in 1:k){
  # for each simulation w_a and w_b are set to a value between 0.5 and 0.95
  w_bV <- round(runif(1, 0.5, 0.95)/5,2)*5 
  w_aV <-  round(runif(1, 0.5, 0.95)/5,2)*5 
  res_val[i,1:2] <- c(w_aV, w_bV) ; res_sansval[i,1:2] <- c(w_aV, w_bV)
  # simulations
  Y <- detect_simulate(N,J,psi,p)
  Z <- Y[,1]
  Y <- Y[,-1]
  W <- apply(Y, c(1,2), function(x){
    if(x==0){
      W <- rbinom(n=1, size=1, prob=(1-w_bV))
    } else{
      W <- rbinom(n=1, size=1, prob=w_aV)
    }
    return(W)
  })
  sumV <- apply(W, 1,sum)
  
  NegLikV <- function(occ){
    psi <- plogis(occ[1])    
    p <- plogis(occ[2]) 
    w_aV <- plogis(occ[3])
    w_bV <- plogis(occ[4])
    cp_a <- (((1-w_bV)*(1-p)+(w_aV*p))^sumV) * ((w_bV*(1-p)+(1-w_aV)*p)^(J-sumV))
    cp_b <- (w_bV^(J-sumV)) * ((1-w_bV)^sumV)
    loglik <- log(cp_a*psi + cp_b*(1-psi))
    return(-sum(loglik))
  }
  
  # Constrained optimisation
  opt_val <- constrOptim(theta=c(psi,p,w_aV+0.01,w_bV), f=NegLikV, NULL, ui=ui, ci=ci)
  res_val[i,3:6] <- plogis(opt_val$par)
  
  # Not-Constrained optimisation
  opt_sansval <- optim(par=c(psi,p,w_aV+0.01,w_bV), fn=NegLikV, hessian=T)
  res_sansval[i,3:6] <- plogis(opt_sansval$par)
}



### plot ----
res_val <- as.data.frame(res_val)
colnames(res_val) <- c("wAinit", "wBinit", "psi", "p", "wA","wB")
res_sansval <- as.data.frame(res_sansval)
colnames(res_sansval) <- c("wAinit", "wBinit", "psi", "p", "wA","wB")

# wB
res_val2 <- res_val %>% 
  group_by(wBinit) %>% # as.factor()
  summarise(bias_wa=median(wA), bias_wa_01=quantile(wA,0.1), bias_wa_09=quantile(wA,0.9),  # (sqrt(mean((wA-wAinit)^2))),
            bias_p=median(p), bias_p_01=quantile(p,0.1), bias_p_09=quantile(p,0.9))# (sqrt(mean(p-0.5)^2))) 
res_val2[,2:4] <- res_val2[,2:4] - res_val2$wBinit
res_val2[,5:7] <- res_val2[,5:7] - 0.5

res_sansval2 <- res_sansval %>% 
  group_by(wBinit) %>% 
  summarise(bias_wa=median(wA), bias_wa_01=quantile(wA,0.1), bias_wa_09=quantile(wA,0.9),  # (sqrt(mean((wA-wAinit)^2))),
            bias_p=median(p), bias_p_01=quantile(p,0.1), bias_p_09=quantile(p,0.9))# (sqrt(mean(p-0.5)^2))) 
res_sansval2[,2:4] <- res_sansval2[,2:4] - res_sansval2$wBinit
res_sansval2[,5:7] <- res_sansval2[,5:7] - 0.5

res3 <- rbind(res_val2, res_sansval2)
res3$cons <- c(rep("with constraint",10), rep("without constaint",10))

#wA
res_val2A <- res_val %>% 
  group_by(wAinit) %>% # as.factor()
  summarise(bias_wa=median(wA), bias_wa_01=quantile(wA,0.1), bias_wa_09=quantile(wA,0.9),  # (sqrt(mean((wA-wAinit)^2))),
            bias_p=median(p), bias_p_01=quantile(p,0.1), bias_p_09=quantile(p,0.9))# (sqrt(mean(p-0.5)^2))) 
res_val2A[,2:4] <- res_val2A[,2:4] - res_val2A$wAinit
res_val2A[,5:7] <- res_val2A[,5:7] - 0.5

res_sansval2A <- res_sansval %>% 
  group_by(wAinit) %>% 
  summarise(bias_wa=median(wA), bias_wa_01=quantile(wA,0.1), bias_wa_09=quantile(wA,0.9),  # (sqrt(mean((wA-wAinit)^2))),
            bias_p=median(p), bias_p_01=quantile(p,0.1), bias_p_09=quantile(p,0.9))# (sqrt(mean(p-0.5)^2))) 
res_sansval2A[,2:4] <- res_sansval2A[,2:4] - res_sansval2A$wAinit
res_sansval2A[,5:7] <- res_sansval2A[,5:7] - 0.5

res3A <- rbind(res_val2A, res_sansval2A)
res3A$cons <- c(rep("with constraint",10), rep("without constaint",10))

plot_grid(
  ggplot(data=res3A)+
    geom_point(aes(x=wAinit, y=bias_p, col=cons), size=2) +
    geom_line(aes(x=wAinit, y=bias_p, col=cons))+
    geom_ribbon(aes(x=wAinit, ymin=bias_p_01, ymax=bias_p_09, fill=cons), alpha=0.2)+
    scale_color_manual(values=c("darkblue", "#808080")) + 
    scale_fill_manual(values=c("darkblue", "#808080")) +
    theme_light(base_size = 12) + xlab(expression(w[A])) + 
    ylab(expression('bias'~hat(p))) + #+ ylim(0,0.24)
    theme(legend.position=c(0.8, 0.15) ,legend.direction = "vertical",
          legend.title =element_blank()),
  ggplot(data=res3A)+
    geom_ribbon(aes(x=wAinit, ymin=bias_wa_01, ymax=bias_wa_09, fill=cons), alpha=0.2)+
    geom_point(aes(x=wAinit, y=bias_wa, col=cons), size=2) +
    geom_line(aes(x=wAinit, y=bias_wa, col=cons))+
    scale_color_manual(values=c("darkblue", "#808080")) + # c("darkblue", "darkgray")
    scale_fill_manual(values=c("darkblue", "#808080")) + # c("#76b6ff", "#808080")
    theme_light(base_size = 12) + xlab(expression(w[A])) + 
    ylab(expression('bias'~hat(w[A]))) + #+ ylim(0,0.24)
    theme(legend.position=c(0.8, 0.85) ,legend.direction = "vertical",
          legend.title =element_blank()),
  
ggplot(data=res3)+
  geom_point(aes(x=wBinit, y=bias_p, col=cons), size=2) +
  geom_line(aes(x=wBinit, y=bias_p, col=cons))+
  geom_ribbon(aes(x=wBinit, ymin=bias_p_01, ymax=bias_p_09, fill=cons), alpha=0.2)+
  scale_color_manual(values=c("darkblue", "#808080")) + 
  scale_fill_manual(values=c("darkblue", "#808080")) +
  theme_light(base_size = 12) + xlab(expression(w[B])) + 
  ylab(expression('bias'~hat(p))) + #+ ylim(0,0.24)
  theme(legend.position=c(0.8, 0.15) ,legend.direction = "vertical",
        legend.title =element_blank()),
ggplot(data=res3)+
  geom_ribbon(aes(x=wBinit, ymin=bias_wa_01, ymax=bias_wa_09, fill=cons), alpha=0.2)+
  geom_point(aes(x=wBinit, y=bias_wa, col=cons), size=2) +
  geom_line(aes(x=wBinit, y=bias_wa, col=cons))+
  scale_color_manual(values=c("darkblue", "#808080")) + # c("darkblue", "darkgray")
  scale_fill_manual(values=c("darkblue", "#808080")) + # c("#76b6ff", "#808080")
  theme_light(base_size = 12) + xlab(expression(w[B])) + 
  ylab(expression('bias'~hat(w[A]))) + #+ ylim(0,0.24)
  theme(legend.position=c(0.8, 0.85) ,legend.direction = "vertical",
        legend.title =element_blank()),
byrow = F
)


## 3. Using an informative prior in Bayesian framework ----

### a. Initialize parameters to simulate data ----
N <- 30  ; J <- 36 ; psi <- 0.8
p <- 0.5 ; w_a <- 0.9 ; w_b <- 0.7 
set.seed(123)

detect_simulate <- function(N, J, psi, p){
  1:N %>%
    map_df(function(i) {
      Z <- rbinom(1, 1, prob=psi)
      Y <- rbinom(J ,1, prob=Z*p)
      return(c(Z=Z,Y=Y))
    })
}

### b. Building model ----
# one for each w_a prior
occuFP_unif01 <- function() {
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
}

occuFP_unif81 <- function() {
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
  wa ~ dunif(0.8,1)
  wb ~ dunif(0,1)
}

occuFP_beta <- function() {
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
  wa ~ dbeta(45,5) # dunif(0.8,1) dunif(0,1)
  wb ~ dunif(0,1)
}

hist(rbeta(seq(0,1,0.001),45,5), main="Shape of beta distribution B(45,5)")


# specify parameters that need to be estimated
parameters <- c("psi","p", "wa", "wb")

# specify nb iterations for burn-in and final inference
nb.burnin <- 1000
nb.iterations <- 5000

# simulations runs
k = 100

unif01 <- matrix(nrow=k, ncol=8)
unif81 <- matrix(nrow=k, ncol=8)
beta455 <- matrix(nrow=k, ncol=8)

psi_exp <- vector() ; p_exp <- vector()

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
  
  psi_exp[i] <- sum(Y[,1])/N
  p_exp[i] <- sum(rowSums(Y[,-1]))/(J*sum(Y[,1]))
  
  # # input data
  # datax <- list(w = W, nsites = N, nvisits = J)
  # 
  # # list of lists of initial values (one for each MCMC chain)
  # init1 <- list(psi = 0.1, p = 0.1, wa = 0.8, wb = 0.5, 
  #               z = as.numeric(apply(Y,1,sum)>0)) # fp = 0.1,
  # init2 <- list(psi = 0.5, p = 0.5, wa = 0.9, wb = 0.7, 
  #               z = as.numeric(apply(Y,1,sum)>0))
  # init3 <- list(psi = 0.9, p = 0.9, wa = 0.95, wb = 0.9, 
  #               z = as.numeric(apply(Y,1,sum)>0))
  # inits <- list(init1, init2, init3)
  # 
  # # run Jags
  # res_unif1 <- jags(data = datax,
  #             inits = inits,
  #             parameters.to.save = parameters,
  #             model.file = occuFP_unif01, 
  #             n.chains = 3,
  #             n.iter = nb.iterations, # includes burn-in!
  #             n.burnin = nb.burnin)
  # 
  # unif01[i,1:4] <- unlist(res_unif1$BUGSoutput$mean)[2:5]
  # unif01[i,5:8] <- unlist(res_unif1$BUGSoutput$median)[2:5]
  # 
  # res_unif2 <- jags(data = datax,
  #                  inits = inits,
  #                  parameters.to.save = parameters,
  #                  model.file = occuFP_unif81, 
  #                  n.chains = 3,
  #                  n.iter = nb.iterations, # includes burn-in!
  #                  n.burnin = nb.burnin)
  # 
  # unif81[i,1:4] <- unlist(res_unif2$BUGSoutput$mean)[2:5]
  # unif81[i,5:8] <- unlist(res_unif2$BUGSoutput$median)[2:5]
  # 
  # res_beta <- jags(data = datax,
  #                  inits = inits,
  #                  parameters.to.save = parameters,
  #                  model.file = occuFP_beta, 
  #                  n.chains = 3,
  #                  n.iter = nb.iterations, # includes burn-in!
  #                  n.burnin = nb.burnin)
  # 
  # beta455[i,1:4] <- unlist(res_beta$BUGSoutput$mean)[2:5]
  # beta455[i,5:8] <- unlist(res_beta$BUGSoutput$median)[2:5]
}

### c. Posterior distribution ----

# palette 
Palette <- c("#cad2c5", "#84a98c", "#52796f", "#354f52")

#### UNIF (0,1) ----
colnames(unif01) <- c("p_mean", "psi_mean", "w_a_mean", "w_b_mean",  
                       "p_median", "psi_median", "w_a_median", "w_b_median") 
res_u01 <- as.data.frame(unif01) %>% 
  relocate(psi_median, .after=w_b_mean) %>% 
  relocate(psi_mean, .before=p_mean)
# bias_u01 <- as.data.frame(t(apply(res_u01, 1, function(x){x-rep(c(0.8,0.5,0.9,0.7),2)})))
bias_u01 <- as.data.frame(t(apply(res_u01, 1, function(x){x-rep(c(0.8,0.5,0.9,0.7),2)})))

bias_u01_mean <- data.frame(param=stack(bias_u01[,1:4])[,1], 
                             gro=c(rep("psi",100),rep("p",100),rep("wa",100),rep("wb",100)))
bias_u01_mean$gro <- factor(bias_u01_mean$gro, levels=c('wb','wa','p','psi'), ordered=T)

bias_u01_med <- data.frame(param=stack(bias_u01[,5:8])[,1], 
                            gro=c(rep("psi",100),rep("p",100),rep("wa",100),rep("wb",100)))
bias_u01_med$gro <- factor(bias_u01_med$gro, levels=c('wb','wa','p','psi'), ordered=T)

#boxplot
u01 <- ggplot() +
  geom_boxplot(data=bias_u01_mean, aes(x=param, y=gro, fill=gro), outlier.alpha = 0.6,
               alpha=0.7) + # fill="#E25E3E", col="#E25E3E", 
  geom_jitter(data=bias_u01_mean, aes(x=param, y=gro, col=gro), size=0.4)+ # color="#E25E3E", 
  geom_vline(xintercept=0) + #xlim(-0.1,0.1) +
  theme_light(base_size = 10) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values=Palette) + scale_colour_manual(values=Palette) + 
  scale_y_discrete(labels=c( expression(hat(w[B])),  expression(hat(w[A])),
                             expression(hat(p)), expression(hat(psi)))) + 
  # scale_x_continuous(breaks=seq(-0.4, 0.3,by=0.1)) +
  xlim(-0.25,0.2) +
  theme(axis.text.y = element_text(family="serif", face="bold", size=20, 
                                   angle=360, vjust = 0.5),
        legend.position = "none")


#### UNIF (0.8,1) ----
colnames(unif81) <- c("p_mean", "psi_mean", "w_a_mean", "w_b_mean",  
                      "p_median", "psi_median", "w_a_median", "w_b_median") 
res_u081 <- as.data.frame(unif81) %>% 
  relocate(psi_median, .after=w_b_mean) %>% 
  relocate(psi_mean, .before=p_mean)
bias_u081 <- as.data.frame(t(apply(res_u081, 1, function(x){x-rep(c(0.8,0.5,0.9,0.7),2)})))

bias_u081_mean <- data.frame(param=stack(bias_u081[,1:4])[,1], 
                            gro=c(rep("psi",100),rep("p",100),rep("wa",100),rep("wb",100)))
bias_u081_mean$gro <- factor(bias_u081_mean$gro, levels=c('wb','wa','p','psi'), ordered=T)

bias_u081_med <- data.frame(param=stack(bias_u081[,5:8])[,1], 
                           gro=c(rep("psi",100),rep("p",100),rep("wa",100),rep("wb",100)))
bias_u081_med$gro <- factor(bias_u081_med$gro, levels=c('wb','wa','p','psi'), ordered=T)

u81 <- ggplot() +
  geom_boxplot(data=bias_u081_mean, aes(x=param, y=gro, fill=gro), outlier.alpha = 0.6,
               alpha=0.7) + # fill="#E25E3E", col="#E25E3E", 
  geom_jitter(data=bias_u081_mean, aes(x=param, y=gro, col=gro), size=0.4) + #color="#E25E3E",
  geom_vline(xintercept=0) + #xlim(-0.1,0.1) +
  theme_light(base_size = 10) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values=Palette) + scale_colour_manual(values=Palette) + 
  scale_y_discrete(labels=c( expression(hat(w[B])),  expression(hat(w[A])),
                             expression(hat(p)), expression(hat(psi)))) + 
  # scale_x_continuous(breaks=seq(-0.4, 0.3,by=0.1)) + 
  xlim(-0.25,0.2) +
  theme(axis.text.y = element_text(family="serif", face="bold", size=20, 
                                   angle=360, vjust = 0.5),
        legend.position = "none")


#### BETA (45,5) ----
colnames(beta455) <- c("p_mean", "psi_mean", "w_a_mean", "w_b_mean",  
                      "p_median", "psi_median", "w_a_median", "w_b_median") 
res_b455 <- as.data.frame(beta455) %>% 
  relocate(psi_median, .after=w_b_mean) %>% 
  relocate(psi_mean, .before=p_mean)
bias_b455 <- as.data.frame(t(apply(res_b455, 1, function(x){x-rep(c(0.8,0.5,0.9,0.7),2)})))

bias_b455_mean <- data.frame(param=stack(bias_b455[,1:4])[,1], 
                             gro=c(rep("psi",100),rep("p",100),rep("wa",100),rep("wb",100)))
bias_b455_mean$gro <- factor(bias_b455_mean$gro, levels=c('wb','wa','p','psi'), ordered=T)

bias_b455_med <- data.frame(param=stack(bias_b455[,5:8])[,1], 
                            gro=c(rep("psi",100),rep("p",100),rep("wa",100),rep("wb",100)))
bias_b455_med$gro <- factor(bias_b455_med$gro, levels=c('wb','wa','p','psi'), ordered=T)

# Boxplot
b455 <- ggplot() +
  geom_boxplot(data=bias_b455_mean, aes(x=param, y=gro, fill=gro), outlier.alpha = 0.6,
               alpha=0.7) + # fill="#E25E3E", col="#E25E3E"
  geom_jitter(data=bias_b455_mean, aes(x=param, y=gro, col=gro), size=0.4) + # color="#E25E3E",
  geom_vline(xintercept=0) + #xlim(-0.1,0.1) +
  theme_light(base_size = 10) + xlab(NULL) + ylab(NULL) +
  scale_fill_manual(values=Palette) + scale_colour_manual(values=Palette) + 
  scale_y_discrete(labels=c( expression(hat(w[B])),  expression(hat(w[A])),
                             expression(hat(p)), expression(hat(psi)))) +
  # scale_x_continuous(breaks=seq(-0.4, 0.3,by=0.1)) + 
  xlim(-0.25,0.2) +
  theme(axis.text.y = element_text(family="serif", face="bold", size=20, 
                                   angle=360, vjust = 0.5),
        legend.position = "none")


plot_grid(u01, u81, b455,nrow=3)
