# Using informative prior to account for identifiability issues 
# in occupancy models with identification errors
#---------------------------------------------------------------#

library(purrr) ; library(dplyr) ; library(ggplot2) ; 
library(cowplot) ; library(ggridges)

## 1. Classical Estimation from the identification layer ----

### a. Frequentist inferences with 12 visits ----
# Initial settings to simulate data
N <- 30       # sites
J <- 12  # visits
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


#### i. Determining initial values for Maximum Likelihood Estimation ----
set.seed(123)
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

# defining initial values
psi0 <- length(sumI[which(sumI!=0)])/N
p0 <- sqrt(w_a*p)
w_a0 <- sqrt(w_a*p)
w_b0 <- 1- mean(sumI[which(Z==0)])/J


#### ii. Estimation from the Minimum of the negative log-Likelihood ----
my_optim <- optim(par=c(psi0,p0,w_a0, w_b0), fn=NegLik, hessian=T)

print(paste("PSI=", plogis(my_optim$par[1]),
            "p=", plogis(my_optim$par[2]),
            "w_a=", plogis(my_optim$par[3]),
            "w_b=", plogis(my_optim$par[4])))


### iii. Simulations for estimates ----
k = 1000
res_simu12 <- matrix(nrow=k, ncol=4)
# res_simu36 <- matrix(nrow=k, ncol=4)

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
  res_simu12[i,] <- plogis(my_optimS$par)
  # res_simu36[i,] <- plogis(my_optimS$par)
}

colnames(res_simu12) <- c("psi", "p", "w_a", "w_b")
# head(res_simu12)
colnames(res_simu36) <- c("psi", "p", "w_a", "w_b")

#### vi. Plot ----
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

### b. Frequentist inferences with 36 visits ----
# Initial settings to simulate data
J <- 36  # visits

#### i. Determining initial values for Maximum Likelihood Estimation ----
set.seed(123)
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

# defining initial values
psi0 <- length(sumI[which(sumI!=0)])/N
p0 <- sqrt(w_a*p)
w_a0 <- sqrt(w_a*p)
w_b0 <- 1- mean(sumI[which(Z==0)])/J


#### ii. Estimation from the Minimum of the negative log-Likelihood ----
my_optim <- optim(par=c(psi0,p0,w_a0, w_b0), fn=NegLik, hessian=T)

print(paste("PSI=", plogis(my_optim$par[1]),
            "p=", plogis(my_optim$par[2]),
            "w_a=", plogis(my_optim$par[3]),
            "w_b=", plogis(my_optim$par[4])))


#### vi. Simulations for estimates ----
k = 1000
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
  res_simu36[i,] <- plogis(my_optimS$par)
}
colnames(res_simu36) <- c("psi", "p", "w_a", "w_b")

#### v. Plot ----
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

### c. plot inferences for 12 and 36 visits
J1236 <- plot_grid(J12,J36, byrow=FALSE, ncol=2)
J1236

#### Standard Deviation ----
sqrt(sum((res_simu12[,1]-mean(res_simu12[,1]))^2)/(k-1))
sqrt(sum((res_simu36[,1]-mean(res_simu36[,1]))^2)/(k-1))

#### Biais ----
mean(res_simu36[,2])-0.5 ; mean(res_simu36[,3])-0.9


### d. Compare to the original model (without identification layer)
# Negative Log Likelihood of the original model
NegLik_org <- function(occ){
  psi <- plogis(occ[1])
  p <- plogis(occ[2])
  cp <- (p^sumI) * ((1-p)^(J-sumI))
  loglik <- log((cp*psi) + (1-avecD)*(1-psi))
  return(-sum(loglik))
}

k = 1000
res_org <- matrix(nrow=k, ncol=2)

set.seed(123)
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
  # Here we consider the identification matrix W like it was the detection matrix Y
  avecD <- vector(length=N)
  for(j in 1:nrow(W)){
    ifelse(sumI[j]==0, avecD[j]<-0, avecD[j]<-1)
  }
  
  # defining initial values
  psi0_org <- length(sumI[which(sumI!=0)])/N
  p0_org <- sum(sumI)/(sum(Y[,1])*J) # == mean(sumI[which(Z==1)]/J)
  
  my_mle <- optim(fn=NegLik_org, par=c(psi0_org,p0_org), hessian=T)
  res_org[i,] <- plogis(my_mle$par)
}

# plot
denspsi_org <- density(res_org[,1], bw=0.02)
densp_org <- density(res_org[,2], bw=0.02)

psi_org <- ggplot()+
  geom_histogram(aes(res_org[,1], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray", binwidth = 0.05) + 
  geom_density(aes(x = res_org[,1], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2,linewidth=0.8, bw=0.02)+ 
  geom_vline(xintercept=psi, col="darkred") +
  geom_point(aes(x=psi, y=mean(denspsi_org$y[which(round(denspsi_org$x, digit=2)==psi)]/max(denspsi_org$y))), 
             col="darkred", size=3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1)) +
  ylab(NULL) + xlab(NULL) + theme_light(base_size = 10) 

p_org <- ggplot()+
  geom_histogram(aes(res_org[,2], y=after_stat(ncount)), fill="#ccc5b9", col="darkgray",binwidth = 0.05) +
  geom_density(aes(x = res_org[,2], y=after_stat(scaled)), fill="#354f52", col="#354f52",alpha=0.2, linewidth=0.8, bw=0.02)+ 
  geom_vline(xintercept=p, col="darkred") +
  geom_point(aes(x=p, y=mean(densp_org$y[which(round(densp_org$x, digit=2)==p)]/max(densp_org$y))), 
             col="darkred", size=3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits=c(0,1)) +
  ylab(NULL) + xlab(NULL) + theme_light(base_size = 10)

plot_grid(psi_org,p_org,ncol=1,
                        label_x=0, label_y=0, label_size = 8, hjust = -0.5, vjust = -0.5)

mean(res_org[,1]) ; quantile(res_org[,1],0.1) ; quantile(res_org[,1],0.9)
