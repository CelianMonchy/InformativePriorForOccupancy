## 2. Addressing identifiability issues with a constrained optimisation wa > 1-wb ----

library(purrr) ; library(dplyr) ; library(ggplot2) ; 
library(cowplot) ; library(ggridges)

# Initial settings to simulate data
N <- 30       # sites
J <- 36  # visits
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


# constraints declaration for constrOprim() function
ui <- c(0,0,1,1)
ci <- 1

k = 1000 # number of simulations
res_val <- matrix(nrow=k, ncol=6) # results with constraint
res_sansval <- matrix(nrow=k, ncol=6) # results without constraint

set.seed(123)
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
