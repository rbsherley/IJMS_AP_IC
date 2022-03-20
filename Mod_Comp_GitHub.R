#install.packages("jagsUI")
#install.packages("rjags")
#install.packages("plotrix")
require(plotrix)
require(jagsUI)
require(rjags)

load("MetaAnalysis.rdata")
dat$icode <- as.numeric(as.factor(dat$Island))-1

# Parameters to trace
params <- c("mu.D","mu.R","beta.island","mu.D.cond","mu.R.cond","beta.island.cond", "sigma","sigma.cond")

# data for JAGS:
jags.data <- list(surv = dat$Surv, se.surv = dat$SurvSE, 
                  cond = dat$Cond, se.surv.cond = dat$CondSE,
                  n = length(dat$Surv), icode = dat$icode)

# Initials:
inits <- function() list(mu.D = rnorm(1,0), beta.island = rnorm(1,0),
                         mu.D.cond = rnorm(1,0), beta.island.cond = rnorm(1,0))

# MCMC settings:
ni <- 120000
nt <- 5
nb <- 20000
nc <- 3


## JAGS model:

cat("
    model {
  # Priors
  mu.D ~ dnorm(0, 0.0001)          # mean survival effect size at Dassen
  beta.island ~ dnorm(0, 0.0001)  # difference Robben to Dassen in survival
  sigma ~ dunif(0, 100)
  tau <- pow(sigma, -2)  	# Residual Standard Error
  
  mu.D.cond ~ dnorm(0, 0.0001)          # mean condition effect size at Dassen
  beta.island.cond ~ dnorm(0, 0.0001)  # difference Robben to Dassen in condition
  sigma.cond ~ dunif(0, 100)
  tau.cond <- pow(sigma.cond, -2)  	# Residual Standard Error
  
  # Likelihood
  for (i in 1:n)  # for each of 10 model estimates
  {
    surv[i] ~ dnorm(mu[i], tau.error[i])    
    mu[i] ~ dnorm(mu.a[i], tau)
    mu.a[i] <- mu.D + beta.island*icode[i]  # linear predictor 
    tau.error[i] <- pow(se.surv[i], -2)
    
    
    cond[i] ~ dnorm(mu.cond[i], tau.error.cond[i])    
    mu.cond[i] ~ dnorm(mu.a.cond[i], tau.cond)
    mu.a.cond[i] <- mu.D.cond + beta.island.cond*icode[i]  # linear predictor 
    tau.error.cond[i] <- pow(se.surv.cond[i], -2)
  }
  
  # Derived quantities
  mu.R <- mu.D + beta.island         # mean effect size at Robben
  mu.R.cond <- mu.D.cond + beta.island.cond
  
} # end model
    ",file ="MA.jags")


## Call JAGS from R 
ma.mod <- jags(jags.data, parallel=T,inits, params, "MA.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=T)

print(ma.mod, digits = 3)



pdf(file = "Figure1.pdf", width = 10, height = 10)
#quartz(width = 10, height = 10)
par(mfrow=c(2,1),mar=
      c(2.2, 3.2, 1.2, 0.4),mgp = c(2, 0.6, 0),cex=1.2)#,family='serif')
##############################
####### Chick Survival ######
##############################
plotrix::plotCI(x=0.75,y=1.20,ui=2.29,li=0.10,pch=15,col="white",xlim=c(0.75,5.25),ylim=c(-1.1,2.5),xaxt='n',lwd=2,cex=1.4,
                xlab='',ylab="% improvement in pop. growth")
# Meta analysis results for Dassen survival effect size
polygon(c((seq(0.75,2.75,0.25)), rev((seq(0.75,2.75,0.25)))), c(rep(ma.mod$q2.5$mu.D,9), rev(rep(ma.mod$q97.5$mu.D,9))),
        col=adjustcolor("orange",alpha.f=0.4), border = NA)
lines(x=seq(0.75,2.75,0.25),y=rep(ma.mod$mean$mu.D,9),col='orange',lwd=2)

# Meta analysis results for Robben survival effect size
polygon(c((seq(3.25,5.25,0.25)), rev((seq(3.25,5.25,0.25)))), c(rep(ma.mod$q2.5$mu.R,9), rev(rep(ma.mod$q97.5$mu.R,9))),
        col=adjustcolor("purple",alpha.f=0.4), border = NA)
lines(x=seq(3.25,5.25,0.25),y=rep(ma.mod$mean$mu.R,9),col='purple',lwd=2)

abline(v=3,lty=1)
abline(h=1,lty=3)
abline(h=0,lty=)

# Dassen survival from Table 1 of Butterworth and Ross-Gillespie (2021; SWG-PEl/39rev) (https://zivahub.uct.ac.za/articles/report/Updated_analysis_of_results_from_data_arising_from_the_Island_Closure_Experiment/15073404)
# Aggregated
plotCI(x=0.75,y=1.20,ui=2.29,li=0.10,pch=15,col="black",add=T,lwd=2,cex=1.4)
# Disggregated
plotCI(x=1.25,y=0.84,ui=1.78,li=-0.10,add=T,pch=16,col="grey60",lwd=2,cex=1.4)
# Dassen from best fitting model with the interaction (5.SHE,I*C,YN) in SM of Sydeman et al. 2021 - Year/NestID random, island*closure fixed effects
plotCI(x=1.75,y=round((0.124*0.09568504)*100,2),
                ui=round((0.124*(0.09568504+(1.96*0.03906169)))*100,2),
                li=round((0.124*(0.09568504-(1.96*0.03906169)))*100,2),add=T,pch=16,col="orange",lwd=2,cex=1.4)

# Dassen survival from best fitting model (D3) in Table SM-3b of Butterworth and Ross-Gillespie 2022
plotCI(x=2.25,y=1.17,ui=1.76,li=0.58,add=T,pch=16,col="grey60",lwd=2,cex=1.4)
# Dassen survival from 2nd best fitting model in Table SM-3b of Butterworth and Ross-Gillespie 2022 (where the Year SD is fixed)
plotCI(x=2.75,y=1.14,ui=2.08,li=0.20,add=T,pch=16,col="grey60",lwd=2,cex=1.4)


# Robben survival from Table 1 of Butterworth and Ross-Gillespie 2021
# Aggregated
plotCI(x=3.25,y=0.10,ui=1.20,li=-0.99,add=T,pch=15,col="black",lwd=2,cex=1.4)
# Disggregated
plotCI(x=3.75,y=0.51,ui=1.45,li=-0.42,add=T,pch=16,col="grey60",lwd=2,cex=1.4)

# Robben from best fitting model with the interaction (5.SHE,I*C,YN) in SM of Sydeman et al. 2021 - Year/NestID random, island*closure fixed effects
plotCI(x=4.25,y=round((0.124*0.09094246)*100,2),
                ui=round((0.124*(0.09094246+(1.96*0.03098658)))*100,2),
                li=round((0.124*(0.09094246-(1.96*0.03098658)))*100,2),add=T,pch=16,col='purple',lwd=2,cex=1.4)

# Robben survival from best fitting model in Table SM-3b of Butterworth and Ross-Gillespie 2022
plotCI(x=4.75,y=1.18,ui=1.64,li=0.72,add=T,pch=16,col="grey60",lwd=2,cex=1.4)
# Robben survival from 2nd best fitting model in Table SM-3b of Butterworth and Ross-Gillespie 2022 (where the Year SD is fixed)
plotCI(x=5.25,y=0.81,ui=1.83,li=-0.21,add=T,pch=16,col="grey60",lwd=2,cex=1.4)

axis(1,at=seq(0.75,5.25,0.5),
     labels=c("RGB21","RGB21","SYD21","BRG22","BRG22",
              "RGB21","RGB21","SYD21","BRG22","BRG22"),cex.axis=0.8)
text(x=seq(0.75,5.25,0.5),y=rep(-1.1,10),c("A","D","D5","D3","D3b","A","D","D5","D3","D3b"))
text(x=1.75,y=2.5,"Dassen Island")
text(x=4.25,y=2.5,"Robben Island")
mtext("A", at = 2.8,side=2, line=2,cex=1.5,las=1)
mtext("Chick Survival",side=3, line=0.1,cex=1.4)

##############################
####### Chick Condition ######
##############################
plotCI(x=0.75,y=-0.06,ui=0.77,li=-0.90,pch=15,col="white",xlim=c(0.75,5.25),ylim=c(-1.1,2),xaxt='n',lwd=2,cex=1.4,
                xlab='',ylab="% improvement in pop. growth")
# Meta analysis results for Dassen survival effect size
polygon(c((seq(0.75,2.75,0.25)), rev((seq(0.75,2.75,0.25)))), c(rep(ma.mod$q2.5$mu.D.cond,9), rev(rep(ma.mod$q97.5$mu.D.cond,9))),
        col=adjustcolor("orange",alpha.f=0.4), border = NA)
lines(x=seq(0.75,2.75,0.25),y=rep(ma.mod$mean$mu.D.cond,9),col='orange',lwd=2)

# Meta analysis results for Robben survival effect size
polygon(c((seq(3.25,5.25,0.25)), rev((seq(3.25,5.25,0.25)))), c(rep(ma.mod$q2.5$mu.R.cond,9), rev(rep(ma.mod$q97.5$mu.R.cond,9))),
        col=adjustcolor("purple",alpha.f=0.4), border = NA)
lines(x=seq(3.25,5.25,0.25),y=rep(ma.mod$mean$mu.R.cond,9),col='purple',lwd=2)

abline(v=3,lty=1)
abline(h=1,lty=3)
abline(h=0,lty=)

# Dassen condition from Table 1 of Butterworth and Ross-Gillespie (2021; SWG-PEl/39rev) (https://zivahub.uct.ac.za/articles/report/Updated_analysis_of_results_from_data_arising_from_the_Island_Closure_Experiment/15073404)
# Aggregated
plotCI(x=0.75,y=-0.06,ui=0.77,li=-0.90,pch=15,col="black",add=T,lwd=2,cex=1.4)
# Disggregated
plotCI(x=1.25,y=0.20,ui=1.05,li=-0.64,add=T,pch=16,col="grey60",lwd=2,cex=1.4)

# Dassen condition from best fitting model 3.D,YM in Figure S1 of Sydeman et al. 2021 - Year/Month random, island*closure fixed effects
plotCI(x=1.75,y=round((0.108*-0.0004112)*100,2),
                ui=round((0.108*(-0.0004112+(1.96*0.02184)))*100,2),
                li=round((0.108*(-0.0004112-(1.96*0.02184)))*100,2),add=T,pch=16,col="orange",lwd=2,cex=1.4)

# Dassen condition from best fitting model in Table SM-3a of Butterworth and Ross-Gillespie 2022
plotCI(x=2.25,y=0.27,ui=0.94,li=-0.39,add=T,pch=16,col="grey60",lwd=2,cex=1.4)

# Dassen condition from model D2 in Table SM-3a of Butterworth and Ross-Gillespie 2022 (which is the model suggested by the 2020 panel)
plotCI(x=2.75,y=0.21,ui=0.95,li=-0.54,add=T,pch=16,col="grey60",lwd=2,cex=1.4)


# Robben condition from Table 1 of Butterworth and Ross-Gillespie (2021; SWG-PEl/39rev) (https://zivahub.uct.ac.za/articles/report/Updated_analysis_of_results_from_data_arising_from_the_Island_Closure_Experiment/15073404)
# Aggregated
plotCI(x=3.25,y=0.92,ui=1.75,li=0.09,pch=16,col="black",add=T,lwd=2,cex=1.4)
# Disggregated
plotCI(x=3.75,y=0.92,ui=1.76,li=0.07,add=T,pch=16,col="grey60",lwd=2,cex=1.4)
# Robben condition from best fitting model 3.D,YM in Figure S1 of Sydeman et al. 2021 - Year/Month random, island*closure fixed effects
plotCI(x=4.25,y=round((0.108*0.09776)*100,2),
                ui=round((0.108*(0.09776+(1.96*0.02298)))*100,2),
                li=round((0.108*(0.09776-(1.96*0.02298)))*100,2),add=T,pch=16,col='purple',lwd=2,cex=1.4)

# Robben condition from best fitting model (D3) in Table SM-3a of Butterworth and Ross-Gillespie 2022
plotCI(x=4.75,y=0.73,ui=1.51,li=-0.04,add=T,pch=16,col="grey60",lwd=2,cex=1.4)


# Robben condition from model D2 in Table SM-3a of Butterworth and Ross-Gillespie 2022 (which is the model suggested by the 2020 panel)
plotCI(x=5.25,y=0.92,ui=1.67,li=0.17,add=T,pch=16,col="grey60",lwd=2,cex=1.4)

mtext("B", at = 2.3,side=2, line=2,cex=1.5,las=1)
axis(1,at=seq(0.75,5.25,0.5),
     labels=c("RGB21","RGB21","SYD21","BRG22","BRG22",
              "RGB21","RGB21","SYD21","BRG22","BRG22"),cex.axis=0.8)
text(x=seq(0.75,5.25,0.5),y=rep(-1.1,10),c("A","D","D3","D3","D2","A","D","D3","D3","D2"))
text(x=1.75,y=2.5,"Dassen Island")
text(x=4.25,y=2.5,"Robben Island")
mtext("Chick Condition",side=3, line=0.1,cex=1.4)
dev.off()
