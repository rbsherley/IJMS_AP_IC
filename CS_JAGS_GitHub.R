#set folder to save models in to keep tidy
model.dir <-"Models/WC_CS/"
dir.create(model.dir,showWarnings = F,recursive=T)

#set folder to save utputs in to keep tidy
output.dir<-"Output/WC_CS/"
dir.create(output.dir,showWarnings = F,recursive=T)

### Read in data and load libraries ###
load("CS_disagg_data.rdata")
require(jagsUI)
require(rjags)
rjags::load.module('dic')
require(loo)

######################################################################################################
######################################################################################################
######              JAGS data prep:               ########################
######################################################################################################
######################################################################################################

# Fish Biomass data
### Double check values against Excel spreadsheet ###
CS$SB <- CS$Sardine_Biomass_m_tonnes
CS$AB <- CS$Anchovy_Biomass_m_tonnes
CS$S <- CS$SB-mean(CS$SB)
CS$A <- CS$AB-mean(na.omit(CS$AB))
mu.A = mean(na.omit(CS$A))
prec.A = sd(na.omit(CS$A))^-2

CS$Y <- CS$Year-2007
CS$Year <- as.factor(CS$Year)
CS$ICodeR <- as.numeric(as.factor(CS$Island)) # "DI" = 1 "RI" = 2 # For the Island Random effect in M6

CS$t <- CS$ChickDays
CS$t.cen <- 0
CS$t[CS$ChickFailure==0] <- NA
CS$t.cen[CS$ChickFailure==0] <- CS$ChickDays[CS$ChickFailure==0]
isCensored <- NA
isCensored[CS$ChickFailure==0] <-T
isCensored[CS$ChickFailure==1] <-F

### Data for prediction
new.t <- seq(0,74, by=1) # this will be used for prediction

# Organise data for censored observations
censorLimit = max(CS$ChickDays[CS$ChickFailure==1])
censorLimitVec = rep(censorLimit, length(CS$t))

CS$t.cen1 <- CS$t.cen
CS$t.cen1[CS$t.cen != 0] <- CS$t.cen[CS$t.cen != 0]+1
CS$t.cen1[CS$t.cen == 0] <- (CS$ChickDays[CS$t.cen == 0]+1)

# Initial values
# intial values of censored data:
tInit = rep(NA, length(CS$t.cen))
censorLimit = max(CS$ChickDays)+2 # needed to be set like this to get this dataset to run
censorLimitVec = rep(censorLimit, length(CS$t))
tInit[isCensored] = censorLimitVec[isCensored]

######################################################################################################
######################################################################################################
######     Model with Year/NestID random effect and island + closure main effect      ################
######################################################################################################
######################################################################################################


# Parameters monitored
params <- c("alpha","beta.island","beta.close","beta.sard","beta.anch",
            "S.DO","S.DC","S.RO","S.RC", "sigma.re","sigma","sigma.gamma","loglik")

## JAGS model:

cat("
    model
    {
    ##############################
    # Likelihood
    ##############################
    for (i in 1:N) { 
    anch[i] ~ dnorm(mu.A, prec.A)
    isCensored[i] ~ dinterval(t[i],t.cen[i])
    t[i] ~ dlnorm(mu[i], tau)
    mu[i] <- b.yn[year[i],nest[i]] + beta.island*island[i] + beta.close*close[i] + beta.sard*sard[i] + beta.anch*anch[i] 
    # this line calculates the loglikelihood, used later to calculate the loo
    k[i] <- ifelse(isCensored[i], t.cen[i], t[i])
    loglik[i] <- log(ifelse(isCensored[i], 1 - plnorm(k[i], mu[i], tau), dlnorm(k[i], mu[i], tau)))
    } 
    
    ##############################
    # Priors and constraints
    ##############################
    alpha ~ dnorm(0.0, 0.000001)
    beta.island ~ dnorm(0.0, 0.000001)
    beta.close ~ dnorm(0.0, 0.000001)
    beta.sard ~ dnorm(0.0, 0.000001)
    beta.anch ~ dnorm(0.0, 0.000001)
    tau ~ dgamma(0.001, 0.001)
    sigma <- sqrt(1/tau)
   
    
    for (y in 1:Y) {
    for (c in 1:C) {
     b.yn[y,c] ~ dnorm(b.y[y], tau.re) # Nested random effect
    }
      b.y[y] ~ dnorm(alpha,tau.gamma)
  }

  num ~ dnorm(0,0.0016)
  denom ~ dnorm(0,1)
  sigma.re <- abs(num/denom) # SD hyperparameter for random month within year effect
  tau.re <- 1/(sigma.re*sigma.re)
    
  num.gamma ~ dnorm(0,0.0016)
  denom.gamma ~ dnorm(0,1)
  sigma.gamma <- abs(num.gamma/denom.gamma) # SD hyperparameter for random month within year effect
  tau.gamma <- 1/(sigma.gamma*sigma.gamma)

    
    ##############################
    # Derived parameters
    ##############################
    S.DO <-   1-phi((log(74) - alpha)*sqrt(tau))
    S.DC <-   1-phi((log(74) - (alpha + beta.close))*sqrt(tau))  
    S.RO <-   1-phi((log(74) - (alpha + beta.island))*sqrt(tau))
    S.RC <-   1-phi((log(74) - (alpha + beta.close + beta.island))*sqrt(tau)) #+ beta.inter
    #bC.bI <- beta.close + beta.inter         # Effect for Island/Closure
    } # end model
    ",file =paste0(model.dir,"/CS_WC_YN.jags")
)

jags.data <- list(t = CS$t, t.cen = CS$t.cen1, isCensored=as.numeric(isCensored), N = dim(CS)[1],
                  island=CS$IslandCode, close=CS$Closure, mu.A = mu.A, prec.A = prec.A,anch=CS$A, sard=CS$S,
                  year=CS$Y, nest=CS$NestCode, Y=length(unique(CS$Y)), C=length(unique(CS$NestCode)))

inits <- function() list(t=tInit, alpha = rnorm(1,0), beta.close = rnorm(1,0), beta.island = rnorm(1,0), beta.sard = rnorm(1,0), beta.anch= rnorm(1,0),
                         num=runif(1,0,25),denom=runif(1,0,1),
                         num.gamma=runif(1,0,25),denom.gamma=runif(1,0,1))

# MCMC settings:
ni <- 120000
nt <- 10
nb <- 20000
nc <- 3

## Call JAGS from R 
cs.mod <- jags(jags.data, parallel=T,inits, params, paste0(model.dir,"CS_WC_YN.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=T)

print(cs.mod, digits = 3)

## Summarize posteriors
jagsM4 <- cs.mod
write.csv(file = paste0(output.dir,"jagsM4.csv"), x = jagsM4$summary)
rm(cs.mod)
save(jagsM4, file = paste0(output.dir,"jagsM4.rda"))


######################################################################################################
######################################################################################################
######      Model with Year/NestID random effect and island x closure interaction  ###################
######################################################################################################
######################################################################################################

# Parameters monitored
params <- c("alpha","beta.island","beta.close","beta.inter","bC.bI","beta.sard","beta.anch",
            "S.DO","S.DC","S.RO","S.RC", "sigma.re","sigma","sigma.gamma","loglik")

## JAGS model:

cat("
    model
    {
    ##############################
    # Likelihood
    ##############################
    for (i in 1:N) { 
    anch[i] ~ dnorm(mu.A, prec.A)
    isCensored[i] ~ dinterval(t[i],t.cen[i])
    t[i] ~ dlnorm(mu[i], tau)
    mu[i] <- b.yn[year[i],nest[i]] + beta.island*island[i] + beta.close*close[i] + beta.sard*sard[i] + beta.anch*anch[i] + beta.inter*island[i]*close[i]
    # this line calculates the loglikelihood, used later to calculate the loo
    k[i] <- ifelse(isCensored[i], t.cen[i], t[i])
    loglik[i] <- log(ifelse(isCensored[i], 1 - plnorm(k[i], mu[i], tau), dlnorm(k[i], mu[i], tau)))
    } 
    
    ##############################
    # Priors and constraints
    ##############################
    alpha ~ dnorm(0.0, 0.000001)
    beta.island ~ dnorm(0.0, 0.000001)
    beta.close ~ dnorm(0.0, 0.000001)
    beta.sard ~ dnorm(0.0, 0.000001)
    beta.anch ~ dnorm(0.0, 0.000001)
    beta.inter ~ dnorm(0, 0.000001)
    tau ~ dgamma(0.001, 0.001)
    sigma <- sqrt(1/tau)
   
    
    for (y in 1:Y) {
    for (c in 1:C) {
     b.yn[y,c] ~ dnorm(b.y[y], tau.re) # Nested random effect
    }
      b.y[y] ~ dnorm(alpha,tau.gamma)
  }

  num ~ dnorm(0,0.0016)
  denom ~ dnorm(0,1)
  sigma.re <- abs(num/denom) # SD hyperparameter for random month within year effect
  tau.re <- 1/(sigma.re*sigma.re)
    
  num.gamma ~ dnorm(0,0.0016)
  denom.gamma ~ dnorm(0,1)
  sigma.gamma <- abs(num.gamma/denom.gamma) # SD hyperparameter for random month within year effect
  tau.gamma <- 1/(sigma.gamma*sigma.gamma)

    
    ##############################
    # Derived parameters
    ##############################
    S.DO <-   1-phi((log(74) - alpha)*sqrt(tau))
    S.DC <-   1-phi((log(74) - (alpha + beta.close))*sqrt(tau))  
    S.RO <-   1-phi((log(74) - (alpha + beta.island))*sqrt(tau))
    S.RC <-   1-phi((log(74) - (alpha + beta.close + beta.island))*sqrt(tau)) #+ beta.inter
    bC.bI <- beta.close + beta.inter         # Effect for Island/Closure
    } # end model
    ",file =paste0(model.dir,"/CS_WC_YN_I.jags")
)

jags.data <- list(t = CS$t, t.cen = CS$t.cen1, isCensored=as.numeric(isCensored), N = dim(CS)[1],
                  island=CS$IslandCode, close=CS$Closure, mu.A = mu.A, prec.A = prec.A,anch=CS$A, sard=CS$S,
                  year=CS$Y, nest=CS$NestCode, Y=length(unique(CS$Y)), C=length(unique(CS$NestCode)))

inits <- function() list(t=tInit, alpha = rnorm(1,0), beta.close = rnorm(1,0), beta.island = rnorm(1,0),beta.inter = rnorm(1,0), beta.sard = rnorm(1,0), beta.anch= rnorm(1,0), num=runif(1,0,25),denom=runif(1,0,1),num.gamma=runif(1,0,25),denom.gamma=runif(1,0,1))


## Call JAGS from R 
cs.mod <- jags(jags.data, parallel=T,inits, params, paste0(model.dir,"CS_WC_YN_I.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=T)

print(cs.mod, digits = 3)

## Summarize posteriors
jagsM5 <- cs.mod
write.csv(file = paste0(output.dir,"jagsM5.csv"), x = jagsM5$summary)
rm(cs.mod)
save(jagsM5, file = paste0(output.dir,"jagsM5.rda"))


######################################################################################################
######################################################################################################
######     Model with Year/Island random effect and island + closure main effect      ################
######################################################################################################
######################################################################################################


# Parameters monitored
params <- c("alpha","beta.island","beta.close","beta.sard","beta.anch",
            "S.DO","S.DC","S.RO","S.RC", "sigma.re","sigma","sigma.gamma","loglik")

## JAGS model:

cat("
    model
    {
    ##############################
    # Likelihood
    ##############################
    for (i in 1:N) { 
    anch[i] ~ dnorm(mu.A, prec.A)
    isCensored[i] ~ dinterval(t[i],t.cen[i])
    t[i] ~ dlnorm(mu[i], tau)
    mu[i] <- b.yn[year[i],islandR[i]] + beta.island*island[i] + beta.close*close[i] + beta.sard*sard[i] + beta.anch*anch[i] 
    # this line calculates the loglikelihood, used later to calculate the loo
    k[i] <- ifelse(isCensored[i], t.cen[i], t[i])
    loglik[i] <- log(ifelse(isCensored[i], 1 - plnorm(k[i], mu[i], tau), dlnorm(k[i], mu[i], tau)))
    } 
    
    ##############################
    # Priors and constraints
    ##############################
    alpha ~ dnorm(0.0, 0.000001)
    beta.island ~ dnorm(0.0, 0.000001)
    beta.close ~ dnorm(0.0, 0.000001)
    beta.sard ~ dnorm(0.0, 0.000001)
    beta.anch ~ dnorm(0.0, 0.000001)
    tau ~ dgamma(0.001, 0.001)
    sigma <- sqrt(1/tau)
   
    
    for (y in 1:Y) {
    for (c in 1:C) {
     b.yn[y,c] ~ dnorm(b.y[y], tau.re) # Nested random effect
    }
      b.y[y] ~ dnorm(alpha,tau.gamma)
  }

  num ~ dnorm(0,0.0016)
  denom ~ dnorm(0,1)
  sigma.re <- abs(num/denom) # SD hyperparameter for random month within year effect
  tau.re <- 1/(sigma.re*sigma.re)
    
  num.gamma ~ dnorm(0,0.0016)
  denom.gamma ~ dnorm(0,1)
  sigma.gamma <- abs(num.gamma/denom.gamma) # SD hyperparameter for random month within year effect
  tau.gamma <- 1/(sigma.gamma*sigma.gamma)

    
    ##############################
    # Derived parameters
    ##############################
    S.DO <-   1-phi((log(74) - alpha)*sqrt(tau))
    S.DC <-   1-phi((log(74) - (alpha + beta.close))*sqrt(tau))  
    S.RO <-   1-phi((log(74) - (alpha + beta.island))*sqrt(tau))
    S.RC <-   1-phi((log(74) - (alpha + beta.close + beta.island))*sqrt(tau)) #+ beta.inter
    #bC.bI <- beta.close + beta.inter         # Effect for Island/Closure
    } # end model
    ",file =paste0(model.dir,"/CS_WC_YI.jags")
)

jags.data <- list(t = CS$t, t.cen = CS$t.cen1, isCensored=as.numeric(isCensored), N = dim(CS)[1],
                  island=CS$IslandCode, close=CS$Closure, mu.A = mu.A, prec.A = prec.A,anch=CS$A, sard=CS$S,
                  year=CS$Y, islandR=CS$ICodeR, Y=length(unique(CS$Y)), C=length(unique(CS$ICodeR)))

inits <- function() list(t=tInit, alpha = rnorm(1,0), beta.close = rnorm(1,0), beta.island = rnorm(1,0), beta.sard = rnorm(1,0), beta.anch= rnorm(1,0),
                         num=runif(1,0,25),denom=runif(1,0,1),
                         num.gamma=runif(1,0,25),denom.gamma=runif(1,0,1))


## Call JAGS from R 
cs.mod <- jags(jags.data, parallel=T,inits, params, paste0(model.dir,"CS_WC_YI.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, DIC=T)

print(cs.mod, digits = 3)

## Summarize posteriors
jagsM6 <- cs.mod
write.csv(file = paste0(output.dir,"jagsM6.csv"), x = jagsM6$summary)
rm(cs.mod)
save(jagsM6, file = paste0(output.dir,"jagsM6.rda"))



######################################################################################################
######################################################################################################
######################                      Model selection:                    ######################
######################################################################################################
######################################################################################################

#load("Output/WC_CS/jagsM4.rda")
log_lik = jagsM4$sims.list$loglik
reff = loo::relative_eff(log_lik, chain_id=rep(1,nrow(log_lik)))
M4.LOOIC = loo::loo(log_lik,r_eff=reff)
rm(reff)
rm(log_lik)

#load("Output/WC_CS/jagsM5.rda")
log_lik = jagsM5$sims.list$loglik
reff = loo::relative_eff(log_lik, chain_id=rep(1,nrow(log_lik)))
M5.LOOIC = loo::loo(log_lik, r_eff=reff)
rm(reff)
rm(log_lik)

#load("Output/WC_CS/jagsM6_int.rda")
log_lik = jagsM6$sims.list$loglik
reff = loo::relative_eff(log_lik, chain_id=rep(1,nrow(log_lik)))
M6.LOOIC = loo::loo(log_lik,r_eff=reff)
rm(reff)
rm(log_lik)

M4.LOOIC
M5.LOOIC
M6.LOOIC
min(M4.LOOIC$estimates[3],M5.LOOIC$estimates[3],M6.LOOIC$estimates[3])
max(M4.LOOIC$estimates[3],M5.LOOIC$estimates[3],M6.LOOIC$estimates[3])-min(M4.LOOIC$estimates[3],M5.LOOIC$estimates[3],M6.LOOIC$estimates[3])


######################################################################################################
######################################################################################################
############             Effect sizes for the plot in CS_GitHub:                             #########
######################################################################################################
######################################################################################################

-mean(jagsM4$sims.list$S.RC-jagsM4$sims.list$S.RO)
sd(jagsM4$sims.list$S.RC-jagsM4$sims.list$S.RO)
-mean(jagsM4$sims.list$S.DC-jagsM4$sims.list$S.DO)
sd(jagsM4$sims.list$S.DC-jagsM4$sims.list$S.DO)

-mean(jagsM5$sims.list$S.RC-jagsM5$sims.list$S.RO)
sd(jagsM5$sims.list$S.RC-jagsM5$sims.list$S.RO)
-mean(jagsM5$sims.list$S.DC-jagsM5$sims.list$S.DO)
sd(jagsM5$sims.list$S.DC-jagsM5$sims.list$S.DO)

-mean(jagsM6$sims.list$S.RC-jagsM6$sims.list$S.RO)
sd(jagsM6$sims.list$S.RC-jagsM6$sims.list$S.RO)
-mean(jagsM6$sims.list$S.DC-jagsM6$sims.list$S.DO)
sd(jagsM6$sims.list$S.DC-jagsM6$sims.list$S.DO)

