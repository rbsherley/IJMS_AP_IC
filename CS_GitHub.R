#install.packages("Rmisc")
#install.packages("ggplot2")
require(Rmisc)
require(ggplot2)
load("CS_annual_mean.rdata")

# extract means, SEs and CIs
dat.ctv <- summarySE(dat, measurevar="Surv", groupvars=c("Island","Closure"))

# plot the raw means and approximate 95% CIs
pd <- position_dodge(0.5)
ggplot(dat.ctv, aes(x=Island, y=Surv, color=Closure)) + 
  geom_errorbar(aes(ymin=Surv-ci, ymax=Surv+ci), width=.1,position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd)
# The data suggest that any closure effect is in the same direction at the two islands.

# This model gives the same results as the "best" model (3rd row) in Table B2 of Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev)
m1<- lmerTest::lmer(Surv~Island/Closure+(1|Year),data=dat,control = lme4::lmerControl(optimizer="bobyqa"),REML=T)
summary(m1) # The main effect of closure (closure at Dassen) is significant at the 5% level
            # The interaction (the closure effect at Robben) is not significant in this model.
plot(emmeans::emmeans(m1, c("Island","Closure")))
# percentage closure effect size at Dassen (14%)
round((abs(m1@beta[3])/(m1@beta[1]+m1@beta[3]))*100,1)
# percentage closure effect size at Robben (1.2%)
round((abs(m1@beta[4])/(m1@beta[1]+m1@beta[2]+m1@beta[4]))*100,1)

# But is there any support for including the interaction effect (the island-specific closure effect)?
# The initial plots of the raw data suggest effects in the same direciton and of similar magnitude
# Use AICc-based model selection (with maximum likelihood estimation) to select the best fitting model
m1<- lmerTest::lmer(Surv~Island/Closure+(1|Year),data=dat,REML=F,control = lme4::lmerControl(optimizer="bobyqa"))
m2<- lmerTest::lmer(Surv~Island+Closure+(1|Year),data=dat,REML=F,control = lme4::lmerControl(optimizer="bobyqa"))
m3<- lmerTest::lmer(Surv~Closure+(1|Year),data=dat,REML=F,control = lme4::lmerControl(optimizer="bobyqa"))
m4<- lmerTest::lmer(Surv~Island+(1|Year),data=dat,REML=F,control = lme4::lmerControl(optimizer="bobyqa"))
m5<- lmerTest::lmer(Surv~1+(1|Year),data=dat,REML=F,control = lme4::lmerControl(optimizer="bobyqa"))

bbmle::ICtab(m1,m2,m3,m4,m5,base=T,weights=T,type="AICc")
ICtable = bbmle::ICtab(m1,m2,m3,m4,m5,base=T,weights=T,type="AICc")
# model 3 (closure only) has 72% of the AICc weight
# model 2 (island + closure) has 15% of the AICc weight and a delta AICc to the interaction model (model 1) of 2.6.
# So, keeping the interaction in the model is not well supported by the data

# refit that model using REML for final inference
m3<- lmerTest::lmer(Surv~Closure+(1|Year),data=dat,REML=T,control = lme4::lmerControl(optimizer="bobyqa"))
summary(m3)
plot(emmeans::emmeans(m3, c("Closure")))
# percentage closure effect size (7.9%) on average across the two islands
round((abs(m1@beta[2])/m1@beta[1])*100,1)

# In summary, if model selection is used, there is no support to retain the island specific closure effects in the model using the aggregated chick survival data. And the overall closure effect size is 7.9%, as opposed to the 14% effect at Dassen Island and the 1.2% effect size at Robben Island reported by Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev)

# And refit model 2 with REML for the plot
m2<- lmerTest::lmer(Surv~Closure+Island+(1|Year),data=dat,REML=T,control = lme4::lmerControl(optimizer="bobyqa"))
summary(m2)
## Covert effects to % population change using the approach in Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev)

# # # # # # # # 
#     Plot    #
# # # # # # # # 

par(mfrow=c(1,1),mar=c(3.5, 3.5, 0.5, 0.5),mgp = c(2, 0.6, 0))#,family='serif')
# Dassen from Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) selected model
plotrix::plotCI(x=0.75,y=round((0.124*-0.096472)*100,2),
li=round((0.124*(-0.096472-(1.96*0.044128)))*100,2),# note that Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) use 2 * SE for approx. 95% CIs rather than 1.96 
ui=round((0.124*(-0.096472+(1.96*0.044128)))*100,2),pch=16,xlim=c(0.6,8.4),ylim=c(-2.5,1),xaxt='n',xlab='Model',ylab="% change in pop. growth rate",col="orange")
abline(h=-1,lty=3)
abline(h=0,lty=)
# Roben from Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) selected model
plotrix::plotCI(x=1.25,y=round((0.124*-0.008454)*100,2),
li=round((0.124*(-0.008454-(1.96*0.044128)))*100,2), # note that Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) use 2 * SE for approx. 95% CIs rather than 1.96 
ui=round((0.124*(-0.008454+(1.96*0.044128)))*100,2),add=T,pch=16,col='purple')
text(print(round(ICtable$AICc[4],2)),x=1,y=-2.5,cex=0.8)

# Aggregated data, closure effect from the additive model (m2)
plotrix::plotCI(x=2,y=round((0.124*-0.05318)*100,2),
                li=round((0.124*(-0.05318-(1.96*0.01832)))*100,2),
                ui=round((0.124*(-0.05318+(1.96*0.01832)))*100,2),add=T,pch=16,col="#80164b")
text(print(round(ICtable$AICc[2],2)),x=2,y=-2.5,cex=0.8)

# Aggregated data, closure effect from best fitting model (m3)
plotrix::plotCI(x=3,y=round((0.124*-0.05297)*100,2),
                li=round((0.124*(-0.05297-(1.96*0.01759)))*100,2),
ui=round((0.124*(-0.05297+(1.96*0.01759)))*100,2),add=T,pch=16,col='#80164b')
text(print(round(ICtable$AICc[1],2)),x=3,y=-2.5,cex=0.8)

# Dassen from M3 of Sherley 2020 corrected - Year/NestID no interaction
# Effect sizes from jagsM4 in CS_JAGS_GitHub.R
plotrix::plotCI(x=3.75,y=round((0.124*-0.07623778)*100,2),
                li=round((0.124*(-0.07623778-(1.96*0.01572658)))*100,2),
                ui=round((0.124*(-0.07623778+(1.96*0.01572658)))*100,2),add=T,pch=16,col="orange")

# Robben from M3 of Sherley 2020 corrected - Year/NestID no interaction
# Effect sizes from jagsM4 in CS_JAGS_GitHub.R
plotrix::plotCI(x=4.25,y=round((0.124*-0.0789118)*100,2),
li=round((0.124*(-0.0789118-(1.96*0.01569493)))*100,2),
ui=round((0.124*(-0.0789118+(1.96*0.01569493)))*100,2),add=T,pch=16,col='purple')
text("16561.8",x=4,y=-2.5,cex=0.8) #LOOIC from jagsM3 in CS_JAGS_GitHub.R

# Dassen from Sherley 2020 corrected - Year/NestID with interaction
# Effect sizes from jagsM5 in CS_JAGS_GitHub.R
plotrix::plotCI(x=4.75,y=round((0.124*-0.09568504)*100,2),
                li=round((0.124*(-0.09568504-(1.96*0.03906169)))*100,2),
                ui=round((0.124*(-0.09568504+(1.96*0.03906169)))*100,2),add=T,pch=16,col="orange")

# Robben from Sherley 2020 corrected - Year/NestID with interaction
# Effect sizes from jagsM5 in CS_JAGS_GitHub.R
plotrix::plotCI(x=5.25,y=round((0.124*-0.09094246)*100,2),
                li=round((0.124*(-0.09094246-(1.96*0.03098658)))*100,2),
                ui=round((0.124*(-0.09094246+(1.96*0.03098658)))*100,2),add=T,pch=16,col='purple')
text("16563.7",x=5,y=-2.5,cex=0.8) #LOOIC from jagsM3_int in CS_JAGS_GitHub.R

# Dassen from Sherley 2020 corrected - Year/Island with no interaction
# Effect sizes from jagsM6 in CS_JAGS_GitHub.R
plotrix::plotCI(x=5.75,y=round((0.124*-0.05412317)*100,2),
                li=round((0.124*(-0.05412317-(1.96*0.01823597)))*100,2),
                ui=round((0.124*(-0.05412317+(1.96*0.01823597)))*100,2),add=T,pch=16,col="orange")

# Robben from Sherley 2020 corrected - Year/Island with no interaction
# Effect sizes from jagsM6 in CS_JAGS_GitHub.R
plotrix::plotCI(x=6.25,y=round((0.124*-0.05473583)*100,2),
                li=round((0.124*(-0.05473583-(1.96*0.01830139)))*100,2),
                ui=round((0.124*(-0.05473583+(1.96*0.01830139)))*100,2),add=T,pch=16,col="purple")
text("17146.9",x=6,y=-2.5,cex=0.8) #LOOIC from jagsM4 in CS_JAGS_GitHub.R

# Dassen from OLSPS 2021 (Gibbs sampling)
plotrix::plotCI(x=6.75,y=-0.942,li=-1.761,ui=-0.025,add=T,pch=16,col="orange")
# Robben from OLSPS 2021 (Gibbs sampling)
plotrix::plotCI(x=7.25,y=-0.397,li=-1.240,ui=0.620,add=T,pch=16,col='purple')

# Raw (unmodelled) mean differences
# Dassen
points(x=7.75,y=round((0.124*-0.0665567)*100,2),pch=3,col="orange",lwd=2)
# Robben
points(x=8.25,y=round((0.124*-0.07796)*100,2),pch=3,col="purple",lwd=2)

axis(1,at=1:8,labels=c("1.A,I*C,Y","2.A,I+C,Y","3.A,C,Y","4.SHE,I+C,YN","5.SHE,I*C,YN","6.SHE,I+C,YI","7.OLS,I*C,Y","8.Raw means"),cex.axis=0.7)

#EOF
