require(bbmle)
require(lme4)
load("CC_disagg_data.rdata")

Cond$Close <- as.factor(Cond$Closure)
Cond$Island <- as.factor(Cond$Island)
Cond$Close <- relevel(Cond$Close,ref="O")

M1 <- lmerTest::lmer(Condition~Island/Close+(1|Year/Month),REML=T, data=Cond,control=lmerControl(optimizer="bobyqa"))
M2 <- lmerTest::lmer(Condition~Island/Close+(1|Year/Island),REML=T, data=Cond,control=lmerControl(optimizer="bobyqa"))
ICtab(M1,M2,base=T,type="AICc",weights=T) # The model with Year/Month is far better supported than the model with Year/Island
summary(M1)
summary(M2)
# The closure effect at Robben Island is statistically significant in both models

# # # # # # # # 
#     Plot    #
# # # # # # # # 

par(mfrow=c(1,1),mar=c(3.5, 3.5, 0.5, 0.5),mgp = c(2, 0.6, 0))#,family='serif')
# Dassen Aggregated from Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) selected model
# note that Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) use 2 * SE for approx. 95% CIs rather than 1.96 
plotrix::plotCI(x=0.75,y=0.06,li=-0.77,ui=0.90,pch=16,xlim=c(0.6,5.4),ylim=c(-2.5,1),xaxt='n',xlab='Model',ylab="% change in pop. growth rate",col="orange")
abline(h=-1,lty=3)
abline(h=0,lty=)
text("0.039",x=0.75,y=-2.5,cex=0.8)
# Robben Aggregated from Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) selected model
# note that Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) use 2 * SE for approx. 95% CIs rather than 1.96 
plotrix::plotCI(x=1.25,y=-0.92,li=-1.75, ui=-0.09,add=T,pch=16,col='purple')
text("0.038",x=1.25,y=-2.5,cex=0.8)

# Dassen Disggregated from Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) selected model
plotrix::plotCI(x=1.75,y=-0.20,li=-1.05,ui=0.64,add=T,pch=16,col="orange")
text("0.039",x=1.75,y=-2.5,cex=0.8)
# Robben Disggregated from Ross-Gillespie and Butterworth (2021; SWG-PEl/39rev) selected model
plotrix::plotCI(x=2.25,y=-0.92,li=-1.76, ui=-0.07,add=T,pch=16,col='purple')
text("0.039",x=2.25,y=-2.5,cex=0.8)

# Dissagregated data, closure effect from Year/Month model
# Dassen
plotrix::plotCI(x=2.75,y=round((0.108*0.0004112)*100,2),
                li=round((0.108*(0.0004112-(1.96*0.02184)))*100,2),
                ui=round((0.108*(0.0004112+(1.96*0.02184)))*100,2),add=T,pch=16,col="orange")
text("0.022",x=2.75,y=-2.5,cex=0.8)
# Robben
plotrix::plotCI(x=3.25,y=round((0.108*-0.09776)*100,2),
                li=round((0.108*(-0.09776-(1.96*0.02298)))*100,2),
                ui=round((0.108*(-0.09776+(1.96*0.02298)))*100,2),add=T,pch=16,col='purple')
text("0.023",x=3.25,y=-2.5,cex=0.8)

# Dissagregated data, closure effect from Year/Island model
# Dassen
plotrix::plotCI(x=3.75,y=round((0.108*-0.01886)*100,2),
                li=round((0.108*(-0.01886-(1.96*0.03908)))*100,2),
                ui=round((0.108*(-0.01886+(1.96*0.03908)))*100,2),add=T,pch=16,col="orange")
text("0.039",x=3.75,y=-2.5,cex=0.8)

# Robben
plotrix::plotCI(x=4.25,y=round((0.108*-0.08483)*100,2),
                li=round((0.108*(-0.08483-(1.96*0.03911)))*100,2),
                ui=round((0.108*(-0.08483+(1.96*0.03911)))*100,2),add=T,pch=16,col='purple')
text("0.039",x=4.25,y=-2.5,cex=0.8)

# Raw (unmodelled) mean differences
# Dassen
points(x=4.75,y=round((0.108*-0.0194378)*100,2),pch=3,col="orange",lwd=2)
# Robben
points(x=5.25,y=round((0.108*-0.1140992)*100,2),pch=3,col="purple",lwd=2)
axis(1,at=1:5,labels=c("1.RGB,A,Y","2.RGB,D,YI","3.D,YM","4.D,YI","Raw means"),cex.axis=1)

#EOF
