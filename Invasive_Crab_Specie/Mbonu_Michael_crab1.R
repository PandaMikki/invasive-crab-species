# Set working directory and read in the data file


Crab<-read.csv(file="Data_for_Dryad.csv", sep=",", header=T)
summary(Crab)

# Create a table with the number of samples clustered by Site & Sex
### 1C
with(Crab, table(Site, Sex))

table(Crab$Site, Crab$Sex)

###  1E

s<-subset(Crab, NumberLimbMissing!=0 & NumberLimbRegen!=0) # Excluding 0 values

cor.test(Crab$NumberLimbMissing,Crab$NumberLimbRegen) # 0.56 inclu. 0
cor.test(s$NumberLimbMissing,s$NumberLimbRegen) # 0.96 exclu. 0

par(mfrow=c(1,1)) #set plotting parameter to single row and single column
plot(s$NumberLimbMissing~s$NumberLimbRegen)+abline(0,1) # Plotting correlation

with(Crab,table(Site, NumberLimbMissing))
with(Crab,table(Site, NumberLimbRegen))

### 1F ----

summary(Crab$PerCrabDryMass) # Found a value of 0

### 1G ----

hist(Crab$Final.oxygen.concentration) # Non-symmetrical

plot(density(Crab$Final.oxygen.concentration))+
  abline(v=mean(Crab$Final.oxygen.concentration),col="red")
## QUESTION 2 ----

### A ----

#  They found interaction between body size & sex t=-2.02; p=0.045
#  and an interaction between temperature and site t=3.10; p=0.002.
#  Step() function used by authors for model selection. 

summary(lm(Metabolic.rate ~
             NumberLimbMissing+NumberLimbRegen+PerCrabDryMass+Sex+Site+Temp,data=Crab))

steptest<-step(lm(Metabolic.rate ~
                    (PerCrabDryMass+Sex+Site+Temp)^2,data=Crab))

summary(steptest)

crabmass.reg<- lm(PerCrabDryMass~Sex, data = Crab)
crabtemp.reg<- lm(Temp~Site, data = Crab)
summary(crabmass.reg)
summary(crabtemp.reg)

anova(crabmass.reg)
anova(crabtemp.reg)
## QUESTION 3 ----

### A ----

#  Sub-setting only the female crabs in New Jersey

NJ_f<-subset(Crab, Site=="NJ"& Sex=="Female")

#  As mentioned in section 3.1 they have looked at the metabolic rate with the
#  predictor variables dry body mass, missing limbs, regenerating limbs & the
#  temperature with all interactions included in the (linear) model.From this I
#  would assume that they looked at all multi-interactions as well. So using * 
#  instead of +. That would look like this;

model_njf<-lm(Metabolic.rate~PerCrabDryMass*NumberLimbMissing*NumberLimbRegen*Temp, data=NJ_f)

#  Then they used the step() function to select the best fitting model from all
#  the interactions that are possible in the above model. This function is a part
#  of the "base" package.

modelnjf<-step(model_njf) # Save the best model from the test

summary(modelnjf) 

#  This does not seem to align with what the authors have written, so maybe it 
#  would make more sense to look at only two way interactions.

c_model_njf<-lm(Metabolic.rate~(PerCrabDryMass+NumberLimbMissing+NumberLimbRegen+Temp)^2, data=NJ_f)

cmodelnjf<-step(c_model_njf) # Save the best model from the test

summary(cmodelnjf)

#  The coefficients are now similar to the ones in the paper. The 0.191 estimate
#  means that for every 1 gram increase in mass, the metabolic rate rises with
#  0.191 ml O2/h-1. We can plot this to look at how it is presented in the paper.

attach(NJ_f)

plot(PerCrabDryMass, Metabolic.rate, pch=NumberLimbRegen,
     ylab="Metabolic rate (ml O2-1 h-1)", xlab= "Dry mass (g)")

detach(NJ_f)

### B ----

#  Sub-setting only the male crabs in New Jersey

NJ_m<-subset(Crab, Site=="NJ" & Sex=="Male")


#  As in A we can look at the model with all the two way interactions.

c_model_njm<-lm(Metabolic.rate~(PerCrabDryMass+NumberLimbMissing+NumberLimbRegen+Temp)^2, data=NJ_m)

modelnjm<-step(c_model_njm) # Save the best model from the test
summary(modelnjm)

#  The 1.106 estimate means that for every 1 gram increase in mass, the metabolic
#  rate rises with 1.106 ml O2/h-1. This corresponds to what is written in the paper

### C ----

#  In question A I think I already did this. But I might not completely understand
#  the question that they ask. Maybe in A it was just supposed to use only the text
#  as an indication of what they have put in the model. 

library(MuMIn)

dredge(lm(Metabolic.rate~(PerCrabDryMass+NumberLimbMissing+NumberLimbRegen+Temp)^2, data=NJ_f, na.action=na.fail))

#  The lowest AIC score is not the model they used, but the model that was used 
#  is within the delta 2 from the top ranked model. The logical thing to do would 
#  be to select the simplest model, which would include only the dry mass & 
#  temperature. So why did it not get selected by the step function. I do not know
#  the answer yet.

## QUESTION 4 ----

### A ----

#  Sub-setting only the female crabs in New Hampshire

NH_f<-subset(Crab, Site=="NH" & Sex=="Female")

#  In the text they mention predictors dry mass, missing limbs, interaction missing
#  limbs & temperature and the separate effect (non-significant) of temperature for female crabs 
#  is mentioned.

modelnhf<-lm(Metabolic.rate~PerCrabDryMass+NumberLimbMissing+NumberLimbMissing:Temp+Temp, data=NH_f)

summary(modelnhf)

#  The coefficients can be interpreted as followed; for every 1 gram increase in
#  dry mass, the metabolic rate rises with 0.172 ml O2/h-1 and for every 1 degrees
#  increase in temperature the metabolic rate rises with 0.009 ml O2/h-1.

#  Sub-setting only the male crabs in New Hampshire

NH_m<-subset(Crab, Site=="NH" & Sex=="Male")

#  In the text the mention dry mass and dry mass:temperature interaction only for
#  male crabs.

modelnhm<-lm(Metabolic.rate~PerCrabDryMass+PerCrabDryMass:Temp, data=NH_m)

summary(modelnhm)

#  This does not show the same results so I assumed that they must have included 
#  the separate effect of temperature once again. Even though it is not stated.

modelnhm2<-lm(Metabolic.rate~PerCrabDryMass+PerCrabDryMass:Temp+Temp, data=NH_m)

summary(modelnhm2)

#  The coefficients can be interpreted as followed; for every 1 gram increase in
#  dry mass, the metabolic rate decreases with 0.502 ml O2/h-1 and for every 1 
#  degrees increase in temperature the metabolic rate decreases with 0.031 ml O2/h-1.

### B ----

#  We can either look at the paper and compare, but we can also look at the
#  coefficients produced by the models used in the paper. We already have those
#  for NJ and NH, but we can make them quite quick for CT as well.

#  Sub-setting only the female crabs in Connecticut.

CT_f<-subset(Crab, Site=="CT" & Sex=="Female")

#  They only included dry mass in the model for female crabs in CT.

modelctf<-lm(Metabolic.rate~PerCrabDryMass, data=CT_f)

summary(modelctf)

#  Sub-setting only the male crabs in Connecticut.

CT_m<-subset(Crab, Site=="CT" & Sex=="Male")

#  For the male crabs they included dry mass, and both regenerating and missing
#  limbs.

modelctm<-lm(Metabolic.rate~PerCrabDryMass+NumberLimbMissing+NumberLimbRegen, data=CT_m)

summary(modelctm)

#  Now we can put all these estimate for the different sites and sexes in a little
#  table and maybe plot it so it is easier to compare them. 

DryMassCoeff <- data.frame(Site=rep(c('NJ', 'NH', 'CT'), each=2),
                           Sex=rep(c('Female', 'Male'), times=3),
                           Estimate=c(0.191,1.106,0.172,-0.502,0.184,0.324))
library(knitr)

kable(DryMassCoeff)

ggplot(data=DryMassCoeff, mapping = aes(x=Site, y=Estimate, colour=Sex))+
  geom_point()

#  From this we can see that the female crabs have similar increase in metabolic rate
#  as a function of mass over the different locations. Male crabs show much higher
#  variability. 

### C ----




## QUESTION 5 ----

### A ----

ECT<-subset(Crab, Site=="NJ"|Site=="NH")

plot(Metabolic.rate, Temp,
     xlab="Metabolic rate (ml O2-1 h-1)", ylab= "Temperature", data=ECT)+
  abline(lm(Metabolic.rate~Temp, data=ECT))

#  As we can see, the data does not look linear at all. We can look at the
#  residual plots for this model.

par(mfrow=c(2,2))

plot(glm(Metabolic.rate~Temp, data=ECT))

#  We can approach this non-linearity by performing....
#  Is there an confounding effect of dry mass?

### B ----

## QUESTION 6 ----