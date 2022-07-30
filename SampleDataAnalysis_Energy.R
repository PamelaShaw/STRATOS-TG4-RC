###################################################################################################################
########################################   Sample Data Analysis   #################################################
######################################## Calibrating Log-Energy ###################################################
###################################################################################################################
# In this example we demonstrate:
# (1) How to obtain an outcome model regression parameter that is adjusted for mediation using Midthune's method of
#     estimating the association parameter when one of the calibration model covariates is a mediator of the 
#     relationship between the exposure and outcome
# (2) How to obtain standard errors that adjust for the uncertainty in the calibration model step 
#     when a calibration model covariate mediates the outcome-exposure relationship
###################################################################################################################

library(dplyr)
library(survey)
library(gee)
library(stringr)
library(sas7bdat)
library(data.table)
library(tidyr)
library(xtable)

options(survey.adjust.domain.lonely=TRUE)
options(survey.lonely.psu="adjust")

#Read in the data
load(file="U:/Baldoni_Sim_June2022/TG4Simulation/SampleData_1.RData")
samp <- as.data.table(samp)

#Regression Calibration Model 
samp.solnas <- samp[(solnas==1),]
lm.lenergy <- lm(c_ln_kcal_bio1 ~ c_ln_kcal_avg + c_age + c_ln_bmi+ high_chol + usborn + female + bkg_pr + bkg_o,data=samp.solnas)

#Create estimated (calibrated) exposure, xhat
samp[,c_ln_kcal_calib := predict(lm.lenergy,newdata=samp)]
summary(predict(lm.lenergy,newdata=samp))

#Create survey design
samp.design = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=samp)

####### Midthune's method
#1. Perform logistic regression of outcome on calibrated log energy, mediator, 
#  and any other variable in the calibration equation.
#  Note the estimated coefficients of calibrated log energy (be) and log BMI (bbmi).

model1 <- svyglm(high_risk  ~ c_ln_kcal_calib + c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o, 
                 design=samp.design,family=stats::quasibinomial())
sum_model1<-summary(model1)

be<-sum_model1$coefficients[c("c_ln_kcal_calib"),1]
se_be<-sum_model1$coefficients[c("c_ln_kcal_calib"),2]

bbmi<-sum_model1$coefficients[c("c_ln_bmi"),1]
se_bbmi<-sum_model1$coefficients[c("c_ln_bmi"),2]

#Model without BMI
model1a <- svyglm(high_risk  ~ c_ln_kcal_calib + c_age  + high_chol + usborn + female + bkg_pr + bkg_o, 
                  design=samp.design,family=stats::quasibinomial())
sum_model1a<-summary(model1a)

be_nobmi<-sum_model1a$coefficients[c("c_ln_kcal_calib"),1]
se_be_nobmi<-sum_model1a$coefficients[c("c_ln_kcal_calib"),2]

#Model without xhat
model1b <- svyglm(high_risk  ~ c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o, 
                  design=samp.design,family=stats::quasibinomial())
sum_model1b<-summary(model1b)

bbmi_noe<-sum_model1b$coefficients[c("c_ln_bmi"),1]
se_bbmi_noe<-sum_model1b$coefficients[c("c_ln_bmi"),2]

#2.	#  (a)	Perform linear regression of biomarker repeat measure on biomarker original measure and other covariates in the outcome regression 
model2a <- lm(c_ln_kcal_bio2   ~ c_ln_kcal_bio1 + c_age + high_chol + usborn + female + bkg_pr + bkg_o, 
              data=samp.solnas)

#Create estimated (calibrated) exposure xhat on subset 
samp.solnas[,c_ln_kcal_calib_subset := predict(model2a,newdata=samp.solnas)]

# (b) Perform linear regression of mediator on calibrated exposure from subset and Z. Record the estimated coefficient, ge.
model2b <- lm(c_ln_bmi ~ c_ln_kcal_calib_subset + c_age + high_chol + usborn + female + bkg_pr + bkg_o,
              data=samp.solnas)
sum_model2b<-summary(model2b)

ge<-sum_model2b$coefficients["c_ln_kcal_calib_subset",1]
se_ge<-sum_model2b$coefficients["c_ln_kcal_calib_subset",2]

#3.	An approximately unbiased estimate of the desired regression coefficient for calibrated log energy
#   i.e., in the logistic regression of outcome on calibrated log energy and Z is given by be + ge bbmi


RegCoeff_CalibratedEnergy<-be+ge*bbmi
RegCoeff_CalibratedEnergy

#20% increase in consumption
exp(RegCoeff_CalibratedEnergy*log(1.2))

#4.	Standard error estimation that adjusts for uncertainty added by estimated exposure
#######################################################################################

#Read in bootstrap replicate weights 
load(file=paste0('U:/Baldoni_Sim_June2022/TG4Simulation/Bootweights.RData'))

#Now start imputation outer loop
mi.list.boot = list()

#Set seed for reproducibility 
set.seed(27567)

#Number of MI bootstrap runs within each sample
M<-25

#Number of bootstrap weights for which we have replicates
NBOOT<-1000

#Create empty object to save results in
bootans<-NULL


for(j in 1:M){
  
  mi.list.boot[[j]] = list()
  
  #Take bootstrap sample
  bID1<-sample(samp[which(samp$solnas_rep==1),]$subid,95,replace=T)
  bID2<-sample(samp[which(samp$solnas==1),]$subid,(450-95),replace=T)
  bID_solnas<-c(bID1,bID2)
  
  samp.solnas.boot<-as.data.frame(setDT(samp.solnas, key = 'subid')[J(bID_solnas)])
  
  #bootstrap calibration equation
  lm.lenergy.boot <- lm(c_ln_kcal_bio1 ~ c_ln_kcal_avg + c_age + c_ln_bmi+ high_chol + usborn + female + bkg_pr + bkg_o,
                        data=samp.solnas.boot)

  #Create estimated (calibrated) exposure 
  samp[,c_ln_kcal_calib_boot := predict(lm.lenergy.boot,newdata=samp)]

  for(b in 1:NBOOT){  cat('Boot MI: ',j,". Bootweight: ",b,".\n")

    
    #select weights for bth bootstrap
    
    if(b==1){
      #2.	#  (a)	Perform linear regression of biomarker repeat measure on biomarker original measure and other covariates in the outcome regression 
      model2a_boot <- lm(c_ln_kcal_bio2   ~ c_ln_kcal_bio1 + c_age + high_chol + usborn + female + bkg_pr + bkg_o,
                         data=samp.solnas.boot)
      samp.solnas$c_ln_kcal_calib_boot<-predict(model2a_boot,newdata=samp.solnas)
      
      # (b) Perform linear regression of mediator on calibrated exposure from subset and Z. Record the estimated coefficient, ge.
      
      model2b_boot <- lm(c_ln_bmi ~ c_ln_kcal_calib_boot +  c_age + high_chol + usborn + female + bkg_pr + bkg_o,
                         data=samp.solnas)
      sum_model2b_boot<-summary(model2b_boot)
      
      ge_boot<-sum_model2b_boot$coefficients["c_ln_kcal_calib_boot",1]
      
      #Create HCHS/SOL Sample Design using PSU ID, strata, and weights
      samp.design.boot <- svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=samp)
      
    }
    else {
    #select weights for bth bootstrap
    
      bIDII_1<-sample(samp.solnas.boot[which(samp.solnas.boot$solnas_rep==1),]$subid,95,replace=T)
      bIDII_2<-sample(samp.solnas.boot[which(samp.solnas.boot$solnas==1),]$subid,(450-95),replace=T)
      bIDII_solnas<-c(bIDII_1,bIDII_2)
      
      #select those for inside bootstrap loop
      SOLNAS_ALL_insideboot<-as.data.frame(setDT(samp.solnas, key = 'subid')[J(bIDII_solnas)])
      
      #2.	#  (a)	Perform linear regression of biomarker repeat measure on biomarker original measure and other covariates in the outcome regression 
      model2a_boot <- lm(c_ln_kcal_bio2   ~ c_ln_kcal_bio1 + c_age + high_chol + usborn + female + bkg_pr + bkg_o,
                         data=SOLNAS_ALL_insideboot)
      SOLNAS_ALL_insideboot$c_ln_kcal_calib_subset_boot<-predict(model2a_boot,newdata=SOLNAS_ALL_insideboot)

      # (b) Perform linear regression of mediator on calibrated exposure from subset and Z. Record the estimated coefficient, ge.
      
      model2b_boot <- lm(c_ln_bmi ~ c_ln_kcal_calib_subset_boot + c_age + high_chol + usborn + female + bkg_pr + bkg_o,
                        data=SOLNAS_ALL_insideboot)
      sum_model2b_boot<-summary(model2b_boot)
      
      ge_boot<-sum_model2b_boot$coefficients["c_ln_kcal_calib_subset_boot",1]
      
      #Now begin bootstrap weight process
      bootweight_b<-bootrepweights[,b] #check 

      #add bth bootstrap weights to dataset - can't use merge function because of duplicate IDs
      #sort both datasets by ID
      samp_boot<-data.frame(samp,BOOTWEIGHT=bootweight_b)
      
      #Create HCHS/SOL Sample Design using PSU ID, strata, and weights
      samp.design.boot <- svydesign(id=~BGid, strata=~strat, weights=~BOOTWEIGHT, data=samp_boot)
      }
      #1. Perform logistic regression of outcome on calibrated log energy, mediator, 
      #  and any other variable in the calibration equation.
      #  Note the estimated coefficients of calibrated log energy (be) and log BMI (bbmi).
      
      model1_boot <- svyglm(high_risk  ~ c_ln_kcal_calib_boot + c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o,
                          design=samp.design.boot,family=stats::quasibinomial())
      sum_model1_boot<-summary(model1_boot)
    
      be_boot<-sum_model1_boot$coefficients[c("c_ln_kcal_calib_boot"),1]
      bbmi_boot<-sum_model1_boot$coefficients[c("c_ln_bmi"),1]
    
      #Model without BMI
      model1a_boot <- svyglm(high_risk  ~ c_ln_kcal_calib_boot + c_age  + high_chol + usborn + female + bkg_pr + bkg_o, 
                             design=samp.design.boot,family=stats::quasibinomial())
      sum_model1a_boot<-summary(model1a_boot)
      
      be_noBMI_boot<-sum_model1a_boot$coefficients["c_ln_kcal_calib_boot","Estimate"]
      
      #Model without Xhat
      model1b_boot <- svyglm(high_risk  ~ c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o,
                             design=samp.design.boot,family=stats::quasibinomial())
      sum_model1b_boot<-summary(model1b_boot)
      
      bbmi_noE_boot<-sum_model1b_boot$coefficients["c_ln_bmi","Estimate"]
      
      
      #3.	An approximately unbiased estimate of the desired regression coefficient for calibrated log energy
      #   i.e., in the logistic regression of outcome on calibrated log energy and Z is given by be + ge bbmi
      
      RegCoeff_CalibratedEnergy_boot<-be_boot+ge_boot*bbmi_boot
      RegCoeff_CalibratedEnergy_boot
    
    diff<-RegCoeff_CalibratedEnergy_boot-be_noBMI_boot
    
    mi.list.boot[[j]][[b]] <-data.frame(Imp=j,BootWeight=b,be_boot,bbmi_boot,
                                        ge_boot,RegCoeff_CalibratedEnergy_boot,be_noBMI_boot,bbmi_noE_boot,diff)
    
    bootans <-rbind(bootans,c(j,b,be_boot,bbmi_boot,ge_boot,RegCoeff_CalibratedEnergy_boot,be_noBMI_boot,bbmi_noE_boot,diff))
    
  }
  
  mi.list.boot[[j]] <- do.call(rbind,mi.list.boot[[j]])
  
  }

#Now take results from iterative procedure above and compute adjusted variance 

#Assign null object
BootVar<-NULL

#Write function that summarizes results for a particular coefficient
MIbootvar_summarize<-function(bootlist,varname){
  for(j in 1:length(bootlist)){
    dat_bootweight<-bootlist[[j]][-1,]
    dat_origweight<-bootlist[[j]][1,]  
    
    V<-var(dat_bootweight[,paste0(varname)])
    MADSQ<-(mad(dat_bootweight[,paste0(varname)]))^2
    MeanBeta<-dat_origweight[,paste0(varname)] 
    
    BootVar<-rbind(BootVar,c(V,MeanBeta,MADSQ))
  }
  
  BootVar<-as.data.frame(BootVar) 
  colnames(BootVar)<-c("bootvar","beta","bootmad")
  return(BootVar)
}

#Write function that actually computes the variance

MIboot_Var<-function(V,betas,mad,M){
  MeanBeta<-mean(betas) 
  MedBeta<-median(betas)
  Diff<-(betas-MeanBeta)^2
  
  Var<-(median(V)+(mad(betas)^2))
  
  return(sqrt(Var))
} 


be_boot_hr1<-MIbootvar_summarize(mi.list.boot,"be_boot")
benoBMI_boot_hr1<-MIbootvar_summarize(mi.list.boot,"be_noBMI_boot")
beTARGET_boot_hr1<-MIbootvar_summarize(mi.list.boot,"RegCoeff_CalibratedEnergy_boot")

#Table with results for a 20% increase 
EnergyBMI_Example_Table<-rbind(cbind(round(exp(be*log(1.2)),2),round(exp(be-1.96*MIboot_Var(be_boot_hr1$bootvar,be_boot_hr1$beta,be_boot_hr1$bootmad,M))^log(1.2),2),
                                               round(exp(be+1.96*MIboot_Var(be_boot_hr1$bootvar,be_boot_hr1$beta,be_boot_hr1$bootmad,M))^log(1.2),2)),
                                         cbind(round(exp(be_nobmi*log(1.2)),2),round(exp(be_nobmi-1.96*MIboot_Var(benoBMI_boot_hr1$bootvar,benoBMI_boot_hr1$beta,benoBMI_boot_hr1$bootmad,M))^log(1.2),2),
                                               round(exp(be_nobmi+1.96*MIboot_Var(benoBMI_boot_hr1$bootvar,benoBMI_boot_hr1$beta,benoBMI_boot_hr1$bootmad,M))^log(1.2),2)),
                                         cbind(round(exp(RegCoeff_CalibratedEnergy*log(1.2)),2),round(exp(RegCoeff_CalibratedEnergy-1.96*MIboot_Var(beTARGET_boot_hr1$bootvar,beTARGET_boot_hr1$beta,beTARGET_boot_hr1$bootmad,M))^log(1.2),2),
                                               round(exp(RegCoeff_CalibratedEnergy+1.96*MIboot_Var(beTARGET_boot_hr1$bootvar,beTARGET_boot_hr1$beta,beTARGET_boot_hr1$bootmad,M))^log(1.2),2)))


#Give row and column names to table
colnames(EnergyBMI_Example_Table)<-c("OR for 20% Increase","Lower 95% CI","Upper 95% CI")
row.names(EnergyBMI_Example_Table)<-c("Including BMI in outcome model","Omitting BMI from outcome model",
                                      "Mudthune's Method")

#Print table results 
EnergyBMI_Example_Table

