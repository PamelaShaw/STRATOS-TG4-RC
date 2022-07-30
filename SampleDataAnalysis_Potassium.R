###################################################################################################################
########################################   Sample Data Analysis   #################################################
###################################### Calibrating Log-Potassium ##################################################
###################################################################################################################
# In this example we demonstrate:
# (1) The different results that can be obtained by leaving a calibration model covariate out of the outcome model.
# (2) Model-based standard errors that do not adjust for the uncertainty added by the calibration model step
#     can occasionally be too small.
###################################################################################################################

library(dplyr)
library(survey)
library(gee)
library(stringr)
library(sas7bdat)
library(data.table)
library(tidyr)

#Read in the data
load(file="U:/Baldoni_Sim_June2022/TG4Simulation/SampleData_1.RData")
samp <- as.data.table(samp)

#Regression calibration model 
samp.solnas <- samp[(solnas==1),]
lm.lpotassium <- lm(c_ln_k_bio1 ~ c_ln_k_avg + c_age + c_ln_bmi+ high_chol + usborn + female + bkg_pr + bkg_o + supp_use,data=samp.solnas)

#Create estimated (calibrated) exposure for potassium
samp[,c_ln_k_calib := predict(lm.lpotassium,newdata=samp)]
summary(predict(lm.lpotassium,newdata=samp))

#Create survey design
samp.design = svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=samp)

#Logistic regression of outcome on calibrated log potassium (Xhat) and all other covariates (Z)
OriginalModel_pota <- svyglm(hypertension  ~ c_ln_k_calib + c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o+ supp_use, 
                             design=samp.design,family=stats::quasibinomial())
sum_origmodel_pota<-summary(OriginalModel_pota)

#Logistic regression of outcome on calibrated log potassium (Xhat) and all other covariates (Z), EXCEPT supplement use
TestModel1_pota <- svyglm(hypertension  ~ c_ln_k_calib + c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o, 
                          design=samp.design,family=stats::quasibinomial())
sum_testmodel1_pota<-summary(TestModel1_pota)

#Standard error calculation using adjusted MI approach of Baldoni et al. (2021)

#Set seed for reproducibility
set.seed(19452)

#Create empty object for saving results 
bootans<-NULL

#Number of bootstrap re-samples
Boots=25

for(j in 1:Boots){cat("Boot: ",j,"\n")

  # Take bootstrap sample
  bID_solnas<-sample(samp[which(samp$solnas==1),]$subid,450,replace=T)
  SOLNAS_boot<-as.data.frame(setDT(samp, key = 'subid')[J(bID_solnas)])

  #Regression calibration model 
  lm.lpotassium.boot <- lm.lpotassium <- lm(c_ln_k_bio1 ~ c_ln_k_avg + c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o+ supp_use,
                                              data=SOLNAS_boot)
  
  #Create estimated (calibrated) exposure for potassium
  samp[,c_ln_k_calib_boot := predict(lm.lpotassium.boot,newdata=samp)]
  
  #Recreate survey design with new calibrated exposure in data
  samp.design.boot <- svydesign(id=~BGid, strata=~strat, weights=~bghhsub_s2, data=samp)

  #Fit models with c_ln_k_calib_boot
  ###################################
  #Logistic regression of outcome on calibrated log potassium (Xhat) and all other covariates (Z)
  OriginalModel_pota_boot <- svyglm(hypertension  ~ c_ln_k_calib_boot + c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o+ supp_use, 
                               design=samp.design.boot,family=stats::quasibinomial())
  sum_origmodel_pota_boot<-summary(OriginalModel_pota_boot)
  
  #Logistic regression of outcome on calibrated log potassium (Xhat) and all other covariates (Z), EXCEPT supplement use
  TestModel1_pota_boot <- svyglm(hypertension  ~ c_ln_k_calib_boot + c_age + c_ln_bmi + high_chol + usborn + female + bkg_pr + bkg_o, 
                            design=samp.design.boot,family=stats::quasibinomial())
  sum_testmodel1_pota_boot<-summary(TestModel1_pota_boot)
  
  #Save results
  bootans <-rbind(bootans,c(j,sum_origmodel_pota_boot$coefficients["c_ln_k_calib_boot",c(1:2)],
                            sum_testmodel1_pota_boot$coefficients["c_ln_k_calib_boot",c(1:2)]))

}

bootans<-data.frame(bootans)
colnames(bootans)<-c("N_Imp","potass_orig","potass_orig_SE","potass_test","potass_test_SE")

#Summarize results
#########################
#Save coefficients from correct model
beta_potass_original<-sum_origmodel_pota$coefficients["c_ln_k_calib",1]

#Save coefficients from incorrect model - without supp use
beta_potass_test<-sum_testmodel1_pota$coefficients["c_ln_k_calib",1]

#Compute MI-adjusted variance for both correct (original) and incorrect (test) models
var_potass_orig<-(median((bootans$potass_orig_SE)^2)+(mad(bootans$potass_orig)^2))
var_potass_test<-(median((bootans$potass_test_SE)^2)+(mad(bootans$potass_test)^2))

#Standard errors from original/test models for showing underestimated SEs
INCORRECT_SE_potass_orig<-sum_origmodel_pota$coefficients["c_ln_k_calib",2]
INCORRECT_SE_potass_test<-sum_testmodel1_pota$coefficients["c_ln_k_calib",2]

#Put all results together in a table
#Table presents results for a 20% increase in potassium
###############################################################
Potassium_Example_Table<-rbind(cbind(round(exp(beta_potass_original*log(1.2)),2),round(exp(beta_potass_original-1.96*(INCORRECT_SE_potass_orig))^log(1.2),2),
                                     round(exp(beta_potass_original+1.96*(INCORRECT_SE_potass_orig))^log(1.2),2)),
                               cbind(round(exp(beta_potass_original*log(1.2)),2),round(exp(beta_potass_original-1.96*sqrt(var_potass_orig))^log(1.2),2),
            round(exp(beta_potass_original+1.96*sqrt(var_potass_orig))^log(1.2),2)),
      cbind(round(exp(beta_potass_test*log(1.2)),2),round(exp(beta_potass_test-1.96*(INCORRECT_SE_potass_test))^log(1.2),2),
            round(exp(beta_potass_test+1.96*(INCORRECT_SE_potass_test))^log(1.2),2)),
      cbind(round(exp(beta_potass_test*log(1.2)),2),round(exp(beta_potass_test-1.96*sqrt(var_potass_test))^log(1.2),2),
            round(exp(beta_potass_test+1.96*sqrt(var_potass_test))^log(1.2),2)))

#Give row and column names to table
colnames(Potassium_Example_Table)<-c("OR for 20% Increase","Lower 95% CI","Upper 95% CI")
row.names(Potassium_Example_Table)<-c("Correct Model (Naive SE)","Correct Model (Adjusted SE)",
                                      "Incorrect Model (Naive SE)","Incorrect Model (Adjusted SE)")

#Print table results 
Potassium_Example_Table

