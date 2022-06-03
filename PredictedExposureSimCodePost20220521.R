
### This code generates results for Table 1 for manuscript:
### "Issues in Implementing Regression Calibration Analyses" by Boe et al. 2022
### Code written by Pamela Shaw and Eunyung Park. 
### Contact: Pamela.A.Shaw@kp.org

#### Read in necessary packages
library(openxlsx)  
library(mvtnorm)
library(boot)

### This program writes out results in Excel. Set PATH to local directory
PATH<-"~/Dropbox/My Documents")
setwd("PATH")

 # Simulation settings for manuscript
  N = 2500
  n = 250
  Nboot=1000
  ME= "SME"   ### measurement error model for main instrument
  cptype="perc"  ### CI type for bootstrap
  CORSET = 2 
  
## Need to Change for correctY and correct RC to be the 4 possible binary pairs and rerun code to get all results in manuscript
  correctY= 1
  correctRC= 1
  
##############Set up Simulation #####################

set.seed(2020124)
NSIM=1000

ans=ansn=rcans=NULL
ste=sten=rcste=Bste=NULL
cpI=cpXhat=cpBZ=cpBU=NULL

expit=function(beta,xvec) {   ### xvec should be (1,x,z, u)
	xb=sum(beta*xvec)
	return( exp(xb)/(1+ exp(xb)) )
}

#### Logistic Outome model 
b0= -1
bx=log(1.5)
bz=-log(1.3)
bu=log(1.75)
beta=c(b0,bx,bz,bu)

### covariance structure of X,Z,U

muX=0
muZ=0
muU=0

sigx=1
sigz=1
sigu=1


# Set up Correlation scenarios 
if (CORSET==1) {
  # High Correlation
   rhoxz=.5
   rhoxu=.5
   rhozu=.5
               }
 
if (CORSET==2) {
 
  # Medium Correlation
   rhoxz=.3
   rhoxu=.3
   rhozu=.3
                }

if (CORSET==3) {

   # Low Correlation
   rhoxz=.1
   rhoxu=.1
   rhozu=.1
               }

if (CORSET==4) {

   # No Correlation with U 
   rhoxz=.1
   rhoxu=0
   rhozu=0
              }

if (CORSET==5) {

   # No Correlation any pair
   rhoxz=0
   rhoxu=0
   rhozu=0
               } 

# SME Q=a0 + a1x+ a2*z + a3*u+sigq 
   a0=0.4
   a1=.5
   a2=.5
   a3=.2
   sigq=0.7
   
### Need subset with W = X + epsilon. variance of epsilon=vareps
   vareps = 0.7
   
   if(ME=="CME") vareps=1.5

####### Start simulation
for (i in 1:NSIM){ 
	
   XZU=rmvnorm(N,c(muX,muZ,muU),cbind(c(sigx^2,rhoxz*sigx*sigz,rhoxu*sigx*sigu), 
          c(rhoxz*sigx*sigz,sigz^2,rhozu*sigu*sigz),c(rhoxu*sigx*sigu,rhozu*sigz*sigu,sigu^2)))

# Intercept and XZU 
   Xmat=cbind(rep(1,N),XZU)

#### calculate probablilty that Y=1 given covariates
   pvec=NULL
   for(j in 1:N) pvec=c(pvec,expit(beta, Xmat[j,]))

   Y=ifelse(runif(N)<pvec,1,0)

   X=XZU[,1]
   Z=XZU[,2]
   U=XZU[,3]

# Error in W- additive error from classical ME
   e1e2=rmvnorm(N,c(0,0),cbind(c(vareps,0), c(0,vareps)))
   e1=e1e2[,1]
   e2=e1e2[,2]

### Classical Error model : biomarker
   W1=X+e1
   W2=X+e2

# Biomarker W1 and W2 only collected on n subjects for either validation or reliability study
   W1[(n+1):N]=NA
   W2[(n+1):N]=NA
  
# SME error prone covariate observed on everyone
   Q=a0+ a1*X + a2*Z + a3*U + rnorm(N,0,sigq)

   data=data.frame(Y, X, W1, W2, Q, Z, U)

   data$Subset = c(rep(TRUE, n), rep(FALSE, (N-n))) 

# Naive model : y on Q and z, U 
   if (ME == "SME") {  if(correctY==0) YmodelN=glm(Y~Q+Z, data=data, family="binomial") 
					if(correctY==1) YmodelN=glm(Y~Q+Z+U, data=data, family="binomial")
                    if(correctRC==0) rc=lm(W1~Q+Z, data=data, subset=Subset)
                    if(correctRC==1) rc=lm(W1~Q+Z+U, data=data, subset=Subset)


    
   } else if (ME == "CME") {  if(correctY==0) YmodelN=glm(Y~W1+Z, data=data, family="binomial") 
                              if(correctY==1) YmodelN=glm(Y~W1+Z+U, data=data, family="binomial") 
                              if(correctRC==0) rc=lm(W2 ~ W1 + Z, data=data, subset=Subset)
                              if(correctRC==1) rc=lm(W2 ~ W1 + Z + U, data=data, subset=Subset)
   }


   ansn=rbind(ansn, YmodelN$coef)  
   sten=rbind(sten, coef(summary(YmodelN))[, "Std. Error"]) 


####################setting Naive CI - (n=1000)
   Lownci=ansn-1.96*sten
   Uppnci=ansn+1.96*sten

   covern = function(x) ifelse(Lownci[,x] < beta[x] & Uppnci[,x] > beta[x], 1, 0)
   CoverNaive = sapply(c(1:3),covern)

# Consider regression calibration for Systematic measurement error model Q using W1 with Z and U
   rcans=rbind(rcans, rc$coef)
   rcste=rbind(rcste, coef(summary(rc))[,"Std. Error"])# Empirical SE
   Xhat=predict(rc,newdata=data) 
  
# Fit outcome model
   if(correctY==0) Ymodel=glm(Y~Xhat+Z, family="binomial") 
   if(correctY==1) Ymodel=glm(Y~Xhat+Z+U, family="binomial") 
   ans=rbind(ans, Ymodel$coef) 
   ste=rbind(ste, coef(summary(Ymodel))[, "Std. Error"])

#############Fit Linear RC SETTING CI FOR Ymodel2 wo boot
   LowRCci= ans-1.96*ste
   UppRCci= ans+1.96*ste

   coverrc = function(x) ifelse(LowRCci[,x] < beta[x] & UppRCci[,x] > beta[x], 1, 0)
   CoverRCM = sapply(c(1:3),coverrc)

# now define bootstrap function
    boot_funR = function(dat, inds) {

          bdataR = dat[inds,]
             
			if (ME == "SME") { 
				if(correctRC==0) brc_modR = lm(W1 ~ Q + Z, data = bdataR, subset = Subset)
				if(correctRC==1) brc_modR = lm(W1 ~ Q + Z + U, data = bdataR, subset = Subset)
  			} else if (ME == "CME") {  
  				if(correctRC==0) brc_modR = lm(W2 ~ W1 + Z, data = bdataR, subset = Subset)
  				if(correctRC==1) brc_modR = lm(W2 ~ W1 + Z + U, data = bdataR, subset = Subset)
 		    }
          
     	 bXhatR = predict(brc_modR, newdata = bdataR)
    
    	 if(correctY==0) bfinal_modelR = glm(Y ~ bXhatR + Z , data = bdataR, family="binomial")
		 if(correctY==1) bfinal_modelR = glm(Y ~ bXhatR + Z + U, data = bdataR, family="binomial")
    
    	 return(bfinal_modelR$coef)
  	 }
  
   rc_bootR = boot(data, boot_funR, strata = data$Subset, R = Nboot)
   Bste= rbind(Bste,apply(rc_bootR$t,2,sd))

################ SETTING CI FOR Ymodel WITH BOOTSTRAP STE
   LowBRCci = ans-1.96*Bste
   UppBRCci = ans+1.96*Bste

   coverB = function(x) ifelse(LowBRCci[,x] < beta[x] & UppBRCci[,x] > beta[x], 1, 0)
   CoverBoot = sapply(c(1:3),coverB)

   ############### Boot.ci
   ################### Intercept 
    Result0=boot.ci(rc_bootR, index=1, conf = 0.95, type = cptype)
    if (cptype=="bca") { BCoverI=ifelse(Result0$bca[1, 4] < b0 & Result0$bca[1, 5] > b0, 1, 0)
    
    } else if (cptype=="perc") { BCoverI=ifelse(Result0$percent[1, 4] < b0 & Result0$percent[1, 5] > b0, 1, 0)
    
    }
    cpI = rbind(cpI, BCoverI) 
    
    ################### Xhat 

    Result1=boot.ci(rc_bootR, index=2, conf = 0.95, type = cptype)

       if (cptype=="bca") { BCoverXhat=ifelse(Result1$bca[1, 4] < bx & Result1$bca[1, 5] > bx, 1, 0)
       
    } else if (cptype=="perc") { BCoverXhat=ifelse(Result1$percent[1, 4] < bx & Result1$percent[1, 5] > bx, 1, 0)
    
    }
    
    cpXhat = rbind(cpXhat, BCoverXhat)

    #################################### Z 
    Result2=boot.ci(rc_bootR, index=3, conf =0.95, type = cptype)
    
        if (cptype=="bca") { BCoverZ=ifelse(Result2$bca[1, 4] < bz & Result2$bca[1, 5] > bz, 1, 0)
        
    } else if (cptype=="perc") { BCoverZ=ifelse(Result2$percent[1, 4] < bz & Result2$percent[1, 5] > bz, 1, 0)
    
    }
    
    cpBZ = rbind(cpBZ, BCoverZ)
    

} ##### End of Simulation loop

######################################################################################################
############################## Summarize Results ###################################################
######################################################################################################
   
################# NAIVE MODEL Y=Q+Z+U ###################################################
nPerc_biasmean = ((apply(ansn[,1:3], 2, mean) - c(b0, bx, bz))/c(b0,bx,bz)) * 100
nPerc_biasmed = ((apply(ansn[,1:3], 2, median) - c(b0, bx, bz))/c(b0,bx,bz)) * 100
nMSEmean = (apply(ansn[,1:3], 2, mean) - c(b0, bx, bz))^2 + apply(sten[,1:3], 2, mean)^2
nMSEmedian = (apply(ansn[,1:3], 2, median) - c(b0, bx, bz))^2 + apply(ansn[,1:3], 2, mad)^2
BiasMSEn=cbind(nPerc_biasmean, nPerc_biasmed, nMSEmean, nMSEmedian) 

Naiven = cbind(c("N_Int","N_Q","N_Z"),apply(ansn[,1:3],2, mean), apply(ansn[,1:3],2, median),apply(sten[,1:3],2, mean), apply(sten[,1:3],2, median), 
                 apply(ansn[,1:3],2, sd), apply(ansn[,1:3],2, mad), apply(CoverNaive, 2, mean) ) 
Naive=cbind(Naiven, BiasMSEn)
colnames(Naive) = c("Coeff","Mean of coeff", "Median of coeff","ASE mean Std.Error of coeff", 
                       "Median Std.Error of coeff","ESE Sd of coeff", "Mad", "CoverProb","PerMeanBias", 
                       "PerMedianBias", "MSEmean", "MSEmedian") 

######## linear RC Model 
Perc_biasmean = ((apply(ans[,1:3], 2, mean) - c(b0, bx, bz))/c(b0,bx,bz)) * 100
Perc_biasmed = ((apply(ans[,1:3], 2, median) - c(b0, bx, bz))/c(b0,bx,bz)) * 100
MSEmean = (apply(ans[,1:3], 2, mean) - c(b0, bx, bz))^2 + apply(ste[,1:3], 2, mean)^2
MSEmedian = (apply(ans[,1:3], 2, median) - c(b0, bx, bz))^2 + apply(ans[,1:3], 2, mad)^2
BiasMSE=cbind(Perc_biasmean, Perc_biasmed, MSEmean, MSEmedian) 

FITRCr = cbind(c("RCN_Int","RCN_X", "RCN_Z"),apply(ans[,1:3],2, mean), apply(ans[,1:3],2, median), apply(ste[,1:3],2, mean), apply(ste[,1:3],2, median), apply(ans[,1:3],2, sd), apply(ans[,1:3],2, mad), apply(CoverRCM, 2, mean) ) 
FITRC = cbind (FITRCr, BiasMSE)
colnames(FITRC) = c("Coeff","Mean of coeff", "Median of coeff","ASE mean Std.Error of coeff", 
                       "Median Std.Error of coeff","ESE Sd of coeff", "Mad", "CoverProb","PerMeanBias", 
                       "PerMedianBias", "MSEmean", "MSEmedian") 

############################## all result Combined
Modelinfo=rbind(Naive, FITRC) 

rownames(Modelinfo)=c("Intercept1","Q","Zn", "Intercept2","Xhat","Z")

Modelinfo=data.frame(Modelinfo) 

############################################### BOOTSTRP ste#########################################################
boot = cbind(apply(Bste[,1:3], 2, mean), apply(Bste[,1:3], 2, median), apply(CoverBoot, 2, mean) )
boot =data.frame(boot)
######## perc boot CP
CP.perc= rbind(apply(cpI, 2, mean) , apply(cpXhat, 2, mean), apply(cpBZ, 2, mean))
CP.perc=data.frame(CP.perc) ######## BCA boot CP

#### Need to add some "empty cells" to paste all results together.
BootM=cbind(
  rbind(boot*0, boot),
  rbind(CP.perc*0, CP.perc)
)

BootM<-cbind(c("RCBN_Int","RCBN_X","RCB2_Z","RCB2_Int","RCB2_X","RCB2_Z"),BootM)
bootname<-paste(cptype,"BootCP",sep="")
colnames(BootM) = c("Coeff","BootASE mean Std.Error of coeff", "BootMedian Std.Error of coeff","SimpleBootCP", bootname) 
rownames(BootM)=c("Intercept1","Q","Zn", "Intercept2","Xhat","Z")

Model=cbind(Modelinfo, BootM)

######################################## Write out results #########################################################
#### This creates the filename, embedding key settigns for the chosen parameters
filename=paste("ME_Y",correctY,"RC",correctRC,"_test",N,"n",n,ME,"COR",CORSET,".xlsx",sep="")
print(filename)

results<-createWorkbook("results")
addWorksheet(results, "estimates")

### estimates
writeData(results,"estimates",Model)

########## prevalence
Yprev=mean(Y) 
addWorksheet(results, "prevalence")
writeData(results,"prevalence",Yprev)

################## simulation data 
Alliter=cbind(ansn, sten, Lownci, Uppnci,ans, ste, LowRCci,UppRCci, Bste, LowBRCci,UppBRCci )
addWorksheet(results, "AllIterations")
writeData(results,"AllIterations",Alliter)
saveWorkbook(results, file=filename, overwrite=TRUE) 
print(date())