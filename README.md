# STRATOS-TG4-RC
Code for examples from the STRATOS-TG4 guidance paper on regression calibration

Files include:

(1) SampleData_1.RData, a simulated dataset representing hypothetical data with similar features to the Hispanic Community Health Study/Study of Latinos (HCHS/SOL) cohort that will be used to illustrate our examples from the paper.

(2) Bootweights.RData, a set of replicate bootstrap weights that will be to account for the complex survey design in the variance estimation stage in  Sample_Data_Analysis_Energy.R, described in (4) below.

(3) Sample_Data_Analysis_Potassium.R, an R code file that performs a sample data analysis to assess the association between potassium intake and hypertension using the simulated data from (1) above. Note that the user will need to change the path in the line of code where the data is read in. This code illustrates the example from Section 3 of our manuscript. First, we use sub-study data to fit the calibration model equation. Then, we perform a logistic regression of hypertension on the estimated exposure from regression calibration and other confounders. Next, we fit another logistic regression outcome model, this time excluding a calibration model covariate (supplement use) from the model. Finally, we use the multiple imputation-based variance approach of Baldoni et al. (2021) to calculate standard error estimates of the regression parameter that account for the complex survey design and are adjusted for the uncertainty added by the calibration model step. In this example, we demonstrate two points: (a) different odds ratios may be obtained when a calibration model covariate is left out of the outcome model and (b) model-based standard errors that do not adjust for the uncertainty added by the calibration model step may occasionally be too small.

(4) Sample_Data_Analysis_Energy.R, an R code file that performs a sample data analysis to assess the association between energy intake and a metabolic syndrome outcome using the simulated data from (1) and replicate bootstrap weights from (2) above. Note that the user will need to change the path in the lines of code where the data and bootstrap weights are read in. This code illustrates the example from Section 5 of our manuscript that shows how to implement Midthune’s method for obtaining an estimate of the regression coefficient for energy that accounts for mediation, when BMI is a mediator of the energy-health outcome relationship. Standard errors for the regression parameters are estimated using the approach outlined in Part 7 of the Appendix.

