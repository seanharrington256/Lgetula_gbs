#### This script reads in multiple bestlhoods files and calculates the AIC and deltaAIC for each

################################################################################################################################
####  NOTE & WARNING: This will only work if you output 1 and exactly 1 parameter to this file for each estimated parameter
####     i.e., if you output a complex parameter AND the estimated parameters, your number of parameters will
####     be wrong. Or if you estimated a parameter but didn't output it, your number of parameters will also be 
####    wrong -- e.g., if you don't care about population sizes and hid these from the output for some reason
################################################################################################################################


model<-list.files(pattern="bestlhoods") # list out the bestlhoods files
output<-lapply(model, read.table, header=TRUE) # read in the actual values of parameters and likelihoods
lhoods<-sapply(output, function(x) x$MaxEstLhood) # get just the estimated likelihood
npars<-sapply(output, function(x) length(x)-2)  # get the number of parameters, again assuming the note at the top of this script
AIC<-(-2*lhoods)+(2*npars) # calculate the AIC score from number of parameters and the LnL
fit_table<-cbind(lhoods, npars, AIC)  # Combine this information together into a single object
rownames(fit_table)<-model
fit_table<-fit_table[order(fit_table[,"AIC"]),]  # sort by AIC
deltaAIC<-fit_table[,"AIC"]-fit_table[1,"AIC"]  ## calculate deltaAIC
rel_lhood<-exp(-0.5*deltaAIC)  # calculate the relative likelihood - step in getting to AIC
sum_rel_lhood<-sum(rel_lhood)  # sum of relative likelihoods, another step in getting to AICw
AICw<-rel_lhood/sum_rel_lhood  #calculate AICw
models_summ_fit<-cbind(round(fit_table, 1), round(deltaAIC,1), round(AICw, 3)) # combine everything relevant together
colnames(models_summ_fit)<-c("LnL", "N_pars", "AIC", "deltaAIC", "AICw") # rename the columns
##Write out a csv with the model fit data
write.csv(models_summ_fit, file="FSC_model_Fits.csv")


