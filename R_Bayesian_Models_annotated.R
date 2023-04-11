# *************************************************************************************************
# TITLE: 			  R CODE FOR RUNNING BAYESIAN MODELS
# LAST REVISED:	08Feb2023
# BACKGROUND: 	Data were obtained from a randomized, prospective clinical trial investigating the performance of 
#					    	two analgesics at an ambulatory surgical center. The purpose of the original trial was to determine 
#					    	if treatment with a single dose (30mcg) of sublingual sufentanil would lead to a decreased length of 
#					     	stay in the postanesthesia care unit (PACU) when compared to the standard of care, IV fentanyl, for 
#				    		those undergoing general anesthesia. The sample includes 66 patients, 33 per drug treatment group.
# REFERENCE:    https://doi.org/10.1155/2022/5237877 ; Bürkner PC (2017). "brms: An R Package for Bayesian Multilevel
#               Models using Stan." Journal of Statistical Software, 80(1), 1--28. doi:10.18637/jss.v080.i01.
# *************************************************************************************************

# CALL LIBRARIES
library(dplyr) #data management
library(here) #directory for where data is stored
library(tidyverse) #data management
library(Rcpp) #required for other packages
library(brms) #Bayesian modeling capabilities
library(bayestestR) #particular Bayesian tests

# READ IN DATA OF DRUG TRIAL FROM LOCAL LOCATION
trial <- here("DataProcessed", "drugtrial.csv") %>%
  read.csv(.)


###########################################################################################################################
###########################################################################################################################
# Linear Regression Models
###########################################################################################################################
###########################################################################################################################


#########################################
# Frequentist simple linear regression
#########################################

# Syntax: <name of model object> <- glm(<outcome variable> ~ <predictor variable>, data = <datasetname>, family=<distribution corresponding to model type>) 
lin_reg <- glm(in_phase_1_to_out_of_phase_2 ~ groupn, 
               data=trial, 
               family='gaussian')

# Syntax: summary(<model object>) - function to show model parameter estimates/results
summary(lin_reg)

# Syntax: confint() - print confidence intervals in console
confint(lin_reg)

# Syntax: plot() - print diagnostic plots in console
plot(lin_reg)



#########################################
# Bayesian simple linear regression
#########################################

# Set initial starting values for chains by creating a list, will be used for all simple linear regressions
# Syntax: list(<model parameter> = <starting value>); be sure to list all parameters
inits <- list(
  Intercept = 0,
  sigma     = 1,
  beta      = 0 )

# Syntax: <new_list> <- list(<initial values list name>) - Create list of all initial values
list_of_inits <- list(inits, inits, inits)


# Syntax: using brm function for Bayesian modeling
#         <name of model object> <- brm(<outcome variable> ~ <predictor variable>, 
#                                   data = <datasetname>, 
#                                   family=<distribution corresponding to model type>,
#                                   prior = c(set_prior("<distribution(mean,SD)", class = "<name>")),
#                                   seed = <value - for reproducibility>,
#                                   init = <name of initial values list>,
#                                   warmup = <sets the # of burn-in iterations (those that will be 'thrown out')>,
#                                   iter = <# of total iterations for each chain including burn-in>
#                                   chains = <# of chains>,
#                                   cores = <#> to use for executing chains in parallel - for processing)
# Note:                             inv_gamma(0.01 0.01) commands that the dispersion parameter is set with shape and scale at 0.01 



########## 
## Pseudo Vague Prior: Variance = 1000 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(1000))
##########

fit_lin_1 <-brm(in_phase_1_to_out_of_phase_2 ~ groupn,
                data=trial,
                family='gaussian',
                prior = c(set_prior("normal(0,31.62278)", class = "b"),
                          set_prior("normal(0,31.62278)", class ="Intercept"),
                          set_prior("inv_gamma(0.01,0.01)", class="sigma")),
                seed= 123,
                init=list_of_inits,
                warmup = 1000, iter = 10000, chains = 2, cores=4,
                save_pars = save_pars(all = TRUE))

# Summarize parameters
summary(fit_lin_1)

# Obtain highest density posterior interval
bayestestR::hdi(fit_lin_1, ci=0.95) 

# Syntax: plot() - print Bayesian diagnostic plots to console, plots in one figure
plot(fit_lin_1)

# Request plots individually 
mcmc_plot(fit_lin_1, type="hist") #histogram
mcmc_plot(fit_lin_1, type="trace") #trace plot
mcmc_plot(fit_lin_1, type="acf") #autocorrelation plot

# Syntax: prior_summary() - print priors used in console
prior_summary(fit_lin_1)

# Extract posterior chains
post_samp <- as_draws(fit_lin_1)

# Combine and extract drug group posterior estimates (can add more list items if more than 2 chains)
xpost <- c(post_samp[[1]]$b_groupn, post_samp[[2]]$b_groupn) 

# Calculate the posterior probability that our group predictor is less than 0
mean(xpost < 0) 



##########
## Vague Prior: Variance = 10,000 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10,000))
##########
fit_lin_2 <- brm(in_phase_1_to_out_of_phase_2 ~ groupn, 
                 data=trial, 
                 family='gaussian', 
                 prior = c(set_prior("normal(0,100)", class = "b"),
                           set_prior("normal(0,100)", class = "Intercept"),
                           set_prior("inv_gamma(0.01,0.01)", class="sigma")),
                 seed= 123,
                 init=list_of_inits,
                 warmup = 1000, iter = 10000, chains = 2, cores=4)

summary(fit_lin_2)
bayestestR::hdi(fit_lin_2, ci=0.95) 
plot(fit_lin_2)
mcmc_plot(fit_lin_2, type="hist") 
mcmc_plot(fit_lin_2, type="trace") 
mcmc_plot(fit_lin_2, type="acf") 
prior_summary(fit_lin_2)


# OPTION 1 for calculating posterior probabilities:
# Extract posterior chains
post_samp2 <- as_draws(fit_lin_2)
xpost2 <- c(post_samp2[[1]]$b_groupn, post_samp2[[2]]$b_groupn) 

# Calculate the posterior probability that our group predictor is less than 0
mean(xpost2 < 0) 



# OPTION 2 for calculating posterior probabilities:
# Extract posterior chains
post_samp2 <- as_draws_df(fit_lin_2)

# Create an indicator for group < 0
post_samp2 <-post_samp2 %>%
  mutate(indicator=ifelse(b_groupn<0,1,0))

# Calculate the posterior probability
summary(post_samp2$indicator)



##########
## Optimistic Prior
## Group: Mean= -30, Var=100 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(100))
## Intercept: Mean=0, Variance = 10000 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10000))
##########
fit_lin_3 <- brm(in_phase_1_to_out_of_phase_2 ~ groupn, 
                 data=trial, 
                 family='gaussian', 
                 prior = c(set_prior("normal(-30,10)", class = "b", coef = "groupn"),
                           set_prior("normal(0,100)", class = "Intercept"),
                           set_prior("inv_gamma(0.01,0.01)", class="sigma")),
                 seed= 123,
                 init=list_of_inits,
                 warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(fit_lin_3)
bayestestR::hdi(fit_lin_3, ci=0.95) #get 95% HDP Credible Intervals
plot(fit_lin_3)
mcmc_plot(fit_lin_3, type="hist") 
mcmc_plot(fit_lin_3, type="trace") 
mcmc_plot(fit_lin_3, type="acf") 
prior_summary(fit_lin_3)

# Extract posterior chains
post_samp3 <- as_draws(fit_lin_3)

xpost3 <- c(post_samp3[[1]]$b_groupn, post_samp3[[2]]$b_groupn) 

# Calculate the posterior probability that our group predictor is less than 0
mean(xpost3 < 0) 


##########
## Skeptical Prior
## Group: Mean= 0, Var=100 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(100))
## Intercept: Mean=0, Variance = 10000 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10000))
##########
fit_lin_4 <- brm(in_phase_1_to_out_of_phase_2 ~ groupn, 
                 data=trial, 
                 family='gaussian', 
                 prior = c(set_prior("normal(0,10)", class = "b", coef = "groupn"),
                           set_prior("normal(0,100)", class = "Intercept"),
                           set_prior("inv_gamma(0.01,0.01)", class="sigma")),
                 seed= 123,
                 init=list_of_inits,
                 warmup = 1000, iter = 10000, chains = 2, cores=4)

summary(fit_lin_4)
bayestestR::hdi(fit_lin_4, ci=0.95) #get 95% HDP Credible Intervals
plot(fit_lin_4)
mcmc_plot(fit_lin_4, type="hist") 
mcmc_plot(fit_lin_4, type="trace") 
mcmc_plot(fit_lin_4, type="acf") 
prior_summary(fit_lin_4)

# Extract posterior chains
post_samp4 <- as_draws(fit_lin_4)

xpost4 <- c(post_samp4[[1]]$b_groupn, post_samp4[[2]]$b_groupn) 

# Calculate the posterior probability that our group predictor is less than 0
mean(xpost4 < 0) 




#########################################
# Frequentist multiple linear regression
#########################################
lin_reg_adj <- glm(in_phase_1_to_out_of_phase_2 ~ groupn + sex_n + proc_length_center, 
                   data=trial, 
                   family='gaussian')
summary(lin_reg_adj)
confint(lin_reg_adj)
plot(lin_reg_adj)







#########################################
# Bayesian multiple linear regression
#########################################


##########
## Pseudo Vague Prior: Variance = 1000 --> sd is defined in prior statement (so converted to correct scale SD=sqrt(1000))
##########
fit_mlr1 <-brm(in_phase_1_to_out_of_phase_2 ~ groupn + sex_n + proc_length_center,
               data=trial,
               family='gaussian',
               prior = c(set_prior("normal(0,31.6227766017)", class = "b"),
                         set_prior("normal(0,31.6227766017)", class="Intercept"),
                         set_prior("inv_gamma(0.01,0.01)", class="sigma")),
               seed= 123,
               init=list_of_inits,
               warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(fit_mlr1)
bayestestR::hdi(fit_mlr1, ci=0.95) #get 95% HDP Credible Intervals
plot(fit_mlr1)
mcmc_plot(fit_mlr1, type="hist") 
mcmc_plot(fit_mlr1, type="trace") 
mcmc_plot(fit_mlr1, type="acf") 
prior_summary(fit_mlr1) 

# Extract posterior chains
post_samp <- as_draws(fit_mlr1)

xpost <- c(post_samp[[1]]$b_groupn, post_samp[[2]]$b_groupn) 

# Calculate the posterior probability that our group predictor is less than 0
mean(xpost < 0) 



##########
## Vague Prior: Variance = 10000 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10000))
##########
fit_mlr2 <-brm(in_phase_1_to_out_of_phase_2 ~ groupn + sex_n + proc_length_center,
               data=trial,
               family='gaussian',
               prior = c(set_prior("normal(0,100)", class = "b"),
                         set_prior("normal(0,100)", class="Intercept"),
                         set_prior("inv_gamma(0.01,0.01)", class="sigma")),
               seed= 123,
               init=list_of_inits,
               warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(fit_mlr2)
bayestestR::hdi(fit_mlr2, ci=0.95) #get 95% HDP Credible Intervals
plot(fit_mlr2)
mcmc_plot(fit_mlr2, type="hist") 
mcmc_plot(fit_mlr2, type="trace") 
mcmc_plot(fit_mlr2, type="acf") 
prior_summary(fit_mlr2) 

# Extract posterior chains
post_samp2 <- as_draws(fit_mlr2)

xpost2 <- c(post_samp2[[1]]$b_groupn, post_samp2[[2]]$b_groupn) 

# Calculate the posterior probability that our group predictor is less than 0
mean(xpost2 < 0) 



##########
## Optimistic Prior
## Group: Mean=-30, Var=100 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(100))
## Intercept/Covariates: Mean=0, Variance = 10000 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10000))
##########

# Syntax: within set_prior(), use 'coef' to specify values for particular variable name, if desired
fit_mlr3 <-brm(in_phase_1_to_out_of_phase_2 ~ groupn + sex_n + proc_length_center,
               data=trial,
               family='gaussian',
               prior = c(set_prior("normal(-30,10)", class = "b", coef = "groupn"),
                         set_prior("normal(0,100)", class = "b"),
                         set_prior("normal(0,100)", class="Intercept"),
                         set_prior("inv_gamma(0.01,0.01)", class="sigma")),
               seed= 123,
               init=list_of_inits,
               warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(fit_mlr3)
bayestestR::hdi(fit_mlr3, ci=0.95) #get 95% HDP Credible Intervals
plot(fit_mlr3)
mcmc_plot(fit_mlr3, type="hist") 
mcmc_plot(fit_mlr3, type="trace") 
mcmc_plot(fit_mlr3, type="acf") 
prior_summary(fit_mlr3) 

# Extract posterior chains
post_samp3 <- as_draws(fit_mlr3)

xpost3 <- c(post_samp3[[1]]$b_groupn, post_samp3[[2]]$b_groupn) 

# Calculate the posterior probability that our group predictor is less than 0
mean(xpost3 < 0) 



##########
## Skeptical Priors
## Group: Mean= 0, Var=100 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(100))
## Intercept/Covariates: Mean=0, Variance = 10000 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10000))
##########

# Syntax: within set_prior(), use 'coef' to specify values for particular variable name, if desired
fit_mlr4 <-brm(in_phase_1_to_out_of_phase_2 ~ groupn + sex_n + proc_length_center,
               data=trial,
               family='gaussian',
               prior = c(set_prior("normal(0,10)", class = "b", coef = "groupn"),
                         set_prior("normal(0,100)", class = "b"),
                         set_prior("normal(0,100)", class="Intercept"),
                         set_prior("inv_gamma(0.01,0.01)", class="sigma")),
               seed= 123,
               init=list_of_inits,
               warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(fit_mlr4)
bayestestR::hdi(fit_mlr4, ci=0.95) #get 95% HDP Credible Intervals
plot(fit_mlr4)
mcmc_plot(fit_mlr4, type="hist") 
mcmc_plot(fit_mlr4, type="trace") 
mcmc_plot(fit_mlr4, type="acf") 
prior_summary(fit_mlr4) 

# Extract posterior chains
post_samp4 <- as_draws(fit_mlr4)

xpost4 <- c(post_samp4[[1]]$b_groupn, post_samp4[[2]]$b_groupn) 

# Calculate the posterior probability that our group predictor is less than 0
mean(xpost4 < 0) 





###########################################################################################################################
###########################################################################################################################
# Logistic Regression Models
###########################################################################################################################
###########################################################################################################################

###############################s##########
# Frequentist logistic regression
#########################################
log_reg <- glm(blockn ~ groupn, 
               data=trial, 
               family='binomial'(link="logit"))
plot(log_reg)

# Results on logOdds scale
summary(log_reg)
confint(log_reg)
plot(log_reg)

# Results on Odds Ratio scale (exponentiate the logOdds results)
exp(log_reg$coefficients)
exp(confint(log_reg))

#########################################
# Bayesian simple logistic regression
#########################################

# set initial starting values for chains
inits <- list(
  Intercept = 0,
  sigma     = 1,
  beta      = 0 )
list_of_inits <- list(inits, inits, inits)



##########
## Vague Prior: Variance = 10 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10))
##########
fit1 <-brm(blockn ~ groupn,
           data=trial,
           family='bernoulli',
           prior = c(set_prior("normal(0,3.16227766017)", class = "b"),
                     set_prior("normal(0,3.16227766017)", class="Intercept")),
           seed= 123,
           init=list_of_inits,
           warmup = 1000, iter = 10000, chains = 2, cores=4)
plot(fit1)
mcmc_plot(fit1, type="hist") 
mcmc_plot(fit1, type="trace") 
mcmc_plot(fit1, type="acf") 
prior_summary(fit1) 

# Extract posterior chains
post_samp <- as_draws_df(fit1)

# Create an indicator for exp(group) < 1, because null is 1 on OR scale
post_samp <-post_samp %>%
  mutate(exp_groupn=exp(b_groupn),
         exp_intercept=exp(b_Intercept),
         indicator=ifelse(exp_groupn<1, 1, 0))

# Calculate the posterior probability (i.e., the mean)
summary(post_samp$indicator)

# Get estimates and 95% HDP Credible Intervals for exponentiated (i.e., on Odds Ratio scale) intercept and group posterior values
summary(post_samp)
bayestestR::hdi(post_samp$exp_intercept, ci=0.95) 
bayestestR::hdi(post_samp$exp_groupn, ci=0.95) 
bayestestR::hdi(post_samp$b_Intercept, ci=0.95) #LogOdds scale
bayestestR::hdi(post_samp$b_groupn, ci=0.95) #LogOdds scale 



##########
## Pseudo Vague Prior: Variance = 1 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(1))
##########
fit2 <- brm(blockn ~ groupn, 
            data=trial, 
            family='bernoulli', 
            prior = c(set_prior("normal(0,1)", class = "b"),
                      set_prior("normal(0,1)", class = "Intercept")),
            seed= 123,
            init=list_of_inits,
            warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(fit2)
bayestestR::hdi(fit2, ci=0.95) 
plot(fit2)
mcmc_plot(fit2, type="hist") 
mcmc_plot(fit2, type="trace") 
mcmc_plot(fit2, type="acf") 
prior_summary(fit2)

# Extract posterior chains
post_samp2 <- as_draws_df(fit2)

# Create an indicator for exp(group) < 1, because null is 1 on OR scale
post_samp2 <-post_samp2 %>%
  mutate(exp_groupn=exp(b_groupn),
         exp_intercept=exp(b_Intercept),
         indicator=ifelse(exp_groupn<1, 1, 0))

# Calculate the posterior probability (i.e., the mean)
summary(post_samp2$indicator)

# Get estimates and 95% HDP Credible Intervals for exponentiated (i.e., on Odds Ratio scale) intercept and group posterior values
summary(post_samp2)
bayestestR::hdi(post_samp2$exp_intercept, ci=0.95) 
bayestestR::hdi(post_samp2$exp_groupn, ci=0.95) 
bayestestR::hdi(post_samp2$b_Intercept, ci=0.95) #LogOdds scale
bayestestR::hdi(post_samp2$b_groupn, ci=0.95) #LogOdds scale 



##########
## Optimistic Prior
## Group: Mean= -0.5, Variance = 2 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(2))
## Intercept: Mean=0, Variance = 10 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10))
##########
fit3 <- brm(blockn ~ groupn, 
            data=trial, 
            family='bernoulli', 
            prior = c(set_prior("normal(-0.5,1.41421356237)", class = "b", coef = "groupn"),
                      set_prior("normal(0,3.16227766017)", class = "Intercept")),
            seed= 123,
            init=list_of_inits,
            warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(fit3)
bayestestR::hdi(fit3, ci=0.95) 
plot(fit3)
mcmc_plot(fit3, type="hist") 
mcmc_plot(fit3, type="trace") 
mcmc_plot(fit3, type="acf") 
prior_summary(fit3)

# Extract posterior chains
post_samp3 <- as_draws_df(fit3)

# Create an indicator for exp(group) < 1, because null is 1 on OR scale
# Exponentiate the intercept and group posterior values
post_samp3 <-post_samp3 %>%
  mutate(exp_groupn=exp(b_groupn),
         exp_intercept=exp(b_Intercept),
         indicator=ifelse(exp_groupn<1, 1, 0))

# Calculate the posterior probability (i.e., the mean)
summary(post_samp3$indicator)

# Get estimates and 95% HDP Credible Intervals for exponentiated (i.e., on Odds Ratio scale) intercept and group posterior values
summary(post_samp3)
bayestestR::hdi(post_samp3$exp_intercept, ci=0.95) 
bayestestR::hdi(post_samp3$exp_groupn, ci=0.95) 
bayestestR::hdi(post_samp3$b_Intercept, ci=0.95) #LogOdds scale
bayestestR::hdi(post_samp3$b_groupn, ci=0.95) #LogOdds scale 



##########
## Skeptical Prior
## Group: Mean= 0, Variance = 2 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(2))
## Intercept: Mean=0, Variance = 10 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10))
##########
fit4 <- brm(blockn ~ groupn, 
            data=trial, 
            family='bernoulli', 
            prior = c(set_prior("normal(0,1.41421356237)", class = "b", coef = "groupn"),
                      set_prior("normal(0,3.16227766017)", class = "Intercept")),
            seed= 123,
            init=list_of_inits,
            warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(fit4)
bayestestR::hdi(fit4, ci=0.95) 
plot(fit4)
mcmc_plot(fit4, type="hist") 
mcmc_plot(fit4, type="trace") 
mcmc_plot(fit4, type="acf") 
prior_summary(fit4)

# Extract posterior chains
post_samp4 <- as_draws_df(fit4)

# Create an indicator for exp(group) < 1, because null is 1 on OR scale
# Exponentiate the intercept and group posterior values
post_samp4 <-post_samp4 %>%
  mutate(exp_groupn=exp(b_groupn),
         exp_intercept=exp(b_Intercept),
         indicator=ifelse(exp_groupn<1, 1, 0))

# Calculate the posterior probability (i.e., the mean)
summary(post_samp4$indicator)

# Get estimates and 95% HDP Credible Intervals for exponentiated (i.e., on Odds Ratio scale) intercept and group posterior values
summary(post_samp4)
bayestestR::hdi(post_samp4$exp_intercept, ci=0.95) 
bayestestR::hdi(post_samp4$exp_groupn, ci=0.95) 
bayestestR::hdi(post_samp4$b_Intercept, ci=0.95) #LogOdds scale
bayestestR::hdi(post_samp4$b_groupn, ci=0.95) #LogOdds scale 





#########################################
# Frequentist multiple logistic regression
#########################################
adj_log_reg <- glm(blockn ~ groupn + sex_n + proc_length_center, 
                   data=trial, 
                   family='binomial'(link="logit"))
              
# Results on logOdds scale
summary(adj_log_reg)
confint(adj_log_reg)
plot(adj_log_reg)

# Results on Odds Ratio scale
exp(adj_log_reg$coefficients)
exp(confint(adj_log_reg))



#########################################
# Bayesian multiple logistic regression
#########################################


##########
## Vague Prior: Variance = 10 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10))
##########
adj_fit1 <-brm(blockn ~ groupn + sex_n + proc_length_center,
               data=trial,
               family='bernoulli',
               prior = c(set_prior("normal(0,3.16227766017)", class = "b"),
                         set_prior("normal(0,3.16227766017)", class="Intercept")),
               seed= 123,
               init=list_of_inits,
               warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(adj_fit1, digits=10)
print(summary(adj_fit1),digits=10) #get posterior estimate for proc_length
bayestestR::hdi(adj_fit1, ci=0.95) 
plot(adj_fit1)
mcmc_plot(adj_fit1, type="hist") 
mcmc_plot(adj_fit1, type="trace") 
mcmc_plot(adj_fit1, type="acf") 
prior_summary(adj_fit1) 

# Extract posterior chains
post_samp <- as_draws_df(adj_fit1)

# Create an indicator for exp(group) < 1, because null is 1 on OR scale
# Exponentiate the intercept and group posterior values
post_samp <-post_samp %>%
  mutate(exp_groupn=exp(b_groupn),
         exp_intercept=exp(b_Intercept),
         exp_sex=exp(b_sex_n),
         exp_proc=exp(b_proc_length_center),
         indicator=ifelse(exp_groupn<1, 1, 0))

# Calculate the posterior probability (i.e., the mean)
summary(post_samp$indicator)

# Get estimates and 95% HDP Credible Intervals for exponentiated (i.e., on Odds Ratio scale) intercept and group posterior values
summary(post_samp)
bayestestR::hdi(post_samp$exp_intercept, ci=0.95) 
bayestestR::hdi(post_samp$exp_groupn, ci=0.95) 
bayestestR::hdi(post_samp$exp_sex, ci=0.95) 
bayestestR::hdi(post_samp$exp_proc, ci=0.95) 



##########
## Pseudo Vague Prior: Variance = 1 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(1))
##########
adj_fit2 <- brm(blockn ~ groupn + sex_n + proc_length_center, 
                data=trial, 
                family='bernoulli', 
                prior = c(set_prior("normal(0,1)", class = "b"),
                          set_prior("normal(0,1)", class = "Intercept")),
                seed= 123,
                init=list_of_inits,
                warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(adj_fit2)
bayestestR::hdi(adj_fit2, ci=0.95) 
plot(adj_fit2)
mcmc_plot(adj_fit2, type="hist") 
mcmc_plot(adj_fit2, type="trace") 
mcmc_plot(adj_fit2, type="acf") 
prior_summary(adj_fit2)

# Extract posterior chains
post_samp2 <- as_draws_df(adj_fit2)

# Create an indicator for exp(group) < 1, because null is 1 on OR scale
# Exponentiate the intercept and group posterior values
post_samp2 <-post_samp2 %>%
  mutate(exp_groupn=exp(b_groupn),
         exp_intercept=exp(b_Intercept),
         exp_sex=exp(b_sex_n),
         exp_proc=exp(b_proc_length_center),
         indicator=ifelse(exp_groupn<1, 1, 0))

# Calculate the posterior probability (i.e., the mean)
summary(post_samp2$indicator)

# Get estimates and 95% HDP Credible Intervals for exponentiated (i.e., on Odds Ratio scale) intercept and group posterior values
summary(post_samp2)
bayestestR::hdi(post_samp2$exp_intercept, ci=0.95) 
bayestestR::hdi(post_samp2$exp_groupn, ci=0.95) 
bayestestR::hdi(post_samp2$exp_sex, ci=0.95) 
bayestestR::hdi(post_samp2$exp_proc, ci=0.95) 


##########
## Optimistic Prior
## Group: Mean= -0.5, Variance = 2 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(2))
## Intercept: Mean=0, Variance = 10 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10))
##########
adj_fit3 <- brm(blockn ~ groupn + sex_n + proc_length_center, 
                data=trial, 
                family='bernoulli', 
                prior = c(set_prior("normal(0,3.16227766017)", class = "b"), 
                          set_prior("normal(-0.5,1.41421356237)", class = "b", coef = "groupn"),
                          set_prior("normal(0,3.16227766017)", class = "Intercept")),
                seed= 123,
                init=list_of_inits,
                warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(adj_fit3)
bayestestR::hdi(adj_fit3, ci=0.95) 
plot(adj_fit3)
mcmc_plot(adj_fit3, type="hist") 
mcmc_plot(adj_fit3, type="trace") 
mcmc_plot(adj_fit3, type="acf") 
prior_summary(adj_fit3)

# Extract posterior chains
post_samp3 <- as_draws_df(adj_fit3)

# Create an indicator for exp(group) < 1, because null is 1 on OR scale
# Exponentiate the intercept and group posterior values
post_samp3 <-post_samp3 %>%
  mutate(exp_groupn=exp(b_groupn),
         exp_intercept=exp(b_Intercept),
         exp_sex=exp(b_sex_n),
         exp_proc=exp(b_proc_length_center),
         indicator=ifelse(exp_groupn<1, 1, 0))

# Calculate the posterior probability (i.e., the mean)
summary(post_samp3$indicator)

# Get estimates and 95% HDP Credible Intervals for exponentiated (i.e., on Odds Ratio scale) intercept and group posterior values
summary(post_samp3)
bayestestR::hdi(post_samp3$exp_intercept, ci=0.95) 
bayestestR::hdi(post_samp3$exp_groupn, ci=0.95) 
bayestestR::hdi(post_samp3$exp_sex, ci=0.95) 
bayestestR::hdi(post_samp3$exp_proc, ci=0.95) 




##########
## Skeptical Prior
## Group: Mean= 0, Variance = 2 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(2))
## Intercept: Mean=0, Variance = 10 --> sd is defined in set_prior statement (so converted to correct scale SD=sqrt(10))
##########
adj_fit4 <- brm(blockn ~ groupn + sex_n + proc_length_center, 
                data=trial, 
                family='bernoulli', 
                prior = c(set_prior("normal(0,3.16227766017)", class = "b"), 
                          set_prior("normal(0,1.41421356237)", class = "b", coef = "groupn"),
                          set_prior("normal(0,3.16227766017)", class = "Intercept")),
                seed= 123,
                init=list_of_inits,
                warmup = 1000, iter = 10000, chains = 2, cores=4)
summary(adj_fit4)
print(summary(adj_fit4),digits=10) #get posterior estimate for proc_length
bayestestR::hdi(adj_fit4, ci=0.95) 
plot(adj_fit4)
mcmc_plot(adj_fit4, type="hist") 
mcmc_plot(adj_fit4, type="trace") 
mcmc_plot(adj_fit4, type="acf") 
prior_summary(adj_fit4)

# Extract posterior chains
post_samp4 <- as_draws_df(adj_fit4)

# Create an indicator for exp(group) < 1, because null is 1 on OR scale
# Exponentiate the intercept and group posterior values
post_samp4 <-post_samp4 %>%
  mutate(exp_groupn=exp(b_groupn),
         exp_intercept=exp(b_Intercept),
         exp_sex=exp(b_sex_n),
         exp_proc=exp(b_proc_length_center),
         indicator=ifelse(exp_groupn<1, 1, 0))

# Calculate the posterior probability (i.e., the mean)
summary(post_samp4$indicator)

# Get estimates and 95% HDP Credible Intervals for exponentiated (i.e., on Odds Ratio scale) intercept and group posterior values
summary(post_samp4)
bayestestR::hdi(post_samp4$exp_intercept, ci=0.95) 
bayestestR::hdi(post_samp4$exp_groupn, ci=0.95) 
bayestestR::hdi(post_samp4$exp_sex, ci=0.95) 
bayestestR::hdi(post_samp4$exp_proc, ci=0.95) 




# END OF PROGRAM 

print(sessionInfo())


#SESSION INFORMATION 
#R version 4.2.1 (2022-06-23 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19044)







