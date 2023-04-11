**************************************************************************************************
* TITLE: 			STATA CODE FOR RUNNING BAYESIAN MODELS
* LAST REVISED:		08FEB2022
* BACKGROUND: 		Data were obtained from a randomized, prospective clinical trial investigating the performance of 
*					two analgesics at an ambulatory surgical center. The purpose of the original trial was to determine 
*					if treatment with a single dose (30mcg) of sublingual sufentanil would lead to a decreased length of 
*					stay in the postanesthesia care unit (PACU) when compared to the standard of care, IV fentanyl, for 
*					those undergoing general anesthesia. The sample includes 66 patients, 33 per drug treatment group.
* NOTE:				Within syntax explanations below, <> indicates that an argument or value must be changed when 
*					running the code with new data sources.
**************************************************************************************************

* Import dataset of drug trial from local location
* Syntax: import delimited "filepath.extension"
import delimited "C:\Users\salauren\Documents\Projects\Bayesian\DataProcessed\drugtrial.csv"


**************************************************************************************************
* LINEAR REGRESSION
**************************************************************************************************

* Run frequentist simple linear model of time to readiness of discharge
* Syntax:'regress' fits the linear regression, followed by <Y-variable-name> then <covariate(s)>
regress in_phase_1_to_out_of_phase_2 groupn


******************************************
******** SIMPLE LINEAR REGRESSION ******** 
******************************************

* Run Bayes simple linear model  of time to readiness of discharge
* Note: '_cons' is the intercept term name within STATA
* Syntax: 	in order
*			'bayes' calls Bayesian framework for modeling
*			'igammaprior(0.01 0.01)' commands that the dispersion parameter is set with shape and scale at 0.01 
*			'prior({outcome-variable:variable-name1} {outcome-variable:variable-name2} {etc}, distribution(x, y))' specify the prior values for the listed variables
*					Note: depending on distribution, x, y are either (mean, variance) [for normal distribution] or (lower, upper) [for uniform distribution]
*					Note: when the distribution and prior values are the same for multiple variables, the entire list of variables is then followed by the dist(x,y)
*			'nchains(x)'  specifies the number of chains to run in MCMC
*			'mcmcsize(x)' specifies the number of iterations for each chain
*			'burnin(x)'   sets the number of burn-in iterations (those that will be 'thrown out')
*			'rseed(x)'    sets the seed for reproducibility 
*			'init1({outcome-variable:} x {sigma2} y)' set the initial values (x and y) for the first chain [repeat the statement with init2 for the second chain]
*					Note: use of the colon specifies that all variables are collectively called in the statement
*			'nomleinitial' suppresses the maximum likelihood estimates (MLE) for the starting values of the model parameters
*			'initsummary'  specifies that the initial values used for simulation be printed in output
*			'hpd' prints the HPD credible intervals rather than default (equal-tailed credible intervals)
*			':'   the colon explains that the following arguments specify the model to be fitted
*			'regress' fits the linear regression, followed by Y-variable-name first and then covariate(s)


** Vague Priors ~ Normal(0, Variance=10,000)
bayes, igammaprior(0.01 0.01) prior({in_phase_1_to_out_of_phase_2:groupn} {in_phase_1_to_out_of_phase_2:_cons}, normal(0,10000)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) init1({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) init2({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) nomleinitial initsummary hpd : regress in_phase_1_to_out_of_phase_2 groupn

* Retrieve parameter estimates with 95% HPD credible intervals in separate statement
* Syntax:	alternative to 'hpd' in BAYES line is to use 'bayesstats summary, hpd clevel(95)' as its own line to call the HPD results
bayesstats summary, hpd clevel(95)

* Look at diagnostic plots (reminder: _cons is the intercept term)
* Syntax:	'bayesgraph' call for graphical commands
*			'diagnostics' calls for trace plots, autocorrelation plots, etc.
*			'_all' requests graphs for all model parameters or {variable-name} can call for individual model parameter graphs
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sigma2}

* Calculate the posterior probability of drug group <0 (Note: Upper() means 'less than hypothesis value' [i.e., 0 upper limit], while lower() is 'greater than')
bayestest interval {in_phase_1_to_out_of_phase_2:groupn}, upper(0)


** Pseudo Vague Priors ~ Normal(0, Variance=1000)
bayes, igammaprior(0.01 0.01) prior({in_phase_1_to_out_of_phase_2:groupn} {in_phase_1_to_out_of_phase_2:_cons}, normal(0,1000)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) init1({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) init2({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) nomleinitial initsummary hpd : regress in_phase_1_to_out_of_phase_2 groupn
bayesgraph diagnostics _all
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sigma2}
bayestest interval {in_phase_1_to_out_of_phase_2:groupn}, upper(0)


** Optimistic Priors: Group ~ Normal(-30, Variance=100), Intercept ~Normal(0,Variance=1000)
*	Syntax: Differing priors for separate model parameters can be specified; here, the primary explanatory variable (groupn) has informative values while the intercept has different prior values
bayes, igammaprior(0.01 0.01) prior({in_phase_1_to_out_of_phase_2:groupn}, normal(-30,100)) prior({in_phase_1_to_out_of_phase_2:_cons}, normal(0,10000)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) init1({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) init2({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) nomleinitial initsummary hpd : regress in_phase_1_to_out_of_phase_2 groupn
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sigma2}
bayestest interval {in_phase_1_to_out_of_phase_2:groupn}, upper(0)


** Skeptical Priors: Group ~ Normal(0, Variance=100)
*	Syntax: Differing priors for separate model parameters can be specified; here, the primary explanatory variable (groupn) has skeptical values while the intercept has different prior values
bayes, igammaprior(0.01 0.01) prior({in_phase_1_to_out_of_phase_2:groupn}, normal(0,100)) prior({in_phase_1_to_out_of_phase_2:_cons}, normal(0,10000)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) init1({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) init2({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) nomleinitial initsummary hpd : regress in_phase_1_to_out_of_phase_2 groupn
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sigma2}
bayestest interval {in_phase_1_to_out_of_phase_2:groupn}, upper(0)



******************************************
******* MULTIPLE LINEAR REGRESSION ******* 
******************************************

* Run frequentist multiple linear model of time to readiness of discharge
regress in_phase_1_to_out_of_phase_2 groupn sex_n proc_length_center


** Pseudo Vague Priors ~ Normal(0, Variance=1000)
bayes, igammaprior(0.01 0.01) prior({in_phase_1_to_out_of_phase_2:}, normal(0,1000)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) init1({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) init2({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) nomleinitial initsummary hpd : regress in_phase_1_to_out_of_phase_2 groupn sex_n proc_length_center 
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sex_n}
bayesgraph diagnostics {proc_length_center}
bayesgraph diagnostics {sigma2}
bayestest interval {in_phase_1_to_out_of_phase_2:groupn}, upper(0)


** Vague Priors ~ Normal(0, Variance=10,000)
bayes, igammaprior(0.01 0.01) prior({in_phase_1_to_out_of_phase_2:}, normal(0,10000)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) init1({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) init2({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) nomleinitial initsummary hpd : regress in_phase_1_to_out_of_phase_2 groupn sex_n proc_length_center
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sex_n}
bayesgraph diagnostics {proc_length_center}
bayesgraph diagnostics {sigma2}
bayestest interval {in_phase_1_to_out_of_phase_2:groupn}, upper(0)


** Optimistic Priors ~ Normal(-30, Variance=100)
bayes, prior({in_phase_1_to_out_of_phase_2:groupn}, normal(-30,100)) prior({in_phase_1_to_out_of_phase_2:_cons} {in_phase_1_to_out_of_phase_2:sex_n} {in_phase_1_to_out_of_phase_2:proc_length_center}, normal(0,10000)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) init1({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) init2({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) nomleinitial initsummary hpd : regress in_phase_1_to_out_of_phase_2 groupn sex_n proc_length_center
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sex_n}
bayesgraph diagnostics {proc_length_center}
bayesgraph diagnostics {sigma2}
bayestest interval {in_phase_1_to_out_of_phase_2:groupn}, upper(0)


** Skeptical Priors ~ Normal(0, Variance=100)
bayes, prior({in_phase_1_to_out_of_phase_2:groupn}, normal(0,100)) prior({in_phase_1_to_out_of_phase_2:_cons} {in_phase_1_to_out_of_phase_2:sex_n} {in_phase_1_to_out_of_phase_2:proc_length_center}, normal(0,10000)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) init1({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) init2({in_phase_1_to_out_of_phase_2:} 0  {sigma2} 1) nomleinitial initsummary hpd : regress in_phase_1_to_out_of_phase_2 groupn sex_n proc_length_center
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sex_n}
bayesgraph diagnostics {proc_length_center}
bayesgraph diagnostics {sigma2}
bayestest interval {in_phase_1_to_out_of_phase_2:groupn}, upper(0)


**************************************************************************************************
* LOGISTIC REGRESSION
**************************************************************************************************


*******************************************
******* SIMPLE LOGISTIC REGRESSION ******** 
*******************************************

* Frequentist simple logistic regression
* Syntax: 	'logistic' displays Odds Ratios (ORs) estimates
* 			'logit' displays logOdds parameter estimates
*			'logit/logistic' fits the logistic regression, followed by Y-variable-name first and then covariate(s)
logistic blockn groupn
logit blockn groupn 


* Run Bayes simple logistic model of time to readiness of discharge
* Note: '_cons' is the intercept term name within STATA
* Syntax: 	in order
*			'bayes' calls Bayesian framework for modeling
*			'igammaprior(0.01 0.01)' commands that the dispersion parameter is set with shape and scale at 0.01 
*			'prior({outcome-variable:variable-name1} {outcome-variable:variable-name2} {etc}, distribution(x, y))' specify the prior values for the listed variables
*					Note: depending on distribution, x, y are either (mean, variance) [for normal distribution] or (lower, upper) [for uniform distribution]
*					Note: when the distribution and prior values are the same for multiple variables, the entire list of variables is then followed by the dist(x,y)
*			'nchains(x)'  specifies the number of chains to run in MCMC
*			'mcmcsize(x)' specifies the number of iterations for each chain
*			'burnin(x)'   sets the number of burn-in iterations (those that will be 'thrown out')
*			'rseed(x)'    sets the seed for reproducibility 
*			'init1({outcome-variable:} x {sigma2} y)' set the initial values (x and y) for the first chain [repeat the statement with init2 for the second chain]
*					Note: use of the colon specifies that all variables are collectively called in the statement
*			'nomleinitial' suppresses the maximum likelihood estimates (MLE) for the starting values of the model parameters
*			'initsummary'  specifies that the initial values used for simulation be printed in output
*			'hpd' prints the HPD credible intervals rather than default (equal-tailed credible intervals)
*			'eform' 	specifies that the results are displayed in exponentiated (Odds ratio) form
*			':'   the colon explains that the following arguments specify the model to be fitted
*			'regress' fits the linear regression, followed by Y-variable-name first and then covariate(s)

** Pseudo Vague Priors ~ Normal(0, Variance=1)
bayes, prior({blockn:}, normal(0,1)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) initall({blockn:} 0) initsummary nomleinitial hpd eform: logit blockn groupn
bayesgraph diagnostics _all
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}

* this is testing the posterior probability on the logOdds scale so keep the hypothesis value at 0
bayestest interval {blockn:groupn}, upper(0)


** Vague Priors ~ Normal(0, Variance=10)
bayes, prior({blockn:}, normal(0,10)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) initall({blockn:} 0) initsummary nomleinitial hpd eform: logit blockn groupn
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayestest interval {blockn:groupn}, upper(0)


** Optimistic Priors ~ Normal(-0.5, Variance=2)
bayes, prior({blockn:groupn}, normal(-0.5,2)) prior({blockn:_cons}, normal(0,10)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) initall({blockn:} 0) initsummary nomleinitial hpd eform: logit blockn groupn
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayestest interval {blockn:groupn}, upper(0)


** Skeptical Priors ~ Normal(0, Variance=2)
bayes, prior({blockn:groupn}, normal(0,2)) prior({blockn:_cons}, normal(0,10)) nchains(2) mcmcsize(9000) burnin(1000) rseed(123) initall({blockn:} 0) initsummary nomleinitial hpd eform: logit blockn groupn
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayestest interval {blockn:groupn}, upper(0)



********************************************
******* MULTIPLE LOGISTIC REGRESSION ******* 
********************************************

* Frequentist simple logistic regression

* logistic block(Y)=groupn(B) - this gives you Odds Ratios estimates
logistic blockn groupn proc_length_center sex_n

* logit displays the logOdds parameter estimates
logit blockn groupn proc_length_center sex_n


** Pseudo Vague Priors ~ Normal(0, Variance=1)
bayes, prior({blockn:}, normal(0,1)) nchains(2) burnin(1000) rseed(123) initall({blockn:} 0) nomleinitial initsummary hpd eform: logit blockn groupn sex_n proc_length_center 
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sex_n}
bayesgraph diagnostics {proc_length_center}
bayestest interval {blockn:groupn}, upper(0)


** Vague Priors ~ Normal(0, Variance=10)
bayes, prior({blockn:}, normal(0,10)) nchains(2) burnin(1000) rseed(123) initall({blockn:} 0) nomleinitial initsummary hpd eform: logit blockn groupn sex_n proc_length_center
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sex_n}
bayesgraph diagnostics {proc_length_center}
bayestest interval {blockn:groupn}, upper(0)


** Optimistic Priors ~ Normal(-0.5, Variance=2)
bayes, prior({blockn:groupn}, normal(-0.5,2)) prior({blockn:_cons} {blockn:sex_n} {blockn:proc_length_center}, normal(0,10)) nchains(2) burnin(1000) rseed(123) initall({blockn:} 0) nomleinitial initsummary hpd eform: logit blockn groupn sex_n proc_length_center 
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sex_n}
bayesgraph diagnostics {proc_length_center}
bayestest interval {blockn:groupn}, upper(0)


** Skeptical Priors ~ Normal(0, Variance=2)
bayes, prior({blockn:groupn}, normal(0,2)) prior({blockn:_cons} {blockn:sex_n} {blockn:proc_length_center}, normal(0,10)) nchains(2) burnin(1000) rseed(123) initall({blockn:} 0) nomleinitial initsummary hpd eform: logit blockn groupn sex_n proc_length_center 
bayesgraph diagnostics {_cons}
bayesgraph diagnostics {groupn}
bayesgraph diagnostics {sex_n}
bayesgraph diagnostics {proc_length_center}
bayestest interval {blockn:groupn}, upper(0)



* END OF PROGRAM