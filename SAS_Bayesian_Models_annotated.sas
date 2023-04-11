*************************************************************************************************
*	TITLE: 			SAS CODE FOR FITTING MODELS IN A BAYESIAN FRAMEWORK
*	LAST REVISED:	09Feb2022
*	BACKGROUND: 	Data were obtained from a randomized, prospective clinical trial investigating the performance of 
*					two analgesics at an ambulatory surgical center. The purpose of the original trial was to determine 
*					if treatment with a single dose (30mcg) of sublingual sufentanil would lead to a decreased length of 
*				 	stay in the postanesthesia care unit (PACU) when compared to the standard of care, IV fentanyl, for 
*	    			those undergoing general anesthesia. The sample includes 66 patients, 33 per drug treatment group.
*					Resource: https://support.sas.com/resources/papers/proceedings14/SAS400-2014.pdf
*	NOTE:			Within syntax explanations below, <> indicates that an argument or value must be changed when 
					running the code with new data sources.
*************************************************************************************************; run;



************************************************************************* ;

*	SET WORKING DIRECTORY AND PROGRAM OPTIONS *;

************************************************************************* ;
%let root   =  C:\Users\USERNAME\Documents\Bayesian;
libname data  		"&root\Data";
options fmtsearch =	(data.formats work) ; 
ods graphics on;

* Reading in csv/excel;
proc import datafile = "&root\Data\drugtrial.csv"
	out = data.trial
	dbms = csv 
		replace;
	getnames = yes;
	guessingrows = 50 ;
run;


************************************************************************* ;

*	FORMATS		*;

************************************************************************* ;
proc format library = data.formats;
	value groupn_f
		1 = "Sublingual Sufentanil"
		0 = "Fentanyl" 	;
	value sex_f
		1 = 'Female'
		0 = 'Male'	;
	value yesno 
		0 = 'No'
		1 = 'Yes';
run;


************************************************************************* ;

*	LINEAR REGRESSION MODELING	*;

************************************************************************* ; run;

*****************************************************************
*	UNADJUSTED MODELS;
*****************************************************************;

* PROC GENMOD used to run frequentist linear and logistic regressions;
* Syntax:	'class' statement lists the categorical variable(s) and specify reference group name
			'model' statement:  <Y-variable-name> = <covariate(s)> / dist = <distribution> link = <link function>
				Note: normal distribution and identity link for linear models, binomial probability distribution and logit for logistic regression;

* Frequentist simple linear regression; 
title "SLR using PROC GENMOD (normal, identity link)";
proc genmod data=data.trial ;
	class groupn(ref="Fentanyl") / param=ref;
	model in_phase_1_to_out_of_phase_2 = groupn / 	dist=normal
													link=identity	;
	format groupn groupn_f.;
run;
title;

*	PROC MCMC:
*	Syntax: 'proc mcmc' statement tells SAS to use the MCMC procedure

			--- within proc mcmc: 
				'diag = <options>' - specify which convergence diagnostic tests to include in output
				'dic' - computes the Deviance Information Criterion (DIC)
				'nbi = <#>'  - specify number of burn-in iterations
				'nmc = <#>'  - specify number of iterations after burn-in
				'seed = <#>' - set the seed for reproducibility
				'plots(smooth)= <option>' - details which plots to display in output
				'statistics = all' - display posterior statistics
				'outpost = dataset' - generates dataset with the iteration numbers and posterior samples of all model parameters

			--- after proc mcmc statement:
				'parms <(var list)> <#>' - specifies the names and starting values of the parameters in the model
					Note: 	1. the variables listed are new names for the variables as model parameters (e.g., "beta_varname", "b_varname"', etc.)
							2. put variables with the same initial value in parentheses to apply value in sweeping statement.
						  	3. for variables with different starting values, include individual 'parms' statements without the parentheses.
							4. every parameter in the 'parms' statement must have a corresponding prior distribution in the 'prior' statement.
				'prior parameter-list ~ <distribution(x, y)>' - specifies the prior distribution for model parameters, can have individual statements for model parameters
					with different specifications on distribution and values for mean/variance or low/upper bounds.
				'mu = model-parameter-list' - defines the model for analysis, 
					Note: intercept term is model parameter name followed by model-param-name*variable-name (similar to regression equation notation)
				'model <outcome variable> ~ <distribution(mu, var=<dispersion-paramter-name>)' - specify the dependent variable and its distribution specifications ;		


* 	Pseudo-Vague Normal Prior N(0, 1000);
title  "Bayesian SLR using PROC MCMC";
title2 "Priors: N(0, 1000)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";
proc mcmc data=data.trial 
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123
	dic
	plots(smooth)=all
	statistics=all
	outpost=post1;
	parms (beta0 b_groupn) 0; /*starting values = 0*/
	parms sigma2 1; /*starting value = 1*/
	prior beta0 b_groupn ~ normal(mean=0, var=1000);
	prior sigma2 ~ igamma(shape = .01, scale = .01);
	mu=beta0 + b_groupn*groupn ;
	model in_phase_1_to_out_of_phase_2 ~ normal(mu, var=sigma2);
	format groupn groupn_f.;
run;

data Prob; *create an indicator for the values of treatment group to calculate posterior probability based on hypothesis;
	set post1; 
	Indicator = (B_groupn < 0); 
	label Indicator ='druggrp < 0';
run;

proc means data=prob(keep=Indicator); *calculate the posterior probability (i.e., the mean of the indicator);
run;
title;


*****************************************************************;
*	EQUIVALENT CODE FOR GENMOD + BAYES OPTIONS;
*****************************************************************;
*	PROC GENMOD wtih BAYES statement:
*	Syntax: 'proc genmod' statement tells SAS to use the genmod procedure and follows the same syntax as frequentist framework
			'BAYES' statement: included in order specify fitting the model within Bayesian framework

			--- Within BAYES statement: 
				'seed = <#>' - set the seed for reproducibility
				'nbi = <#>'  - specify number of burn-in iterations
				'nmc = <#>'  - specify number of iterations after burn-in
				'diag = <options>' - specify which convergence diagnostic tests to include in output
				'coeffprior=<distribution>(input=<dataset>)' - call the dataset where priors are specified for model parameters, in this case dataset = 'Prior'
				'plots= <option>' - details which plots to display in output
				'dispersionprior=<distribution>(shape = <a>, scale = <b>)' - set the shape and scale of the distribution for the variance parameter
				'initial=<dataset>'	- set the initial values for the chain by calling the dataset previously specified, in this case dataset = 'initial_values'
				'outpost = dataset' - generates dataset with the iteration numbers and posterior samples of all model parameters

			--- Set initial chain values and priors:
				1). initial values: must create a dataset with the starting chain values for each parameter in the model
						columns are 'chain' and list of the model parameters
						row values are the desired starting values for the chain
				2). priors: must create a dataset with the priors specified for each parameter
						columns are _TYPE_ (i.e., mean and variance), and list of the model parameters
						row values are the desired prior values for respective mean and variance for each parameter


* Vague Normal prior (0, 1000) using PROC GENMOD + BAYES; 
title "Bayesian SLR using PROC GENMOD + BAYES (normal, link-identity)";
title2 "Pseudo Vague Priors: N(0, 1000)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";

data initial_values; /*set initial values of the Markov chain in dataset to call in BAYES statement*/
	input chain $ intercept groupnSublingual_Sufentanil Dispersion scale;
datalines;
chain1 0 0 1 1
;
run;

data Prior; /*set priors in dataset to call in BAYES statement*/
	input _TYPE_ $ intercept groupnSublingual_Sufentanil; 
datalines; 
Mean 0 0
Var 1000 1000
; 
run;

ods graphics on; *must specify ODS graphics in order to generate plots;
proc genmod data=data.trial;
	class 	groupn(ref="Fentanyl") / param=ref;
	model 	in_phase_1_to_out_of_phase_2 = groupn / dist=normal link=identity	;
	bayes 	seed=123 nbi=1000 nmc=9000 diag=all coeffprior=normal(input=Prior) plots=all
			dispersionprior=igamma(shape = .01, scale = .01) initial=initial_values outpost=posterior;
	format 	groupn groupn_f.;
run;

* end of PROC GENMOD + BAYES example;



* 	Vague Normal Prior N(0, 10,000); 
title  "Bayesian SLR using PROC MCMC";
title2 "Priors: N(0, 10,000)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";
proc mcmc data=data.trial 
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123
	dic
	plots(smooth)=all 
	statistics=all
	outpost=post2;
	parms (beta0 b_groupn) 0;
	parms sigma2 1;
	prior beta0 b_groupn ~ normal(mean=0, var=10000);
	prior sigma2 ~ igamma(shape = .01, scale = .01);
	mu=beta0 + b_groupn*groupn ;
	model in_phase_1_to_out_of_phase_2 ~ normal(mu, var=sigma2);
	format groupn groupn_f.;
run;

data Prob; 
	set post2; 
	Indicator = (B_groupn < 0); 
	label Indicator ='druggrp < 0';
run;

proc means data=prob(keep=Indicator); 
run;
title;

* 	Optimistic Informative Prior N(-30, 100);
title "Bayesian SLR using PROC MCMC";
title2 "Prior: Normal(-30,100)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";
proc mcmc data=data.trial 
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123
	dic
	plots(smooth)=all 
	statistics=all
	outpost=post3;
	parms (beta0 b_groupn) 0;
	parms sigma2 1;
	prior beta0 ~ normal(mean=0, var=10000);
	prior b_groupn ~ normal(mean=-30, var=100);
	prior sigma2 ~ igamma(shape = .01, scale = .01);
	mu=beta0 + b_groupn*groupn;
	model in_phase_1_to_out_of_phase_2 ~ normal(mu, var=sigma2);	
	format groupn groupn_f. ;
run;

data Prob; 
	set post3; 
	Indicator = (B_groupn < 0); 
	label Indicator ='druggrp < 0';
run;

proc means data=prob(keep=Indicator); 
run;
title;

* 	Skeptical Prior N(0, 100);
title "Bayesian SLR using PROC MCMC";
title2 "Prior: Normal(0,100)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";
proc mcmc data=data.trial 
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123
	dic
	plots(smooth)=all 
	statistics=all
	outpost=post4;
	parms (beta0 b_groupn) 0;
	parms sigma2 1;
	prior beta0 ~ normal(mean=0, var=10000);
	prior b_groupn ~ normal(mean=0, var=100);
	prior sigma2 ~ igamma(shape = .01, scale = .01);
	mu=beta0 + b_groupn*groupn;
	model in_phase_1_to_out_of_phase_2 ~ normal(mu, var=sigma2);	
	format groupn groupn_f. ;
run;

data Prob; 
	set post4; 
	Indicator = (B_groupn < 0); 
	label Indicator ='druggrp < 0';
run;

proc means data=prob(keep=Indicator); 
run;
title;



*****************************************************************
*	ADJUSTED MODELS;
*****************************************************************;

* 	Frequentist Multiple Linear Regression;
title "MLR using PROC GENMOD (normal, identity link)";
proc genmod data=data.trial ;
	class groupn(ref="Fentanyl") sex_n(ref="Male") /param=ref;
	model in_phase_1_to_out_of_phase_2 = groupn sex_n proc_length_center 
							/ 	dist=normal
								link=identity	;
	format groupn groupn_f. sex_n sex_F.;
run;

* 	Pseudo-Vague Normal Prior N(0, 1000);
title  "Bayesian MLR using PROC MCMC";
title2 "Prior: Normal(0, 1000)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";
proc mcmc data=data.trial 
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123
	dic
	plots(smooth)=all 
	statistics=all
	outpost=post1;
	parms (beta0 b_groupn b_sex b_procleng) 0;
	parms sigma2 1;
	prior beta0 b_sex b_procleng ~ normal(mean=0, var=1000);
	prior b_groupn  ~ normal(mean=0, var=1000);
	prior sigma2 ~ igamma(shape =  .01, scale = .01);
	mu=beta0 + b_groupn*groupn + b_sex*sex_n + b_procleng*proc_length_center;
	model in_phase_1_to_out_of_phase_2 ~ normal(mu, var=sigma2);	
	format groupn groupn_f. sex_n sex_F.;
run;

data Prob; 
	set post1; 
	Indicator = (B_groupn < 0); 
	label Indicator ='druggrp < 0';
run;

proc means data=prob(keep=Indicator); 
run;
title;

* 	Vague Normal Prior N(0, 10000);
title "Bayesian MLR using PROC MCMC";
title2 "Prior: Normal(0, 10,000)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";
proc mcmc data=data.trial 
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123
	dic 
	plots(smooth)=all 
	statistics=all
	outpost=post2;
	parms (beta0 b_groupn b_sex  b_procleng ) 0;
	parms sigma2 1;
	prior sigma2 ~ igamma(shape = .01, scale = .01);
	prior beta0 b_groupn  b_sex b_procleng ~ normal(mean=0, var=10000);
	mu=beta0 + b_groupn*groupn + b_sex*sex_n  + b_procleng*proc_length_center;
	model in_phase_1_to_out_of_phase_2 ~ normal(mu, var=sigma2);	
	format groupn groupn_f. sex_n sex_F.;
run;

data Prob;
	set post2; 
	Indicator = (B_groupn < 0); 
	label Indicator ='druggrp < 0';
run;

proc means data=prob(keep=Indicator); 
run;
title;


*	Optimistic Informative Prior N(-30, 100);
title "Bayesian MLR using PROC MCMC";
title2 "Prior: Normal(-30,100)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";
proc mcmc data=data.trial 
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123
	dic 
	plots(smooth)=all 
	statistics=all
	outpost=post3;
	parms (beta0 b_groupn b_sex b_procleng) 0;
	parms sigma2 1;
	prior beta0 b_sex b_procleng ~ normal(mean=0, var=10000);
	prior b_groupn ~ normal(mean=-30, var=100);
	prior sigma2 ~ igamma(shape =  .01, scale = .01);
	mu=beta0 + b_groupn*groupn + b_sex*sex_n + b_procleng*proc_length_center;
	model in_phase_1_to_out_of_phase_2 ~ normal(mu, var=sigma2);	
	format groupn groupn_f. sex_n sex_F.;
run;

data Prob; 
	set post3; 
	Indicator = (B_groupn < 0); 
	label Indicator ='druggrp < 0';
run;

proc means data=prob(keep=Indicator);
run;
title;


* 	Skeptical Prior N(0, 100);
title "Bayesian MLR using PROC MCMC";
title2 "Prior: Normal(0,100)";
title3 "Defined Sigma^2: Igamma(shape = .01, scale = .01)";
proc mcmc data=data.trial 
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123
	dic 
	plots(smooth)=all 
	statistics=all
	outpost=post4;
	parms (beta0 b_groupn b_sex b_procleng) 0;
	parms sigma2 1;
	prior beta0 b_sex b_procleng ~ normal(mean=0, var=10000);
	prior b_groupn ~ normal(mean=0, var=100);
	prior sigma2 ~ igamma(shape =  .01, scale = .01);
	mu=beta0 + b_groupn*groupn + b_sex*sex_n + b_procleng*proc_length_center;
	model in_phase_1_to_out_of_phase_2 ~ normal(mu, var=sigma2);	
	format groupn groupn_f. sex_n sex_F.;
run;

data Prob; 
	set post4; 
	Indicator = (B_groupn < 0); 
	label Indicator ='druggrp < 0';
run;

proc means data=prob(keep=Indicator); 
run;
title;



************************************************************************* ;

*	LOGISTIC REGRESSION MODELING	*;

************************************************************************* ; run;

*****************************************************************
*	UNADJUSTED MODELS;
*****************************************************************;

* Syntax:	'class' statement: lists the categorical variable(s) and specify reference group name
			'model' statement:  <Y-variable-name> = <covariate(s)> / dist = <distribution> link = <link function>
				Note: binomial probability distribution and logit for logistic regression
			'estimate' statement: exponentiate log-odds results to Odds Ratio scale;

* Frequentist simple logistic regression; 
title1 "Logistic Unadjusted Model of Block";
proc genmod data=data.trial  ;
	class groupn(ref="Fentanyl")/ param=ref;
	model blockn(event="Yes") = groupn / 	dist = bin 
											link = logit ;
	format groupn groupn_f. blockn yesno.; 
	estimate 'Beta0' intercept 1 / exp;
	estimate 'B_group' groupn 1 / exp;
run;

*	Syntax: 'proc mcmc' statement tells SAS to use the MCMC procedure

			--- within proc mcmc: 
				'diag = <options>' - specify which convergence diagnostic tests to include in output
				'dic' - computes the Deviance Information Criterion (DIC)
				'nbi = <#>'  - specify number of burn-in iterations
				'nmc = <#>'  - specify number of iterations after burn-in
				'seed = <#>' - set the seed for reproducibility
				'plots(smooth)= <option>' - details which plots to display in output
				'statistics = all' - display posterior statistics
				'monitor = (<var list>)' - specify variables that are computed in subsequent beginnodata/endnodata block
				'outpost = dataset' - generates dataset with the iteration numbers and posterior samples of all model parameters

			--- after proc mcmc statement:
				'parms <(var list)> <#>' - specifies the names and starting values of the parameters in the model
					Note: 	1. the variables listed are new names for the variables as model parameters (e.g., "beta_varname", "b_varname"', etc.)
							2. put variables with the same initial value in parentheses to apply value in sweeping statement.
						  	3. for variables with different starting values, include individual 'parms' statements without the parentheses.
							4. every parameter in the 'parms' statement must have a corresponding prior distribution in the 'prior' statement.
				'prior parameter-list ~ <distribution(x, y)>' - specifies the prior distribution for model parameters, can have individual statements for model parameters
					with different specifications on distribution and values for mean/variance or low/upper bounds.
				'beginnodata / endnodata' - for calculations that relate to parameters only, with logistic regression you want to exponentiate values from logOdds scale to obtain Odds ratios
							within these statements: any computations that are identical to every observation, such as transformation of parameters, can be executed
				'p = logitic(<model-parameter-list>)' - defines the model for analysis, 
					Note: intercept term is model parameter name followed by model-param-name*variable-name (similar to regression equation notation)
				'model <outcome variable> ~ <distribution(p)>' - specify the dependent variable and its distribution probability ;		

* 	Pseudo-Vague Prior N(0,1);
title1 "Bayesian Unadjusted Model of Block";
title2 "Normal (0, 1) Priors";
title3 "Results on the logOdds scale"; *****results maintained on the logOdds Scale;
proc mcmc data=data.trial 
	plots=all
	nbi=1000 
	nmc=9000 
    seed=123 
	dic 
	plots(smooth)=all 
	statistics=all 
	outpost=post1;
	parms (beta0 B_groupn) 0;
	prior beta0 ~ normal(0,var=1);
	prior B_groupn ~ normal(0,var=1);
	p = logistic(beta0 + B_groupn*groupn);
	model blockn ~ binary(p);
	format groupn groupn_f. blockn yesno.;
run;

title1 "Bayesian Unadjusted Model of Block";
title2 "Normal (0, 1) Priors";
title3 "Results on the Odds Ratio scale";*****results transformed to Odds Ratio Scale;
proc mcmc data=data.trial 
	plots=all
	nbi=1000 
	nmc=9000 
    seed=123 
	dic 
	plots(smooth)=all 
	statistics=all 
	monitor=(or_beta0 or_groupn) 
	outpost=post1;
	parms (beta0 B_groupn) 0;
	prior beta0 ~ normal(0,var=1);
	prior B_groupn ~ normal(0,var=1);
		beginnodata;
		or_beta0 = exp(beta0);
		or_groupn = exp(B_groupn);
		endnodata;  	
	p = logistic(beta0 + B_groupn*groupn);
	model blockn ~ binary(p);
	format groupn groupn_f. blockn yesno.;
run;

data Prob; *create an indicator for the values of treatment group to calculate posterior probability based on hypothesis;
	set post1; 
	Indicator = (or_groupn < 1); *posterior probability that odds of group <1;
	label Indicator ='exp(druggrp)<1';
run;

proc means data=prob; *calculate the posterior probability (i.e., mean of indicator) and obtain OR estimates for parameters;
	var Indicator or_groupn or_beta0 ;
run;
title;

*****************************************************************;
*	EQUIVALENT CODE FOR GENMOD + BAYES OPTIONS;
*****************************************************************;

* Vague Normal prior (0, 1) using PROC GENMOD + BAYES; 
title "Bayesian SLR using PROC GENMOD + BAYES (normal, link-identity)";
title2 "Pseudo Vague Priors: N(0, 1)";
data initial_values; /* initial values of the Markov chains */
	input chain $ intercept groupnSublingual_Sufentanil scale;
datalines;
chain1 0 0 1
;
run;

data Prior; /*set priors in dataset to call in BAYES statement*/
	input _TYPE_ $ intercept groupnSublingual_Sufentanil; 
datalines; 
Mean 0 0
Var 1 1
; 
run;

ods graphics on; *must specify ODS graphics in order to generate plots;
proc genmod data=data.trial ;
	class groupn(ref="Fentanyl") / param=ref;
	model blockn(event="Yes") = groupn /   dist = bin 
							link = logit ;
	format groupn groupn_f. blockn yesno. ;
	bayes seed=123 nbi=1000 nmc=9000 DIAG=all coeffprior=normal(input=Prior) initial=initial_values outpost=genmod_post1;
run;

* end of PROC GENMOD + BAYES example;



* 	Vague Prior N(0,10);
title1 "Bayesian Unadjusted Model of Block";
title2 "Normal (0, 10) Priors";
title3 "Results on the Odds Ratio scale";
proc mcmc data=data.trial
	plots=all
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123 
	dic
	plots(smooth)=all 
	statistics=all
	monitor=(or_beta0 or_groupn) 
	outpost=post2;
   	parms (beta0 b_groupn) 0;
   	prior beta0 b_groupn ~ normal(mean=0, var=10);
		beginnodata;
		or_beta0 = exp(beta0);
		or_groupn = exp(B_groupn);
		endnodata;  	
   	p = logistic(beta0 + b_groupn*groupn) ;
   	model blockn ~ binary(p);
	format groupn groupn_f. blockn yesno.;
run;

data Prob; 
	set post2; 
	Indicator = (or_groupn < 1); *posterior probability that odds of group <1;
	label Indicator ='exp(druggrp)<1';
run;

proc means data=prob; *calculate the posterior probability (i.e., mean of indicator) and obtain OR estimates for parameters;
	var Indicator or_groupn or_beta0 ;
run;
title;


*	Optimistic Informative Prior N(-0.5, 2);
title1 "Bayesian Adjusted Model";
title2 "Normal (-0.5, 2) for 'Group' and N(0,1) Priors";
title3 "Results on the Odds Ratio scale";
proc mcmc data=data.trial 
	plots=all
	nbi=1000 
    nmc=9000
	seed=123 
	dic 	 
	plots(smooth)=all 
	statistics=all
	monitor=(or_beta0 or_groupn) 
	outpost=post3;
	parms (beta0 b_groupn) 0;
	prior beta0 ~ normal(mean=0, var=10);
	prior b_groupn ~ normal(mean=-0.5, var=2);
		beginnodata;
		or_beta0 = exp(beta0);
		or_groupn = exp(B_groupn);
		endnodata;  	
	p = logistic(beta0 + b_groupn*groupn);
   	model blockn ~ binary(p);
	format groupn groupn_f. blockn yesno.;
run;

data Prob; 
	set post3; 
	Indicator = (or_groupn < 1); *posterior probability that odds of group <1;
	label Indicator ='exp(druggrp)<1';
run;

proc means data=prob; *calculate the posterior probability (i.e., mean of indicator) and obtain OR estimates for parameters;
	var Indicator or_groupn or_beta0 ;
run;
title;

* 	Skeptical Prior N(0, 2);
title1 "Bayesian Adjusted Model";
title2 "Normal (0, 2) for 'Group' and N(0,1) Priors";
title3 "Results on the logOdds scale";
proc mcmc data=data.trial 
	plots=all
	nbi=1000 
    nmc=9000
	seed=123 
	dic 	 
	plots(smooth)=all 
	statistics=all
	outpost=post4;
	parms (beta0 b_groupn) 0;
	prior beta0 ~ normal(mean=0, var=10);
	prior b_groupn ~ normal(mean=0, var=2);
	p = logistic(beta0 + b_groupn*groupn);
   	model blockn ~ binary(p);
	format groupn groupn_f. blockn yesno.;
run;

title1 "Bayesian Adjusted Model";
title2 "Normal (0, 2) for 'Group' and N(0,1) Priors";
title3 "Results on the Odds Ratio scale";
proc mcmc data=data.trial 
	plots=all
	nbi=1000 
    nmc=9000
	seed=123 
	dic 	 
	plots(smooth)=all 
	statistics=all
	monitor=(or_beta0 or_groupn) 
	outpost=post4;
	parms (beta0 b_groupn) 0;
	prior beta0 ~ normal(mean=0, var=10);
	prior b_groupn ~ normal(mean=0, var=2);
		beginnodata;
		or_beta0 = exp(beta0);
		or_groupn = exp(B_groupn);
		endnodata;  	
	p = logistic(beta0 + b_groupn*groupn);
   	model blockn ~ binary(p);
	format groupn groupn_f. blockn yesno.;
run;

data Prob; 
	set post4; 
	Indicator = (or_groupn < 1); *posterior probability that odds of group <1;
	label Indicator ='exp(druggrp)<1';
run;

proc means data=prob; *calculate the posterior probability (i.e., mean of indicator) and obtain OR estimates for parameters;
	var Indicator or_groupn or_beta0 ;
run;
title;



*****************************************************************
*	ADJUSTED MODELS;
*****************************************************************; 

* 	Frequentist multiple logistic regression; 
title1 "Frequentist Adjusted Model of Block";
proc genmod data=data.trial  ;
	class groupn(ref="Fentanyl") sex_n(ref="Male")/ param =ref;
	model blockn(event="Yes") = groupn sex_n proc_length_center 
						  / dist = bin 
							link = logit ;
	format groupn groupn_f. sex_n sex_F. blockn yesno.;
	estimate 'Beta0' intercept 1 / exp;
	estimate 'B_group' groupn 1 / exp;
	estimate 'B_sex' sex_n 1 / exp;
	estimate 'B_proc_length' proc_length_center 1 / exp;
run;
title;

* 	Pseudo-Vague Prior N(0,1);
title1 "Bayesian Unadjusted Model of Block";
title2 "Normal (0, 1) Priors";
title3 "Results on the Odds Ratio scale";
proc mcmc data=data.trial
	plots=all
	diag=none   
	nbi=1000 
	nmc=9000 
	seed=123 
	dic
	plots(smooth)=all 
	statistics=all
	monitor=(or_beta0 or_groupn or_sex or_proc_length) 
	outpost=post1;
	parms (beta0 b_groupn b_sex B_proc_length_center) 0;
	prior beta0 b_sex B_proc_length_center ~ normal(0, var=1);
	prior b_groupn ~ normal(0, var=1);	
		beginnodata;
		or_beta0 = exp(beta0);
		or_groupn = exp(B_groupn);
		or_sex = exp(b_sex);
		or_proc_length = exp(B_proc_length_center);
		endnodata;  	
	p = logistic(beta0 + b_groupn*groupn + b_sex*sex_n + B_proc_length_center*proc_length_center);
   	model blockn ~ binary(p);
run;

data Prob; 
	set post1; 
	Indicator = (or_groupn < 1); *posterior probability that odds of group <1;
	label Indicator ='exp(druggrp)<1';
run;

proc means data=prob; *calculate the posterior probability (i.e., mean of indicator) and obtain OR estimates for parameters;
	var Indicator or_groupn or_beta0 or_sex or_proc_length;
run;
title;



* 	Vague Prior N(0,10);
title1 "Bayesian Adjusted Model of Block";
title2 "Normal (0,10) Priors";
title3 "Results on the Odds Ratio scale";
proc mcmc data=data.trial
	diag=none   
	nbi=1000
	nmc=9000 
	seed=123 
	dic
	plots(smooth)=all 
	statistics=all
	monitor=(or_beta0 or_groupn or_sex or_proc_length) 
	outpost=post2;
	parms (beta0 b_groupn b_sex B_proc_length_center) 0;
	prior beta0 b_groupn b_sex B_proc_length_center ~ normal(0, var=10);
		beginnodata;
		or_beta0 = exp(beta0);
		or_groupn = exp(B_groupn);
		or_sex = exp(b_sex);
		or_proc_length = exp(B_proc_length_center);
		endnodata;  	
	p = logistic(beta0 + b_groupn*groupn + b_sex*sex_n + B_proc_length_center*proc_length_center);
   	model blockn ~ binary(p);
run;

data Prob; 
	set post2; 
	Indicator = (or_groupn < 1); *posterior probability that odds of group <1;
	label Indicator ='exp(druggrp)<1';
run;

proc means data=prob; *calculate the posterior probability (i.e., mean of indicator) and obtain OR estimates for parameters;
	var Indicator or_groupn or_beta0 or_sex or_proc_length;
run;
title;



*	Optimistic Informative Prior N(-0.5, 2);
title1 "Bayesian Adjusted Model";
title2 "Normal (-0.5, 2) for 'Group' and N(0,10) Priors";
title3 "Results on the Odds Ratio scale";
proc mcmc data=data.trial
	plots=all
	nbi=1000 
    nmc=9000
	seed=123 
	dic 	
	plots(smooth)=all 
	statistics=all
	monitor=(or_beta0 or_groupn or_sex or_proc_length) 
	outpost=post3;
	parms (beta0 b_groupn b_sex B_proc_length_center) 0;
	prior beta0 b_sex B_proc_length_center ~ normal(mean=0, var=10);
	prior b_groupn ~ normal(mean=-0.5, var=2);
		beginnodata;
		or_beta0 = exp(beta0);
		or_groupn = exp(B_groupn);
		or_sex = exp(b_sex);
		or_proc_length = exp(B_proc_length_center);
		endnodata;  	
	p = logistic(beta0 + b_groupn*groupn + b_sex*sex_n + B_proc_length_center*proc_length_center);
   	model blockn ~ binary(p);
	format groupn groupn_f. sex_n sex_F. blockn yesno.;
run;

data Prob;
	set post3; 
	Indicator = (or_groupn < 1); *posterior probability that odds of group <1;
	label Indicator ='exp(druggrp)<1';
run;

proc means data=prob; *calculate the posterior probability (i.e., mean of indicator) and obtain OR estimates for parameters;
	var Indicator or_groupn or_beta0 or_sex or_proc_length;
run;
title;


* 	Skeptical Prior N(0, 2);
title1 "Bayesian Adjusted Model";
title2 "Normal (0, 2) for 'Group' and N(0,10) Priors";
title3 "Results on the Odds Ratio scale";
proc mcmc data=data.trial
	plots=all
	nbi=1000 
    nmc=9000
	seed=123 
	dic 	 
	plots(smooth)=all 
	statistics=all
	monitor=(or_beta0 or_groupn or_sex or_proc_length) 
	outpost=post4;
	parms (beta0 b_groupn b_sex B_proc_length_center) 0;
	prior beta0 b_sex B_proc_length_center ~ normal(mean=0, var=10);
	prior b_groupn ~ normal(mean=0, var=2);
		beginnodata;
		or_beta0 = exp(beta0);
		or_groupn = exp(B_groupn);
		or_sex = exp(b_sex);
		or_proc_length = exp(B_proc_length_center);
		endnodata;  	
	p = logistic(beta0 + b_groupn*groupn + b_sex*sex_n + B_proc_length_center*proc_length_center);
   	model blockn ~ binary(p);
	format groupn groupn_f. sex_n sex_F. blockn yesno.;
run;

data Prob; 
	set post4; 
	Indicator = (or_groupn < 1); *posterior probability that odds of group <1;
	label Indicator ='exp(druggrp)<1';
run;

proc means data=prob; *calculate the posterior probability (i.e., mean of indicator) and obtain OR estimates for parameters;
	var Indicator or_groupn or_beta0 or_sex or_proc_length;
run;
title;


ods graphics off;
run; quit;
/* END OF PROGRAM */ ;
