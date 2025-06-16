######################################################################################################
#  R-code to replicate the simulation study outlined in Section 3 of the paper (Scenario III: NPH and nonlinear)
#
#   Penalised spline estimation of covariate-specific time-dependent ROC curves
#
#  by Maria Xose Rodriguez Alvarez and Vanda Inacio
#
#  Contact: mxrodriguez@uvigo.gal
#			vanda.inacio@ed.ac.uk			
#
######################################################################################

rm(list = ls())

##########################################################################
# Load the needed packages
##########################################################################
library(CondTimeROC)	# Available at https://github.com/mxrodriguezUVIGO/CondTimeROC
library(ggplot2)		# Graphical library

##########################################################################
# Load the needed functions
##########################################################################
source("r-code-functions-simulations-cTimeROC.R")	# Available at https://github.com/mxrodriguezUVIGO/CondTimeROC-paper

#########################################################################################################
# Simulation Scenario III (NPH and nonlinear): Only one covariate value and one time point
##########################################################################################################

	#####################################################
	# Input parameters
	#####################################################
		# Number of simulations
		n.sim <- 500

		# Cutoffs at which to compute the FPR and TPR
		cutoffs <- seq(-3, 5, l = 200)

		# FPR at  at which to compute the ROC curve
		seq.fpr <- seq(0, 1, l = 101)  # 

		# Formulas for the models
			# Song & Zhou
			formula.sz <- list(coxph = "Surv(t, delta) ~ y + x", lm = "y ~ x")

			# New approach
			formula.ps <- list(hazard = "Surv(t, delta) ~ ti(t, k = 8, bs = 'cr') + ti(x, k = 8, bs = 'cr') + ti(y, k = 8, bs = 'cr') + 
					ti(x, y, k = 8, bs = 'cr') +
					ti(x, t, k = 8, bs = 'cr') +
					ti(y, t, k = 8, bs = 'cr')", 
					biomarker = list(mean = "y ~ s(x, k = 13, bs = 'cr')", var = "~ s(x, k = 13, bs = 'cr')"))

	#####################################################
	# n = 300
	#####################################################
		m3.300.x0 <- cTimeROC_simulation(n.sim = n.sim, 
			n = 300, 
			cov.value = 0, 
			predict.time = 0.24,
			scenario = "m3", 
			cutoffs = cutoffs, 
			seq.fpr = seq.fpr, 
			formula.sz = formula.sz, 
			formula.ps = formula.ps,
			select = TRUE,  
			smooth = "ns", 
			bwsel = "npbw",
			seed = NULL)

		# Plot results
		plot(m3.300.x0)

		# If the goal is to estimate at different covariate values and time points,
		# this can be done as follows. It saves a lot of time since the models do not need to be refitted.
		m3.300.x2 	<- cTimeROC_simulation(n.sim = n.sim, 
			n = 300, 
			cov.value = 2, 
			predict.time = 0.24,
			scenario = "m3", 
			cutoffs = cutoffs, 
			seq.fpr = seq.fpr, 
			formula.sz = formula.sz, 
			formula.ps = formula.ps,
			select = TRUE,  
			smooth = "ns", bwsel = "npbw",
			fitted.models.ps = m3.300.x0$fitted.models.ps, # Fitted models (hazard and biomarker) for the proposed approach
			lbd = m3.300.x0$lbd, 			# Bandwidths for the non-parametric approach				
			data.sim = m3.300.x0$data.sim) 	# Data used in the original simulation


	#####################################################
	# n = 600
	#####################################################
		m3.600.x0 <- cTimeROC_simulation(n.sim = n.sim, 
			n = 600, 
			cov.value = 0, 
			predict.time = 0.24,
			scenario = "m3", 
			cutoffs = cutoffs, 
			seq.fpr = seq.fpr, 
			formula.sz = formula.sz, 
			formula.ps = formula.ps,
			select = TRUE,  
			smooth = "ns", bwsel = "npbw",
			seed = NULL)