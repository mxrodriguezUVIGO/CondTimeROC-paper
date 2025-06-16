#######################################################################################################
#  R-code containing the necessary functions to replicate the simulation study outlined in Section 3 of the paper.
#
#   Penalised spline estimation of covariate-specific time-dependent ROC curves
#
#  by Maria Xose Rodriguez Alvarez and Vanda Inacio
#
#  Contact: mxrodriguez@uvigo.gal
#			vanda.inacio@ed.ac.uk			
#
######################################################################################
################################################################################
# Some auxiliary functions
################################################################################
	# Time-varying effects
	time.effects <- function(t, x, y, scenario) {
		if(scenario == "m1") { # Linear on x and y and PH
			res <- 2 + log(t + 0.2)
		} else if(scenario == "m2") { # Non-linear on x and y and PH
			res <- 2 + log(t + 0.2)
		} else if(scenario == "m3") { # Non-linear on x and y and NPH with marker
			res <- (2 + log(t + 0.2))*(y^3/20) # NPH with marker
		} else {
			error("Scenario not implemented")
		}
		res
	}

	# No time-varying effects
	notime.effects <- function(x, y, scenario) {
		if(scenario == "m1") { # Linear on x and y and PH
			res <- y + 0.1*x
		} else if(scenario == "m2") { # Non-linear on x and y and PH
			res <- y^3/20 + 0.5*sin(2*(x + 1.5))
		} else if(scenario == "m3") { # Non-linear on x and y and NPH with marker
			res <- 0.5*sin(2*(x + 1.5))
		} else {
			error("Scenario not implemented")
		}
		res
	}

	obtain.capital.lambda <- function(t ,x, y, scenario){
		res <- sapply(t, function(tt, x, y, scenario) {
				seq.t.integral <- seq(0, tt, l = 101)
				seq.inter.time.cov.mar  <- exp(time.effects(seq.t.integral, x, y, scenario))
				exp(notime.effects(x, y, scenario))*CondTimeROC:::simpson(seq.inter.time.cov.mar, seq.t.integral)
			}, x = x, y = y, scenario = scenario)
		res <- c(res)
		res
	}

	root.finding.function <- function(log_t, x, y, u, scenario) {
		res <- exp(-obtain.capital.lambda(exp(log_t), x, y, scenario)) - u
		res

	}

	# Function to delete unnecesary info from a mgcv object (save space)
	remove_info <- function(x) {
		obj <- c("linear.predictors",
			"fitted.values",
			"edf",
			"edf1",
			"Vp",
			"Ve",
			"Vc",
			"wt",
			"y",
			"prior.weights",
			"weights",
			"residuals", "dinfo", "db.drho", "offset", "R")

		for(j in 1:length(obj)) {		
			x[[obj[j]]] <- NULL
		}
		# Keep only some elements of the model
		x$model <- x$model[1:2,]
		x
	}

################################################################################
# Generate simulated data according to the desired scenario (m1, m2 and m3)
################################################################################
	sim.data <- function(n, scenario = "m1") {
		x <- rnorm(n, 1, 1)
		y <- rnorm(n, x, 1)
		u <- runif(n)

		t <- sapply(1:n, function(i, x, y, u, scenario) {
			exp(uniroot(root.finding.function, interval = c(-20, 20), extendInt = "yes", x = x[i], y = y[i], u = u[i], scenario = scenario)$root)
		}, x = x, y = y, u = u, scenario = scenario)

		if(scenario == "m1") {
			c <- pmin(rexp(n, rate = 1/(0.01 + 0.2*abs(x))), 20)
		} else if(scenario == "m2") {
			c <- pmin(rexp(n, rate = 1/(0.01 + 0.4*abs(x))), 20)
		} else if(scenario == "m3") {
			c <- pmin(rexp(n, rate = 1/(0.01 + abs(x))), 20)
		} else {
			error("Scenario not implemented")
		}	

		delta <- as.numeric(t <= c)
		tc <- pmin(t,c)
		data <- data.frame(x = x, y = y, t = tc, t_orig = t, c = c, delta = delta)
		data
	}

################################################################################
# Obtain theoretical quantities (FPF, TPF, ROC curve and AUC)
################################################################################
	theor.data <- function(cutoffs = NULL, seq.fpr = seq(0, 1, l = 100), cov.value, predict.time, scenario = 1) {
		f.tpr <- function(u, t, x, scenario) {
			# uu - x since Y = X + epsilon
			res <- sapply(u, function(uu, x, t, scenario) {
				aux <- obtain.capital.lambda(t, x, uu, scenario)
				(1 - exp(-aux))*dnorm(uu-x)
			}, x = x, t = t, scenario = scenario)
			res <- c(res)
			res
		}
		f.fpr <- function(u, t, x, scenario) {
			# uu - x since Y = X + epsilon
			res <- sapply(u, function(uu, x, t,scenario) {
				aux <- obtain.capital.lambda(t, x, uu,scenario)
				exp(-aux)*dnorm(uu-x)
			}, x = x, t = t, scenario = scenario)
			res <- c(res)
			res
		}
		AUC.computation <- function(u, fpr, tpr) {
			res <- approxfun(fpr, tpr, yleft = min(tpr), yright = max(tpr))(u)
		}
		
		cutoffs <- yy <- if(is.null(cutoffs)) {
			seq(-3, 5, l = 50)
		} else {
			cutoffs
		}
		tpr.teor <- fpr.teor <- vector(l = length(cutoffs))
		ROC <- vector(l = length(seq.fpr))

		fpr.den <- integrate(f.fpr, -Inf, Inf, t = predict.time, x = cov.value, scenario = scenario)
		tpr.den <- integrate(f.tpr, -Inf, Inf, t = predict.time, x = cov.value, scenario = scenario)
			
		for(i in 1:length(cutoffs)) {
			tpr.teor[i] <- integrate(f.tpr, yy[i], Inf, t = predict.time, x = cov.value, scenario = scenario)$value/tpr.den$value
			fpr.teor[i] <- integrate(f.fpr, yy[i], Inf, t = predict.time, x = cov.value, scenario = scenario)$value/fpr.den$value
		}
		ROC <- approxfun(fpr.teor, tpr.teor, yleft = min(tpr.teor), yright = max(tpr.teor))(seq.fpr)
		yy  <- seq(-3,5,l = 500)
		tpr.AUC <- fpr.AUC <- vector(l = length(yy))
		for(i in 1:length(yy)) {
			tpr.AUC[i] <- integrate(f.tpr, yy[i], Inf, t = predict.time, x = cov.value, scenario = scenario)$value/tpr.den$value
			fpr.AUC[i] <- integrate(f.fpr, yy[i], Inf, t = predict.time, x = cov.value, scenario = scenario)$value/fpr.den$value
		}
		AUC <- integrate(AUC.computation, 0, 1, fpr = fpr.AUC, tpr = tpr.AUC)
		res <- list(FPR = fpr.teor, TPR = tpr.teor, ROC = ROC, AUC = AUC)
		res
	}

################################################################################
# Function for the simulations associated with the paper:
#
# Arguments:
#   n.sim              Number of runs.
#   n                  Sample size.
#   cov.value          Covariate value at which to compute the time-dependent ROC curve.
#   predict.time       Time at which to compute the time-dependent ROC curve.
#   scenarios          Simulation scenario: m1 (I), m2 (II), m3 (III).
#   cutoffs            Set of cutoffs where the FPR and TPR will be computed. 
#                      If NULL, the cutoffs are the unique values of the marker in the sample.
#   seq.fpr            Set of FPR values at which the time-dependent ROC curve will be computed.
#   formula.ps         For the new approach, a named list with two elements: 
#                         - "hazard": formula for the model of the conditional hazard function
#                         - "biomarker": named list with formulas for the location-scale regression 
#                           model of the biomarker (mean and variance)
#   n.timepoints       For the new approach and the piecewise exponential model, number of time points.
#   select             For the new approach, a logical indicating if the double penalty approach for model building/selection described in Marra and Wood (2011) is implemented. 
#   fitted.models.ps   For the new approach, a list of length n.sim, where each element is a named list 
#                      containing the fitted models for the conditional hazard function ("hazard") 
#                      and the biomarker ("biomarker"). 
#                      Used to save time if the aim is to estimate at different covariate values and/or time points.
#   formula.sz         For the semiparametric approach, a named list with two elements:
#                         - "coxph": formula for the Cox model
#                         - "lm": formula for the linear model
#   smooth             For the nonparametric approach, estimator:
#                         - "nu": unsmoothed estimator
#                         - "ns": smoothed estimator
#   bwsel              For the nonparametric approach, bandwidth selector:
#                         - "ALbw": Altman and Leger (1995)
#                         - "dpik": Sheather and Jones (1991)
#                         - "npbw": Li et al. (2013)
#                         - "NNE": Nearest neighbour
#   ipcw               For the nonparametric approach, inverse probability weights computed 
#                      using the Beran estimator (1981) or the Akritas estimator (1994).
#   span               For the nonparametric approach, span for NNE (also needed when using 
#                      IPCW based on Akritas' approach).
#   lambda             For the nonparametric approach, smoothing parameter for NNE 
#                      (also needed when using IPCW based on Akritas' approach). 
#                      Either lambda or span is required for NNE.
#   kernel             For the nonparametric approach, kernel function: either "gaussian" or "epanechnikov".
#   lbd                For the nonparametric approach, bandwidths for the simulations. 
#                      By default, they are computed during the simulation, but it is also possible 
#                      to provide their values. In this case, supply a matrix of dimension 2 x n.sim 
#                      with the bandwidths for the covariate and the marker for each run. 
#                      Used to save time if the aim is to estimate at different covariate values and/or time points.
#   data.sim           A list of length n.sim, where each element is a data frame. 
#                      Used to save time if the aim is to estimate at different covariate values and/or time points.
#   seed               Seed for the simulations. Default is NULL (then set to 123).
################################################################################
	cTimeROC_simulation <- function(n.sim = 500, 
		n = 300, 
		cov.value, 
		predict.time, 
		scenario = "m1", 
		cutoffs = NULL, 
		seq.fpr = seq(0, 1, length = 101),			
		formula.ps = list(hazard = "Surv(t, delta) ~ ti(t, k = 8, bs = 'cr') + ti(x, k = 8, bs = 'cr') + ti(y, k = 8, bs = 'cr') + 
				ti(x, y, k = 8, bs = 'cr') +
				ti(x, t, k = 8, bs = 'cr') +
				ti(y, t, k = 8, bs = 'cr')", 
				biomarker = list(mean = "y ~ s(x, k = 13, bs = 'cr')", var = "~ s(x, k = 13, bs = 'cr')")),
		n.timepoints = NULL,
		select = TRUE,
		fitted.models.ps = NULL,
		formula.sz = list(coxph = "Surv(t, delta) ~ y + x", lm = "y ~ x"), 
		smooth = c("nu","ns"), 
		bwsel = c("ALbw","dpik","npbw", "NNE"), 
		ipcw = c("Beran", "Akritas"), 
		span = NULL, 
		lambda = NULL, 
		kernel = c("gaussian", "epanechnikov"), 
		lbd = NULL, 	
		data.sim = NULL, 
		seed = NULL) {
		
		form_biomarker  <- formula.ps$biomarker
		form_hazard  <- as.formula(formula.ps$hazard)

		if(is.null(seed)) seed <- 123
		set.seed(seed)

		smooth <- match.arg(smooth)
		bwsel <- match.arg(bwsel)
		kernel <- match.arg(kernel)
		ipcw <- match.arg(ipcw)
		
		if(!is.null(cutoffs)) {	
			FPR.ps <- TPR.ps <- FPR.sz <- TPR.sz <- FPR.np <- TPR.np <- matrix(NA, ncol = n.sim, nrow = length(cutoffs) + 2)
		}
		ROC.ps <- ROC.sz <- ROC.np <- matrix(NA, ncol = n.sim, nrow = length(seq.fpr))
		AUC.ps <- AUC.sz <- AUC.np <- vector(l = n.sim)

		conv.ps <- vector(l = n.sim)
		
		# Smoothing parameters non-parametric approach
		lbd.temp <- if(is.null(lbd)) {
			matrix(ncol = n.sim, nrow = 2)
		} else {
			lbd
		}

		# Fitted models penalised-spline based method
		fitted.models.ps.temp <- if(is.null(fitted.models.ps)) {
			list()
		} else {
			fitted.models.ps
		}

		theor.model <- theor.data(cutoffs, seq.fpr, cov.value, predict.time, scenario = scenario)
		
		if(is.null(data.sim)) {
			data.sim <- list()
			for(i in 1:n.sim) {
				data.sim[[i]] <- sim.data(n = n, scenario = scenario)
			}
		}
		for(i in 1:n.sim) {
			print(i)
			data <- data.sim[[i]]
			# Song and Zhou
				op <- options(warn = 1)
				roc.sz <- ctimeROC.sp(formula.coxph = formula.sz$coxph, 
					formula.lm = formula.sz$lm, 
					data = data, predict.time = predict.time, cov.values = data.frame(x = cov.value), cutoffs = cutoffs)
				AUC.sz[i]<- roc.sz$AUC
				ROC.sz[,i] <- CondTimeROC:::obtain.ROCCurve(seq.fpr, rev(roc.sz$FPR), rev(roc.sz$TPR))$ROC 

			# Penalised-based approach (new)
				# Warning to error (convergence)
				op <- options(warn = 2)
				if(is.null(fitted.models.ps)) {
					roc.ps <- try(ctimeROC.ps(formula.hazard = form_hazard, formula.biomarker = form_biomarker, 
						data = data, n.timepoints = n.timepoints, predict.time = predict.time, 
						cov.values = data.frame(x = cov.value), cutoffs = cutoffs, select = select))
					
					if(inherits(roc.ps, "try-error")) {
						conv.ps[i] <- 0
						fitted.models.ps.temp[[i]] <- NULL
					} else {
						conv.ps[i] <- 1
						fitted.models.ps.temp[[i]] <- list(hazard = remove_info(roc.ps$fitted.models$hazard), 
							biomarker = roc.ps$fitted.models$biomarker)
					}
					
				} else {
					if(!is.null(fitted.models.ps.temp[[i]])) {
						roc.ps <- ctimeROC.ps(formula.hazard = form_hazard, formula.biomarker = form_biomarker, 
							data = data, n.timepoints = n.timepoints, predict.time = predict.time, 
							cov.values = data.frame(x = cov.value), cutoffs = cutoffs, select = select,
							fitted.models = fitted.models.ps.temp[[i]])
						conv.ps[i] <- 1
					} else {
						conv.ps[i] <- 0
					}
				}
				if(conv.ps[i]) {
					AUC.ps[i]<- roc.ps$AUC
					ROC.ps[,i] <- CondTimeROC:::obtain.ROCCurve(seq.fpr, rev(roc.ps$FPR), rev(roc.ps$TPR))$ROC
				} else {
					AUC.ps[i] <- NA
					ROC.ps[,i] <- NA
				}
			
			# Non parametric
				if(is.null(lbd)) {							
					roc.np <- ctimeROC.np(marker = "y", covariate = "x", time = "t", status = "delta", 
						data = data, predict.time = predict.time, cov.value = cov.value, cutoffs = cutoffs, 
						smooth = smooth, bwsel = bwsel, span = span, lambda = lambda, kernel = kernel, ipcw = ipcw, 
						AUC.calc = "cexpr")
					lbd.temp[,i] <- roc.np$lbd
				} else {
					roc.np <- ctimeROC.np(marker = "y", covariate = "x", time = "t", status = "delta", 
						data = data, predict.time = predict.time, cov.value = cov.value, cutoffs = cutoffs, 
						smooth = smooth, bwsel = bwsel, lbd.x = lbd[1,i], lbd.y = lbd[2,i], span = span, 
						lambda = lambda, kernel = kernel, ipcw = ipcw, AUC.calc = "cexpr")
				}
				AUC.np[i]<- roc.np$AUC
				ROC.np[,i] <- CondTimeROC:::obtain.ROCCurve(seq.fpr, rev(roc.np$FPR), rev(roc.np$TPR))$ROC

			if(!is.null(cutoffs)) {
				if(conv.ps[i]) {	
					FPR.ps[,i] <- roc.ps$FPR
					TPR.ps[,i] <- roc.ps$TPR
				} else {
					FPR.ps[,i] <- NA
					TPR.ps[,i] <- NA
				}
				FPR.sz[,i] <- roc.sz$FPR
				TPR.sz[,i] <- roc.sz$TPR
				FPR.np[,i] <- roc.np$FPR
				TPR.np[,i] <- roc.np$TPR
			}
			if(conv.ps[i]) rm(roc.ps)
			rm(roc.sz)
			rm(roc.np)
			gc()
		}
		res <- if(!is.null(cutoffs)) {
			list(n.sim = n.sim, 
				cov.value = cov.value, 
				predict.time = predict.time,
				cutoffs = cutoffs, 
				seq.fpr = seq.fpr, 
				theor = theor.model, 
				FPR = list(ps = FPR.ps, sz = FPR.sz, np = FPR.np), 
				TPR = list(ps = TPR.ps, sz = TPR.sz, np = TPR.np), 
				ROC = list(ps = ROC.ps, sz = ROC.sz, np = ROC.np), 
				AUC = list(ps = AUC.ps, sz = AUC.sz, np = AUC.np), 
				lbd = lbd.temp, 
				fitted.models.ps = fitted.models.ps.temp,
				conv = list(ps = conv.ps), 
				data.sim = data.sim)
		} else {
			list(n.sim = n.sim,
				cov.value = cov.value, 
				predict.time = predict.time, 
				seq.fpr = seq.fpr, 
				theor = theor.model, 
				ROC = list(ps = ROC.ps, sz = ROC.sz, np = ROC.np), 
				AUC = list(ps = AUC.ps, sz = AUC.sz, np = AUC.np), 
				lbd = lbd.temp, 
				fitted.models.ps = fitted.models.ps.temp,
				conv = list(ps = conv.ps), 
				data.sim = data.sim)
		}
		class(res) <- "cTimeROC.Simulation"
		op <- options(warn = 1)
		res
		
	}
################################################################################
# Plot function for the simulations
################################################################################
	plot.cTimeROC.Simulation <- function(x, ...) {
		n.sim <- x$n.sim
		cov.value <- x$cov.value
		predict.time <- x$predict.time
		n.fpf <- length(x$seq.fpr)
		seq.fpr <- x$seq.fpr

		############################################
		# ROC curve (RMSE)
		############################################
			RMSE.ps <- apply(x$ROC$ps, 2, function(y) {sqrt(mean((y - x$theor$ROC)^2, na.rm = TRUE))})
			RMSE.sz <- apply(x$ROC$sz, 2, function(y) {sqrt(mean((y - x$theor$ROC)^2, na.rm = TRUE))})
			RMSE.np <- apply(x$ROC$np, 2, function(y) {sqrt(mean((y - x$theor$ROC)^2, na.rm = TRUE))})

			df <- data.frame(RMSE = c(RMSE.ps, RMSE.sz, RMSE.np),
				predict.time = rep(rep(predict.time, each = n.sim), 3),
				cov.value = rep(rep(cov.value, each = n.sim), 3),
				Approach = rep(c("New approach","Semi-parametric","Non-parametric"), each = n.sim))

			df$Approach <- factor(df$Approach, levels = c("New approach","Semi-parametric", "Non-parametric"), 
				labels = c("New approach","Semi-parametric", "Non-parametric"), ordered = TRUE)

			g_roc_rmse <- ggplot(df, aes(x = Approach, y = RMSE)) +  
				geom_violin(aes(colour = Approach, fill = Approach)) +  
				geom_boxplot(width = 0.2) +
				scale_fill_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +
				scale_color_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +			
				labs(title = "Empirical Root Mean Squared Error (ERMSE) \n Time-dependent ROC curve", x = "", y = "ERMSE")  + 
				theme_bw() +
				theme(legend.position = 'bottom',
					strip.text.y = element_text(size = 15), 
					strip.text.x = element_text(size = 15), 
					axis.text.x = element_blank(), 
					plot.title = element_text(hjust = 0.5, size = 15), 
					axis.text = element_text(size = 10), 
					axis.title = element_text(size = 15), 
					legend.title = element_text(size = 15), 
					legend.text = element_text(size = 15))
			print(g_roc_rmse)
		
		############################################
		# ROC curve
		############################################
			ROC.ps.m <- apply(x$ROC$ps, 1, mean, na.rm = TRUE)
			ROC.ps.ql <- apply(x$ROC$ps, 1, quantile, 0.025, na.rm = TRUE)
			ROC.ps.qh <- apply(x$ROC$ps, 1, quantile, 0.975, na.rm = TRUE)

			ROC.sz.m <- apply(x$ROC$sz, 1, mean, na.rm = TRUE)
			ROC.sz.ql <- apply(x$ROC$sz, 1, quantile, 0.025, na.rm = TRUE)
			ROC.sz.qh <- apply(x$ROC$sz, 1, quantile, 0.975, na.rm = TRUE)

			ROC.np.m <- apply(x$ROC$np, 1, mean, na.rm = TRUE)
			ROC.np.ql <- apply(x$ROC$np, 1, quantile, 0.025, na.rm = TRUE)
			ROC.np.qh <- apply(x$ROC$np, 1, quantile, 0.975, na.rm = TRUE)

			df <- data.frame(theor = x$theor$ROC,
				ROC.m = c(ROC.ps.m, ROC.sz.m, ROC.np.m),
				ROC.ql = c(ROC.ps.ql, ROC.sz.ql, ROC.np.ql),
				ROC.qh = c(ROC.ps.qh, ROC.sz.qh, ROC.np.qh),
				seq.fpr = rep(seq.fpr, 3),
				Approach = rep(c("New approach","Semi-parametric","Non-parametric"), each = n.fpf))

			df$Approach <- factor(df$Approach, levels = c("New approach","Semi-parametric", "Non-parametric"), 
				labels = c("New approach","Semi-parametric", "Non-parametric"), ordered = TRUE)
		
			g_roc <- ggplot(df, aes(x = seq.fpr, y = theor)) + 
				geom_line(linewidth = 1) +
				geom_line(aes(x = seq.fpr, y = ROC.m, color = Approach), linewidth = 1, linetype = "dashed") +
				scale_color_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +
				geom_ribbon(aes(x = seq.fpr, y = ROC.m, ymin = ROC.ql, ymax = ROC.qh, fill = Approach, color = Approach), linewidth = 0.2, alpha = 0.2) +
				scale_fill_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +
				labs(title = "Time-dependent ROC curve", x = "1-Specificity", y = "Sensitivity") +
				ylim(0,1) + 
				theme_bw() +			
				theme(legend.position = 'bottom',
					strip.text.y = element_text(size = 15), 
					strip.text.x = element_text(size = 15), 
					axis.text.x = element_blank(), 
					plot.title = element_text(hjust = 0.5, size = 15), 
					axis.text = element_text(size = 10), 
					axis.title = element_text(size = 15), 
					legend.title = element_text(size = 15), 
					legend.text = element_text(size = 15))
			
			x11()	
			print(g_roc)

		############################################
		# Bias AUC (Violin plot)
		############################################
			bias.ps <- x$AUC$ps - x$theor$AUC$value
			bias.sz <- x$AUC$sz - x$theor$AUC$value
			bias.np <- x$AUC$np - x$theor$AUC$value

			df <- data.frame(bias = c(bias.ps, bias.sz, bias.np),
				predict.time = rep(rep(predict.time, each = n.sim), 3),
				cov.value = rep(rep(cov.value, each = n.sim), 3),
				Approach = rep(c("New approach","Semi-parametric","Non-parametric"), each = n.sim))

			df$Approach <- factor(df$Approach, levels = c("New approach","Semi-parametric", "Non-parametric"), 
				labels = c("New approach","Semi-parametric", "Non-parametric"), ordered = TRUE)

			g_auc <- ggplot(df, aes(x = Approach, y = bias)) +  
				geom_violin(aes(colour = Approach, fill = Approach)) +  
				geom_boxplot(width = 0.2) +
				scale_fill_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +
				scale_color_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +
				geom_hline(yintercept = 0, col = "black", linetype = 2) +
				labs(title = "Bias time-dependent AUC", x = "", y = "Bias")  + 
				theme_bw() +
				theme(legend.position = 'bottom',
					strip.text.y = element_text(size = 15), 
					strip.text.x = element_text(size = 15), 
					axis.text.x = element_blank(), 
					plot.title = element_text(hjust = 0.5, size = 15), 
					axis.text = element_text(size = 10), 
					axis.title = element_text(size = 15), 
					legend.title = element_text(size = 15), 
					legend.text = element_text(size = 15))
			x11()
			print(g_auc)

		#######################################################
		# TPF (Sensitivity) and FPF (1-Specificty): RMSE
		#######################################################	
			if(!is.null(x$cutoffs)) {
		
				RMSE.TPR.ps <- apply(x$TPR$ps[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y - x$theor$TPR)^2, na.rm = TRUE))})
				RMSE.TPR.sz <- apply(x$TPR$sz[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y - x$theor$TPR)^2, na.rm = TRUE))})
				RMSE.TPR.np <- apply(x$TPR$np[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y - x$theor$TPR)^2, na.rm = TRUE))})

				df <- data.frame(RMSE = c(RMSE.TPR.ps, RMSE.TPR.sz, RMSE.TPR.np),
					predict.time = rep(rep(predict.time, each = n.sim), 3),
					cov.value = rep(rep(cov.value, each = n.sim), 3),
					Approach = rep(c("New approach","Semi-parametric","Non-parametric"), each = n.sim))

				df$Approach <- factor(df$Approach, levels = c("New approach","Semi-parametric", "Non-parametric"), 
					labels = c("New approach","Semi-parametric", "Non-parametric"), ordered = TRUE)

				g_TPR_rmse <- ggplot(df, aes(x = Approach, y = RMSE)) +  
					geom_violin(aes(colour = Approach, fill = Approach)) +  
					geom_boxplot(width = 0.2) +
					scale_fill_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +
					scale_color_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +			
					labs(title = "Empirical Root Mean Squared Error (ERMSE) \n Time-dependent Sensitivity", x = "", y = "ERMSE")  + 
					theme_bw() +
					theme(legend.position = 'bottom',
						strip.text.y = element_text(size = 15), 
						strip.text.x = element_text(size = 15), 
						axis.text.x = element_blank(), 
						plot.title = element_text(hjust = 0.5, size = 15), 
						axis.text = element_text(size = 10), 
						axis.title = element_text(size = 15), 
						legend.title = element_text(size = 15), 
						legend.text = element_text(size = 15))
				x11()	
				print(g_TPR_rmse)

				RMSE.FPR.ps <- apply(x$FPR$ps[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y - x$theor$FPR)^2, na.rm = TRUE))})
				RMSE.FPR.sz <- apply(x$FPR$sz[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y - x$theor$FPR)^2, na.rm = TRUE))})
				RMSE.FPR.np <- apply(x$FPR$np[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y - x$theor$FPR)^2, na.rm = TRUE))})

				df <- data.frame(RMSE = c(RMSE.FPR.ps, RMSE.FPR.sz, RMSE.FPR.np),
					predict.time = rep(rep(predict.time, each = n.sim), 3),
					cov.value = rep(rep(cov.value, each = n.sim), 3),
					Approach = rep(c("New approach","Semi-parametric","Non-parametric"), each = n.sim))

				df$Approach <- factor(df$Approach, levels = c("New approach","Semi-parametric", "Non-parametric"), 
					labels = c("New approach","Semi-parametric", "Non-parametric"), ordered = TRUE)

				g_FPR_rmse <- ggplot(df, aes(x = Approach, y = RMSE)) +  
					geom_violin(aes(colour = Approach, fill = Approach)) +  
					geom_boxplot(width = 0.2) +
					scale_fill_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +
					scale_color_manual(values = c("#D55E00", "#009E73", "#56B4E9")) +			
					labs(title = "Empirical Root Mean Squared Error (ERMSE) \n Time-dependent 1-Specificity", x = "", y = "ERMSE")  + 
					theme_bw() +
					theme(legend.position = 'bottom',
						strip.text.y = element_text(size = 15), 
						strip.text.x = element_text(size = 15), 
						axis.text.x = element_blank(), 
						plot.title = element_text(hjust = 0.5, size = 15), 
						axis.text = element_text(size = 10), 
						axis.title = element_text(size = 15), 
						legend.title = element_text(size = 15), 
						legend.text = element_text(size = 15))
				x11()
				print(g_FPR_rmse)
		 	}
	}
################################################################################
# Print function for the simulations
################################################################################
	print.cTimeROC.Simulation <- function(x, ...) {
		mse.ROC.ps <- apply(x$ROC$ps, 2, function(y) {sqrt(mean((y-x$theor$ROC)^2))})
		mse.ROC.sz <- apply(x$ROC$sz, 2, function(y) {sqrt(mean((y-x$theor$ROC)^2))})
		mse.ROC.np <- apply(x$ROC$np, 2, function(y) {sqrt(mean((y-x$theor$ROC)^2))})

		if(!is.null(x$cutoffs)) {
			mse.FPR.ps <- apply(x$FPR$ps[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y-x$theor$FPR)^2))})
			mse.FPR.sz <- apply(x$FPR$sz[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y-x$theor$FPR)^2))})
			mse.FPR.np <- apply(x$FPR$np[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y-x$theor$FPR)^2))})
		
			mse.TPR.ps <- apply(x$TPR$ps[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y-x$theor$TPR)^2))})
			mse.TPR.sz <- apply(x$TPR$sz[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y-x$theor$TPR)^2))})
			mse.TPR.np <- apply(x$TPR$np[2:(length(x$cutoffs)+1),], 2, function(y) {sqrt(mean((y-x$theor$TPR)^2))})
		}
		AUC.ps <- x$AUC$ps - x$theor$AUC$value
		AUC.sz <- x$AUC$sz - x$theor$AUC$value
		AUC.np <- x$AUC$np - x$theor$AUC$value
		
		cat("PS approach \n")
		if(!is.null(x$cutoffs)) {
			cat(paste("RMSE - FPR:", round(mean(mse.FPR.ps, na.rm = TRUE)*100, 3), " (",round(sd(mse.FPR.ps, na.rm = TRUE)*100, 3),")\n", sep = ""))
			cat(paste("RMSE - TPR:", round(mean(mse.TPR.ps, na.rm = TRUE)*100, 3), " (",round(sd(mse.TPR.ps, na.rm = TRUE)*100, 3),")\n", sep = ""))
			cat(paste("RMSE - ROC:", round(mean(mse.ROC.ps, na.rm = TRUE)*100, 3), " (",round(sd(mse.ROC.ps, na.rm = TRUE)*100, 3),")\n", sep = ""))
		}
		cat(paste("AUC (Bias):", round(mean(AUC.ps, na.rm = TRUE)*100, 3), " (",round(sd(AUC.ps, na.rm = TRUE)*100, 3),")\n", sep = ""))


		cat("Song & Zhou approach \n")
		if(!is.null(x$cutoffs)) {
			cat(paste("RMSE - FPR:", round(mean(mse.FPR.sz, na.rm = TRUE)*100, 3), " (",round(sd(mse.FPR.sz, na.rm = TRUE)*100, 3),")\n", sep = ""))
			cat(paste("RMSE - TPR:", round(mean(mse.TPR.sz, na.rm = TRUE)*100, 3), " (",round(sd(mse.TPR.sz, na.rm = TRUE)*100, 3),")\n", sep = ""))
			cat(paste("RMSE - ROC:", round(mean(mse.ROC.sz, na.rm = TRUE)*100, 3), " (",round(sd(mse.ROC.sz, na.rm = TRUE)*100, 3),")\n", sep = ""))
		}
		cat(paste("AUC (Bias):", round(mean(AUC.sz, na.rm = TRUE)*100, 3), " (",round(sd(AUC.sz, na.rm = TRUE)*100, 3),")\n", sep = ""))

		cat("Non parametric \n")
		if(!is.null(x$cutoffs)) {
			cat(paste("RMSE - FPR:", round(mean(mse.FPR.np, na.rm = TRUE)*100, 3), " (",round(sd(mse.FPR.np, na.rm = TRUE)*100, 3),")\n", sep = ""))
			cat(paste("RMSE - TPR:", round(mean(mse.TPR.np, na.rm = TRUE)*100, 3), " (",round(sd(mse.TPR.np, na.rm = TRUE)*100, 3),")\n", sep = ""))
			cat(paste("RMSE - ROC:", round(mean(mse.ROC.np, na.rm = TRUE)*100, 3), " (",round(sd(mse.ROC.np, na.rm = TRUE)*100, 3),")\n", sep = ""))
		}
		cat(paste("AUC (Bias):", round(mean(AUC.np, na.rm = TRUE)*100, 3), " (",round(sd(AUC.np, na.rm = TRUE)*100, 3),")\n", sep = ""))
	}