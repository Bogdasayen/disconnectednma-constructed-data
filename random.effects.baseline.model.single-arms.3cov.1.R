# Random effects on baseline models with single-arm studies
# Uses 3 covariates to fit baseline models and incorporate single arms studies
# Covariates not used in the connecteded/pure RCT NMA (subsumed into nuisance parameters mu[i] anyway)
# Version currently uses a fixed intercept and only 1 covariate for baseline model but up to 3 may be included.
# Howard Thom 14-October-2017

# Changes from previous versions
# Uses only reference treatment arms to fit the baseline model
# Uses 'cut()' to prevent feedback from single-arm studies to baseline model

# Data are same as TSD format: ns, nt, na, r, n, t
# The data on baseline arms (from RCTs) are
# ns.base : number of baseline arms
# r.base : number of events in baseline arms
# n.base : number of patients in baseline arms
# x.base : Matrix of covariates with one row per covariate
# The data on single arm studies are:
# ns.single : number of single arm studies
# r.single : number of events in single arm studies
# n.single : number of patients in single arm studies
# t.single : treatment in single arm study
# x.single : matrix of covariates with one row per covariate
# cov.index : vector of numbers indicating which row of x.single/x.base to use for each covariate

# NOTE: All x.base and x.single must be defined. Suggest using mean covariate value when not reported.
# NOTE: Naively pools single arm and RCT evidence on treatment effects (not sure enough evidence to fit hierarchical model)



# Model with baseline prediction (fixed intercept, up to 3 covariates)
# Binomial likelihood, logit link 
# Simultaneous baseline and treatment effects model for multi-arm trials 
# Includes single-arm studies by predicting baseline response using baseline model
model.random.effects.baseline.single.3cov.re<-function()
{ # *** PROGRAM STARTS 
	# Model for RCTs ##############################################
	for(i in 1:ns){ # LOOP THROUGH STUDIES 
		w[i,1] <- 0 # adjustment for multi-arm trials is zero for control arm 
		delta[i,1] <- 0 # treatment effect is zero for control arm 
		mu[i] ~ dnorm(0,0.0001) # Baseline are a nuisance for RCT evidence
		for (k in 1:na[i]) { # LOOP THROUGH ARMS 
			r[i,k] ~ dbin(p[i,k],n[i,k]) # binomial likelihood 
			logit(p[i,k]) <- mu[i] + delta[i,k] # model for linear predictor 
			rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
			dev.NA[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) #Deviance contribution including NAs 
			+ (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k]))) 
			dev[i,k] <- dev.NA[i,k]*(1-equals(n[i,1],1)) # Deviance contribution with correction for NAs 
		} 
		resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
		for (k in 2:na[i]) { # LOOP THROUGH ARMS 
			delta[i,k] ~ dnorm(md[i,k],taud[i,k]) # trial-specific LOR distributions 
			md[i,k] <- d[t[i,k]] - d[t[i,1]] + sw[i,k] # mean of LOR distributions (with multi-arm trial correction) 
			taud[i,k] <- tau *2*(k-1)/k # precision of LOR distributions (with multi-arm trial correction) 
			w[i,k] <- (delta[i,k] - d[t[i,k]] + d[t[i,1]]) # adjustment for multi-arm RCTs 
			sw[i,k] <- sum(w[i,1:k-1])/(k-1) # cumulative adjustment for multi-arm trials 
		} 
	}

	# Model for baseline effects ###################################
	# Adapted from Program 1 of NICE DSU TSD 5
	for (i in 1:ns.base){ # LOOP THROUGH STUDIES
		r.base[i] ~ dbin(p.base[i],n.base[i]) # Likelihood
		# Linear predictor options
		# Fixed intercept (default)
		logit(p.base[i]) <- m + x.base[cov.index[1],i]*beta.base[1] # + x.base[cov.index[2],i]*beta.base[2] #+ x.base[cov.index[3],i]*beta.base[3] # Log-odds of response
		# Random intercept
		#logit(p.base[i]) <- mu.base[i] + x.base[cov.index[1],i]*beta.base[1] # + x.base[cov.index[2],i]*beta.base[2] + x.base[cov.index[3],i]*beta.base[3] # Log-odds of response
		mu.base[i] ~ dnorm(m,tau.m) # Random effects model
	}
	for(i in 1:3){
		beta.base[i]~dnorm(0,0.0001) #  vague prior for covariate effects
	}
	m ~ dnorm(0,.0001) # vague prior for mean
	var.m <- 1/tau.m # between-trial variance
	tau.m <- pow(sd.m,-2) # between-trial precision = (1/between-trial variance)
	sd.m ~ dunif(0,5) # vague prior for between-trial SD

	# Model for single-arm studies #######################################
	# Prevent feedback to baseline model (these two may not be necessary as mu.single is cut below
	m.cut<-cut(m)
	tau.m.cut<-cut(tau.m)
	for(i in 1:3){
		beta.base.cut[i]<-cut(beta.base[i])
	}
	for(i in 1:ns.single)
	{
		mu.single[i]~dnorm(m.cut,tau.m.cut) # Sampled from baseline random effects model
		mu.single.cut[i]<-cut(mu.single[i]) # Prevent feedback to baseline model
		r.single[i]~dbin(p.single[i],n.single[i])
		# Linear predictor options 
		# Fixed intercept (default)
		logit(p.single[i])<-m.cut + delta.single[i] + x.single[cov.index[1],i]*beta.base.cut[1] #+ x.single[cov.index[2],i]*beta.base.cut[2] #+ x.single[cov.index[3],i]*beta.base.cut[3] 
		# Random intercept	
		#logit(p.single[i])<-mu.single.cut[i] + delta.single[i] + x.single[cov.index[1],i]*beta.base.cut[1] #+ x.single[cov.index[2],i]*beta.base.cut[2] #+ x.single[cov.index[3],i]*beta.base.cut[3] 
		delta.single[i]~dnorm(d[t.single[i]],tau) # Treatment effect relative to reference
		rhat.single[i] <- p.single[i] * n.single[i] # expected value of the numerators 
		dev.single[i] <- 2 * (r.single[i] * (log(r.single[i])-log(rhat.single[i])) 
		+ (n.single[i]-r.single[i]) * (log(n.single[i]-r.single[i]) - log(n.single[i]-rhat.single[i]))) #Deviance contribution 
	}

	# Calculate deviance ##############################################
	totresdev.single<-sum(dev.single[]) # Total residual deviance for single-arm studies
	totresdev.rct <- sum(resdev[]) # Total residual deviance for RCTs
	totresdev<-totresdev.single+totresdev.rct #Total Residual Deviance 
	
	# Specify remaining priors
	d[1]<-0 # treatment effect is zero for reference treatment 
	for (k in 2:nt){ d[k] ~ dnorm(0,0.0001) } # vague priors for treatment effects 
	sd ~ dunif(0,5) # vague prior for between-trial SD 
	tau <- pow(sd,-2) # between-trial precision = (1/between-trial variance) 
}


# The random effects on baseline and fixed treatment effect model
# Including single arm studies by baseline prediction (fixed intercept, up to 3 covariates)
model.random.effects.baseline.single.3cov.fe<-function()
{
	# Model for RCTs ###################################
	for(i in 1:ns){ # LOOP THROUGH STUDIES 
		mu[i] ~ dnorm(0,0.0001) # random effect on baselines
		for (k in 1:na[i]) { # LOOP THROUGH ARMS 
			r[i,k] ~ dbin(p[i,k],n[i,k]) # Binomial likelihood 
			logit(p[i,k]) <- mu[i] + delta[i,k]
			delta[i,k]<-d[t[i,k]] - d[t[i,1]] # model for linear predictor 
			rhat[i,k] <- p[i,k] * n[i,k] # expected value of the numerators 
			dev[i,k] <- 2 * (r[i,k] * (log(r[i,k])-log(rhat[i,k])) 
			+ (n[i,k]-r[i,k]) * (log(n[i,k]-r[i,k]) - log(n[i,k]-rhat[i,k]))) #Deviance contribution 
		} 
	resdev[i] <- sum(dev[i,1:na[i]]) # summed residual deviance contribution for this trial 
	} 

	# Model for baseline effects ###################################
	# Adapted from Program 1 of NICE DSU TSD 5
	for (i in 1:ns.base){ # LOOP THROUGH STUDIES
		r.base[i] ~ dbin(p.base[i],n.base[i]) # Likelihood
		# Linear predictor options
		# Fixed intercept (default)
		logit(p.base[i]) <- m + x.base[cov.index[1],i]*beta.base[1]  #+ x.base[cov.index[2],i]*beta.base[2] #+ x.base[cov.index[3],i]*beta.base[3] # Log-odds of response
		# Random intercept
		#logit(p.base[i]) <- mu.base[i] + x.base[cov.index[1],i]*beta.base[1] #+ x.base[cov.index[2],i]*beta.base[2] #+ x.base[cov.index[3],i]*beta.base[3] # Log-odds of response
		mu.base[i] ~ dnorm(m,tau.m) # Random effects model
	}
	for(i in 1:3){
		beta.base[i]~dnorm(0,0.0001) #  vague prior for covariate effects
	}
	m ~ dnorm(0,.0001) # vague prior for mean
	var.m <- 1/tau.m # between-trial variance
	tau.m <- pow(sd.m,-2) # between-trial precision = (1/between-trial variance)
	sd.m ~ dunif(0,5) # vague prior for between-trial SD


	# Model for single-arm studies ###################################
	# Prevent feedback to baseline model (these two may not be necessary as mu.single is cut below
	m.cut<-cut(m)
	tau.m.cut<-cut(tau.m)
	for(i in 1:3){
		beta.base.cut[i]<-cut(beta.base[i])
	}
	for(i in 1:ns.single)
	{
		mu.single[i]~dnorm(m.cut,tau.m.cut) # Sampled from baseline random effects model	
		mu.single.cut[i]<-cut(mu.single[i]) # Prevent feedback to baseline model
		r.single[i]~dbin(p.single[i],n.single[i])
		# Linear predictor options
		# Fixed intercept (default)
		logit(p.single[i])<-m.cut  + delta.single[i] + x.single[cov.index[1],i]*beta.base.cut[1] #+ x.single[cov.index[2],i]*beta.base.cut[2] #+ x.single[cov.index[3],i]*beta.base.cut[3]
		# Random intercept
		#logit(p.single[i])<-mu.single.cut[i]  + delta.single[i] + x.single[cov.index[1],i]*beta.base.cut[1] #+ x.single[cov.index[2],i]*beta.base.cut[2] + x.single[cov.index[3],i]*beta.base.cut[3]
		delta.single[i]<-d[t.single[i]] # Treatment effect relative to reference

		rhat.single[i] <- p.single[i] * n.single[i] # expected value of the numerators 
		dev.single[i] <- 2 * (r.single[i] * (log(r.single[i])-log(rhat.single[i])) 
		+ (n.single[i]-r.single[i]) * (log(n.single[i]-r.single[i]) - log(n.single[i]-rhat.single[i]))) #Deviance contribution 
	}

	# Calculate deviance ###################################
	totresdev.single<-sum(dev.single[]) # Total residual deviance for single-arm studies
	totresdev.rct <- sum(resdev[]) # Total residual deviance for RCTs
	totresdev<-totresdev.single+totresdev.rct #Total Residual Deviance 
	
	# Priors for remaining parameters
	d[1]<-0 # treatment effect is zero for reference treatment 
	for (k in 2:nt){ d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects 
}



