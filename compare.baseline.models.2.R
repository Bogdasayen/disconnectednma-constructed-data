# Script to perform outputs for baseline model selection
# Howard Thom 21-October-2017
# Relies on data generated by test.re.nma.3.R, so is run in that script, 
# just before the NMAs are run

# Version 2 produces various plots

source("baseline.models.2.R")

bugs.inits.base<-bugs.inits.base.re<-list()
bugs.baseline.fe.nocov<-bugs.baseline.fe.age<-bugs.baseline.fe.male<-bugs.baseline.fe.chads2<-bugs.baseline.fe.age.male<-bugs.baseline.fe.age.chads2<-bugs.baseline.fe.male.chads2<-bugs.baseline.fe.age.male.chads2<-list()
bugs.baseline.re.age.male<-bugs.baseline.re.age.chads2<-bugs.baseline.re.male.chads2<-bugs.baseline.re.age.male.chads2<-list()
bugs.baseline.re.nocov<-bugs.baseline.re.age<-bugs.baseline.re.male<-bugs.baseline.re.chads2<-list()

n.chains<-2			 # 2
num.sims=30000*n.chains  # 150000*n.chains
burn.in=30000*n.chains	 # 100000*n.chains


# Inverse of logistic link function
expit<-function(x)
{
	return(exp(x)/(1+exp(x)))
}
# Logistic link function
logit<-function(x)
{
	return(log(x/(1-x)))
}


# Look at correlation between covariate and event probability	
for(outcome.i in 1:length(outcome.names))
{
	# Age and CHADS2 look correlated, but gender seems independent.
	jpeg(file=paste(baseline.directory,"/results/correlation plots/",outcome.names[outcome.i],".vs.",dimnames(bugs.data.single[[outcome.i]]$x)[[1]][1],".jpg",sep=""))
	plot(x=bugs.data.single[[outcome.i]]$x.base[1,],y=logit(bugs.data.single[[outcome.i]]$r.base/bugs.data.single[[outcome.i]]$n.base),xlab=dimnames(bugs.data.single[[outcome.i]]$x)[[1]][1],main=paste(outcome.names[outcome.i],"vs age"),ylab=paste("Log odds of",outcome.names[outcome.i]))
	dev.off()
	jpeg(file=paste(baseline.directory,"/results/correlation plots/",outcome.names[outcome.i],".vs.",dimnames(bugs.data.single[[outcome.i]]$x)[[1]][2],".jpg",sep=""))
	plot(x=bugs.data.single[[outcome.i]]$x.base[2,],y=logit(bugs.data.single[[outcome.i]]$r.base/bugs.data.single[[outcome.i]]$n.base),xlab=dimnames(bugs.data.single[[outcome.i]]$x)[[1]][2],main=paste(outcome.names[outcome.i],"vs proportion male"),ylab=paste("Log odds of",outcome.names[outcome.i]))
	dev.off()
	jpeg(file=paste(baseline.directory,"/results/correlation plots/",outcome.names[outcome.i],".vs.",dimnames(bugs.data.single[[outcome.i]]$x)[[1]][3],".jpg",sep=""))
	plot(x=bugs.data.single[[outcome.i]]$x.base[3,],y=logit(bugs.data.single[[outcome.i]]$r.base/bugs.data.single[[outcome.i]]$n.base),xlab=dimnames(bugs.data.single[[outcome.i]]$x)[[1]][3],main=paste(outcome.names[outcome.i],"vs CHADS2 stroke risk"),ylab=paste("Log odds of",outcome.names[outcome.i]))
	dev.off()
}

# Conduct baseline meta-regressions
for(outcome.i in 1:length(outcome.names))
{
	# Initial values for baseline models (use same set for all models, so specify some parameters that are not in the model)
	inits1<-list(m=0.1,sd.m=1,mu.base=rep(0.1,bugs.data.single[[outcome.i]]$ns.base),beta.base=rep(0.5,3))
	inits2<-list(m=0.1,sd.m=1,mu.base=rep(0.1,bugs.data.single[[outcome.i]]$ns.base),beta.base=rep(0.5,3))
	inits1.re<-list(m=0.1,sd.m=1,mu.base=rep(0.1,bugs.data.single[[outcome.i]]$ns.base),sd=1,beta.base=rep(0.5,3))
	inits2.re<-list(m=0.5,sd.m=2,mu.base=rep(0.5,bugs.data.single[[outcome.i]]$ns.base),sd=0.5,beta.base=rep(-0.5,3))
	bugs.inits.base[[outcome.i]]<-list(inits1,inits2)		
	bugs.inits.base.re[[outcome.i]]<-list(inits1.re,inits2.re)
	#bugs.inits.single.3cov.re[[outcome.i]]<-list(inits1.re,inits2.re)

	# Fit the various baseline models.
	# No covariates
	bugs.baseline.fe.nocov[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","totresdev"),model=model.baseline.fe.nocov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	bugs.baseline.re.nocov[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base.re[[outcome.i]],parameters.to.save=c("m","sd.m","mu.new","totresdev"),model=model.baseline.re.nocov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)

	# With age as covariate
	bugs.data.single[[outcome.i]]$cov.index[1]<-1
	bugs.baseline.fe.age[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","beta.base","mu.new","totresdev","dev"),model=model.baseline.fe.1cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	bugs.baseline.re.age[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base.re[[outcome.i]],parameters.to.save=c("m","sd.m","beta.base","totresdev"),model=model.baseline.re.1cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=TRUE)


	# With prop male as covariate
	bugs.data.single[[outcome.i]]$cov.index[1]<-2
	bugs.baseline.fe.male[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","beta.base","mu.new","totresdev"),model=model.baseline.fe.1cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	bugs.baseline.re.male[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base.re[[outcome.i]],parameters.to.save=c("m","sd.m","beta.base","totresdev"),model=model.baseline.re.1cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=TRUE)
	# With CHADS2 as covariate
	bugs.data.single[[outcome.i]]$cov.index[1]<-3
	bugs.baseline.fe.chads2[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","beta.base","mu.new","totresdev"),model=model.baseline.fe.1cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	bugs.baseline.re.chads2[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","sd.m","beta.base","totresdev"),model=model.baseline.re.1cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)

	# With age and prop male 
	bugs.data.single[[outcome.i]]$cov.index[1:2]<-c(1,2)
	bugs.baseline.fe.age.male[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","beta.base","totresdev"),model=model.baseline.fe.2cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	bugs.baseline.re.age.male[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","sd.m","beta.base","totresdev"),model=model.baseline.re.2cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	# With age and chads2 
	bugs.data.single[[outcome.i]]$cov.index[1:2]<-c(1,3)
	bugs.baseline.fe.age.chads2[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","beta.base","totresdev"),model=model.baseline.fe.2cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	bugs.baseline.re.age.chads2[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","sd.m","beta.base","totresdev"),model=model.baseline.re.2cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	# With prop male and chads2 
	bugs.data.single[[outcome.i]]$cov.index[1:2]<-c(2,3)
	bugs.baseline.fe.male.chads2[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","beta.base","totresdev"),model=model.baseline.fe.2cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	bugs.baseline.re.male.chads2[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","sd.m","beta.base","totresdev"),model=model.baseline.re.2cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	# With age, prop male, and chads2
	bugs.data.single[[outcome.i]]$cov.index[1:3]<-c(1,2,3)
	bugs.baseline.fe.age.male.chads2[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","beta.base","totresdev"),model=model.baseline.fe.3cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)
	bugs.baseline.re.age.male.chads2[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.base[[outcome.i]],parameters.to.save=c("m","sd.m","beta.base","totresdev"),model=model.baseline.re.3cov,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,debug=FALSE)


	# Reset the covariate indices as they are used elsewhere
	bugs.data.single[[outcome.i]]$cov.index[1:3]<-c(1:3)	

	# Build a table comparing the residual deviance, DIC and estimate of m.base
	model.comparison<-matrix(nrow=16,ncol=6)
	rownames(model.comparison)<-c("Random effects","Fixed effects","Fixed effects with age","Fixed effects with prop male","Fixed effects with CHADS2",
			"Fixed effects with age and prop male","Fixed effects with age and CHADS2","Fixed effects with prop male and CHADS2","Fixed effects with all 3",
			"Random effects with age","Random effects with prop male","Random effects with CHADS2",
			"Random effects with age and prop male","Random effects with age and CHADS2","Random effects with prop male and CHADS2","Random effects with all 3")
	colnames(model.comparison)<-c("Residual deviance","DIC","m","sd.m","beta.base","Number data points")
	model.comparison[,"Number data points"]<-bugs.data.single[[outcome.i]]$ns.base
	model.comparison[1,c("Residual deviance","DIC","m","sd.m")]<-c(format.results(bugs.baseline.re.nocov[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.re.nocov[[outcome.i]]$DIC,format.results(bugs.baseline.re.nocov[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.re.nocov[[outcome.i]]$summary["sd.m",c(1,3,7)]))
	model.comparison[2,c("Residual deviance","DIC","m")]<-c(format.results(bugs.baseline.fe.nocov[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.fe.nocov[[outcome.i]]$DIC,format.results(bugs.baseline.fe.nocov[[outcome.i]]$summary["m",c(1,3,7)]))
	model.comparison[3,c("Residual deviance","DIC","m","beta.base")]<-c(format.results(bugs.baseline.fe.age[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.fe.age[[outcome.i]]$DIC,format.results(bugs.baseline.fe.age[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.fe.age[[outcome.i]]$summary["beta.base[1]",c(1,3,7)]))
	model.comparison[4,c("Residual deviance","DIC","m","beta.base")]<-c(format.results(bugs.baseline.fe.male[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.fe.male[[outcome.i]]$DIC,format.results(bugs.baseline.fe.male[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.fe.male[[outcome.i]]$summary["beta.base[1]",c(1,3,7)]))
	model.comparison[5,c("Residual deviance","DIC","m","beta.base")]<-c(format.results(bugs.baseline.fe.chads2[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.fe.chads2[[outcome.i]]$DIC,format.results(bugs.baseline.fe.chads2[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.fe.chads2[[outcome.i]]$summary["beta.base[1]",c(1,3,7)]))

	# Only report model fit and m for multiple covariate models (they don't converge)
	model.comparison[6,c("Residual deviance","DIC","m")]<-c(format.results(bugs.baseline.fe.age.male[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.fe.age.male[[outcome.i]]$DIC,format.results(bugs.baseline.fe.age.male[[outcome.i]]$summary["m",c(1,3,7)]))
	model.comparison[7,c("Residual deviance","DIC","m")]<-c(format.results(bugs.baseline.fe.age.chads2[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.fe.age.chads2[[outcome.i]]$DIC,format.results(bugs.baseline.fe.age.chads2[[outcome.i]]$summary["m",c(1,3,7)]))
	model.comparison[8,c("Residual deviance","DIC","m")]<-c(format.results(bugs.baseline.fe.male.chads2[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.fe.male.chads2[[outcome.i]]$DIC,format.results(bugs.baseline.fe.male.chads2[[outcome.i]]$summary["m",c(1,3,7)]))
	model.comparison[9,c("Residual deviance","DIC","m")]<-c(format.results(bugs.baseline.fe.age.male.chads2[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.fe.age.male.chads2[[outcome.i]]$DIC,format.results(bugs.baseline.fe.age.male.chads2[[outcome.i]]$summary["m",c(1,3,7)]))


	model.comparison[10,c("Residual deviance","DIC","m","sd.m","beta.base")]<-c(format.results(bugs.baseline.re.age[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.re.age[[outcome.i]]$DIC,format.results(bugs.baseline.re.age[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.re.age[[outcome.i]]$summary["sd.m",c(1,3,7)]),format.results(bugs.baseline.re.age[[outcome.i]]$summary["beta.base[1]",c(1,3,7)]))
	model.comparison[11,c("Residual deviance","DIC","m","sd.m","beta.base")]<-c(format.results(bugs.baseline.re.male[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.re.male[[outcome.i]]$DIC,format.results(bugs.baseline.re.male[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.re.male[[outcome.i]]$summary["sd.m",c(1,3,7)]),format.results(bugs.baseline.re.male[[outcome.i]]$summary["beta.base[1]",c(1,3,7)]))
	model.comparison[12,c("Residual deviance","DIC","m","sd.m","beta.base")]<-c(format.results(bugs.baseline.re.chads2[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.re.chads2[[outcome.i]]$DIC,format.results(bugs.baseline.re.chads2[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.re.chads2[[outcome.i]]$summary["sd.m",c(1,3,7)]),format.results(bugs.baseline.re.chads2[[outcome.i]]$summary["beta.base[1]",c(1,3,7)]))

	# Only report model fit and m for multiple covariate models (they don't converge)
	model.comparison[13,c("Residual deviance","DIC","m","sd.m")]<-c(format.results(bugs.baseline.re.age.male[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.re.age.male[[outcome.i]]$DIC,format.results(bugs.baseline.re.age.male[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.re.age.male[[outcome.i]]$summary["sd.m",c(1,3,7)]))
	model.comparison[14,c("Residual deviance","DIC","m","sd.m")]<-c(format.results(bugs.baseline.re.age.chads2[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.re.age.chads2[[outcome.i]]$DIC,format.results(bugs.baseline.re.age.chads2[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.re.age.chads2[[outcome.i]]$summary["sd.m",c(1,3,7)]))
	model.comparison[15,c("Residual deviance","DIC","m","sd.m")]<-c(format.results(bugs.baseline.re.male.chads2[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.re.male.chads2[[outcome.i]]$DIC,format.results(bugs.baseline.re.male.chads2[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.re.male.chads2[[outcome.i]]$summary["sd.m",c(1,3,7)]))
	model.comparison[16,c("Residual deviance","DIC","m","sd.m")]<-c(format.results(bugs.baseline.re.age.male.chads2[[outcome.i]]$summary["totresdev",c(1,3,7)]),bugs.baseline.re.age.male.chads2[[outcome.i]]$DIC,format.results(bugs.baseline.re.age.male.chads2[[outcome.i]]$summary["m",c(1,3,7)]),format.results(bugs.baseline.re.age.male.chads2[[outcome.i]]$summary["sd.m",c(1,3,7)]))

	# Export the model comparison results for this outcome
	write.csv(model.comparison,file=paste(baseline.directory,"/results/baseline models/model.comparison.",outcome.names[outcome.i],".2.csv",sep=""))

}

# Save the bugs objects in case they're needed
#save(file=paste(baseline.directory,"/results/baseline models/bugs objects/bugs.baseline.rda",sep=""),
	bugs.baseline.re.nocov,bugs.baseline.fe.age,bugs.baseline.fe.male,bugs.baseline.fe.chads2,
	bugs.baseline.fe.age.male,bugs.baseline.fe.male.chads2,bugs.baseline.fe.age.chads2,bugs.baseline.fe.age.male.chads2,
	bugs.baseline.re.age,bugs.baseline.re.male,bugs.baseline.re.chads2,
	bugs.baseline.re.age.male,bugs.baseline.re.male.chads2,bugs.baseline.re.age.chads2,bugs.baseline.re.age.male.chads2)
load(file=paste(baseline.directory,"/results/baseline models/bugs objects/bugs.baseline.rda",sep=""))

# Plot the baseline response against covariates
# Include predictions
# Look at correlation between covariate and event probability	
for(outcome.i in 1:length(outcome.names))
{
	cov.names<-c("Age","Male","CHADS2")
	for(i.cov in 1:3)
	{
		jpeg(file=paste(baseline.directory,"/results/baseline models/regression comparison/",outcome.names[outcome.i]," vs ",cov.names[i.cov],".jpeg",sep=""))
		
		# Define the limits of the plot
		x.min=min(bugs.data.single[[outcome.i]]$x.base[i.cov,],bugs.data.single[[outcome.i]]$x.single[i.cov,])
		x.max=max(bugs.data.single[[outcome.i]]$x.base[i.cov,],bugs.data.single[[outcome.i]]$x.single[i.cov,])
		y.min=-1+min(bugs.baseline.re.nocov[[outcome.i]]$summary["mu.new",3],logit(bugs.data.single[[outcome.i]]$r.base/bugs.data.single[[outcome.i]]$n.base),logit(bugs.data.single[[outcome.i]]$r.single/bugs.data.single[[outcome.i]]$n.single)[bugs.data.single[[outcome.i]]$t.single==1])
		y.max=max(bugs.baseline.re.nocov[[outcome.i]]$summary["mu.new",7],logit(bugs.data.single[[outcome.i]]$r.base/bugs.data.single[[outcome.i]]$n.base),logit(bugs.data.single[[outcome.i]]$r.single/bugs.data.single[[outcome.i]]$n.single)[bugs.data.single[[outcome.i]]$t.single==1])

		# Plot the points observed
		plot(c(0,0),col=0,xlim=c(x.min,x.max),ylim=c(y.min,y.max),xlab=dimnames(bugs.data.single[[outcome.i]]$x)[[1]][i.cov],ylab=paste("Log odds of",outcome.names[outcome.i]))
		title(paste(outcome.names[outcome.i],"vs",cov.names[i.cov]))
		points(x=bugs.data.single[[outcome.i]]$x.base[i.cov,],y=logit(bugs.data.single[[outcome.i]]$r.base/bugs.data.single[[outcome.i]]$n.base),pch=15,cex=4*sqrt(bugs.data.single[[outcome.i]]$n.base)/max(sqrt(bugs.data.single[[outcome.i]]$n.base)))

		# Plot observed and predicted response from single arm studies not used in fitting baseline model
		points(x=bugs.data.single[[outcome.i]]$x.single[i.cov,bugs.data.single[[outcome.i]]$t.single==1],y=logit(bugs.data.single[[outcome.i]]$r.single/bugs.data.single[[outcome.i]]$n.single)[bugs.data.single[[outcome.i]]$t.single==1],pch=15,
			cex=4*sqrt(bugs.data.single[[outcome.i]]$n.single[bugs.data.single[[outcome.i]]$t.single==1])/max(sqrt(bugs.data.single[[outcome.i]]$n.base)),col=2)
		# Add a star as the studies are too small for boxes to be visible
		points(x=bugs.data.single[[outcome.i]]$x.single[i.cov,bugs.data.single[[outcome.i]]$t.single==1],y=logit(bugs.data.single[[outcome.i]]$r.single/bugs.data.single[[outcome.i]]$n.single)[bugs.data.single[[outcome.i]]$t.single==1],pch=4,col=2)

		# Now plot the predictions
		if(i.cov==1){bugs.summary<-bugs.baseline.fe.age[[outcome.i]]$summary}
		if(i.cov==2){bugs.summary<-bugs.baseline.fe.male[[outcome.i]]$summary}
		if(i.cov==3){bugs.summary<-bugs.baseline.fe.chads2[[outcome.i]]$summary}
		r.names<-rownames(bugs.summary)
		for(i.single in which(bugs.data.single[[outcome.i]]$t.single==1))
		{
			# Plot the posterior predictive means
			points(x=bugs.data.single[[outcome.i]]$x.single[i.cov,i.single],bugs.summary[grep("mu.new",r.names),1][i.single])
			# Plot the credible intervals
			lines(x=rep(bugs.data.single[[outcome.i]]$x.single[i.cov,i.single],2),y=bugs.summary[grep("mu.new",r.names),c(3,7)][i.single,])
		}

		
		# Plot the random effects on baseline distributions (same for all covariates)
		lines(x=c(x.min,x.max),y=rep(bugs.baseline.re.nocov[[outcome.i]]$summary["mu.new",1],2))
		lines(x=c(x.min,x.max),y=rep(bugs.baseline.re.nocov[[outcome.i]]$summary["m",3],2),lty=2)
		lines(x=c(x.min,x.max),y=rep(bugs.baseline.re.nocov[[outcome.i]]$summary["m",7],2),lty=2)
		lines(x=c(x.min,x.max),y=rep(bugs.baseline.re.nocov[[outcome.i]]$summary["mu.new",3],2),lty=3)
		lines(x=c(x.min,x.max),y=rep(bugs.baseline.re.nocov[[outcome.i]]$summary["mu.new",7],2),lty=3)

		legend("bottomright",legend=c("Observed","New (single arms)","Predicted"),pch=c(15,4,1),col=c(1,2,1))
	
		#polygon(x=c(x.min,x.max,x.min,x.max),y=c(rep(bugs.baseline.re.nocov[[outcome.i]]$summary["mu.new",7],2),rep(bugs.baseline.re.nocov[[outcome.i]]$summary["mu.new",3],2)))
		dev.off()
	} # End loop over covariates
}

# Plot baseline response against predictions
# Only need include the chosen model (RE on baseline, no covariates)
# Include r.base/n.base and r.single/n.single (for warfarin arms) and distributions of m and mu.new
for(outcome.i in 1:n.outcomes)
{
	mu.posterior<-bugs.baseline.re.nocov[[outcome.i]]$summary["m",c(1,3,7)]
	mu.posterior.predictive<-bugs.baseline.re.nocov[[outcome.i]]$summary["mu.new",c(1,3,7)]
	# Mean observed log odds of event
	mean.base<-logit(bugs.data.single[[outcome.i]]$r.base/bugs.data.single[[outcome.i]]$n.base)
	# Mean observed log odds of event for new data points (those in single arm trials)
	mean.single<-logit(bugs.data.single[[outcome.i]]$r.single/bugs.data.single[[outcome.i]]$n.single)[bugs.data.single[[outcome.i]]$t.single==1]

	n.lines<-2+length(mean.base)+length(mean.single)
	x.min<-min(c(mean.base,mean.single,mu.posterior.predictive))-3
	x.max<-max(c(mean.base,mean.single,mu.posterior.predictive))+1

	# Export a plot comparing observed, new, and predicted baselines
	jpeg(file=paste(baseline.directory,"/results/baseline models/forest plots/forest.",outcome.names[[outcome.i]],".jpeg",sep=""))
	plot(c(0,0),xlim=c(x.min,x.max),ylim=c(0,n.lines),axes=FALSE,xlab="Log odds",ylab="")
	axis(side=1)
	title(paste("Observed vs predicted log odds of ",outcome.names[[outcome.i]],sep=""))
	text("Posterior predictive",x=x.min,y=0,pos=4)
	points(x=mu.posterior.predictive[1],y=0)
	lines(x=mu.posterior.predictive[c(2,3)],y=c(0,0))
	text("Posterior",x=x.min,y=1,pos=4)
	points(x=mu.posterior[1],y=1)
	lines(x=mu.posterior[c(2,3)],y=c(1,1))
	lines(x=rep(mu.posterior[1],2),y=c(-1,n.lines))
	text("New baseline(s)",x=x.min,pos=4,y=mean(c(1:length(mean.single))+1))
	points(x=mean.single,y=c(1:length(mean.single))+1,pch=4) 		#,cex=4*sqrt(bugs.data.single[[outcome.i]]$n.single)/max(sqrt(bugs.data.single[[outcome.i]]$n.base)))
	text("Observed baselines",x=x.min,pos=4,y=mean(c((1+length(mean.single)):(length(mean.single)+length(mean.base)))+1))
	points(x=mean.base,y=c((1+length(mean.single)):(length(mean.single)+length(mean.base)))+1,pch=15,cex=4*sqrt(bugs.data.single[[outcome.i]]$n.base)/max(sqrt(bugs.data.single[[outcome.i]]$n.base)))
	dev.off()	
}
