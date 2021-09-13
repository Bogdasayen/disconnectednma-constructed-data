# Script to test code for random effects on baseline and
# inclusion of single-arm studies in an NMA
# Version 3 updated to apply to the NOACs AF data
# Howard Thom 25-September-2017

# V4 updated to use plug-in estimator for ALM

# V5
# Use the standard deviation from connected RCTs as informative prior for disconnected RCTs and single-arms
# Centre the covariate NMA models (baseline models already fixed) and put random effect
# around intercept with covariates to improve fit and stability


library(R2OpenBUGS)

# Function to draw the network plot
source("mtm.networkplot.fun.R")
# The following are the NICE TSD2 1c (RE) and 1d (FE) functions
source("independent.baselines.model.R")
# The following are modifications of the NICE TSD2 1c (RE) and 1d (FE) functions
source("random.effects.baseline.model.3.R")


# The following uses random effects on baseline to include single arm studies and keeps RCT network otherwise separate from single-arm studies
source("random.effects.baseline.model.single-arms.3.R")
# The following use fixed effects with (up to) 3 covariates on baseline to include single arm studies and keeps RCT network otherwise separate from single-arm studies
source("random.effects.baseline.model.single-arms.3cov.1.R")

# The following uses random effects on baseline to include disconnected RCTs but keeps RCT network (with reference) otherwise separate
source("random.effects.baseline.model.disconnected.5.R")
# The following uses fixed effects with (up to) 3 covariates on baseline to include disconnected RCTs but keeps RCT network (with reference) otherwise separate
source("random.effects.baseline.model.disonnected.3cov.2.R")

# Load Joy Leahy's code
source("matching.model.hierarchical.1.R")

# Plugin estimator models
source("plugin.model.single.3.R")
source("plugin.model.disconnected.3.R")

load("af.nma.cov.rda")
source("clean.af.data.1.R")

# Data on atrial fibrillation
cov.names<-dimnames(bugs.data[[1]]$x)[[1]][1]
outcome.names<-names(bugs.data)
n.outcomes<-length(outcome.names)
#study.names<-rownames(bugs.data[[1]]$r) # This depends on outcome so don't use

# Define the covariate indices for use in the models
for(outcome.i in 1:length(outcome.names))
{
	bugs.data[[outcome.i]]$cov.index<-c(1,2,3)
}

# Export a simple data summary table for publication
for(outcome.i in c(1,4)) {
  treatment.matrix.temp <- matrix("-", nrow = bugs.data[[outcome.i]]$ns, ncol = max(bugs.data[[outcome.i]]$na))
  events.matrix.temp <- matrix("-", nrow = bugs.data[[outcome.i]]$ns, ncol = max(bugs.data[[outcome.i]]$na))
  rownames(treatment.matrix.temp) <- rownames(events.matrix.temp) <- rownames(bugs.data[[outcome.i]]$r)
  colnames(treatment.matrix.temp) <- colnames(events.matrix.temp) <- paste("Arm", 1:max(bugs.data[[outcome.i]]$na))
  for(study.i in 1:bugs.data[[outcome.i]]$ns) {
    treatment.matrix.temp[study.i, 1:bugs.data[[outcome.i]]$na[study.i]] <- treatment.names[[outcome.i]][bugs.data[[outcome.i]]$t[study.i, 1:bugs.data[[outcome.i]]$na[study.i]]]
    events.matrix.temp[study.i, 1:bugs.data[[outcome.i]]$na[study.i]] <- paste0(bugs.data[[outcome.i]]$r[study.i, ], "/", bugs.data[[outcome.i]]$n[study.i, ])[1:bugs.data[[outcome.i]]$na[study.i]]
  }
  write.csv(events.matrix.temp, file = paste0(baseline.directory, "/results/data summaries/events.", outcome.names[outcome.i], ".csv"))
  write.csv(treatment.matrix.temp, file = paste0(baseline.directory, "/results/data summaries/treatments.", outcome.names[outcome.i], ".csv"))
}



# Set up the baseline data
for(outcome.i in 1:length(outcome.names))
{
	# Now extract arms from single-arm studies or RCTs that are on reference (tr=1)
	tr<-bugs.data[[outcome.i]]$t
	ns.base<-sum(tr==1,na.rm=TRUE)
	r.base<-bugs.data[[outcome.i]]$r[tr==1][1:ns.base]
	n.base<-bugs.data[[outcome.i]]$n[tr==1][1:ns.base]
	x.base<-matrix(NA,ncol=ns.base,nrow=3)
	for(cov.i in 1:3){
		x.base[cov.i,]<-bugs.data[[outcome.i]]$x[cov.i,,][tr==1][1:ns.base]

		# Set the covariate to the trial mean if not reported (code a bit convoluted)
		trial.means<-rowMeans(bugs.data[[outcome.i]]$x[cov.i,,],na.rm=TRUE)
		x.means<-bugs.data[[outcome.i]]$x[cov.i,,]
		for(trial.i in 1:bugs.data[[outcome.i]]$ns)
		{
			x.means[trial.i,1:bugs.data[[outcome.i]]$na[trial.i]]<-rep(trial.means[trial.i],bugs.data[[outcome.i]]$na[trial.i])
		}
		# First set to trial means
		for(i.baseline.arm in 1:ns.base)
		{
			if(is.na(x.base[cov.i,i.baseline.arm])){
				x.base[cov.i,i.baseline.arm]<-x.means[tr==1][1:ns.base][i.baseline.arm]
			}
		}
	}

	
	# Look at correlation between covariate and event probability	
	# Age and CHADS2 look correlated, but gender seems independent.
	jpeg(file=paste(baseline.directory,"/results/correlation plots/",outcome.names[outcome.i],".vs.",dimnames(bugs.data[[outcome.i]]$x)[[1]][1],".jpg",sep=""))
	#plot(x=bugs.data[[outcome.i]]$x[1,,],y=bugs.data[[outcome.i]]$r/bugs.data[[outcome.i]]$n,ylab=names(bugs.data)[outcome.i],xlab=dimnames(bugs.data[[outcome.i]]$x)[[1]][1])
	plot(x=x.base[1,],y=r.base/n.base,ylab=names(bugs.data)[outcome.i],xlab=dimnames(bugs.data[[outcome.i]]$x)[[1]][1],main="Correlation reference arms")
	dev.off()
	jpeg(file=paste(baseline.directory,"/results/correlation plots - all data/",outcome.names[outcome.i],".vs.",dimnames(bugs.data[[outcome.i]]$x)[[1]][2],".jpg",sep=""))
	plot(x=x.base[2,],y=r.base/n.base,ylab=names(bugs.data)[outcome.i],xlab=dimnames(bugs.data[[outcome.i]]$x)[[1]][2],main="Correlation reference arms")
	dev.off()
	jpeg(file=paste(baseline.directory,"/results/correlation plots/",outcome.names[outcome.i],".vs.",dimnames(bugs.data[[outcome.i]]$x)[[1]][3],".jpg",sep=""))
	plot(x=x.base[3,],y=r.base/n.base,ylab=names(bugs.data)[outcome.i],xlab=dimnames(bugs.data[[outcome.i]]$x)[[1]][3],main="Correlation reference arms")
	dev.off()

	# Need treatment matrix in special vector format
	tr.temp<-tr<-bugs.data[[outcome.i]]$t
	for(i in 1:dim(tr.temp)[1])
	{
		tr.temp[i,]<-i
	}
	t1<-c(tr.temp[!is.na(tr)])
	t2<-c(tr[!is.na(tr)])

	# These network plots don't make sense in terms of treatment ordering
	jpeg(file=paste(baseline.directory,"/results/network plots/",outcome.names[outcome.i],".jpg",sep=""),width=960,height=480)
	mtm.networkplot.fun(t1=t1,t2=t2,percomparison=FALSE,nameoftreatments=paste(1:bugs.data[[outcome.i]]$nt,treatment.names[[outcome.i]]))
	title(outcome.names[outcome.i])
	dev.off()
}

# Construct two alternative data scenarios (both disconnected but one uses single-arm studies, the other RCTs)
# 1. "single": Convert dabigatran 110mg and 150mg arms into single-arm studies 
# Leave the warfarin arms from RE-LY and PETRO as single-arm studies; comparing fit of predicted mu to them helps assess the model
# 2. "disconnect": Remove common comparator (warfarin) arms from the dabigatran 110mg and 150mg RCTs giving a disconnected network of
# RE-LY (dab 110 vs dab 150), PETRO (9 dabig doses with/without aspirin), AF-DABIG-VKA-JAPAN (no change, dab 110 vs dab 150).
# PETRO and AF-DABIG-VKA-JAPAN are only in the Clinically relevant bleeding network.

bugs.data.single<-list()
bugs.data.disconnect<-list()

for(outcome.i in 1:n.outcomes)
{
	bugs.data.single[[outcome.i]]<-bugs.data.disconnect[[outcome.i]]<-bugs.data[[outcome.i]]

	tr<-bugs.data[[outcome.i]]$t
	# Which are the dabigatran treamtent codes?
	remove.treat<-which(treatment.names[[outcome.i]]=="Dabigatran (110mg bd)" | treatment.names[[outcome.i]]=="Dabigatran (150mg bd)")
	# A boolean vector indicating whether or not to remove the trial (if it contains dabigatran
	remove.trial<-rowSums(bugs.data[[outcome.i]]$t==remove.treat[1] | bugs.data[[outcome.i]]$t==remove.treat[2],na.rm=TRUE)>0

	# Extract the (unchanged) non-dabigatran RCTs
	bugs.data.single[[outcome.i]]$t<-bugs.data.disconnect[[outcome.i]]$t<-bugs.data[[outcome.i]]$t[!remove.trial,]	
	bugs.data.single[[outcome.i]]$r<-bugs.data.disconnect[[outcome.i]]$r<-bugs.data[[outcome.i]]$r[!remove.trial,]	
	bugs.data.single[[outcome.i]]$n<-bugs.data.disconnect[[outcome.i]]$n<-bugs.data[[outcome.i]]$n[!remove.trial,]	
	bugs.data.single[[outcome.i]]$x<-bugs.data.disconnect[[outcome.i]]$x<-bugs.data[[outcome.i]]$x[,!remove.trial,]	
	bugs.data.single[[outcome.i]]$na<-bugs.data.disconnect[[outcome.i]]$na<-bugs.data[[outcome.i]]$na[!remove.trial]	
	bugs.data.single[[outcome.i]]$ns<-bugs.data.disconnect[[outcome.i]]$ns<-bugs.data[[outcome.i]]$ns-sum(remove.trial)
	# nt doesn't change as it is a characteristic of all included treatments

	# Build single arm studies
	# Remove treatment=NA arms
	t.single<-as.vector(bugs.data[[outcome.i]]$t[remove.trial,])[!is.na(as.vector(bugs.data[[outcome.i]]$t[remove.trial,]))]
	ns.single<-length(t.single)
	x.single<-matrix(NA,nrow=3,ncol=ns.single)
	rownames(x.single)<-dimnames(bugs.data[[outcome.i]]$x)[[1]]
	r.single<-as.vector(bugs.data[[outcome.i]]$r[remove.trial,])[!is.na(as.vector(bugs.data[[outcome.i]]$t[remove.trial,]))]
	n.single<-as.vector(bugs.data[[outcome.i]]$n[remove.trial,])[!is.na(as.vector(bugs.data[[outcome.i]]$t[remove.trial,]))]
	for(i.cov in 1:dim(x.single)[1]){
		x.single[i.cov,]<-as.vector(bugs.data[[outcome.i]]$x[i.cov,remove.trial,])[!is.na(as.vector(bugs.data[[outcome.i]]$t[remove.trial,]))]
		# Set the covarite to the mean if not observed
		x.single[i.cov,is.na(x.single[i.cov,])]<-mean(x.single[i.cov,],na.rm=TRUE)
	}	


	# Build the disconnected RCTs
	t.disc<-bugs.data[[outcome.i]]$t[remove.trial,]
	ns.disc<-sum(remove.trial)
	r.disc<-bugs.data[[outcome.i]]$r[remove.trial,]
	n.disc<-bugs.data[[outcome.i]]$n[remove.trial,]
	na.disc<-bugs.data[[outcome.i]]$na[remove.trial]
	x.disc<-bugs.data[[outcome.i]]$x[,remove.trial,]
	if(ns.disc==1)
	{
		r.disc<-matrix(r.disc,nrow=1,ncol=length(r.disc))
		n.disc<-matrix(n.disc,nrow=1,ncol=length(n.disc))
		t.disc<-matrix(t.disc,nrow=1,ncol=length(t.disc))
		x.disc<-array(data=x.disc,dim=c(dim(x.disc)[1],1,dim(x.disc)[2]))
	}
	# Remove the reference (warfarin) t=1 from the RCTs
	for(i.trial in 1:ns.disc)
	{
		reference.arm<-t.disc[i.trial,]==1
		na.disc[i.trial]<-na.disc[i.trial]-sum(reference.arm,na.rm=TRUE)	
		r.disc[i.trial,]<-c(r.disc[i.trial,!reference.arm][1:na.disc[i.trial]],rep(NA,dim(bugs.data[[outcome.i]]$r)[2]-na.disc[i.trial]))
		n.disc[i.trial,]<-c(n.disc[i.trial,!reference.arm][1:na.disc[i.trial]],rep(NA,dim(bugs.data[[outcome.i]]$r)[2]-na.disc[i.trial]))
		t.disc[i.trial,]<-c(t.disc[i.trial,!reference.arm][1:na.disc[i.trial]],rep(NA,dim(bugs.data[[outcome.i]]$r)[2]-na.disc[i.trial]))
		for(i.cov in 1:dim(x.disc)[1]){
			x.disc[i.cov,i.trial,]<-c(x.disc[i.cov,i.trial,!reference.arm][1:na.disc[i.trial]],rep(NA,dim(bugs.data[[outcome.i]]$r)[2]-na.disc[i.trial]))
			# Set the baseline covariate to the mean if it is not reported
			x.disc[i.cov,i.trial,is.na(x.disc[i.cov,i.trial,])]<-mean(x.disc[i.cov,,],na.rm=TRUE)
		}
	}
	# Avoid the "expected collation operator" error by including an unused element to make sure na.disc is a vectorand r.disc, n.disc, t.disc matrices
	if(ns.disc==1)
	{
		na.disc<-c(na.disc,NA)
		r.disc<-rbind(r.disc,rep(NA,length(r.disc)))
		n.disc<-rbind(n.disc,rep(NA,length(n.disc)))
		t.disc<-rbind(t.disc,rep(NA,length(t.disc)))
		x.disc.new<-array(NA,dim=c(dim(x.disc)[1],2,dim(x.disc)[3]))
		x.disc.new[,1,]<-x.disc
		x.disc<-x.disc.new
	}
	
	
	bugs.data.single[[outcome.i]]<-c(bugs.data.single[[outcome.i]],list("ns.single"=ns.single,
			"r.single"=r.single,"n.single"=n.single,"t.single"=t.single,"x.single"=x.single))
	bugs.data.disconnect[[outcome.i]]<-c(bugs.data.disconnect[[outcome.i]],list("ns.disc"=ns.disc,
			"r.disc"=r.disc,"n.disc"=n.disc,"t.disc"=t.disc,"x.disc"=x.disc,"na.disc"=na.disc))
}


# Re-plot the (now disconnected) evidence network
# This should plot the connected and disconnected components together
for(outcome.i in 1:n.outcomes)
{
	# Need treatment matrix in special vector format
	tr.temp<-bugs.data.disconnect[[outcome.i]]$t
	for(i in 1:dim(tr.temp)[1])
	{
	tr.temp[i,]<-i
	}
	# List of trials
	t1<-c(tr.temp[!is.na(bugs.data.disconnect[[outcome.i]]$t)])
	# List of treatments
	t2<-c(bugs.data.disconnect[[outcome.i]]$t[!is.na(bugs.data.disconnect[[outcome.i]]$t)])

	# Plot the connected network with dabigatran nodes removed
	jpeg(file=paste(baseline.directory,"/results/network plots/subnetwork.connected.",outcome.names[outcome.i],".jpg",sep=""),width=960,height=480)
	mtm.networkplot.fun(t1=t1,t2=t2,percomparison=FALSE,nameoftreatments=treatment.names[[outcome.i]][unique(t2)])
	#title(paste("Disconnected",outcome.names[outcome.i]))	
	dev.off()

	tr.temp.disc<-bugs.data.disconnect[[outcome.i]]$t.disc
	for(i in (dim(tr.temp)[1]+1):(dim(tr.temp)[1]+dim(tr.temp.disc)[1]))
	{
	tr.temp.disc[i-dim(tr.temp)[1],]<-i
	}
	t1.disc<-tr.temp.disc[!is.na(bugs.data.disconnect[[outcome.i]]$t.disc)]
	t2.disc<-bugs.data.disconnect[[outcome.i]]$t.disc[!is.na(bugs.data.disconnect[[outcome.i]]$t.disc)]
	t1<-c(t1,tr.temp.disc[!is.na(bugs.data.disconnect[[outcome.i]]$t.disc)])
	t2<-c(t2,bugs.data.disconnect[[outcome.i]]$t.disc[!is.na(bugs.data.disconnect[[outcome.i]]$t.disc)])

	# Plot the network of dabigatran nodes
	jpeg(file=paste(baseline.directory,"/results/network plots/subnetwork.disconnected.",outcome.names[outcome.i],".jpg",sep=""),width=960,height=480)
	mtm.networkplot.fun(t1=t1.disc,t2=t2.disc,percomparison=FALSE,nameoftreatments=treatment.names[[outcome.i]][unique(t2.disc)])
	#title(paste("Disconnected",outcome.names[outcome.i]))	
	dev.off()


	# Plot the disconnected network with all (connected and disconneced) treatments
	jpeg(file=paste(baseline.directory,"/results/network plots/disconnected",outcome.names[outcome.i],".jpg",sep=""),width=960,height=480)
	mtm.networkplot.fun(t1=t1,t2=t2,percomparison=FALSE,nameoftreatments=paste(1:bugs.data[[outcome.i]]$nt,treatment.names[[outcome.i]]))
	title(paste("Disconnected",outcome.names[outcome.i]))	
	dev.off()
}


# Now extract arms from single-arm studies or RCTs that are on reference (tr=1)
for(outcome.i in 1:n.outcomes)
{
	# Data may be used twice - once in NMA, once in baseline model
	# Single and disconected are the same so use either to extract baseline/reference arm data
	r.base<-bugs.data.single[[outcome.i]]$r[which(bugs.data.single[[outcome.i]]$t==1)]
	n.base<-bugs.data.single[[outcome.i]]$n[which(bugs.data.single[[outcome.i]]$t==1)]
	x.base<-matrix(NA,nrow=dim(bugs.data.single[[outcome.i]]$x)[1],ncol=length(n.base))
	rownames(x.base)<-dimnames(bugs.data[[outcome.i]]$x)[[1]]
	for(i.cov in 1:3){
		x.base[i.cov,]<-bugs.data.single[[outcome.i]]$x[i.cov,,][which(bugs.data.single[[outcome.i]]$t==1)]
		# Set the baseline covariate to the mean if it is not reported
		x.base[i.cov,is.na(x.base[i.cov,])]<-mean(x.base[i.cov,],na.rm=TRUE)
	}
	ns.base<-length(r.base)

	bugs.data.single[[outcome.i]]<-c(bugs.data.single[[outcome.i]],list("ns.base"=ns.base,
			"r.base"=r.base,"n.base"=n.base,"x.base"=x.base))

	# Add average for baseline natural history model
	bugs.data.single[[outcome.i]]$x.base.mean<-rowMeans(bugs.data.single[[outcome.i]]$x.base)

	bugs.data.disconnect[[outcome.i]]<-c(bugs.data.disconnect[[outcome.i]],list("ns.base"=ns.base,
			"r.base"=r.base,"n.base"=n.base,"x.base"=x.base))

	# Add average for baseline natural history model
	bugs.data.disconnect[[outcome.i]]$x.base.mean<-rowMeans(bugs.data.disconnect[[outcome.i]]$x.base)

}


# Build RCTs matching a single arm from an existing RCT to the single arm studies
bugs.data.matched<-list()
# Data for the plugin estimators
# This is the same as the single-arm and disconnected RCT dataset but with matched RCTs recorded
bugs.data.disconnect.plugin<-bugs.data.single.plugin<-list()
matched.rct.disc<-matched.rct.single<-list()

# First find the best matches to the single-arm studies
for(outcome.i in 1:n.outcomes)
{
	matched.rct.single[[outcome.i]]<-rep(NA,bugs.data.single[[outcome.i]]$ns.single)
	MATCHEDr<-MATCHEDn<-MATCHEDt<-matrix(NA,nrow=bugs.data.single[[outcome.i]]$ns.single,ncol=2)
	MATCHEDna<-rep(2,bugs.data.single[[outcome.i]]$ns.single)
	MATCHEDns<-bugs.data.single[[outcome.i]]$ns.single
	t.matched<-rep(NA,bugs.data.single[[outcome.i]]$ns.single)
	for(i.single in 1:bugs.data.single[[outcome.i]]$ns.single)
	{
		# Use euclidean distance to choose closest arm match
		distance<-bugs.data.single[[outcome.i]]$t
		for(i.rct in 1:bugs.data.single[[outcome.i]]$ns)
		{
			for(i.arm in 1:bugs.data.single[[outcome.i]]$na[i.rct])
			{
				distance[i.rct,i.arm]<-dist(rbind(x.single[,i.single],
					bugs.data.single[[outcome.i]]$x[,i.rct,i.arm]))
			}
		}
		min.rct<-which(rowSums(distance==min(distance,na.rm=TRUE),na.rm=TRUE)>0)
		matched.rct.single[[outcome.i]][i.single]<-min.rct
		min.arm<-which(colSums(distance==min(distance,na.rm=TRUE),na.rm=TRUE)>0)
		# Which treatment is going to be added? For diagnostics as we need at least one treatment 1 arm...
		t.matched[i.single]<-bugs.data.single[[outcome.i]]$t[min.rct,min.arm]
		MATCHEDt[i.single,]<-c(bugs.data.single[[outcome.i]]$t.single[i.single],bugs.data.single[[outcome.i]]$t[min.rct,min.arm])
		MATCHEDr[i.single,]<-c(bugs.data.single[[outcome.i]]$r.single[i.single],bugs.data.single[[outcome.i]]$r[min.rct,min.arm])
		MATCHEDn[i.single,]<-c(bugs.data.single[[outcome.i]]$n.single[i.single],bugs.data.single[[outcome.i]]$n[min.rct,min.arm])
		# Order the constructed RCT by treatment
		MATCHEDr[i.single,]<-MATCHEDr[i.single,order(MATCHEDt[i.single,])]
		MATCHEDn[i.single,]<-MATCHEDn[i.single,order(MATCHEDt[i.single,])]
		MATCHEDt[i.single,]<-MATCHEDt[i.single,order(MATCHEDt[i.single,])]
	}
	# Define number of treatments in the RCTs and matched (constructed) RCTs
	RCTnt<-length(unique(c(bugs.data.single[[outcome.i]]$t[!is.na(bugs.data.single[[outcome.i]]$t)])))
	MATCHEDnt<-length(unique(c(MATCHEDt[!is.na(MATCHEDt)])))

	# Which treatments are in only matched (1), both matched and RCT (2), or only RCT (3)
	RCTorMatched<-rep(NA, bugs.data.single[[outcome.i]]$nt)
	for(i.treat in 1:bugs.data.single[[outcome.i]]$nt)
	{
		# Matched only
		if(sum(c(i.treat==bugs.data.single[[outcome.i]]$t),na.rm=TRUE)==0 & sum(c(i.treat==MATCHEDt),na.rm=TRUE)>0){ 
		RCTorMatched[i.treat]<-1}
		# Both
		if(sum(c(i.treat==bugs.data.single[[outcome.i]]$t),na.rm=TRUE)>0 & sum(c(i.treat==MATCHEDt),na.rm=TRUE)>0){ 
		RCTorMatched[i.treat]<-2}
		# RCT only
		if(sum(c(i.treat==bugs.data.single[[outcome.i]]$t),na.rm=TRUE)>0 & sum(c(i.treat==MATCHEDt),na.rm=TRUE)==0){ 
		RCTorMatched[i.treat]<-3}
	}

	# Build the data for Joy Leahy's matchined BUGS script
	bugs.data.matched[[outcome.i]]<-list("RCTns"=bugs.data.single[[outcome.i]]$ns,"RCTnt"=RCTnt,"RCTna"=bugs.data.single[[outcome.i]]$na,
			"MATCHEDns"=MATCHEDns,"COnt"=bugs.data.single[[outcome.i]]$nt,"MATCHEDna"=MATCHEDna,"MATCHEDnt"=MATCHEDnt,
			"RCTn"=bugs.data.single[[outcome.i]]$n,"RCTr"=bugs.data.single[[outcome.i]]$r,"RCTt"=bugs.data.single[[outcome.i]]$t,
			"MATCHEDn"=MATCHEDn,"MATCHEDr"=MATCHEDr,"MATCHEDt"=MATCHEDt,
			"RCTMultiplier"=1,"MatchedMultiplier"=1,"RCTorMatched"=RCTorMatched)

	bugs.data.single.plugin[[outcome.i]]<-bugs.data.single[[outcome.i]]
	# Add a trailing NA to avoid BUGS getting confused between scalars and vectors
	bugs.data.single.plugin[[outcome.i]]$matched.rct<-c(matched.rct.single[[outcome.i]],NA)
}
names(matched.rct.single)<-names(bugs.data.single.plugin)<-names(bugs.data.matched)<-outcome.names

# Now find the best matches to the disconnected RCTs
for(outcome.i in 1:n.outcomes)
{
	matched.rct.disc[[outcome.i]]<-rep(NA,bugs.data.disconnect[[outcome.i]]$ns.disc)
	for(i.disc in 1:bugs.data.disconnect[[outcome.i]]$ns.disc)
	{
		# Use euclidean distance to choose closest arm match
		distance<-rep(NA,bugs.data.disconnect[[outcome.i]]$ns)
		for(i.rct in 1:bugs.data.disconnect[[outcome.i]]$ns)
		{
			# Distance between unweighted average of arms of each trial
			distance[i.rct]<-dist(rbind(
				rowMeans(x.disc[,i.disc,1:bugs.data.disconnect[[outcome.i]]$na.disc[i.disc]]),
				rowMeans(bugs.data.disconnect[[outcome.i]]$x[,i.rct,1:bugs.data.disconnect[[outcome.i]]$na[i.rct]])
				))
		}
		min.rct<-which.min(distance)
		matched.rct.disc[[outcome.i]][i.disc]<-min.rct
	}

	bugs.data.disconnect.plugin[[outcome.i]]<-bugs.data.disconnect[[outcome.i]]
	bugs.data.disconnect.plugin[[outcome.i]]$matched.rct<-c(matched.rct.disc[[outcome.i]],NA)
}
names(matched.rct.disc)<-names(bugs.data.disconnect.plugin)<-outcome.names

# Function to format mean and 95% credible interval with fewer digits
format.results<-function(x)
{
	return(paste(format(x[1],digits=3)," (",format(x[2],digits=3),", ",format(x[3],digits=3),")",sep=""))
}



# Compare deviance and DIC of baseline models
# Only needed to run this once.
if(0)
{
	source("compare.baseline.models.2.R")
}

n.chains<-2			 # 2
num.sims=30000*n.chains  # 150000*n.chains
burn.in=30000*n.chains	 # 100000*n.chains

# Lists to store the bugs inits and outputs
bugs.inits<-bugs.inits.re<-bugs.inits.single<-bugs.inits.single.re<-bugs.inits.single.3cov<-bugs.inits.single.3cov.re<-list()
bugs.inits.disc<-bugs.inits.disc.re<-bugs.inits.disc.3cov<-bugs.inits.disc.3cov.re<-list()
bugs.object.fe<-bugs.object.re<-bugs.object.rebase.single.fe<-bugs.object.rebase.single.re<-list()
bugs.object.single.rct.only.re<-bugs.object.single.rct.only.fe<-list()
bugs.object.rebase.single.chads2.fe<-bugs.object.rebase.single.chads2.re<-list()
bugs.object.rebase.disc.fe<-bugs.object.rebase.disc.re<-list()
bugs.object.rebase.disc.chads2.fe<-bugs.object.rebase.disc.chads2.re<-list()
bugs.object.single.rct.only.fe<-bugs.object.single.rct.only.re<-list()
bugs.object.disc.rct.only.fe<-bugs.object.disc.rct.only.re<-list()
bugs.object.plugin.single.fe<-bugs.object.plugin.single.re<-list()
bugs.object.plugin.disc.fe<-bugs.object.plugin.disc.re<-list()
bugs.inits.matched<-bugs.object.matched<-list()

for(outcome.i in 1:n.outcomes)
{
	print(outcome.names[outcome.i])
	# Initial values for models using all (only) RCT data
	inits1<-list(d=c(NA,rep(0.5,bugs.data[[outcome.i]]$nt-1)),mu=rep(0.5,bugs.data[[outcome.i]]$ns))
	inits2<-list(d=c(NA,rep(-0.5,bugs.data[[outcome.i]]$nt-1)),mu=rep(1,bugs.data[[outcome.i]]$ns))
	inits1.re<-list(d=c(NA,rep(0.5,bugs.data[[outcome.i]]$nt-1)),mu=rep(0.5,bugs.data[[outcome.i]]$ns),sd=1)
	inits2.re<-list(d=c(NA,rep(-0.5,bugs.data[[outcome.i]]$nt-1)),mu=rep(1,bugs.data[[outcome.i]]$ns),sd=0.5)
	bugs.inits[[outcome.i]]<-list(inits1,inits2)
	bugs.inits.re[[outcome.i]]<-list(inits1.re,inits2.re)

	# Initial values for models using mixture of RCT and single-arm data
	inits1<-list(d=c(NA,rep(1,bugs.data.single[[outcome.i]]$nt-1)),mu=rep(-0.5,bugs.data.single[[outcome.i]]$ns),mu.single=rep(-0.5,bugs.data.single[[outcome.i]]$ns.single),m=0.1,sd.m=1)
	inits2<-list(d=c(NA,rep(-1,bugs.data.single[[outcome.i]]$nt-1)),mu=rep(-1,bugs.data.single[[outcome.i]]$ns),mu.single=rep(-1,bugs.data.single[[outcome.i]]$ns.single),m=0.1,sd.m=1)
	inits1.re<-list(d=c(NA,rep(1,bugs.data.single[[outcome.i]]$nt-1)),mu=rep(-0.5,bugs.data.single[[outcome.i]]$ns),mu.single=rep(-0.5,bugs.data.single[[outcome.i]]$ns.single),sd=1,m=0.1,sd.m=1)
	inits2.re<-list(d=c(NA,rep(-1,bugs.data.single[[outcome.i]]$nt-1)),mu=rep(-1,bugs.data.single[[outcome.i]]$ns),mu.single=rep(-1,bugs.data.single[[outcome.i]]$ns.single),sd=0.5,m=0.1,sd.m=1)
	bugs.inits.single[[outcome.i]]<-list(inits1,inits2)
	bugs.inits.single.re[[outcome.i]]<-list(inits1.re,inits2.re)

	# Initial values for models using all RCT data single-arm studies and fixed effects on baseline with (up to) 3 covariates
	inits1<-list(d=c(NA,rep(1,bugs.data.single[[outcome.i]]$nt-1)),mu=rep(-0.5,bugs.data.single[[outcome.i]]$ns),m=0.1,sd.m=1,mu.base=rep(0.1,bugs.data.single[[outcome.i]]$ns.base),beta.base=rep(0.5,3))
	inits2<-list(d=c(NA,rep(-1,bugs.data.single[[outcome.i]]$nt-1)),mu=rep(-1,bugs.data.single[[outcome.i]]$ns),m=0.5,sd.m=2,mu.base=rep(0.5,bugs.data.single[[outcome.i]]$ns.base),beta.base=rep(-0.5,3))
	inits1.re<-list(d=c(NA,rep(1,bugs.data.single[[outcome.i]]$nt-1)),mu=rep(-0.5,bugs.data.single[[outcome.i]]$ns),m=0.1,sd.m=1,mu.base=rep(0.1,bugs.data.single[[outcome.i]]$ns.base),sd=1,beta.base=rep(0.5,3))
	inits2.re<-list(d=c(NA,rep(-1,bugs.data.single[[outcome.i]]$nt-1)),mu=rep(-1,bugs.data.single[[outcome.i]]$ns),m=0.5,sd.m=2,mu.base=rep(0.5,bugs.data.single[[outcome.i]]$ns.base),sd=0.5,beta.base=rep(-0.5,3))
	bugs.inits.single.3cov[[outcome.i]]<-list(inits1,inits2)
	bugs.inits.single.3cov.re[[outcome.i]]<-list(inits1.re,inits2.re)

	# Initial values for models using at least two disconnected networks of RCTs and random effects on baseline
	inits1<-list(d=c(NA,rep(.1,bugs.data.disconnect[[outcome.i]]$nt-1)),mu=rep(0.5,bugs.data.disconnect[[outcome.i]]$ns),mu.disc=rep(0.5,bugs.data.disconnect[[outcome.i]]$ns.disc),m=0.1,sd.m=1)
	inits2<-list(d=c(NA,rep(0.05,bugs.data.disconnect[[outcome.i]]$nt-1)),mu=rep(0.1,bugs.data.disconnect[[outcome.i]]$ns),mu.disc=rep(0.1,bugs.data.disconnect[[outcome.i]]$ns.disc),m=0.1,sd.m=1)
	inits1.re<-list(d=c(NA,rep(1,bugs.data.disconnect[[outcome.i]]$nt-1)),mu=rep(-0.5,bugs.data.disconnect[[outcome.i]]$ns),mu.disc=rep(-0.5,bugs.data.disconnect[[outcome.i]]$ns.disc),sd=1,m=0.1,sd.m=1)
	inits2.re<-list(d=c(NA,rep(-1,bugs.data.disconnect[[outcome.i]]$nt-1)),mu=rep(-1,bugs.data.disconnect[[outcome.i]]$ns),mu.disc=rep(-1,bugs.data.disconnect[[outcome.i]]$ns.disc),sd=0.5,m=0.1,sd.m=1)
	bugs.inits.disc[[outcome.i]]<-list(inits1,inits2)
	bugs.inits.disc.re[[outcome.i]]<-list(inits1.re,inits2.re)

	# Initial values for models using at least two disconnected networks of RCTs and fixed effects on baseline with (up to) 3 covariates
	inits1<-list(d=c(NA,rep(1,bugs.data.disconnect[[outcome.i]]$nt-1)),mu=rep(-0.5,bugs.data.disconnect[[outcome.i]]$ns),m=0.1,sd.m=1,mu.base=rep(0.1,bugs.data.disconnect[[outcome.i]]$ns.base),beta.base=rep(0.5,3))
	inits2<-list(d=c(NA,rep(-1,bugs.data.disconnect[[outcome.i]]$nt-1)),mu=rep(-1,bugs.data.disconnect[[outcome.i]]$ns),m=0.5,sd.m=2,mu.base=rep(0.5,bugs.data.disconnect[[outcome.i]]$ns.base),beta.base=rep(-0.5,3))
	inits1.re<-list(d=c(NA,rep(1,bugs.data.disconnect[[outcome.i]]$nt-1)),mu=rep(-0.5,bugs.data.disconnect[[outcome.i]]$ns),m=0.1,sd.m=1,mu.base=rep(0.1,bugs.data.disconnect[[outcome.i]]$ns.base),sd=1,beta.base=rep(0.5,3))
	inits2.re<-list(d=c(NA,rep(-1,bugs.data.disconnect[[outcome.i]]$nt-1)),mu=rep(-1,bugs.data.disconnect[[outcome.i]]$ns),m=0.5,sd.m=2,mu.base=rep(0.5,bugs.data.disconnect[[outcome.i]]$ns.base),sd=0.5,beta.base=rep(-0.5,3))
	bugs.inits.disc.3cov[[outcome.i]]<-list(inits1,inits2)
	bugs.inits.disc.3cov.re[[outcome.i]]<-list(inits1.re,inits2.re)

	# Initial values for aggregate level matching
	inits1<-list(d=c(NA,rep(.1,bugs.data.matched[[outcome.i]]$COnt-1)),RCTmu=rep(0.05,bugs.data.matched[[outcome.i]]$RCTns),RCTsd=.25,MATCHEDmu=rep(0.05,bugs.data.matched[[outcome.i]]$MATCHEDns),MATCHEDsd=.25)
	inits2<-list(d=c(NA,rep(0.05,bugs.data.matched[[outcome.i]]$COnt-1)),RCTmu=rep(0.1,bugs.data.matched[[outcome.i]]$RCTns),RCTsd=0.5,MATCHEDmu=rep(0.1,bugs.data.matched[[outcome.i]]$MATCHEDns),MATCHEDsd=0.5)
	bugs.inits.matched[[outcome.i]]<-list(inits1,inits2)


	# Independent baselines (All RCT data). This is standard NICE approach
	bugs.object.fe[[outcome.i]]<-bugs(data=bugs.data[[outcome.i]],inits=bugs.inits[[outcome.i]],parameters.to.save=c("mu","d","totresdev"),model=model.independent.baseline.fe,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)
	bugs.object.re[[outcome.i]]<-bugs(data=bugs.data[[outcome.i]],inits=bugs.inits.re[[outcome.i]],parameters.to.save=c("mu","d","sd","totresdev"),model=model.independent.baseline.re,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=FALSE)

	# Independent baselines (Only RCT data after single studies removed). This is standard NICE approach
	bugs.object.single.rct.only.fe[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.single[[outcome.i]],parameters.to.save=c("mu","d","totresdev"),model=model.independent.baseline.fe,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)
	bugs.object.single.rct.only.re[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.single.re[[outcome.i]],parameters.to.save=c("mu","d","sd","totresdev"),model=model.independent.baseline.re,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=FALSE)

	# Independent baselines (Only RCT data after disconnected studies removed). This is standard NICE approach
	bugs.object.disc.rct.only.fe[[outcome.i]]<-bugs(data=bugs.data.disconnect[[outcome.i]],inits=bugs.inits.disc[[outcome.i]],parameters.to.save=c("mu","d","totresdev"),model=model.independent.baseline.fe,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)
	bugs.object.disc.rct.only.re[[outcome.i]]<-bugs(data=bugs.data.disconnect[[outcome.i]],inits=bugs.inits.disc.re[[outcome.i]],parameters.to.save=c("mu","d","sd","totresdev"),model=model.independent.baseline.re,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=1,debug=FALSE)

	# Take the mu from the matched RCT
	# Adding NA at end to avoid confusion between vectors and scalars
	bugs.data.single.plugin[[outcome.i]]$mu.plugin.mean<-
		c(bugs.object.single.rct.only.fe[[outcome.i]]$summary[bugs.data.single.plugin[[outcome.i]]$matched.rct,"mean"],NA)
	bugs.data.single.plugin[[outcome.i]]$mu.plugin.prec<-
		1/c(bugs.object.single.rct.only.fe[[outcome.i]]$summary[bugs.data.single.plugin[[outcome.i]]$matched.rct,"sd"],NA)^2
	bugs.data.disconnect.plugin[[outcome.i]]$mu.plugin.mean<-
		c(bugs.object.disc.rct.only.fe[[outcome.i]]$summary[bugs.data.disconnect.plugin[[outcome.i]]$matched.rct,"mean"],NA)
	bugs.data.disconnect.plugin[[outcome.i]]$mu.plugin.prec<-
		1/c(bugs.object.disc.rct.only.fe[[outcome.i]]$summary[bugs.data.disconnect.plugin[[outcome.i]]$matched.rct,"sd"],NA)^2

	# And informative priors for sd.disc from the connected RCTs
	# This is for both plug-in and reference prediction models
	bugs.data.disconnect.plugin[[outcome.i]]$sd.connected.mean<-bugs.data.disconnect[[outcome.i]]$sd.connected.mean<-
		bugs.object.single.rct.only.re[[outcome.i]]$summary["sd","mean"]
	bugs.data.disconnect.plugin[[outcome.i]]$sd.connected.tau<-bugs.data.disconnect[[outcome.i]]$sd.connected.tau<-
		1/bugs.object.single.rct.only.re[[outcome.i]]$summary["sd","sd"]^2
	bugs.data.single.plugin[[outcome.i]]$sd.connected.mean<-bugs.data.single[[outcome.i]]$sd.connected.mean<-
		bugs.object.single.rct.only.re[[outcome.i]]$summary["sd","mean"]
	bugs.data.single.plugin[[outcome.i]]$sd.connected.tau<-bugs.data.single[[outcome.i]]$sd.connected.tau<-
		1/bugs.object.single.rct.only.re[[outcome.i]]$summary["sd","sd"]^2


	# Plug-in estimator (single-arm studies)
	print("Single arms plug-in estimator")
	bugs.object.plugin.single.fe[[outcome.i]]<-bugs(data=bugs.data.single.plugin[[outcome.i]],inits=bugs.inits.single[[outcome.i]],parameters.to.save=c("mu","d","totresdev"),model=model.plugin.single.fe,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)
	# This worked for outcome.i 1 and 4 (although sd uniform prior had to be bounded above by 2)
	bugs.object.plugin.single.re[[outcome.i]]<-bugs(data=bugs.data.single.plugin[[outcome.i]],inits=bugs.inits.single.re[[outcome.i]],parameters.to.save=c("mu","d","sd","sd.disc","totresdev"),model.file=paste0(baseline.directory,"/code/plugin.model.single.re.3.txt"),	
		clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)


	# Plug-in estimator (disconnected RCTs)
	print("Plug-in estimator for disconnected RCTs")
	bugs.object.plugin.disc.fe[[outcome.i]]<-bugs(data=bugs.data.disconnect.plugin[[outcome.i]],inits=bugs.inits.disc[[outcome.i]],parameters.to.save=c("mu","d","totresdev.disc","totresdev"),model=model.plugin.disc.fe,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)
	# Note that bugs.seed was 10 for stroke but changed to 5 for bleed as it crashed with 10. 
	bugs.object.plugin.disc.re[[outcome.i]]<-bugs(data=bugs.data.disconnect.plugin[[outcome.i]],inits=bugs.inits.disc.re[[outcome.i]],parameters.to.save=c("mu","d","sd","sd.disc","totresdev.disc","totresdev"),model.file=paste0(baseline.directory,"/code/plugin.model.disconnected.re.3.txt"),	
		clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=5,n.thin=1,debug=FALSE)


	# Random effects baseline (some single-arm studies)
	print("Single arms rebase")
	bugs.object.rebase.single.fe[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.single[[outcome.i]],parameters.to.save=c("mu","d","m","sd.m","mu.single","totresdev"),model=model.random.effects.baseline.single.fe,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)
	# This worked for outcome.i 1 and 4 (although sd uniform prior had to be bounded above by 2)
	bugs.object.rebase.single.re[[outcome.i]]<-bugs(data=bugs.data.single[[outcome.i]],inits=bugs.inits.single.re[[outcome.i]],parameters.to.save=c("mu","d","m","sd.m","mu.single","sd","sd.disc","totresdev"),model.file=paste0(baseline.directory,"/code/random.effects.baseline.model.single-arms.re.3.txt"),	
		clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)


	# Random effects baseline (with disconnected RCTs)
	print("Disconnected rebase")

	# Initial values need collation operator error
	if(length(bugs.inits.disc.re[[outcome.i]][[1]]$mu.disc)==1){
		bugs.inits.disc[[outcome.i]][[1]]$mu.disc<-bugs.inits.disc[[outcome.i]][[2]]$mu.disc<-NULL
		bugs.inits.disc.re[[outcome.i]][[1]]$mu.disc<-bugs.inits.disc.re[[outcome.i]][[2]]$mu.disc<-NULL
	}
	bugs.object.rebase.disc.fe[[outcome.i]]<-bugs(data=bugs.data.disconnect[[outcome.i]],inits=bugs.inits.disc[[outcome.i]],parameters.to.save=c("mu","d","m","sd.m","mu.disc","totresdev.disc","totresdev"),model=model.random.effects.baseline.disc.fe,clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)
	# The below worked for both outcome.i 1
	# Ran into update error for node <delta[4,2]>  for outcome.i 4
	bugs.object.rebase.disc.re[[outcome.i]]<-bugs(data=bugs.data.disconnect[[outcome.i]],inits=bugs.inits.disc.re[[outcome.i]],parameters.to.save=c("mu","d","m","sd.m","mu.disc","sd","sd.disc","totresdev.disc","totresdev"),model.file=paste0(baseline.directory,"/code/random.effects.baseline.model.disconnected.re.4.txt"),	
		clearWD=TRUE,summary.only=FALSE,n.iter=(num.sims+burn.in),n.burnin=burn.in,n.chains=n.chains,bugs.seed=1,n.thin=10,debug=FALSE)
}

save(file=paste(baseline.directory,"/results/bugs objects/bugs.objects.4.rda",sep=""),bugs.object.rebase.disc.chads2.re,bugs.object.matched,
		bugs.object.disc.rct.only.fe, bugs.object.disc.rct.only.re, bugs.object.single.rct.only.fe, bugs.object.single.rct.only.re,
		bugs.object.rebase.disc.re,bugs.object.rebase.disc.fe,
		bugs.object.rebase.single.chads2.re,
		bugs.object.rebase.single.re,bugs.object.rebase.single.fe,
		bugs.object.re,bugs.object.fe,
		bugs.object.plugin.disc.fe,bugs.object.plugin.disc.re,
		bugs.object.plugin.single.fe,bugs.object.plugin.single.re,
		bugs.data.disconnect.plugin, 	bugs.data.single.plugin)
#load(file=paste(baseline.directory,"/results/bugs objects/bugs.objects.4.rda",sep=""))

# Need some comparison
comparison.matrix<-matrix(nrow=bugs.data[[outcome.i]]$nt,ncol=5)
colnames(comparison.matrix)<-c("Single arm?","Fixed effects","Random effects","Rbase fixed effects","Rbase random effects")
for(i in 1:bugs.data[[outcome.i]]$nt)
{
	comparison.matrix[i,"Single arm?"]<-sum(tr.single==i)>0
	comparison.matrix[i,"Random effects"]<-format.results(bugs.object.re[[outcome.i]]$summary[i,c(1,3,7)])
	comparison.matrix[i,"Rbase random effects"]<-format.results(bugs.object.rebase.single.re[[outcome.i]]$summary[i,c(1,3,7)])
	comparison.matrix[i,"Fixed effects"]<-format.results(bugs.object.fe[[outcome.i]]$summary[i,c(1,3,7)])
	comparison.matrix[i,"Rbase fixed effects"]<-format.results(bugs.object.rebase.single.fe[[outcome.i]]$summary[i,c(1,3,7)])
}

# Generate forest plots comparing treatment effects under different models
source("plot.model.comparison.4.R")
source("plot.model.comparison.re.2.R")

