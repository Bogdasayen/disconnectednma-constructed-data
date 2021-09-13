# v2 splits the RCT (connected) only and single-arm (disconnected) into
# two separate (more readable) plots.
# After V3 created a separate script for random study effects

format.effect.summary<-function(x,n.digits=2)
{
  return(paste0(format(x[1],digits=n.digits, nsmall = n.digits)," (",format(x[2],digits=n.digits, nsmall = n.digits),", ",format(x[3],digits=n.digits, nsmall = n.digits),")"))
}


# Single-arm NMA first #################################################
cex.effect.estimate <- 0.45

for(outcome.i in c(1,4))
{
	# Number of treatments
	n.treat<-bugs.data.single[[outcome.i]]$nt

	# Categorise treatments by evidence available
	single.treats<-unique(bugs.data.single[[outcome.i]]$t.single)
	rct.treats<-unique(c(bugs.data.single[[outcome.i]]$t[!is.na(bugs.data.single[[outcome.i]]$t)]))
	# Treatments in both single and rct studies
	both.treats<-single.treats[which(is.element(single.treats,rct.treats))]
	# Treatments only in single arm studies
	single.only.treats<-single.treats[!is.element(single.treats,both.treats)]
	# Treatments only in rct studies
	rct.only.treats<-rct.treats[!is.element(rct.treats,both.treats)]

	# Treatment lists are shared across outputs
	r.names<-rownames(bugs.object.re[[outcome.i]]$summary)
	
	maxima<-c(bugs.object.re[[outcome.i]]$summary[grep("d",r.names),7][1:n.treat-1],
		bugs.object.rebase.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.single.re[[outcome.i]]$summary)),7][1:n.treat-1],
		bugs.object.plugin.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.single.re[[outcome.i]]$summary)),7][1:n.treat-1])
	minima<-c(bugs.object.re[[outcome.i]]$summary[grep("d",r.names),3][1:n.treat-1],
		bugs.object.rebase.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.single.re[[outcome.i]]$summary)),3][1:n.treat-1],
		bugs.object.plugin.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.single.re[[outcome.i]]$summary)),3][1:n.treat-1])

	# Ignore treatments with no data
	x.max<-max(maxima[maxima<50])
	x.min<-min(minima[minima>-50])


	# Open a pdf file for output
	pdf(file=paste(baseline.directory,"/results/forest plots/v2.single.arms.rcts.",outcome.names[outcome.i],".forest.re.pdf",sep=""))

	# This depends on what is being plotted
	if(outcome.i!=4){cex.text<-0.9}
	if(outcome.i==4){cex.text<-0.9}
	
	y.min=0; y.max=3*length(rct.only.treats)+3
	
	
	# Store the old plot parameters
	old.mar <- par()$mar
	old.fig <- par()$fig
	old.mfrow <- par()$mfrow
	old.scipen <- options()$scipen
	# Turn off scientific notation
	options(scipen = 999)
	
	# Create 3 plots
	# One for treatment names, one for the forest plot, and one for effect summaries
	par(mfrow = c(1,3))
	par(fig = c(0, 0.3, 0, 1))
	par(mar=c(5.1,1,4.1,0.1))
	
	# Plot treatment names
	plot(c(0,0), xlab = "", col = 0, xlim = c(0, 4), ylim = c(y.min, y.max), ylab = "", axes = FALSE)
	
	# Plot the RCT only treatment effects
	for(i.treat in rct.only.treats)
	{
	  text(treatment.names[[outcome.i]][i.treat],x=0,y=3*(which(rct.only.treats==i.treat)),pos=4,cex=cex.text)
	}
	
	
  # Plot lines and points
	# Build plot area
	par(fig = c(0.3, 0.7, 0, 1), new = TRUE)
	par(mar = c(5.1, 0.1, 4.1, 0.1))
	
	plot(c(0,0),col=0,xlim=c(x.min,x.max),ylim=c(y.min,y.max),ylab="",xlab=paste("Log odds ratio vs",treatment.names[[1]][1]),axes=FALSE)
	axis(side=1)	
	title(paste(outcome.names[outcome.i],"True (all RCTs)"), cex.main = cex.text)

	# Plot the RCT only treatment effects
	for(i.treat in rct.only.treats)
	{
	
		# All RCTs ('true' effects)
		points(y=3*(which(rct.only.treats==i.treat)),x=bugs.object.single.rct.only.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.single.rct.only.re[[outcome.i]]$summary)),1][i.treat-1])
		lines(y=rep(3*(which(rct.only.treats==i.treat)),2),x=bugs.object.single.rct.only.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.single.rct.only.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,])

		# RE on baseline
		points(y=-0.5+3*(which(rct.only.treats==i.treat)),x=bugs.object.rebase.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.single.re[[outcome.i]]$summary)),1][i.treat-1],col=2)
		lines(y=rep(-0.5+3*(which(rct.only.treats==i.treat)),2),x=bugs.object.rebase.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.single.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,],col=2)

		# Plug-in method
		points(y=-1+3*(which(rct.only.treats==i.treat)),x=bugs.object.plugin.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.single.re[[outcome.i]]$summary)),1][i.treat-1],col=3)
		lines(y=rep(-1+3*(which(rct.only.treats==i.treat)),2),x=bugs.object.plugin.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.single.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,],col=3)
	}
	legend("topright",legend=c("True (all RCTs)","Reference prediction","Aggregate level matching"),lty=1,col=c(1,2,3),cex=cex.text)

	
	par(fig = c(0.7, 1, 0, 1), new = TRUE)
	
	par(mar=c(5.1,1,4.1,1))
	text.plot.x <- c(0,6)
	plot(c(0, 0), xlab = "", col = 0, xlim = text.plot.x, ylim = c(y.min, y.max), ylab = "",
	     axes=FALSE)
	
	#lines(x = rep(1.1*text.plot.x[2]/3, 2), y = c(y.min, y.max))
	#lines(x = rep(2*text.plot.x[2]/3, 2), y = c(y.min, y.max))
	text(x=text.plot.x[2]/6, y = y.max, pos = 1, "True (all RCTs)", cex = cex.effect.estimate, col = 1)
	text(x = 3*text.plot.x[2]/6, y = y.max, pos = 1, "RP", cex = cex.effect.estimate, col = 2)
	text(x = 5*text.plot.x[2]/6, y = y.max, pos = 1,"ALM", cex = cex.effect.estimate, col = 3)
	for(i.treat in rct.only.treats) {
	  text(y=3*(which(rct.only.treats==i.treat)), x = 0, format.effect.summary(bugs.object.single.rct.only.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.single.rct.only.re[[outcome.i]]$summary)),c(1, 3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 1)
	  text(y=3*(which(rct.only.treats==i.treat)), x = text.plot.x[2]/3, format.effect.summary(bugs.object.rebase.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.single.re[[outcome.i]]$summary)),c(1, 3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 2)
    text(y=3*(which(rct.only.treats==i.treat)), x = 2*text.plot.x[2]/3, format.effect.summary(bugs.object.plugin.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.single.re[[outcome.i]]$summary)),c(1,3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 3)
	}
	
	# End export of plot
	dev.off()
	
	
	# Reset parameters
	par(mar = old.mar)
	par(fig = old.fig)
	par(mfrow = old.mfrow)
	options(scipen = old.scipen)
	


	# Open a pdf file for output
	pdf(file=paste(baseline.directory,"/results/forest plots/v2.single.arms.single.",outcome.names[outcome.i],".forest.re.pdf",sep=""))

	# This depends on what is being plotted
	if(outcome.i!=4){cex.text<-0.9}
	if(outcome.i==4){cex.text<-0.65}
	
	# Redefine the y limits
	y.min=0; y.max=3*length(single.only.treats)+3

	
	# Store the old plot parameters
	old.mar <- par()$mar
	old.fig <- par()$fig
	old.mfrow <- par()$mfrow
	old.scipen <- options()$scipen
	# Turn off scientific notation
	options(scipen = 999)
	
	
	# Create 3 plots
	# One for treatment names, one for the forest plot, and one for effect summaries
	par(mfrow = c(1,3))
	par(fig = c(0, 0.3, 0, 1))
	par(mar=c(5.1,1,4.1,0.1))
	
	# Plot treatment names
	plot(c(0,0), xlab = "", col = 0, xlim = c(0, 4), ylim = c(y.min, y.max), ylab = "", axes = FALSE)
	
	# Plot the RCT only treatment effects
	for(i.treat in single.only.treats)
	{
	  text(treatment.names[[outcome.i]][i.treat],x=0,y=3*(which(single.only.treats==i.treat)),pos=4,cex=cex.text)
	}
	
	
	# Plot lines and points
	# Build plot area
	par(fig = c(0.3, 0.7, 0, 1), new = TRUE)
	par(mar = c(5.1, 0.1, 4.1, 0.1))
	
	# Build plot area
	plot(c(0,0),col=0,xlim=c(x.min,x.max),ylim=c(y.min,y.max),axes=FALSE,ylab="",xlab=paste("Log odds ratio vs",treatment.names[[1]][1]))
	axis(side=1)	
	title(paste(outcome.names[outcome.i],"single-arm only"), cex.main = cex.text)

	# Plot the Single-arm only treatment effects
	for(i.treat in single.only.treats)
	{

		# All RCTs ('true' effects)
		points(y=3*(which(single.only.treats==i.treat)),x=bugs.object.re[[outcome.i]]$summary[grep("d",r.names),1][i.treat-1])
		lines(y=rep(3*(which(single.only.treats==i.treat)),2),x=bugs.object.re[[outcome.i]]$summary[grep("d",r.names),c(3,7)][i.treat-1,])

		# RE on baseline
		points(y=-0.5+3*(which(single.only.treats==i.treat)),x=bugs.object.rebase.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.single.re[[outcome.i]]$summary)),1][i.treat-1],col=2)
		lines(y=rep(-0.5+3*(which(single.only.treats==i.treat)),2),x=bugs.object.rebase.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.single.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,],col=2)

		# Plug-in method
		points(y=-1+3*(which(single.only.treats==i.treat)),x=bugs.object.plugin.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.single.re[[outcome.i]]$summary)),1][i.treat-1],col=3)
		lines(y=rep(-1+3*(which(single.only.treats==i.treat)),2),x=bugs.object.plugin.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.single.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,],col=3)
	}
	
	legend("topright",legend=c("True (all RCTs)","Reference prediction","Aggregate level matching"),lty=1,col=c(1,2,3),cex=cex.text)
	
	
	par(fig = c(0.7, 1, 0, 1), new = TRUE)
	
	par(mar=c(5.1,1,4.1,1))
	text.plot.x <- c(0,6)
	plot(c(0, 0), xlab = "", col = 0, xlim = text.plot.x, ylim = c(y.min, y.max), ylab = "",
	     axes=FALSE)
	
	#lines(x = rep(1.1*text.plot.x[2]/3, 2), y = c(y.min, y.max))
	#lines(x = rep(2*text.plot.x[2]/3, 2), y = c(y.min, y.max))
	text(x=text.plot.x[2]/6, y = y.max, pos = 1, "True (all RCTs)", cex = cex.effect.estimate, col = 1)
	text(x = 3*text.plot.x[2]/6, y = y.max, pos = 1, "RP", cex = cex.effect.estimate, col = 2)
	text(x = 5*text.plot.x[2]/6, y = y.max, pos = 1,"ALM", cex = cex.effect.estimate, col = 3)
	for(i.treat in single.only.treats) {
	  text(y=3*(which(single.only.treats==i.treat)), x = 0, 
	       format.effect.summary(bugs.object.re[[outcome.i]]$summary[grep("d",r.names),c(1,3,7)][i.treat-1,]),
	       pos = 4, cex = cex.effect.estimate, col = 1)
	  text(y=3*(which(single.only.treats==i.treat)), x = text.plot.x[2]/3, 
	       format.effect.summary(bugs.object.rebase.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.single.re[[outcome.i]]$summary)),c(1,3,7)][i.treat-1,]), 
	       pos = 4, cex = cex.effect.estimate, col = 2)
	  text(y=3*(which(single.only.treats==i.treat)), x = 2*text.plot.x[2]/3, format.effect.summary(bugs.object.plugin.single.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.single.re[[outcome.i]]$summary)),c(1,3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 3)
	}
	
	
	dev.off()
	# Reset parameters
	par(mar = old.mar)
	par(fig = old.fig)
	par(mfrow = old.mfrow)
	options(scipen = old.scipen)
	
}

# Then disconnected RCT networks

for(outcome.i in c(1,4))
{
	# Number of treatments
	n.treat<-bugs.data.disconnect[[outcome.i]]$nt

	# Categorise treatments by evidence available
	disconnect.treats<-unique(c(bugs.data.disconnect[[outcome.i]]$t.disc[!is.na(bugs.data.disconnect[[outcome.i]]$t.disc)]))
	rct.treats<-unique(c(bugs.data.disconnect[[outcome.i]]$t[!is.na(bugs.data.disconnect[[outcome.i]]$t)]))
	# Treatments in both single and rct studies
	both.treats<-disconnect.treats[which(is.element(disconnect.treats,rct.treats))]
	# Treatments only in single arm studies
	disconnect.only.treats<-disconnect.treats[!is.element(disconnect.treats,both.treats)]
	# Treatments only in rct studies
	rct.only.treats<-rct.treats[!is.element(rct.treats,both.treats)]

	# Treatment lists are shared across outputs
	r.names<-rownames(bugs.object.re[[outcome.i]]$summary)
	
	maxima<-c(bugs.object.re[[outcome.i]]$summary[grep("d",r.names),7][1:n.treat-1],
		bugs.object.rebase.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.disc.re[[outcome.i]]$summary)),7][1:n.treat-1],
		bugs.object.plugin.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.disc.re[[outcome.i]]$summary)),7][1:n.treat-1])
	minima<-c(bugs.object.re[[outcome.i]]$summary[grep("d",r.names),3][1:n.treat-1],
		bugs.object.rebase.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.disc.re[[outcome.i]]$summary)),3][1:n.treat-1],
		bugs.object.plugin.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.disc.re[[outcome.i]]$summary)),3][1:n.treat-1])
	# Ignore treatments with no data
	x.max<-max(maxima[maxima<10])
	x.min<-min(minima[minima>-10])


	# Open a pdf file for output
	pdf(file=paste(baseline.directory,"/results/forest plots/v2.disconnected.connected",outcome.names[outcome.i],".forest.re.pdf",sep=""))

	# This depends on what is being plotted
	if(outcome.i!=4){cex.text<-.9}
	if(outcome.i==4){cex.text<-0.75}
	
	
	y.min=0; y.max=2*length(rct.only.treats)+3
	
	# Store the old plot parameters
	old.mar <- par()$mar
	old.fig <- par()$fig
	old.mfrow <- par()$mfrow
	old.scipen <- options()$scipen
	# Turn off scientific notation
	options(scipen = 999)
	
	par(mfrow = c(1,3))
	par(fig = c(0, 0.3, 0, 1))
	par(mar=c(5.1,1,4.1,0.1))
	
	# Plot treatment names
	plot(c(0,0), xlab = "", col = 0, xlim = c(0, 4), ylim = c(y.min, y.max), ylab = "", axes = FALSE)
	
	# Plot the RCT only treatment effects
	for(i.treat in rct.only.treats[-1])
	{
	  text(treatment.names[[outcome.i]][i.treat],x=0,y=2*(which(rct.only.treats==i.treat)),pos=4,cex=cex.text)
	}
	

	# Plot lines and points
	# Build plot area
	par(fig = c(0.3, 0.7, 0, 1), new = TRUE)
	par(mar = c(5.1, 0.1, 4.1, 0.1))
	
	# Build plot area
	plot(c(0,0),col=0,xlim=c(x.min,x.max),ylim=c(y.min,y.max),axes=FALSE,ylab="",xlab=paste("Log odds ratio vs",treatment.names[[1]][1]))
	axis(side=1)	
	title(paste(outcome.names[outcome.i],"connected RCTs only"), cex.main = cex.text)

	#text(x=x.min,y=2.25+2*length(disconnect.only.treats),"Disconnected from reference",pos=4,cex=cex.text)
	#lines(x=c(x.min,x.max),y=rep(1.5+2*length(disconnect.only.treats),2))
	#text(x=x.min,y=3.25+2*(length(disconnect.only.treats)+length(rct.only.treats)),"Connected to reference",pos=4,cex=cex.text)
	#lines(x=c(x.min,x.max),y=rep(2.5+2*(length(disconnect.only.treats)+length(rct.only.treats)),2))

	# Plot the connected to reference treatment effects
	for(i.treat in rct.only.treats[-1])
	{
	
		# All RCTs ('true' effects)
		points(y=2*(which(rct.only.treats==i.treat)),x=bugs.object.disc.rct.only.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.disc.rct.only.re[[outcome.i]]$summary)),1][i.treat-1])
		lines(y=rep(2*(which(rct.only.treats==i.treat)),2),x=bugs.object.disc.rct.only.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.disc.rct.only.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,])

		# RE on baseline
		points(y=-0.5+2*(which(rct.only.treats==i.treat)),x=bugs.object.rebase.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.disc.re[[outcome.i]]$summary)),1][i.treat-1],col=2)
		lines(y=rep(-0.5+2*(which(rct.only.treats==i.treat)),2),x=bugs.object.rebase.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.disc.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,],col=2)

		# FE on baseline but with CHADS2 as covariate
		points(y=-1+2*(which(rct.only.treats==i.treat)),x=bugs.object.plugin.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.disc.re[[outcome.i]]$summary)),1][i.treat-1],col=3)
		lines(y=rep(-1+2*(which(rct.only.treats==i.treat)),2),x=bugs.object.plugin.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.disc.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,],col=3)
	}
	legend("topright",legend=c("RCT only","Reference prediction","Aggregate level matching"),lty=1,col=c(1:3),cex=cex.text)
	
	
	
	par(fig = c(0.7, 1, 0, 1), new = TRUE)
	
	par(mar=c(5.1,1,4.1,1))
	text.plot.x <- c(0,6)
	plot(c(0, 0), xlab = "", col = 0, xlim = text.plot.x, ylim = c(y.min, y.max), ylab = "",
	     axes=FALSE)
	cex.effect.estimate <- 0.45
	text(x=text.plot.x[2]/6, y = y.max, pos = 1, "True (all RCTs)", cex = cex.effect.estimate, col = 1)
	text(x = 3*text.plot.x[2]/6, y = y.max, pos = 1, "RP", cex = cex.effect.estimate, col = 2)
	text(x = 5*text.plot.x[2]/6, y = y.max, pos = 1,"ALM", cex = cex.effect.estimate, col = 3)
	for(i.treat in rct.only.treats[-1]) {
	  text(y=2*(which(rct.only.treats==i.treat)), x = 0, format.effect.summary(bugs.object.disc.rct.only.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.disc.rct.only.re[[outcome.i]]$summary)),c(1,3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 1)
	  text(y=2*(which(rct.only.treats==i.treat)), x = text.plot.x[2]/3, format.effect.summary(bugs.object.rebase.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.disc.re[[outcome.i]]$summary)),c(1,3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 2)
	  text(y=2*(which(rct.only.treats==i.treat)), x = 2*text.plot.x[2]/3, format.effect.summary(bugs.object.plugin.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.disc.re[[outcome.i]]$summary)),c(1,3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 3)
	}
	dev.off()
	# Reset parameters
	par(mar = old.mar)
	par(fig = old.fig)
	par(mfrow = old.mfrow)
	options(scipen = old.scipen)
	
	# Open a pdf file for output
	pdf(file=paste(baseline.directory,"/results/forest plots/v2.disconnected.disconnected",outcome.names[outcome.i],".forest.re.pdf",sep=""))

	
	# This depends on what is being plotted
	if(outcome.i!=4){cex.text<-.9}
	if(outcome.i==4){cex.text<-0.65}
	
	
	y.min=0; y.max=2*length(disconnect.only.treats)+1
	
	
	# Store the old plot parameters
	old.mar <- par()$mar
	old.fig <- par()$fig
	old.mfrow <- par()$mfrow
	old.scipen <- options()$scipen
	# Turn off scientific notation
	options(scipen = 999)

	# Create 3 plots
	# One for treatment names, one for the forest plot, and one for effect summaries
	par(mfrow = c(1,3))
	par(fig = c(0, 0.3, 0, 1))
	par(mar=c(5.1,1,4.1,0.1))
	
	
	# Plot treatment names
	plot(c(0,0), xlab = "", col = 0, xlim = c(0, 4), ylim = c(y.min, y.max), ylab = "", axes = FALSE)
	
	# Plot the RCT only treatment effects
	for(i.treat in disconnect.only.treats)
	{
	  text(treatment.names[[outcome.i]][i.treat],x=0,y=2*(which(disconnect.only.treats==i.treat)),pos=4,cex=cex.text)
	}
	
	# Plot lines and points
	# Build plot area
	par(fig = c(0.3, 0.7, 0, 1), new = TRUE)
	par(mar = c(5.1, 0.1, 4.1, 0.1))
	
	plot(c(0,0),col=0,xlim=c(x.min,x.max),ylim=c(y.min,y.max),axes=FALSE,ylab="",xlab=paste("Log odds ratio vs",treatment.names[[1]][1]))
	axis(side=1)	
	title(paste(outcome.names[outcome.i],"disconnected RCTs"), cex.main = cex.text)


	# Plot the disconnected from reference treatment effects
	for(i.treat in disconnect.only.treats)
	{
		# All RCTs ('true' effects)
		points(y=2*(which(disconnect.only.treats==i.treat)),x=bugs.object.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.re[[outcome.i]]$summary)),1][i.treat-1])
		lines(y=rep(2*(which(disconnect.only.treats==i.treat)),2),x=bugs.object.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,])

		# RE on baseline
		points(y=-0.5+2*(which(disconnect.only.treats==i.treat)),x=bugs.object.rebase.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.disc.re[[outcome.i]]$summary)),1][i.treat-1],col=2)
		lines(y=rep(-0.5+2*(which(disconnect.only.treats==i.treat)),2),x=bugs.object.rebase.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.disc.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,],col=2)

		# FE on baseline but with CHADS2 as covariate
		points(y=-1+2*(which(disconnect.only.treats==i.treat)),x=bugs.object.plugin.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.disc.re[[outcome.i]]$summary)),1][i.treat-1],col=3)
		lines(y=rep(-1+2*(which(disconnect.only.treats==i.treat)),2),x=bugs.object.plugin.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.disc.re[[outcome.i]]$summary)),c(3,7)][i.treat-1,],col=3)
	}
	legend("topright",legend=c("True (all RCTs)","Reference prediction","Aggregate level matching"),lty=1,col=c(1,2,3),cex=cex.text)
	
	par(fig = c(0.7, 1, 0, 1), new = TRUE)
	
	par(mar=c(5.1,1,4.1,1))
	text.plot.x <- c(0,6)
	plot(c(0, 0), xlab = "", col = 0, xlim = text.plot.x, ylim = c(y.min, y.max), ylab = "",
	     axes=FALSE)
	cex.effect.estimate <- 0.45
	text(x=text.plot.x[2]/6, y = y.max, pos = 1, "True (all RCTs)", cex = cex.effect.estimate, col = 1)
	text(x = 3*text.plot.x[2]/6, y = y.max, pos = 1, "RP", cex = cex.effect.estimate, col = 2)
	text(x = 5*text.plot.x[2]/6, y = y.max, pos = 1,"ALM", cex = cex.effect.estimate, col = 3)
	for(i.treat in disconnect.only.treats) {
	  text(y=2*(which(disconnect.only.treats==i.treat)), x = 0, format.effect.summary(bugs.object.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.re[[outcome.i]]$summary)),c(1,3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 1)
	  text(y=2*(which(disconnect.only.treats==i.treat)), x = text.plot.x[2]/3, format.effect.summary(bugs.object.rebase.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.rebase.disc.re[[outcome.i]]$summary)),c(1,3, 7)][i.treat-1]), pos = 4, cex = cex.effect.estimate, col = 2)
	  text(y=2*(which(disconnect.only.treats==i.treat)), x = 2*text.plot.x[2]/3, format.effect.summary(bugs.object.plugin.disc.re[[outcome.i]]$summary[grep("d",rownames(bugs.object.plugin.disc.re[[outcome.i]]$summary)),c(1,3,7)][i.treat-1,]), pos = 4, cex = cex.effect.estimate, col = 3)
	}
	
	dev.off()
	# Reset parameters
	par(mar = old.mar)
	par(fig = old.fig)
	par(mfrow = old.mfrow)
	options(scipen = old.scipen)
	
}





