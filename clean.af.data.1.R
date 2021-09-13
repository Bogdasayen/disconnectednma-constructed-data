# Script to remove NA arms from the AF data


# Build a list of treatments just included for each outcome
treatment.names<-list()

# Clean the data so that NA arms are removed
for(outcome.i in 1:length(outcomes.to.save))
{
	# Remove trials with NA or 0 in all arms 
	trials.to.keep<-rowSums(!is.na(bugs.data[[outcome.i]]$r))>0 & rowSums(bugs.data[[outcome.i]]$r,na.rm=TRUE)>0
	bugs.data[[outcome.i]]$r<-bugs.data[[outcome.i]]$r[trials.to.keep,]
	bugs.data[[outcome.i]]$n<-bugs.data[[outcome.i]]$n[trials.to.keep,]
	bugs.data[[outcome.i]]$t<-bugs.data[[outcome.i]]$t[trials.to.keep,]
	bugs.data[[outcome.i]]$x<-bugs.data[[outcome.i]]$x[,trials.to.keep,]
	bugs.data[[outcome.i]]$na<-bugs.data[[outcome.i]]$na[trials.to.keep]
	bugs.data[[outcome.i]]$ns<-sum(trials.to.keep)
	bugs.data[[outcome.i]]$nt<-length(unique(c(bugs.data[[outcome.i]]$t[!is.na(bugs.data[[outcome.i]]$nt)])))-1 # Minus 1 for the NA
	
	# Apply a continutity correction (this was done in main DOAC NMA)
	bugs.data[[outcome.i]]$r[bugs.data[[outcome.i]]$r==0]<-0.5

	# Fix treatments to be consecutive
	tr<-bugs.data[[outcome.i]]$t
	t.pre<-tr # For debugging
	t.list.premap<-unique(c(tr[!is.na(tr)]))
	
	# Map the treatment names to consecutive integers
	tr<-tr+10000
	t.nonconsec<-unique(c(tr[!is.na(tr)]))
	for(i in 1:bugs.data[[outcome.i]]$nt)
	{
		tr[tr==t.nonconsec[i]]<-i
	}
	# Save the renamed treatments
	treatment.names[[outcome.i]]<-t.names[t.list.premap]

	# Make sure the treatment are listed in numerical order
	for(i in 1:bugs.data[[outcome.i]]$ns)
	{
		bugs.data[[outcome.i]]$n[i,]<-bugs.data[[outcome.i]]$n[i,order(tr[i,])]
		bugs.data[[outcome.i]]$r[i,]<-bugs.data[[outcome.i]]$r[i,order(tr[i,])]
		bugs.data[[outcome.i]]$x[,i,]<-bugs.data[[outcome.i]]$x[,i,order(tr[i,])]
		tr[i,]<-tr[i,order(tr[i,])]
	}
	bugs.data[[outcome.i]]$t<-tr
}