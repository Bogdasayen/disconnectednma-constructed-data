# Code provided by Joy Leahy October 2017
# RCtMultiplier and MatchedMultiplier dictate the amount of weight each type of evidence receives. RCTMultiplier usually stays at 1, but you can down-weight the matched by choosing a lower number. 
# For each treatment in the vector "RCTorMatched" put in 1 if it's only in a matched study, 2 if it's in both and 3 if it's only in an RCT.
# Treatment 1 should be in both.

model.aggregate.level.matching<-function(){                   # *** PROGRAM STARTS

### Overall Level ###

d[1]<-0
RCTd[1]<-0
MATCHEDd[1]<-0

for (k in 2:COnt){
RCTd[k]~dnorm(d[k],COtauRCT)
MATCHEDd[k]~dnorm(d[k],COtauMATCHED)
d[k]~dnorm(0,0.298)}



for (k in 1:COnt){
d_branch[k, 3]<-RCTd[k]
d_branch[k, 1]<-MATCHEDd[k]
d_branch[k, 2]<-d[k]
Used_d[k]<-d_branch[k, RCTorMatched[k]] }

	
COsd~dunif(0,2)
COvar<-pow(COsd,2)
#COtau<-1/COvar
COtauMATCHED<-MatchedMultiplier/COvar
COtauRCT<-RCTMultiplier/COvar

for (c in 1:(COnt-1)) {
for (k in (c+1):COnt) {
or[k,c] <- exp(Used_d[k] - Used_d[c])
lor[k,c] <- (Used_d[k]-Used_d[c])
}
}
	

### RCT bit ###

for(i in 1:RCTns){          # LOOP THROUGH STUDIES
   RCTw[i,1] <- 0  # adjustment for multi-arm trials is zero for control arm
   RCTdelta[i,1] <- 0      # treatment effect is zero for control arm
   RCTmu[i] ~ dnorm(0,0.298)       # vague priors for all trial baselines
   for (k in 1:RCTna[i]) {     # LOOP THROUGH ARMS
      RCTr[i,k] ~ dbin(RCTp[i,k],RCTn[i,k]) # binomial likelihood
      logit(RCTp[i,k]) <- RCTmu[i] + RCTdelta[i,k] # model for linear predictor
}

  for (k in 2:RCTna[i]) {     # LOOP THROUGH ARMS
# trial-specific LOR distributions
      RCTdelta[i,k] ~ dnorm(RCTmd[i,k],RCTtaud[i,k])
	
# mean of LOR distributions (with multi-arm trial correction)
      RCTmd[i,k] <- RCTd[RCTt[i,k]] - RCTd[RCTt[i,1]]+ RCTsw[i,k]
# precision of LOR distributions (with multi-arm trial correction)
	  RCTtaud[i,k] <- RCTtau*2*(k-1)/k
# adjustment for multi-arm s
      RCTw[i,k] <- (RCTdelta[i,k] - RCTd[RCTt[i,k]] + RCTd[RCTt[i,1]])
# cumulative adjustment for multi-arm trials
      RCTsw[i,k] <- sum(RCTw[i,1:k-1])/(k-1)
   }
 }
# vague priors for treatment effects
RCTsd ~ dunif(0,2) # vague prior for between-trial SD
RCTtau <- pow(RCTsd,-2) # between-trial precision = (1/between-trial variance)



### MATCHED bit ###

for(i in 1:MATCHEDns){            # LOOP THROUGH STUDIES
   MATCHEDw[i,1] <- 0  # adjustment for multi-arm trials is zero for control arm
   MATCHEDdelta[i,1] <- 0      # treatment effect is zero for control arm
   MATCHEDmu[i] ~ dnorm(0,0.298)       # vague priors for all trial baselines
   for (k in 1:MATCHEDna[i]) {     # LOOP THROUGH ARMS
      MATCHEDr[i,k] ~ dbin(MATCHEDp[i,k],MATCHEDn[i,k]) # binomial likelihood
      logit(MATCHEDp[i,k]) <- MATCHEDmu[i] + MATCHEDdelta[i,k] # model for linear predictor
}

  for (k in 2:MATCHEDna[i]) {     # LOOP THROUGH ARMS
# trial-specific LOR distributions
      MATCHEDdelta[i,k] ~ dnorm(MATCHEDmd[i,k],MATCHEDtaud[i,k])
	
# mean of LOR distributions (with multi-arm trial correction)
      MATCHEDmd[i,k] <- MATCHEDd[MATCHEDt[i,k]] - MATCHEDd[MATCHEDt[i,1]]+ MATCHEDsw[i,k]
# precision of LOR distributions (with multi-arm trial correction)
	  MATCHEDtaud[i,k] <- MATCHEDtau*2*(k-1)/k
# adjustment for multi-arm s
      MATCHEDw[i,k] <- (MATCHEDdelta[i,k] - MATCHEDd[MATCHEDt[i,k]] + MATCHEDd[MATCHEDt[i,1]])
# cumulative adjustment for multi-arm trials
      MATCHEDsw[i,k] <- sum(MATCHEDw[i,1:k-1])/(k-1)
   }
 }
# vague priors for treatment effects
MATCHEDsd ~ dunif(0,2) # vague prior for between-trial SD
MATCHEDtau <- pow(MATCHEDsd,-2) # between-trial precision = (1/between-trial variance)




}
# *** PROGRAM ENDS











