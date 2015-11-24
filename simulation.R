### Introduction
# This demo reproduces the ordered hypothesis testing simulation in the paper:
# Ang Li and Rina Foygel Barber, "Accumulation tests for FDR control in ordered hypothesis testing" (2015). Available from http://arxiv.org/abs/1505.07352
# The demo was run using R version 3.2.0.

# The demo considers an ordered hypothesis testing problem with n hypotheses, of which n_1 contain true signals. We produce an ordered list of p-values as follows:
# * First, generate (Z_prior)_i ~ N(mu_i,1) where mu_i=0 for index i corresponding to a null hypothesis, and mu_i=mu_intermix for a nonnull (true signal) index i. This z-score represents prior information, e.g. data from a prior experiment.
# * Then, reorder the hypotheses according to their z-scores, with the largest z-scores (in absolute value) first. If mu_intermix is large, then this yields good separation: the n_1 true signals will mostly be placed early in the list. If mu_intermix is small, however, then this yields poor separation: many true signals will be mixed in with the nulls throughout the list.
# * Next, for each hypothesis, generate a new z-score Z_i~N(mu_i,1) where now mu_i=0 or mu_i=mu_signal, for either a strong signal strength (mu_signal large) or a weak signal strength (mu_signal small). Produce p-values p_i with a two-sided z-test. See paper for details.

# We then run accumulation tests, specifically, we text the HingeExp method (with parameter C=2), the SeqStep and SeqStep+ (Barber and Candès 2014) (each with parameters C=2), and the ForwardStop method (G'Sell et al 2013).

### Setup
# We first define some functions for running the simulation and for plotting the results.
Setup_simulation_functions = function(){
  # PSeq: returns the simulated ordered hypotheses -- the P-value sequence, and labels of nonnulls and nulls
  # n (num of predictors), pr (prop of Hyps that is nonnulls), mu1 (controls the separation of nonnulls and nulls), mu2 (controls signal strength)
  PSeq <<- function(n, pr, mu1, mu2, seed){
    set.seed(seed)
    pr=floor(n*pr)/n # so that n*pr is a whole number
    z = rnorm(n*pr, mean = mu1, sd = 1)  # determines the intermix
    NonNullP = 2*(1-pnorm(abs(z), mean = 0, sd = 1))
    NullP = runif(n*(1-pr), min = 0, max = 1)
    PriorP=cbind(c(NonNullP,NullP),c(rep(1,n*pr),rep(0,n*(1-pr))))
    NonNullIndex=sort(rank(c(NonNullP,NullP),ties.method='random')[1:(n*pr)])
    NullIndex = setdiff(seq(1, n, 1), NonNullIndex)
    zNew = rnorm(n*pr, mean = mu2, sd = 1) # determines the signal strength
    NewNonNullP = 2*(1-pnorm(abs(zNew), mean = 0, sd = 1))
    NewNullP = runif(n*(1-pr))
    PSeq = rep(0, n)
    PSeq[NonNullIndex] = NewNonNullP   # nonNullP
    PSeq[NullIndex] = NewNullP
    return(cbind(PSeq, is.element(1:n,NonNullIndex)))
  }
  
 # CompEffectFunc: returns the average power and observed FDR for accumulation test with different h() functions
  CompEffectFunc <<- function(n, pr, mu1, mu2, hfuns, numerator_pluses,denominator_pluses, seeds, set_alpha){
  numSimu=length(seeds)
    Power = matrix(0,length(set_alpha),length(hfuns)) 
	  FDR = matrix(0,length(set_alpha),length(hfuns)) 
    for (i in 1:numSimu){
		  Pvals = PSeq(n, pr, mu1, mu2, seeds[i])     
		  for(j in 1:length(hfuns)){
		  	kh = AccumulationTest(Pvals[,1],hfuns[[j]],alpha=set_alpha,numerator_plus=numerator_pluses[j],denominator_plus=denominator_pluses[j],output_type='khat')
	  		Power[,j] = Power[,j] + unlist(lapply(kh, function(x) sum(Pvals[seq_len(x),2])/(n*pr))) # Pvals[,2] is an indicator of nonnulls
  			FDR[,j] = FDR[,j] + unlist(lapply(kh, function(x) sum(1-Pvals[seq_len(x),2])/max(x, 1))) 
		  }	
	  }
    Power=Power/numSimu
    FDR=FDR/numSimu
    return(list(Power, FDR))
  }

# FDPSeqFunc: returns the estimated FDP by h() functions, and the actual FDP
  FDPSeqFunc <<- function(n,pr,mu1,mu2,hfuns,numerator_pluses, denominator_pluses,seeds){
	  numSimu=length(seeds)
	  FDPtrue=rep(0,n)
	  FDPhat=matrix(0,n,length(hfuns))
	  for(i in 1:numSimu){
      Pvals = PSeq(n, pr, mu1, mu2, seeds[i])     
		  for(j in 1:length(hfuns)){
	      FDPhat[,j]=FDPhat[,j]+AccumulationTest(Pvals[,1],hfuns[[j]],numerator_plus=numerator_pluses[j],denominator_plus=denominator_pluses[j],output_type='FDPest')
      }
		  FDPtrue=FDPtrue+cumsum(1-Pvals[,2])/(1:n)
	  }
	  FDPhat=FDPhat/numSimu;FDPtrue=FDPtrue/numSimu
	  return(list(FDPhat,FDPtrue))
  }
  print('Simulation functions defined.')
}

Setup_simulation_plotting_functions = function(){
    FDR_and_power_plot <<- function(FDRmat, powermat, names, cols, pchs, alpha, alpha_display_div=1, title='',leg_coords=integer(0),show_legend=FALSE){
    alpha_display=alpha[seq(1,length(alpha),by=alpha_display_div)]
	  if(length(leg_coords)==0){leg_x=min(alpha);leg_y=1}else{leg_x=leg_coords[1];leg_y=leg_coords[2]}
	  xzero=min(alpha)-.05*(max(alpha)-min(alpha))
	  xzero1=min(alpha)-.17*(max(alpha)-min(alpha))
	  xzero2=min(alpha)-.07*(max(alpha)-min(alpha))
	  xzero3=min(alpha)-.1*(max(alpha)-min(alpha))
	  plot(0:1,0:1,type='n',xlab='Target FDR',ylab='',xlim=c(xzero1,max(alpha)),ylim=c(-.55,1.05),axes=FALSE,main=title,cex.main=2,cex.lab=1.5)
	  segments(xzero2,0,max(alpha),0)
	  axis(side=1,at=alpha_display,cex.axis=1.5)
  	segments(xzero,-.55,xzero,-.05)
	  segments(xzero,.05,xzero,1.05)
	  for(i in 0:5){
		  segments(xzero,-.55+i/10,xzero2,-.55+i/10)
		  text(xzero3,-.55+i/10,i/10,cex=1.5)
	  }
	  for(i in 0:10){
		  segments(xzero,.05+i/10,xzero2,0.05+i/10)
		  text(xzero3,.05+i/10,i/10,cex=1.5)
	  }
	  text(xzero1,-.3,'Avg. observed FDR',srt=90,cex=1.5)
	  text(xzero1,.55,'Average power',srt=90,cex=1.5)
	  for(i in 1:length(names)){
		  points(alpha,FDRmat[,i]-.55,type='l',col=cols[i])
	  	points(alpha,powermat[,i]+.05,type='l',col=cols[i])		
  		points(alpha,FDRmat[,i]-.55,pch=pchs[i],col=cols[i],cex=1.5)
		  points(alpha,powermat[,i]+.05,pch=pchs[i],col=cols[i],cex=1.5)		
	  }
	  points(alpha,alpha-.55,type='l',lty='dotted',col='gray50',lwd=2)
	  if(show_legend){
		  legend(leg_x,leg_y,legend=names,col=cols,pch=pchs,lty='solid',cex=1.5)
	  }
  }

  FDP_vs_k_plot <<- function(FDPhat_mat, FDPtrue, names, cols, pchs, title='',xlab='k',ylab='Estimated FDP(k)', kmax=0,show_legend=FALSE,leg_coords=integer(0)){
    n=length(FDPtrue)
  	if(length(leg_coords)==0){leg_x=kmax/10;leg_y=1}else{leg_x=leg_coords[1];leg_y=leg_coords[2]}
    if (kmax==0){kmax=n}
	  plot(0:1,0:1,type='n',xlim=c(0,kmax),ylim=c(0,1),main=title,xlab=xlab,ylab=ylab,cex.main=2,cex.axis=1.5,cex.lab=1.5)
	  points(1:kmax,FDPtrue[1:kmax],type='l',col='gray50',lty='dotted', lwd=2)
	  for(i in 1:length(names)){
		  points(1:kmax,FDPhat_mat[1:kmax,i],type='l',col=cols[i])
		  points((1:kmax)[(1:(floor(kmax)/20))*20],FDPhat_mat[(1:kmax)[(1:(floor(kmax)/20))*20],i],col=cols[i],pch=pchs[i],cex=1.5)
	  }
	  if(show_legend){
		  legend(leg_x,leg_y,legend=c(names,'True FDP'),col=c(cols,'gray50'),pch=c(pchs,0),lty=c(rep('solid',length(names)),'dotted'),lwd=c(rep(1,length(names)),2),cex=1.5,pt.cex=c(rep(1.5,length(names)),0))
	  }
  }
  
  get_titles <<- function(mu_low,mu_high){
    matrix(c(substitute(paste('Poor separation & weak signal (',mu[1]==mu1,', ',mu[2]==mu2,')'),list(mu1=toString(mu_low),mu2=toString(mu_low))),
    substitute(paste('Good separation & weak signal (',mu[1]==mu1,', ',mu[2]==mu2,')'),list(mu1=toString(mu_high),mu2=toString(mu_low))),
    substitute(paste('Poor separation & strong signal (',mu[1]==mu1,', ',mu[2]==mu2,')'),list(mu1=toString(mu_low),mu2=toString(mu_high))),
    substitute(paste('Good separation & strong signal (',mu[1]==mu1,', ',mu[2]==mu2,')'),list(mu1=toString(mu_high),mu2=toString(mu_high)))),2,2)
  }
  print('Simulation plotting functions defined.')
}

Setup_simulation_functions()
Setup_simulation_plotting_functions()
source('accumulation_test_functions.R')


### Parameters for the simulations:
hfuns=c(create_SeqStep_function(C=2),create_SeqStep_function(C=2),create_HingeExp_function(C=2),create_ForwardStop_function())
numerator_pluses=c(2,0,0,0)
denominator_pluses=c(1,0,0,0)
names=c('SeqStep+ (C=2)','SeqStep (C=2)','HingeExp (C=2)','ForwardStop')
mu_low=2;mu_high=3; 
n=1000; ntrials=100; 
pr_list=c(0.2,0.1) # proportion of true signals
alphaseq=(2:10)/40
seeds=1:ntrials


### Simulations part 1
# Settings: n=1000 hypotheses, k*=200 true signals

# First, we run each method with target FDR level alpha (ranging in alpha={0.05,0.075,0.1,...,0.25}). For each method and each alpha, we record the actual FDR attained, and the power to detect signals, averaged over 100 trials.
pr=pr_list[1]
#pdf(paste0('FDR_vs_power_plot_n_',n,'_kstar_',n*pr,'.pdf'),14,14)
par(mfrow=c(2,2))
for(i in 1:2){for(j in 1:2){
  mu1=c(mu_low,mu_high)[i];mu2=c(mu_low,mu_high)[j]
	temp=CompEffectFunc(n,pr,mu1,mu2,hfuns,numerator_pluses,denominator_pluses,seeds, alphaseq)
FDR_and_power_plot(temp[[2]],temp[[1]],names,cols=c('black','blue','purple','red'),pchs=15:18,alpha=alphaseq, alpha_display_div=2, title=get_titles(mu_low,mu_high)[i,j],show_legend=(i==1 && j==1))
}};#dev.off()

# Now we look at how each method estimates the FDP, along the sequence of p-values. At each k, we compare the true FDP among the first k p-values, with the estimated FDP at k for each method. (All results averaged over 100 trials).
#pdf(paste0('FDP_vs_k_plot_n_',n,'_kstar_',n*pr,'.pdf'),14,14)
par(mfrow=c(2,2))
for(i in 1:2){for(j in 1:2){
	mu1=c(mu_low,mu_high)[i];mu2=c(mu_low,mu_high)[j]
	temp=FDPSeqFunc(n,pr,mu1,mu2,hfuns,numerator_pluses,denominator_pluses,seeds)
	FDP_vs_k_plot(temp[[1]],temp[[2]],names,cols=c('black','blue','purple','red'),pchs=15:18,title=get_titles(mu_low,mu_high)[i,j],kmax=300,show_legend=(i==1 && j==1))
}};#dev.off()


### Simulations part 2
# Settings: n=1000 hypotheses, k*=100 true signals

# First, we run each method with target FDR level alpha (ranging in alpha={0.05,0.075,0.1,...,0.25}). For each method and each alpha, we record the actual FDR attained, and the power to detect signals, averaged over 100 trials.
pr=pr_list[2]
#pdf(paste0('FDR_vs_power_plot_n_',n,'_kstar_',n*pr,'.pdf'),14,14)
par(mfrow=c(2,2))
for(i in 1:2){for(j in 1:2){
  mu1=c(mu_low,mu_high)[i];mu2=c(mu_low,mu_high)[j]
  temp=CompEffectFunc(n,pr,mu1,mu2,hfuns,numerator_pluses,denominator_pluses,seeds, alphaseq)
FDR_and_power_plot(temp[[2]],temp[[1]],names,cols=c('black','blue','purple','red'),pchs=15:18,alpha=alphaseq, alpha_display_div=2, title=get_titles(mu_low,mu_high)[i,j],show_legend=(i==1 && j==1))
}};#dev.off()

# Now we look at how each method estimates the FDP, along the sequence of p-values. At each k, we compare the true FDP among the first k p-values, with the estimated FDP at k for each method. (All results averaged over `r ntrials` trials).
#pdf(paste0('FDP_vs_k_plot_n_',n,'_kstar_',n*pr,'.pdf'),14,14)
par(mfrow=c(2,2))
for(i in 1:2){for(j in 1:2){
  mu1=c(mu_low,mu_high)[i];mu2=c(mu_low,mu_high)[j]
	temp=FDPSeqFunc(n,pr,mu1,mu2,hfuns,numerator_pluses,denominator_pluses,seeds)
	FDP_vs_k_plot(temp[[1]],temp[[2]],names,cols=c('black','blue','purple','red'),pchs=15:18,title=get_titles(mu_low,mu_high)[i,j],kmax=300,show_legend=(i==1 && j==1))
}};#dev.off()


### References 

# Barber, Rina Foygel, and Emmanuel Candès. 2014. “Controlling the False Discovery Rate via Knockoffs.” ArXiv Preprint ArXiv:1404.5609.

# G’Sell, Max Grazier, Stefan Wager, Alexandra Chouldechova, and Robert Tibshirani. 2013. “False Discovery Rate Control for Sequential Selection Procedures, with Application to the Lasso.” ArXiv Preprint ArXiv:1309.5352.