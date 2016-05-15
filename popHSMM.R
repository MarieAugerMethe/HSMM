#!#
#########################################################
# Based on the code from Langrock et al. 2012
# Fit the different models to the bear population
# M0: All bears are the same (all parameters are the same accros diet)
# M1: Bears of different diet differ in their Step length (scale of Weibull)
  # M1I: Differ in intensive behaviour (potentially foraging)
  # M1E: Differ in extensive behaviour (potentially travelling)
  # M1B: Differ in both intensive and extensive behaviours
# M2: Bears of different diet differ in their the number of steps in a behaviour (poisson lambda)
  # M2I: Differ in intensive behaviour (potentially foraging)
  # M2E: Differ in extensive behaviour (potentially travelling)
  # M2B: Differ in both intensive and extensive behaviours
# M3: Bears of different diets differ in how correlated their movement is (rho from wrapped Cauchy)
  # M3I: Differ in intensive behaviour (potentially foraging)
  # M3E: Differ in extensive behaviour (potentially travelling)
  # M3B: Differ in both intensive and extensive behaviours
# M4: Bears of different diets differ in their mean turning angle (mu from wrapped Cauchy)
  # M4I: Differ in intensive behaviour (potentially foraging)
  # M4E: Differ in extensive behaviour (potentially travelling)
  # M4B: Differ in both intensive and extensive behaviours

#####
# Load packages
library(CircStats) # for wrapped Cauchy

#####
# Get data for the bears

##
# Movement data
# observations ("OBS") need to be given in an n x (2*p) matrix, 
# with p being the number of individuals, 
# and the observed step lengths and turning angles for individual i 
# given in columns 2*(i-1)+1 and 2*(i-1)+2, respectively
OBS <- read.csv("beardata.csv")

##
# Covariate data - diet!
bearsCap <- read.csv("grizzlyIs.csv")

# Same order
OBS[1:2,1:8]
bearsCap[1:4,]

#####################################################
# Define functions common to all models

###
# function that derives the t.p.m. of the HMM that represents the HSMM
# using poisson instead of neg binom
genGamma <- function(m,lamb){
  Gamma <- diag(m[1]+m[2])*0
  probs1 <- dpois(0:(m[1]-1),lambda=lamb[1])
  probs2 <- dpois(0:(m[2]-1),lambda=lamb[2])
  
  # Denominator of c(r): 1 - sum_{k=1}^{r-1}p(k), so for r=1 -> 0 (because empty sum equals 0)
  # Use cumulative distribution function because it is the sum of the prob
  den1 <- 1 - c(0,ppois(0:(m[1]-2),lambda=lamb[1]))
  den2 <- 1 - c(0,ppois(0:(m[2]-2),lambda=lamb[2]))
  
  # To remove the chance of getting Inf
  probs1[which(den1<1e-12)] <- 1
  den1[which(den1<1e-12)] <- 1
  probs2[which(den2<1e-12)] <- 1
  den2[which(den2<1e-12)] <- 1
  
  # state aggregate 1
  Gamma[1:m[1],m[1]+1] <- probs1/den1 # c_1(r) for r=1,2,...,N_1* in first column of Beh 2
  diag(Gamma[1:(m[1]-1),2:m[1]]) <- 1-Gamma[1:(m[1]-1),m[1]+1] # 1-c_1(r), for r=1,2,...,N_1*-1
  Gamma[m[1],m[1]] <- 1 - Gamma[m[1],m[1]+1] # 1-c_1(N_1*)
  
  # state aggregate 2
  Gamma[m[1]+(1:m[2]),1] <- probs2/den2 # c_2(r) for r=1,2,...,N_2* in first column of Beh 1
  diag(Gamma[m[1]+1:(m[2]-1),m[1]+2:m[2]]) <- 1 - Gamma[m[1]+1:(m[2]-1),1] # 1-c_2(r), for r=1,2,...,N_2*-1
  Gamma[m[1]+m[2],m[1]+m[2]] <- 1 - Gamma[m[1]+m[2],1] # 1-c_2(N_2*)
  return(Gamma)
}

## function that transforms each of the (possibly constrained) parameters to the real line
move.HSMM.pn2pw <- function(a,b,kapp,gamLam,co){
  ta <- log(a) # Weibull scale
  tb <- log(b)
  tkappa <- logit(kapp)
  tgamLam <- log(gamLam) # Mean number of step poisson
  #tco <- c(log((co[1])/(2*pi-co[1])),log((pi+co[2])/(pi-co[2])))
  tco <- co
  parvect <- c(ta,tb,tkappa,tgamLam,tco)
  return(parvect)
}

## inverse transformation back to the natural parameter space            
move.HSMM.pw2pn <- function(parvect,M){
  epar <- exp(parvect)
  aE <- 2
  rE <- 2
  gE <- 2
  #"M0"
  a <- rep(epar[1:2], each=nDiet)
  if(M == "M1B"){
    aE <- 2*nDiet
    a <- epar[1:aE]
  }
  if(M == "M1I"){
    aE <- nDiet + 1
    a <- c(epar[1:nDiet], rep(epar[aE], nDiet))
  }
  if(M == "M1E"){
    aE <- nDiet + 1
    a <- c(rep(epar[1], nDiet), epar[2:aE])
  }
  b <- epar[aE+(1:2)]

  kapp <- rep(inv.logit(parvect[aE+(3:4)]),each=nDiet)
  if(M == "M3B"){
    rE <- 2*nDiet
    kapp <- epar[aE+(1:rE)]
  }
  if(M == "M3I"){
    rE <- nDiet + 1
    kapp <- c(epar[aE+(1:nDiet)], rep(epar[aE+rE], nDiet))
  }
  if(M == "M3E"){
    rE <- nDiet + 1
    kapp <- c(rep(epar[aE+1], nDiet), epar[aE+(2:rE)])
  }

  gamLam <- c(rep(epar[aE+rE+3],nDiet), rep(epar[aE+rE+4],nDiet))
  if(M == "M2B"){
    gE <- 2*nDiet
    gamLam <- epar[aE+rE+(2+(1:gE))]
  }
  if(M == "M2I"){
    gE <- nDiet + 1
    gamLam <- c(epar[aE+rE+(2+(1:nDiet))], rep(epar[aE+rE+(2+(gE))],nDiet)) 
  }
  if(M == "M2E"){
    gE <- nDiet + 1
    gamLam <- c(rep(epar[aE+rE+3],nDiet), epar[aE+rE+(2+(2:gE))])
  }
  
#   co <- c(rep(2*pi*inv.logit(parvect[aE+rE+gE+3]),nDiet),
#           rep(pi*(exp(parvect[aE+rE+gE+4])-1)/(exp(parvect[aE+rE+gE+4])+1),nDiet))
#   if(M == "M4B"){
#     co <- c(2*pi*inv.logit(parvect[aE+rE+gE+2+(1:nDiet)]),
#             pi*(exp(parvect[aE+rE+gE+2+nDiet+(1:nDiet)])-1)/(exp(parvect[aE+rE+gE+2+nDiet+(1:nDiet)])+1))
#   }
#   if(M == "M4I"){
#     co <- c(2*pi*inv.logit(parvect[aE+rE+gE+2+(1:nDiet)]),
#             rep(pi*(exp(parvect[aE+rE+gE+3+nDiet])-1)/(exp(parvect[aE+rE+gE+3+nDiet])+1),nDiet))
#   }
#   if(M == "M4E"){
#     co <- c(rep(2*pi*inv.logit(parvect[aE+rE+gE+3]),nDiet),
#             pi*(exp(parvect[aE+rE+gE+3+(1:nDiet)])-1)/(exp(parvect[aE+rE+gE+3+(1:nDiet)])+1))
#   }
  
  co <- c(rep(parvect[aE+rE+gE+3],nDiet),
          rep(parvect[aE+rE+gE+4],nDiet))
  if(M == "M4B"){
    co <- c(parvect[aE+rE+gE+2+(1:nDiet)],
            parvect[aE+rE+gE+2+nDiet+(1:nDiet)])
  }
  if(M == "M4I"){
    co <- c(parvect[aE+rE+gE+2+(1:nDiet)],
            parvect[aE+rE+gE+3+nDiet],nDiet)
  }
  if(M == "M4E"){
    co <- c(rep(parvect[aE+rE+gE+3],nDiet),
            parvect[aE+rE+gE+3+(1:nDiet)])
  }
  
  return(list(a=a,b=b,kappa=kapp,gamLam=gamLam,co=co))
}

###
# Function that runs the numerical maximization of the above likelihood function and returns the results
HSMMmle <- function(OBS,a0,b0,kappa0,gammaLam0,co0,M){
  parvect0 <- move.HSMM.pn2pw(a0,b0,kappa0,gammaLam0,co0)
  mod <- nlm(HSMMmllk,parvect0,OBS,M,print.level=2,hessian=TRUE,stepmax=stepm,iterlim=4000) ## hessian=TRUE only for confidence intervals 
  pn <- move.HSMM.pw2pn(mod$estimate, M)
  modAIC <- mod$minimum*2 + 2*length(parvect0)
  list(a=pn$a,b=pn$b,kappa=pn$kappa,gamLam=pn$gamLam,co=pn$co,H=mod$hessian,
       mllk=mod$minimum, AIC=modAIC,OptCode=mod$code)
}


#########################################################
# Neg. log likelihood functions

HSMMmllk <- function(parvect,OBS, M){
  n.ind <- ncol(OBS)/2
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  mllk.all <- 0 # Starting the likelihood
  for (ani in 1:n.ind){
    gamma <- genGamma(m,lpn$gamLam[c(1,nDiet+1)+(bearsCap$F_Group[ani]-1)]) # Creating transition probility matrix
    delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    allprobs <- matrix(rep(1,sum(m)*n),nrow=n)
    
    # For behaviour 1
    # Step length probability
    allprobs[!is.na(obs[,1]),1:m[1]] <- dweibull(obs[!is.na(obs[,1]),1],
                                                 shape=lpn$b[1],
                                                 scale=lpn$a[bearsCap$F_Group[ani]])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),1:m[1]] <- dwrpcauchy(obs[!is.na(obs[,2]),2],
                                                   mu=lpn$co[bearsCap$F_Group[ani]],
                                                   rho=lpn$kappa[bearsCap$F_Group[ani]])*allprobs[!is.na(obs[,2]),1:m[1]]
    
    # For behaviour 2
    # Step length probability
    allprobs[!is.na(obs[,1]),(m[1]+1):sum(m)] <- dweibull(obs[!is.na(obs[,1]),1],
                                                          shape=lpn$b[2],
                                                          scale=lpn$a[nDiet+bearsCap$F_Group[ani]])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)] <- dwrpcauchy(obs[!is.na(obs[,2]),2],
                                                            mu=lpn$co[nDiet+bearsCap$F_Group[ani]],
                                                            rho=lpn$kappa[nDiet+bearsCap$F_Group[ani]])*allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)]
    
    foo <- delta  
    lscale <- 0
    for (i in 1:n){
      foo <- foo%*%gamma*allprobs[i,]  
      sumfoo <- sum(foo)
      lscale <- lscale+log(sumfoo)
      foo <- foo/sumfoo
    }
    mllk.all <- mllk.all-lscale  
  } 
  return(mllk.all)
}

########################################
# Parameters used accros models
## size of state aggregates
m <- c(20,20)
# Number of diets in dataset
nDiet <- length(unique(bearsCap$F_Group))
#n.ind <- ncol(OBS)/2
# Parameters that we do not change
b0 <- c(0.8,0.85) # Weibull shape parameters
kappa0 <- c(0.2,0.4) # wrapped Cauchy concentration parameters
co0 <- c(pi,0) # wrapped Cauchy mean parameters
gammaLam0 <- c(1,2) # poisson lambda state dwell-time number of step
stepm <- 500

############################################
# Fitting models

###
# Fitting M0
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8,0.9) # Weibull scale parameters - for all individual

## run the numerical maximization
m0HSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0,"M0")
m0HSMM

# means (note that I reversed a & b compared to dweibull)
m0HSMM$a[1]*gamma(1+1/m0HSMM$b[1])
m0HSMM$a[2]*gamma(1+1/m0HSMM$b[2])

###
# Fitting M1B
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(rep(0.8,nDiet),rep(0.9,nDiet)) # Weibull scale parameters - for all individual

## run the numerical maximization
m1bHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0,"M1B")
m1bHSMM

###
# Fitting M1I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(rep(0.8,nDiet), 0.9) # Weibull scale parameters - for all individual

## run the numerical maximization
m1iHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0,"M1I")
m1iHSMM

###
# Fitting M1E
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, rep(0.9,nDiet)) # Weibull scale parameters - for all individual

## run the numerical maximization
m1eHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0,"M1E")
m1eHSMM

###
# Fitting M2B
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, 0.9) # Weibull scale parameters - for all individual
gammaLam0 <- c(rep(3,nDiet),rep(1,nDiet))  # state dwell-time parameters

## run the numerical maximization
m2bHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0, "M2B")
m2bHSMM

###
# Fitting M2I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
gammaLam0 <- c(rep(0.39,nDiet),3.75)  # negative binomial state dwell-time size parameters

## run the numerical maximization
m2iHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0, "M2I")
m2iHSMM

###
# Fitting M2E
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
gammaLam0 <- c(0.39,rep(3.75,nDiet))  # negative binomial state dwell-time size parameters

## run the numerical maximization
m2eHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0,"M2E")
m2eHSMM

###
# Fitting M3B
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, 0.9) # Weibull scale parameters - for all individual
gammaLam0 <- c(1,2) # poisson lambda state dwell-time number of step
kappa0 <- c(rep(0.2,nDiet),rep(0.4,nDiet)) # wrapped Cauchy concentration parameters

## run the numerical maximization
m3bHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0, "M3B")
m3bHSMM


###
# Fitting M3I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
kappa0 <- c(rep(0.2,nDiet),0.4) # wrapped Cauchy concentration parameters

m3iHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0, "M3I")
m3iHSMM

###
# Fitting M3E
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
kappa0 <- c(0.2,rep(0.4,nDiet)) # wrapped Cauchy concentration parameters

m3eHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0,"M3E")
m3eHSMM

###
# Fitting M4B
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
kappa0 <- c(0.2,0.4) # wrapped Cauchy concentration parameters
co0 <- rep(c(pi,0),each=nDiet) # wrapped Cauchy mean parameters

m4bHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0, "M4B")
m4bHSMM

###
# Fitting M4I
co0 <- c(rep(pi,nDiet),0) # wrapped Cauchy mean parameters
m4iHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0, "M4I")
m4iHSMM

###
# Fitting M4E
co0 <- c(pi,rep(0,nDiet)) # wrapped Cauchy mean parameters
m4eHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaLam0,co0, "M4E")
m4eHSMM


# Compare AIC
modAIC <- c(M0=m0HSMM$AIC,M1B=m1bHSMM$AIC,M1I=m1iHSMM$AIC,M1E=m1eHSMM$AIC,M2B=m2bHSMM$AIC,M2I=m2iHSMM$AIC,M2E=m2eHSMM$AIC)
cbind(AIC=modAIC, deltaAIC=modAIC - min(modAIC), 
      OptCode=c(m0HSMM$OptCode,
                m1bHSMM$OptCode,m1iHSMM$OptCode,m1eHSMM$OptCode,
                m2bHSMM$OptCode,m2iHSMM$OptCode,m2eHSMM$OptCode))

m1bHSMM$a
m1iHSMM$a
m1eHSMM$a
m2bHSMM$a
m2iHSMM$a
m2eHSMM$a

m1bHSMM$gamLam
m1iHSMM$gamLam
m1eHSMM$gamLam
m2bHSMM$gamLam
m2iHSMM$gamLam
m2eHSMM$gamLam
