#########################################################
# Based on the code from Langrock et al. 2012
# Fit the different models to the bear population
# M0: All bears are the same (all parameters are the same accros diet)
# M1: Bears of different diet differ in their Step length (scale of Weibull)
  # M1I: Differ in intensive behaviour (potentially foraging)
  # M1E: Differ in extensive behaviour (potentially travelling)
  # M1B: Differ in both intensive and extensive behaviours
# M2: Bears of different diet differ in their the number of steps in a behaviour (negative binomial size)
  # M2I: Differ in intensive behaviour (potentially foraging)
  # M2E: Differ in extensive behaviour (potentially travelling)
  # M2B: Differ in both intensive and extensive behaviours

#####
# Load packages
library(CircStats) # for wrapped Cauchy
#library(boot) # not sure for what?

#####
# Get data for the bears

##
# Movement data
# observations ("OBS") need to be given in an n x (2*p) matrix, 
# with p being the number of individuals, 
# and the observed step lengths and turning angles for individual i 
# given in columns 2*(i-1)+1 and 2*(i-1)+2, respectively
#OBS <- read.csv("minidata.csv")
OBS <- read.csv("beardata.csv")

##
# Covariate data - diet!
bearsSI <- read.csv("GB_SI.csv")

# Link bear ID (diet) with PTT ID (movement)
bearsCap <- read.csv("MackenzieCollarDeploy.csv") 
# Sort by Year (so we get first capture)
bearsCap <- bearsCap[order(bearsCap$Date),]
bearsCap <- bearsCap[match(unique(substring(colnames(OBS),4,8)), bearsCap$PTT_ID),]
bearsCap <- cbind(bearsCap[,c("BearID","PTT_ID")], bearsSI[match(bearsCap$BearID, bearsSI$BearID),
                                    c("Sex", "deltaN", "dealtaC", "Group")])
bearsCap$GroupFL <- as.numeric(as.factor(bearsCap$Group))

#####################################################
# Define functions common to all models
## function that derives the t.p.m. of the HMM that represents the HSMM (see Langrock and Zucchini, 2011) 
gen.GammaO <- function(m,pSize,pSP){
  Gamma <- diag(sum(m))*0
  ## state aggregate 1
  Gamma[1,m[1]+1] <- dnbinom(0,size=pSize[1],prob=pSP[1])
  Gamma[1,2] <- 1-Gamma[1,m[1]+1]
  for (i in 2:(m[1]-1)){
    cc <- rep(1,sum(m))
    #cc[1:(i-1)] <- diag(matrix(Gamma[1:(i-1),2:i]))
    for (k in 1:(i-1)) {cc[k] <- Gamma[k,k+1]}
    dd <- prod(cc)
    if (dd>1e-12) Gamma[i,m[1]+1] <- dnbinom(i-1,size=pSize[1],prob=pSP[1])/dd
    if (dd<1e-12) Gamma[i,m[1]+1] <- 1
    Gamma[i,i+1] <- 1-Gamma[i,m[1]+1]
  } 
  cc <- rep(1,sum(m))
  for (k in 1:(m[1]-1)){cc[k] <- Gamma[k,k+1]}
  dd <- prod(cc)
  if (dd>1e-12) Gamma[m[1],m[1]+1] <- dnbinom(m[1]-1,size=pSize[1],prob=pSP[1])/dd
  if (dd<1e-12) Gamma[m[1],m[1]+1] <- 1
  Gamma[m[1],m[1]] <- 1-Gamma[m[1],m[1]+1] 
  
  ## state aggregate 2
  Gamma[m[1]+1,1] <- dnbinom(0,size=pSize[2],prob=pSP[2])
  Gamma[m[1]+1,m[1]+2] <- 1-Gamma[m[1]+1,1]
  for (i in 2:(m[2]-1)){
    cc <- rep(1,sum(m))
    for (k in 1:(i-1)) {cc[k] <- Gamma[m[1]+k,m[1]+k+1]}
    dd <- prod(cc)
    if (dd>1e-12) Gamma[m[1]+i,1] <- dnbinom(i-1,size=pSize[2],prob=pSP[2])/dd
    if (dd<1e-12) Gamma[m[1]+i,1] <- 1
    Gamma[m[1]+i,m[1]+i+1] <- 1-Gamma[m[1]+i,1]
  }
  cc<-rep(1,sum(m))
  for (k in 1:(m[2]-1)) {cc[k] <- Gamma[m[1]+k,m[1]+k+1]}
  dd<-prod(cc)
  if (dd>1e-12) Gamma[m[1]+m[2],1] <- dnbinom(m[2]-1,size=pSize[2],prob=pSP[2])/dd
  if (dd<1e-12) Gamma[m[1]+m[2],1] <- 1
  Gamma[m[1]+m[2],m[1]+m[2]] <- 1-Gamma[m[1]+m[2],1] 
  Gamma 
}

gen.Gamma <- function(m,pSize,pSP){
  Gamma <- diag(m[1]+m[2])*0
  # p(r), for r=1,2,...N* (note that r=1 is 0 in dnbinom)
  probs1 <- dnbinom(0:(m[1]-1),size=pSize[1],prob=pSP[1]) 
  probs2 <- dnbinom(0:(m[2]-1),size=pSize[2],prob=pSP[2])
  
  # Denominator of c(r): 1 - sum_{k=1}^{r-1}p(k), so for r=1 -> 0 (because empty sum equals 0)
  # Use cumulative distribution function because it is the sum of the prob
  den1 <- 1 - c(0,pnbinom(0:(m[1]-2),size=pSize[1],prob=pSP[1]))
  den2 <- 1 - c(0,pnbinom(0:(m[2]-2),size=pSize[2],prob=pSP[2]))
  
#   den1 <- 1 - c(0,cumsum(dnbinom(0:(m[1]-2),size=pSize[1],prob=pSP[1])))
#   den2 <- 1 - c(0,cumsum(dnbinom(0:(m[2]-2),size=pSize[2],prob=pSP[2])))
  
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

gen.Gamma.repar <- function(m,pSize,pMu){
  Gamma <- diag(sum(m))*0
  ## state aggregate 1
  Gamma[1,m[1]+1] <- dnbinom(0,size=pSize[1],mu=pMu[1])
  Gamma[1,2] <- 1-Gamma[1,m[1]+1]
  for (i in 2:(m[1]-1)){
    cc <- rep(1,sum(m))
    #cc[1:(i-1)] <- diag(matrix(Gamma[1:(i-1),2:i]))
    for (k in 1:(i-1)) {cc[k] <- Gamma[k,k+1]}
    dd <- prod(cc)
    if (dd>1e-12) Gamma[i,m[1]+1] <- dnbinom(i-1,size=pSize[1],mu=pMu[1])/dd
    if (dd<1e-12) Gamma[i,m[1]+1] <- 1
    Gamma[i,i+1] <- 1-Gamma[i,m[1]+1]
  } 
  cc <- rep(1,sum(m))
  for (k in 1:(m[1]-1)){cc[k] <- Gamma[k,k+1]}
  dd <- prod(cc)
  if (dd>1e-12) Gamma[m[1],m[1]+1] <- dnbinom(m[1]-1,size=pSize[1],mu=pMu[1])/dd
  if (dd<1e-12) Gamma[m[1],m[1]+1] <- 1
  Gamma[m[1],m[1]] <- 1-Gamma[m[1],m[1]+1] 
  
  ## state aggregate 2
  Gamma[m[1]+1,1] <- dnbinom(0,size=pSize[2],mu=pMu[2])
  Gamma[m[1]+1,m[1]+2] <- 1-Gamma[m[1]+1,1]
  for (i in 2:(m[2]-1)){
    cc <- rep(1,sum(m))
    for (k in 1:(i-1)) {cc[k] <- Gamma[m[1]+k,m[1]+k+1]}
    dd <- prod(cc)
    if (dd>1e-12) Gamma[m[1]+i,1] <- dnbinom(i-1,size=pSize[2],mu=pMu[2])/dd
    if (dd<1e-12) Gamma[m[1]+i,1] <- 1
    Gamma[m[1]+i,m[1]+i+1] <- 1-Gamma[m[1]+i,1]
  }
  cc<-rep(1,sum(m))
  for (k in 1:(m[2]-1)) {cc[k] <- Gamma[m[1]+k,m[1]+k+1]}
  dd<-prod(cc)
  if (dd>1e-12) Gamma[m[1]+m[2],1] <- dnbinom(m[2]-1,size=pSize[2],mu=pMu[2])/dd
  if (dd<1e-12) Gamma[m[1]+m[2],1] <- 1
  Gamma[m[1]+m[2],m[1]+m[2]] <- 1-Gamma[m[1]+m[2],1] 
  Gamma 
}


## function that transforms each of the (possibly constrained) parameters to the real line
# Function to transform parameters
#move.HSMM.pn2pw <- function(a,b,kappa,gamSize,gamSP,co){
move.HSMM.pn2pw <- function(a,b,kappa,gamSize,gamMu,co){
  ta <- log(a) # Weibull scale
  tb <- log(b)
  tkappa <- logit(kappa)
  tgamSize <- log(gamSize) # Size of neg. binom.
  #tgamSP <- log(gamSP/(1-gamSP)) # Succes probility of neg. binom.
  tgamMu <- log(gamMu) # Mean number of neg. binom.
  tco <- c(log((co[1])/(2*pi-co[1])),log((pi+co[2])/(pi-co[2])))
  #parvect <- c(ta,tb,tkappa,tgamSize,tgamSP,tco)
  parvect <- c(ta,tb,tkappa,tgamSize,tgamMu,tco)
  return(parvect)
}

## inverse transformation back to the natural parameter space            
move.HSMM.pw2pn <- function(parvect,M){
  epar <- exp(parvect)
  aE <- 2
  if(M == "M1B"){aE <- 2*nDiet}
  if(M == "M1I" | M == "M1E"){aE <- nDiet + 1}
  gE <- 2
  if(M == "M2B"){gE <- 2*nDiet}
  if(M == "M2I" | M == "M2E"){gE <- nDiet + 1}
  a <- epar[1:aE]
  b <- epar[aE+(1:2)]
  kappa <- inv.logit(parvect[aE+(3:4)])
  #gamSize <- exp(parvect[aE+(4+(1:gE))])
  gamSize <- exp(parvect[aE+5:6])
  #gamSP <- exp(parvect[aE+gE+(5:6)])/(exp(parvect[aE+gE+(5:6)])+1)
  gamMu <- exp(parvect[aE+6+(1:gE)])
  co <- c(2*pi*inv.logit(parvect[aE+gE+7]),pi*(exp(parvect[aE+gE+8])-1)/(exp(parvect[aE+gE+8])+1))
  #return(list(a=a,b=b,kappa=kappa,gamSize=gamSize,gamSP=gamSP,co=co))
  return(list(a=a,b=b,kappa=kappa,gamSize=gamSize,gamMu=gamMu,co=co))
}

###
# Function that runs the numerical maximization of the above likelihood function and returns the results
#HSMMmle <- function(OBS,a0,b0,kappa0,gammaSize0,gammaSP0,co0,mllk,M){
HSMMmle <- function(OBS,a0,b0,kappa0,gammaSize0,gammaMu0,co0,mllk,M){
  #parvect0 <- move.HSMM.pn2pw(a0,b0,kappa0,gammaSize0,gammaSP0,co0)
  parvect0 <- move.HSMM.pn2pw(a0,b0,kappa0,gammaSize0,gammaMu0,co0)
  mod <- nlm(mllk,parvect0,OBS,print.level=2,hessian=TRUE,stepmax=stepm,iterlim=4000) ## hessian=TRUE only for confidence intervals 
  mllk <- mod$minimum
  pn <- move.HSMM.pw2pn(mod$estimate, M)
  modAIC <- mod$minimum*2 + 2*length(parvect0)
  list(a=pn$a,b=pn$b,kappa=pn$kappa,H=mod$hessian,gamSize=pn$gamSize,gamMu=pn$gamMu,co=pn$co,mllk=mllk, AIC=modAIC,OptCode=mod$code)
  #list(a=pn$a,b=pn$b,kappa=pn$kappa,H=mod$hessian,gamSize=pn$gamSize,gamSP=pn$gamSP,co=pn$co,mllk=mllk, AIC=modAIC,OptCode=mod$code)
}


#########################################################
# Neg. log likelihood functions

# M0: All bears are the same (all parameters are the same accros diet)
m0.HSMM.mllk <- function(parvect,OBS){
  n.ind <- ncol(OBS)/2
  M <- "M0" # Define the model used
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  #gamma <- gen.Gamma(m,lpn$gamSize,lpn$gamSP) # Creating transition probility matrix
  gamma <- gen.Gamma.repar(m,lpn$gamSize,lpn$gamMu) # Creating transition probility matrix - reparametrised
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  mllk.all <- 0 # Starting the likelihood
  for (ani in 1:n.ind){  
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    allprobs <- matrix(rep(1,sum(m)*n),nrow=n)
    
    # For behaviour 1
    # Step length probability
    allprobs[!is.na(obs[,1]),1:m[1]] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[1],scale=lpn$a[1])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),1:m[1]] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[1],rho=lpn$kappa[1])*allprobs[!is.na(obs[,2]),1:m[1]]
    
    # For behaviour 2
    # Step length probability
    allprobs[!is.na(obs[,1]),(m[1]+1):sum(m)] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[2],scale=lpn$a[2])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[2],rho=lpn$kappa[2])*allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)]
     
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

#############
# M1: Bears of different diet differ in their Step length (scale of Weibull)
# M1B: Differ in both intensive and extensive behaviours
m1b.HSMM.mllk <- function(parvect,OBS){
  n.ind <- ncol(OBS)/2
  M <- "M1B" # Define the model used
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  #gamma <- gen.Gamma(m,lpn$gamSize,lpn$gamSP) # Creating transition probility matrix
  gamma <- gen.Gamma.repar(m,lpn$gamSize,lpn$gamMu) # Creating transition probility matrix - reparametrised
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  mllk.all <- 0 # Starting the likelihood
  for (ani in 1:n.ind){  
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    allprobs <- matrix(rep(1,sum(m)*n),nrow=n)
    
    # For behaviour 1
    # Step length probability
    allprobs[!is.na(obs[,1]),1:m[1]] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[1],scale=lpn$a[bearsCap$GroupFL[ani]])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),1:m[1]] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[1],rho=lpn$kappa[1])*allprobs[!is.na(obs[,2]),1:m[1]]
    
    # For behaviour 2
    # Step length probability
    allprobs[!is.na(obs[,1]),(m[1]+1):sum(m)] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[2],scale=lpn$a[nDiet+bearsCap$GroupFL[ani]])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[2],rho=lpn$kappa[2])*allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)]
  
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

# M1I: Differ in intensive behaviour (potentially foraging)
m1i.HSMM.mllk <- function(parvect,OBS){
  n.ind <- ncol(OBS)/2
  M <- "M1I" # Define the model used
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  #gamma <- gen.Gamma(m,lpn$gamSize,lpn$gamSP) # Creating transition probility matrix
  gamma <- gen.Gamma.repar(m,lpn$gamSize,lpn$gamMu) # Creating transition probility matrix - reparametrised
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  mllk.all <- 0 # Starting the likelihood
  for (ani in 1:n.ind){  
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    allprobs <- matrix(rep(1,sum(m)*n),nrow=n)
    
    # For behaviour 1
    # Step length probability
    allprobs[!is.na(obs[,1]),1:m[1]] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[1],scale=lpn$a[bearsCap$GroupFL[ani]])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),1:m[1]] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[1],rho=lpn$kappa[1])*allprobs[!is.na(obs[,2]),1:m[1]]
    
    # For behaviour 2
    # Step length probability
    allprobs[!is.na(obs[,1]),(m[1]+1):sum(m)] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[2],scale=lpn$a[nDiet+1])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[2],rho=lpn$kappa[2])*allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)]
    
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

# M1E: Differ in extensive behaviour (potentially travelling)
m1e.HSMM.mllk <- function(parvect,OBS){
  n.ind <- ncol(OBS)/2
  M <- "M1E" # Define the model used
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  #gamma <- gen.Gamma(m,lpn$gamSize,lpn$gamSP) # Creating transition probility matrix
  gamma <- gen.Gamma.repar(m,lpn$gamSize,lpn$gamMu) # Creating transition probility matrix - reparametrised
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  mllk.all <- 0 # Starting the likelihood
  for (ani in 1:n.ind){  
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    allprobs <- matrix(rep(1,sum(m)*n),nrow=n)
    
    # For behaviour 1
    # Step length probability
    allprobs[!is.na(obs[,1]),1:m[1]] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[1],scale=lpn$a[1])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),1:m[1]] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[1],rho=lpn$kappa[1])*allprobs[!is.na(obs[,2]),1:m[1]]
    
    # For behaviour 2
    # Step length probability
    allprobs[!is.na(obs[,1]),(m[1]+1):sum(m)] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[2],scale=lpn$a[1+bearsCap$GroupFL[ani]])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[2],rho=lpn$kappa[2])*allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)]
      
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

# M2: Bears of different diet differ in their the number of steps in a behaviour (negative binomial size)
# M2B: Differ in both intensive and extensive behaviours
m2b.HSMM.mllk <- function(parvect,OBS){
  n.ind <- ncol(OBS)/2
  M <- "M2B" # Define the model used
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  mllk.all <- 0 # Starting the likelihood
  for (ani in 1:n.ind){
    #gamma <- gen.Gamma(m,lpn$gamSize[c(1,nDiet+1)+(bearsCap$GroupFL[ani]-1)],lpn$gamSP) # Creating transition probility matrix
    gamma <- gen.Gamma.repar(m,lpn$gamSize,lpn$gamMu[c(1,nDiet+1)+(bearsCap$GroupFL[ani]-1)]) # Creating transition probility matrix - reparametrised
    delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    allprobs <- matrix(rep(1,sum(m)*n),nrow=n)
    
    # For behaviour 1
    # Step length probability
    allprobs[!is.na(obs[,1]),1:m[1]] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[1],scale=lpn$a[1])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),1:m[1]] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[1],rho=lpn$kappa[1])*allprobs[!is.na(obs[,2]),1:m[1]]
    
    # For behaviour 2
    # Step length probability
    allprobs[!is.na(obs[,1]),(m[1]+1):sum(m)] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[2],scale=lpn$a[2])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[2],rho=lpn$kappa[2])*allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)]
    
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

# M2I: Differ in intensive behaviour (potentially foraging)
m2i.HSMM.mllk <- function(parvect,OBS){
  n.ind <- ncol(OBS)/2
  M <- "M2I" # Define the model used
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  mllk.all <- 0 # Starting the likelihood
  for (ani in 1:n.ind){
    #gamma <- gen.Gamma(m,lpn$gamSize[c(bearsCap$GroupFL[ani], nDiet+1)],lpn$gamSP) # Creating transition probility matrix
    gamma <- gen.Gamma.repar(m,lpn$gamSize,lpn$gamMu[c(bearsCap$GroupFL[ani], nDiet+1)]) # Creating transition probility matrix - reparametrised
    delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    allprobs <- matrix(rep(1,sum(m)*n),nrow=n)
    
    # For behaviour 1
    # Step length probability
    allprobs[!is.na(obs[,1]),1:m[1]] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[1],scale=lpn$a[1])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),1:m[1]] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[1],rho=lpn$kappa[1])*allprobs[!is.na(obs[,2]),1:m[1]]
    
    # For behaviour 2
    # Step length probability
    allprobs[!is.na(obs[,1]),(m[1]+1):sum(m)] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[2],scale=lpn$a[2])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[2],rho=lpn$kappa[2])*allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)]

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

# M2E: Differ in extensive behaviour (potentially travelling)
m2e.HSMM.mllk <- function(parvect,OBS){
  n.ind <- ncol(OBS)/2
  M <- "M2E" # Define the model used
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  mllk.all <- 0 # Starting the likelihood
  for (ani in 1:n.ind){
    #gamma <- gen.Gamma(m,lpn$gamSize[c(1, 1+bearsCap$GroupFL[ani])],lpn$gamSP) # Creating transition probility matrix
    gamma <- gen.Gamma.repar(m,lpn$gamSize,lpn$gamMu[c(1, 1+bearsCap$GroupFL[ani])]) # Creating transition probility matrix - reparametrised
    delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    allprobs <- matrix(rep(1,sum(m)*n),nrow=n)
    
    # For behaviour 1
    # Step length probability
    allprobs[!is.na(obs[,1]),1:m[1]] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[1],scale=lpn$a[1])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),1:m[1]] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[1],rho=lpn$kappa[1])*allprobs[!is.na(obs[,2]),1:m[1]]
    
    # For behaviour 2
    # Step length probability
    allprobs[!is.na(obs[,1]),(m[1]+1):sum(m)] <- dweibull(obs[!is.na(obs[,1]),1],shape=lpn$b[2],scale=lpn$a[2])
    # Turn angle probability
    allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)] <- dwrpcauchy(obs[!is.na(obs[,2]),2],mu=lpn$co[2],rho=lpn$kappa[2])*allprobs[!is.na(obs[,2]),(m[1]+1):sum(m)]
    
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
nDiet <- length(unique(bearsCap$GroupFL))
#n.ind <- ncol(OBS)/2
# Parameters that we do not change
b0 <- c(0.8,0.85) # Weibull shape parameters
kappa0 <- c(0.2,0.4) # wrapped Cauchy concentration parameters
co0 <- c(pi,0) # wrapped Cauchy mean parameters
gammaSP0 <- c(0.14,0.59) # negative binomial state dwell-time success prob. parameters

stepm <- 500

############################################
# Fitting models

###
# Fitting M0
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8,0.9) # Weibull scale parameters - for all individual
gammaSize0 <- c(0.39,1)#3.75)  # negative binomial state dwell-time size parameters

## run the numerical maximization
m0HSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaSP0,co0,m0.HSMM.mllk,"M0")
m0HSMM

###
# Fitting M1B
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(rep(0.8,nDiet),rep(0.9,nDiet)) # Weibull scale parameters - for all individual
gammaSize0 <- c(0.39,3.75)  # negative binomial state dwell-time size parameters

## run the numerical maximization
m1bHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaSP0,co0,m1b.HSMM.mllk,"M1B")
m1bHSMM

###
# Fitting M1I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(rep(0.8,nDiet), 0.9) # Weibull scale parameters - for all individual
gammaSize0 <- c(0.39,3.75)  # negative binomial state dwell-time size parameters

## run the numerical maximization
m1iHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaSP0,co0,m1i.HSMM.mllk,"M1I")
m1iHSMM

###
# Fitting M1I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, rep(0.9,nDiet)) # Weibull scale parameters - for all individual
gammaSize0 <- c(0.39,3.75)  # negative binomial state dwell-time size parameters

## run the numerical maximization
m1eHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaSP0,co0,m1e.HSMM.mllk,"M1E")
m1eHSMM

###
# Fitting M2B
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, 0.9) # Weibull scale parameters - for all individual
gammaSize0 <- c(rep(0.39,nDiet),rep(3.75,nDiet))  # negative binomial state dwell-time size parameters

## run the numerical maximization
m2bHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaSP0,co0,m2b.HSMM.mllk,"M2B")
m2bHSMM

###
# Fitting M2I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, 0.9) # Weibull scale parameters - for all individual
gammaSize0 <- c(rep(0.39,nDiet),3.75)  # negative binomial state dwell-time size parameters

## run the numerical maximization
m2iHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaSP0,co0,m2i.HSMM.mllk,"M2I")
m2iHSMM

###
# Fitting M2E
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, 0.9) # Weibull scale parameters - for all individual
gammaSize0 <- c(0.39,rep(3.75,nDiet))  # negative binomial state dwell-time size parameters

## run the numerical maximization
m2eHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaSP0,co0,m2e.HSMM.mllk,"M2E")
m2eHSMM


# Compare AIC
modAIC <- c(M0=m0HSMM$AIC,M1B=m1bHSMM$AIC,M1I=m1iHSMM$AIC,M1E=m1eHSMM$AIC,M2B=m2bHSMM$AIC,M2I=m2iHSMM$AIC,M2E=m2eHSMM$AIC)
cbind(AIC=modAIC, deltaAIC=modAIC - min(modAIC), 
      OptCode=c(m0HSMM$OptCode,
                m1bHSMM$OptCode,m1iHSMM$OptCode,m1eHSMM$OptCode,
                m2bHSMM$OptCode,m2iHSMM$OptCode,m2eHSMM$OptCode))


############# With new parametrisation of neg binom

########################################
# Parameters used accros models
## size of state aggregates
m <- c(20,20)
# Number of diets in dataset
nDiet <- length(unique(bearsCap$GroupFL))
#n.ind <- ncol(OBS)/2
# Parameters that we do not change
b0 <- c(0.8,0.85) # Weibull shape parameters
kappa0 <- c(0.2,0.4) # wrapped Cauchy concentration parameters
co0 <- c(pi,0) # wrapped Cauchy mean parameters
#gammaSP0 <- c(0.14,0.59) # negative binomial state dwell-time success prob. parameters
gammaSize0 <- c(0.39,1)#3.75)  # negative binomial state dwell-time size parameters

stepm <- 500

############################################
# Fitting models

###
# Fitting M0
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8,0.9) # Weibull scale parameters - for all individual
gammaMu0 <- c(1,3)  # negative binomial state dwell-time mu parameters

## run the numerical maximization
m0HSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaMu0,co0,m0.HSMM.mllk,"M0")
m0HSMM

###
# Fitting M1B
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(rep(0.8,nDiet),rep(0.9,nDiet)) # Weibull scale parameters - for all individual
gammaMu0 <- c(1,3)  # negative binomial state dwell-time mu parameters

## run the numerical maximization
m1bHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaMu0,co0,m1b.HSMM.mllk,"M1B")
m1bHSMM

###
# Fitting M1I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(rep(0.8,nDiet), 0.9) # Weibull scale parameters - for all individual
gammaMu0 <- c(1,3)  # negative binomial state dwell-time mu parameters

## run the numerical maximization
m1iHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaMu0,co0,m1i.HSMM.mllk,"M1I")
m1iHSMM

###
# Fitting M1I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, rep(0.9,nDiet)) # Weibull scale parameters - for all individual
gammaMu0 <- c(1,3)  # negative binomial state dwell-time mu parameters

## run the numerical maximization
m1eHSMM <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaMu0,co0,m1e.HSMM.mllk,"M1E")
m1eHSMM

###
# Fitting M2B
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, 0.9) # Weibull scale parameters - for all individual
gammaMu0 <- c(rep(1,nDiet),rep(3,nDiet))  # negative binomial state dwell-time mu parameters

## run the numerical maximization
m2bHSMMr <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaMu0,co0,m2b.HSMM.mllk,"M2B")
m2bHSMMr

###
# Fitting M2I
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, 0.9) # Weibull scale parameters - for all individual
gammaMu0 <- c(rep(1,nDiet),3)  # negative binomial state dwell-time mu parameters


## run the numerical maximization
m2iHSMMr <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaMu0,co0,m2i.HSMM.mllk,"M2I")
m2iHSMMr

###
# Fitting M2E
## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8, 0.9) # Weibull scale parameters - for all individual
gammaMu0 <- c(1,rep(3,nDiet))  # negative binomial state dwell-time mu parameters

## run the numerical maximization
m2eHSMMr <- HSMMmle(OBS,a0,b0,kappa0,gammaSize0,gammaMu0,co0,m2e.HSMM.mllk,"M2E")
m2eHSMMr


# Compare AIC
modAIC <- c(M0=m0HSMM$AIC,M1B=m1bHSMM$AIC,M1I=m1iHSMM$AIC,M1E=m1eHSMM$AIC,M2B=m2bHSMMr$AIC,M2I=m2iHSMMr$AIC,M2E=m2eHSMMr$AIC)
cbind(AIC=modAIC, deltaAIC=modAIC - min(modAIC), 
      OptCode=c(m0HSMM$OptCode,
                m1bHSMM$OptCode,m1iHSMM$OptCode,m1eHSMM$OptCode,
                m2bHSMMr$OptCode,m2iHSMMr$OptCode,m2eHSMMr$OptCode))

modAIC <-  c(M2Bo=m2bHSMM$AIC,M2Io=m2iHSMM$AIC,M2B=m2bHSMMr$AIC,M2I=m2iHSMMr$AIC,M2E=m2eHSMMr$AIC) 

###########
bearModels<- read.csv("HSMMini0.csv")
2*sum(bearModels[1:4,"mllk"]) + 2*(14*4)
moveHSMM$mllk*2 + 2*18
