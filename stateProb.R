HSMMlalb <- function(parvect, OBS, M, ani){

  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters

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
  
  ###
  # Initialising the alpha and beta vectors and creating a matrix will holds alpha and beta values for all t.
  # Note that to limit the underflow problem we are evaluating them as log values. 
  lalpha <- lbeta <- matrix(NA,sum(m),n)
  
  foo <- delta  
  lscale <- 0
  for (i in 1:n){
    foo <- foo%*%gamma*allprobs[i,]  
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
    lalpha[,i] <- log(foo) + lscale
  }
  
  lbeta[,n] <- rep(0,sum(m)) # beta_T = 1 so log(beta_T)=0 
  foo <- rep(0.5,sum(m)) # beta_T/w_T =psi_T 
  lscale <- log(sum(m)) # log(w_T)
  for (i in (n-1):1){
    foo <- gamma %*% (allprobs[i+1,]*foo) # gamma%*%P(x_{t+1}))*psi_{t+1} = beta_t/w_{t+1}
    lbeta[,i] <- log(foo) + lscale # log(beta_t) - log(w_{t+1} + log(w_{t+1}) = log(beta_t)
    sumfoo <- sum(foo) # w_t/w_{t+1}
    foo <- foo/sumfoo # (beta_t/w_{t+1})/(w_t/w_{t+1}) = psi_t
    lscale <- lscale + log(sumfoo) # log(w_{t+1}) + log(w_t) - log(w_{t+1}) = log(w_t)
  }
  
  return(list(la=lalpha, lb=lbeta))
}

# for hsm with poisson
HSMMw <- function(parvect, OBS, M, ani){
  
  lpn <- move.HSMM.pw2pn(parvect,M) # Transforming the parameters
  
  gamma <- genGamma(m,lpn$gamLam[c(1,nDiet+1)+(bearsCap$F_Group[ani]-1)]) # Creating transition probility matrix
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m))) # Getting the probility of the first step - stationary distribution
  
  # This calculates the weights for the pseudo-residuals
  obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
  n <- max(which(!is.na(obs[,1])))
  fb <- HSMMlalb(parvect, OBS, M, ani)
  la <- fb$la
  lb <- fb$lb
  la <- cbind(log(delta),la)
  lafact <- apply(la,2,max)
  lbfact <- apply(lb,2,max)
  w <- matrix(NA,sum(m),n)
  for (i in 1:n){
    foo <- (exp(la[,i]-lafact[i])%*%gamma)*
      exp(lb[,i]-lbfact[i])
    w[,i] <- foo/sum(foo)
  }
  w <- rbind(colSums(w[1:m[1],]),colSums(w[(m[1]+1):sum(m),]))
  return(w)
}

library(zoo) # for moving window

parMLE <- move.HSMM.pn2pw(m1bHSMM$a,m1bHSMM$b,m1bHSMM$kappa,m1bHSMM$gamLam,m1bHSMM$co)
nInd <- ncol(OBS)/2
ww <- matrix(NA, ncol=nInd, nrow=nrow(OBS))
for(ani in 1:nInd){
  obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
  n <- max(which(!is.na(obs[,1])))
  ww[1:n,ani] <- HSMMw(parMLE, OBS, "M1B", ani)[1,]
}

wSize <- 6*2
wi <- rowMeans(ww, na.rm=TRUE)
wi <- zoo(wi)
wim <- rollapply(wi, wSize, mean, na.rm=TRUE, partial=TRUE)
wisd <- apply(ww, MARGIN=1, function(x) sd(x, na.rm = TRUE))
wisdm <- rollapply(wisd, wSize, mean, na.rm=TRUE, partial=TRUE)
plot(wim, ylim=c(0,1))
lines(wim+wisd, col="red")
lines(wim-wisd, col="red")

# Just FG
FG1 <- bearsCap$F_Group == 1
wi <- rowMeans(ww[,FG1], na.rm=TRUE)
wi <- zoo(wi)
wim <- rollapply(wi, wSize, mean, na.rm=TRUE, partial=TRUE)
wisd <- apply(ww[,FG1], MARGIN=1, function(x) sd(x, na.rm = TRUE))
wisdm <- rollapply(wisd, wSize, mean, na.rm=TRUE, partial=TRUE)
plot(wim, ylim=c(0,1), ty="l")
lines(wim+wisd, col="red")
lines(wim-wisd, col="red")
