library(psbcGroup)
#データ数
n <- 100
#共変量の数
p <- 20
#目的変数に関連する共変量の数
psss<-10
#betaベクトルの真値
beta.true <- c(rep(log(2), psss), rep(0, (p-psss)))
beta.true
#外れ値の割合
outrate<-0.2
#重みベクトルの真値
W<- c(rep(1, n*(1-outrate)), rep(1/exp(5),n*(outrate/2)),rep(1/exp(5),n*(outrate/2)))
w.true<-log(W)
correct<-c(rep(1, psss), rep(0, (p-psss)))
correct
CovX<- matrix(0, nrow=p, ncol=p)
for(i in 1:p){
  for(j in 1:p){
    CovX[i,j] <- 0^abs(i-j)
  }
}
diag(CovX) <- 1
mu<-rep(0, p)
survObj <- list()
#KLはプログラムの試行回数
KL<-1
#プログラムごとにパラメータを保存するための準備
dat<-matrix(c(rep(0, p)),nrow=KL,ncol=p)
meanbeta<-matrix(c(rep(0, p)),nrow=KL,ncol=p)
medianbeta<-matrix(c(rep(0, p)),nrow=KL,ncol=p)
meanlambda<-matrix(c(rep(0, 1)),nrow=KL,ncol=1)
medianlambda<-matrix(c(rep(0, 1)),nrow=KL,ncol=1)
time<-matrix(c(rep(0, 5)),nrow=KL,ncol=5)
sum<-matrix(c(rep(0, KL)),nrow=KL,ncol=1)
acceptbeta<-matrix(c(rep(0, p)),nrow=KL,ncol=p)
MSE<-matrix(c(rep(0, 1)),nrow=KL,ncol=1)

for(i in 1:KL){
  start_time<-proc.time()
  survObj$x <-rmvnorm(n,mean=mu,sigma=CovX)
  pred<-as.vector(exp(rowSums(survObj$x%*%beta.true)+w.true))
  hist(rowSums(survObj$x%*%beta.true))
  t<-rexp(n,rate=pred)
  t1<-t[1:(n*(1-outrate))]
  t2<-t[(n*(1-outrate)+1):n]
  par(mfrow=c(1,2))
  boxplot(t1,ylim=c(0,50));boxplot(t2,ylim=c(0,50))
  x=quantile(rexp(n*(1-outrate), rate = pred[1:(n*(1-outrate))]), 1.0, type = 1)
  x2=quantile(rexp(n*outrate, rate = pred[(n*(1-outrate)+1):n])
              , 1.0, type = 1)
  cen1 <- runif((n*(1-outrate)), 0, 2*x)
  cen2<- runif((n*outrate), 0, 2*x2)
  #tとcenの様子をそれぞれ比較して、小さい要素を採用
  survObj$t1<- pmin(t1, cen1)
  survObj$di1 <- as.numeric(t1 <= cen1)
  survObj$di1
  1-sum(survObj$di1)/(n*(1-outrate))
  #tとcenの様子をそれぞれ比較して、小さい要素を採用
  survObj$t2<- pmin(t2, cen2)
  survObj$t <-c(survObj$t1, survObj$t2)
  survObj$di2 <- as.numeric(t2 <= cen2)
  survObj$di2
  1-sum(survObj$di2)/(n*(outrate))
  #打ち切りかどうかを示すラベル
  survObj$di <-c(survObj$di1, survObj$di2)
  survObj$di
  hist(survObj$di)
  sum[i,]<-sum(survObj$di)/100
  ##############################################
  priorPara <- list()
  #ハザード関数のパラメータ設定
  priorPara$eta0 <- 1
  priorPara$kappa0 <- 1
  priorPara$c0 <- 2
  #λ二乗の分布のパラメータ
  priorPara$r <- 1.5
  priorPara$delta <- 0.0001
  
  priorPara$s <- sort(survObj$t[survObj$di == 1])
  priorPara$s
  
  ###質問
  priorPara$s <- c(priorPara$s, 2*max(survObj$t)
                   -max(survObj$t[-which(survObj$t==max(survObj$t))]))
  priorPara$s
  priorPara$J <- length(priorPara$s)
  priorPara$groupInd <- c(rep(1,p))
  priorPara$J
  priorPara$groupInd
  
  mcmcPara <- list()
  #1回の反復で更新されるβの成分数
  mcmcPara$numBeta <- p
  mcmcPara$numBeta <- n
  
  initial <- list()
  #βの初期値
  #initial$beta.ini <- c(rep(0.5,p))
  
  initial$beta.ini<-c(rep(0.5, p))
  #λ二乗の初期値
  #initial$lambdaSq <- 1
  initial$lambdaSq <- 1
  #σの初期値
  initial$sigmaSq <- runif(1, 0.1, 10)
  initial$sigmaSq
  #質問
  #λ2乗の分布、τの初期値
  initial$tauSq <- rexp(length(unique(priorPara$groupInd)),
                        rate = initial$lambdaSq/2)
  initial$tauSq
  #γ分布の乱数、ハザード関数の初期値
  initial$h <- rgamma(priorPara$J, 1, 1)
  initial$h
  
  setting.interval <-
    function(y, delta, s, J){
      
      n <- length(y)
      
      smax	<- max(s)
      
      case0 <- which(delta == 0)
      case1 <- which(delta == 1)	
      
      case0yleq <- which(delta == 0 & y <= smax)
      case0ygeq <- which(delta == 0 & y > smax)
      case1yleq <- which(delta == 1 & y <= smax)
      case1ygeq <- which(delta == 1 & y > smax)
      
      
      ind.d <- ind.r <- matrix(0, n, J)
      
      for(i in case1yleq){
        d.mat.ind	<- min(which(s - y[i] >=0))
        ind.d[i, d.mat.ind]	<- 1
        ind.r[i, 1:d.mat.ind] <- 1		
      }
      
      for(i in case0yleq){
        cen.j <- min(which(s - y[i] >=0))		
        ind.r[i, 1:cen.j]	<- 1			
      }	
      
      if(length(union(case1ygeq, case0ygeq)) > 0){
        ind.r[union(case1ygeq, case0ygeq),]	<- 1
      }		
      
      ind.r_d	<- ind.r - ind.d;
      
      d	<- colSums(ind.d)
      
      list(ind.r = ind.r, ind.d = ind.d, d = d, ind.r_d = ind.r_d)
    }
  #βの更新、MH法を利用  
  UpdateRPrw <-
    function(survObj, priorPara, mcmcPara, ini){
      
      n 				<- survObj$n
      p				<- survObj$p
      x				<- survObj$x	
      
      J				<- priorPara$J
      ind.r			<- priorPara$ind.r
      ind.d			<- priorPara$ind.d
      ind.r_d			<- priorPara$ind.r_d		
      numBeta			<- mcmcPara$numBeta	
      beta.prop.me	<- mcmcPara$beta.prop.me
      beta.prop.var	<- mcmcPara$beta.prop.var
      
      xbeta			<- ini$xbeta
      wxbeta    <-ini$wxbeta
      be.ini			<- ini$beta.ini
      h				<- ini$h
      sd.bee		<- ini$sd.be
      
      updatej	<- c(1:p)
      accept 	<- rep(0, p)
      
      
      for(j in updatej){
        
        be.prop <- be.ini
        
        xbeta[xbeta > 700] <- 700
        exp.wxbeta		<- exp(wxbeta)
        exp.wxbeta.mat	<- matrix(rep(exp.wxbeta, J), n, J)
        
        first.sum 	<- colSums(exp.wxbeta.mat*ind.r_d)
        
        h.mat					<- matrix(rep(h, n), n, J, byrow = T)
        h.exp.xbeta.mat 	<- - h.mat * exp.wxbeta.mat
        h.exp.xbeta.mat[h.exp.xbeta.mat > -10^(-7)] <- -10^(-7)
        second.sum 	<- colSums(log(1 - exp(h.exp.xbeta.mat))*ind.d)
        #対数尤度関数
        loglh.ini <- sum(-h*first.sum + second.sum)
        
        be.prop[j] <- rnorm(1, mean =be.ini[j] , sd = 1)
        
        wxbeta.prop			<- wxbeta - x[,j]*be.ini[j] + x[,j]*be.prop[j] 
        wxbeta.prop[wxbeta.prop > 700] <- 700
        exp.wxbeta.prop		<- exp(wxbeta.prop)
        exp.wxbeta.mat.prop	<- matrix(rep(exp.wxbeta.prop, J), n, J)
        
        first.sum.prop 	<- colSums(exp.wxbeta.mat.prop*ind.r_d)
        
        h.exp.xbeta.mat.prop 	<- - h.mat * exp.wxbeta.mat.prop
        h.exp.xbeta.mat.prop[h.exp.xbeta.mat.prop > -10^(-7)] <- -10^(-7)
        second.sum.prop 	<- colSums(log(1 - exp(h.exp.xbeta.mat.prop))*ind.d)
        
        loglh.prop <- sum(-h*first.sum.prop + second.sum.prop)
        
        #対数尤度関数、確認必要
        logprior.prop   <- dnorm(be.prop[j] , mean = 0.5 , sd = sd.bee[j], log = TRUE)
        logprior.ini    <- dnorm(be.ini[j]  , mean = 0.5 , sd = sd.bee[j], log = TRUE) 
        logprop.prop    <- dnorm(be.prop[j] , mean = be.ini[j], sd = 1, log = TRUE)
        logprop.ini     <- dnorm(be.ini[j]  , mean = be.prop[j], sd = 1, log = TRUE)
        
        logR  <- loglh.prop - loglh.ini + logprior.prop - logprior.ini + logprop.ini - logprop.prop
        
        u = log(runif(1)) < logR
        
        if(u == 1){
          be.ini[j] <- be.prop[j]
          wxbeta	<- wxbeta.prop
        }
        
        accept[j]<-accept[j] + u
        
      } # end of for loop for j
      
      list(beta.ini = be.ini, accept = accept, wxbeta = wxbeta)
      
    }
  
  #hの更新
  UpdateBH <-
    function(survObj, priorPara, ini){
      
      n <- survObj$n
      p <- survObj$p	
      
      xbeta		<- ini$xbeta
      wxbeta		<- ini$wxbeta
      ind.r_d		<- priorPara$ind.r_d
      c0			<- priorPara$c0
      hPriorSh	<- priorPara$hPriorSh	
      d			<- priorPara$d	
      J			<- priorPara$J
      
      exp.wxbeta	<- exp(wxbeta)
      exp.wxbeta.mat	<- matrix(rep(exp.wxbeta, J), n, J)
      
      h.rate	<-  colSums(exp.wxbeta.mat * ind.r_d) + c0
      h		<- rgamma(J, shape = hPriorSh + d, rate = h.rate)
      
      return(h)
      
    }
  #Τの更新
  UpdateTau.GL <-
    function(survObj, priorPara, ini){
      
      n 				<- survObj$n
      p				<- survObj$p
      lambdaSq	<- ini$lambdaSq
      sigmaSq		<- ini$sigmaSq
      tauSq		<- ini$tauSq
      be.normSq	<- ini$be.normSq
      
      K			<- priorPara$K<-1
      groupInd	<- priorPara$groupInd
      groupNo		<- priorPara$groupNo
      
      nu.ind<-NULL
      nu=sqrt(lambdaSq * sigmaSq/be.normSq)
      nu.ind <- which(nu == Inf)
      if(length(nu.ind) > 0){nu[nu.ind] <- max(nu[-nu.ind]) + 10}
      
      gam <- c()
      
      for (j in 1:K){
        repeat{
          gam[j]  <- rinvGauss(1, nu = nu[j], lambda = lambdaSq)
          if (gam[j] > 0) break    	
        }
        tauSq[j] <- 1/gam[j]
      }
      
      return(tauSq)	
      
    }
  
  #σの更新
  UpdateSigma.GL <-
    function(survObj, priorPara, ini){
      
      p			<- survObj$p
      be.normSq	<- ini$be.normSq
      tauSq		<- ini$tauSq 
      
      sh.sig     <- p/2
      rate.sig   <- 1/2*sum(be.normSq/tauSq)
      
      sigmaSq     <- rigamma(1, a = sh.sig, b = rate.sig)
      
      return(sigmaSq)
      
    }
  #λの更新
  UpdateLambda.GL <-
    function(survObj, priorPara, ini){
      
      p		<- survObj$p
      K		<- priorPara$K
      tauSq	<- ini$tauSq 
      
      r		<- priorPara$r
      delta	<- priorPara$delta
      #確認必要   
      lambdaSq	<- rgamma(1, shape = (p + K)/2 + r, rate = delta + tauSq/2)
      
      return(lambdaSq)
      
    }
  
  #全体のパラメータを更新するためのアルゴリズム  
  psbcGLL<-function (survObj, priorPara, initial, rw , mcmcPara, 
                     num.reps, thin, chain, save) 
  {
    survObj$n <- n <- length(survObj$t)
    survObj$p <- p <- dim(survObj$x)[2]
    eta0 <- priorPara$eta0
    kappa0 <- priorPara$kappa0
    c0 <- priorPara$c0
    r <- priorPara$r
    delta <- priorPara$delta
    s <- priorPara$s
    J <- priorPara$J <- length(priorPara$s)
    groupInd <- priorPara$groupInd
    groupNo <- priorPara$groupNo <- unique(priorPara$groupInd)
    K <- priorPara$K <- length(groupNo)
    m_k <- priorPara$m_k
    m_k <- rep(NA, K)
    for (i in 1:K) {
      m_k[i] <- sum(groupInd == groupNo[i])
    }
    priorPara$m_k <- m_k
    intv <- setting.interval(survObj$t, survObj$di, priorPara$s, 
                             priorPara$J)
    priorPara$ind.r <- intv$ind.r
    priorPara$ind.d <- intv$ind.d
    priorPara$ind.r_d <- intv$ind.r_d
    priorPara$d <- intv$d
    ini <- initial
    beta.ini <- ini$beta.ini
    lambdaSq <- ini$lambdaSq
    sigmaSq <- ini$sigmaSq
    tauSq <- ini$tauSq
    h <- ini$h
    mcmcPara$beta.prop.me <- beta.ini
    
    tauSq.exp <- rep(NA, p)
    for (i in 1:K) {
      tauSq.exp[groupInd == groupNo[i]] <- tauSq[i]
    }
    
    ini$sd.be <- sqrt(sigmaSq * tauSq.exp)
    ini$xbeta <- as.vector(survObj$x %*% beta.ini)
    ini$wxbeta <- w.true+ini$xbeta
    be.normSq <- c()
    for (i in 1:K) {
      be.normSq[i] <- sum(beta.ini[which(groupInd == groupNo[i])]^2)
    }
    ini$be.normSq <- be.normSq
    
    H.star <- alpha0 <- c()
    for (j in 1:J) {
      H.star[j] <- eta0 * s[j]^kappa0
      alpha0[j] <- c0 * H.star[j]
    }
    priorPara$hPriorSh <- diff(c(0, alpha0))
    mcmcOutcome <- list()
    mcmcOutcome$initial <- initial
    mcmcOutcome$priorPara <- priorPara
    beta.p <- beta.ini
    h.p <- h
    tauSq.p <- tauSq
    sigmaSq.p <- sigmaSq
    lambdaSq.p<-lambdaSq
    tauSq.exp.p<-tauSq.exp
    mcmcOutcome$sigmaSq.p <- sigmaSq
    mcmcOutcome$lambdaSq.p <- lambdaSq
    mcmcOutcome$accept.beta <- c(rep(0, p))
    
    outcomeSum <- list()
    dir.create("mcmcOutcome", showWarnings = FALSE)
    for (M in 1:num.reps) {
      if (M%%1000 == 0) {
        cat("Chain", chain, "Iteration", M, fill = TRUE)
      }
      if (rw == FALSE) {
        sampleRP <- 1
      }
      if (rw == TRUE) {
        sampleRP <- UpdateRPrw(survObj, priorPara, mcmcPara, 
                               ini)
      }
      beta.ini <- ini$beta.ini <- sampleRP$beta.ini
      wxbeta <- ini$wxbeta <- sampleRP$wxbeta
      mcmcOutcome$accept.beta <- mcmcOutcome$accept.beta + 
        sampleRP$accept
      for (i in 1:K) {
        be.normSq <- ini$be.normSq[i] <- sum(beta.ini[which(groupInd == 
                                                              groupNo[i])]^2)
      }
      
      
      h <- ini$h <- UpdateBH(survObj, priorPara, ini)
      
      tauSq <- ini$tauSq <- UpdateTau.GL(survObj, priorPara, 
                                         ini)
      for (i in 1:K) {
        tauSq.exp[groupInd == groupNo[i]] <- tauSq[i]
      }
      sigmaSq <- ini$sigmaSq <- UpdateSigma.GL(survObj, priorPara, 
                                               ini)
      lambdaSq <- ini$lambdaSq <- UpdateLambda.GL(survObj, 
                                                  priorPara, ini)
      ini$sd.be <- sqrt(sigmaSq * tauSq.exp)
      
      
      
      
      if (M%%thin == 0) {
        beta.p <- rbind(beta.p, beta.ini, deparse.level = 0)
        h.p <- rbind(h.p, h, deparse.level = 0)
        tauSq.p <- rbind(tauSq.p, tauSq, deparse.level = 0)
        lambdaSq.p <- rbind(lambdaSq.p, lambdaSq, deparse.level = 0)
        sigmaSq.p <- rbind(sigmaSq.p, sigmaSq, deparse.level = 0)
        tauSq.exp.p<-rbind(tauSq.exp.p,tauSq.exp , deparse.level = 0)
        
        mcmcOutcome$sigmaSq.p <- c(mcmcOutcome$sigmaSq.p, 
                                   sigmaSq)
        mcmcOutcome$lambdaSq.p <- c(mcmcOutcome$lambdaSq.p, 
                                    lambdaSq)
        mcmcOutcome$ini <- ini
      }
      for (j in 1:survObj$p) {
        if (M%/%thin > (20%/%thin)) {
          if (beta.ini[j] == beta.p[(M%/%thin + 1 - (20%/%thin)), 
                                    j]) {
            mcmcPara$beta.prop.me[j] <- beta.p[(M%/%thin + 
                                                  1), j]
          }
        }
      }
      
      
      
      if (M%%save == 0 | M == num.reps) {
        save(mcmcOutcome, file = paste("mcmcOutcome/otherAll.ch", 
                                       chain, ".Rdata", sep = ""))
        save(beta.p, file = paste("mcmcOutcome/betaAll.ch", 
                                  chain, ".Rdata", sep = ""))
        save(tauSq.p, file = paste("mcmcOutcome/tauSqAll.ch", 
                                   chain, ".Rdata", sep = ""))
        save(h.p, file = paste("mcmcOutcome/hAll.ch", 
                               chain, ".Rdata", sep = ""))
      }
    }
    ret <- list(beta.p = beta.p[1251:num.reps,], h.p = h.p[1251:num.reps,], tauSq.p = tauSq.p[1251:num.reps,],
                lambdaSq.p=lambdaSq.p[1251:num.reps,],sigmaSq.p = sigmaSq.p[1251:num.reps,],mcmcOutcome = mcmcOutcome,
                tauSq.exp.p=tauSq.exp.p[1251:num.reps,],
                t = survObj$t, di = survObj$di)
    class(ret) <- "psbcGL"
    return(ret)
  }
  
  #試行回数
  num.reps =2000
  chain = 1
  thin = 1
  save = 5
  
  fitGL2 <- psbcGLL(survObj, priorPara, initial, rw=TRUE, mcmcPara,
                    num.reps, thin, chain, save)
  #変数選択を行うための関数
  VSS <- function(fit, X)
  {
    p <- dim(fit$beta.p)[2]
    beta.me <- apply(fit$beta.p, 2, mean)
    w.me <- w.true
    beta.me2<-abs(beta.me)
    sorted_indices <- order(abs(beta.me), decreasing = TRUE)
    beta.me3 <- beta.me[sorted_indices]
    sorted<- order(beta.me2, decreasing = TRUE)
    vec <- 1:p
    result_list <- lapply(seq_along(vec), function(i) vec[1:i])
    
    bicVec <- rep(NA, p)
    
    if(class(fit)=="aftGL")
    {
      w.me <- apply(fit$w.p, 2, mean)
      alpha.me <- mean(fit$alpha.p)
      sigSq.me <- mean(fit$sigSq.p)
      tauSq.me <- apply(fit$tauSq.p, 2, mean)
      
      nu0 <- fit$hyperParams$nu0
      sigSq0 <- fit$hyperParams$sigSq0
      alpha0 <- fit$hyperParams$alpha0
      h0 <- fit$hyperParams$h0
      
      for(i in 1:length(psiVec))
      {
        if(i %% 100 == 0){
          cat("BIC: i = ", i, "out of ", length(psiVec), sep = " ", fill = T)
        }
        betaHat <- beta.me
        
        ind0 <- c(1:p)[-selList[[i]]]
        betaHat[ind0] <- 0
        if(is.na(ind0[1]) & length(ind0) >0)
        {
          betaHat <- rep(0, p)
        }
        
        gam <- coef(lm(w.me ~ as.matrix(X)%*%betaHat))[2]
        if(is.na(gam))
        {
          gam = 1
        }
        
        bicVec[i] <- bic4(t=exp(w.me), x=X, alpha=alpha.me, beta=gam*betaHat, beta_all = beta.me, sigmaSq=sigSq.me, tauSq=tauSq.me, nu0, sigSq0, alpha0, h0)
      }
      
    }else if(class(fit)=="psbcEN" | class(fit)=="psbcFL" | class(fit)=="psbcGL")
    {
      for(i in 1:p)
      {
        betaHat <- beta.me3
        
        ind0 <- c(1:p)[-result_list[[i]]]
        betaHat[ind0] <- 0
        gam <- as.vector(coef(coxph(Surv(fit$t, fit$di) ~ as.matrix(X)%*%betaHat)))
        
        if(is.na(gam))
        {
          gam = 1
        }
        
        a.star <- coxph(Surv(fit$t, fit$di) ~ as.matrix(X)%*%(betaHat*gam)+w.me)
        a.all <- coxph(Surv(fit$t, fit$di) ~ as.matrix(X)%*%beta.me+w.me)
        
        bicVec[i] <- -2 * (a.star$loglik[2] - a.all$loglik[2])+(sum(betaHat != 0)-p)*log(length(fit$t))
      }
    }
    Sel.indd <- result_list[[which(bicVec==min(bicVec))[1]]]
    sorted_indices <- order(sorted)
    
    Sel.ind <-sorted_indices[1:length(Sel.indd)]
    Sel.ind <- sort(Sel.ind)
    list(Sel.ind=Sel.ind)
    
    ##
    cat("\nIndicators for variables selected based on SNC-BIC thresholding\n")
    print(Sel.ind)
  }
  
  VSS(fitGL2, X=survObj$x)
  
  
  a<-VSS(fitGL2, X=survObj$x)
  aa<-c(rep(0,p-length(a)))
  d<-c(a,aa)
  dat[i,]<-d
  
  beta.mean<- matrix(c(rep(0, p)),nrow=1,ncol=p)
  for(j in 1:p){
    beta.mean[,j]<- mean(fitGL2$beta.p[,j])
  }
  beta.mean
  meanbeta[i,]<-beta.mean
  beta.true
  mse <- ( meanbeta-beta.true )^2
  mse<-sum(mse)
  MSE[i]<-mse/p
  
  beta.median<- matrix(c(rep(0, p)),nrow=1,ncol=p)
  for(j in 1:p){
    beta.median[,j]<- median(fitGL2$beta.p[,j])
  }
  beta.median
  medianbeta[i,]<-beta.median
  
  
  lambda.mean<- matrix(c(rep(0, 1)),nrow=1,ncol=1)
  lambda.mean<- mean(fitGL2$lambdaSq.p)
  meanlambda[i]<-lambda.mean
  
  lambda.median<- matrix(c(rep(0, 1)),nrow=1,ncol=1)
  lambda.median<- median(fitGL2$lambdaSq.p)
  medianlambda[i]<-lambda.median
  
  acceptbeta[i,]<- fitGL2$mcmcOutcome$accept.beta
  end_time<-proc.time()
  runtime<-end_time-start_time
  runtime
  runtime <- as.vector(runtime)
  time[i,]<-runtime
}

library(dplyr)
TPR<- matrix(c(rep(0, KL)),nrow=1,ncol=KL)
TNR<- matrix(c(rep(0, KL)),nrow=1,ncol=KL)
PPV<- matrix(c(rep(0, KL)),nrow=1,ncol=KL)
NPV<- matrix(c(rep(0, KL)),nrow=1,ncol=KL)
dat <- dat[complete.cases(dat), ]
for(i in 1:KL){
  b<-dat
  b<- b[b != 0]
  b
  ##本番用
  x<-1:p
  y<-c(rep(0, p))
  d<-data.frame(x,y)
  colnames(d)<-c("変数","zero")
  px<-b
  #px<-c("x1","x2","X3","X4","X5")
  py<-c(rep(1,length(px)))
  pds<-data.frame(px,py)
  colnames(pds)<-c("変数","one")
  
  
  library(dplyr)
  # データフレームをidをキーとして結合
  merged_df <- full_join(d,pds, by = "変数")
  
  merged_df$Pre<-merged_df$zero+merged_df$one
  
  merged_df[is.na(merged_df)] <- 0
  
  
  final<-data.frame(x,merged_df$Pre)
  colnames(final)<-c("変数","予測値")
  
  #正解用
  y<-c(rep(1,psss), rep(0, (p-psss)))
  sin<-data.frame(x,y)
  colnames(sin)<-c("変数","真値")
  
  datdat<-data.frame(final$予測値,sin$真値)
  colnames(datdat)<-c("予測","真値")
  
  new_row <- data.frame(
    予測 = 0,
    真値 = 0
  )
  datdat <- rbind(datdat, new_row)
  table(datdat)
  
  TPR[,i]=(table(datdat)[2,2]/psss)*100
  TPR
  TNR[,i]=((table(datdat)[1,1]-1)/(p-psss))*100
  TNR
  PPV[,i]=(table(datdat)[2,2]/(table(datdat)[2,1]+table(datdat)[2,2]))*100
  PPV
  NPV[,i]=((table(datdat)[1,1]-1)/((table(datdat)[1,1]-1)+table(datdat)[1,2]))*100
  NPV
}
mean(TPR)
mean(TNR)
mean(PPV)
NPV <-NPV[!is.na(NPV)]
mean(NPV)
