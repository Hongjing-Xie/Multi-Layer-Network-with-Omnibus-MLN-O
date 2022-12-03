rm(list=ls())
library(stats); library(MASS); library(Matrix); library(SPAtest)
get_community <- function(Weight,A)
{
  m=ncol(Weight)
  B=list()
  for(i in 2:(m-1)) B[[i-1]]=matrix(0,m,i)
  Sigma=Weight
  dist=1-Sigma
  hc=hclust(as.dist(dist), method = "ward.D2")
  Q <- c()
  for(i in 2:(m-1))
  {
    index=cutree(hc,i)
    B[[i-1]]=sapply(1:i, function(t) as.numeric(index==t))
    Q[i-1] <- get_modularity(A, B[[i-1]])
  }
  
  IND=which(Q==max(Q))
  L=length(IND)
  if (L>1) IND=IND[L]
  khat=IND+1
  Bb=B[[khat-1]]
  return(Bb)
}
# Get modularity
get_modularity <- function(Weight, B)
{
  if (dim(Weight)[1] == 1){
    Q <- 0
  } else {
    W_pos <- Weight * (Weight > 0)
    W_neg <- Weight * (Weight < 0)
    N <- dim(Weight)[1]
    K_pos <- colSums(W_pos)
    K_neg <- colSums(W_neg)
    m_pos <- sum(K_pos)
    m_neg <- sum(K_neg)
    m <- m_pos + m_neg
    cate <- B %*% t(B)
    if (m_pos == 0 & m_neg == 0){
      Q <- 0
    } else {
      if (m_pos == 0){
        Q_positive <- 0
        Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
      } else if (m_neg == 0){
        Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
        Q_negative <- 0
      } else {
        Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
        Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
      }
    }
    Q <- m_pos / m * Q_positive - m_neg / m * Q_negative
  }
}
score=function(x,y)
{
  n=nrow(y)
  yvar=apply(y,2,var)
  Tstat=n*cov(y,x)/sqrt(n*yvar*var(x))
  return(Tstat)
}
SCORE=function(x,y)
{
  n=nrow(y)
  yvar=apply(y,2,var)
  Tstat=n*cov(y,x)/sqrt(n*yvar*var(x))
  pv=1-pchisq(Tstat^2,1)
  return(cbind(pv,Tstat^2))
}
ACAT1=function(pv)
{
  ACAT=sum(tan((0.5-pv)*pi))/length(pv)
  pvalue=0.5-atan(ACAT)/pi
  return(pvalue)
}

T_Omnibus=function(x,y)
{
  K=ncol(y)
  n=nrow(y)
  yvar=apply(y,2,var)
  T=n*cov(y,x)/sqrt(n*yvar*var(x))
  U=t(T)%*%solve(cor(y))%*%T
  p=1-pchisq(U,K)
}
TATES=function(x,y)
{
  m=ncol(y)
  n=nrow(y)
  pv_j=sort(SCORE(x,y)[,1],index.return = TRUE)
  idx=pv_j$ix
  y=y[,idx]
  r=cor(y)
  rho=-0.0008-0.0023*r+0.6226*r^2+0.0149*r^3+0.1095*r^4-0.0219*r^5+0.2179*r^6
  me=numeric(m)
  for(i in 1:m)
  {
    rhoi=rho[1:i,1:i]
    eigeni=eigen(rhoi)$values
    me[i]=i-sum(ifelse(eigeni>1,eigeni,1)-1)
  }
  pv=min(me[m]*pv_j$x/me)
  return(pv)
}
pv.community.l=function(y,G)
{
  y<- as.matrix(y);  G <- as.matrix(G); 
  zscore <- sapply(1:ncol(y), FUN = function(t){
    b <- ScoreTest_SPA(genos = t(as.matrix(G)), pheno = as.numeric(y[,t]),beta.out=TRUE,beta.Cutoff=1)
    if (b$p.value == 1){
      score_values <- (y[,t] - mean(y[,t])) %*% G
      return(sign(score_values)*qnorm(1-10e-17))
    } else {
      return(b$beta/b$SEbeta)
    }
  })
  one=as.matrix(rep(1,ncol(y)))
  Si=ginv(cor(y))
  T_l=(t(one))%*%Si%*%zscore
  pv=1-pchisq((T_l)^2/((t(one))%*%Si%*%one),1)
  return(pv)
}
multi.layer.pv=function(y,G,r=ratio)
{
  A.temp <- lapply(1:nrow(y), function(t1) return((y[t1,] == 1) %*% t(y[t1,] == 1)))
  A <- Reduce("+", A.temp)
  W1=A/(nrow(y)*r)
  Sigma =cor(y)
  B=get_community(W1,A)
  Sigma_inv=ginv(Sigma)
  Sigma_inv_Kl=lapply(apply(B, 2, function(t) which(t==1)), function(z) Sigma_inv[,z])
  y.community=lapply(apply(B,2,FUN= function(z) which(z==1)), function(t) as.matrix(y[,t]) )
  z.ls=Sigma.ls=G.ls=list() # list form
  G.ls=rep(list(as.matrix(G)), length(y.community))
  Sigma.ls=rep(list(Sigma), length(y.community))
  z <- sapply(1:ncol(y), FUN = function(t){
    b <- ScoreTest_SPA(genos = t(as.matrix(G)), pheno = as.numeric(y[,t]),beta.out=TRUE,beta.Cutoff=1)
    if (b$p.value == 1){
      score_values <- (y[,t] - mean(y[,t])) %*% G
      return(sign(score_values)*qnorm(1-10e-17))
    } else {
      return(b$beta/b$SEbeta)
    }
  })
  z.ls=rep(list(as.matrix(z)), length(y.community))
  pv0=ACAT1(mapply(pv.community.l, y=y.community,G=G.ls))
  
  T_l=t(B)%*%Sigma_inv%*%z
  R=t(B)%*%Sigma_inv%*%B
  lambda=eigen(R)$values;Q=eigen(R)$vectors
  P1=diag(lambda^(-1/2))%*%t(Q); P2=Q%*%diag(lambda^(-1/2))%*%t(Q)
  Ts1=P1%*%T_l; Ts2=P2%*%T_l; 
  pv1=ACAT1(1-pchisq(Ts1^2,1));   pv2=ACAT1(1-pchisq(Ts2^2,1))
  
  W=t(B)%*%ginv(Sigma)
  U=t(W)%*%ginv(W%*%Sigma%*%t(W))%*%W
  CLC=t(z)%*%U%*%z
  pv3=1-pchisq(CLC,ncol(B))
  pv=c(pv0,pv1,pv2,pv3)
  return(pv)
}
Binarytrans=function(z) ifelse(rowSums(as.matrix(z))==0,0,1)
MLN_O.pv=function(y,G,r=ratio)
{
  A.temp <- lapply(1:nrow(y), function(t1) return((y[t1,] == 1) %*% t(y[t1,] == 1)))
  A <- Reduce("+", A.temp)
  W1=A/(nrow(y)*r)
  B=get_community(W1,A);L=ncol(B)
  n=nrow(y)
  By=matrix(0,n,L); for(l in 1:L) {By[,l]=Binarytrans(y[,which(B[,l]==1)])}
  Tscore=score(G,By)
  T=t(Tscore)%*%ginv(cor(By))%*%Tscore
  pv1=1-pchisq(T,L)
  pspa <- sapply(1:ncol(By), FUN = function(t){
    ScoreTest_SPA(genos = t(as.matrix(G)), pheno = as.numeric(By[,t]),beta.out=TRUE,beta.Cutoff=1)$p.value
  })
  pv2=ACAT1(pspa)
  pv=c(pv1,pv2)
  return(pv)
}
