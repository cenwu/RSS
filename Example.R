library(glmnet)
library(MASS)
library(regnet)
library(mltools)
library(rmutil)

################################################################## Data Generation
#used for RLSS
lam1=c(0.0618263,
      0.069378263673,
      0.078388735211,
      0.08590394335972,
      0.09757575757,
      0.109,
      0.115470137382,
      0.130700188739,
      0.146850259294,
      0.1699050356225,
      0.1892202048939,
      0.2160672336,
      0.2280923671,
      0.2491268961,
      0.25901743329,
      0.2652395027,
      0.2713290345,
      0.2884520354,
      0.296210169,
      0.301531679,
      0.3061721023,
      0.31210262,
      0.3282122163,
      0.33391954,
      0.348753189,
      0.357361525,
      0.368804628
)



n=500
p=1000
lens = length(lam)
a = matrix(0,nrow = 2*p + 1, ncol = lens)
#replicates
quant=0.5
repc=1
#predictors
k=20
rho=0
slaprop=sllaprop=smcpprop=matrix(0,p+1,repc)
slacoef=NULL
sllacoef=NULL
ladcoef=NULL
lassocoef=NULL
plassocoef=NULL
pladcoef=NULL

#variable collection
ww=NULL
#x,y collection
xx=NULL
yy=NULL

#set.seed(120)
for (q in 1:repc) {
  
  mu = rep(0, p)
  sigma = matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      sigma[i,j]=rho
    }}
  for(i in 1:nrow(sigma)){sigma[i,i] = 1}
  ep=1
  X=mvrnorm(n, mu, sigma, tol = 1e-6)
  xx=c(xx,list(X))
  beta=rep(0, p)
  bbb = sample(1:p,k)
  beta[bbb]= 0.3*c(rep(1, k))*(2*rbinom(k, 1, 0.5)-1) 
  ww=c(ww,list(bbb))
  y = 1+ X %*% beta + rt(n,2)#c(rnorm(0.3*n,0,4),rnorm(0.7*n,0,1)) -quantile(c(rnorm(0.3*10000,0,4),rnorm(0.7*10000,0,1)),probs = quant)#rnorm(n,0,1)#rlnorm(n,0,1) -quantile(rlnorm(10000,0,1),probs = quant) #c(rnorm(0.3*n,0,4),rnorm(0.7*n,0,1)) -quantile(c(rnorm(0.3*10000,0,4),rnorm(0.7*10000,0,1)),probs = quant)# rnorm(n,0,1) #rt(n,2)  #error1 = rt(n,2)-quantile(rnorm(10000),probs = quant)
  yy=c(yy,list(y))
}  


####################################################  Functions
# RLSS: Robust (LAD) LASSO with stability selection
LADL_SS <- function(x,y,t,cutoff)
{
  stab_coef=matrix(0,nrow = p+1, ncol = n)
  freq=matrix(0,p+1) #selected variable matrix
  
  for(s in 1:t) 
  {
    ### randomly split the sample in two sets
    index <- sample(n)
    i1 <- index[1:as.integer(n*0.7)]
    
    #### Normal model  
    tt <- cv.regnet(x[i1,], y[i1,], "continuous", "lasso", lamb.1 = lam1, folds = 5,NULL, clv =NULL, alpha.i=1, robust = TRUE)   
    coef <- regnet(x[i1,], y[i1,], "continuous", "lasso", lamb.1 = tt$lambda[1], NULL, clv =NULL, alpha.i=1, robust = TRUE)
    freq <- freq + matrix(abs(sign(coef$coeff)))
    stab_coef[,s] <- coef$coeff
  }
  
  prop <- freq/t # proportions
  s_result <- (prop > cutoff)*1
  cbind(s_result,prop)
}

#RLP: Robust (LAD) LASSO with assisted tuning
LADLP_SS=function(X,y,pB,SS){
  n=nrow(X)
  p=ncol(X)
  Spi=NULL
  for(tt in 1:pB){
    X_ko1=X[c(sample(nrow(X))), ]
    
    for(i in 1:lens) {
      oo=regnet(cbind(X, X_ko1), y, "continuous", "lasso", lamb.1 = lam1[i], NULL, clv =NULL, alpha.i=1, robust = TRUE)
      a[,i]=as.vector(oo$coeff)
    }
    lasso1=as.matrix(a)
    ii=1
    while (max(abs(lasso1[(p+2):(2*p+1),ii]))>0 & ii<dim(lasso1)[2]){
      ii <- ii+1 
    }
    selected_lasso = which(abs(lasso1[,ii])>0)
    Spi=c(Spi,selected_lasso[-1]-1)
  }
  freq=tabulate(Spi)/pB
  out=which(freq>SS | freq==SS)
  return(out)
}

##LSS: LASSO with stability selection
LASSO_SS <- function(x,y,t,cutoff)
{
  stab_coef=matrix(0,nrow = p + 1, ncol = n)
  freq=matrix(0,p + 1) #selected variable matrix
  
  for(s in 1:t) 
  {
    ### randomly split the sample in two sets
    index <- sample(n)
    i1 <- index[1:as.integer(n/2)]
    
    #### Normal model  
    lasso.cv_s <- cv.glmnet(x[i1,],y[i1,],alpha=1,nfolds=5)   
    lambda_s <- lasso.cv_s$lambda.min  
    lasso.fit_s <- glmnet(x[i1,],y[i1,],family="gaussian",alpha=1,nlambda=5)
    coef <- as.vector(predict(lasso.fit_s, s=lambda_s, type="coefficients"))
    freq <- freq + matrix(abs(sign(coef)))
    stab_coef[,s] <- coef
  }
  
  prop <- freq/t # proportions
  s_result <- (prop > cutoff)*1
  cbind(s_result,prop)
}

#LP: LASSO with permutation-assisted tuning
plassoc=function(X,y,pB,SS){
  n=nrow(X)
  p=ncol(X)
  #X=scale(X)
  Spi=NULL
  for(tt in 1:pB){
    X_ko1=X[c(sample(nrow(X))), ]
    b=glmnet(cbind(X, X_ko1),y)
    lasso=as.matrix(b$beta)
    ii=1
    while (max(abs(lasso[(p+1):(2*p),ii]))==0 & ii<dim(lasso)[2]){
      ii=ii+1 
    }
    selected_lasso = which(abs(lasso[,ii-1])>0)
    Spi=c(Spi,selected_lasso)
  }
  freq=tabulate(Spi)/pB
  out=which(freq>SS | freq==SS)
  return(out)
}

################################################################################# Example 

# RLSS: Robust (LAD) LASSO with stability selection
# The first element corresponds to the intercept. 
features=NULL
for (q in 1:repc) {
  slla=LADL_SS(xx[[q]],yy[[q]],80,0.9)
  sllaprop[,q]=slla[,2]
  #features=c(features,list(sllaprop[,q]))
} 

#RLP: Robust (LAD) LASSO with assisted tuning
# The first element corresponds to the intercept. 
for (q in 1:repc) {
  pB=10; SS=0.9
  select_qplasso=LADLP_SS(xx[[q]],yy[[q]],pB,SS)
  #pladcoef=c(pladcoef,list(select_qplasso))
} 


#RL: Robust (LAD) LASSO using regnet
for (q in 1:repc) {
  tt=cv.regnet(X, y, "continuous", "lasso", lamb.1 = lam1, folds = 5,NULL, clv =NULL, alpha.i=1, robust = TRUE)
  ttt=regnet(X, y, "continuous", "lasso", lamb.1 = tt$lambda[1], NULL, clv =NULL, alpha.i=1, robust = TRUE)
  ladcoef=c(ladcoef,list(which(ttt$coeff[-1]!=0)))
  #save(ladcoef,file="ladcoef.RData")
} 

#LASSO
for (q in 1:repc) {
  cc=cv.glmnet(xx[[q]],yy[[q]],alpha=1,nfolds=5) 
  lasso.fit_s <- glmnet(xx[[q]],yy[[q]],family="gaussian",alpha=1)
  ccc <- as.vector(predict(lasso.fit_s, s=cc$lambda.min, type="coefficients"))
  lassocoef=c(lassocoef,list(which(ccc[-1]!=0)))
  #save(lassocoef,file="lassocoef.RData")
} 

#LSS: LASSO with stability selection
for (q in 1:repc) {
  sla=LASSO_SS(xx[[q]],yy[[q]],80,0.9)
  slaprop[,q]=sla[,2]
  #save(slaprop,file="slaprop.RData")
} 

#LP: Lasso with permutation-assisted tuning
for (q in 1:repc) {
  pB=10; SS=0.9
  select_plasso=plassoc(xx[[q]],yy[[q]],pB,SS)
  plassocoef=c(plassocoef,list(select_plasso))
  #save(plassocoef,file="plassocoef.RData")
}