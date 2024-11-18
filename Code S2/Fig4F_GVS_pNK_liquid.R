
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(R.matlab)
rm(list=ls())

data_total<-readMat("nk_data_230622.mat")
MFI<-data_total$mfi.all

cyto<-lapply(data_total[c('pNK.025.1','pNK.05.1','pNK.1.1')], as.numeric)

leuk_ind=2:13;
x=MFI[leuk_ind,1:5]-1
x[x<0]=0
# x=cbind(x,rowSums(x[1:4]))

ligand_min=apply(x, 2, min)

x_zero_min=sweep(x, 2, ligand_min)
x_normalized=sweep(x_zero_min, 2, apply(x_zero_min,2,max),FUN="/")
x_ten=10*x_normalized
x_ten_threefold= matrix(rep(t(x_ten),3),ncol=5,byrow=TRUE)
y=c(cyto$pNK.025.1[leuk_ind],cyto$pNK.05.1[leuk_ind],cyto$pNK.1.1[leuk_ind])

# y=(y-min(y))/(max(y)-min(y))
x=x_ten_threefold

tau=100
s0=1
n_iter=1e6

draw_beta = draw_gamma = draw_beta_accept = draw_gamma_accept = array(0,c(n_iter,6))
draw_c_vec = draw_c_vec_accept = array(0,c(n_iter,2))
draw_K = draw_sigma = draw_K_accept = draw_sigma_accept = rep(0,n_iter)
# draw_sigma = draw_sigma_accept = rep(0,n_iter)

log_likelihood<-function(x, y, beta, gamma, c_vec, sigma, K){
  
  scaling=rep(c(1,c_vec), each=12)
  x_K = x^2/ (x^2+K^2)
  x = cbind(rep(1,36),x_K)
  err = y-scaling*(as.matrix(x %*% (beta*gamma)))
  log_likelihood=dnorm(err, mean=0, sd=sigma, log=T)
  
  return(sum(log_likelihood))
}

# log_likelihood_with_data<-log_likelihood(x_ten_threefold,y)

log_prior<-function(beta_j,gamma_j){
  if (gamma_j==1) #gamma_j : variable is included in the model
  {
    if (beta_j>=0){
      return(dnorm(beta_j, mean=0,sd=tau,log=T)+log(2))
    }else{
      return (-Inf)
    }
  }
  else
  {
    if (beta_j>=0){
      return(dnorm(beta_j, mean=0,sd=s0,log=T)+log(2))
    }else{
      return (-Inf)
    }
  }
}


##### Initiation
# coef_init=c(2,2,2,2,2,2,2,2,2,2)
# gamma_init=c(0,0,0,0,0)
# 
# beta_init=c(coef_init, gamma_init)
# beta=beta_init  

beta = 10*c(2,2,2,2,2,2) ; gamma = c(1,1,1,1,1,1) ; c_vec = c(1.5,2.0)
K = 5 ;
sigma = 15

stepsize_beta = 2*c(0.1,0.1,0.1,0.1,0.1,0.1)
stepsize_c_vec = 0.3*c(0.01,0.01)

set.seed(3)

##### Gibbs Variable Selection step 
for (i_iter in 1:n_iter){
  
  for (j in 1:6){
    beta_q=beta
    beta_q[j]=rnorm(1,beta[j],stepsize_beta[j])
    if (beta_q[j]>0){
      logProb = log_likelihood(x, y, beta_q, gamma, c_vec,  sigma ,K) - log_likelihood(x, y, beta, gamma, c_vec, sigma ,K)
      logProb = logProb + log_prior(beta_q[j],gamma[j]) - log_prior(beta[j],gamma[j]) 
      if (runif(1)<exp(logProb)){
        draw_beta_accept[i_iter, j]=1
        beta=beta_q
      }
    }
  } # for (j)
  
  ### c_k & K update ###
  for (j in 1:2){
    c_vec_q = c_vec
    c_vec_q[j] = rnorm(1,c_vec[j],stepsize_c_vec[j])
    logProb = log_likelihood(x, y, beta, gamma, c_vec_q,  sigma ,K) - log_likelihood(x, y, beta, gamma, c_vec, sigma, K)
    if (runif(1)<exp(logProb)){
      draw_c_vec_accept[i_iter, j]=1
      c_vec=c_vec_q
    }
  }
  c_vec=pmax(0, c_vec)
  c_vec=pmin(10, c_vec)
  
  K_q = K + 0.1*(2*(runif(1)>0.5) - 1)
  if (K_q >=1){
    logProb = log_likelihood(x, y, beta, gamma, c_vec,  sigma, K_q) - log_likelihood(x, y, beta, gamma, c_vec, sigma, K)
    if (runif(1)<exp(logProb)){
      draw_K_accept[i_iter]=1
      K = K_q
    }
  }
  K=max(1, K)
  K=min(100,K)
  
  # ### Gamma update ###
  for (j in 2:6){
    gamma_1=gamma; gamma_1[j]=1
    gamma_0=gamma; gamma_0[j]=0
    
    
    logProb = log_likelihood(x, y, beta, gamma_1, c_vec, sigma, K) - log_likelihood(x, y, beta, gamma_0, c_vec, sigma,K)
    logProb = logProb + log_prior(beta[j],1) - log_prior(beta[j],0)
    
    prob=1/(1/exp(logProb)+1)
    gamma[j]=(runif(1)<prob)
    
  }
  
  ### Sigma update ###
  scaling=rep(c(1,c_vec), each=12)
  # x_K = x / (x+K)
  x_aug = cbind(rep(1,36),x)
  x_aug_K = x_aug^2/(x_aug^2+K^2)
  
  err = y-scaling*(as.matrix(x_aug_K %*% (beta*gamma)))
  
  sigma=sqrt(1/rgamma(1,1+36/2,1+0.5*sum(err^2)))
  
  if (is.na(sigma)) print(i_iter) ; 
  
  ### Update summary
  draw_beta[i_iter,]=beta
  draw_gamma[i_iter,]=gamma
  draw_c_vec[i_iter,]=c_vec
  draw_K[i_iter]=K
  draw_sigma[i_iter]=sigma
  if (i_iter %% 5000==0){print(i_iter)}
  
} # for (i_iter)



# par(mfrow=c(2,3))
# # plot(draw_gamma[seq(from=5000000,to=10000000,by=1000),1], type="l")
# plot(draw_gamma[seq(from=1,to=n_iter,by=100),2], type="l")
# plot(draw_gamma[seq(from=1,to=n_iter,by=100),3], type="l")
# plot(draw_gamma[seq(from=1,to=n_iter,by=100),4], type="l")
# plot(draw_gamma[seq(from=1,to=n_iter,by=100),5], type="l")
# plot(draw_gamma[seq(from=1,to=n_iter,by=100),6], type="l")


par(mfrow=c(2,3))
# plot(draw_beta[seq(from=5000000,to=10000000,by=1000),1], type="l")
plot(draw_beta[seq(from=1,to=n_iter,by=100),1], type="l")
plot(draw_beta[seq(from=1,to=n_iter,by=100),2]*draw_gamma[seq(from=1,to=n_iter,by=100),2], type="l")
plot(draw_beta[seq(from=1,to=n_iter,by=100),3]*draw_gamma[seq(from=1,to=n_iter,by=100),3], type="l")
plot(draw_beta[seq(from=1,to=n_iter,by=100),4]*draw_gamma[seq(from=1,to=n_iter,by=100),4], type="l")
plot(draw_beta[seq(from=1,to=n_iter,by=100),5]*draw_gamma[seq(from=1,to=n_iter,by=100),5], type="l")
plot(draw_beta[seq(from=1,to=n_iter,by=100),6]*draw_gamma[seq(from=1,to=n_iter,by=100),6], type="l")

par(mfrow=c(2,3))
# plot(draw_gamma[seq(from=5000000,to=10000000,by=1000),1], type="l")
hist(draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=100),2])
hist(draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=100),3])
hist(draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=100),4])
hist(draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=100),5])
hist(draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=100),6])


par(mfrow=c(2,3))
# plot(draw_beta[seq(from=5000000,to=10000000,by=1000),1], type="l")
hist(draw_beta[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),1],freq=F)
hist(draw_beta[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),2]*draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),2],freq=F)
hist(draw_beta[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),3]*draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),3],freq=F)
hist(draw_beta[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),4]*draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),4],freq=F)
hist(draw_beta[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),5]*draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),5],freq=F)
hist(draw_beta[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),6]*draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),6],freq=F)

round( apply( draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),]==0, 2, mean ) , 3 )
round( apply( draw_gamma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),]==1, 2, mean ) , 3 )

# par(mfrow=c(2,3))
# plot(draw_c_vec[seq(from=1,to=n_iter,by=1),1])
# plot(draw_c_vec[seq(from=1,to=n_iter,by=1),2])
# plot(draw_sigma[seq(from=1,to=n_iter,by=1)])

par(mfrow=c(2,3))
hist(draw_c_vec[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),1])
hist(draw_c_vec[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),2])
hist(draw_sigma[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1)])
hist(draw_K[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1)])


mean( draw_c_vec[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),1] < draw_c_vec[seq(from=round(n_iter*0.2)+1,to=n_iter,by=1),2] )

library(R.matlab)
writeMat("draw_beta_liquid_230622.mat",draw_c_vec=draw_c_vec,draw_beta=draw_beta, draw_gamma=draw_gamma,draw_sigma=draw_sigma, draw_K=draw_K)
writeMat("draw_beta_accept_liquid_230622.mat",draw_c_vec_accept=draw_c_vec_accept,draw_beta_accept=draw_beta_accept,draw_gamma_accept=draw_gamma_accept, draw_sigma_accept=draw_sigma_accept,draw_K_accept=draw_K_accept)
