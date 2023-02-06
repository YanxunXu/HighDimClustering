####### mcmc update of mean vectors mu
update_M = function(y,M,phi,z,xi,n,p,lambda0,lambda1){
  for (c in 1:length(unique(z))){
    y_c = as.matrix(y[,z==c])
    var = 1/(sum(z==c)+(lambda0^2*(1-xi)+lambda1^2*xi)/(phi[,c]))
    M[,c] = rnorm(p,mean = rowSums(y_c) * var, sd = sqrt(var))
  }
  return(M)
}

####### mcmc update of auxiliary variable phi
update_phi = function(M,phi,z,xi,p,lambda0,lambda1){
  for (c in 1:length(unique(z))){
    for (j in 1:p){
      ### generalized inverse gaussian
      phi[j,c] = rgig(1,lambda = 0.5,chi = M[j,c]^2*(lambda0^2*(1-xi[j])+lambda1^2*xi[j]),psi = 1)
    }
  }
  return(phi)
}

####### mcmc update of sparsity indicators xi
update_xi = function(M,phi,z,p,lambda0,lambda1,theta){
  t = length(unique(z))
  xi_prob = 1/(1 + (lambda0/lambda1)^t * matrixStats::rowProds(exp(-0.5*(lambda0^2-lambda1^2)*M^2/phi)) * (1-theta)/theta )
  xi = rbinom(p,1,prob = xi_prob)
  return(xi)
}

####### mcmc update of clustering memberships z
update_z = function(y,M,phi,z,xi,n,p,lambda0,lambda1,lambda,alpha){
  #### constants that will be used in MFM model
  log_t = log(seq(1,20))
  log_n = log(n)
  log_dpois = log(lambda/(seq(1,20)+1) )
  #### MFM model
  for (j in 1:n){
    zz = z[-j]
    t = length(unique(zz))
    idx = z[j]
    
    if (idx %in% zz){
      
      new_phi = rexp(p,rate = 0.5)
      sd_new_mu = sqrt(new_phi)/(lambda0*(1-xi)+lambda1*xi)
      new_mu = rnorm(p, mean = 0, sd = sd_new_mu)
      M = as.matrix(M)
      z_prob = log((table(zz)) + alpha) - 0.5 * colSums((y[,j]-M)^2)
      
      
      z_new = -0.5*sum((y[,j]-new_mu)^2) - log_n + 2*log_t[t] + log_dpois[t] 
      z_prob = c(z_prob,z_new)
      bb = max(z_prob)
      tmp = exp(z_prob-bb)
      z_prob = tmp/sum(tmp)
      z[j] = sample(1:(t+1),size = 1,prob = z_prob)
      
      if(z[j]==(t+1)){
        M = cbind(M, new_mu)
        phi = cbind(phi, new_phi)
      }
      
    }else{
      new_phi = rexp(p,rate = 0.5)
      sd_new_mu = sqrt(new_phi)/(lambda0*(1-xi)+lambda1*xi)
      new_mu = rnorm(p, mean = 0, sd = sd_new_mu)
      #### relabeling
      pos = which(z>idx)
      z[pos] = z[pos] - 1
      M = M[,-idx]
      phi = phi[,-idx]
      M = as.matrix(M)
      z_prob = log((table(zz)) + alpha) - 0.5 * colSums((y[,j]-M)^2)
      z_new = -0.5*sum((y[,j]-new_mu)^2) - log_n + 2*log_t[t] + log_dpois[t] 
      z_prob = c(z_prob,z_new)
      
      bb = max(z_prob)
      tmp = exp(z_prob-bb)
      z_prob = tmp/sum(tmp)
      z[j] = sample(1:(t+1),size = 1,prob = z_prob)
      
      if(z[j]==(t+1)){
        M = cbind(M, new_mu)
        phi = cbind(phi, new_phi)
      }
    }
    
    M = as.matrix(M)
    phi = as.matrix(phi)
  }
  
  
  result = list("z" = z,"M" = M,"phi" = phi)
  return(result)
}

####### mcmc update of hyperparameter theta
update_theta = function(xi,beta_theta,p){
  theta = rbeta(1,1+sum(xi),beta_theta+p-sum(xi))
  return(theta)
}


