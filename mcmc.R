####### main function of proposed bayesian mcmc iteration
main_mcmc = function(y,ini_m,ini_z,N,burnin,ini_xi=rep(0,dim(y)[1])){
  n = dim(y)[2]
  p = dim(y)[1]
  burnin = burnin
  N = N #num of iterations
  effsamp = N - burnin
  kmax = 20
  
  ### parameters ###
  z = ini_z
  xi = ini_xi
  M = ini_m
  phi = matrix(1,p,length(unique(z)))
  
  ### output ###
  zout = matrix(0,n,effsamp)
  xiout = matrix(0,p,effsamp)
  Mout = NULL
  phiout = NULL
  zzout = vector("list",effsamp)
  thetaout = rep(0,effsamp)
  
  ### hyperparameters ###
  lambda0 = 100                 #spike
  lambda1 = 1                   #slab
  alpha = 1                     #dirichlet
  theta = 1/p                   #bernoulli
  beta_theta = log(p)*p^(1+0.1) #beta
  lambda = 2                    #poisson
  
  
  ### MCMC ###
  for (i in 1:N){
    
    ### update z ###
    
    result_z = update_z(y,M,phi,z,xi,n,p,lambda0,lambda1,lambda,alpha)
    z = result_z$z
    M = result_z$M
    phi = result_z$phi
    
    ### update mu ###
    M = update_M(y,M,phi,z,xi,n,p,lambda0,lambda1)
    
    ### update phi ###
    phi = update_phi(M,phi,z,xi,p,lambda0,lambda1)
    
    ### update xi ###
    xi = update_xi(M,phi,z,p,lambda0,lambda1,theta)
    
    ### update theta ###
    theta = update_theta(xi,beta_theta,p)
    
    ### store MCMC outputs ###
    if(i>burnin){
      zout[,(i-burnin)] = z
      xiout[,(i-burnin)] = xi
      Mout = append(Mout,list(M))
      phiout = append(phiout,list(phi))
      thetaout[i] = theta
    }
    
    if(i%%100 == 1){
      print(paste0("Iteration ",i))
    }
  }
  
  result = list("ZSamples"=zout,
                "XiSamples"=xiout,
                "MSamples"=Mout,
                "PhiSamples"=phiout,
                "ThetaSamples"=thetaout)
  return(result)
}

####### get initial values for MCMC
####### including mean vectors m, sparsity indicators xi, and labels z
get_initial = function(y){
  sim_mc = Mclust(data=t(y),modelNames = "VII")
  ini_m = sim_mc$parameters$mean
  ini_z = sim_mc$classification
  sim = main_mcmc(y,ini_m,ini_z,500,0)
  ini_xi = sim$XiSamples[,500]
  result = list("ini_m" = ini_m,
                "ini_z" = ini_z,
                "ini_xi" = ini_xi)
  return(result)
}