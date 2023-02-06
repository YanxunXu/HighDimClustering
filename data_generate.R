######### data sampling process
######### return a p*n dim data matrix yy and n dim cluster membership vector zz
sample_y_z = function(k,s,scen,seed,MM=1){
  ########## scen = 1 refers to equal cluster size case
  ########## scen = 2 refers to small cluster case
  ########## scen = 3 refers to the mixture of t case
  set.seed(seed)
  if (scen == 1){
    if (k == 3){
      n = 200
      p = 400
      
      #### set true mean vectors
      M1_true = rep(0,p)
      M2_true = rep(0,p)
      M3_true = rep(0,p)
      M1_true[(1:s)] = rep(3,s)
      M2_true[(1:s)] = rep(-1.5,s)
      M3_true[(1:s)] = rep(0,s)
      M_true = cbind(cbind(M1_true,M2_true),M3_true)
      z_true = rep(0,n)
      y = matrix(0,p,n)
      yy = vector(mode = "list", length = MM)
      zz = vector(mode = "list", length = MM)
      r = rep(0,n)
      
      for (j in 1:MM){
        M1_true[(1:s)] = rep(3,s)
        M2_true[(1:s)] = rep(-1.5,s)
        M3_true[(1:s)] = rep(0,s)
        for (i in 1:n){
          r[i] = runif(1)
          if (r[i]<0.3){
            y[,i] = rnorm(p,mean = M1_true,sd = rep(1,p))
            z_true[i] = 1
          }
          else if (r[i]>=0.3 && r[i]<0.6){
            cov2 = c(rep(1,s),rep(1,length(M2_true)-s))
            y[,i] = rnorm(p,mean = M2_true,sd = sqrt(cov2))
            z_true[i] = 2
          }
          else{
            y[,i] = rnorm(p,mean = M3_true,sd = rep(1,p))
            z_true[i] = 3
          }
        }
        yy[[j]] = y
        zz[[j]] = z_true
      }
    }else if (k == 5){
      n = 200
      p = 400
      #### set true mean vectors
      M1_true = rep(0,p)
      M2_true = rep(0,p)
      M3_true = rep(0,p)
      M4_true = rep(0,p)
      M5_true = rep(0,p)
      
      M1_true[(1:s)] = rep(4,s)
      M2_true[(1:s)] = rep(-4,s)
      M3_true[(1:s)] = rep(0,s)
      M4_true[(1:s)] = rep(c(-4,4),s/2)
      M5_true[(1:s)] = rep(c(1.5,-1.5),s/2)
      
      M_true = cbind(M1_true,M2_true,M3_true,M4_true,M5_true)
      
      z_true = rep(0,n)
      y = matrix(0,p,n)
      yy = vector(mode = "list", length = MM)
      zz = vector(mode = "list", length = MM)
      r = rep(0,n)
      
      for (j in 1:MM){
        M1_true[(1:s)] = rep(4,s)
        M2_true[(1:s)] = rep(-4,s)
        M3_true[(1:s)] = rep(0,s)
        M4_true[(1:s)] = rep(c(-4,4),s/2)
        M5_true[(1:s)] = rep(c(1.5,-1.5),s/2)
        
        for (i in 1:n){
          r[i] = runif(1)
          if (r[i]<0.2){
            y[,i] = rnorm(p,mean = M1_true,sd = rep(1,p))
            z_true[i] = 1
          }
          else if (r[i]>=0.2 && r[i]<0.4){
            cov2 = c(rep(1,s),rep(1,length(M2_true)-s))
            y[,i] = rnorm(p,mean = M2_true,sd = sqrt(cov2))
            z_true[i] = 2
          }
          else if (r[i]>=0.4 && r[i]<0.6){
            y[,i] = rnorm(p,mean = M3_true,sd = rep(1,p))
            z_true[i] = 3
          }
          else if (r[i]>=0.6 && r[i]<0.8){
            y[,i] = rnorm(p,mean = M4_true,sd = rep(1,p))
            z_true[i] = 4
          }
          else{
            y[,i] = rnorm(p,mean = M5_true,sd = rep(1,p))
            z_true[i] = 5
          }
        }
        yy[[j]] = y
        zz[[j]] = z_true
      }
    
  }else if (scen == 2){
    n = 200
    p = 400
    
    #### set true mean vectors
    M1_true = rep(0,p)
    M2_true = rep(0,p)
    M3_true = rep(0,p)
    M1_true[(1:s)] = rep(1,s)
    M2_true[(1:s)] = rep(-1,s)
    M3_true[(1:s)] = rep(0,s)
    M_true = cbind(cbind(M1_true,M2_true),M3_true)
    z_true = rep(0,n)
    y = matrix(0,p,n)
    yy = vector(mode = "list", length = MM)
    zz = vector(mode = "list", length = MM)
    r = rep(0,n)
    
    for (j in 1:MM){
      M1_true[(1:s)] = rep(c(5,2),s/2)
      M2_true[(1:s)] = rep(c(10,5),s/2)
      M3_true[(1:s)] = rep(c(15,2),s/2)
      for (i in 1:n){
        r[i] = runif(1)
        if (r[i]<0.02){
          y[,i] = rnorm(p,mean = M1_true,sd = rep(1,p))
          z_true[i] = 1
        }
        else if (r[i]>=0.01 && r[i]<0.5){
          cov2 = c(4*rep(1,4),rep(1,p-4))
          y[,i] = rnorm(p,mean = M2_true,sd = sqrt(cov2))
          z_true[i] = 2
        }
        else{
          y[,i] = rnorm(p,mean = M3_true,sd = rep(1,p))
          z_true[i] = 3
        }
      }
      yy[[j]] = y
      zz[[j]] = z_true
    }
    
    
    }
    
    else if (scen==3){
      M1_true = rep(0,p)
      M2_true = rep(0,p)
      M3_true = rep(0,p)
      M1_true[c(1,51,101,151,201,251,301,351)] = c(5,2,5,2,5,2,5,2)
      M2_true[c(1,51,101,151,201,251,301,351)] = c(10,5,10,5,10,5,10,5)
      M3_true[c(1,51,101,151,201,251,301,351)] = c(15,2,15,2,15,2,15,2)
      M_true = cbind(cbind(M1_true,M2_true),M3_true)
      z_true = rep(0,n)
      y = matrix(0,p,n)
      r = rep(0,n)
      yy = vector(mode = "list",length = MM)
      zz = vector(mode = "list",length = MM)
      MM = MM
      for (j in 1:1){
        for (i in 1:n){
          r[i] = runif(1)
          if (r[i]<0.2){
            cov1 = diag(length(M1_true))
            y[,i] = mvtnorm::rmvt(1, df = 5,delta = M1_true,sigma = cov1)
            z_true[i] = 1
          }
          else if (r[i]>=0.2 && r[i]<0.6){
            cov2 = diag(length(M2_true))
            cov2[1,1] = 4
            cov2[51,51] = 4
            cov2[101,101] = 4
            cov2[151,151] = 4
            y[,i] = mvtnorm::rmvt(1,df = 5,delta = M2_true,sigma = cov2)
            z_true[i] = 2
          }
          else{
            cov3 = diag(length(M1_true))
            y[,i] = mvtnorm::rmvt(1,df = 5,delta = M3_true,sigma = cov3)
            z_true[i] = 3
          }
        }
        yy[[j]] = y
        zz[[j]] = z_true
      }
      
    }
    
  }
  return(list("y"=y,"z_true"=z_true))
}



