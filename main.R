#### This file is a demo of implementation of proposed Bayesian 
#### method in simulation setup with k=3, s=12 in Scenario I    


## First specify the packages of interest
packages = c("extraDistr", "GIGrvg", "ggplot2", 
             "gridExtra","mclust", "sparcl","coda",
             "stats","aricode","reshape2","cluster")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

source("mcmc.R")
source("mcmc_update.R")
source("data_generate.R")

n_exp = 1        #number of experiment replicates
Niter = 5000     #number of mcmc iterations
burnin = 1000    #number of burn-in iterations

data = vector("list",n_exp)     #store p*n dimensional data with true label
sim_b = vector("list",n_exp)    #store results of bayesian method
ari_b = rep(0,n_exp)            #store results of ARI

k = 3     #number of true clusters
s = 12    #number of non-zero coordinates

for (i in 1:n_exp){
  print(paste0("Experiment ",i,":"))
  seed = i
  ########## sample data including labels
  data[[i]] = sample_y_z(k=k,s=s,scen=1,seed)
  y_true = data[[i]]$y
  
  ############# proposed Bayesian method ##################
  set.seed(seed)
  ## get initial values for MCMC updates
  print("Getting initial values ... usually takes a minute or two")
  invisible(capture.output(ini = get_initial(y_true))) 
  ini_m = ini$ini_m
  ini_xi = ini$ini_xi
  ini_z = ini$ini_z
  print("done")
  ## run MCMC
  print("Main MCMC iterations ...")
  sim_b[[i]] = main_mcmc(y_true,ini_m,ini_z,Niter,burnin,ini_xi)
  print("done")

}


####### calculate ARI
for (m in 1:n_exp){
  ##### get posterior mode of z
  z_unique = unique(sim_b[[m]]$ZSamples,MARGIN=2)
  z_pos_freq = rep(0,ncol(z_unique) )
  for (i in 1:ncol(z_unique)){
    for (j in 1:(Niter-burnin)){
      if (sum(sim_b[[m]]$ZSamples[,j]-z_unique[,i])==0){z_pos_freq[i] = z_pos_freq[i] + 1}
    }
  }
  ##### ari
  ari_b[m] = adjustedRandIndex(z_unique[,which.max(z_pos_freq)],data[[m]]$z_true)
}





####### scatter plots of truth and clustering result
i = 1
plt1 = vector(mode = "list", length = 2)
plt1[[1]] = qplot(data[[i]]$y[1,],data[[i]]$y[2,],colour = factor(data[[i]]$z_true)) + xlab("") + ylab("") + ggtitle("Truth") + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
plt1[[2]] = qplot(data[[i]]$y[1,],data[[i]]$y[2,],colour = factor(sim_b[[i]]$ZSamples[,ncol(sim_b[[i]]$ZSamples)])) + xlab("") + ylab("") + ggtitle("Bayesian") + theme(legend.position="none") + theme(plot.title = element_text(hjust = 0.5))
grid.arrange(plt1[[1]], plt1[[2]])


