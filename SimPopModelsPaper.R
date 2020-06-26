
#Function for obtaining alpha and beta ofa Beta-distribution from given mean (mu) and variance (var)
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}



#Harvest strategy 1 ‘ordinary’ (with threshold of maximum sex ratio): 
#Number of adult males is specified by the harvest rate (H1[6]), but constrained so that there is
#a maximum number of adult famales to each adult male.
#If mRatio=20, at least 5% of post-harvest adults should be male.
#Harvest strategy starts from first year (2018)

SimPop18Ktot_h  <- function(P=11, T_Adf=100,K=Ktot,mRatio=10,phi3_m=phi3.m,phi3_sd=phi3.sd, f_m=f.m, f_sd=f.sd, phi1_m = phi1.m, phi1_sd=phi1.sd,H1=h1,hadf.max=hadf_max,hadf.m=hadf_m,Nmean=N.mean,Nsd=N.sd){
  # P: Number of time steps to simulate
  # T_Adf: threshold of adult females. Do not hunt additional females if N.adf<= T_Adf
  # H (h): scenario of harvest (harvest rates) for P time steps
  # mRatio: set the upper limit of number of adult females per adult male
  # K: In order to stabilize the population size, harvest rate of adult females (H1[5]) was for each year
  # deterimined by the total population size compared to K (carrying capacity)
  # hadf.max and hadf.m: paramters that decides the rate of change of adult female harvest in relation to population size and K
  # Demographic rates are given as mean (_m) and standard deviation (_sd)
  # PHI3: Adult annual survival probability 
  # PHI1: Juvenile summer survival probability
  # f: fertility rate
  # Nmean: mean pre-harvest population size 
  # Nsd: sd of pre-harvest population size

  ############################################################
  # Define the priors for the parameters
  ############################################################ 
  
  ## POPULATION VECTORS
  N <- matrix(ncol=P, nrow=6)               ## Pre harvest pop. vector. No monitoring
  X <- matrix(ncol=P, nrow=6)               ## Post harvest pop. vector. No monitoring
  H <- matrix(ncol=P, nrow=6)               ## Harvest numbers
  
  N_tot <- matrix(ncol=P, nrow=1)
  X_tot <- matrix(ncol=P, nrow=1)
  H_tot <- matrix(ncol=P, nrow=1)
  HU_tot <- matrix(ncol=P, nrow=1)
  HAd_tot <- matrix(ncol=P, nrow=1)
  
  # Initial pre-harvest population sizes
  N[1,1] <- max(round(rnorm(1,Nmean[1], Nsd[1]),0),10)      # calves females
  N[2,1] <- max(round(rnorm(1,Nmean[2], Nsd[2]),0),10)       # calves males
  N[3,1] <- max(round(rnorm(1,Nmean[3], Nsd[3]),0),10)        # yearling females
  N[4,1] <- max(round(rnorm(1,Nmean[4], Nsd[4]),0),10)        # yearling males
  N[5,1] <- max(round(rnorm(1,Nmean[5], Nsd[5]),0),10)      # adult females
  N[6,1] <- max(round(rnorm(1,Nmean[6], Nsd[6]),0),10)	    # adult males
  
  ## DEMOGRAPHIC PARAMETERS
  # fecundity;   
  f <- matrix(ncol=P, nrow=1) 
  # Juvenile summer survival
  PHI1 <- matrix(ncol=P, nrow=1)
  
  phi3_var=phi3_sd*phi3_sd
  dp3<-estBetaParams(phi3_m,phi3_var)  
  PHI3<-rbeta(1,dp3$alpha,dp3$beta)
  
  phi1_var=phi1_sd*phi1_sd
  f_var=f_sd*f_sd  
  
  for(i in 1:P){
    
    dp1<-estBetaParams(phi1_m,phi1_var)
    df<-estBetaParams(f_m,f_var)
    
    PHI1[i]<-rbeta(1,dp1$alpha,dp1$beta)
    f[i]<-rbeta(1,df$alpha,df$beta)
    
  }
  
  h<-H1
  
  #############################
  # SYSTEM PROCESS
  #############################
  
  for (t in 1:(P-1)){ 
    ###########################################################
    # STATE PROCESS;
    # PRE-HARVEST POPULATION VECTORS IN T+1  
    hadf=hadf.max-hadf.m*(K-(N[1,t]+N[2,t]+N[3,t]+N[4,t]+N[5,t]+N[6,t]))/K  
    H5<-round(N[5,t]*hadf)   #B?r legge p? variasjon rundt hadf?
    hlow5<-ifelse(N[5,t]-T_Adf>0,N[5,t]-T_Adf,0)
    H[5,t] <- ifelse(N[5,t]-H5>T_Adf,H5,hlow5)
    T_Adm<-round((N[5,t]-H[5,t])/mRatio)   #mratio=20 , dvs 5 %
    
    #B?r vi legge inn minimum number of males required?
    #H[6,t] <- ifelse(N[6,t]-T_Adm>1,N[6,t]-T_Adm,0)
    H[6,t] <- ifelse(N[6,t]*(1-h[6])>T_Adm,round(N[6,t]*h[6]),max((N[6,t]-T_Adm),0))
    
    N[3,t+1] <- rbinom(1, round(N[1,t]*(1-h[1])), PHI3)
    N[4,t+1] <- rbinom(1, round(N[2,t]*(1-h[2])), PHI3)
    
    N[5,t+1] <- rbinom(1, round(N[3,t]*(1-h[3])+N[5,t]-H[5,t]), PHI3)
    N[6,t+1] <- rbinom(1, round(N[4,t]*(1-h[4])+N[6,t]-H[6,t]), PHI3)
    
    
    N[1,t+1] <- rbinom(1, N[5, t+1], PHI1[t]*f[t]/2)   
    N[2,t+1] <- rbinom(1, N[5, t+1], PHI1[t]*f[t]/2)   
    
    H[1,t]  <-  round(N[1,t]*h[1])
    H[2,t]  <-  round(N[2,t]*h[2])
    H[3,t]  <-  round(N[3,t]*h[3])
    H[4,t]  <-  round(N[4,t]*h[4])
    
  }
  
  hadf=hadf.max-hadf.m*(K-(N[1,P]+N[2,P]+N[3,P]+N[4,P]+N[5,P]+N[6,P]))/K  
  H5<-round(N[5,P]*hadf) 
  hlow5<-ifelse(N[5,P]-T_Adf>0,N[5,P]-T_Adf,0)
  H[5,P] <- ifelse(N[5,P]-H5>T_Adf,H5,hlow5)
  
  T_Adm<-round((N[5,P]-H[5,P])/mRatio)   #mratio=20 , dvs 5 %    
  #H[6,P] <- ifelse(N[6,P]-T_Adm>1,N[6,P]-T_Adm,0)
  #H[6,P] <- ifelse(N[6,P]-H[6,P]>T_Adm,H[6,P],max((N[6,P]-T_Adm),0))
  H[6,P] <- ifelse(N[6,P]*(1-h[6])>T_Adm,round(N[6,P]*h[6]),max((N[6,P]-T_Adm),0))
  
  H[1,P]  <-  round(N[1,P]*h[1])
  H[2,P]  <-  round(N[2,P]*h[2])
  H[3,P]  <-  round(N[3,P]*h[3])
  H[4,P]  <-  round(N[4,P]*h[4])
  
  for (t in 1:P){ 
    
    #############################################################
    # POST-HARVEST POPULATION VECTORS IN T+1
    X[1:4,t] <- (N[1:4,t]-H[1:4,t])
    X[5:6,t] <- (N[5:6,t]-H[5:6,t])
    
    
    #############################################################
    # DERIVED HARVEST NUMBERS  
    #H[,t] <- round(N[,t]*h[,t])
    
    X_tot[t] <- sum(X[,t])  # POST-HARVEST POPULATION size
    N_tot[t] <- sum(N[,t])  # summing up population vector to population size	
    H_tot[t] <- sum(H[,t])
    HU_tot[t] <- sum(H[3:6,t])
    HAd_tot[t] <- sum(H[5:6,t])
    
  }
  
  out <- list(N, X, H, N_tot, X_tot, H_tot,HU_tot, HAd_tot,
              f, PHI1, PHI3)
  names(out) <- c("N", "X","H", "N_tot","X_tot","H_tot","HU_tot",
                  "HAd_tot",
                  "f", "phi1", "phi3")
  out
  
}



#Harvest strategy 2 ‘proactive’ (with operational sex ratio): 
#Number of adult males harvested (Hadm[t]) in year t is set to obtain a specified operational 
#sex ratio (SR = m:f = 1:mRatio) after harvest. The aim for number of adult males 
#after harvest are then determined by mRatio and the number of post-harvest number of adult females.
#Harvest strategy starts from first year (2018)

SimPop18RKtot_h  <- function(P=11, T_Adf=100,K=Ktot,mRatio=10,phi3_m=phi3.m,phi3_sd=phi3.sd, f_m=f.m, f_sd=f.sd, phi1_m = phi1.m, phi1_sd=phi1.sd,H1=h1,hadf.max=hadf_max,hadf.m=hadf_m,Nmean=N.mean,Nsd=N.sd){
  # P: Number of time steps to simulate
  # H (h): scenario of harvest (harvest rates) for P time steps
  # mRatio: the operational sex ratio, the number of adult females per adult male to be obtained after harvest
  # T_Adf: threshold of adult females. Do not hunt additional females if N.adf<= T_Adf
  # K: In order to stabilize the population size, harvest rate of adult females (H1[5]) was for each year
  # deterimined by the total population size compared to K (carrying capacity)
  # hadf.max and hadf.m: paramters that decides the rate of change of adult female harvest in relation to population size and K
  
  # Demographic rates are given as mean (_m) and standard deviation (_sd)
  # PHI3: Adult annual survival probability 
  # PHI1: Juvenile summer survival probability
  # f: fertility rate
  # Nmean: mean pre-harvest population size 
  # Nsd: sd of pre-harvest population size
  
  ############################################################
  # Define the priors for the parameters
  ############################################################ 
  
  ## POPULATION VECTORS
  N <- matrix(ncol=P, nrow=6)               ## Pre harvest pop. vector. No monitoring
  X <- matrix(ncol=P, nrow=6)               ## Post harvest pop. vector. No monitoring
  H <- matrix(ncol=P, nrow=6)               ## Harvest numbers
  
  N_tot <- matrix(ncol=P, nrow=1)
  X_tot <- matrix(ncol=P, nrow=1)
  H_tot <- matrix(ncol=P, nrow=1)
  HU_tot <- matrix(ncol=P, nrow=1)
  HAd_tot <- matrix(ncol=P, nrow=1)
  
  # Initial pre-harvest population sizes
  N[1,1] <- max(round(rnorm(1,Nmean[1], Nsd[1]),0),10)      # calves females
  N[2,1] <- max(round(rnorm(1,Nmean[2], Nsd[2]),0),10)      # calves males
  N[3,1] <- max(round(rnorm(1,Nmean[3], Nsd[3]),0),10)      # yearling females
  N[4,1] <- max(round(rnorm(1,Nmean[4], Nsd[4]),0),10)      # yearling males
  N[5,1] <- max(round(rnorm(1,Nmean[5], Nsd[5]),0),10)      # adult females
  N[6,1] <- max(round(rnorm(1,Nmean[6], Nsd[6]),0),10)	    # adult males
  
  ## DEMOGRAPHIC PARAMETERS
  # fecundity;   
  f <- matrix(ncol=P, nrow=1) 
  # Juvenile summer survival
  PHI1 <- matrix(ncol=P, nrow=1)
  
  phi3_var=phi3_sd*phi3_sd
  dp3<-estBetaParams(phi3_m,phi3_var)  
  PHI3<-rbeta(1,dp3$alpha,dp3$beta)
  
  phi1_var=phi1_sd*phi1_sd
  f_var=f_sd*f_sd  
  
  for(i in 1:P){
    
    dp1<-estBetaParams(phi1_m,phi1_var)
    df<-estBetaParams(f_m,f_var)
    
    PHI1[i]<-rbeta(1,dp1$alpha,dp1$beta)
    f[i]<-rbeta(1,df$alpha,df$beta)
    
  }
  
  h<-H1
  
  #############################
  # SYSTEM PROCESS
  # STATE PROCESS;
  # PRE-HARVEST POPULATION VECTORS IN T+1  
  #############################
 
  for (t in 1:(P-1)){ 
    ###########################################################
    # STATE PROCESS;
    # PRE-HARVEST POPULATION VECTORS IN T+1  
    hadf=hadf.max-hadf.m*(K-(N[1,t]+N[2,t]+N[3,t]+N[4,t]+N[5,t]+N[6,t]))/K  
    H5<-round(N[5,t]*hadf)   
    hlow5<-ifelse(N[5,t]-T_Adf>0,N[5,t]-T_Adf,0)
    H[5,t] <- ifelse(N[5,t]-H5>T_Adf,H5,hlow5)
    
    T_Adm<-round((N[5,t]-H[5,t])/mRatio)   #mratio=20 -> 5 %
    
    H[6,t] <- ifelse(N[6,t]-T_Adm>1,N[6,t]-T_Adm,0)
    
    N[3,t+1] <- rbinom(1, round(N[1,t]*(1-h[1])), PHI3)
    N[4,t+1] <- rbinom(1, round(N[2,t]*(1-h[2])), PHI3)
    
    N[5,t+1] <- rbinom(1, round(N[3,t]*(1-h[3])+N[5,t]-H[5,t]), PHI3)
    N[6,t+1] <- rbinom(1, round(N[4,t]*(1-h[4])+N[6,t]-H[6,t]), PHI3)
    
    
    N[1,t+1] <- rbinom(1, N[5, t+1], PHI1[t]*f[t]/2)   
    N[2,t+1] <- rbinom(1, N[5, t+1], PHI1[t]*f[t]/2)   
    
    H[1,t]  <-  round(N[1,t]*h[1])
    H[2,t]  <-  round(N[2,t]*h[2])
    H[3,t]  <-  round(N[3,t]*h[3])
    H[4,t]  <-  round(N[4,t]*h[4])
    
  }
  
  hadf=hadf.max-hadf.m*(K-(N[1,t]+N[2,t]+N[3,t]+N[4,t]+N[5,t]+N[6,t]))/K   
  H5<-round(N[5,P]*hadf) 
  hlow5<-ifelse(N[5,P]-T_Adf>0,N[5,P]-T_Adf,0)
  H[5,P] <- ifelse(N[5,P]-H5>T_Adf,H5,hlow5)
  
  T_Adm<-round((N[5,P]-H[5,P])/mRatio)   
  
  H[6,P] <- ifelse(N[6,P]-T_Adm>1,N[6,P]-T_Adm,0)
  
  H[1,P]  <-  round(N[1,P]*h[1])
  H[2,P]  <-  round(N[2,P]*h[2])
  H[3,P]  <-  round(N[3,P]*h[3])
  H[4,P]  <-  round(N[4,P]*h[4])
  
  for (t in 1:P){ 
    
    #############################################################
    # POST-HARVEST POPULATION VECTORS IN T+1
    X[1:4,t] <- (N[1:4,t]-H[1:4,t])
    X[5:6,t] <- (N[5:6,t]-H[5:6,t])
    
    
    #############################################################
    # DERIVED HARVEST NUMBERS  
    #H[,t] <- round(N[,t]*h[,t])
    
    X_tot[t] <- sum(X[,t])  # POST-HARVEST POPULATION size
    N_tot[t] <- sum(N[,t])  # summing up population vector to population size	
    H_tot[t] <- sum(H[,t])
    HU_tot[t] <- sum(H[3:6,t])
    HAd_tot[t] <- sum(H[5:6,t])
    
  }
  
  out <- list(N, X, H, N_tot, X_tot, H_tot,HU_tot, HAd_tot,
              f, PHI1, PHI3)
  names(out) <- c("N", "X","H", "N_tot","X_tot","H_tot","HU_tot",
                  "HAd_tot",
                  "f", "phi1", "phi3")
  out
  
}

#Harvest strategy 2 ‘proactive’ (with operational sex ratio).
#Same as above ("SimPop18RKtot_h"), but number of harvested calves is determined as a proportion of harvested adult females
SimPop18RKtot_hcalf  <- function(P=11, T_Adf=100,K=Ktot,mRatio=10,phi3_m=phi3.m,phi3_sd=phi3.sd, f_m=f.m, f_sd=f.sd, phi1_m = phi1.m, phi1_sd=phi1.sd,H1=h1,hc=hcalf,hadf.max=hadf_max,hadf.m=hadf_m,Nmean=N.mean,Nsd=N.sd){
  # P: Number of time steps to simulate
  # H (h): scenario of harvest (harvest rates) for P time steps
  # mRatio: the operational sex ratio, the number of adult females per adult male to be obtained after harvest
  # T_Adf: threshold of adult females. Do not hunt additional females if N.adf<= T_Adf
  # K: In order to stabilize the population size, harvest rate of adult females (H1[5]) was for each year
  # deterimined by the total population size compared to K (carrying capacity)
  # hadf.max and hadf.m: paramters that decides the rate of change of adult female harvest in relation to population size and K
  # hc=hcalf: proportion of adult females harvested for which there is also harvested a calf
  # Demographic rates are given as mean (_m) and standard deviation (_sd)
  # PHI3: Adult annual survival probability 
  # PHI1: Juvenile summer survival probability
  # f: fertility rate
  # N.mean: mean pre-harvest population size 
  # N.sd: sd of pre-harvest population size
  # h: scenario of harvest rates for T time steps
  # Do not hunt additional adult males if N.adm<= minAdm
  
  ############################################################
  # Define the priors for the parameters
  ############################################################ 
  
  ## POPULATION VECTORS
  N <- matrix(ncol=P, nrow=6)               ## Pre harvest pop. vector. No monitoring
  X <- matrix(ncol=P, nrow=6)               ## Post harvest pop. vector. No monitoring
  H <- matrix(ncol=P, nrow=6)               ## Harvest numbers
  
  N_tot <- matrix(ncol=P, nrow=1)
  X_tot <- matrix(ncol=P, nrow=1)
  H_tot <- matrix(ncol=P, nrow=1)
  HU_tot <- matrix(ncol=P, nrow=1)
  HAd_tot <- matrix(ncol=P, nrow=1)
  
  # Initial pre-harvest population sizes
  N[1,1] <- max(round(rnorm(1,Nmean[1], Nsd[1]),0),10)      # calves females
  N[2,1] <- max(round(rnorm(1,Nmean[2], Nsd[2]),0),10)      # calves males
  N[3,1] <- max(round(rnorm(1,Nmean[3], Nsd[3]),0),10)      # yearling females
  N[4,1] <- max(round(rnorm(1,Nmean[4], Nsd[4]),0),10)      # yearling males
  N[5,1] <- max(round(rnorm(1,Nmean[5], Nsd[5]),0),10)      # adult females
  N[6,1] <- max(round(rnorm(1,Nmean[6], Nsd[6]),0),10)	    # adult males
  
  ## DEMOGRAPHIC PARAMETERS
  # fecundity;   
  f <- matrix(ncol=P, nrow=1) 
  # Juvenile summer survival
  PHI1 <- matrix(ncol=P, nrow=1)
  
  phi3_var=phi3_sd*phi3_sd
  dp3<-estBetaParams(phi3_m,phi3_var)  
  PHI3<-rbeta(1,dp3$alpha,dp3$beta)
  
  phi1_var=phi1_sd*phi1_sd
  f_var=f_sd*f_sd  
  
  for(i in 1:P){
    
    dp1<-estBetaParams(phi1_m,phi1_var)
    df<-estBetaParams(f_m,f_var)
    
    PHI1[i]<-rbeta(1,dp1$alpha,dp1$beta)
    f[i]<-rbeta(1,df$alpha,df$beta)
    
  }
  
  #############################
  # SYSTEM PROCESS
  #############################

  h=H1
  
  for (t in 1:(P-1)){ 
    ###########################################################
    # STATE PROCESS;
    # PRE-HARVEST POPULATION VECTORS IN T+1  
    hadf=hadf.max-hadf.m*(K-(N[1,t]+N[2,t]+N[3,t]+N[4,t]+N[5,t]+N[6,t]))/K  
    H5<-round(N[5,t]*hadf)   
    hlow5<-ifelse(N[5,t]-T_Adf>0,N[5,t]-T_Adf,0)
    H[5,t] <- ifelse(N[5,t]-H5>T_Adf,H5,hlow5)
    
    T_Adm<-round((N[5,t]-H[5,t])/mRatio)   
    
    H[6,t] <- ifelse(N[6,t]-T_Adm>1,N[6,t]-T_Adm,0)
    
    H[1,t]  <-  ifelse(round(H[5,t]*hc/2)<N[1,t],round(H[5,t]*hc/2),0)
    H[2,t]  <-  ifelse(round(H[5,t]*hc/2)<N[2,t],round(H[5,t]*hc/2),0)
    
    N[3,t+1] <- rbinom(1, N[1,t]-H[1,t], PHI3)
    N[4,t+1] <- rbinom(1, N[2,t]-H[2,t], PHI3)
    
    N[5,t+1] <- rbinom(1, round(N[3,t]*(1-h[3])+N[5,t]-H[5,t]), PHI3)
    N[6,t+1] <- rbinom(1, round(N[4,t]*(1-h[4])+N[6,t]-H[6,t]), PHI3)
    
    
    N[1,t+1] <- rbinom(1, N[5, t+1], PHI1[t]*f[t]/2)   
    N[2,t+1] <- rbinom(1, N[5, t+1], PHI1[t]*f[t]/2)   
    
    H[3,t]  <-  round(N[3,t]*h[3])
    H[4,t]  <-  round(N[4,t]*h[4])
    
  }
  
  hadf=hadf.max-hadf.m*(K-(N[1,t]+N[2,t]+N[3,t]+N[4,t]+N[5,t]+N[6,t]))/K   
  H5<-round(N[5,P]*hadf) 
  hlow5<-ifelse(N[5,P]-T_Adf>0,N[5,P]-T_Adf,0)
  H[5,P] <- ifelse(N[5,P]-H5>T_Adf,H5,hlow5)
  
  T_Adm<-round((N[5,P]-H[5,P])/mRatio)   
  
  H[6,P] <- ifelse(N[6,P]-T_Adm>1,N[6,P]-T_Adm,0)
  
  H[1,P]  <-  ifelse(round(H[5,P]*hc/2)<N[1,P],round(H[5,P]*hc/2),0)
  H[2,P]  <-  ifelse(round(H[5,P]*hc/2)<N[2,P],round(H[5,P]*hc/2),0)
  H[3,P]  <-  round(N[3,P]*h[3])
  H[4,P]  <-  round(N[4,P]*h[4])
  
  for (t in 1:P){ 
    
    #############################################################
    # POST-HARVEST POPULATION VECTORS IN T+1
    X[1:4,t] <- (N[1:4,t]-H[1:4,t])
    X[5:6,t] <- (N[5:6,t]-H[5:6,t])
    
    
    #############################################################
    # DERIVED HARVEST NUMBERS  
    #H[,t] <- round(N[,t]*h[,t])
    
    X_tot[t] <- sum(X[,t])  # POST-HARVEST POPULATION size
    N_tot[t] <- sum(N[,t])  # summing up population vector to population size	
    H_tot[t] <- sum(H[,t])
    HU_tot[t] <- sum(H[3:6,t])
    HAd_tot[t] <- sum(H[5:6,t])
    
  }
  
  out <- list(N, X, H, N_tot, X_tot, H_tot,HU_tot, HAd_tot,
              f, PHI1, PHI3)
  names(out) <- c("N", "X","H", "N_tot","X_tot","H_tot","HU_tot",
                  "HAd_tot",
                  "f", "phi1", "phi3")
  out
  
}






