
#Posterior probability is calculated by updating the prior probability by using Bayes' theorem. 
#Assuming Sp=1
#Assuming prior = Pr(population (country) is free)
#NPV = Pr(D-|T-) = pfree
#NPV = (1-prior)/((1-prior)+prior(1-SSe)) = (1-prior)/(1-prior*Se)
#SSe = population/system sensitivity

#Discounting prior probability by risk of introduction
pfree.est<-function(SSe,p.intro,priorPFree){
  PostPInf_tl1 = (1-priorPFree)
  priorPInf = PostPInf_tl1 + p.intro - (PostPInf_tl1*p.intro)
  PFree = (1-priorPInf)/(1-priorPInf*SSe)
  as.data.frame(cbind(SSe,priorPFree,p.intro,priorPInf,PFree))
}


#Example
#pp<-pfree.est(SSe=0.45,p.intro=0.01,priorPFree=0.5)


#############################################################################
#Application 2 Estimating probability of detecting disease
#############################################################################

source("Functions_modelling_dSe1.R")

#Calculate dSeTab1 - table with expected Se (estimated by rescaling infection curves) 
#for testing RLN and high and low-quality obex samples

#IncludeSeVar - function to add stocasticity around dSe given by infection curve
#x1A - number of max time steps since infection for adults
#x1Y - number of max time steps since infection for yearlings
#8 time steps = 1 month


##############################################################################
#simSSe3L - function to estimate Surveilance Sensitivity
#Assuming 3 risk levels, separating between yearlings, adult females and adult males
################################################################################
#Time since infection are randomly drawn from time line

#N = c(mean population size yearlings, population size adult females and adult males)
#N_sd = c(sd population size yearlings, population size adult females and adult males)
#n = c(number of yearlings sampled, number of adult females and males sampled)
#RR = c(1,2,4) - relative risk of yearlings, adult females and adult males 
#pstarN = design prevalencegiven as number infected
#PrRLN =  Proportion of samples inlcuding RLN
#dSeTab = table with expected Se (estimated by rescaling infection curves) for testing RLN and high and low-quality obex samples  
#nsim = nr of simulations
#x1Y = vector with possible infection times for infected yearlings
#x1A = vector with possible infection times for infected adults

#rpert-function are utilized to random draw values from a distribution defined by min, max and mode
#PrLQ.min=0.02,
#PrLQ.max=0.60,
#PrLQ.mode=0.22,

simSSe3L<-function(N=N,N_sd=N_sd, nn=nn,RR=RR1,pstarN=pstarN,dSE=dSe1,PrRLN=PrRLN,x1Y=x1Y,x1A=x1A,nsim=nsim){

    #n1 = nn[1] +nn[2] +nn[3]
        
    #PrAf<-nn[2]/n1
    #PrAm<-nn[3]/n1 
     
    SSe<-rep(NA,nsim)
    Pr_AllTestingNeg<-SSe
    

    #Draw random proportion of LQ obex sample from given pert-distribution
    PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
    
    N1<-rnorm(n=nsim,mean=N[1],sd=N_sd[1])
    N2<-rnorm(n=nsim,mean=N[2],sd=N_sd[2])
    N3<-rnorm(n=nsim,mean=N[3],sd=N_sd[3])

    N1[N1<=0]<-1
    N2[N2<=0]<-1
    N3[N3<=0]<-1

    for (i in c(1:nsim)){
    
    Ntot=N1[i] + N2[i] + N3[i]
    
    pstar=round(pstarN/Ntot,4)    
    
    #Adjusted risk         
    ARy <- 1/ (RR[3]*(N3[i]/Ntot)+RR[2]*(N2[i]/Ntot) + (N1[i]/Ntot))  #yearlings
    ARadf <- RR[2]*ARy                            #adult females
    ARadm <- RR[3]*ARy                            #adult males

    nnv1<-ifelse(nn[1]<N1[i],nn[1],N1[i]-1)
    nnv2<-ifelse(nn[2]<N2[i],nn[2],N2[i]-1) 
    nnv3<-ifelse(nn[3]<N3[i],nn[3],N3[i]-1) 
 
    #Draw random time since infection    
    tpY<-sample(x1Y,nnv1,replace=T)          #yearlings
    tpAdf<-sample(x1A,nnv2,replace=T)        #adult females
    tpAdm<-sample(x1A,nnv3,replace=T)        #adult males

    #dSe for yearlings       
    dSE_RLNobexY<-dSE[tpY,"yRLN"]
    dSE_obexYLQ<-dSE[tpY,"yObexLowq"]
    dSE_obexYHQ<-dSE[tpY,"yObexLowq"]
    
    #dSe for adult females       
    dSE_RLNobexAdf<-dSE[tpAdf,"yRLN"]
    dSE_obexAdfLQ<-dSE[tpAdf,"yObexLowq"]
    dSE_obexAdfHQ<-dSE[tpAdf,"yObexHighq"]
    
    #dSe for adult males       
    dSE_RLNobexAdm<-dSE[tpAdm,"yRLN"]
    dSE_obexAdmLQ<-dSE[tpAdm,"yObexLowq"]
    dSE_obexAdmHQ<-dSE[tpAdm,"yObexHighq"]
 
    #dSe for yearlings with stochasticy      
    dSE_RLNobexYi<-IncludeSeVar(Y=dSE_RLNobexY)
    dSE_obexYLQi<-IncludeSeVar(Y=dSE_obexYLQ)
    dSE_obexYHQi<-IncludeSeVar(Y=dSE_obexYHQ)
    
    #dSe for adult females with stochasticity       
    dSE_RLNobexAdfi<-IncludeSeVar(Y=dSE_RLNobexAdf)
    dSE_obexAdfLQi<-IncludeSeVar(Y=dSE_obexAdfLQ)
    dSE_obexAdfHQi<-IncludeSeVar(Y=dSE_obexAdfHQ)
    
    #dSe for adult males with stochasticity      
    dSE_RLNobexAdmi<-IncludeSeVar(Y=dSE_RLNobexAdm)
    dSE_obexAdmLQi<-IncludeSeVar(Y=dSE_obexAdmLQ)
    dSE_obexAdmHQi<-IncludeSeVar(Y=dSE_obexAdmHQ)
    
    #Probability of positive test result given Age class, Sex and tissue type
    PrYRLN<-ARy*PrRLN*dSE_RLNobexYi
    PrAdfRLN<-ARadf*PrRLN*dSE_RLNobexAdfi
    PrAdmRLN<-ARadm*PrRLN*dSE_RLNobexAdmi

    PrYobex<-ARy*(1-PrRLN)*(PrLQv[i]*dSE_obexYLQi + (1-PrLQv[i])* dSE_obexYHQi)
    PrAdfobex<-ARadf*(1-PrRLN)*(PrLQv[i]*dSE_obexAdfLQi + (1-PrLQv[i])* dSE_obexAdfHQi)
    PrAdmobex<-ARadm*(1-PrRLN)*(PrLQv[i]*dSE_obexAdmLQi + (1-PrLQv[i])* dSE_obexAdmHQi)
   
    SeY<-PrYRLN+PrYobex
    SeAdf<-PrAdfRLN+PrAdfobex
    SeAdm<-PrAdmRLN+PrAdmobex

    #Pr_NonTestingPos
    Pr_AllTestingNeg[i]<- (1 - sum(SeY,SeAdf,SeAdm)/Ntot)^(pstar*Ntot)
    
    #SSe = Surveillance sensitivity
    SSe[i] <- 1 - Pr_AllTestingNeg[i]

    }
                
    SSe
    
    }


#############################
#############################
#N and nn is given as a matrix of simulated data (nrow=nsim)
simSSe3LmeanDSE_i<-function(N=N,nn=nn,k1=k,RR=RR1,pstarN=pstarN,dSE=dSe1,PrRLN=PrRLN,x1Y=x1Y,x1A=x1A,nsim=nsim){
  
  SSe<-rep(NA,nsim)
  Pr_AllTestingNeg<-SSe
  
  #Draw random proportion of LQ obex sample from given pert-distribution
  PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
  
  for (i in c(1:nsim)){
    N1=N[3,k1,i]+N[4,k1,i]
    N2=N[5,k1,i]
    N3=N[6,k1,i]
    Ntot=N1 + N2 + N3
    
    pstar=round(pstarN/Ntot,4)    
    
    #Adjusted risk         
    ARy <- 1/ (RR[3]*(N3/Ntot)+RR[2]*(N2/Ntot) + (N1/Ntot))  #yearlings
    ARadf <- RR[2]*ARy                            #adult females
    ARadm <- RR[3]*ARy                            #adult males
    
    nn1=nn[3,k1,i]+nn[4,k1,i]
    nn2=nn[5,k1,i]
    nn3=nn[6,k1,i]
    
    nnv1<-ifelse(nn1<N1,nn1,N1-1)
    nnv2<-ifelse(nn2<N2,nn2,N2-1) 
    nnv3<-ifelse(nn3<N3,nn3,N3-1) 
    
    #Draw random time since infection    
    tpY<-sample(x1Y,nnv1,replace=T)          #yearlings
    tpAdf<-sample(x1A,nnv2,replace=T)        #adult females
    tpAdm<-sample(x1A,nnv3,replace=T)        #adult males
    
    #dSe for yearlings       
    dSE_RLNobexY<-dSE[tpY,"yRLN"]
    dSE_obexYLQ<-dSE[tpY,"yObexLowq"]
    dSE_obexYHQ<-dSE[tpY,"yObexLowq"]
    
    #dSe for adult females       
    dSE_RLNobexAdf<-dSE[tpAdf,"yRLN"]
    dSE_obexAdfLQ<-dSE[tpAdf,"yObexLowq"]
    dSE_obexAdfHQ<-dSE[tpAdf,"yObexHighq"]
    
    #dSe for adult males       
    dSE_RLNobexAdm<-dSE[tpAdm,"yRLN"]
    dSE_obexAdmLQ<-dSE[tpAdm,"yObexLowq"]
    dSE_obexAdmHQ<-dSE[tpAdm,"yObexHighq"]
    
    #dSe for yearlings without stochasticy      
    dSE_RLNobexYi<-dSE_RLNobexY
    dSE_obexYLQi<-dSE_obexYLQ
    dSE_obexYHQi<- dSE_obexYHQ
    
    #dSe for adult females without stochasticity       
    dSE_RLNobexAdfi<-dSE_RLNobexAdf
    dSE_obexAdfLQi<- dSE_obexAdfLQ
    dSE_obexAdfHQi<-dSE_obexAdfHQ
    
    #dSe for adult males without stochasticity      
    dSE_RLNobexAdmi<-dSE_RLNobexAdm
    dSE_obexAdmLQi<-dSE_obexAdmLQ
    dSE_obexAdmHQi<-dSE_obexAdmHQ
    
    #Probability of positive test result given Age class, Sex and tissue type
    PrYRLN<-ARy*PrRLN*dSE_RLNobexYi
    PrAdfRLN<-ARadf*PrRLN*dSE_RLNobexAdfi
    PrAdmRLN<-ARadm*PrRLN*dSE_RLNobexAdmi
    
    PrYobex<-ARy*(1-PrRLN)*(PrLQv[i]*dSE_obexYLQi + (1-PrLQv[i])* dSE_obexYHQi)
    PrAdfobex<-ARadf*(1-PrRLN)*(PrLQv[i]*dSE_obexAdfLQi + (1-PrLQv[i])* dSE_obexAdfHQi)
    PrAdmobex<-ARadm*(1-PrRLN)*(PrLQv[i]*dSE_obexAdmLQi + (1-PrLQv[i])* dSE_obexAdmHQi)
    
    SeY<-PrYRLN+PrYobex
    SeAdf<-PrAdfRLN+PrAdfobex
    SeAdm<-PrAdmRLN+PrAdmobex
    
    #Pr_NonTestingPos
    Pr_AllTestingNeg[i]<- (1 - sum(SeY,SeAdf,SeAdm)/Ntot)^(pstar*Ntot)
    
    #SSe = Surveillance sensitivity
    SSe[i] <- 1 - Pr_AllTestingNeg[i]
    
  }
  
  SSe
  
}


#N and nn is given as a matrix of simulated data (nrow=nsim)
simSSe3L_i<-function(N=N,nn=nn,k1=k,RR=RR1,pstarN=pstarN,dSE=dSe1,PrRLN=PrRLN,x1Y=x1Y,x1A=x1A,nsim=nsim){
  
  SSe<-rep(NA,nsim)
  Pr_AllTestingNeg<-SSe
  
  #Draw random proportion of LQ obex sample from given pert-distribution
  PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
  
  for (i in c(1:nsim)){
    N1=N[3,k1,i]+N[4,k1,i]
    N2=N[5,k1,i]
    N3=N[6,k1,i]
    Ntot=N1 + N2 + N3
    
    pstar=round(pstarN/Ntot,4)    
    
    #Adjusted risk         
    ARy <- 1/ (RR[3]*(N3/Ntot)+RR[2]*(N2/Ntot) + (N1/Ntot))  #yearlings
    ARadf <- RR[2]*ARy                            #adult females
    ARadm <- RR[3]*ARy                            #adult males
    
    nn1=nn[3,k1,i]+nn[4,k1,i]
    nn2=nn[5,k1,i]
    nn3=nn[6,k1,i]
    
    nnv1<-ifelse(nn1<N1,nn1,N1-1)
    nnv2<-ifelse(nn2<N2,nn2,N2-1) 
    nnv3<-ifelse(nn3<N3,nn3,N3-1) 
    
    #Draw random time since infection    
    tpY<-sample(x1Y,nnv1,replace=T)          #yearlings
    tpAdf<-sample(x1A,nnv2,replace=T)        #adult females
    tpAdm<-sample(x1A,nnv3,replace=T)        #adult males
    
    #dSe for yearlings       
    dSE_RLNobexY<-dSE[tpY,"yRLN"]
    dSE_obexYLQ<-dSE[tpY,"yObexLowq"]
    dSE_obexYHQ<-dSE[tpY,"yObexLowq"]
    
    #dSe for adult females       
    dSE_RLNobexAdf<-dSE[tpAdf,"yRLN"]
    dSE_obexAdfLQ<-dSE[tpAdf,"yObexLowq"]
    dSE_obexAdfHQ<-dSE[tpAdf,"yObexHighq"]
    
    #dSe for adult males       
    dSE_RLNobexAdm<-dSE[tpAdm,"yRLN"]
    dSE_obexAdmLQ<-dSE[tpAdm,"yObexLowq"]
    dSE_obexAdmHQ<-dSE[tpAdm,"yObexHighq"]
    
    #dSe for yearlings with stochasticy      
    dSE_RLNobexYi<-IncludeSeVar(Y=dSE_RLNobexY)
    dSE_obexYLQi<-IncludeSeVar(Y=dSE_obexYLQ)
    dSE_obexYHQi<-IncludeSeVar(Y=dSE_obexYHQ)
    
    #dSe for adult females with stochasticity       
    dSE_RLNobexAdfi<-IncludeSeVar(Y=dSE_RLNobexAdf)
    dSE_obexAdfLQi<-IncludeSeVar(Y=dSE_obexAdfLQ)
    dSE_obexAdfHQi<-IncludeSeVar(Y=dSE_obexAdfHQ)
    
    #dSe for adult males with stochasticity      
    dSE_RLNobexAdmi<-IncludeSeVar(Y=dSE_RLNobexAdm)
    dSE_obexAdmLQi<-IncludeSeVar(Y=dSE_obexAdmLQ)
    dSE_obexAdmHQi<-IncludeSeVar(Y=dSE_obexAdmHQ)
    
    #Probability of positive test result given Age class, Sex and tissue type
    PrYRLN<-ARy*PrRLN*dSE_RLNobexYi
    PrAdfRLN<-ARadf*PrRLN*dSE_RLNobexAdfi
    PrAdmRLN<-ARadm*PrRLN*dSE_RLNobexAdmi
    
    PrYobex<-ARy*(1-PrRLN)*(PrLQv[i]*dSE_obexYLQi + (1-PrLQv[i])* dSE_obexYHQi)
    PrAdfobex<-ARadf*(1-PrRLN)*(PrLQv[i]*dSE_obexAdfLQi + (1-PrLQv[i])* dSE_obexAdfHQi)
    PrAdmobex<-ARadm*(1-PrRLN)*(PrLQv[i]*dSE_obexAdmLQi + (1-PrLQv[i])* dSE_obexAdmHQi)
    
    SeY<-PrYRLN+PrYobex
    SeAdf<-PrAdfRLN+PrAdfobex
    SeAdm<-PrAdmRLN+PrAdmobex
    
    #Pr_NonTestingPos
    Pr_AllTestingNeg[i]<- (1 - sum(SeY,SeAdf,SeAdm)/Ntot)^(pstar*Ntot)
    
    #SSe = Surveillance sensitivity
    SSe[i] <- 1 - Pr_AllTestingNeg[i]
    
  }
  
  SSe
  
}


#Include stochasticity in RR
simSSe3L_RRi<-function(N=N,N_sd=N_sd, nn=nn,RR=RR1,pstarN=pstarN,dSE=dSe1,PrRLN=PrRLN,x1Y=x1Y,x1A=x1A,nsim=nsim){
  
  #n1 = nn[1] +nn[2] +nn[3]
  
  #PrAf<-nn[2]/n1
  #PrAm<-nn[3]/n1 
  
  SSe<-rep(NA,nsim)
  Pr_AllTestingNeg<-SSe
  
  
  #Draw random proportion of LQ obex sample from given pert-distribution
  PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
  
  N1<-rnorm(n=nsim,mean=N[1],sd=N_sd[1])
  N2<-rnorm(n=nsim,mean=N[2],sd=N_sd[2])
  N3<-rnorm(n=nsim,mean=N[3],sd=N_sd[3])
  
  N1[N1<=0]<-1
  N2[N2<=0]<-1
  N3[N3<=0]<-1
  
  for (i in c(1:nsim)){
    
    Ntot=N1[i] + N2[i] + N3[i]
    
    pstar=round(pstarN/Ntot,4)    
    
    #Adjusted risk         
    ARy <- 1/ (RR[i,3]*(N3[i]/Ntot)+RR[i,2]*(N2[i]/Ntot) + (N1[i]/Ntot))  #yearlings
    ARadf <- RR[i,2]*ARy                            #adult females
    ARadm <- RR[i,3]*ARy                            #adult males
    
    nnv1<-ifelse(nn[1]<N1[i],nn[1],N1[i]-1)
    nnv2<-ifelse(nn[2]<N2[i],nn[2],N2[i]-1) 
    nnv3<-ifelse(nn[3]<N3[i],nn[3],N3[i]-1) 
    
    #Draw random time since infection    
    tpY<-sample(x1Y,nnv1,replace=T)          #yearlings
    tpAdf<-sample(x1A,nnv2,replace=T)        #adult females
    tpAdm<-sample(x1A,nnv3,replace=T)        #adult males
    
    #dSe for yearlings       
    dSE_RLNobexY<-dSE[tpY,"yRLN"]
    dSE_obexYLQ<-dSE[tpY,"yObexLowq"]
    dSE_obexYHQ<-dSE[tpY,"yObexLowq"]
    
    #dSe for adult females       
    dSE_RLNobexAdf<-dSE[tpAdf,"yRLN"]
    dSE_obexAdfLQ<-dSE[tpAdf,"yObexLowq"]
    dSE_obexAdfHQ<-dSE[tpAdf,"yObexHighq"]
    
    #dSe for adult males       
    dSE_RLNobexAdm<-dSE[tpAdm,"yRLN"]
    dSE_obexAdmLQ<-dSE[tpAdm,"yObexLowq"]
    dSE_obexAdmHQ<-dSE[tpAdm,"yObexHighq"]
    
    #dSe for yearlings with stochasticy      
    dSE_RLNobexYi<-IncludeSeVar(Y=dSE_RLNobexY)
    dSE_obexYLQi<-IncludeSeVar(Y=dSE_obexYLQ)
    dSE_obexYHQi<-IncludeSeVar(Y=dSE_obexYHQ)
    
    #dSe for adult females with stochasticity       
    dSE_RLNobexAdfi<-IncludeSeVar(Y=dSE_RLNobexAdf)
    dSE_obexAdfLQi<-IncludeSeVar(Y=dSE_obexAdfLQ)
    dSE_obexAdfHQi<-IncludeSeVar(Y=dSE_obexAdfHQ)
    
    #dSe for adult males with stochasticity      
    dSE_RLNobexAdmi<-IncludeSeVar(Y=dSE_RLNobexAdm)
    dSE_obexAdmLQi<-IncludeSeVar(Y=dSE_obexAdmLQ)
    dSE_obexAdmHQi<-IncludeSeVar(Y=dSE_obexAdmHQ)
    
    #Probability of positive test result given Age class, Sex and tissue type
    PrYRLN<-ARy*PrRLN*dSE_RLNobexYi
    PrAdfRLN<-ARadf*PrRLN*dSE_RLNobexAdfi
    PrAdmRLN<-ARadm*PrRLN*dSE_RLNobexAdmi
    
    PrYobex<-ARy*(1-PrRLN)*(PrLQv[i]*dSE_obexYLQi + (1-PrLQv[i])* dSE_obexYHQi)
    PrAdfobex<-ARadf*(1-PrRLN)*(PrLQv[i]*dSE_obexAdfLQi + (1-PrLQv[i])* dSE_obexAdfHQi)
    PrAdmobex<-ARadm*(1-PrRLN)*(PrLQv[i]*dSE_obexAdmLQi + (1-PrLQv[i])* dSE_obexAdmHQi)
    
    SeY<-PrYRLN+PrYobex
    SeAdf<-PrAdfRLN+PrAdfobex
    SeAdm<-PrAdmRLN+PrAdmobex
    
    #Pr_NonTestingPos
    Pr_AllTestingNeg[i]<- (1 - sum(SeY,SeAdf,SeAdm)/Ntot)^(pstar*Ntot)
    
    #SSe = Surveillance sensitivity
    SSe[i] <- 1 - Pr_AllTestingNeg[i]
    
  }
  
  SSe
  
}

#Include stochasticity in RR
#N and nn is given as a matrix of simulated data (nrow=nsim)
simSSe3LmeanDSE_RRi<-function(N=N,nn=nn,k1=k,RR=RR1,pstarN=pstarN,dSE=dSe1,PrRLN=PrRLN,x1Y=x1Y,x1A=x1A,nsim=nsim){
#RR1 is a matrix with nsim rows and 3 columns with relativ risk of Yearlings, adult females and adult males versus yearlings. 
  SSe<-rep(NA,nsim)
  Pr_AllTestingNeg<-SSe
  
  #Draw random proportion of LQ obex sample from given pert-distribution
  PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
  
  for (i in c(1:nsim)){
    N1=N[3,k1,i]+N[4,k1,i]
    N2=N[5,k1,i]
    N3=N[6,k1,i]
    Ntot=N1 + N2 + N3
    
    pstar=round(pstarN/Ntot,4)    
    
    #Adjusted risk         
    ARy <- 1/ (RR[i,3]*(N3/Ntot)+RR[i,2]*(N2/Ntot) + (N1/Ntot))  #yearlings
    ARadf <- RR[i,2]*ARy                            #adult females
    ARadm <- RR[i,3]*ARy                            #adult males
    
    nn1=nn[3,k1,i]+nn[4,k1,i]
    nn2=nn[5,k1,i]
    nn3=nn[6,k1,i]
    
    nnv1<-ifelse(nn1<N1,nn1,N1-1)
    nnv2<-ifelse(nn2<N2,nn2,N2-1) 
    nnv3<-ifelse(nn3<N3,nn3,N3-1) 
    
    #Draw random time since infection    
    tpY<-sample(x1Y,nnv1,replace=T)          #yearlings
    tpAdf<-sample(x1A,nnv2,replace=T)        #adult females
    tpAdm<-sample(x1A,nnv3,replace=T)        #adult males
    
    #dSe for yearlings       
    dSE_RLNobexY<-dSE[tpY,"yRLN"]
    dSE_obexYLQ<-dSE[tpY,"yObexLowq"]
    dSE_obexYHQ<-dSE[tpY,"yObexLowq"]
    
    #dSe for adult females       
    dSE_RLNobexAdf<-dSE[tpAdf,"yRLN"]
    dSE_obexAdfLQ<-dSE[tpAdf,"yObexLowq"]
    dSE_obexAdfHQ<-dSE[tpAdf,"yObexHighq"]
    
    #dSe for adult males       
    dSE_RLNobexAdm<-dSE[tpAdm,"yRLN"]
    dSE_obexAdmLQ<-dSE[tpAdm,"yObexLowq"]
    dSE_obexAdmHQ<-dSE[tpAdm,"yObexHighq"]
    
    #dSe for yearlings without stochasticy      
    dSE_RLNobexYi<-dSE_RLNobexY
    dSE_obexYLQi<-dSE_obexYLQ
    dSE_obexYHQi<- dSE_obexYHQ
    
    #dSe for adult females without stochasticity       
    dSE_RLNobexAdfi<-dSE_RLNobexAdf
    dSE_obexAdfLQi<- dSE_obexAdfLQ
    dSE_obexAdfHQi<-dSE_obexAdfHQ
    
    #dSe for adult males without stochasticity      
    dSE_RLNobexAdmi<-dSE_RLNobexAdm
    dSE_obexAdmLQi<-dSE_obexAdmLQ
    dSE_obexAdmHQi<-dSE_obexAdmHQ
    
    #Probability of positive test result given Age class, Sex and tissue type
    PrYRLN<-ARy*PrRLN*dSE_RLNobexYi
    PrAdfRLN<-ARadf*PrRLN*dSE_RLNobexAdfi
    PrAdmRLN<-ARadm*PrRLN*dSE_RLNobexAdmi
    
    PrYobex<-ARy*(1-PrRLN)*(PrLQv[i]*dSE_obexYLQi + (1-PrLQv[i])* dSE_obexYHQi)
    PrAdfobex<-ARadf*(1-PrRLN)*(PrLQv[i]*dSE_obexAdfLQi + (1-PrLQv[i])* dSE_obexAdfHQi)
    PrAdmobex<-ARadm*(1-PrRLN)*(PrLQv[i]*dSE_obexAdmLQi + (1-PrLQv[i])* dSE_obexAdmHQi)
    
    SeY<-PrYRLN+PrYobex
    SeAdf<-PrAdfRLN+PrAdfobex
    SeAdm<-PrAdmRLN+PrAdmobex
    
    #Pr_NonTestingPos
    Pr_AllTestingNeg[i]<- (1 - sum(SeY,SeAdf,SeAdm)/Ntot)^(pstar*Ntot)
    
    #SSe = Surveillance sensitivity
    SSe[i] <- 1 - Pr_AllTestingNeg[i]
    
  }
  
  SSe
  
}

