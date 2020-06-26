#######################################
#Run population simulation model
#######################################
MakeSimPopTabRKh18<-function(Simff=SimPopK_h,Nsim=1000,H1tab=h1,mfRatio=10,Tadf=2500)
{
  simdatH<-SimPop18Ktot_h(T_Adf=Tadf,mRatio=mfRatio,H1=H1tab)
  
  Ntab<-array(NA,c(dim(simdatH$N),Nsim))
  Htab<-array(NA,c(dim(simdatH$H),Nsim))
  Xtab<-array(NA,c(dim(simdatH$X),Nsim))

  Ntot_tab<-Xtot_tab<-Htot_tab<-f_tab<-phi1_tab<-phi3_tab<-array(NA,c(Nsim,11))
  
  for (i in 1:Nsim){
    
    simdati<-Simff(T_Adf=Tadf,mRatio=mfRatio,H1=H1tab)
    Ntab[,,i]<-simdati$N
    Htab[,,i]<-simdati$H
    Xtab[,,i]<-simdati$X
    Ntot_tab[i,]<-simdati$N_tot
    Xtot_tab[i,]<-simdati$X_tot
    Htot_tab[i,]<-simdati$H_tot
    f_tab[i,]<-simdati$f
    phi1_tab[i,]<-simdati$phi1
    phi3_tab[i,]<-simdati$phi3
  }
  
  Nmean<-apply(Ntab[,,],c(1,2),mean)
  Nlow<-round(apply(Ntab[,,],c(1,2),quantile,prob=0.025))
  Nhigh<-round(apply(Ntab[,,],c(1,2),quantile,prob=0.975))
  
  Xmean<-apply(Xtab[,,],c(1,2),mean)
  Xmean1 <-round(Xmean)
  Xlow<-round(apply(Xtab[,,],c(1,2),quantile,prob=0.025))
  Xhigh<-round(apply(Xtab[,,],c(1,2),quantile,prob=0.975))
  
  Htab1<-round(apply(Htab[,,],c(1,2),mean))
  Htab1low<-round(apply(Htab[,,],c(1,2),quantile,prob=0.025))
  Htab1high<-round(apply(Htab[,,],c(1,2),quantile,prob=0.975))
  
  Nmean1 <-round(Nmean)
  
  H_SD<-round(apply(Htab[,,],c(1,2),sd))
  N_SD<-apply(Ntab[,,],c(1,2),sd)
  
  out <- list(Ntab,Htab,Xtab,Nmean1,Nlow,Nhigh,Xmean1,Xlow,Xhigh, Htab1,Htab1low,Htab1high,N_SD,H_SD)
  names(out) <- c("Ntab","Htab","Xtab","Nmean1","Nlow","Nhigh","Xmean1","Xlow","Xhigh", "Htab1","Htab1low",
                  "Htab1high","N_SD","H_SD")
  
  out
  
}


#############################################################
#############################################################
#Posterior probability is calculated by updating the prior probability by using Bayes' theorem. 
#Assuming Sp=1
#Assuming prior = Pr(population (country) is free)
#NPV = Pr(D-|T-) = pfree
#NPV = (1-prior)/((1-prior)+prior(1-SSe)) = (1-prior)/(1-prior*Se)

#SSe = population/system sensitivity (one column for each of P years)
#p.intro = vector with risk of introduction for each of P years

pfree.estPYr<-function(SSe,p.intro,priorPFree){
  
  PFreeTab<-array(NA,dim(SSe))
  PostPInf_tl1 = (1-priorPFree)
  priorPInf = PostPInf_tl1 + p.intro[1] - (PostPInf_tl1*p.intro[1])
  PFree = (1-priorPInf)/(1-priorPInf*SSe[,1])
  PFreeTab[,1]<-PFree
  P=ncol(SSe)
  for (i in 2:P)  {
    #p.intro1=ifelse(length(p.intro==1),p.intro,p.intro[i])
    p.intro1 = p.intro[i]
    PostPInf_tl1<- (1-PFreeTab[,(i-1)])
    priorPInf = PostPInf_tl1 + p.intro1 - (PostPInf_tl1*p.intro1)
    PFree = (1-priorPInf)/(1-priorPInf*SSe[,i])
    PFreeTab[,i]<-PFree
  } 
  PFreeTab 
}




#New harvest scenario from 2018
#mod1 = output from population simulation
#pstarV = vector of design prevalence for each of k years
calcPFreedom10meanDSE_i<-function(mod=mod1, PopHist=cbind(PSe16,PSe17),simSSe_i=simSSe3LmeanDSE_i,pstarI=pstarV)
{                                                  
  PopSemod1<-array(0,c(nsim1,11))
  for (k in 1:11){
    pstarS<-pstarI[k]
    #Use iterations to account or variation in nn-number harvested/tested
    PopSeV<-simSSe_i(N=mod$Ntab,nn=mod$Htab,k1=k,RR=RR1,pstarN=pstarS,PrRLN=PrRLN,dSE=dSeTab1,
                   nsim=nsim1,x1Y=x1Y,x1A=x1A)
    PopSemod1[,k]=PopSeV
    
  }
  
  SSe1=cbind(PopHist,PopSemod1)
  pfreemod1<- pfree.estPYr(SSe=SSe1,p.intro=c(pintro1,pintro1,rep(pintro,11)),priorPFree=0.5)
  colnames(pfreemod1)<-c("y2016","y2017","y2018","y2019",
                         "y2020","y2021","y2022","y2023","y2024","y2025","y2026","y2027","y2028")
  pfreemod1
  
}


#New harvest scenario from 2018
#mod1 = output from population simulation
#pstarV = vector of design prevalence for each of k years
calcPFreedom10_i<-function(mod=mod1, PopHist=cbind(PSe16,PSe17),simSSe_i=simSSe3L_i,pstarI=pstarV)
{                                                  
  PopSemod1<-array(0,c(nsim1,11))
  for (k in 1:11){
    pstarS<-pstarI[k]
    #Use iterations (1:nsim1) to account or variation in nn-number harvested/tested
    PopSeV<-simSSe_i(N=mod$Ntab,nn=mod$Htab,k1=k,RR=RR1,pstarN=pstarS,PrRLN=PrRLN,dSE=dSeTab1,
                     nsim=nsim1,x1Y=x1Y,x1A=x1A)
    PopSemod1[,k]=PopSeV
    
  }
  
  SSe1=cbind(PopHist,PopSemod1)
  pfreemod1<- pfree.estPYr(SSe=SSe1,p.intro=c(pintro1,pintro1,rep(pintro,11)),priorPFree=0.5)
  colnames(pfreemod1)<-c("y2016","y2017","y2018","y2019",
                         "y2020","y2021","y2022","y2023","y2024","y2025","y2026","y2027","y2028")
  pfreemod1
  
}


calcPFreedom10_RRi<-function(mod=mod1, PopHist=cbind(PSe16,PSe17),simSSe_i=simSSe3LmeanDSE_RRi,RR1matrix=RRmat,pstarI=pstarV)
{                                                  
  PopSemod1<-array(0,c(nsim1,11))
  for (k in 1:11){
    pstarS<-pstarI[k]
    #Use iterations (1:nsim1) to account or variation in nn-number harvested/tested
    PopSeV<-simSSe_i(N=mod$Ntab,nn=mod$Htab,k1=k,RR=RR1matrix,pstarN=pstarS,PrRLN=PrRLN,dSE=dSeTab1,
                     nsim=nsim1,x1Y=x1Y,x1A=x1A)
    PopSemod1[,k]=PopSeV
    
  }
  
  SSe1=cbind(PopHist,PopSemod1)
  pfreemod1<- pfree.estPYr(SSe=SSe1,p.intro=c(pintro1,pintro1,rep(pintro,11)),priorPFree=0.5)
  colnames(pfreemod1)<-c("y2016","y2017","y2018","y2019",
                         "y2020","y2021","y2022","y2023","y2024","y2025","y2026","y2027","y2028")
  pfreemod1
  
}



assembleDatA<-function(mod1){
  Nmean1<-mod1$Nmean1 [,c(1:8)]
  colnames(Nmean1)<-c(2018,2019,2020,2021,2022,2023,2024,2025)
  rownames(Nmean1)<-c("N0f","N0fm","N1f","N1m","Nadf","Nadm")
  Nmean1<-as.data.frame(Nmean1)
  Nmean1$Var<-"Nmean1"
  
  Nsd1<-round(mod1$N_SD,1)[,c(1:8)]
  colnames(Nsd1)<-c(2018,2019,2020,2021,2022,2023,2024,2025)
  rownames(Nsd1)<-c("N0f","N0fm","N1f","N1m","Nadf","Nadm")
  Nsd1<-as.data.frame(Nsd1)
  Nsd1$Var<-"N_SD"
  
  Htab1<-mod1$Htab1 [,c(1:8)]
  colnames(Htab1)<-c(2018,2019,2020,2021,2022,2023,2024,2025)
  rownames(Htab1)<-c("N0f","N0fm","N1f","N1m","Nadf","Nadm")
  Htab1<-as.data.frame(Htab1)
  Htab1$Var<-"Htab1"
  
  Sum_exclCalves<-apply(mod1$Nmean1[-c(1,2),c(1:8)],2,sum)
  Var<-"Sum_exclCalves"
  Sum_exclCalves<-c(Sum_exclCalves,Var)
  
  Nlow1<-mod1$Nlow[,c(1:8)]
  colnames(Nlow1)<-c(2018,2019,2020,2021,2022,2023,2024,2025)
  rownames(Nlow1)<-c("N0f","N0fm","N1f","N1m","Nadf","Nadm")
  Nlow1<-as.data.frame(Nlow1)
  Nlow1$Var<-"Nlow1"
  
  Nhigh1<-mod1$Nhigh[,c(1:8)]
  colnames(Nhigh1)<-c(2018,2019,2020,2021,2022,2023,2024,2025)
  rownames(Nhigh1)<-c("N0f","N0fm","N1f","N1m","Nadf","Nadm")
  Nhigh1<-as.data.frame(Nhigh1)
  Nhigh1$Var<-"Nhigh1"
  
  Hlow1<-mod1$Htab1low[,c(1:8)]
  colnames(Hlow1)<-c(2018,2019,2020,2021,2022,2023,2024,2025)
  rownames(Hlow1)<-c("N0f","N0fm","N1f","N1m","Nadf","Nadm")
  Hlow1<-as.data.frame(Hlow1)
  Hlow1$Var<-"Hlow1"
  
  Hhigh1<-mod1$Htab1high[,c(1:8)]
  colnames(Hhigh1)<-c(2018,2019,2020,2021,2022,2023,2024,2025)
  rownames(Hhigh1)<-c("N0f","N0fm","N1f","N1m","Nadf","Nadm")
  Hhigh1<-as.data.frame(Hhigh1)
  Hhigh1$Var<-"Hhigh1"
  
  
  samleT<-rbind(Nmean1,Sum_exclCalves,Nsd1,Htab1,Nlow1,Nhigh1,Hlow1,Hhigh1)
  return(samleT)
}



MakeSummaryPstarScenarioVar<-function(mod1){
  
  pfreemodS2=calcPFreedom10(mod=mod1,PopHist=cbind(PopSe16_Rpstar2,PopSe17_Rpstar2),pstarI=rep(pstar2,11))

  pfreemod=calcPFreedom10(mod=mod1,PopHist=cbind(PopSe16_Rpstar4,PopSe17_Rpstar4),pstarI=rep(pstar4,11))

  pfreemodS6=calcPFreedom10(mod=mod1,PopHist=cbind(PopSe16_Rpstar6,PopSe17_Rpstar6),pstarI=rep(pstar6,11))

  pfreemodS10=calcPFreedom10(mod=mod1,PopHist=cbind(PopSe16_Rpstar10,PopSe17_Rpstar10),pstarI=rep(pstar10,11))

  pstarVA<-c(3,3,rep(4,9))
  pfreemodSA=calcPFreedom10(mod=mod1,PopHist=cbind(PopSe16_Rpstar2,PopSe17_Rpstar2),pstarI=pstarVA)

  Pfree.pstar_scenario<-rbind(
    apply(pfreemodS2,2,mean),
    apply(pfreemod,2,mean),
    apply(pfreemodS6,2,mean),
    apply(pfreemodS10,2,mean),
    apply(pfreemodSA,2,mean))
  rownames(Pfree.pstar_scenario)<-c("pstar2","pstar4","pstar6","pstar10","pstarSA")
  
  Pfree.pstar_scenarioLow<-rbind(
    apply(pfreemodS2,2,quantile,prob=0.025),
    apply(pfreemod,2,quantile,prob=0.025),
    apply(pfreemodS6,2,quantile,prob=0.025),
    apply(pfreemodS10,2,quantile,prob=0.025),
    apply(pfreemodSA,2,quantile,prob=0.025))
  rownames(Pfree.pstar_scenarioLow)<-c("pstar2","pstar4","pstar6","pstar10","pstarSA")
  
  Pfree.pstar_scenarioHigh<-rbind(
    apply(pfreemodS2,2,quantile,prob=0.975),
    apply(pfreemod,2,quantile,prob=0.975),
    apply(pfreemodS6,2,quantile,prob=0.975),
    apply(pfreemodS10,2,quantile,prob=0.975),
    apply(pfreemodSA,2,quantile,prob=0.975))
  rownames(Pfree.pstar_scenarioHigh)<-c("pstar2","pstar4","pstar6","pstar10","pstarSA")
  
  Pfree.pstarscenario <- list("Pfree.pstarS"=Pfree.pstar_scenario, "Pfree.pstarSLow"=Pfree.pstar_scenarioLow, 
                              "Pfree.pstarSHigh"=Pfree.pstar_scenarioHigh)
  return(Pfree.pstarscenario)
}


#Plot function used in PlotScenarioNf2
MakeProbPlotNf<-function(pfreeS){
  plot(seq(2015,2028,1),c(0.5,pfreeS["pstar4",])*100,type="b",lwd=2,ylab="Probability of freedom (%)",xlab="Year",ylim=c(50,100))
  abline(h=99,lty=2)
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstar2",])*100,type="b",lwd=2,col="skyblue")
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstar6",])*100,type="b",col="blue",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstar10",])*100,type="b",col="purple",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstarSA",])*100,type="b",col="green",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstar4",])*100,type="b",col="black",lwd=2)
  
  legend("bottomright",
         c("pstar=2","pstar=4","pstar=6","pstar=10","pstar= 2 - 4"),bty="n",lwd=2,lty=1,
         col=c("skyblue","black","blue","purple","green"),pch=19,cex=1)
}


#Plot function used in PlotScenarioHV
MakeProbPlot<-function(pfreeS){
  plot(seq(2015,2028,1),c(0.5,pfreeS["pstar4",])*100,type="b",lwd=2,ylab="Probability of freedom (%)",xlab="Year",ylim=c(50,100))
  abline(h=99,lty=2)
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstar2",])*100,type="b",lwd=2,col="skyblue")
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstar6",])*100,type="b",col="blue",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstar10",])*100,type="b",col="purple",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeS["pstar4",])*100,type="b",col="black",lwd=2)
  
  legend("bottomright",
         c("pstar=2","pstar=4","pstar=6","pstar=10"),bty="n",lwd=2,lty=1,
         col=c("skyblue","black","blue","purple"),pch=19,cex=1)
}


