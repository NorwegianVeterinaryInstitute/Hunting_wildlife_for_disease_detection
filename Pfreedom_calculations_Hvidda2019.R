################
#Hardangervidda#
################

#RR1 #c(1,2,6), vector with RR for yearlings, adult females and adult males relative to yearlings, or if stochasstic, 
    #a matrix with one row for each simulation
pintro1 #0.01     1%
pintro #0.001     0.1%
prior #0.5
nsim1 #1000

S16
S17
S18
S19

pRLN2016
pRLN2017
pRLN2018
pRLN2019

############################
############################
############################

#Design prevalence
pstar1 =1
pstar2 =2
pstar3 =3
pstar4 =4
pstar5  =5
pstar6  =6
pstar10 =10

#Design prevalence:
#2016 -> 4 individuals

rownames(Nmat) <-c("N15","N16","N17","N18","N19")
rownames(Nmat_sd) <-c("N15","N16","N17","N18","N19")


####################
#2016
####################
N16_yad2=c(Nmat["N16","N1f"]+Nmat["N16","N1m"],Nmat["N16","Nadf"],Nmat["N16","Nadm"])
Nmat_sdW16<-(Nmat_sd["N16","N1f"]*Nmat["N16","N1f"]+Nmat_sd["N16","N1m"]*Nmat["N16","N1m"])/(Nmat["N16","N1f"]+Nmat["N16","N1m"])
N16sd_yad2=c(Nmat_sdW16,Nmat_sd["N16","Nadf"],Nmat_sd["N16","Nadm"]) 

N=N16_yad2
N_sd=N16sd_yad2

PrRLN=pRLN2016   

nn=S16    
PopSe16_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
summary(PopSe16_Rpstar4)

PopSe16_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)


#######################
#2017
#######################

N17_yad2=c(Nmat["N17","N1f"]+Nmat["N17","N1m"],Nmat["N17","Nadf"],Nmat["N17","Nadm"])
Nmat_sdC17<-(Nmat_sd["N17","N1f"]*Nmat["N17","N1f"]+Nmat_sd["N17","N1m"]*Nmat["N17","N1m"])/(Nmat["N17","N1f"]+Nmat["N17","N1m"])
N17sd_yad2=c(Nmat_sdC17,Nmat_sd["N17","Nadf"],Nmat_sd["N17","Nadm"]) 

N=N17_yad2
N_sd=N17sd_yad2

PrRLN=pRLN2017  #0.095

nn=S17   
PopSe17_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)

#######################
#2018
#######################

N18_yad2=c(Nmat["N18","N1f"]+Nmat["N18","N1m"],Nmat["N18","Nadf"],Nmat["N18","Nadm"])
Nmat_sdC18<-(Nmat_sd["N18","N1f"]*Nmat["N18","N1f"]+Nmat_sd["N18","N1m"]*Nmat["N18","N1m"])/(Nmat["N18","N1f"]+Nmat["N18","N1m"])
N18sd_yad2=c(Nmat_sdC18,Nmat_sd["N18","Nadf"],Nmat_sd["N18","Nadm"]) 

N=N18_yad2
N_sd=N18sd_yad2

PrRLN=pRLN2018   #0.74

nn=S18   
PopSe18_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)

#######################
#2019
#######################

N19_yad2=c((Nmat["N19","N1f"]+Nmat["N19","N1m"]),Nmat["N19","Nadf"],Nmat["N19","Nadm"])
Nmat_sdC19<-(Nmat_sd["N19","N1f"]*Nmat["N19","N1f"]+Nmat_sd["N19","N1m"]*Nmat["N19","N1m"])/(Nmat["N19","N1f"]+Nmat["N19","N1m"])
N19sd_yad2=c(Nmat_sdC19,Nmat_sd["N19","Nadf"],Nmat_sd["N19","Nadm"]) 

N=N19_yad2
N_sd=N19sd_yad2

PrRLN=pRLN2019   #0.76

nn=S19
PopSe19_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
summary(PopSe19_Rpstar4)
PopSe19_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)


####################
#Combine years
####################

#Used for harvest scenario simulations
PrRLN=pRLN2018

pintro1   #0.01
pintro    #0.001

pH16<-pfree.est(SSe=PopSe16_Rpstar4,p.intro=pintro1,priorPFree=0.5)
pH17<-pfree.est(SSe=PopSe17_Rpstar4,p.intro=pintro1,priorPFree=pH16$PFree)
pH18<-pfree.est(SSe=PopSe18_Rpstar4,p.intro=pintro,priorPFree=pH17$PFree)
pH19<-pfree.est(SSe=PopSe19_Rpstar4,p.intro=pintro,priorPFree=pH18$PFree)

########################################################################################
resultPFreeHV<-rbind(c(mean(pH16$PFree),quantile(pH16$PFree,prob=c(0.025,0.5,0.975))),
                   c(mean(pH17$PFree),quantile(pH17$PFree,prob=c(0.025,0.5,0.975))),
                   c(mean(pH18$PFree),quantile(pH18$PFree,prob=c(0.025,0.5,0.975))),
                   c(mean(pH19$PFree),quantile(pH19$PFree,prob=c(0.025,0.5,0.975))))
resultPFreeHV
########################################################################################


