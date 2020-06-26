###################
#Nordfjella zone 2#
###################
#source("Functions_PrDetectingDiseasePopSizeVar.R")

#RR! =c(1,2,6) 

#pintro1=0.05     #5%
#pintro=0.001     #0.1%
#prior=0.5

#nsim1=1000

#Design prevalence:
#2016-2017: 2 individuals
#2018-2019: 3 individuals
#2020 -> 4 individuals

####################
#2016
####################
N16_yad2=c(Nmat2["N16","N1f"]+Nmat2["N16","N1m"],Nmat2["N16","Nadf"],Nmat2["N16","Nadm"])
Nmat_sdW16<-(Nmatnf2_sd["N16","N1f"]*Nmat2["N16","N1f"]+Nmatnf2_sd["N16","N1m"]*Nmat2["N16","N1m"])/(Nmat2["N16","N1f"]+Nmat2["N16","N1m"])
N16sd_yad2=c(Nmat_sdW16,Nmatnf2_sd["N16","Nadf"],Nmatnf2_sd["N16","Nadm"]) 

N=N16_yad2
N_sd=N16sd_yad2

nn=S16    #c(6,15,19) 
nn
Nt=sum(N)
PrRLN=pRLN2016   #0.92 #  0.845
PrRLN
pstar2 =2
pstar3 =3
pstar4 =4
pstar5  =5
pstar6  =6
pstar10 =10


PrLQ.min #0.02
PrLQ.max #0.6
PrLQ.mode #0.22

#x1A=192
PopSe16_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe16_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)

#summary(PopSe16_Rpstar2)

#######################
#2017
#######################

N17_yad2=c(Nmat2["N17","N1f"]+Nmat2["N17","N1m"],Nmat2["N17","Nadf"],Nmat2["N17","Nadm"])
Nmat_sdC17<-(Nmatnf2_sd["N17","N1f"]*Nmat2["N17","N1f"]+Nmatnf2_sd["N17","N1m"]*Nmat2["N17","N1m"])/(Nmat2["N17","N1f"]+Nmat2["N17","N1m"])
N17sd_yad2=c(Nmat_sdC17,Nmatnf2_sd["N17","Nadf"],Nmatnf2_sd["N17","Nadm"]) 

N=N17_yad2
N_sd=N17sd_yad2

nn=S17   #c(5,19,20)
nn

Nt=sum(N)
PrRLN=pRLN2017  #0.93
PrRLN
     
PopSe17_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe17_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)

#summary(PopSe17_Rpstar2)

########################

#2018
N18_yad2=c(Nmat2["N18","N1f"]+Nmat2["N18","N1m"],Nmat2["N18","Nadf"],Nmat2["N18","Nadm"])
Nmat_sdC18<-(Nmatnf2_sd["N18","N1f"]*Nmat2["N18","N1f"]+Nmatnf2_sd["N18","N1m"]*Nmat2["N18","N1m"])/(Nmat2["N18","N1f"]+Nmat2["N18","N1m"])
N18sd_yad2=c(Nmat_sdC18,Nmatnf2_sd["N18","Nadf"],Nmatnf2_sd["N18","Nadm"]) 

N=N18_yad2
N_sd=N18sd_yad2
Nt=sum(N)
Nt
S18

PrRLN=pRLN2018   #0.975

nn=S18     
PopSe18_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)

#summary(PopSe18_Rpstar3)

nn=S18w19 
nn
PopSe18w_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18w_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18w_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18w_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18w_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe18w_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)

#summary(PopSe18w_Rpstar3)


#2019
N19_yad2=c(Nmat2["N19","N1f"]+Nmat2["N19","N1m"],Nmat2["N19","Nadf"],Nmat2["N19","Nadm"])
Nmat_sdC19<-(Nmatnf2_sd["N19","N1f"]*Nmat2["N19","N1f"]+Nmatnf2_sd["N19","N1m"]*Nmat2["N19","N1m"])/(Nmat2["N19","N1f"]+Nmat2["N19","N1m"])
N19sd_yad2=c(Nmat_sdC19,Nmatnf2_sd["N19","Nadf"],Nmatnf2_sd["N19","Nadm"]) 

N
N_sd
S19
#S19F<-c(0,40,70) 

Nt=sum(N)
PrRLN=pRLN2019   #0.975
PrRLN            #0.83

nn=S19     
PopSe19_Rpstar2<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar2,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar3<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar3,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar4<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar4,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar5<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar6<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar6,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)
PopSe19_Rpstar10<-simSSe_variant(N=N,N_sd=N_sd,nn=nn,RR=RR1,pstarN=pstar10,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim1,x1Y=x1Y,x1A=x1A)

#summary(PopSe19_Rpstar3)

####################
#Combine years
####################
#Used for harvest scenario simulations
PrRLN=pRLN2018

pintro1 #[1] 0.05
pintro  #0.001

pH16<-pfree.est(SSe=PopSe16_Rpstar2,p.intro=pintro1,priorPFree=0.5)
apply(pH16,2,mean)
quantile(pH16$PFree,prob=c(0.025,0.5,0.975))

pH17<-pfree.est(SSe=PopSe17_Rpstar2,p.intro=pintro1,priorPFree=pH16$PFree)
apply(pH17,2,mean)
quantile(pH17$PFree,prob=c(0.025,0.5,0.975))

#Designprevalens:
#3 dyr
pH18<-pfree.est(SSe=PopSe18_Rpstar3,p.intro=pintro,priorPFree=pH17$PFree)
apply(pH18,2,mean)  #65 %
quantile(pH18$PFree,prob=c(0.025,0.5,0.975))

pH18w<-pfree.est(SSe=PopSe18w_Rpstar3,p.intro=pintro,priorPFree=pH17$PFree)
apply(pH18w,2,mean)  #76
quantile(pH18w$PFree,prob=c(0.025,0.5,0.975))


pH19<-pfree.est(SSe=PopSe19_Rpstar3,p.intro=pintro,priorPFree=pH18w$PFree)
apply(pH19,2,mean)  #87%
quantile(pH19$PFree,prob=c(0.025,0.5,0.975))


########################################################################################
######################################################################
########################################################################################
resultPFree<-rbind(c(mean(pH16$PFree),quantile(pH16$PFree,prob=c(0.025,0.5,0.975))),
                    c(mean(pH17$PFree),quantile(pH17$PFree,prob=c(0.025,0.5,0.975))),
                    c(mean(pH18$PFree),quantile(pH18$PFree,prob=c(0.025,0.5,0.975))),
                   c(mean(pH18w$PFree),quantile(pH18w$PFree,prob=c(0.025,0.5,0.975))),
                    c(mean(pH19$PFree),quantile(pH19$PFree,prob=c(0.025,0.5,0.975))))
resultPFree

############################

plot(c(2015,2016,2017,2018,2018.4,2019),c(50,100*resultPFree[,1]),type="b",pch=19,ylim=c(50,100),
     ylab="Probability of freedom from CWD (%)",xlab="Year",lwd=2,col="blue",xlim=c(2015,2025),las=1)
abline(h=99,lty=2,col="skyblue3",lwd=1)
abline(h=90,lty=2,col="darkblue",lwd=1)

segments(x0=c(2016,2017,2018,2018.4,2019), y0=resultPFree[,2]*100,y1=resultPFree[,4]*100,col="blue",lwd=2)
title("Nordfjella zone 2")

###########################



