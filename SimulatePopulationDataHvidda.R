######################################################
#Read data
######################################################

  #Estimated population size

  Nmat=read.csv2("Nmat_HV.csv")
  rownames(Nmat)=as.character(Nmat[,1])
  Nmat<-Nmat[,-1]
  
  Nmat<-as.data.frame(round(Nmat))
  Nmat$Total<-round(apply(Nmat,1,sum))
  names(Nmat)<-c("N0f","N0m","N1f","N1m","Nadf","Nadm","Total")

  N.sd=read.csv2("Nmatsd_HV.csv")
  rownames(N.sd)=as.character(N.sd[,1])
  N.sd<-N.sd[,-1]
  N.sd<-as.data.frame(round(N.sd))
  names(N.sd)<-c("N0f","N0m","N1f","N1m","Nadf","Nadm")
  Nmat_sd<-N.sd

  Nmat
  Nmat_sd

  N.sd15=N.sd["N15",]
  N.sd16=N.sd["N16",]
  N.sd17=N.sd["N17",]
  N.sd18=N.sd["N18",]
  
  N.mean15=Nmat["N15",][1:6]
  N.mean16=Nmat["N16",][1:6]
  N.mean17=Nmat["N17",][1:6]
  N.mean18=Nmat["N18",][1:6]
  
  #Demographic rates (estimated from Bayesian population model)
   demoratesHV<-read.csv2("demorates_HV.csv")
   rownames(demoratesHV)=as.character(demoratesHV[,1])
   demoratesHV<-demoratesHV[,-1]
   demoratesHV
   
   #Number harvested
   NhuntedHV<-read.csv2("nHarvested_HV.csv")
   rownames(NhuntedHV)=as.character(NhuntedHV[,1])
   NhuntedHV<-NhuntedHV[,-1]
   NhuntedHV
   
   H15=NhuntedHV["H15",][1:6]
   H16=NhuntedHV["H16",][1:6]
   H17=NhuntedHV["H17",][1:6]
   
  #####
   demorates<-read.csv2("demorates_HV.csv")
   f_m <- demorates[1,"mean"]
   phi1_m<- demorates[2,"mean"]
   phi3_m<- demorates[3,"mean"]
   
   f_sd<- demorates[1,"sd"]
   phi1_sd<- demorates[2,"sd"]
   phi3_sd<- demorates[3,"sd"]
   
   f.m   =  f_m      #0.649 
   phi1.m  = phi1_m  #0.942 
   phi3.m  =  phi3_m #0.934 
   
   f.sd   =  f_sd
   phi1.sd  = phi1_sd 
   phi3.sd = phi3_sd
   
   #Number of samples tested, 
   #In this data: NA age class is distributed according to distribution of known age class
   Stab<-read.csv2(file=paste("Stab_HV.csv",sep=""))
   S16 <-as.numeric(Stab[Stab[,1]=="S16",c(2:4)])
   S17 <-as.numeric(Stab[Stab[,1]=="S17",c(2:4)])
   S18 <-as.numeric(Stab[Stab[,1]=="S18",c(2:4)])
   S19 <-as.numeric(Stab[Stab[,1]=="S19",c(2:4)])
   
   #Proportion of samples including RLN in addition to obex 
   pRLN2016=as.numeric(Stab[Stab[,1]=="S16","pRLN"])
   pRLN2017=as.numeric(Stab[Stab[,1]=="S17","pRLN"])
   pRLN2018=as.numeric(Stab[Stab[,1]=="S18","pRLN"]) 
   pRLN2019=as.numeric(Stab[Stab[,1]=="S19","pRLN"]) 
   
   #harvest rate
   hratesHV<- read.csv2("hrates_HV.csv")
   rownames(hratesHV)=as.character(hratesHV[,1])
   hratesHV<-hratesHV[,-1]
   hratesHV
   
   #Define "ordinary" hunting rate as mean of 3 years
   h13<-round(apply(hratesHV,2,mean),3)
   h13
   
   #Strategy where 0 calves and 0 yearlings are harvested
   h0<-c(0,0,0,0,0,0.10,0.10)
   
   ######################################################
   #Set parameter values
   ######################################################
   
pstar4=4

RR1=c(1,2,6)
  
pintro1=0.01
pintro=0.001     #0.1%
prior=0.5

nsim1=1000

######################################################
#Source files needed for simulations 
######################################################

source("Functions_PrDetectingDisease.R")
#defines the functions sim3L etc used for simulating the diagnostic test sensitivity

source("SimFunctionsPaper.R")
source("SimPopModelsPaper.R")


#############################################################################
#Estimate probability of freedom obtained after each year of harvesting
#and testing for CWD (2016-2019)
#############################################################################

simSSe_variant <- simSSe3L
source("Pfreedom_calculations_Hvidda2019.R")

HPopSe16=PopSe16_Rpstar4
HPopSe17=PopSe17_Rpstar4
HPopSe18=PopSe18_Rpstar4
summary(HPopSe16) #0.11
summary(HPopSe17) #0.31
summary(HPopSe18) #0.26

quantile(pH18$PFree,prob=c(0.025,0.5,0.975))    #66.6 #67.9%
quantile(pH19$PFree,prob=c(0.025,0.5,0.975))    #82.0

pfree16<-round(mean(pH16$PFree),3)
pfree16Low<-round(quantile(pH16$PFree,prob=0.025),3)
pfree16High<-round(quantile(pH16$PFree,prob=0.975),3)

pfree17<-round(mean(pH17$PFree),3)
pfree17Low<-round(quantile(pH17$PFree,prob=0.025),3)
pfree17High<-round(quantile(pH17$PFree,prob=0.975),3)


pfree18hv<-round(mean(pH18$PFree),3)
pfree18hvLow<-round(quantile(pH18$PFree,prob=0.025),3)
pfree18hvHigh<-round(quantile(pH18$PFree,prob=0.975),3)

pfree19<-round(mean(pH19$PFree),3)
pfree19Low<-round(quantile(pH19$PFree,prob=0.025),3)
pfree19High<-round(quantile(pH19$PFree,prob=0.975),3)

pfree19
pfree19Low
pfree19High

testedHV<-c(0.5,pfree16,pfree17,pfree18hv,pfree19)
testedHVLow<-c(0.5,pfree16Low,pfree17Low,pfree18hvLow,pfree19Low)
testedHVHigh<-c(0.5,pfree16High,pfree17High,pfree18hvHigh,pfree19High)

rbind(testedHV,testedHVLow,testedHVHigh)

###################################
##################################
#Run simulations 
##################################

#Specify function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioHV.R 
calcPFreedom10<-calcPFreedom10meanDSE_i


P=11   #number of years

f_m      #0.649
phi1_m   #0.942
phi3_m   #0.934

f_sd
phi1_sd
phi3_sd

f.m   =  f_m      
phi1.m  = phi1_m  
phi3.m  =  phi3_m 

f.sd   =  f_sd
phi1.sd  = phi1_sd 
phi3.sd = phi3_sd


N.mean18<-as.numeric(Nmat["N18",][1:6])
N.sd18<-as.numeric(Nmat_sd["N18",][1:6])

N.mean=N.mean18
N.sd=N.sd18

N.mean
N.sd

pstar4s=pstar4 #pstar for year 2018 with data

#############################
#############################

PrRLN=pRLN2018

#Harvest rate of adult females (hadf) adjusted in relation to carrying capacity (Ktot)
Ktot=10000 
hadf_max=0.14+0.033
hadf_m=hadf_max
hadf=hadf_max-hadf_m*(Ktot-sum(N.mean18))/Ktot 
hadf #0.143

####
#Population simulation model -  functions from SimPopModels.R:
#SimPop18RKtot_h
#SimPop18Ktot_h

#TEST
mr5test.2500<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h13,mfRatio=5,Tadf=2500)
#dim(mr5test.2500$Htab)

##############
#Scenario 1
#############
#hadf at level of ordinary harvest rate h13["Hadf"]
hadf
h13["Hadf"]

h11=h13

mmKh13_3000<-MakeSimPopTabRKh18(Simff=SimPop18Ktot_h,Nsim=1000,H1tab=h11,mfRatio=5,Tadf=3000)
objects(mmKh13_3000)


#Scenario No 1 #
modS<-mmKh13_3000
mod1=modS

source("PlotScenarioHV.R")
#dtab
source("MakeSummaryVar.R")
#tabV
S2hv_tabV<-tabV

probfreeS
S2hv_probfreeS<-probfreeS
S2hv_probfreeSLow<-probfreeSLow
S2hv_probfreeSHigh<-probfreeSHigh

#Strategy with double harvest rate
#Scenario No 2 #
Ktot=10000 
hadf_max=0.173*2
hadf_m=hadf_max
h11=h13
hadf=hadf_max-hadf_m*(Ktot-sum(N.mean18))/Ktot 
hadf 
h11*2

mmKh13D_3000<-MakeSimPopTabRKh18(Simff=SimPop18Ktot_h,Nsim=1000,H1tab=h11*2,mfRatio=5,Tadf=3000)

modS<-mmKh13D_3000
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
S2Dhv_tabV<-tabV
probfreeS
S2Dhv_probfreeS<-probfreeS
S2Dhv_probfreeSLow<-probfreeSLow
S2Dhv_probfreeSHigh<-probfreeSHigh


#Scenario No 3 #

mmKh13D_2800<-MakeSimPopTabRKh18(Simff=SimPop18Ktot_h,Nsim=1000,H1tab=h11*2,mfRatio=5,Tadf=2800)

modS<-mmKh13D_2800
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
S2Dtadf2800hv_tabV<-tabV
probfreeS
S2Dtadf2800hv_probfreeS<-probfreeS
S2Dtadf2800hv_probfreeSLow<-probfreeSLow
S2Dtadf2800hv_probfreeSHigh<-probfreeSHigh

#Scenario No 4 #
mmKh13D_2500<-MakeSimPopTabRKh18(Simff=SimPop18Ktot_h,Nsim=1000,H1tab=h11*2,mfRatio=5,Tadf=2500)

modS<-mmKh13D_2500
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
S2Dtadf2500hv_tabV<-tabV
probfreeS
S2Dtadf2500hv_probfreeS<-probfreeS
S2Dtadf2500hv_probfreeSLow<-probfreeSLow
S2Dtadf2500hv_probfreeSHigh<-probfreeSHigh



#hadf back to basic level:
Ktot=10000 
hadf_max=0.14 + 0.033
hadf_m=hadf_max
hadf=hadf_max-hadf_m*(Ktot-sum(N.mean18))/Ktot 
hadf 
h13[5]


###########################################################################################
##########################################################################################
#mratio 1:3, 1:5, 1:10, 1:20
#hadf
#########################################################################################
###########################################################################################

mr5test.3000<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h13,mfRatio=5,Tadf=3000)

hcalf=0
mr5test.3000_c00<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h0,mfRatio=5,Tadf=3000)

#hcalf=0.4
#mr5test.3000_c04<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h13,mfRatio=5,Tadf=3000)

#Scenario No. 7 #
modS<-mr5test.3000
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
mr5hv_tabV<-tabV
probfreeS
mr5hv_probfreeS<-probfreeS
mr5hv_probfreeSLow<-probfreeSLow
mr5hv_probfreeSHigh<-probfreeSHigh


#Scenario No. 8 #
modS<-mr5test.3000_c00
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
mr5c00hv_tabV<-tabV
probfreeS
mr5c00hv_probfreeS<-probfreeS
mr5c00hv_probfreeSLow<-probfreeSLow
mr5c00hv_probfreeSHigh<-probfreeSHigh


mr3test.3000<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h13,mfRatio=3,Tadf=3000)
hcalf=0
mr3test.3000_c00<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h0,mfRatio=3,Tadf=3000)


#Scenario 5#
modS<-mr3test.3000
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
mr3hv_tabV<-tabV
probfreeS
mr3hv_probfreeS<-probfreeS
mr3hv_probfreeSLow<-probfreeSLow
mr3hv_probfreeSHigh<-probfreeSHigh

#Scenario 6#
modS<-mr3test.3000_c00
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
mr3c00hv_tabV<-tabV
probfreeS
mr3c00hv_probfreeS<-probfreeS
mr3c00hv_probfreeSLow<-probfreeSLow
mr3c00hv_probfreeSHigh<-probfreeSHigh


mr10test.3000<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h13,mfRatio=10,Tadf=3000)
hcalf=0
mr10test.3000_c00<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h0,mfRatio=10,Tadf=3000)

#Scenario 9#
modS<-mr10test.3000
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
mr10hv_tabV<-tabV
probfreeS
mr10hv_probfreeS<-probfreeS
mr10hv_probfreeSLow<-probfreeSLow
mr10hv_probfreeSHigh<-probfreeSHigh

#Scenario 10#
modS<-mr10test.3000_c00
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
mr10c00hv_tabV<-tabV
probfreeS
mr10c00hv_probfreeS<-probfreeS
mr10c00hv_probfreeSLow<-probfreeSLow
mr10c00hv_probfreeSHigh<-probfreeSHigh


mr20test.3000<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h13,mfRatio=20,Tadf=3000)
hcalf=0
mr20test.3000_c00<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h0,mfRatio=20,Tadf=3000)

#Scenario 11#
modS<-mr20test.3000
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
mr20hv_tabV<-tabV
probfreeS
mr20hv_probfreeS<-probfreeS
mr20hv_probfreeSLow<-probfreeSLow
mr20hv_probfreeSHigh<-probfreeSHigh

#Scenario 12#
modS<-mr20test.3000_c00
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
mr20c00hv_tabV<-tabV
probfreeS
mr20c00hv_probfreeS<-probfreeS
mr20c00hv_probfreeSLow<-probfreeSLow
mr20c00hv_probfreeSHigh<-probfreeSHigh

######################
#Effect RR   #########
######################
modS1= mmKh13_3000 
modS = mr5test.3000_c00

#################
RR1=c(1,1,1)
#################
source("Pfreedom_calculations_Hvidda2019.R")

#Scenario 13#
mod1=modS1
source("PlotScenarioHV.R")
probfreeS
Hv_RR111_probfreeS<-probfreeS
Hv_RR111_probfreeSLow<-probfreeSLow
Hv_RR111_probfreeSHigh<-probfreeSHigh
source("MakeSummaryVar.R")
Hv_RR111_tabV<-tabV


#Scenario 18#
mod1=modS
source("PlotScenarioHV.R")
dtab
probfreeS
mr5Hv_RR111_probfreeS<-probfreeS
mr5Hv_RR111_probfreeSLow<-probfreeSLow
mr5Hv_RR111_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
tabV
mr5Hv_RR111_tabV<-tabV

#################
RR1=c(1,2,2)
#################
source("Pfreedom_calculations_Hvidda2019.R")
#PopSe16=PopSe16_Rpstar4
#PopSe17=PopSe17_Rpstar4
#PopSe18=PopSe18_Rpstar4
#mean(PopSe16)
#mean(PopSe17)
#mean(PopSe18)

#Scenario 14#
mod1=modS1
source("PlotScenarioHV.R")
probfreeS
Hv_RR122_probfreeS<-probfreeS
Hv_RR122_probfreeSLow<-probfreeSLow
Hv_RR122_probfreeSHigh<-probfreeSHigh
source("MakeSummaryVar.R")
tabV
Hv_RR122_tabV<-tabV

#Scenario 19#
mod1=modS
source("PlotScenarioHV.R")
probfreeS
mr5Hv_RR122_probfreeS<-probfreeS
mr5Hv_RR122_probfreeSLow<-probfreeSLow
mr5Hv_RR122_probfreeSHigh<-probfreeSHigh
source("MakeSummaryVar.R")
tabV
mr5Hv_RR122_tabV<-tabV



#################
RR1=c(1,2,4)
#################
source("Pfreedom_calculations_Hvidda2019.R")

#Scenario 15#
modS1<-mmKh13_3000
mod1=modS1
source("PlotScenarioHV.R")
probfreeS
Hv_RR124_probfreeS<-probfreeS
Hv_RR124_probfreeSLow<-probfreeSLow
Hv_RR124_probfreeSHigh<-probfreeSHigh
source("MakeSummaryVar.R")
Hv_RR124_tabV<-tabV

#Scenario 20#
modS<-mr5test.3000_c00
mod1=modS
source("PlotScenarioHV.R")
probfreeS
mr5Hv_RR124_probfreeS<-probfreeS
mr5Hv_RR124_probfreeSLow<-probfreeSLow
mr5Hv_RR124_probfreeSHigh<-probfreeSHigh
source("MakeSummaryVar.R")
mr5Hv_RR124_tabV<-tabV
mod1=modS1

#################
RR1=c(1,2,6)
#################
simSSe_variant <- simSSe3L
source("Pfreedom_calculations_Hvidda2019.R")

#Scenario 1#
modS1<-mmKh13_3000
mod1=modS1
source("PlotScenarioHV.R")
probfreeS
Hv_RR126_probfreeS<-probfreeS
Hv_RR126_probfreeSLow<-probfreeSLow
Hv_RR126_probfreeSHigh<-probfreeSHigh
source("MakeSummaryVar.R")
Hv_RR126_tabV<-tabV

#Scenario 8#
modS<-mr5test.3000_c00
mod1=modS
source("PlotScenarioHV.R")
probfreeS
mr5Hv_RR126_probfreeS<-probfreeS
mr5Hv_RR126_probfreeSLow<-probfreeSLow
mr5Hv_RR126_probfreeSHigh<-probfreeSHigh
source("MakeSummaryVar.R")
mr5Hv_RR126_tabV<-tabV
mod1=modS1


#######################################
#Test stochastic RR    ################
#######################################
rradf<-rpert(nsim1,x.min=1,x.max=2.5,x.mode=2)
hist(rradf)
rradm<-rpert(nsim1,x.min=2,x.max=6.5,x.mode=6)
hist(rradm)
rrY<-rep(1,nsim1)
RRmat=cbind(rrY,rradf,rradm)

RR1=RRmat
simSSe_variant <- simSSe3L_RRi
source("Pfreedom_calculations_Hvidda2019.R")

resultPFreeHV_Stoch1<-resultPFreeHV
resultPFreeHV_Stoch1 #Prob free for emp. data assuming stochastic RR version 1

PrRLN=pRLN2018
PrRLN
N.mean=as.numeric(N.mean18[1:6])
N.sd=as.numeric(N.sd18)
h13
calcPFreedom10<-calcPFreedom10_RRi 

#Scenario 16 (harvest strategy No 1 with stochastic RR)
modS<-mmKh13_3000
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")

probfreeS
S2hv_probfreeSRRi<-probfreeS
S2hv_probfreeSRRiLow<-probfreeSLow
S2hv_probfreeSRRiHigh<-probfreeSHigh

#Scenario 21 (harvest strategy No 8 with stochastic RR)
mod0<-mr5test.3000_c00
mod1=mod0
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")

probfreeS
mr5c00hv_probfreeSRRi<-probfreeS
mr5c00hv_probfreeSRRiLow<-probfreeSLow
mr5c00hv_probfreeSRRiHigh<-probfreeSHigh




################################
#Stochastic RR Version 2
################################
rradf<-rpert(nsim1,x.min=1,x.max=2.5,x.mode=1.5)
hist(rradf)
rradm<-rpert(nsim1,x.min=2,x.max=6.5,x.mode=4.5)
hist(rradm)
rrY<-rep(1,nsim1)
RRmat2=cbind(rrY,rradf,rradm)

RR1=RRmat2
simSSe_variant <- simSSe3L_RRi
source("Pfreedom_calculations_Hvidda2019.R")
resultPFreeHV_Stoch2<-resultPFreeHV
resultPFreeHV_Stoch2 ##Prob free for emp. data assuming stochastic RR version 1

PrRLN=pRLN2018
N.mean=as.numeric(N.mean18[1:6])
N.sd=as.numeric(N.sd18)

calcPFreedom10<-calcPFreedom10_RRi 

#Scenario 17 (harvest strategy No 1 with stochastic RR)
modS<-mmKh13_3000
mod1=modS
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")
probfreeS
S2hv_probfreeSRRi2<-probfreeS
S2hv_probfreeSRRi2Low<-probfreeSLow
S2hv_probfreeSRRi2High<-probfreeSHigh

S2hv_probfreeSRRi2[5,]
S2hv_probfreeS[5,]

#Scenario 22 (harvest strategy No 8 with stochastic RR)
mod0<-mr5test.3000_c00
mod1=mod0
source("PlotScenarioHV.R")
source("MakeSummaryVar.R")

probfreeS
mr5c00hv_probfreeSRRi2<-probfreeS
mr5c00hv_probfreeSRRi2Low<-probfreeSLow
mr5c00hv_probfreeSRRi2High<-probfreeSHigh

mr5c00hv_probfreeSRRi2[5,]
mr5c00hv_probfreeS[5,]

#############################
#############################
RR1=c(1,2,6) #Back to basic level
simSSe_variant <- simSSe3L
source("Pfreedom_calculations_Hvidda2019.R")
calcPFreedom10<-calcPFreedom10meanDSE_i
#############################

prior=0.5
nsim1=1000
pintro1=0.05
pintro=0.05 #Increase from basic 0.001
nsim1=1000


modS1= mmKh13_3000 
modS = mr5test.3000_c00

#Scenario No 23
mod1=modS1
source("PlotScenarioHV.R")
#dtab
#probfreeS
Hv_probfreeS05<-probfreeS
Hv_probfreeS05Low<-probfreeSLow
Hv_probfreeS05High<-probfreeSHigh


#Scenario No 28
mod1=modS
source("PlotScenarioHV.R")
#dtab
#probfreeS
mr5Hv_probfreeS05<-probfreeS
mr5Hv_probfreeS05Low<-probfreeSLow
mr5Hv_probfreeS05High<-probfreeSHigh


########################################################################
#Update with other combinations of pintro1 and pintro
prior=0.5
nsim1=1000
pintro1=0.05 #Increase from basic 0.01
pintro=0.01  #Increase from basic 0.001
nsim1=1000

###########################################################################

modS1= mmKh13_3000 
modS = mr5test.3000_c00

#Scenario No 24
mod1=modS1
source("PlotScenarioHV.R")
#dtab
#probfreeS
Hv_probfreeS0501<-probfreeS
Hv_probfreeS0501Low<-probfreeSLow
Hv_probfreeS0501High<-probfreeSHigh


#Scenario No 29
mod1=modS
source("PlotScenarioHV.R")
#dtab
#probfreeS
mr5Hv_probfreeS0501<-probfreeS
mr5Hv_probfreeS0501Low<-probfreeSLow
mr5Hv_probfreeS0501High<-probfreeSHigh

########################################################################
#Update with other combinations of pintro1 and pintro
prior=0.5
nsim1=1000
pintro1=0.05 #Increase from basic 0.01
pintro=0.005  #Increase from basic 0.001
nsim1=1000

#########################################################################

modS1= mmKh13_3000 
modS = mr5test.3000_c00

#Scenario No 25
mod1=modS1
source("PlotScenarioHV.R")
#dtab
#probfreeS
Hv_probfreeS05005<-probfreeS
Hv_probfreeS05005Low<-probfreeSLow
Hv_probfreeS05005High<-probfreeSHigh


#Scenario No 30
mod1=modS
source("PlotScenarioHV.R")
#dtab
#probfreeS
mr5Hv_probfreeS05005<-probfreeS
mr5Hv_probfreeS05005Low<-probfreeSLow
mr5Hv_probfreeS05005High<-probfreeSHigh

########################################################################
#Update with other combinations of pintro1 and pintro
prior=0.5
nsim1=1000
pintro1=0.05 #Increase from basic 0.01
pintro=0.001 #basic 0.001
nsim1=1000


########################################################################################

modS1= mmKh13_3000 
modS = mr5test.3000_c00

#Scenario No 26
mod1=modS1
source("PlotScenarioHV.R")
#dtab
#probfreeS
Hv_probfreeS05001<-probfreeS
Hv_probfreeS05001Low<-probfreeSLow
Hv_probfreeS05001High<-probfreeSHigh


#Scenario No 31
mod1=modS
source("PlotScenarioHV.R")
#dtab
#probfreeS
mr5Hv_probfreeS05001<-probfreeS
mr5Hv_probfreeS05001Low<-probfreeSLow
mr5Hv_probfreeS05001High<-probfreeSHigh

########################################################################
#Update with other combinations of pintro1 and pintro
prior=0.5
nsim1=1000
pintro1=0.01
pintro=0.01 #Increase from basic 0.001
nsim1=1000

#########################################################################

modS1= mmKh13_3000 
modS = mr5test.3000_c00

#Scenario No 27
mod1=modS1
source("PlotScenarioHV.R")
#dtab
#probfreeS
Hv_probfreeS01<-probfreeS
Hv_probfreeS01Low<-probfreeSLow
Hv_probfreeS01High<-probfreeSHigh


#Scenario No 32
mod1=modS
source("PlotScenarioHV.R")
#dtab
#probfreeS
mr5Hv_probfreeS01<-probfreeS
mr5Hv_probfreeS01Low<-probfreeSLow
mr5Hv_probfreeS01High<-probfreeSHigh




################################
#Back to standard level
################################
pintro1=0.01
pintro=0.001     #0.1%
#simSSe_variant <- simSSe3L
#source("Pfreedom_calculations_Hvidda2019.R")
#Function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioHV.R 
#calcPFreedom10<-calcPFreedom10meanDSE_i

#####################################################
#source(SummaryResTabHV.R)
#####################################################

