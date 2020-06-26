
######################################################
#Read data
#From 2018: Data for Nf2 includes a small reindeer flock in a neighbouring area
######################################################
Nhunted<-read.csv2("nHarvested_Nf2.csv")
rownames(Nhunted)=as.character(Nhunted[,1])
Nhunted<-Nhunted[,-1]

#Estimated population size
N.mean=read.csv2("Nmat_Nf2.csv")
rownames(N.mean)=as.character(N.mean[,1])
N.mean<-N.mean[,-1]

N.sd=read.csv2("Nmatsd_Nf2.csv")
rownames(N.sd)=as.character(N.sd[,1])
N.sd<-N.sd[,-1]

Nmat2=N.mean
Nmatnf2_sd=N.sd

N.sd15=Nmatnf2_sd["N15",]
N.sd16=Nmatnf2_sd["N16",]
N.sd17=Nmatnf2_sd["N17",]
N.sd18=Nmatnf2_sd["N18",]

N.mean15=Nmat2["N15",][1:6]
N.mean16=Nmat2["N16",][1:6]
N.mean17=Nmat2["N17",][1:6]
N.mean18=Nmat2["N18",][1:6]

H15=Nhunted["H15",][1:6]
H16=Nhunted["H16",][1:6]
H17=Nhunted["H17",][1:6]

###
demorates<-read.csv2("demorates_Nf2.csv")
f_m <- demorates[1,"mean"]
phi1_m<- demorates[2,"mean"]
phi3_m<- demorates[3,"mean"]

f_sd<- demorates[1,"sd"]
phi1_sd<- demorates[2,"sd"]
phi3_sd<- demorates[3,"sd"]

#Number of samples tested
StabNf2<-read.csv2(file="Stab_Nf2.csv")
S16 <-as.numeric(StabNf2[StabNf2[,1]=="S16",c(2:4)])
S17 <-as.numeric(StabNf2[StabNf2[,1]=="S17",c(2:4)])
S18 <-as.numeric(StabNf2[StabNf2[,1]=="S18",c(2:4)])
#S18F <-c(2,26,107)   #suggested quoata before hunting 2018

#Including extra harvested winter 2019
S18w19<-S18 +as.numeric(StabNf2[StabNf2[,1]=="S1819W",c(2:4)])

S19 <-as.numeric(StabNf2[StabNf2[,1]=="S19",c(2:4)])

#Set number of NA ageclass to yearling
S19[1]<-S19[1]+as.numeric(StabNf2[StabNf2[,1]=="S19","ageNA"])
S19

#Proportion of samples including RLN in addition to obex 
pRLN2016=as.numeric(StabNf2[StabNf2[,1]=="S16","pRLN"])
pRLN2017=as.numeric(StabNf2[StabNf2[,1]=="S17","pRLN"])
pRLN2018=as.numeric(StabNf2[StabNf2[,1]=="S18","pRLN"]) 
pRLN2019=as.numeric(StabNf2[StabNf2[,1]=="S19","pRLN"]) 

#Define "ordinary" hunting rate as mean of 3 years
hrates<-read.csv2("hrates_nf2.csv")
h13<-round(apply(hrates[,-1],2,mean),3)
h13

#Strategy where 0 calves and 0 yearlings are harvested
h0<-c(0,0,0,0,0.10,0.10)

######################################################
#Set parameter values
######################################################
pstar4=4
RR1=c(1,2,6)

pintro1=0.05
pintro=0.001     #0.1%
prior=0.5

nsim1=1000

######################################################
#Source files needed for simulations 
######################################################

source("Functions_PrDetectingDisease.R")
#defines the functions sim3L etc used for simulating the diagnostic test sensitivity

source("SimPopModelsPaper.R")
source("SimFunctionsPaper.R")

#############################################################################
#Estimate probability of freedom obtained after each year of harvesting
#and testing for CWD (2016-2019)
#############################################################################
simSSe_variant <- simSSe3L
source("Pfreedom_calculations_Nf2.R")

quantile(pH18$PFree,prob=c(0.025,0.5,0.975))    
quantile(pH18w$PFree,prob=c(0.025,0.5,0.975))    
quantile(pH19$PFree,prob=c(0.025,0.5,0.975))   
mean(pH19$PFree)    

#Save outcome - used for making figures 
pfree16<-round(mean(pH16$PFree),3)
pfree16Low<-round(quantile(pH16$PFree,prob=0.025),3)
pfree16High<-round(quantile(pH16$PFree,prob=0.975),3)

pfree17<-round(mean(pH17$PFree),3)
pfree17Low<-round(quantile(pH17$PFree,prob=0.025),3)
pfree17High<-round(quantile(pH17$PFree,prob=0.975),3)

pfree18<-round(mean(pH18$PFree),3)
pfree18Low<-round(quantile(pH18$PFree,prob=0.025),3)
pfree18High<-round(quantile(pH18$PFree,prob=0.975),3)

pfree18W<-round(mean(pH18w$PFree),3)
pfree18WLow<-round(quantile(pH18w$PFree,prob=0.025),3)
pfree18WHigh<-round(quantile(pH18w$PFree,prob=0.975),3)

pfree19<-round(mean(pH19$PFree),3)
pfree19Low<-round(quantile(pH19$PFree,prob=0.025),3)
pfree19High<-round(quantile(pH19$PFree,prob=0.975),3)

pfree19
pfree19Low
pfree19High

testedNf<-c(0.5,pfree16,pfree17,pfree18,pfree18W,pfree19)
testedNfLow<-c(0.5,pfree16Low,pfree17Low,pfree18Low,pfree18WLow,pfree19Low)
testedNfHigh<-c(0.5,pfree16High,pfree17High,pfree18High,pfree18WHigh,pfree19High)

#The updated probability of freedom after surveillence 2019:
pfree19Nf2=pfree19
pfree19Nf2Low=pfree19Low
pfree19Nf2High=pfree19High

#The updated probability of freedom for each year of surveillance:
testedNf #mean
testedNfLow #2.5th percentile
testedNfHigh #97.5th percentile


###################################
##################################
#Run simulation scenario
##################################
###################################

#Specify function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioNf2.R 
calcPFreedom10<-calcPFreedom10meanDSE_i

P=11   #number of years

PrRLN=pRLN2018
PrRLN

N.mean=as.numeric(N.mean18[1:6])
N.sd=as.numeric(N.sd18)

f.m   =  f_m      #0.632 
phi1.m  = phi1_m  #0.948 
phi3.m  =  phi3_m #0.971

f.sd   =  f_sd
phi1.sd  = phi1_sd 
phi3.sd = phi3_sd

h13

#Harvest rate of adult females (hadf) adjusted in relation to carrying capacity (Ktot)
Ktot=800 
hadf_max=0.141
hadf_m=hadf_max
hadf=hadf_max-hadf_m*(Ktot-Nmat2["N18","Total"])/Ktot 
hadf #0.113

#Population simulation model (functions from SimPopModels.R):
#SimPop18Ktot_h - Harvest strategy 1 ‘ordinary’ (with threshold of maximum sex ratio)
#SimPop18RKtot_h - Harvest strategy 2 ‘proactive’ (with operational sex ratio)

#TEST
mr5test.150<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h13,mfRatio=5,Tadf=150)
#dim(mr5test.150$Htab)

##################################
#Scenario No. 1
#################################
#RR1=c(1,2,6)
RR1

#hadf at level of ordinary harvest rate h13["Hadf"]
hadf
h13["Hadf"]

mmKh13_150<-MakeSimPopTabRKh18(Simff=SimPop18Ktot_h,Nsim=1000,H1tab=h13,mfRatio=5,Tadf=150)

#Strategy with double harvest rate
hadf_max=0.141*2
hadf_m=hadf_max
hadf=hadf_max-hadf_m*(Ktot-Nmat2["N18","Total"])/Ktot 
hadf #0.225

mmKh13D_150<-MakeSimPopTabRKh18(Simff=SimPop18Ktot_h,Nsim=1000,H1tab=h13*2,mfRatio=5,Tadf=150)

#harvest rates back to basic level
hadf_max=0.141
hadf_m=hadf_max
hadf=hadf_max-hadf_m*(Ktot-Nmat2["N18","Total"])/Ktot 
hadf #0.112

#Scenario No 1 #
modS<-mmKh13_150
mod1=modS
source("PlotScenarioNf2.R")
dtab
source("MakeSummaryVar.R")
tabV
S1nf2_tabV<-tabV

probfreeS
S1nf2_probfreeS<-probfreeS
S1nf2_probfreeSLow<-probfreeSLow
S1nf2_probfreeSHigh<-probfreeSHigh

#Scenario No 2 #
modS<-mmKh13D_150
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
S1Dnf2_tabV<-tabV

probfreeS
S1Dnf2_probfreeS<-probfreeS
S1Dnf2_probfreeSLow<-probfreeSLow
S1Dnf2_probfreeSHigh<-probfreeSHigh

##########################################################################################
##########################################################################################
#Scenario for harvest strategies B
#mratio 1:3, 1:5, 1:10, 1:20
#Kadf
#hadf
#########################################################################################

Ktot=800 
hadf_max=0.141
hadf_m=hadf_max

h13
h0
h11=h13

####################################################################################
#mf Ratio 1:3
####################################################################################
mr3test.150<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h11,mfRatio=3,Tadf=150)
hcalf=0 #hyearling=0
mr3test.150_c00<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h0,mfRatio=3,Tadf=150)

#Scenario No 3
modS<-mr3test.150
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
mr3nf2_tabV<-tabV
probfreeS

mr3nf2_probfreeS<-probfreeS
mr3nf2_probfreeSLow<-probfreeSLow
mr3nf2_probfreeSHigh<-probfreeSHigh

#Scenario No 4
modS<-mr3test.150_c00
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
mr3c00nf2_tabV<-tabV

probfreeS
mr3c00nf2_probfreeS<-probfreeS
mr3c00nf2_probfreeSLow<-probfreeSLow
mr3c00nf2_probfreeSHigh<-probfreeSHigh

####################################################################################
#mf Ratio 1:5
####################################################################################
mr5test.150<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h11,mfRatio=5,Tadf=150)
hcalf=0 #hyearling=0
mr5test.150_c00<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h0,mfRatio=5,Tadf=150)

#Scenario No 5
modS<-mr5test.150
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
mr5nf2_tabV<-tabV
probfreeS

mr5nf2_probfreeS<-probfreeS
mr5nf2_probfreeSLow<-probfreeSLow
mr5nf2_probfreeSHigh<-probfreeSHigh

#Scenario no 6
modS<-mr5test.150_c00
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
mr5c00nf2_tabV<-tabV

probfreeS
mr5c00nf2_probfreeS<-probfreeS
mr5c00nf2_probfreeSLow<-probfreeSLow
mr5c00nf2_probfreeSHigh<-probfreeSHigh

####################################################################################
#mfRatio 1:10
####################################################################################
mr10test.150<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h11,mfRatio=10,Tadf=150)
hcalf=0
mr10test.150_c00<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h0,mfRatio=10,Tadf=150)

#Scenario No 7
modS<-mr10test.150
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
mr10nf2_tabV<-tabV

probfreeS
mr10nf2_probfreeS<-probfreeS
mr10nf2_probfreeSLow<-probfreeSLow
mr10nf2_probfreeSHigh<-probfreeSHigh


#Scenario No 8
modS<-mr10test.150_c00
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
mr10c00nf2_tabV<-tabV

probfreeS
mr10c00nf2_probfreeS<-probfreeS
mr10c00nf2_probfreeSLow<-probfreeSLow
mr10c00nf2_probfreeSHigh<-probfreeSHigh


####################################################################################
#mf Ratio 1:20
####################################################################################
mr20test.150<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_h,Nsim=1000,H1tab=h11,mfRatio=20,Tadf=150)
hcalf=0
mr20test.150_c00<-MakeSimPopTabRKh18(Simff=SimPop18RKtot_hcalf,Nsim=1000,H1tab=h0,mfRatio=20,Tadf=150)

#Scenario9
modS<-mr20test.150
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
mr20nf2_tabV<-tabV

probfreeS
mr20nf2_probfreeS<-probfreeS
mr20nf2_probfreeSLow<-probfreeSLow
mr20nf2_probfreeSHigh<-probfreeSHigh

#Scenario No 10
modS<-mr20test.150_c00
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
tabV
mr20c00nf2_tabV<-tabV

probfreeS
mr20c00nf2_probfreeS<-probfreeS
mr20c00nf2_probfreeSLow<-probfreeSLow
mr20c00nf2_probfreeSHigh<-probfreeSHigh


###################
#Effect RR #########
###################
pintro1=0.05    #0.1%
prior=0.5
pintro=0.001
nsim1=1000
####

modS0<-mr5test.150_c00
modS<-mmKh13_150 

#######################################
RR1=c(1,1,1)
source("Pfreedom_calculations_Nf2.R")
#######################################
#Scenario No 11
mod1=modS
source("PlotScenarioNf2.R")
probfreeS
S1Nf2_RR111_probfreeS<-probfreeS
S1Nf2_RR111_probfreeSLow<-probfreeSLow
S1Nf2_RR111_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
S1Nf2_RR111_tabV<-tabV


#Scenario No 16
mod1=modS0
source("PlotScenarioNf2.R")
mr5c00Nf2_RR111_probfreeS<-probfreeS
mr5c00Nf2_RR111_probfreeSLow<-probfreeSLow
mr5c00Nf2_RR111_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
mr5c00Nf2_RR111_tabV<-tabV

#######################################
RR1=c(1,2,2)
source("Pfreedom_calculations_Nf2.R")
######################################
#Scenario No 12
mod1=modS
source("PlotScenarioNf2.R")
#dtab
probfreeS
S1Nf2_RR122_probfreeS<-probfreeS
S1Nf2_RR122_probfreeSLow<-probfreeSLow
S1Nf2_RR122_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
S1Nf2_RR122_tabV<-tabV


#Scenario No 17
mod1=modS0
source("PlotScenarioNf2.R")
#dtab
probfreeS
mr5c00Nf2_RR122_probfreeS<-probfreeS
mr5c00Nf2_RR122_probfreeSLow<-probfreeSLow
mr5c00Nf2_RR122_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
mr5c00Nf2_RR122_tabV<-tabV

#######################################
RR1=c(1,2,4)
source("Pfreedom_calculations_Nf2.R")
######################################
#Scenario No 13
mod1=modS
source("PlotScenarioNf2.R")
#dtab
probfreeS
S1Nf2_RR124_probfreeS<-probfreeS
S1Nf2_RR124_probfreeSLow<-probfreeSLow
S1Nf2_RR124_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
S1Nf2_RR124_tabV<-tabV


#Scenario No 18
mod1=modS0
source("PlotScenarioNf2.R")
#dtab
probfreeS
mr5c00Nf2_RR124_probfreeS<-probfreeS
mr5c00Nf2_RR124_probfreeSLow<-probfreeSLow
mr5c00Nf2_RR124_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
#tabV
mr5c00Nf2_RR124_tabV<-tabV

#######################################
RR1=c(1,2,6)
simSSe_variant <- simSSe3L
source("Pfreedom_calculations_Nf2.R")
######################################
#Scenario No 1
mod1=modS
source("PlotScenarioNf2.R")
#dtab
probfreeS
S1Nf2_RR126_probfreeS<-probfreeS
S1Nf2_RR126_probfreeSLow<-probfreeSLow
S1Nf2_RR126_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
S1Nf2_RR126_tabV<-tabV


#Scenario No 6
mod1=modS0
source("PlotScenarioNf2.R")
#dtab
probfreeS
mr5c00Nf2_RR126_probfreeS<-probfreeS
mr5c00Nf2_RR126_probfreeSLow<-probfreeSLow
mr5c00Nf2_RR126_probfreeSHigh<-probfreeSHigh

source("MakeSummaryVar.R")
#tabV
mr5c00Nf2_RR126_tabV<-tabV


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
source("Pfreedom_calculations_Nf2.R")

resultPFree_Stoch1<-resultPFree
resultPFree_Stoch1 #Prob free for emp. data assuming stochastic RR version 1

PrRLN=pRLN2018
PrRLN
N.mean=as.numeric(N.mean18[1:6])
N.sd=as.numeric(N.sd18)
h13
#Function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioNf2.R 
calcPFreedom10<-calcPFreedom10_RRi 

#Scenario 14 (Harvest strategy No 1 with stochastic RR)#
modS<-mmKh13_150
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")

probfreeS
S1nf2_probfreeSRRi<-probfreeS
S1nf2_probfreeSRRiLow<-probfreeSLow
S1nf2_probfreeSRRiHigh<-probfreeSHigh


#Scenario 19 (harvest strategy No 6 with stochastic RR)
mod0<-mr5test.150_c00
mod1=mod0
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")

probfreeS
mr5c00nf2_probfreeSRRi<-probfreeS
mr5c00nf2_probfreeSRRiLow<-probfreeSLow
mr5c00nf2_probfreeSRRiHigh<-probfreeSHigh



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
source("Pfreedom_calculations_Nf2.R")
resultPFree_Stoch2<-resultPFree
resultPFree_Stoch2 ##Prob free for emp. data assuming stochastic RR version 1

PrRLN=pRLN2018
N.mean=as.numeric(N.mean18[1:6])
N.sd=as.numeric(N.sd18)

calcPFreedom10<-calcPFreedom10_RRi #calcPFreedom10_RRi

#Scenario 15 (Harvest strategy No 1 with stochastic RR)#
modS<-mmKh13_150
mod1=modS
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")
probfreeS
S1nf2_probfreeSRRi2<-probfreeS
S1nf2_probfreeSRRi2Low<-probfreeSLow
S1nf2_probfreeSRRi2High<-probfreeSHigh

S1nf2_probfreeSRRi2[5,]
S1nf2_probfreeS[5,]

#Scenario 20 (harvest strategy No 6 with stochastic RR)
mod0<-mr5test.150_c00
mod1=mod0
source("PlotScenarioNf2.R")
source("MakeSummaryVar.R")

probfreeS
mr5c00nf2_probfreeSRRi2<-probfreeS
mr5c00nf2_probfreeSRRi2Low<-probfreeSLow
mr5c00nf2_probfreeSRRi2High<-probfreeSHigh

mr5c00nf2_probfreeSRRi2[5,]
mr5c00nf2_probfreeS[5,]

#############################
##############################
#Back to standard level 
pintro1=0.05
pintro=0.001     #0.1%
RR1<-c(1,2,6)

simSSe_variant <- simSSe3L
source("Pfreedom_calculations_Nf2.R")

calcPFreedom10<-calcPFreedom10meanDSE_i #calcPFreedom10_RRi
PrRLN=pRLN2018

#############################
##Change import rates 
#############################
prior=0.5
nsim1=1000
pintro1=0.05
pintro=0.05 #Increase from basic 0.001

########################################################################################
#Function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioNf2.R 
calcPFreedom10<-calcPFreedom10meanDSE_i

modS<-mmKh13_150
modS0<-mr5test.150_c00

#Scenario No 21
mod1=modS
source("PlotScenarioNf2.R")
dtab
probfreeS
S1Nf2_probfreeS05<-probfreeS
S1Nf2_probfreeS05Low<-probfreeSLow
S1Nf2_probfreeS05High<-probfreeSHigh

#Scenario No 26
mod1=modS0
source("PlotScenarioNf2.R")
dtab
probfreeS
mr5c00Nf2_probfreeS05<-probfreeS
mr5c00Nf2_probfreeS05Low<-probfreeSLow
mr5c00Nf2_probfreeS05High<-probfreeSHigh

#Update with other combinations of pintro1 and pintro
prior=0.5
nsim1=1000
pintro1=0.05 #basic 
pintro=0.01 #basic 0.001
nsim1=1000


########################################################################################
#Function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioNf2.R 
calcPFreedom10<-calcPFreedom10meanDSE_i

modS<-mmKh13_150
modS0<-mr5test.150_c00

#Scenario No 22
mod1=modS
source("PlotScenarioNf2.R")
dtab
probfreeS
S1Nf2_probfreeS0501<-probfreeS
S1Nf2_probfreeS0501Low<-probfreeSLow
S1Nf2_probfreeS0501High<-probfreeSHigh

#Scenario No 27
mod1=modS0
source("PlotScenarioNf2.R")
dtab
probfreeS
mr5c00Nf2_probfreeS0501<-probfreeS
mr5c00Nf2_probfreeS0501Low<-probfreeSLow
mr5c00Nf2_probfreeS0501High<-probfreeSHigh



#Update with other combinations of pintro1 and pintro
prior=0.5
nsim1=1000
pintro1=0.05 #basic 
pintro=0.005 #basic 0.001
nsim1=1000

########################################################################################
#Function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioNf2.R 
calcPFreedom10<-calcPFreedom10meanDSE_i

modS<-mmKh13_150
modS0<-mr5test.150_c00

#Scenario No 23
mod1=modS
source("PlotScenarioNf2.R")
dtab
probfreeS
S1Nf2_probfreeS05005<-probfreeS
S1Nf2_probfreeS05005Low<-probfreeSLow
S1Nf2_probfreeS05005High<-probfreeSHigh

#Scenario No 28
mod1=modS0
source("PlotScenarioNf2.R")
dtab
probfreeS
mr5c00Nf2_probfreeS05005<-probfreeS
mr5c00Nf2_probfreeS05005Low<-probfreeSLow
mr5c00Nf2_probfreeS05005High<-probfreeSHigh



########################################################################
#Update with other combinations of pintro1 and pintro
prior=0.5
nsim1=1000
pintro1=0.01 #basic 
pintro=0.01 #basic 0.001
nsim1=1000

########################################################################################
#Function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioNf2.R 
calcPFreedom10<-calcPFreedom10meanDSE_i

modS<-mmKh13_150
modS0<-mr5test.150_c00

#Scenario No 24
mod1=modS
source("PlotScenarioNf2.R")
dtab
probfreeS
S1Nf2_probfreeS01<-probfreeS
S1Nf2_probfreeS01Low<-probfreeSLow
S1Nf2_probfreeS01High<-probfreeSHigh

#Scenario No 29
mod1=modS0
source("PlotScenarioNf2.R")
dtab
probfreeS
mr5c00Nf2_probfreeS01<-probfreeS
mr5c00Nf2_probfreeS01Low<-probfreeSLow
mr5c00Nf2_probfreeS01High<-probfreeSHigh



#Update with other combinations of pintro1 and pintro
prior=0.5
nsim1=1000
pintro1=0.01 #basic 
pintro=0.001 #basic 0.001
nsim1=1000


########################################################################################
#Function to be used in function MakeSummaryPstarScenarioVar in PlotScenarioNf2.R 
calcPFreedom10<-calcPFreedom10meanDSE_i

modS<-mmKh13_150
modS0<-mr5test.150_c00

#Scenario No 25
mod1=modS
source("PlotScenarioNf2.R")
dtab
probfreeS
S1Nf2_probfreeS01001<-probfreeS
S1Nf2_probfreeS01001Low<-probfreeSLow
S1Nf2_probfreeS01001High<-probfreeSHigh

#Scenario No 30
mod1=modS0
source("PlotScenarioNf2.R")
dtab
probfreeS
mr5c00Nf2_probfreeS01001<-probfreeS
mr5c00Nf2_probfreeS01001Low<-probfreeSLow
mr5c00Nf2_probfreeS01001High<-probfreeSHigh




######################################################
################################
#Back to standard level
################################
pintro1=0.05
pintro=0.001     #0.1%
prior=0.5
nsim1=1000


#####################################################
#source(SummaryResTabNf2.R)
#####################################################
