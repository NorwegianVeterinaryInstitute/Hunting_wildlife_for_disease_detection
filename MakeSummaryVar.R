#summarize output data from dtab, which is a table made by 
#"PlotScenario.R" from mod1 (output of population simulation model)

#Number of adult males in pop after hunt:
postHarvest_Nmales<-as.numeric(dtab["Nadm",c(1:8)])-as.numeric(dtab["Nadm2",c(1:8)])
postHarvest_Nmales
postHarvest_Nfemales<-as.numeric(dtab["Nadf",c(1:8)])-as.numeric(dtab["Nadf2",c(1:8)])
postHarvest_Nfemales

#Tot pop-harvest population size (included calves)
Npop<-dtab[dtab$Var=="Nmean1",c(1:8)]
Hpop<-dtab[dtab$Var=="Htab1",c(1:8)]
postHarvestPop<-apply(apply(Npop,2,as.numeric),2,sum)-apply(apply(Hpop,2,as.numeric),2,sum)
postHarvest_Nmales/postHarvestPop

#Number of adult males hunted:
dtab["Nadm2",]

#Poportion of adult males in postharves pop - excluding calves
postHarvestPop_c<-apply(apply(Npop[-c(1:2),],2,as.numeric),2,sum)-apply(apply(Hpop[-c(1:2),],2,as.numeric),2,sum)
postHarvest_Nmales/postHarvestPop_c
postHarvestpAdm<-round(postHarvest_Nmales/postHarvestPop,3)
postHarvestpAdf<-round(postHarvest_Nfemales/postHarvestPop,3)

postHarvestCalvesN<-as.numeric(dtab["N0fm",c(1:8)])+as.numeric(dtab["N0f",c(1:8)])-as.numeric(dtab["N0f2",c(1:8)])-as.numeric(dtab["N0fm2",c(1:8)])
postHarvestCalves<-round(postHarvestCalvesN/postHarvestPop,3)

PopSize<-apply(apply(Npop,2,as.numeric),2,sum)
PopSize_c<-apply(apply(Npop[-c(1:2),],2,as.numeric),2,sum)
Nharvested<-apply(apply(Hpop,2,as.numeric),2,sum)
Nharvested_c<-apply(apply(Hpop[-c(1:2),],2,as.numeric),2,sum)

tabV<-rbind(PopSize,PopSize_c,Nharvested,Nharvested_c,postHarvest_Nfemales,postHarvest_Nmales,postHarvestCalvesN,postHarvestCalves,postHarvestpAdf,postHarvestpAdm)
#tabV

