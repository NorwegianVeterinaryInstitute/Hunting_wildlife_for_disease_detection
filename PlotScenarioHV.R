#mod1: output from population simulation
#assembleDatA: function specified in  SimFunctionsPaper.R
#calcPFreedom10meanDSE_i: function specified in SimFunctionsPaper.R
#MakeSummaryPstarScenarioVa: function from "SimFunctionsPaper.R"
#MakeProbPlot: function from SimFunctionsPaper.R

#dtab<-samleDatA(mod1)
#Assemble output data (mod1) from population simulation model
dtab<-assembleDatA(mod1)

#Specify function used in MakeSummaryPstarScenarioVar.R
#calcPFreedom10<-calcPFreedom10meanDSE_i
#calcPFreedom10<-calcPFreedom10_i #this alternative includes more variation in dSe and takes longer time

#MakeSummaryPstarScenarioVar - function from "SimFunctionsPaper.R"
#Calculate probability of freedom for each simulated year and
#Summarize the frequency ditstributions by mean, 2.5 and 97.5 percentiles
pfreemodSV<-MakeSummaryPstarScenarioVar(mod1)

probfreeS=round(pfreemodSV$Pfree.pstarS,3)
probfreeSLow=round(pfreemodSV$Pfree.pstarSLow,3)
probfreeSHigh=round(pfreemodSV$Pfree.pstarSHigh,3)

par(mfrow=c(1,1))
MakeProbPlot(probfreeS) #Function in SimFunctionsPaper.R



