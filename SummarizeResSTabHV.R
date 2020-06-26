######################################
#Area HV         #####################
######################################

##########################################################################################
#MakeRowSTab: Function for selecting summary output of each scenario
#Output from scenarios are saved as tables (design prevelence in rows and year in columns): 
#probtab = xxxx_probfreeS: Table with mean probability of freedom after each year of testing
#probtabLow = xxxx_probfreeSLow: Table with 2.5 percentile
#probtabHigh = xxxx_probfreeSHigh: Table with 97.5 percentile
#pstar4: rowname of the row with results for the design prevalence used in area HV
#sdtab = xxxx_tabv: Table with outcome (mean numbers) from population model simulations as 
#specified in rownames for each year (columns)
##########################################################################################
MakeRowSTab4<-function(probtab,probtabLow,probtabHigh,sdtab)
  {
  c(probtab["pstar4",c("y2018")],probtabLow["pstar4",c("y2018")],probtabHigh["pstar4",c("y2018")],
  probtab["pstar4",c("y2022")],probtabLow["pstar4",c("y2022")],probtabHigh["pstar4",c("y2022")],
  min(which(probtab["pstar4",]>0.900))-2, #start at year 2018 
  min(which(probtab["pstar4",]>0.990))-2, #start at year 2018 
  min(which(round(probtab["pstar4",],2)>=0.99))-2, #start at year 2018 
  sdtab["Nharvested_c","2018"],
  sum(sdtab["Nharvested_c",1:5]),
  sdtab["postHarvest_Nfemales","2018"],
  sdtab["postHarvestpAdf","2018"],
  sdtab["postHarvestpAdm","2018"],
  sdtab["postHarvestCalves","2018"],
  sdtab["PopSize","2018"]-sdtab["Nharvested","2018"],
  sdtab["PopSize","2022"]-sdtab["Nharvested","2022"],
  (sdtab["PopSize_c","2018"]-sdtab["Nharvested_c","2018"]),
  (sdtab["PopSize_c","2022"]-sdtab["Nharvested_c","2022"]))
}


#####

S2hv_probfreeS
name="S2hv_"
#probtab=S1nf2_probfreeS
#probtabLow=S1nf2_probfreeSLow
#probtabHigh=S1nf2_probfreeSHigh
#sdtab<-S1nf2_tabV
sumdatS2hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS2hv

S2Dhv_probfreeS
name="S2Dhv_"
sumdatS2Dhv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS2Dhv

S2Dtadf2800hv_probfreeS
name="S2Dtadf2800hv_"
sumdatS2Dhv2800<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS2Dhv2800


S2Dtadf2500hv_probfreeS
name="S2Dtadf2500hv_"
sumdatS2Dhv2500<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS2Dhv2500


mr3hv_probfreeS
name="mr3hv_"
sumdatmr3hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr3hv

mr3c00hv_probfreeS
name="mr3c00hv_"
sumdatmr3c00hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr3c00hv

mr5hv_probfreeS
name="mr5hv_"
sumdatmr5hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5hv

mr5c00hv_probfreeS
name="mr5c00hv_"
sumdatmr5c00hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5c00hv

mr10hv_probfreeS
name="mr10hv_"
sumdatmr10hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr10hv

mr10c00hv_probfreeS
name="mr10c00hv_"
sumdatmr10c00hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr10c00hv

mr20hv_probfreeS
name="mr20hv_"
sumdatmr20hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr20hv

mr20c00hv_probfreeS
name="mr20c00hv_"
sumdatmr20c00hv<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr20c00hv


#####

Hv_RR111_probfreeS
name="Hv_RR111_"

sumdatHvRR111<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatHvRR111

Hv_RR122_probfreeS
name="Hv_RR122_"

sumdatHvRR122<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatHvRR122

Hv_RR124_probfreeS
name="Hv_RR124_"

sumdatHvRR124<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatHvRR124

###

mr5Hv_RR111_probfreeS
name="mr5Hv_RR111_"

sumdatmr5hvRR111<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5hvRR111

mr5Hv_RR122_probfreeS
name="mr5Hv_RR122_"

sumdatmr5hvRR122<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5hvRR122

mr5Hv_RR124_probfreeS
name="mr5Hv_RR124_"

sumdatmr5hvRR124<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5hvRR124

S2hv_probfreeSRRi
name="S2hv_"
tabV=S2hv_tabV
sumdatS2hv_RRi<-MakeRowSTab4(probtab=get(paste0(name,"probfreeSRRi")),probtabLow=get(paste0(name,"probfreeSRRiLow")),probtabHigh=get(paste0(name,"probfreeSRRiHigh")),sdtab=tabV)
sumdatS2hv_RRi

mr5c00hv_probfreeSRRi
name="mr5c00hv_"
tabV=mr5c00hv_tabV
sumdatmr5c00hv_RRi<-MakeRowSTab4(probtab=get(paste0(name,"probfreeSRRi")),probtabLow=get(paste0(name,"probfreeSRRiLow")),probtabHigh=get(paste0(name,"probfreeSRRiHigh")),sdtab=tabV)
sumdatmr5c00hv_RRi

S2hv_probfreeSRRi2
name="S2hv_"
tabV=S2hv_tabV
sumdatS2hv_RRi2<-MakeRowSTab4(probtab=get(paste0(name,"probfreeSRRi2")),probtabLow=get(paste0(name,"probfreeSRRi2Low")),probtabHigh=get(paste0(name,"probfreeSRRi2High")),sdtab=tabV)
sumdatS2hv_RRi2

mr5c00hv_probfreeSRRi2
name="mr5c00hv_"
tabV=mr5c00hv_tabV
sumdatmr5c00hv_RRi2<-MakeRowSTab4(probtab=get(paste0(name,"probfreeSRRi2")),probtabLow=get(paste0(name,"probfreeSRRi2Low")),probtabHigh=get(paste0(name,"probfreeSRRi2High")),sdtab=tabV)
sumdatmr5c00hv_RRi2


Hv_probfreeS05001
name="Hv_"
tabV=S2hv_tabV
sumdatS2hv_05001<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS05001")),probtabLow=get(paste0(name,"probfreeS05001Low")),probtabHigh=get(paste0(name,"probfreeS05001High")),sdtab=tabV)
sumdatS2hv_05001

mr5Hv_probfreeS05001
name="mr5Hv_"
tabV=mr5c00hv_tabV
sumdatmr5c00hv_05001<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS05001")),probtabLow=get(paste0(name,"probfreeS05001Low")),probtabHigh=get(paste0(name,"probfreeS05001High")),sdtab=tabV)
sumdatmr5c00hv_05001

Hv_probfreeS0501
name="Hv_"
tabV=S2hv_tabV
sumdatS2hv_0501<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS0501")),probtabLow=get(paste0(name,"probfreeS0501Low")),probtabHigh=get(paste0(name,"probfreeS0501High")),sdtab=tabV)
sumdatS2hv_0501

mr5Hv_probfreeS0501
name="mr5Hv_"
tabV=mr5c00hv_tabV
sumdatmr5c00hv_0501<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS0501")),probtabLow=get(paste0(name,"probfreeS0501Low")),probtabHigh=get(paste0(name,"probfreeS0501High")),sdtab=tabV)
sumdatmr5c00hv_0501

Hv_probfreeS01
name="Hv_"
tabV=S2hv_tabV
sumdatS2hv_01<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS01")),probtabLow=get(paste0(name,"probfreeS01Low")),probtabHigh=get(paste0(name,"probfreeS01High")),sdtab=tabV)
sumdatS2hv_01

mr5Hv_probfreeS01
name="mr5Hv_"
tabV=mr5c00hv_tabV
sumdatmr5c00hv_01<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS01")),probtabLow=get(paste0(name,"probfreeS01Low")),probtabHigh=get(paste0(name,"probfreeS01High")),sdtab=tabV)
sumdatmr5c00hv_01

Hv_probfreeS05005
name="Hv_"
tabV=S2hv_tabV
sumdatS2hv_05005<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS05005")),probtabLow=get(paste0(name,"probfreeS05005Low")),probtabHigh=get(paste0(name,"probfreeS05005High")),sdtab=tabV)
sumdatS2hv_05005

mr5Hv_probfreeS05005
name="mr5Hv_"
tabV=mr5c00hv_tabV
sumdatmr5c00hv_05005<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS05005")),probtabLow=get(paste0(name,"probfreeS05005Low")),probtabHigh=get(paste0(name,"probfreeS05005High")),sdtab=tabV)
sumdatmr5c00hv_05005


Hv_probfreeS05
name="Hv_"
tabV=S2hv_tabV
sumdatS2hv_05<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS05")),probtabLow=get(paste0(name,"probfreeS05Low")),probtabHigh=get(paste0(name,"probfreeS05High")),sdtab=tabV)
sumdatS2hv_05

mr5Hv_probfreeS05
name="mr5Hv_"
tabV=mr5c00hv_tabV
sumdatmr5c00hv_05<-MakeRowSTab4(probtab=get(paste0(name,"probfreeS05")),probtabLow=get(paste0(name,"probfreeS05Low")),probtabHigh=get(paste0(name,"probfreeS05High")),sdtab=tabV)
sumdatmr5c00hv_05

#####

#Table with one row for each harvest scenario /strategy
HVSUMTAB<-rbind(sumdatS2hv,sumdatS2Dhv,sumdatS2Dhv2800,sumdatS2Dhv2500,
                sumdatmr3hv,sumdatmr3c00hv,
                sumdatmr5hv,sumdatmr5c00hv,
                sumdatmr10hv,sumdatmr10c00hv,
                sumdatmr20hv,sumdatmr20c00hv,
                sumdatHvRR111,sumdatHvRR122,sumdatHvRR124,
                sumdatS2hv_RRi,
                sumdatS2hv_RRi2,
                sumdatmr5hvRR111,sumdatmr5hvRR122,sumdatmr5hvRR124,
                sumdatmr5c00hv_RRi,
                sumdatmr5c00hv_RRi2,
                sumdatS2hv_05,
                sumdatS2hv_0501,
                sumdatS2hv_05005,
                sumdatS2hv_05001,
                sumdatS2hv_01,
                sumdatmr5c00hv_05,
                sumdatmr5c00hv_0501,
                sumdatmr5c00hv_05005,
                sumdatmr5c00hv_05001,
                sumdatmr5c00hv_01
) 


colnames(HVSUMTAB)<-c("PropFree1","PFree1_2.5%","PFree1_97.5%",
                      "PropFree5","PFree5_2.5%","PFree5_97.5%","NrYear_90.0%","NrYear_99.0%","NrYear_99%",
                      "NrTested1yr","NrTested5yr","NAdf","pAdf","pAdm","pCalves","Popsize1","Popsize5",
                      "PopSize1_c","PopSize5_c")
#HVSUMTAB
#write.csv2(HVSUMTAB,file="HV_summarySTab2&3.csv")

