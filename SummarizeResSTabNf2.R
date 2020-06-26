######################################
#Nf2          #####################
######################################

##########################################################################################
#MakeRowSTab: Function for selecting summary output of each scenario
#Output from scenarios are saved as tables (design prevelence in rows and year in columns): 
#probtab = xxxx_probfreeS: Table with mean probability of freedom after each year of testing
#probtabLow = xxxx_probfreeSLow: Table with 2.5 percentile
#probtabHigh = xxxx_probfreeSHigh: Table with 97.5 percentile
#pstarSA: rowname of the row with results for the design prevalence used in area Nf2
#sdtab = xxxx_tabv: Table with outcome (mean numbers) from population model simulations as 
#specified in rownames for each year (columns)
##########################################################################################
MakeRowSTab<-function(probtab,probtabLow,probtabHigh,sdtab)
  {
  c(probtab["pstarSA",c("y2018")],probtabLow["pstarSA",c("y2018")],probtabHigh["pstarSA",c("y2018")],
  probtab["pstarSA",c("y2022")],probtabLow["pstarSA",c("y2022")],probtabHigh["pstarSA",c("y2022")],
  min(which(probtab["pstarSA",]>0.900))-2, #start at year 2018 
  min(which(probtab["pstarSA",]>0.990))-2, #start at year 2018 
  min(which(round(probtab["pstarSA",],2)>=0.99))-2, #start at year 2018 
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

S1nf2_probfreeS
name="S1nf2_"
#probtab=S1nf2_probfreeS
#probtabLow=S1nf2_probfreeSLow
#probtabHigh=S1nf2_probfreeSHigh
#sdtab<-S1nf2_tabV
sumdatS1nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS1nf2

S1Dnf2_probfreeS
name="S1Dnf2_"
sumdatSD1nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatSD1nf2

mr3nf2_probfreeS
name="mr3nf2_"
sumdatmr3nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr3nf2

mr3c00nf2_probfreeS
name="mr3c00nf2_"
sumdatmr3c00nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr3c00nf2

mr5nf2_probfreeS
name="mr5nf2_"
sumdatmr5nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5nf2

mr5c00nf2_probfreeS
name="mr5c00nf2_"
sumdatmr5c00nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5c00nf2

mr10nf2_probfreeS
name="mr10nf2_"
sumdatmr10nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr10nf2

mr10c00nf2_probfreeS
name="mr10c00nf2_"
sumdatmr10c00nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr10c00nf2

mr20nf2_probfreeS
name="mr20nf2_"
sumdatmr20nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr20nf2

mr20c00nf2_probfreeS
name="mr20c00nf2_"
sumdatmr20c00nf2<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr20c00nf2


#####

S1Nf2_RR111_probfreeS
name="S1Nf2_RR111_"

sumdatS1nf2RR111<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS1nf2RR111

S1Nf2_RR122_probfreeS
name="S1Nf2_RR122_"

sumdatS1nf2RR122<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS1nf2RR122

S1Nf2_RR124_probfreeS
name="S1Nf2_RR124_"

sumdatS1nf2RR124<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS1nf2RR124

S1Nf2_RR126_probfreeS
name="S1Nf2_RR126_"

sumdatS1nf2RR126<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatS1nf2RR126


mr5c00Nf2_RR111_probfreeS
name="mr5c00Nf2_RR111_"

sumdatmr5c00nf2RR111<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5c00nf2RR111

mr5c00Nf2_RR122_probfreeS
name="mr5c00Nf2_RR122_"

sumdatmr5c00nf2RR122<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5c00nf2RR122

mr5c00Nf2_RR124_probfreeS
name="mr5c00Nf2_RR124_"

sumdatmr5c00nf2RR124<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5c00nf2RR124

mr5c00Nf2_RR126_probfreeS
name="mr5c00Nf2_RR126_"

sumdatmr5c00nf2RR126<-MakeRowSTab(probtab=get(paste0(name,"probfreeS")),probtabLow=get(paste0(name,"probfreeSLow")),probtabHigh=get(paste0(name,"probfreeSHigh")),sdtab=get(paste0(name,"tabV")))
sumdatmr5c00nf2RR126

#####

S1Nf2_probfreeS05
name="S1Nf2_"
tabV=S1nf2_tabV
sumdatS1nf2_05<-MakeRowSTab(probtab=get(paste0(name,"probfreeS05")),probtabLow=get(paste0(name,"probfreeS05Low")),probtabHigh=get(paste0(name,"probfreeS05High")),sdtab=tabV)
sumdatS1nf2_05

mr5c00Nf2_probfreeS05
name="mr5c00Nf2_"
tabV=mr5c00nf2_tabV
sumdatmr5c00nf2_05<-MakeRowSTab(probtab=get(paste0(name,"probfreeS05")),probtabLow=get(paste0(name,"probfreeS05Low")),probtabHigh=get(paste0(name,"probfreeS05High")),sdtab=tabV)
sumdatmr5c00nf2_05

#####


S1Nf2_probfreeS01
name="S1Nf2_"
tabV=S1nf2_tabV
sumdatS1nf2_01<-MakeRowSTab(probtab=get(paste0(name,"probfreeS01")),probtabLow=get(paste0(name,"probfreeS01Low")),probtabHigh=get(paste0(name,"probfreeS01High")),sdtab=tabV)
sumdatS1nf2_01

mr5c00Nf2_probfreeS01
name="mr5c00Nf2_"
tabV=mr5c00nf2_tabV
sumdatmr5c00nf2_01<-MakeRowSTab(probtab=get(paste0(name,"probfreeS01")),probtabLow=get(paste0(name,"probfreeS01Low")),probtabHigh=get(paste0(name,"probfreeS01High")),sdtab=tabV)
sumdatmr5c00nf2_01


S1Nf2_probfreeS01001
name="S1Nf2_"
tabV=S1nf2_tabV
sumdatS1nf2_01001<-MakeRowSTab(probtab=get(paste0(name,"probfreeS01001")),probtabLow=get(paste0(name,"probfreeS01001Low")),probtabHigh=get(paste0(name,"probfreeS01001High")),sdtab=tabV)
sumdatS1nf2_01001

mr5c00Nf2_probfreeS01001
name="mr5c00Nf2_"
tabV=mr5c00nf2_tabV
sumdatmr5c00nf2_01001<-MakeRowSTab(probtab=get(paste0(name,"probfreeS01001")),probtabLow=get(paste0(name,"probfreeS01001Low")),probtabHigh=get(paste0(name,"probfreeS01001High")),sdtab=tabV)
sumdatmr5c00nf2_01001


S1Nf2_probfreeS0501
name="S1Nf2_"
tabV=S1nf2_tabV
sumdatS1nf2_0501<-MakeRowSTab(probtab=get(paste0(name,"probfreeS0501")),probtabLow=get(paste0(name,"probfreeS0501Low")),probtabHigh=get(paste0(name,"probfreeS0501High")),sdtab=tabV)
sumdatS1nf2_0501

mr5c00Nf2_probfreeS0501
name="mr5c00Nf2_"
tabV=mr5c00nf2_tabV
sumdatmr5c00nf2_0501<-MakeRowSTab(probtab=get(paste0(name,"probfreeS0501")),probtabLow=get(paste0(name,"probfreeS0501Low")),probtabHigh=get(paste0(name,"probfreeS0501High")),sdtab=tabV)
sumdatmr5c00nf2_0501


S1Nf2_probfreeS05005
name="S1Nf2_"
tabV=S1nf2_tabV
sumdatS1nf2_05005<-MakeRowSTab(probtab=get(paste0(name,"probfreeS05005")),probtabLow=get(paste0(name,"probfreeS05005Low")),probtabHigh=get(paste0(name,"probfreeS05005High")),sdtab=tabV)
sumdatS1nf2_05005

mr5c00Nf2_probfreeS05005
name="mr5c00Nf2_"
tabV=mr5c00nf2_tabV
sumdatmr5c00nf2_05005<-MakeRowSTab(probtab=get(paste0(name,"probfreeS05005")),probtabLow=get(paste0(name,"probfreeS05005Low")),probtabHigh=get(paste0(name,"probfreeS05005High")),sdtab=tabV)
sumdatmr5c00nf2_05005


S1nf2_probfreeSRRi
name="S1nf2_"
tabV=S1nf2_tabV
sumdatS1nf2_RRi<-MakeRowSTab(probtab=get(paste0(name,"probfreeSRRi")),probtabLow=get(paste0(name,"probfreeSRRiLow")),probtabHigh=get(paste0(name,"probfreeSRRiHigh")),sdtab=tabV)
sumdatS1nf2_RRi

mr5c00nf2_probfreeSRRi
name="mr5c00nf2_"
tabV=mr5c00nf2_tabV
sumdatmr5c00nf2_RRi<-MakeRowSTab(probtab=get(paste0(name,"probfreeSRRi")),probtabLow=get(paste0(name,"probfreeSRRiLow")),probtabHigh=get(paste0(name,"probfreeSRRiHigh")),sdtab=tabV)
sumdatmr5c00nf2_RRi


S1nf2_probfreeSRRi2
name="S1nf2_"
tabV=S1nf2_tabV
sumdatS1nf2_RRi2<-MakeRowSTab(probtab=get(paste0(name,"probfreeSRRi2")),probtabLow=get(paste0(name,"probfreeSRRi2Low")),probtabHigh=get(paste0(name,"probfreeSRRi2High")),sdtab=tabV)
sumdatS1nf2_RRi2

mr5c00nf2_probfreeSRRi2
name="mr5c00nf2_"
tabV=mr5c00nf2_tabV
sumdatmr5c00nf2_RRi2<-MakeRowSTab(probtab=get(paste0(name,"probfreeSRRi2")),probtabLow=get(paste0(name,"probfreeSRRi2Low")),probtabHigh=get(paste0(name,"probfreeSRRi2High")),sdtab=tabV)
sumdatmr5c00nf2_RRi2

#####
#Table with one row for each harvest scenario /strategy
Nf2SUMTAB<-rbind(sumdatS1nf2,sumdatSD1nf2,
                 sumdatmr3nf2,sumdatmr3c00nf2,
                 sumdatmr5nf2,sumdatmr5c00nf2,
                 sumdatmr10nf2,sumdatmr10c00nf2,
                 sumdatmr20nf2,sumdatmr20c00nf2,
                 sumdatS1nf2RR111,sumdatS1nf2RR122,sumdatS1nf2RR124,#sumdatS1nf2RR126,
                 sumdatS1nf2_RRi,sumdatS1nf2_RRi2,
                 sumdatmr5c00nf2RR111,sumdatmr5c00nf2RR122,sumdatmr5c00nf2RR124,#sumdatmr5c00nf2RR126,
                 sumdatmr5c00nf2_RRi,sumdatmr5c00nf2_RRi2,
                 sumdatS1nf2_05,
                 sumdatS1nf2_0501,
                 sumdatS1nf2_05005,
                 sumdatS1nf2_01,
                 sumdatS1nf2_01001,
                 sumdatmr5c00nf2_05,
                 sumdatmr5c00nf2_0501,
                 sumdatmr5c00nf2_05005,
                 sumdatmr5c00nf2_01,
                 sumdatmr5c00nf2_01001
                 ) 



colnames(Nf2SUMTAB)<-c("PropFree1","PFree1_2.5%","PFree1_97.5%",
                    "PropFree5","PFree5_2.5%","PFree5_97.5%","NrYear_90.0%","NrYear_99.0%","NrYear_99%",
                    "NrTested1yr","NrTested5yr","NAdf","pAdf","pAdm","pCalves","Popsize1","Popsize5",
                    "PopSize1_c","PopSize5_c")

#write.csv2(Nf2SUMTAB,file="Nf2_summarySTab2&3.csv")



