
probfreeTabSA5<-rbind(S1nf2_probfreeS["pstarSA",],
                      mr3c00nf2_probfreeS["pstarSA",],
                      mr5c00nf2_probfreeS["pstarSA",],
                      mr10c00nf2_probfreeS["pstarSA",],
                      mr20c00nf2_probfreeS["pstarSA",]
)

probfreeTabS4HV<-rbind(S2hv_probfreeS["pstar4",],
                       mr3c00hv_probfreeS["pstar4",],
                       mr5c00hv_probfreeS["pstar4",],
                       mr10c00hv_probfreeS["pstar4",],
                       mr20c00hv_probfreeS["pstar4",]
)




MakeProbPlotNfS4t5<-function(pfreeSTab,ptested1=testedNf[1:4],ptested2=testedNf[5:6]){
  plot(seq(2015,2028,1),c(0.5,pfreeSTab[1,])*100,las=1,type="l",col="black",lwd=2,ylab="Probability of freedom (%)",
       xlab="Year",ylim=c(50,100))
  abline(h=99,lty=2)
  lines(seq(2015,2028,1),c(0.5,pfreeSTab[2,])*100,type="l",lwd=2,col="green")
  lines(seq(2015,2028,1),c(0.5,pfreeSTab[3,])*100,type="l",col="skyblue2",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeSTab[4,])*100,type="l",col="blue",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeSTab[5,])*100,type="l",col="purple",lwd=2)
  lines(seq(2015,2018,1),c(ptested1[1:4])*100,type="l",lwd=1.5,col="red",lty=1)
  lines(c(2018,2018.4),c(ptested1[4],ptested2[1])*100,type="l",lwd=1.5,col="red",lty=1)
  lines(c(2018,2018.4,2019),c(ptested1[4],ptested2[1:2])*100,type="l",lwd=1.5,col="red",lty=1)
  
  points(seq(2016,2018,1),c(ptested1[-1])*100,lwd=2,pch=19,col="red")
  points(c(2018,2018.4),c(ptested1[4],ptested2[1])*100,col="red",lwd=2,pch=19)
  
  points(c(2018.4,2019),c(ptested2)*100,col="red",lwd=2,pch=19)
  lines(seq(2015,2017,1),c(0.5,pfreeSTab[1,1:2])*100,type="l",lwd=2,col="black")
  legend("bottomright",
         c("emp. data","ordinary","m:f=1:3","m:f=1:5","m:f=1:10","m:f=1:20"),bty="n",#lwd=2,lty=1,
         lwd=c(1,2,2,2,2,2),lty=c(1,1,1,1,1,1),
         col=c("red","black","green","skyblue2","blue","purple"),pch=c(16,16,16,16,16,16),cex=1.0)#0.9)
}



MakeProbPlotS4t<-function(pfreeSTab){
  plot(seq(2015,2028,1),c(0.5,pfreeSTab[1,])*100,type="l",col="black",lwd=2,las=1,ylab="Probability of freedom (%)",xlab="Year",ylim=c(50,100))
  abline(h=99,lty=2)
  lines(seq(2015,2028,1),c(0.5,pfreeSTab[2,])*100,type="l",lwd=2,col="green")
  lines(seq(2015,2028,1),c(0.5,pfreeSTab[3,])*100,type="l",col="skyblue2",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeSTab[4,])*100,type="l",col="blue",lwd=2)
  lines(seq(2015,2028,1),c(0.5,pfreeSTab[5,])*100,type="l",col="purple",lwd=2)
  lines(seq(2015,2017,1),c(0.5,pfreeSTab[1,1:2])*100,type="l",col="black",lwd=2)
  
  
  points(seq(2016,2019,1),c(testedHV[2:5])*100,pch=19,lwd=2,col="red")
  lines(seq(2017,2019,1),c(testedHV[3:5])*100,lwd=1.5,col="red",lty=1)
  
    legend("bottomright",
         c("emp. data","ordinary","m:f=1:3","m:f=1:5","m:f=1:10","m:f=1:20"),bty="n",lwd=c(1,2,2,2,2,2),lty=c(1,1,1,1,1,1),
         col=c("red","black","green","skyblue2","blue","purple"),pch=c(16,16,16,16,16,16),cex=1)#0.9)
}




######################
#Figure 3AB ###########
######################
par(mfrow=c(2,1),mai = c(0.6, 1, 0.1, 0.6))

#Figur 1 Paper
MakeProbPlotNfS4t5
MakeProbPlotNfS4t5(probfreeTabSA5)
segments(x0=c(2015,2016,2017,2018,2018.4,2019), y0=testedNfLow*100,y1=testedNfHigh*100,col="red",lwd=2)
#segments(x0=c(2016,2017,2018,2018.4,2019), y0=resultPFree_Stoch1[,2]*100,y1=resultPFree_Stoch1[,4]*100,col="red",lwd=1)

resultPFree_Stoch1


legend(2019.7,81,"Nordfjella zone 2",bty="n")
MakeProbPlotS4t(probfreeTabS4HV)
legend(2019.8,81,"Hardangervidda",bty="n")
segments(x0=c(2015,2016,2017,2018,2019), y0=testedHVLow*100,y1=testedHVHigh*100,col="red",lwd=2)
#segments(x0=c(2016,2017,2018,2019), y0=resultPFreeHV_Stoch1[,2]*100,y1=resultPFreeHV_Stoch1[,4]*100,col="red",lwd=1)


######################
#Figure 4AB ########### 
######################

#Data Nordfjella zone 2
Nmat2 #pre-harvest population
Nhunted #number hunted
StabNf2 #number tested

S18 <-as.numeric(StabNf2[StabNf2[,1]=="S18",c(2:4)])
#Including extra harvested winter 2019
S18w19<-S18 +as.numeric(StabNf2[StabNf2[,1]=="S1819W",c(2:4)])
S19 <-as.numeric(StabNf2[StabNf2[,1]=="S19",c(2:4)])
#Set number of NA ageclass to yearling
S19[1]<-S19[1]+as.numeric(StabNf2[StabNf2[,1]=="S19","ageNA"])
S19
#Number harvested = number sampled + 29 calves
H19=sum(S19)+29

#Proportion of adult males in post-harvest population
padm15=round((Nmat2["N15","Nadm"]-Nhunted["H15","Hadm"])/(Nmat2["N15","Total"]-sum(Nhunted["H15",])),3)
padm16=round((Nmat2["N16","Nadm"]-Nhunted["H16","Hadm"])/(Nmat2["N16","Total"]-sum(Nhunted["H16",])),3)
padm17=round((Nmat2["N17","Nadm"]-Nhunted["H17","Hadm"])/(Nmat2["N17","Total"]-sum(Nhunted["H17",])),3)
padm18=round((Nmat2["N18","Nadm"]-Nhunted["H18","Hadm"])/(Nmat2["N18","Total"]-sum(Nhunted["H18",])),3)
padm18w=round((Nmat2["N18","Nadm"]-S18w19[3])/(Nmat2["N18","Total"]-(sum(Nhunted["H18",])+52)),3)
padm19=round((Nmat2["N19","Nadm"]-S19[3])/(Nmat2["N19","Total"]-H19),3)
padmNf2<-c(padm15,padm16,padm17,padm18,padm18w,padm19)

#Number of adult females in post-harvest population
adf15=Nmat2["N15","Nadf"]-Nhunted["H15","Hadf"]
adf16=Nmat2["N16","Nadf"]-Nhunted["H16","Hadf"]
adf17=Nmat2["N17","Nadf"]-Nhunted["H17","Hadf"]
adf18=Nmat2["N18","Nadf"]-Nhunted["H18","Hadf"]
adf18w=adf18-2
adf19=Nmat2["N19","Nadf"]-22
adfNf2<-c(adf15,adf16,adf17,adf18,adf18w,adf19)

#N pop post-harvest 
pN15=Nmat2["N15","Total"]-sum(Nhunted["H15",])
pN16=Nmat2["N16","Total"]-sum(Nhunted["H16",])
pN17=Nmat2["N17","Total"]-sum(Nhunted["H17",])
pN18=Nmat2["N18","Total"]-sum(Nhunted["H18",])
pN18w=Nmat2["N18","Total"]-(sum(Nhunted["H18",])+52)
pN19=Nmat2["N19","Total"]-H19
pN_Nf2<-c(pN15,pN16,pN17,pN18,pN18w,pN19)


#For Hvidda:
Nmat #pre-harvest population
NhuntedHV #number hunted

#Proportion of adult males in post-harvest population
padm15HV=round((Nmat["N15","Nadm"]-NhuntedHV["H15","Hadm"])/(Nmat["N15","Total"]-sum(NhuntedHV["H15",])),3)
padm16HV=round((Nmat["N16","Nadm"]-NhuntedHV["H16","Hadm"])/(Nmat["N16","Total"]-sum(NhuntedHV["H16",])),3)
padm17HV=round((Nmat["N17","Nadm"]-NhuntedHV["H17","Hadm"])/(Nmat["N17","Total"]-sum(NhuntedHV["H17",])),3)
padm18HV=round((Nmat["N18","Nadm"]-NhuntedHV["H18","Hadm"])/(Nmat["N18","Total"]-sum(NhuntedHV["H18",])),3)
padm19HV=round((Nmat["N19","Nadm"]-NhuntedHV["H19","Hadm"])/(Nmat["N19","Total"]-sum(NhuntedHV["H19",])),3)
padmHV<-c(padm15HV,padm16HV,padm17HV,padm18HV,padm19HV)

#Number of adult females in post-harvest population
adf15HV=Nmat["N15","Nadf"]-NhuntedHV["H15","Hadf"]
adf16HV=Nmat["N16","Nadf"]-NhuntedHV["H16","Hadf"]
adf17HV=Nmat["N17","Nadf"]-NhuntedHV["H17","Hadf"]
adf18HV=Nmat["N18","Nadf"]-NhuntedHV["H18","Hadf"]
adf19HV=Nmat["N19","Nadf"]-NhuntedHV["H19","Hadf"]

adfHV<-c(adf15HV,adf16HV,adf17HV,adf18HV,adf19HV)

#N pop post-harvest 
pN15HV=Nmat["N15","Total"]-sum(NhuntedHV["H15",])
pN16HV=Nmat["N16","Total"]-sum(NhuntedHV["H16",])
pN17HV=Nmat["N17","Total"]-sum(NhuntedHV["H17",])
pN18HV=Nmat["N18","Total"]-sum(NhuntedHV["H18",])
pN19HV=Nmat["N19","Total"]-sum(NhuntedHV["H19",])
pN_HV<-c(pN15HV,pN16HV,pN17HV,pN18HV,pN19HV)



#tiff(filename = "Fig4.tif",width = 365, height = 645)
#pdf(file = "Fig4.pdf",width = 5.25, height = 8.70)

par(mfrow=c(2,1),mai = c(0.6, 1.2, 0.6, 0.6))

plot(seq(2017,2020,1),c(pN17,mr20nf2_tabV["PopSize",1:3]-mr20nf2_tabV["Nharvested",1:3]),type="l",pch=19,col="purple",
     las=1,lwd=2,ylab="Post-harvest population size",xlab="Year",
     ylim=c(100,700),xlim=c(2014.5,2020))
lines(seq(2017,2020,1),c(pN17,mr10nf2_tabV["PopSize",1:3]-mr10nf2_tabV["Nharvested",1:3]),type="l",pch=19,lwd=2,col="blue")
lines(seq(2017,2020,1),c(pN17,mr5nf2_tabV["PopSize",1:3]-mr5nf2_tabV["Nharvested",1:3]),type="l",pch=19,col="skyblue2",lwd=2)
lines(seq(2017,2020,1),c(pN17,mr3nf2_tabV["PopSize",1:3]-mr3nf2_tabV["Nharvested",1:3]),type="l",pch=19,col="green",lwd=2)
lines(c(2017,2018,2018.4,2019),pN_Nf2[3:6],type="l",lty=1,lwd=1,col="red")
points(c(2015,2016,2017,2018,2018.4,2019),pN_Nf2,pch=19,lwd=2,col="red")
lines(seq(2015,2017,1),pN_Nf2[1:3],type="l",pch=19,lwd=2,col="black")
lines(seq(2017,2020,1),c(pN17,S1nf2_tabV["PopSize",1:3]-S1nf2_tabV["Nharvested",1:3]),type="l",pch=19,lwd=2,col="black")

lines(seq(2015,2020,1),c(adfNf2[1:3],S1nf2_tabV["postHarvest_Nfemales",1:3]),type="l",lwd=3,col="darkgrey")


lines(c(2017,2018,2018.4,2019),adfNf2[3:6],type="l",lwd=1,col="darkgrey")
points(c(2015,2016,2017,2018,2018.4,2019),c(adfNf2),pch=19,lwd=2,col="darkgrey")

legend("topleft",
       c("emp. data","ordinary","m:f=1:3","m:f=1:5","m:f=1:10","m:f=1:20"),bty="n",lwd=c(1,2,2,2,2,2),lty=1,
       col=c("red","black","green","skyblue2","blue","purple"),pch=c(19,NA,NA,NA,NA,NA),cex=1)
legend("bottomleft",
       c("Adult females"),col=c("darkgrey"),pch=NA,lwd=c(3),lty=c(1),bty="n")
text(2016.5,650,"Nordfjella zone 2")

plot(seq(2017,2020,1),c(pN17HV,mr20hv_tabV["PopSize",1:3]-mr20hv_tabV["Nharvested",1:3]),type="l",pch=19,col="purple",
     las=1,lwd=2,ylab="Post-harvest population size",xlab="Year",
     ylim=c(2000,9000),xlim=c(2014.5,2020))
lines(seq(2017,2020,1),c(pN17HV,mr10hv_tabV["PopSize",1:3]-mr10hv_tabV["Nharvested",1:3]),type="l",pch=19,lwd=2,col="blue")
lines(seq(2017,2020,1),c(pN17HV,mr5hv_tabV["PopSize",1:3]-mr5hv_tabV["Nharvested",1:3]),type="l",pch=19,col="skyblue2",lwd=2)
lines(seq(2017,2020,1),c(pN17HV,mr3hv_tabV["PopSize",1:3]-mr3hv_tabV["Nharvested",1:3]),type="l",pch=19,col="green",lwd=2)
lines(c(2017,2018,2019),pN_HV[3:5],type="l",lty=1,lwd=1,col="red")
points(c(2015,2016,2017,2018,2019),pN_HV,pch=19,lwd=2,col="red")
lines(seq(2015,2017,1),pN_HV[1:3],type="l",pch=19,lwd=2,col="black")
lines(seq(2017,2020,1),c(pN17HV,S2hv_tabV["PopSize",1:3]-S2hv_tabV["Nharvested",1:3]),type="l",pch=19,lwd=2,col="black")

lines(seq(2015,2020,1),c(adfHV[1:3],S2hv_tabV["postHarvest_Nfemales",1:3]),type="l",lwd=3,col="darkgrey")


lines(c(2017,2018,2019),adfHV[3:5],type="l",lwd=1,col="darkgrey")
points(c(2015,2016,2017,2018,2019),c(adfHV),pch=19,lwd=2,col="darkgrey")

legend("left",
       c("emp. data","ordinary","m:f=1:3","m:f=1:5","m:f=1:10","m:f=1:20"),bty="n",lwd=c(1,2,2,2,2,2),lty=1,
       col=c("red","black","green","skyblue2","blue","purple"),pch=c(19,NA,NA,NA,NA,NA),cex=1)
legend("bottomleft",
       c("Adult females"),col=c("darkgrey"),pch=NA,lwd=c(3),lty=c(1),bty="n")
text(2016.5,8500,"Hardangervidda")


#dev.off()



##################################
#Figure 5AB (saved as tiff: 667 x 645, PDF: 6.95 inches * 6.72 inches)
##################################

par(mfrow=c(2,1),mai = c(0.6, 1, 0.1, 0.6))

plot(seq(2015,2028,1),c(0.5,mr5c00Nf2_RR124_probfreeS["pstarSA",])*100,type="l",col="blue",las=1,lwd=4,ylab="Probability of freedom (%)",xlab="Year",ylim=c(50,100))
abline(h=99,lty=2)

lines(seq(2015,2028,1),c(0.5,mr5c00Nf2_RR122_probfreeS["pstarSA",])*100,type="l",lwd=2,col="dodgerblue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5c00Nf2_RR126_probfreeS["pstarSA",])*100,type="l",lwd=4,col="darkmagenta",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5c00Nf2_RR111_probfreeS["pstarSA",])*100,type="l",lwd=1,col="skyblue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5c00Nf2_RR124_probfreeS["pstarSA",])*100,type="l",lwd=3,col="blue",lty=1)

lines(seq(2015,2028,1),c(0.5,S1Nf2_RR126_probfreeS["pstarSA",])*100,type="l",lwd=4,col="black",lty=1)
lines(seq(2015,2028,1),c(0.5,S1Nf2_RR124_probfreeS["pstarSA",])*100,type="l",lwd=3,col="azure4",lty=1)
lines(seq(2015,2028,1),c(0.5,S1Nf2_RR122_probfreeS["pstarSA",])*100,type="l",lwd=2,col="azure3",lty=1)
lines(seq(2015,2028,1),c(0.5,S1Nf2_RR111_probfreeS["pstarSA",])*100,type="l",lwd=1,col="black",lty=1)

legend("bottom",
       c("proactive","RR=1:1:1","RR=1:2:2","RR=1:2:4","RR=1:2:6"),bty="n",
       #col=c(NA,"dodgerblue","blue","blue3","darkblue"),
       col=c(NA,"skyblue","dodgerblue","blue","darkmagenta"),
       lty=c(NA,1,1,1,1),lwd=c(NA,1,2,3,4),pch=NA,cex=1)
legend("bottomright",
       c("ordinary","RR=1:1:1","RR=1:2:2","RR=1:2:4","RR=1:2:6"),bty="n",
       col=c(NA,"black","azure3","azure4","black"),
       lty=c(NA,1,1,1,1),lwd=c(NA,1,2,3,4),pch=NA,cex=1)


mr5Nf2_probfreeS05001=mr5c00nf2_probfreeS
mr5Nf2_probfreeS05=mr5c00Nf2_probfreeS05

plot(seq(2015,2028,1),c(0.5,mr5Nf2_probfreeS05001["pstarSA",])*100,las=1,type="l",pch=19,col="dodgerblue",lwd=2,ylab="Probability of freedom (%)",xlab="Year",ylim=c(50,100))
abline(h=99,lty=2)

lines(seq(2015,2028,1),c(0.5,mr5Nf2_probfreeS05001["pstar2",])*100,type="l",lwd=2,col="skyblue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Nf2_probfreeS05001["pstar4",])*100,type="l",lwd=2,col="blue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Nf2_probfreeS05001["pstar6",])*100,type="l",lwd=2,col="purple",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Nf2_probfreeS05001["pstar10",])*100,type="l",lwd=2,col="darkblue",lty=1)

legend("bottomright",
       c("p*=2","p*=2-4","p*=4","p*=6","p*=10"),bty="n",
       col=c("skyblue","dodgerblue","blue","purple","darkblue"),lty=1,lwd=c(2,2,2,2,2),pch=NA,cex=1)

lines(seq(2015,2028,1),c(0.5,mr5Nf2_probfreeS05["pstarSA",])*100,type="l",lwd=1,col="red",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Nf2_probfreeS05001["pstarSA",])*100,type="l",lwd=2,col="dodgerblue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Nf2_probfreeS05001["pstar2",])*100,type="l",lwd=2,col="skyblue",lty=1)


legend("bottom",
       c("p*=2-4, pIntro=5%"),bty="n",
       col=c("red"),lty=1,lwd=c(1),pch=c(NA),cex=1)



##################################
#Figure Supplementary 1 AB (saved as tiff: 667 x 645, PDF: 6.95 inches * 6.72 inches)
##################################
par(mfrow=c(2,1),mai = c(0.6, 1, 0.1, 0.6))

plot(seq(2015,2028,1),c(0.5,mr5Hv_RR124_probfreeS["pstar4",])*100,type="l",col="blue",las=1,lwd=4,ylab="Probability of freedom (%)",xlab="Year",ylim=c(50,100))
abline(h=99,lty=2)

lines(seq(2015,2028,1),c(0.5,mr5Hv_RR122_probfreeS["pstar4",])*100,type="l",lwd=2,col="dodgerblue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Hv_RR126_probfreeS["pstar4",])*100,type="l",lwd=4,col="darkmagenta",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Hv_RR111_probfreeS["pstar4",])*100,type="l",lwd=1,col="skyblue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Hv_RR124_probfreeS["pstar4",])*100,type="l",lwd=3,col="blue",lty=1)

lines(seq(2015,2028,1),c(0.5,Hv_RR126_probfreeS["pstar4",])*100,type="l",lwd=4,col="black",lty=1)
lines(seq(2015,2028,1),c(0.5,Hv_RR124_probfreeS["pstar4",])*100,type="l",lwd=3,col="azure4",lty=1)
lines(seq(2015,2028,1),c(0.5,Hv_RR122_probfreeS["pstar4",])*100,type="l",lwd=2,col="azure3",lty=1)
lines(seq(2015,2028,1),c(0.5,Hv_RR111_probfreeS["pstar4",])*100,type="l",lwd=1,col="black",lty=1)

legend("bottom",
       c("proactive","RR=1:1:1","RR=1:2:2","RR=1:2:4","RR=1:2:6"),bty="n",
       #col=c(NA,"dodgerblue","blue","blue3","darkblue"),
       col=c(NA,"skyblue","dodgerblue","blue","darkmagenta"),
       lty=c(NA,1,1,1,1),lwd=c(NA,1,2,3,4),pch=NA,cex=1)
legend("bottomright",
       c("ordinary","RR=1:1:1","RR=1:2:2","RR=1:2:4","RR=1:2:6"),bty="n",
       col=c(NA,"black","azure3","azure4","black"),
       lty=c(NA,1,1,1,1),lwd=c(NA,1,2,3,4),pch=NA,cex=1)


mr5Hv_probfreeS05001=mr5Hv_RR126_probfreeS #01001
mr5Hv_probfreeS05

plot(seq(2015,2028,1),c(0.5,mr5Hv_probfreeS05001["pstar4",])*100,las=1,type="l",pch=19,col="dodgerblue",lwd=2,ylab="Probability of freedom (%)",xlab="Year",ylim=c(50,100))
abline(h=99,lty=2)

lines(seq(2015,2028,1),c(0.5,mr5Hv_probfreeS05001["pstar2",])*100,type="l",lwd=2,col="skyblue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Hv_probfreeS05001["pstar6",])*100,type="l",lwd=2,col="purple",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Hv_probfreeS05001["pstar10",])*100,type="l",lwd=2,col="darkblue",lty=1)

legend("bottomright",
       c("p*=2","p*=4","p*=6","p*=10"),bty="n",
       col=c("skyblue","dodgerblue","purple","darkblue"),lty=1,lwd=c(2,2,2,2),pch=NA,cex=1)

lines(seq(2015,2028,1),c(0.5,mr5Hv_probfreeS05["pstar4",])*100,type="l",lwd=1,col="red",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Hv_probfreeS05001["pstar4",])*100,type="l",lwd=2,col="dodgerblue",lty=1)
lines(seq(2015,2028,1),c(0.5,mr5Hv_probfreeS05001["pstar2",])*100,type="l",lwd=2,col="skyblue",lty=1)


legend("bottom",
       c("p*=4, pIntro=5%"),bty="n",
       col=c("red"),lty=1,lwd=c(1),pch=c(NA),cex=1)



