
wd = "D:/ISmalia_DSSAT"

Obs_data = read.csv(paste(wd,"/OBS_Calib-Data_.csv",sep=""))

#-------------------------SC10 cultivar------------------------------

OBStotalDW <- round(Obs_data$Value[Obs_data$Vari=="totaldryweight" & Obs_data$Nit=="190" & Obs_data$DAP %in% c("60","75") & Obs_data$Trt %in% c(1:20)],digits=0)
OBSGY <- round(Obs_data$Value[Obs_data$Vari=="GrainYield" & Obs_data$Nit=="190" & Obs_data$DAP=="110" & Obs_data$Trt %in% c(1:20)],digits=0)


par <- read.csv(paste(wd,"/Corn_calibration.csv",sep=""))

X <- par$Calib_SC10

#-----edit the CUL&ECO files----------

out_string <- paste0("!MAIZE CULTIVAR COEFFICIENTS: MZCER047 MODEL\n",
                     "!AHMED ATTIA (MAR-2019)\n",
                     "@VAR#  VRNAME.......... EXPNO   ECO#    P1    P2    P5    G2    G3 PHINT")

CalibP <- c(formatC(X[1],format="f",digits=1),formatC(X[2],format="f",digits=3),formatC(X[3],format="f",digits=1),
            formatC(X[4],format="f",digits=1), formatC(X[5],format="f",digits=2),
            formatC(X[6],format="f",digits=2))
CalibP[5] <- paste(" ",CalibP[5],sep="")

cat(out_string,
    "990001 LONG SEASON          . IB0001",CalibP,file="C:/DSSAT47/Genotype/MZCER047.CUL",fill=T,append = F)

out_string2 <- paste0("*MAIZE ECOTYPE COEFFICIENTS: MZCER047 MODEL\n",
                      "@ECO#  ECONAME.........  TBASE  TOPT ROPT   P20  DJTI  GDDE  DSGFT  RUE   KCAN  TSEN  CDAY")


CalibP2 <- c(formatC(X[7],format="f",digits=1),formatC(X[8],format="f",digits=1),formatC(X[9],format="f",digits=1),
             formatC(X[10],format="f",digits=1),formatC(X[11],format="f",digits=1),
             formatC(X[12],format="f",digits=1),formatC(X[13],format="f",digits=1),formatC(X[14],format="f",digits=1),
             formatC(X[15],format="f",digits=2))
CalibP2[2] <- paste(" ",CalibP2[2],sep="")
CalibP2[3] <- paste("",CalibP2[3],sep="")
CalibP2[4] <- paste(" ",CalibP2[4],sep="")
CalibP2[5] <- paste("  ",CalibP2[5],sep="")
CalibP2[6] <- paste("  ",CalibP2[6],sep="")
CalibP2[7] <- paste(" ",CalibP2[7],sep="")
CalibP2[8] <- paste(" ",CalibP2[8],sep="")
CalibP2[9] <- paste("  ",CalibP2[9],sep="")

cat(out_string2,
    "IB0001 GENERIC MIDWEST1   ",CalibP2,file="C:/DSSAT47/Genotype/MZCER047.ECO",fill=T,append = F)


setwd(paste("C:/DSSAT47/Maize",sep = ""))

#--- write paramters used on the screen
message("")
message("Running DSSAT-MaizeCERES...")

#--- Call DSSAT047.exe and run X files list within DSSBatch.v47
system("C:/DSSAT47/DSCSM047.EXE MZCER047 B DSSBatch.v47",show.output.on.console = F)


plantgro <- read.dssat("C:/DSSAT47/Maize/PlantGro.OUT")


SIMtotalDW60 <- 0 
SIMtotalDW75 <- 0 
SIMtotalDW <- 0
SIMGY <- 0 

for(i in 1:length(plantgro)){
  
  data=as.data.frame(plantgro[[i]])
  
  SIMtotalDW60[i] <- data$CWAD[data$DAP==60]
  SIMtotalDW75[i] <- data$CWAD[data$DAP==75]
  SIMGY[i] <- tail(data$GWAD,n=1)
  
}


simtotalDW=c(SIMtotalDW60,SIMtotalDW75)
simGY=c(SIMGY)

#-------------------plot--------------------

par(mfrow=c(2,2),mar=c(4,4,3,2)+0.1,mgp=c(2.5,0.7,0))

plot(OBSGY/1000,simGY/1000,xlim=c(0,5),ylim=c(0,5),
     ylab=expression(paste("Simulated grain yield (Mg ",ha^-1,")")),xlab=expression(paste("Observed grain yield (Mg ",ha^-1,")")),
     main="SC10",font.main=2,pch=21,yaxs="i",bg="darkgray",xaxs="i",cex.axis=0.9,cex=1.4,las=1)

SimregGY <- simGY/1000
ObsregGY <- OBSGY/1000
reg1<- lm(SimregGY~ObsregGY)
abline(reg1,pch=4,col=2,lwd=2, lty=2)
abline(0:1)
modeleval_GY <- modeval(simGY,OBSGY)
text(1,4.8,label=bquote("R"^2~":" ~ .(round(modeleval_GY$R2[[1]],digits=2))),cex=0.8) 
text(1.4,4.4,label=noquote(paste0("RMSE : ",round(modeleval_GY$RMSE[[1]]/1000,digits=3)," Mg/ha")),cex=0.8)
text(3.5,1,label=noquote(paste0("RRMSE : ",round(modeleval_GY$RRMSE[[1]],digits=2),"%")),cex=0.8) 
text(3.3,0.6,label=noquote(paste0("NSE : ",round(modeleval_GY$EF[[1]],digits=2))),cex=0.8)

plot(OBStotalDW/1000,simtotalDW/1000,xlim=c(0,8),ylim=c(0,8),
     ylab=expression(paste("Simulated total dry matter (Mg ",ha^-1,")")),xlab=expression(paste("Observed total dry matter (Mg ",ha^-1,")")),
     main="SC10",font.main=2,pch=21,yaxs="i",bg="darkgray",xaxs="i",cex.axis=0.9,cex=1.5,las=1)

SimregtotalDW <- simtotalDW/1000
ObsregtotalDW <- OBStotalDW/1000
reg1<- lm(SimregtotalDW~ObsregtotalDW)
abline(reg1,pch=4,col=2,lwd=2, lty=2)
abline(0:1)
modeleval_totDM <- modeval(SimregtotalDW,ObsregtotalDW)
text(2,7.6,label=bquote("R"^2~":" ~ .(round(modeleval_totDM$R2[[1]],digits=2))),cex=0.8) 
text(2.7,6.8,label=noquote(paste0("RMSE : ",round(modeleval_totDM$RMSE[[1]],digits=3)," Mg/ha")),cex=0.8)
text(6.0,1.5,label=noquote(paste0("RRMSE : ",round(modeleval_totDM$RRMSE[[1]],digits=2),"%")),cex=0.8) 
text(6.0,0.8,label=noquote(paste0("NSE : ",round(modeleval_totDM$EF[[1]],digits=2))),cex=0.8)

Trt105 <- as.data.frame(plantgro[[1]])
Trt205 <- as.data.frame(plantgro[[2]])
Trt305 <- as.data.frame(plantgro[[3]])
Trt405 <- as.data.frame(plantgro[[4]])
Trt505 <- as.data.frame(plantgro[[5]])
Trt106 <- as.data.frame(plantgro[[6]])
Trt206 <- as.data.frame(plantgro[[7]])
Trt306 <- as.data.frame(plantgro[[8]])
Trt406 <- as.data.frame(plantgro[[9]])
Trt506 <- as.data.frame(plantgro[[10]])

plot(Trt105$DAP[1:95],Trt105$CWAD[1:95]/1000,type="l",ylim=c(0,9),xlim=c(0,110), ylab=expression(paste("Total dry matter (Mg  ",ha^-1,")")),
     xlab="Days after planting",xaxs="i",cex.axis=0.9,cex=1.5,las=1,col=1,main="SC10 2005")
points(60,Obs_data$Value[Obs_data$Trt=="1" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=1,cex=1.5,col=1)
points(75,Obs_data$Value[Obs_data$Trt=="1" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=1,cex=1.5,col=1)
lines(Trt205$DAP[1:95],Trt205$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=2)
points(60,Obs_data$Value[Obs_data$Trt=="2" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=2,cex=1.5,col=2)
points(75,Obs_data$Value[Obs_data$Trt=="2" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=2,cex=1.5,col=2)
lines(Trt305$DAP[1:95],Trt305$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=3)
points(60,Obs_data$Value[Obs_data$Trt=="3" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=3,cex=1.5,col=3)
points(75,Obs_data$Value[Obs_data$Trt=="3" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=3,cex=1.5,col=3)
lines(Trt405$DAP[1:95],Trt405$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=4)
points(60,Obs_data$Value[Obs_data$Trt=="4" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=4,cex=1.5,col=4)
points(75,Obs_data$Value[Obs_data$Trt=="4" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=4,cex=1.5,col=4)
lines(Trt505$DAP[1:95],Trt505$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=5)
points(60,Obs_data$Value[Obs_data$Trt=="5" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=5,cex=1.5,col=5)
points(75,Obs_data$Value[Obs_data$Trt=="5" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=5,cex=1.5,col=5)
legend("topleft",c("F2","F3","F4","F5","F6"),pch=19,lty=c(1,2,3),col=c(1:5),bty="n",cex=0.9)

plot(Trt106$DAP[1:95],Trt106$CWAD[1:95]/1000,type="l",ylim=c(0,9),xlim=c(0,110), ylab=expression(paste("Total dry matter (Mg  ",ha^-1,")")),
     xlab="Days after planting",xaxs="i",cex.axis=0.9,cex=1.5,las=1,col=1,main="SC10 2006")
points(60,Obs_data$Value[Obs_data$Trt=="11" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=1,cex=1.5,col=1)
points(75,Obs_data$Value[Obs_data$Trt=="11" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=1,cex=1.5,col=1)
lines(Trt206$DAP[1:95],Trt206$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=2)
points(60,Obs_data$Value[Obs_data$Trt=="12" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=2,cex=1.5,col=2)
points(75,Obs_data$Value[Obs_data$Trt=="12" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=2,cex=1.5,col=2)
lines(Trt306$DAP[1:95],Trt306$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=3)
points(60,Obs_data$Value[Obs_data$Trt=="13" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=3,cex=1.5,col=3)
points(75,Obs_data$Value[Obs_data$Trt=="13" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=3,cex=1.5,col=3)
lines(Trt406$DAP[1:95],Trt406$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=4)
points(60,Obs_data$Value[Obs_data$Trt=="14" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=4,cex=1.5,col=4)
points(75,Obs_data$Value[Obs_data$Trt=="14" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=4,cex=1.5,col=4)
lines(Trt506$DAP[1:95],Trt506$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=5)
points(60,Obs_data$Value[Obs_data$Trt=="15" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=5,cex=1.5,col=5)
points(75,Obs_data$Value[Obs_data$Trt=="15" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=5,cex=1.5,col=5)
legend("topleft",c("F2","F3","F4","F5","F6"),pch=19,lty=c(1,2,3),col=c(1:5),bty="n",cex=0.9)

#------------------------------TC310-------------------

OBStotalDW <- round(Obs_data$Value[Obs_data$Vari=="totaldryweight" & Obs_data$Nit=="190" & Obs_data$DAP %in% c("60","75") & Obs_data$Trt %in% c(21:40)],digits=0)
OBSGY <- round(Obs_data$Value[Obs_data$Vari=="GrainYield" & Obs_data$Nit=="190" & Obs_data$DAP=="110" & Obs_data$Trt %in% c(21:40)],digits=0)

X <- par$Calib_TW310

#-----edit the CUL&ECO files----------

out_string <- paste0("!MAIZE CULTIVAR COEFFICIENTS: MZCER047 MODEL\n",
                     "!AHMED ATTIA (MAR-2019)\n",
                     "@VAR#  VRNAME.......... EXPNO   ECO#    P1    P2    P5    G2    G3 PHINT")

CalibP <- c(formatC(X[1],format="f",digits=1),formatC(X[2],format="f",digits=3),formatC(X[3],format="f",digits=1),
            formatC(X[4],format="f",digits=1), formatC(X[5],format="f",digits=2),
            formatC(X[6],format="f",digits=2))
CalibP[5] <- paste(" ",CalibP[5],sep="")

cat(out_string,
    "990001 LONG SEASON          . IB0001",CalibP,file="C:/DSSAT47/Genotype/MZCER047.CUL",fill=T,append = F)

out_string2 <- paste0("*MAIZE ECOTYPE COEFFICIENTS: MZCER047 MODEL\n",
                      "@ECO#  ECONAME.........  TBASE  TOPT ROPT   P20  DJTI  GDDE  DSGFT  RUE   KCAN  TSEN  CDAY")


CalibP2 <- c(formatC(X[7],format="f",digits=1),formatC(X[8],format="f",digits=1),formatC(X[9],format="f",digits=1),
             formatC(X[10],format="f",digits=1),formatC(X[11],format="f",digits=1),
             formatC(X[12],format="f",digits=1),formatC(X[13],format="f",digits=1),formatC(X[14],format="f",digits=1),
             formatC(X[15],format="f",digits=2))
CalibP2[2] <- paste(" ",CalibP2[2],sep="")
CalibP2[3] <- paste("",CalibP2[3],sep="")
CalibP2[4] <- paste(" ",CalibP2[4],sep="")
CalibP2[5] <- paste("  ",CalibP2[5],sep="")
CalibP2[6] <- paste("  ",CalibP2[6],sep="")
CalibP2[7] <- paste(" ",CalibP2[7],sep="")
CalibP2[8] <- paste(" ",CalibP2[8],sep="")
CalibP2[9] <- paste("  ",CalibP2[9],sep="")

cat(out_string2,
    "IB0001 GENERIC MIDWEST1   ",CalibP2,file="C:/DSSAT47/Genotype/MZCER047.ECO",fill=T,append = F)


setwd(paste("C:/DSSAT47/Maize",sep = ""))

#--- write paramters used on the screen
message("")
message("Running DSSAT-MaizeCERES...")

#--- Call DSSAT047.exe and run X files list within DSSBatch.v47
system("C:/DSSAT47/DSCSM047.EXE MZCER047 B DSSBatch.v47",show.output.on.console = F)


plantgro <- read.dssat("C:/DSSAT47/Maize/PlantGro.OUT")


SIMtotalDW60 <- 0 
SIMtotalDW75 <- 0 
SIMtotalDW <- 0
SIMGY <- 0 

for(i in 1:length(plantgro)){
  
  data=as.data.frame(plantgro[[i]])
  
  SIMtotalDW60[i] <- data$CWAD[data$DAP==60]
  SIMtotalDW75[i] <- data$CWAD[data$DAP==75]
  SIMGY[i] <- tail(data$GWAD,n=1)
  
}


simtotalDW=c(SIMtotalDW60,SIMtotalDW75)
simGY=c(SIMGY)

#-------------------plot--------------------


plot(OBSGY/1000,simGY/1000,xlim=c(0,5),ylim=c(0,5),
     ylab=expression(paste("Simulated grain yield (Mg ",ha^-1,")")),xlab=expression(paste("Observed grain yield (Mg ",ha^-1,")")),
     main="TC310",font.main=2,pch=21,yaxs="i",bg="darkgray",xaxs="i",cex.axis=0.9,cex=1.4,las=1)

SimregGY <- simGY/1000
ObsregGY <- OBSGY/1000
reg1<- lm(SimregGY~ObsregGY)
abline(reg1,pch=4,col=2,lwd=2, lty=2)
abline(0:1)
modeleval_GY <- modeval(simGY,OBSGY)
text(1,4.8,label=bquote("R"^2~":" ~ .(round(modeleval_GY$R2[[1]],digits=2))),cex=0.8) 
text(1.4,4.4,label=noquote(paste0("RMSE : ",round(modeleval_GY$RMSE[[1]]/1000,digits=3)," Mg/ha")),cex=0.8)
text(3.5,1,label=noquote(paste0("RRMSE : ",round(modeleval_GY$RRMSE[[1]],digits=2),"%")),cex=0.8) 
text(3.3,0.6,label=noquote(paste0("NSE : ",round(modeleval_GY$EF[[1]],digits=2))),cex=0.8)

plot(OBStotalDW/1000,simtotalDW/1000,xlim=c(0,8),ylim=c(0,8),
     ylab=expression(paste("Simulated total dry matter (Mg ",ha^-1,")")),xlab=expression(paste("Observed total dry matter (Mg ",ha^-1,")")),
     main="TC310",font.main=2,pch=21,yaxs="i",bg="darkgray",xaxs="i",cex.axis=0.9,cex=1.5,las=1)

SimregtotalDW <- simtotalDW/1000
ObsregtotalDW <- OBStotalDW/1000
reg1<- lm(SimregtotalDW~ObsregtotalDW)
abline(reg1,pch=4,col=2,lwd=2, lty=2)
abline(0:1)
modeleval_totDM <- modeval(SimregtotalDW,ObsregtotalDW)
text(2,7.6,label=bquote("R"^2~":" ~ .(round(modeleval_totDM$R2[[1]],digits=2))),cex=0.8) 
text(2.7,6.8,label=noquote(paste0("RMSE : ",round(modeleval_totDM$RMSE[[1]],digits=3)," Mg/ha")),cex=0.8)
text(6.0,1.5,label=noquote(paste0("RRMSE : ",round(modeleval_totDM$RRMSE[[1]],digits=2),"%")),cex=0.8) 
text(6.0,0.8,label=noquote(paste0("NSE : ",round(modeleval_totDM$EF[[1]],digits=2))),cex=0.8)

Trt105 <- as.data.frame(plantgro[[1]])
Trt205 <- as.data.frame(plantgro[[2]])
Trt305 <- as.data.frame(plantgro[[3]])
Trt405 <- as.data.frame(plantgro[[4]])
Trt505 <- as.data.frame(plantgro[[5]])
Trt106 <- as.data.frame(plantgro[[6]])
Trt206 <- as.data.frame(plantgro[[7]])
Trt306 <- as.data.frame(plantgro[[8]])
Trt406 <- as.data.frame(plantgro[[9]])
Trt506 <- as.data.frame(plantgro[[10]])

plot(Trt105$DAP[1:95],Trt105$CWAD[1:95]/1000,type="l",ylim=c(0,9),xlim=c(0,110), ylab=expression(paste("Total dry matter (Mg  ",ha^-1,")")),
     xlab="Days after planting",xaxs="i",cex.axis=0.9,cex=1.5,las=1,col=1,main="TC310 2005")
points(60,Obs_data$Value[Obs_data$Trt=="21" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=1,cex=1.5,col=1)
points(75,Obs_data$Value[Obs_data$Trt=="21" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=1,cex=1.5,col=1)
lines(Trt205$DAP[1:95],Trt205$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=2)
points(60,Obs_data$Value[Obs_data$Trt=="22" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=2,cex=1.5,col=2)
points(75,Obs_data$Value[Obs_data$Trt=="22" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=2,cex=1.5,col=2)
lines(Trt305$DAP[1:95],Trt305$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=3)
points(60,Obs_data$Value[Obs_data$Trt=="23" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=3,cex=1.5,col=3)
points(75,Obs_data$Value[Obs_data$Trt=="23" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=3,cex=1.5,col=3)
lines(Trt405$DAP[1:95],Trt405$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=4)
points(60,Obs_data$Value[Obs_data$Trt=="24" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=4,cex=1.5,col=4)
points(75,Obs_data$Value[Obs_data$Trt=="24" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=4,cex=1.5,col=4)
lines(Trt505$DAP[1:95],Trt505$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=5)
points(60,Obs_data$Value[Obs_data$Trt=="25" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=5,cex=1.5,col=5)
points(75,Obs_data$Value[Obs_data$Trt=="25" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=5,cex=1.5,col=5)
legend("topleft",c("F2","F3","F4","F5","F6"),pch=19,lty=c(1,2,3),col=c(1:5),bty="n",cex=0.9)

plot(Trt106$DAP[1:95],Trt106$CWAD[1:95]/1000,type="l",ylim=c(0,9),xlim=c(0,110), ylab=expression(paste("Total dry matter (Mg  ",ha^-1,")")),
     xlab="Days after planting",xaxs="i",cex.axis=0.9,cex=1.5,las=1,col=1,main="TC310 2006")
points(60,Obs_data$Value[Obs_data$Trt=="31" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=1,cex=1.5,col=1)
points(75,Obs_data$Value[Obs_data$Trt=="31" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=1,cex=1.5,col=1)
lines(Trt206$DAP[1:95],Trt206$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=2)
points(60,Obs_data$Value[Obs_data$Trt=="32" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=2,cex=1.5,col=2)
points(75,Obs_data$Value[Obs_data$Trt=="32" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=2,cex=1.5,col=2)
lines(Trt306$DAP[1:95],Trt306$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=3)
points(60,Obs_data$Value[Obs_data$Trt=="33" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=3,cex=1.5,col=3)
points(75,Obs_data$Value[Obs_data$Trt=="33" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=3,cex=1.5,col=3)
lines(Trt406$DAP[1:95],Trt406$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=4)
points(60,Obs_data$Value[Obs_data$Trt=="34" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=4,cex=1.5,col=4)
points(75,Obs_data$Value[Obs_data$Trt=="34" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=4,cex=1.5,col=4)
lines(Trt506$DAP[1:95],Trt506$CWAD[1:95]/1000,cex.axis=0.9,cex=1.5,las=1,col=5)
points(60,Obs_data$Value[Obs_data$Trt=="35" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="60"]/1000,pch=21,bg=5,cex=1.5,col=5)
points(75,Obs_data$Value[Obs_data$Trt=="35" & Obs_data$Vari=="totaldryweight" & Obs_data$DAP=="75"]/1000,pch=21,bg=5,cex=1.5,col=5)
legend("topleft",c("F2","F3","F4","F5","F6"),pch=19,lty=c(1,2,3),col=c(1:5),bty="n",cex=0.9)


