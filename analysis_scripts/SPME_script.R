library(reshape2)
setwd("Dropbox/R projects/Scent")

spme<-read.csv("SPME_data_no_blanks.csv",header=TRUE)

spme_wide<-dcast(spme, Sample~Name, sum, value.var="Area")
spme_ID<-read.csv("SPME_info.csv",header=TRUE)

spme_wide<-cbind(spme_wide, spme_ID[,2:4])

CB323<-spme_wide[which(spme_wide$Site=="CB323"),]
CB323.U.dose<-CB323[which(CB323$Type=="3.Fl"|CB323$Type=="6.Fl"),]
CB323.U.dose<-CB323.U.dose[order(CB323.U.dose$Type),]
CB323.U.dose.3<-CB323[which(CB323$Type=="3.Fl"),]
CB323.U.dose.6<-CB323[which(CB323$Type=="6.Fl"),]

ts <- sapply(2:48,function(x) (t.test((CB323.U.dose[1:3,x]), CB323.U.dose[4:6,x])))
Pvalue <- t(as.matrix(ts[3,]))

CB323<-spme_wide[which(spme_wide$Site=="CB323"),]
CB323.U.petals<-CB323[which(CB323$Type=="Petals"|CB323$Type=="Not Petals"),]
CB323.U.petals<-CB323.U.petals[order(CB323.U.petals$Type),]

ts <- sapply(2:48,function(x) (t.test((CB323.U.petals[1:3,x]), CB323.U.petals[4:6,x])))
Pvalue <- t(as.matrix(ts[3,]))


SC.U<-spme_wide[which(spme_wide$Site=="SC"&spme_wide$Species=="Ung"),]
SC.U.petals<-SC.U[which(SC.U$Type=="Petals"|SC.U$Type=="Not Petals"),]
SC.U.petals<-SC.U.petals[order(SC.U.petals$Type),]

ts <- sapply(2:48,function(x) (t.test((SC.U.petals[1:3,x]), SC.U.petals[4:6,x])))
Pvalue <- t(as.matrix(ts[3,]))

SC.C<-spme_wide[which(spme_wide$Site=="SC"&spme_wide$Species=="Cyl"),]
SC.C.petals<-SC.C[which(SC.C$Type=="Petals"|SC.C$Type=="Not Petals"),]
SC.C.petals<-SC.C.petals[order(SC.C.petals$Type),]

ts <- sapply(2:48,function(x) (t.test((SC.C.petals[1:3,x]), SC.C.petals[4:6,x],paired=TRUE)))
Pvalue <- t(as.matrix(ts[3,]))

MHG.C<-spme_wide[which(spme_wide$Site=="MHG"&spme_wide$Species=="Cyl"),]
MHG.C.petals<-MHG.C[which(MHG.C$Type=="Petals"|MHG.C$Type=="Not Petals"),]
MHG.C.petals<-MHG.C.petals[order(MHG.C.petals$Type),]

ts <- sapply(2:48,function(x) (t.test((MHG.C.petals[1:3,x]), MHG.C.petals[4:6,x],paired=TRUE)))
Pvalue <- t(as.matrix(ts[3,]))
