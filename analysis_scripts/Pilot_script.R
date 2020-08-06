#packages

library(tidyverse)
library(readxl)
library(vegan)

pilot<-readxl::read_excel("Data/Pilot_data.xlsx", na=".")
pilot<-pilot[-56,]

pilot.sqrt<-as.data.frame(sapply(7:19, function(x) sqrt(pilot[,x])))
pilot.log<-as.data.frame(sapply(7:19, function(x) log(pilot[,x])))
  
pilot.log[is.na(pilot.log)]<-0
pilot.sqrt[is.na(pilot.sqrt)]<-0
log.matrix<-as.matrix(pilot.log)
sqrt.matrix<-as.matrix(pilot.sqrt)
  
 nmds<- metaMDS(sqrt.matrix, k=2, trymax=50)
 
 nmds.spp.fit <- envfit(nmds, sqrt.matrix, permutations = 999) #this gives us which compounds are significantly contributing to the separation
 
 #colvec <- c("tan1", "plum3", "lightskyblue3", "hotpink1")
 colvec <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")
 
 ordiplot(nmds, type = "n")
 #ordiellipse(nmds, groups = pilot$Sp, draw = "polygon", lty = 1, col = "grey", kind="se", conf=0.95, alpha=0.5, label=TRUE)
 orditorp(nmds, display="sites", labels=F, pch=21, col='black', bg =colvec[as.factor(pilot$Sp)], cex=2)
 plot(nmds.spp.fit, p.max=0.01, col="black",cex=1.2)
 legend(-2.2,1.45, legend=levels(as.factor(pilot$Sp)), col="black", pch=21, pt.bg=colvec, cex=1)
 legend(1.2,1.45, "stress = 0.139", bty = "n", cex = 1)
 
 
 