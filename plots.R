#packages:
library(ggplot2)
library(plotrix)

#load data, generate subsets etc.:
mass.data<-read.csv("data/all_data_standardized_mass.csv", header=TRUE)
mass.data$Site.Type<-factor(mass.data$Site.Type, levels=c("One", "Two", "Four"))

mass.data$total<-rowSums(cbind(mass.data[,c(12:65)]))

mass.data$mono<-(mass.data$er.2.6.dimethyl.1.3.5.7.octatetraene..trans.+mass.data$er.6.methyl.5.hepten.2.one +mass.data$er.alpha.pinene +mass.data$er.alpha.terpineol + mass.data$er.alpha.terpinene+ mass.data$er.alpha.terpinolene+ mass.data$er.beta.myrcene + mass.data$er.beta.phellandrene + mass.data$er.beta.pinene+mass.data$er.borneol + mass.data$er.cis.beta.ocimene + mass.data$er.gamma.terpinene + mass.data$er.geraniol + mass.data$er.limonene + mass.data$er.linalool + mass.data$er.myroxide + mass.data$er.para.cymene + mass.data$er.sabinene+ mass.data$er.sabinene.hydrate + mass.data$er.terpinen.4.ol + mass.data$er.trans.beta.ocimene + mass.data$er.verbenone)
mass.data$glv<-(mass.data$er.1.hexanol+mass.data$er.cis.3.hexen.1.ol+mass.data$er.cis.3.hexenyl.acetate+mass.data$er.cis.jasmone+ mass.data$er.trans.2.hexen.1.ol )
mass.data$arom<-(mass.data$er.2.amino.phenyl.ethenone+ mass.data$er.2.phenyl.ethanol+ mass.data$er.benzyl.acetate+ mass.data$er.benzyl.alcohol+ mass.data$er.cinnamic.alcohol+ mass.data$er.methyl.nicotinate+ mass.data$er.methyl.salicylate+ mass.data$er.trans.cinnamic.aldehyde+ mass.data$er.veratrole)
mass.data$sesqui<-(mass.data$er.alloaromadendrene+ mass.data$er.alpha.bergamotene + mass.data$er.alpha.humulene + mass.data$er.beta.cadinene + mass.data$er.beta.farnesene + mass.data$er.beta.longipinene + mass.data$er.caryophyllene.oxide + mass.data$er.DMNT + mass.data$er.gamma.cadinene + mass.data$er.germacrene.D + mass.data$er.intermedeol + mass.data$er.patchoulane + mass.data$er.trans.beta.caryophyllene + mass.data$er.trans.DMNT + mass.data$er.trans.trans.alpha.farnesene + mass.data$er.unknown.C15H24+mass.data$er.unknown.C15H24.2 + mass.data$er.unknown.C15H24.3)


mass.data.compounds<-mass.data[,12:65]
mass.data$Code<-as.factor(paste(mass.data$Species,mass.data$Site.Type))

##Generate capscale:

scents.cap <-capscale(mass.data.compounds ~ Code, mass.data, dist="bray") 

##Visualizing the capscale:

pl<-plot(scents.cap, type="n", correlation=TRUE, xlim=c(-0.5,0.5))
with(mass.data, points(pl, "site",pch=c(21,22)[mass.data$Species], bg= color[mass.data$Code],cex=1.5 ))
text(pl, "sp", scaling=1, arrow=TRUE, length=0.05, col=4, cex=0.6, xpd=TRUE)


##plot means of the four compound classes by species and community types
means<-as.data.frame(sapply(66:70,function(x) tapply(mass.data[,x], list(mass.data$Species, mass.data$Site.Type), mean)))

colnames(means)<-c("mono", "glv", "arom", "sesqui")
Species<-c("Cyl", "Ung","Cyl", "Ung","Cyl", "Ung")
Site.Type<-c("One","One","Two","Two","Four","Four")
means$Species<-as.factor(Species)
means$Site.Type<-as.factor(Site.Type)
means$Site.Type<-factor(means$Site.Type, levels=c("One", "Two", "Four"))

ses<-as.data.frame(sapply(66:70,function(x) tapply(mass.data[,x], list(mass.data$Species, mass.data$Site.Type), std.error)))
colnames(ses)<-c("mono_se", "glv_se", "arom_se", "sesqui_se")

dat<-cbind(means, ses)

#now can make this plot for each of the four compound classes in the object above:

ggplot(aes(x=Site.Type, y=arom, color=Species), data=dat)+geom_point(size=10, position=position_dodge(width = 0.50))+geom_errorbar(aes(ymin= arom-arom_se, ymax= arom+arom_se), width=0, size=1, position=position_dodge(width = 0.50))+theme_classic(base_size=20)+ylab("")+xlab("Community Type")+ theme(legend.position = "none")


##plot means of the nine compounds that showed CD patterns by species and community types

means<-as.data.frame(sapply(c(20,31,36,13,41,54,47,64,59), function(x) tapply(mass.data[,x], list(mass.data$Species, mass.data$Site.Type), mean)))

colnames(means)<-c("a.pinene", "b.pinene", "cis.hex", "amino.phenyl","gamma.terp","sab.hyd","meth.nic","ver", "cinn" )
Species<-c("Cyl", "Ung","Cyl", "Ung","Cyl", "Ung")
Site.Type<-c("One","One","Two","Two","Four","Four")
means$Species<-as.factor(Species)
means$Site.Type<-as.factor(Site.Type)
means$Site.Type<-factor(means$Site.Type, levels=c("One", "Two", "Four"))

ses<-as.data.frame(sapply(c(20,31,36,13,41,54,47,64,59),function(x) tapply(mass.data[,x], list(mass.data$Species, mass.data$Site.Type), std.error)))
colnames(ses)<-c("a.pinene_se", "b.pinene_se", "cis.hex_se", "amino.phenyl_se","gamma.terp_se","sab.hyd_se","meth.nic_se","ver_se", "cinn_se" )

dat<-cbind(means, ses)


##can now run this plot code for each of the nine compounds in the above object:
ggplot(aes(x=Site.Type, y=ver, color=Species), data=dat)+geom_point(size=10, position=position_dodge(width = 0.50))+geom_errorbar(aes(ymin= ver-ver_se, ymax= ver+ ver_se), width=0, size=1, position=position_dodge(width = 0.50))+theme_classic(base_size=20)+ylab("")+xlab("Community type")+ theme(legend.position = "none")