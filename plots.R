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