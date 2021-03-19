###Outline for Floral Scent Analysis

#basics:

library(vegan)
library(multcomp)
library(lme4)
library(lmerTest)
library(emmeans)

#the entire, mass-standardized, data file:

mass.data<-read.csv("data/all_data_standardized_mass.csv", header=TRUE)
mass.data$Site.Type<-factor(mass.data$Site.Type, levels=c("One", "Two", "Four"))

#adding in total scent:

mass.data$total<-rowSums(cbind(mass.data[,c(12:65)]))

#adding in type of compounds:

mass.data$mono<-(mass.data$er.2.6.dimethyl.1.3.5.7.octatetraene..trans.+mass.data$er.6.methyl.5.hepten.2.one +mass.data$er.alpha.pinene +mass.data$er.alpha.terpineol + mass.data$er.alpha.terpinene+ mass.data$er.alpha.terpinolene+ mass.data$er.beta.myrcene + mass.data$er.beta.phellandrene + mass.data$er.beta.pinene+mass.data$er.borneol + mass.data$er.cis.beta.ocimene + mass.data$er.gamma.terpinene + mass.data$er.geraniol + mass.data$er.limonene + mass.data$er.linalool + mass.data$er.myroxide + mass.data$er.para.cymene + mass.data$er.sabinene+ mass.data$er.sabinene.hydrate + mass.data$er.terpinen.4.ol + mass.data$er.trans.beta.ocimene + mass.data$er.verbenone)
mass.data$glv<-(mass.data$er.1.hexanol+mass.data$er.cis.3.hexen.1.ol+mass.data$er.cis.3.hexenyl.acetate+mass.data$er.cis.jasmone+ mass.data$er.trans.2.hexen.1.ol )
mass.data$arom<-(mass.data$er.2.amino.phenyl.ethenone+ mass.data$er.2.phenyl.ethanol+ mass.data$er.benzyl.acetate+ mass.data$er.benzyl.alcohol+ mass.data$er.cinnamic.alcohol+ mass.data$er.methyl.nicotinate+ mass.data$er.methyl.salicylate+ mass.data$er.trans.cinnamic.aldehyde+ mass.data$er.veratrole)
mass.data$sesqui<-(mass.data$er.alloaromadendrene+ mass.data$er.alpha.bergamotene + mass.data$er.alpha.humulene + mass.data$er.beta.cadinene + mass.data$er.beta.farnesene + mass.data$er.beta.longipinene + mass.data$er.caryophyllene.oxide + mass.data$er.DMNT + mass.data$er.gamma.cadinene + mass.data$er.germacrene.D + mass.data$er.intermedeol + mass.data$er.patchoulane + mass.data$er.trans.beta.caryophyllene + mass.data$er.trans.DMNT + mass.data$er.trans.trans.alpha.farnesene + mass.data$er.unknown.C15H24+mass.data$er.unknown.C15H24.2 + mass.data$er.unknown.C15H24.3)

#compounds only:

mass.data.compounds<-mass.data[,12:65]


#A. Multivariate analyses:
#first, adonis:
scents.permanova<-adonis(mass.data.compounds~Site.Type*Species, data=mass.data, strata=mass.data$Site.Type, permutations=999, method="bray")

scents.permanova

#second, capscale:


mass.data$Code<-as.factor(paste(mass.data$Species,mass.data$Site.Type))

scents.cap <-capscale(mass.data.compounds ~ Code, mass.data, dist="bray") 
scents.cap
summary(scents.cap)

#adding the capscale scores to the mass.data.compounds object, to be able to assess how correlated each compound is with the cap axes:

scores<-as.data.frame(scores(scents.cap, display="sites"))


mass.data.compounds.scores<-cbind(mass.data.compounds, scores)

cors <- as.data.frame(sapply(1:54,function(x) cor.test(mass.data.compounds.scores[,x],mass.data.compounds.scores$CAP1)))

cors2 <- as.data.frame(sapply(1:54,function(x) cor.test(mass.data.compounds.scores[,x],mass.data.compounds.scores$CAP2)))
pvals<-sapply(1:54,function(x) cors[,x]$p.value)
pvals2<-sapply(1:54,function(x) cors2[,x]$p.value)
pvals.a<-p.adjust(pvals, method="BH")
pvals.a2<-p.adjust(pvals2, method="BH")


co1<-sapply(1:54,function(x) cors[,x]$estimate)
co2<-sapply(1:54,function(x) cors2[,x]$estimate)


compounds<-colnames(mass.data.compounds)

c.table<-as.data.frame(cbind(compounds,cor1=co1,cap1=pvals.a,cor2=co2, cap2=pvals.a2))
c.table$cap1<-as.numeric(as.character(c.table$cap1))
c.table$cap2<-as.numeric(as.character(c.table$cap2))
c.table$sig1<-rep("No", times=dim(c.table)[1])
c.table$sig2<-rep("No", times=dim(c.table)[1])
c.table$sig1<-ifelse(c.table$cap1 < 0.01, "Yes", c.table$sig1)
c.table$sig2<-ifelse(c.table$cap2 < 0.01, "Yes", c.table$sig2)


## B. Univariate analyses
#first, total scent:


total<-lmer(sqrt(total)~Site.Type*Species+ (1|Site), data=mass.data)

hist(resid(total))
plot(predict(total),resid(total)) 
summary(total)
anova(total)


emmeans(total, "Species",type="response" )
contrast<-emmeans(total, pairwise~Species|Site.Type, type="response")

#second, compound classes:

ses<-lmer(log(glv+0.0001)~Site.Type*Species+ (1|Site), data=mass.data)
hist(resid(ses))
plot(predict(ses),resid(ses)) 

anova(ses)

contrast<-emmeans(ses, pairwise~Site.Type|Species, type="response")
contrast



contrast.matrix3 <- rbind(
  `C-U: one vs two` = c(0, 0, 0, 0, 1, 0),
  `C-U: one vs four` = c(0, 0, 0, 0, 0, 1),
  `C-U: two vs four` = c(0, 0, 0, 0, -1, 1))

comps3 <- glht(ses, contrast.matrix3)
summary(comps3)

ses<-lmer(log(sesqui+0.0001)~Site.Type*Species+ (1|Site), data=mass.data)
hist(resid(ses))
plot(predict(ses),resid(ses)) 

anova(ses)

contrast<-emmeans(ses, pairwise~Site.Type|Species, type="response")
contrast


contrast.matrix3 <- rbind(
  `C-U: one vs two` = c(0, 0, 0, 0, 1, 0),
  `C-U: one vs four` = c(0, 0, 0, 0, 0, 1),
  `C-U: two vs four` = c(0, 0, 0, 0, -1, 1))

comps3 <- glht(ses, contrast.matrix3)
summary(comps3)

ses<-lmer(sqrt(mono)~Site.Type*Species+ (1|Site), data=mass.data)
hist(resid(ses))
plot(predict(ses),resid(ses)) 

anova(ses)

contrast<-emmeans(ses, pairwise~Site.Type|Species, type="response")
contrast



contrast.matrix3 <- rbind(
  `C-U: one vs two` = c(0, 0, 0, 0, 1, 0),
  `C-U: one vs four` = c(0, 0, 0, 0, 0, 1),
  `C-U: two vs four` = c(0, 0, 0, 0, -1, 1))

comps3 <- glht(ses, contrast.matrix3)
summary(comps3)

ses<-lmer(sqrt(arom)~Site.Type*Species+ (1|Site), data=mass.data)
hist(resid(ses))
plot(predict(ses),resid(ses)) 

anova(ses)

contrast<-emmeans(ses, pairwise~Site.Type|Species, type="response")
contrast



contrast.matrix3 <- rbind(
  `C-U: one vs two` = c(0, 0, 0, 0, 1, 0),
  `C-U: one vs four` = c(0, 0, 0, 0, 0, 1),
  `C-U: two vs four` = c(0, 0, 0, 0, -1, 1))

comps3 <- glht(ses, contrast.matrix3)
summary(comps3)


#third, the nine compounds that were significantly correlated with the cap axes:

ses<-lmer(sqrt(er.2.amino.phenyl.ethenone)~Site.Type*Species+ (1|Site), data=mass.data)
ses<-lmer(log(er.alpha.pinene+0.001)~Site.Type*Species+ (1|Site), data=mass.data)
ses<-lmer(log(er.beta.pinene+0.001)~Site.Type*Species+ (1|Site), data=mass.data)
ses<-lmer(log(er.gamma.terpinene+0.001)~Site.Type*Species+ (1|Site), data=mass.data)
ses<-lmer(log(er.sabinene.hydrate+0.001)~Site.Type*Species+ (1|Site), data=mass.data)
ses<-lmer(log(er.methyl.nicotinate+0.001)~Site.Type*Species+ (1|Site), data=mass.data)
ses<-lmer(sqrt(er.trans.cinnamic.aldehyde)~Site.Type*Species+ (1|Site), data=mass.data)
ses<-lmer(sqrt(er.veratrole)~Site.Type*Species+ (1|Site), data=mass.data)
ses<-lmer(log(er.cis.3.hexenyl.acetate+0.001)~Site.Type*Species+ (1|Site), data=mass.data)
anova(ses)

contrast<-emmeans(ses, pairwise~Site.Type|Species, type="response")
contrast

  
  
  contrast.matrix3 <- rbind(
  `C-U: one vs two` = c(0, 0, 0, 0, 1, 0),
  `C-U: one vs four` = c(0, 0, 0, 0, 0, 1),
  `C-U: two vs four` = c(0, 0, 0, 0, -1, 1))
  
  comps3 <- glht(ses, contrast.matrix3)
summary(comps3)

