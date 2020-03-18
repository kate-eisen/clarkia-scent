###Outline for Floral Scent Analysis

#basics:

library(vegan)
library(multcomp)
library(lme4)
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

#Frameworks:

#A. Can analyze data with all 54 compounds by species x site type
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


co1<-sapply(1:54,function(x) print(cors[,x]$estimate))
co2<-sapply(1:54,function(x) print(cors2[,x]$estimate))


compounds<-colnames(mass.data.compounds)

c.table<-as.data.frame(cbind(compounds,cor1=co1,cap1=pvals.a,cor2=co2, cap2=pvals.a2))
c.table$cap1<-as.numeric(as.character(c.table$cap1))
c.table$cap2<-as.numeric(as.character(c.table$cap2))
c.table$sig1<-rep("No", times=dim(c.table)[1])
c.table$sig2<-rep("No", times=dim(c.table)[1])
c.table$sig1<-ifelse(c.table$cap1 < 0.01, "Yes", c.table$sig1)
c.table$sig2<-ifelse(c.table$cap2 < 0.01, "Yes", c.table$sig2)


ggplot(aes(x=Site.Type, y=er.veratrole, fill=Species), data=mass.data)+geom_boxplot()

model<-lmer(er.veratrole~Site.Type*Species+(1|Site), data=mass.data)


##
Code<-as.factor(mass.data$Code)
#Y<-cbind(mass.data.compounds[1:54])
manova.scents<-manova(as.matrix(mass.data.compounds)~Code)
scents.can <- candisc(manova.scents, data=mass.data, scores=TRUE)
coefs<-as.data.frame(scents.can$coeffs.std)



greedy.wilks(mass.data.compounds, Code, data=mass.data)
stepclass(mass.data.compounds, grouping=mass.data$Code, method="lda")



#B. Can analyze data with a subset of the 54 compounds via variable selection by species x site type

#C. Can analyze total scent emission by species x site type

total<-lmer(sqrt(total)~Site.Type*Species+ (1|Site), data=mass.data)

hist(resid(total))
plot(predict(total),resid(total)) 
summary(total)
anova(total)

ggplot(aes(x=Site.Type, y=er.alpha.pinene, fill=Species), data=mass.data)+geom_boxplot()+xlab("")+ylab("Aromatics scent emission")+theme_classic()

emmeans(total, "Species",type="response" )
contrast<-emmeans(total, pairwise~Species|Site.Type, type="response")
contrast<-emmeans(total, pairwise~Site.Type|Species, type="response")



#D. Can analyze scent emission in categories (e.g. monoterpenes, sesequiterpenes, etc.) by species x site types

ses<-lmer(log(glv+0.0001)~Site.Type*Species+ (1|Site), data=mass.data)
hist(resid(ses))
plot(predict(ses),resid(ses)) 

ses<-lmer(log(sesqui+0.0001)~Site.Type*Species+ (1|Site), data=mass.data)
hist(resid(ses))
plot(predict(ses),resid(ses)) 

ses<-lmer(sqrt(mono)~Site.Type*Species+ (1|Site), data=mass.data)
hist(resid(ses))
plot(predict(ses),resid(ses)) 

ses<-lmer(sqrt(arom)~Site.Type*Species+ (1|Site), data=mass.data)
hist(resid(ses))
plot(predict(ses),resid(ses)) 

mass.data$c.comps<-mass.data$er.veratrole+mass.data$er.trans.cinnamic.aldehyde
mass.data$u.comps<-mass.data$er.2.amino.phenyl.ethenone+
mass.data$er.alpha.pinene+mass.data$er.beta.pinene+
mass.data$er.gamma.terpinene+mass.data$er.sabinene.hydrate+
mass.data$er.methyl.nicotinate

size.data<-read.csv("../GH CG 2015/CG_All_Traits_Feb_2016.csv", header=TRUE)
size.data.cu<-size.data[which(size.data$Sp=="Cyl"|size.data$Sp=="Ung"),]
size.data.cu<-size.data.cu[!is.na(size.data.cu$Petal.area.mm),]

tapply(size.data.cu$Petal.area.mm/10, list(size.data.cu$Sp, size.data.cu$Site.Type), mean)


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

contrast<-emmeans(ses, pairwise~Species|Site.Type, type="response")
contrast

#contrast.matrix3 <- rbind(
 # `C vs. U: Four vs One` = c(0, 0, 0, 0, -1, 0),
  #`C vs. U: Four vs Two` = c(0, 0, 0, 0, 0, -1),
 # `C vs. U: One vs Two` = c(0, 0, 0, 0, 1, -1))
  
  
  contrast.matrix3 <- rbind(
  `C-U: alone vs together` = c(0, 0, 0, 0, 1, 0),
  `C-U: alone vs with friends` = c(0, 0, 0, 0, 0, 1),
  `C-U: together vs with friends` = c(0, 0, 0, 0, -1, 1))
  
  comps3 <- glht(ses, contrast.matrix3)
summary(comps3)





#Main Questions:

#1. Are there differences across species or site types in emission rates?

#1A. What is responsible for this variation?
#1B. Do the most important/most emitted compounds vary across the site types?
#1C. Do the compounds that are important in the profiles of both species change across site types?
#1D. Do the compounds that are unique to each species change across site types?


#Other questions:

#1. Where are compounds made?

#1A. Are the same compounds made in the same floral parts across the species?

#2. Does the correlation between compounds change across species and site types?

#use package "cocor"

library(cocor)
#comparision a: alpha humulene and caryophyllene
cocor(~er.alpha.humulene +er.trans.beta.caryophyllene|er.alpha.humulene +er.trans.beta.caryophyllene, data=list(mass.data.u.four, mass.data.u.two))

#comparison b: GLVs
cocor(~er.cis.3.hexen.1.ol +er.cis.3.hexenyl.acetate |er.cis.3.hexen.1.ol +er.cis.3.hexenyl.acetate, data=list(mass.data.c.one, mass.data.c.four))

cocor(~mono +sesqui |mono +sesqui, data=list(mass.data.u.one, mass.data.u.four))

#glvs are positively correlated in U but not in C

  

#questions:


##for adonis, how do you determine what variables are driving the observed results?

##which compounds do we expect to be strongly correlated?
#a. alpha humulene and caryophyllene
#b. the GLVs
#c. classes of compounds, eg. monoterpenes and sesquiterpenes?


