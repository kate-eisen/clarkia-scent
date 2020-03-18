###Outline for Floral Scent Analysis

#basics:




library(multcomp)
library(lme4)
library(emmeans)

#the entire, mass-standardized, data file:

mass.data<-read.csv("data/all_data_standardized_mass.csv", header=TRUE)

##take out columns that don't need to be in the file; then adjust column number selections as needed throughout
mass.data$Site.Type<-factor(mass.data$Site.Type, levels=c("One", "Two", "Four"))

#adding in total scent:

mass.data$total<-rowSums(cbind(mass.data[,c(17:70)]))

#adding in type of compounds:

mass.data$mono<-(mass.data$er.2.6.dimethyl.1.3.5.7.octatetraene..trans.+mass.data$er.6.methyl.5.hepten.2.one +mass.data$er.alpha.pinene +mass.data$er.alpha.terpineol + mass.data$er.alpha.terpinene+ mass.data$er.alpha.terpinolene+ mass.data$er.beta.myrcene + mass.data$er.beta.phellandrene + mass.data$er.beta.pinene+mass.data$er.borneol + mass.data$er.cis.beta.ocimene + mass.data$er.gamma.terpinene + mass.data$er.geraniol + mass.data$er.limonene + mass.data$er.linalool + mass.data$er.myroxide + mass.data$er.para.cymene + mass.data$er.sabinene+ mass.data$er.sabinene.hydrate + mass.data$er.terpinen.4.ol + mass.data$er.trans.beta.ocimene + mass.data$er.verbenone)
mass.data$glv<-(mass.data$er.1.hexanol+mass.data$er.cis.3.hexen.1.ol+mass.data$er.cis.3.hexenyl.acetate+mass.data$er.cis.jasmone+ mass.data$er.trans.2.hexen.1.ol )
mass.data$arom<-(mass.data$er.2.amino.phenyl.ethenone+ mass.data$er.2.phenyl.ethanol+ mass.data$er.benzyl.acetate+ mass.data$er.benzyl.alcohol+ mass.data$er.cinnamic.alcohol+ mass.data$er.methyl.nicotinate+ mass.data$er.methyl.salicylate+ mass.data$er.trans.cinnamic.aldehyde+ mass.data$er.veratrole)
mass.data$sesqui<-(mass.data$er.alloaromadendrene+ mass.data$er.alpha.bergamotene + mass.data$er.alpha.humulene + mass.data$er.beta.cadinene + mass.data$er.beta.farnesene + mass.data$er.beta.longipinene + mass.data$er.caryophyllene.oxide + mass.data$er.DMNT + mass.data$er.gamma.cadinene + mass.data$er.germacrene.D + mass.data$er.intermedeol + mass.data$er.patchoulane + mass.data$er.trans.beta.caryophyllene + mass.data$er.trans.DMNT + mass.data$er.trans.trans.alpha.farnesene + mass.data$er.unknown.C15H24+mass.data$er.unknown.C15H24.2 + mass.data$er.unknown.C15H24.3)

#plotting by group means

means<-as.data.frame(sapply(72:75,function(x) tapply(mass.data[,x], list(mass.data$Species, mass.data$Site.Type), mean)))

colnames(means)<-c("mono", "glv", "arom", "sesqui")
Species<-c("Cyl", "Ung","Cyl", "Ung","Cyl", "Ung")
Site.Type<-c("One","One","Two","Two","Four","Four")
means$Species<-as.factor(Species)
means$Site.Type<-as.factor(Site.Type)
means$Site.Type<-factor(means$Site.Type, levels=c("One", "Two", "Four"))

ses<-as.data.frame(sapply(72:75,function(x) tapply(mass.data[,x], list(mass.data$Species, mass.data$Site.Type), std.error)))
colnames(ses)<-c("mono_se", "glv_se", "arom_se", "sesqui_se")

dat<-cbind(means, ses)



ggplot(aes(x=Site.Type, y=arom, color=Species), data=dat)+geom_point(size=10, position=position_dodge(width = 0.50))+geom_errorbar(aes(ymin= arom-arom_se, ymax= arom+arom_se), width=0, size=1, position=position_dodge(width = 0.50))+theme_classic(base_size=20)+ylab("")+xlab("Community Type")+ theme(legend.position = "none")

###

means<-as.data.frame(sapply(c(25,36,41,18,46,59,52,69,64), function(x) tapply(mass.data[,x], list(mass.data$Species, mass.data$Site.Type), mean)))

colnames(means)<-c("a.pinene", "b.pinene", "cis.hex", "amino.phenyl","gamma.terp","sab.hyd","meth.nic","ver", "cinn" )
Species<-c("Cyl", "Ung","Cyl", "Ung","Cyl", "Ung")
Site.Type<-c("One","One","Two","Two","Four","Four")
means$Species<-as.factor(Species)
means$Site.Type<-as.factor(Site.Type)
means$Site.Type<-factor(means$Site.Type, levels=c("One", "Two", "Four"))

ses<-as.data.frame(sapply(c(25,36,41,18,46,59,52,69,64),function(x) tapply(mass.data[,x], list(mass.data$Species, mass.data$Site.Type), std.error)))
colnames(ses)<-c("a.pinene_se", "b.pinene_se", "cis.hex_se", "amino.phenyl_se","gamma.terp_se","sab.hyd_se","meth.nic_se","ver_se", "cinn_se" )

dat<-cbind(means, ses)



ggplot(aes(x=Site.Type, y=ver, color=Species), data=dat)+geom_point(size=10, position=position_dodge(width = 0.50))+geom_errorbar(aes(ymin= ver-ver_se, ymax= ver+ ver_se), width=0, size=1, position=position_dodge(width = 0.50))+theme_classic(base_size=20)+ylab("")+xlab("Community type")+ theme(legend.position = "none")






mode<-lmer(er.trans.beta.ocimene~Site.Type*Species+(1|Site), data=mass.data)
hist(resid(mode))
plot(predict(mode),resid(mode))
anova(mode)


emmeans(mode, "Species",type="response" )
contrast<-emmeans(mode, pairwise~Species|Site.Type, type="response")
contrast
contrast<-emmeans(mode, pairwise~Site.Type|Species, type="response")

contrast.matrix3 <- rbind(
  `C vs. U: One vs Two` = c(0, 0, 0, 0, -1, 0),
  `C vs. U: One vs Four` = c(0, 0, 0, 0, 0, -1),
  `C vs. U: Two vs Four` = c(0, 0, 0, 0, 1, -1))
  
  comps3 <- glht(mode, contrast.matrix3)
summary(comps3)

peas<-read.csv("P_vals.csv", header=TRUE)
peas.adj<-p.adjust(peas$p, method="BH")
peas$p.adj<-peas.adj


tapply(mass.data$er.gamma.terpinene, list( mass.data$Site.Type, mass.data$Species), mean)

plotting<-read.csv("Sig_means_plotting.csv", header=TRUE)
plotting<-plotting[13:36,]
plotting$Site.Type<-factor(plotting$Site.Type, levels=c("One", "Two", "Four"))

ggplot(aes(x=Site.Type, y=Mean, color=Species, fill=Species, color=Species, shape=Compound), data=plotting)+geom_point(size=10, position=position_dodge(width=1))+theme_classic()+geom_errorbar(aes(x=Site.Type, ymin=Mean-SE, ymax=Mean+SE),data=plotting, position=position_dodge(width=1), color="#B8A7A4")+scale_shape_manual(values=c(15:18))+scale_color_manual(values=c("#0E2057","#A51647"))+scale_fill_manual(values=c("#0E2057","#A51647"))


mod<-hurdle(er.trans.beta.ocimene~Site.Type*Species, data=mass.data)


ocimene<-mass.data[which(mass.data$er.trans.beta.ocimene!=0),]

km1 <- mixed_model(er.trans.beta.ocimene ~ Site.Type*Species, random = ~ 1 | Site, data = mass.data, 
                  family = hurdle.lognormal(), n_phis = 1,
                  zi_fixed = ~ Species)
                  
                  km2 <- mixed_model(er.trans.beta.ocimene ~ Site.Type+Species, random = ~ 1 | Site, data = mass.data, 
                  family = hurdle.lognormal(), n_phis = 1,
                  zi_fixed = ~ Species)

#need to figure out how to run tweedie

#also, is there some kind of hurdle model that is better for this data?
hist(resid(mode))
plot(predict(mode),resid(mode)) 

#compounds only:

mass.data.compounds<-mass.data[,17:70]

#by species:

mass.data.sp<-mass.data[,-c(1:14,16,71)]
mass.data.sp$Species<-as.factor(mass.data.sp$Species)

#working with unguiculata only:

mass.data.u<-mass.data[which(mass.data$Species=="Ung"),]
mass.data.u.comp<-mass.data.u[,c(17:70)]


counts <- sapply(17:70,function(x) length(which(mass.data.u[,x]!=0)))
counts.u<-data.frame(counts, compounds=colnames(mass.data.u)[17:70])



mass.data.u.one<-mass.data.u[which(mass.data.u$Site.Type=="One"),]
mass.data.u.two<-mass.data.u[which(mass.data.u$Site.Type=="Two"),]
mass.data.u.four<-mass.data.u[which(mass.data.u$Site.Type=="Four"),]

p<-pint.p(mass.data.u.one[,c(17:38,40:52,54:63,65:68,70)], tails=2)

p2<-pint.p(mass.data.u.two[,c(17:38,40:56, 58:63,65:68,70)], tails=2)
p3<-pint.p(mass.data.u.four[,c(17:38,40:56, 58,60:68)], tails=2)

#working with cylindrica only:

mass.data.c<-mass.data[which(mass.data$Species=="Cyl"),]
mass.data.c.comp<-mass.data.c[,c(17:70)]

mass.data.c.one<-mass.data.c[which(mass.data.c$Site.Type=="One"),]
mass.data.c.two<-mass.data.c[which(mass.data.c$Site.Type=="Two"),]
mass.data.c.four<-mass.data.c[which(mass.data.c$Site.Type=="Four"),]


p<-pint.p(mass.data.c.one[,c(17:18,20:27, 29:51,53:55,58,61:64,66,68,69)], tails=2)

p2<-pint.p(mass.data.c.two[,c(17,19:27, 29:51,53:55,57,58:59,61:64,65:66,68,69)],tails=2)
p3<-pint.p(mass.data.c.four[,c(17,19:43,46:52,54:55, 57:58,61:64,66:69)],tails=2)


counts <- sapply(17:70,function(x) length(which(mass.data.c[,x]!=0)))

sum<-sapply(17:70,function(x) sum(mass.data.u.one[,x]))
counts.c<-data.frame(counts, compounds=colnames(mass.data.c)[17:70])


#by the site types:

mass.data.one<-mass.data[which(mass.data$Site.Type=="One"),]
mass.data.one.comp<-mass.data.one[,c(17:70)]

mass.data.two<-mass.data[which(mass.data$Site.Type=="Two"),]
mass.data.two.comp<-mass.data.two[,c(17:70)]

mass.data.four<-mass.data[which(mass.data$Site.Type=="Four"),]
mass.data.four.comp<-mass.data.four[,c(17:70)]


#Frameworks:

#A. Can analyze data with all 54 compounds by species x site type

mass.data.comp<-read.csv("all_data_edit_standardized_er_mass_vc_edit_F_composite_compounds.csv", header=TRUE)

mass.data.comp.compounds<-mass.data.comp[,c(17:47)]

scents.permanova.u<-adonis(mass.data.u.comp~Site.Type/Site, data=mass.data.u, strata=mass.data.u$Site, permutations=999, method="bray")
scents.permanova.c<-adonis(mass.data.c.comp~Site.Type/Site, data=mass.data.c, strata=mass.data.c$Site, permutations=999, method="bray")

scents.permanova.one<-adonis(mass.data.one.comp~Species, data=mass.data.one, strata=mass.data.one$Site, permutations=999, method="bray")
scents.permanova.two<-adonis(mass.data.two.comp~Species, data=mass.data.two, strata=mass.data.two$Site, permutations=999, method="bray")
scents.permanova.four<-adonis(mass.data.four.comp~Species, data=mass.data.four, strata=mass.data.four$Site, permutations=999, method="bray")

scents.permanova<-adonis(mass.data.compounds~Site.Type*Species, data=mass.data, strata=mass.data$Site.Type, permutations=999, method="bray")


#trying this with z scores

zscores<-read.csv("all_data_edit_standardized_er_mass_vc_edit_F_z_scores_final.csv", header=TRUE)
zscores.compounds<-zscores[,c(17:70)]
scents.permanova<-adonis(zscores.compounds~Site.Type*Species, data=zscores, strata=zscores$Site.Type, permutations=999, method="euclidean")



#if worried about site to site variability, could standardize the volatiles at the site level with z scores



scents.permanova<-adonis(mass.data.comp.compounds~Site.Type*Species, data=mass.data.comp, strata=mass.data.comp$Site, permutations=999, method="bray")


pairwise.adonis2(mass.data.compounds~Species*Site.Type/Site, data=mass.data, strata="Site")

mass.data$Code<-as.factor(paste(mass.data$Species,mass.data$Site.Type))
mass.data.comp$Code<-as.factor(paste(mass.data.comp$Species,mass.data.comp$Site.Type))
mass.data$Code.2<-as.factor(paste(mass.data$Species, mass.data$Site))

#can run this with either the composites or the full dataset:
scents.cap <-capscale(mass.data.compounds ~ Code, mass.data, dist="bray") 
scents.cap
plot(scents.cap, display=c("sp", "sites","cn"),type="text") 
plot(scents.cap, display=c("sp", "cn"),type="text")
plot(scents.cap, display=c("sp"),type="text")
anova(scents.cap) 

##trying without trans beta ocimene and cis 3 hexenyl acetate b/c of scale issues:

mass.data.c.sub<-mass.data.compounds[,-c(25,47)]

scents.permanova<-adonis(mass.data.c.sub~Site.Type*Species, data=mass.data, strata=mass.data$Site, permutations=999, method="bray")

mass.data$Code<-as.factor(paste(mass.data$Species,mass.data$Site.Type))

scents.cap <-capscale(mass.data.compounds ~ Code, mass.data, dist="bray") 
scents.cap
plot(scents.cap, type="n", xlim=c(-0.5,0.5))

ordiellipse(scents.cap, groups=mass.data$Code, col=color, kind="se",conf=0.95 , draw="polygon", alpha=50, border=color)
ordibar(scents.cap, groups=mass.data$Code,  col = color, kind="se" )
ordispider(scents.cap, groups=mass.data$Code,  col = color )

plot(scents.cap, type="n", xlim=c(-0.5,0.5))

points(scents.cap, pch=c(21,22)[mass.data$Species], bg= color[mass.data$Code],cex=1.5)

#need to figure out how to better define legend
legend(-3, -2, legend=levels(mass.data$Code), pch=c(21,21,21,22,22,22), col=color, cex=1)

color<-c("red","purple", "green", "blue", "pink", "grey")

color<-c("darkorchid3", "orchid1", "mediumorchid", "deeppink2", "hotpink1", "violetred3")

color<-c("lightslateblue", "firebrick1", "goldenrod1", "lightslateblue", "firebrick1", "goldenrod1")

color<-c("dodgerblue1", "firebrick1", "goldenrod1", "cadetblue1", "indianred1", "lightgoldenrod1")

color.sites<-c("firebrick1", "dodgerblue1", "cadetblue1", "goldenrod1", "indianred1", "springgreen2", "royalblue2", "seagreen3", "lightgoldenrod1", "brown2", "gold", "darkolivegreen3")

pl<-plot(scents.cap, type="n", correlation=TRUE, xlim=c(-0.5,0.5))
#with(mass.data, points(pl, "site",pch=c(21,22)[mass.data$Species], bg= color[mass.data$Code],cex=1.5 ))
with(mass.data, points(pl, "site",pch=c(21,22)[mass.data$Species], bg= color.sites[mass.data$Site],cex=1.5 ))

text(pl, "sp", scaling=1, arrow=TRUE, length=0.05, col=4, cex=0.6, xpd=TRUE)


points(scents.cap, pch=c(21,22)[mass.data$Species], bg= color[mass.data$Code],cex=1.5)
text(scents.cap, "sp", arrow=TRUE, length=0.05, col=4, cex=0.6, xpd=TRUE)







points(scents.cap, type="points", display="sites", col=color[mass.data$Code]) 
plot(scents.cap, display=c("cn"),type="text")
anova(scents.cap) 

scores<-read.csv("CAP_scores.csv", header=TRUE)
#or
scores<-read.csv("CAP_scores_composites.csv", header=TRUE)

#or
scores<-read.csv("CAP_scores_mDist.csv", header=TRUE)


mass.data.compounds.scores<-cbind(mass.data.compounds, scores)
#or
mass.data.compounds.scores<-cbind(mass.data.comp.compounds, scores)



##figure out how to make a matrix into a dataframe that works for this!!


#for the composite compounds, there are now 31 compound collumns
cors <- as.data.frame(sapply(1:31,function(x) cor.test(mass.data.compounds.scores[,x],mass.data.compounds.scores$CAP1)))

cors2 <- as.data.frame(sapply(1:31,function(x) cor.test(mass.data.compounds.scores[,x],mass.data.compounds.scores$CAP2)))
pvals<-sapply(1:31,function(x) print(cors[,x]$p.value))
pvals2<-sapply(1:31,function(x) print(cors2[,x]$p.value))
pvals.a<-p.adjust(pvals, method="BH")
pvals.a2<-p.adjust(pvals2, method="BH")


co1<-sapply(1:31,function(x) print(cors[,x]$estimate))
co2<-sapply(1:31,function(x) print(cors2[,x]$estimate))

compounds<-colnames(mass.data.comp.compounds)


c.table<-as.data.frame(cbind(compounds,cor1=co1,cap1=pvals.a,cor2=co2, cap2=pvals.a2))
c.table$cap1<-as.numeric(as.character(c.table$cap1))
c.table$cap2<-as.numeric(as.character(c.table$cap2))
c.table$sig1<-rep("No", times=dim(c.table)[1])
c.table$sig2<-rep("No", times=dim(c.table)[1])
c.table$sig1<-ifelse(c.table$cap1 < 0.01, "Yes", c.table$sig1)
c.table$sig2<-ifelse(c.table$cap2 < 0.01, "Yes", c.table$sig2)


###


cors <- as.data.frame(sapply(1:54,function(x) cor.test(mass.data.compounds.scores[,x],mass.data.compounds.scores$CAP1)))

cors2 <- as.data.frame(sapply(1:54,function(x) cor.test(mass.data.compounds.scores[,x],mass.data.compounds.scores$CAP2)))
pvals<-sapply(1:54,function(x) print(cors[,x]$p.value))
pvals2<-sapply(1:54,function(x) print(cors2[,x]$p.value))
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


