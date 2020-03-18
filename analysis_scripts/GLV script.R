setwd("../Scent/Follow up")

conc.r<-read.csv("concentrations_ratio.csv", header=TRUE)

ggplot(aes(x=Ratio.N, y=Ratio.W, color=Site, shape=Site), data=conc.r)+geom_point()+facet_wrap(~Species)

conc<-read.csv("concentrations.csv", header=TRUE)
conc<-conc[1:120,]

ggplot(aes(x=Site, y=Con, group=interaction(Type,Site), color=Type), data=conc.p.ol)+geom_boxplot()+facet_wrap(~Species)


conc.ac<-conc[61:120,]
conc.ac.c<-conc.ac[which(conc.ac$Species=="Cyl"),]
conc.ac.c.MCK<-conc.ac.c[which(conc.ac.c$Site=="MCK"),]
conc.ac.c.SC<-conc.ac.c[which(conc.ac.c$Site=="SC"),]
conc.ac.c.MHG<-conc.ac.c[which(conc.ac.c$Site=="MHG"),]
conc.ac.u<-conc.ac[which(conc.ac$Species=="Ung"),]
conc.ac.u.MCK<-conc.ac.u[which(conc.ac.u$Site=="MCK"),]
conc.ac.u.SC<-conc.ac.u[which(conc.ac.u$Site=="SC"),]
conc.ac.u.GRCO<-conc.ac.u[which(conc.ac.u$Site=="GRCO"),]





conc.ol<-conc[1:60,]
conc.ol.c<-conc.ol[which(conc.ol$Species=="Cyl"),]
conc.ol.c.MCK<-conc.ol.c[which(conc.ol.c$Site=="MCK"),]
conc.ol.c.SC<-conc.ol.c[which(conc.ol.c$Site=="SC"),]
conc.ol.c.MHG<-conc.ol.c[which(conc.ol.c$Site=="MHG"),]
conc.ol.u<-conc.ol[which(conc.ol$Species=="Ung"),]
conc.ol.u.MCK<-conc.ol.u[which(conc.ol.u$Site=="MCK"),]
conc.ol.u.GRCO<-conc.ol.u[which(conc.ol.u$Site=="GRCO"),]
conc.ol.u.SC<-conc.ol.u[which(conc.ol.u$Site=="SC"),]



t.test(conc.ac$Con.W.t.mass, conc.ac$Con.N.t.mass, paired=TRUE)


t.test(conc.ol.c$Con.W.t.mass, conc.ol.c$Con.N.t.mass, paired=TRUE)
t.test(conc.ol.u$Con.W.t.mass, conc.ol.u$Con.N.t.mass, paired=TRUE)

t.test(conc.ac.c$Con.W.t.mass, conc.ac.c$Con.N.t.mass, paired=TRUE)
t.test(conc.ac.u$Con.W.t.mass, conc.ac.u$Con.N.t.mass, paired=TRUE)


t.test(conc.ac.u.GRCO$Con.W.t.mass, conc.ac.u.GRCO$Con.N.t.mass, paired=TRUE)


conc.p<-read.csv("concentrations_plot2.csv", header=TRUE)




conc.p.c<-conc.p[which(conc.p$Species=="Cyl"),]
conc.p.c.ac<-conc.p.c[which(conc.p.c$Compound=="cis 3 hexenyl acetate"),]
conc.p.c.ac.N<-conc.p.c.ac[which(conc.p.c.ac$Type=="Con.N.t.mass"),]
conc.p.c.ac.W<-conc.p.c.ac[which(conc.p.c.ac$Type=="Con.W.t.mass"),]
conc.p.c.ac.OG<-conc.p.c.ac[which(conc.p.c.ac$Type=="First CG"),]


conc.p.c.ol<-conc.p.c[which(conc.p.c$Compound=="cis 3 hexen 1 ol"),]
conc.p.c.ol.N<-conc.p.c.ol[which(conc.p.c.ol$Type=="Con.N.t.mass"),]
conc.p.c.ol.W<-conc.p.c.ol[which(conc.p.c.ol$Type=="Con.W.t.mass"),]
conc.p.c.ol.OG<-conc.p.c.ol[which(conc.p.c.ol$Type=="First CG"),]


conc.p.u<-conc.p[which(conc.p$Species=="Ung"),]
conc.p.u.ac<-conc.p.u[which(conc.p.u$Compound=="cis 3 hexenyl acetate"),]
conc.p.u.ac.N<-conc.p.u.ac[which(conc.p.u.ac$Type=="Con.N.t.mass"),]
conc.p.u.ac.W<-conc.p.u.ac[which(conc.p.u.ac$Type=="Con.W.t.mass"),]
conc.p.u.ac.OG<-conc.p.u.ac[which(conc.p.u.ac$Type=="First CG"),]

conc.p.u.ol<-conc.p.u[which(conc.p.u$Compound=="cis 3 hexen 1 ol"),]
conc.p.u.ol.N<-conc.p.u.ol[which(conc.p.u.ol$Type=="Con.N.t.mass"),]
conc.p.u.ol.W<-conc.p.u.ol[which(conc.p.u.ol$Type=="Con.W.t.mass"),]
conc.p.u.ol.OG<-conc.p.u.ol[which(conc.p.u.ol$Type=="First CG"),]

var.test(conc.p.c.ol.OG$Con,conc.p.c.ol.W$Con)

t.test(conc.p.c.ol.OG$Con,conc.p.c.ol.W$Con, var.equal=FALSE)

a<-aov(lm(Con~Type, data=conc.p.c.ac))

plot(predict(a), residuals(a)) #this plot looks good
hist(residuals(a))


ggplot(aes(x=Site, y=Con, group=interaction(Type,Site), color=Type), data=conc.p.c)+geom_boxplot()+facet_wrap(~Compound)

conc.p.ac<-conc.p[which(conc.p$Compound=="cis 3 hexenyl acetate"),]
conc.p.ol<-conc.p[which(conc.p$Compound=="cis 3 hexen 1 ol"),]

tapply(conc.p$Con, list(conc.p$Type, conc.p$Site, conc.p$Species,conc.p$Compound), mean)


ggplot(aes(x=Site, y=Con, group=interaction(Type,Site), color=Type), data=conc.p.ol)+geom_boxplot()+facet_wrap(~Species)


conc.mean<-read.csv("Means.csv", header=TRUE)

conc.mean.ol<-conc.mean[which(conc.mean$Compound=="cis 3 hexen 1 ol"),]
conc.mean.ac<-conc.mean[which(conc.mean$Compound=="cis 3 hexenyl acetate"),]


ggplot(aes(x=Site, y=Con, group=interaction(Type,Site), color=Type), data=conc.mean.ac)+geom_point(position = position_dodge2(width=0.2))+facet_wrap(~Species)+geom_errorbar(aes(x=Site, ymin=Con-SE, ymax=Con+SE), width=0.2, position = position_dodge2(width=0.1))
