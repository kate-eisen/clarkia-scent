#First, testing for differences between the two treatments in the 2019 wounding study

conc<-read.csv("data/glvs_2019_wide.csv", header=TRUE)


#separating out the rows for cis 3 hexenyl acetate, and then further separating by species

conc.ac<-conc[61:120,]
conc.ac.c<-conc.ac[which(conc.ac$Species=="Cyl"),]

conc.ac.u<-conc.ac[which(conc.ac$Species=="Ung"),]



conc.ol<-conc[1:60,]
conc.ol.c<-conc.ol[which(conc.ol$Species=="Cyl"),]
conc.ol.u<-conc.ol[which(conc.ol$Species=="Ung"),]


#tests for cis 3 hexen 1 ol for both species:
t.test(conc.ol.c$Con.W.t.mass, conc.ol.c$Con.N.t.mass, paired=TRUE)
t.test(conc.ol.u$Con.W.t.mass, conc.ol.u$Con.N.t.mass, paired=TRUE)

#tests for cis 3 hexenyl acetate for both species:
t.test(conc.ac.c$Con.W.t.mass, conc.ac.c$Con.N.t.mass, paired=TRUE)
t.test(conc.ac.u$Con.W.t.mass, conc.ac.u$Con.N.t.mass, paired=TRUE)



#Now comparing the 2019 control and wounded values to the 2018 values from the same species and site combinations:

conc.p<-read.csv("data/glvs_2018_2019.csv", header=TRUE)


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

#checking whether we need tests with unequal variance. And we do need that for cis 3 hexen 1 ol
var.test(conc.p.c.ol.OG$Con,conc.p.c.ol.W$Con)

#tests for cis 3 hexen ol for unguiculata
t.test(conc.p.u.ol.OG$Con,conc.p.u.ol.N$Con, var.equal=FALSE)
t.test(conc.p.u.ol.OG$Con,conc.p.u.ol.W$Con, var.equal=FALSE)

#tests for cis 3 hexen ol for cylindrica
t.test(conc.p.c.ol.OG$Con,conc.p.c.ol.N$Con, var.equal=FALSE)
t.test(conc.p.c.ol.OG$Con,conc.p.c.ol.W$Con, var.equal=FALSE)

#tests for cis 3 hexenyl acetate for unguiculata
t.test(conc.p.u.ac.OG$Con,conc.p.u.ac.N$Con, var.equal=FALSE)
t.test(conc.p.u.ac.OG$Con,conc.p.u.ac.W$Con, var.equal=FALSE)

#tests for cis 3 hexenyl acetate for cylindrica
t.test(conc.p.c.ac.OG$Con,conc.p.c.ac.N$Con, var.equal=FALSE)
t.test(conc.p.c.ac.OG$Con,conc.p.c.ac.W$Con, var.equal=FALSE)
