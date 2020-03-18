library(reshape2)

spme<-read.csv("data/SPME_data_no_blanks.csv",header=TRUE)

spme_wide<-dcast(spme, Sample~Name, sum, value.var="Area")
spme_ID<-read.csv("data/SPME_info.csv",header=TRUE)

#adding species and community info:
spme_wide<-cbind(spme_wide, spme_ID[,2:4])

#removing two compounds that are putative contaminants
spme_wide<-spme_wide[,c(-29, -32)]
#removing problematic characters from column name
colnames(spme_wide)[colnames(spme_wide)=="1.hexanol"] <- "er.1.hexanol"
colnames(spme_wide)[colnames(spme_wide)=="2.6.dimethyl.1.3.5.7.octatetraene.(trans)"] <- "er.2.6.dimethyl.1.3.5.7.octatetraene"
colnames(spme_wide)[colnames(spme_wide)=="2.amino.phenyl.ethenone"] <- "er.2.amino.phenyl.ethenone"
colnames(spme_wide)[colnames(spme_wide)=="6.Methyl.5.heptene.2.one"] <- "er.6.Methyl.5.heptene.2.one"
colnames(spme_wide)[colnames(spme_wide)=="unknown.C15H24.(1)"] <- "unknown.C15H24.1"

#removing two more contaminants
spme_wide<-spme_wide[,c(-6,-38)]

#this is the spme dataset. Out of R, I counted the number of compounds that each sample had in total, and in each of the four compound classes:
three.six<-read.csv("3_6_spmes.csv", header=TRUE)

#total:
wilcoxon(three.six[1:12,1], three.six[13:24,1], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#glvs:
test<-wilcox.test(three.six[1:12,7], three.six[13:24,7], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#aromatics:
test<-wilcox.test(three.six[1:12,6], three.six[13:24,6], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#sesquis:
wilcoxon(three.six[1:12,5], three.six[13:24,5], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#monos: 
wilcoxon(three.six[1:12,4], three.six[13:24,4], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#now generating a count of how many compounds are in each sample. Subtracting four because of the four ID columns
counts<-rowSums(spme_wide!=0)-4

test<-data.frame(counts, as.factor(spme_wide$Sample), spme_wide$Type)
test.3.6<-test[which(test$spme_wide.Type=="3.Fl"|test$spme_wide.Type=="6.Fl"),]
test.3.6<-test.3.6[order(test.3.6$spme_wide.Type),]

test<-wilcox.test(test.3.6[1:12,1], test.3.6[13:24,1], paired=TRUE)




wilcoxon(test.3.6[1:12,1], test.3.6[13:24,1], paired=TRUE)


#test[33:38,1]<-c(11,5,3,12,10,7)








petals<-spme_wide[which(spme_wide$Type=="Petals"|spme_wide$Type=="Not Petals"),]


spme_wide$mono<-(spme_wide$er.2.6.dimethyl.1.3.5.7.octatetraene+spme_wide$er.6.Methyl.5.heptene.2.one+spme_wide$alpha.pinene+spme_wide$alpha.terpinene+spme_wide$alpha.terpineol+spme_wide$alpha.terpinolene+spme_wide$beta.myrcene+spme_wide$beta.phellandrene + spme_wide$beta.pinene+spme_wide$cis.beta.ocimene+spme_wide$gamma.terpinene+spme_wide$limonene+
spme_wide$linalool+spme_wide$myroxide+spme_wide$para.cymene+spme_wide$sabine.hydrate+spme_wide$sabinene+spme_wide$terpinene.4.ol+spme_wide$trans.beta.ocimene+spme_wide$verbenone)
spme_wide$glv<-(spme_wide$er.1.hexanol+spme_wide$cis.3.hexen.1.ol +spme_wide$cis.3.hexenyl.acetate +spme_wide$trans.2.hexen.1.ol+spme_wide$trans.2.hexenal)
spme_wide$arom<-(spme_wide$er.2.amino.phenyl.ethenone+spme_wide$acetoin+spme_wide$benzyl.alcohol+spme_wide$methyl.salicylate+spme_wide$veratrole)
spme_wide$sesqui<-(spme_wide$alloaromadendrene+spme_wide$alpha.humulene+spme_wide$bergamotene+spme_wide$beta.cadinene+spme_wide$beta.longipinene+spme_wide$beta.sesquiphellandrene+spme_wide$cis.DMNT+spme_wide$gamma.cadinene+spme_wide$germacrene.D+spme_wide$intermedeol+spme_wide$trans.beta.caryophyllene+spme_wide$trans.trans.alpha.farnesene+spme_wide$unknown.C15H24.1+spme_wide$unknown.C15H24.2+spme_wide$unknown.C15H24.3)

CB323<-spme_wide[which(spme_wide$Site=="CB323"),]
#CB323<-CB323[,c(-29, -32)]
CB323.U.dose<-CB323[which(CB323$Type=="3.Fl"|CB323$Type=="6.Fl"),]
CB323.U.dose<-CB323.U.dose[order(CB323.U.dose$Type),]
CB323.U.dose.3<-CB323[which(CB323$Type=="3.Fl"),]
CB323.U.dose.6<-CB323[which(CB323$Type=="6.Fl"),]

ts <- sapply(2:46,function(x) (t.test(CB323.U.dose[1:3,x], CB323.U.dose[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
p.adjust(Pvalue)

ts <- sapply(50:53,function(x) (t.test(CB323.U.dose[1:3,x], CB323.U.dose[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)


SC.U<-spme_wide[which(spme_wide$Site=="SC"&spme_wide$Species=="Ung"),]
#SC.U<-SC.U[,c(-29, -32)]
SC.U.dose<-SC.U[which(SC.U$Type=="3.Fl"|SC.U$Type=="6.Fl"),]
SC.U.dose<-SC.U.dose[order(SC.U.dose$Type),]


ts <- sapply(2:46,function(x) (t.test(SC.U.dose[1:3,x], SC.U.dose[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
p.adjust(Pvalue)

ts <- sapply(50:53,function(x) (t.test(SC.U.dose[1:3,x], SC.U.dose[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)


SC.C<-spme_wide[which(spme_wide$Site=="SC"&spme_wide$Species=="Cyl"),]
#SC.C<-SC.C[,c(-29, -32)]
SC.C.dose<-SC.C[which(SC.C$Type=="3.Fl"|SC.C$Type=="6.Fl"),]
SC.C.dose<-SC.C.dose[order(SC.C.dose$Type),]


ts <- sapply(2:46,function(x) (t.test(SC.C.dose[1:3,x], SC.C.dose[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
p.adjust(Pvalue)

ts <- sapply(50:53,function(x) (t.test(SC.C.dose[1:3,x], SC.C.dose[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)

MHG.C<-spme_wide[which(spme_wide$Site=="MHG"&spme_wide$Species=="Cyl"),]
#MHG.C<-MHG.C[,c(-29, -32)]
MHG.C.dose<-MHG.C[which(MHG.C$Type=="3.Fl"|MHG.C$Type=="6.Fl"),]
MHG.C.dose<-MHG.C.dose[order(MHG.C.dose$Type),]


ts <- sapply(2:46,function(x) (t.test(MHG.C.dose[1:3,x], MHG.C.dose[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
p.adjust(Pvalue)


ts <- sapply(50:53,function(x) (t.test(MHG.C.dose[1:3,x], MHG.C.dose[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)


CB323<-spme_wide[which(spme_wide$Site=="CB323"),]
CB323<-CB323[,c(-29, -32)]
CB323.U.petals<-CB323[which(CB323$Type=="Petals"|CB323$Type=="Not Petals"),]
CB323.U.petals<-CB323.U.petals[order(CB323.U.petals$Type),]

ts <- sapply(2:46,function(x) (t.test((CB323.U.petals[1:3,x]), CB323.U.petals[4:6,x])))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
p.adjust(Pvalue)

ts <- sapply(50:53,function(x) (t.test(CB323.U.petals[1:3,x], CB323.U.petals[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)


SC.U<-spme_wide[which(spme_wide$Site=="SC"&spme_wide$Species=="Ung"),]
SC.U<-SC.U[,c(-29, -32)]
SC.U.petals<-SC.U[which(SC.U$Type=="Petals"|SC.U$Type=="Not Petals"),]
SC.U.petals<-SC.U.petals[order(SC.U.petals$Type),]

ts <- sapply(2:46,function(x) (t.test((SC.U.petals[1:3,x]), SC.U.petals[4:6,x])))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
p.adjust(Pvalue)

ts <- sapply(50:53,function(x) (t.test(SC.U.petals[1:3,x], SC.U.petals[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)


SC.C<-spme_wide[which(spme_wide$Site=="SC"&spme_wide$Species=="Cyl"),]
SC.C<-SC.C[,c(-29, -32)]
SC.C.petals<-SC.C[which(SC.C$Type=="Petals"|SC.C$Type=="Not Petals"),]
SC.C.petals<-SC.C.petals[order(SC.C.petals$Type),]

ts <- sapply(2:46,function(x) (t.test((SC.C.petals[1:3,x]), SC.C.petals[4:6,x],paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
p.adjust(Pvalue)

ts <- sapply(50:53,function(x) (t.test(SC.C.petals[1:3,x], SC.C.petals[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)

MHG.C<-spme_wide[which(spme_wide$Site=="MHG"&spme_wide$Species=="Cyl"),]
MHG.C<-MHG.C[,c(-29, -32)]
MHG.C.petals<-MHG.C[which(MHG.C$Type=="Petals"|MHG.C$Type=="Not Petals"),]
MHG.C.petals<-MHG.C.petals[order(MHG.C.petals$Type),]

ts <- sapply(2:46,function(x) (t.test((MHG.C.petals[1:3,x]), MHG.C.petals[4:6,x],paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
p.adjust(Pvalue)

ts <- sapply(50:53,function(x) (t.test(MHG.C.petals[1:3,x], MHG.C.petals[4:6,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)


C<-spme_wide[which(spme_wide$Species=="Cyl"),]
C.dose<-C[which(C$Type=="3.Fl"|C$Type=="6.Fl"),]
C.dose<-C.dose[order(C.dose$Type, C.dose$Site),]
C.petals<-C[which(C$Type=="Petals"|C$Type=="Not Petals"),]
C.petals<-C.petals[order(C.petals$Type, C.dose$Site),]


ts <- sapply(50:53,function(x) (t.test(C.petals[1:6,x], C.petals[7:12,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)


U<-spme_wide[which(spme_wide$Species=="Ung"),]
U.dose<-U[which(U$Type=="3.Fl"|U$Type=="6.Fl"),]
U.dose<-U.dose[order(U.dose$Type, U.dose$Site),]
U.petals<-U[which(U$Type=="Petals"|U$Type=="Not Petals"),]
U.petals<-U.petals[order(U.petals$Type, U.dose$Site),]


ts <- sapply(50:53,function(x) (t.test(U.dose[1:6,x], U.dose[7:12,x], paired=TRUE)))
stat<-t(as.matrix(ts[1,]))
Pvalue <- t(as.matrix(ts[3,]))
as.data.frame(t(stat))
as.data.frame(t(Pvalue))
p.adjust(Pvalue)



veg<-spme_wide[55:58,]
veg.1<-veg[1,]
veg.2<-veg[2,]
veg.3<-veg[3,]
veg.4<-veg[4,]

x<-t(as.data.frame(sapply(names(veg.1)[2:46], function(x) {
 veg[paste(x, "_pct")] <<- (veg.1[x] / sum(veg.1[2:46]))*100
 }
 )))
 
x2<-t(as.data.frame(sapply(names(veg.2)[2:46], function(x) {
 veg[paste(x, "_pct")] <<- (veg.2[x] / sum(veg.2[2:46]))*100
 }
 )))

x3<-t(as.data.frame(sapply(names(veg.3)[2:46], function(x) {
 veg[paste(x, "_pct")] <<- (veg.3[x] / sum(veg.3[2:46]))*100
 }
 )))
 x4<-t(as.data.frame(sapply(names(veg.4)[2:46], function(x) {
 veg[paste(x, "_pct")] <<- (veg.4[x] / sum(veg.4[2:46]))*100
 }
 )))
 
 x<-as.data.frame(t(x))
 x2<-as.data.frame(t(x2))
 x3<-as.data.frame(t(x3))
 x4<-as.data.frame(t(x4))
 
percents<-cbind(x,x2,x3,x4)
percents.t<-t(percents)
percents<-cbind(percents.t, veg$Site, veg$Species)

#what package is "gather" from?

per.cast<-gather(percents, compound, percent, 1:45,factor_key=TRUE)
colnames(per.cast)<-c("Site", "Species", "compound","percent")

per.cast.n<-per.cast[which(per.cast$percent!=0),]
per.cast.n$type<-c("acetoin", "mono", "mono", "mono", "sesqui", "mono", "mono", "mono", "cis 3 hexen 1 ol", "cis 3 hexen 1 ol","cis 3 hexen 1 ol","cis 3 hexenyl acetate","cis 3 hexenyl acetate","cis 3 hexenyl acetate","mono","sesqui","sesqui","sesqui","sesqui", "mono","mono", "mono", "GLV","GLV","mono","sesqui")

per.cast$ID<-c(as.character(per.cast$Site), as.character(per.cast$Species))

ggplot(aes(x="", y=percent, fill=type),data=per.cast.n)+
geom_bar(stat="identity")+coord_polar("y", start=0)+facet_grid(Species~Site)


#for petals:

matrix<-read.csv("SPME PA Matrix.csv", header=TRUE)
matrix.c<-matrix[which(matrix$Species=="Cyl",)]
scent.matrix<-matrix[,c(-1,-2)]
petals<-matrix$Sample.Type

scent.dist<-vegdist(scent.matrix, method='jaccard')

scent.div<-adonis2(scent.dist~petals, data=matrix, permutations = 999, method="jaccard", strata="PLOT")
scent.div

wilcoxon(matrix$Total[1:12], matrix$Total[13:24], paired=TRUE)
