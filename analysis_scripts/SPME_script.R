
library(dplyr)

#for comparing the three vs six flower samples, I counted the number of compounds that each sample had in total, and in each of the four compound classes:
three.six<-read.csv("data/3_6_spmes.csv", header=TRUE)

#getting means:
three.six %>% group_by(Type) %>% 
  summarise(m.Total=mean(Total.N), m.Mono=mean(N.mono), m.Sesqui=mean(N.sesqui), m.arom=mean(N.arom), m.GLV=mean(N.glv))

#total:
test <- wilcox.test(three.six[1:12,1], three.six[13:24,1], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#glvs:
test<-wilcox.test(three.six[1:12,7], three.six[13:24,7], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#aromatics:
test<-wilcox.test(three.six[1:12,6], three.six[13:24,6], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#sesquis:
test <- wilcox.test(three.six[1:12,5], three.six[13:24,5], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#monos: 
wilcoxon(three.six[1:12,4], three.six[13:24,4], paired=TRUE)
Zstat<-qnorm(test$p.value/2)

#for petals:

matrix<-read.csv("data/spme_petals_matrix.csv", header=TRUE)

#removing totals and sample ID info
scent.matrix<-matrix[,c(10:50)]
petals<-matrix$Sample.Type

scent.dist<-vegdist(scent.matrix, method='jaccard')

scent.div<-adonis2(scent.dist~petals, data=matrix, permutations = 999, method="jaccard", strata="PLOT")
scent.div

matrix<-matrix[order(matrix$Sample.Type),]

matrix %>% group_by(Sample.Type) %>% 
  summarise(m.Total=mean(Total))

test<-wilcox.test(matrix$Total[1:12], matrix$Total[13:24], paired=TRUE)
Zstat<-qnorm(test$p.value/2)