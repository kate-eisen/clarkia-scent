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
