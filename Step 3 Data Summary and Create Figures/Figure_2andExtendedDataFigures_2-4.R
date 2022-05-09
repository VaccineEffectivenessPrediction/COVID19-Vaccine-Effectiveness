df <- read.csv("VE and whole genome mismatch.csv", header = TRUE, fill = TRUE)
df$type<- factor(df$type,c("efficacy","effectiveness"))
df$platforms <- factor(df$platforms,c("mRNA","Protein subunit","Inactivated","Viral vector"))

# scatter plot for each vaccine platform
cols_lmsp <- c("#00BFC4","#F8766D","#FFB90F","#C77CFF")  
lm.sp <- function(df1,x1,y1,group1,type1,title) {
  splot <- ggplot() + 
    stat_smooth(data=df1,aes(x = x1, y = y1,colour=group1,fill=group1),size=1,alpha=0.2) + 
    stat_smooth(data=df1, aes(x = x1, y = y1),method=lm,size=0.5,colour="grey50",fullrange = TRUE,
                se=FALSE,linetype="dashed") + 
    geom_point(data=df1, aes(x = x1, y = y1,color=group1,shape=type1),size=3.6) 
  
  splot.edited <- splot +
    scale_y_continuous(name="COVID-19 VE (%)",breaks=seq(0, 100, 25)) + scale_x_continuous(name = "Genetic distance") + 
    scale_colour_manual(values = cols_lmsp) + theme_classic() + 
    scale_fill_manual(values = cols_lmsp) +
    labs(title = title) + 
    theme(plot.title = element_text(size=12, face="bold")) + 
    theme(legend.position = "none") +
    expand_limits(y = 0)
  print(splot.edited)
  ggsave(filename = paste(title,".png",sep=""), plot = splot.edited,width = 4.28, height =3)
}

lm.sp(df,df$RBD,df$VE,df$platforms,df$type,"(a) RBD")
lm.sp(df,df$NTD,df$VE,df$platforms,df$type,"(a) NTD")
lm.sp(df,df$All,df$VE,df$platforms,df$type,"(b) S protein")

# scatter plot for each vaccine product
sp <- function(df1,x1,y1,type1,title) {
  splot <- ggplot() +
    geom_point(data=df1, aes(x = x1, y = y1,shape=type1),color="#00C091",size=3.6) +
    stat_smooth(aes(x = x1, y = y1),color="#00C091",fill="#00C091",size=1,alpha=0.2)
  splot.edited <- splot + 
    scale_y_continuous(name="COVID-19 VE (%)",breaks=seq(0, 100, 25)) + scale_x_continuous(name = "Genetic distance",breaks=seq(0, 3, 1)) +
    theme_classic() + labs(title = title) +
    theme(plot.title = element_text(size=12, face="bold")) +
    theme(legend.position = "none") +
    expand_limits(y = c(0,100),x=c(0,3))
  print(splot.edited)
  ggsave(filename = paste(title,".png",sep=""), plot = splot.edited,width = 4.28, height =3)
}

sp.nse <- function(df1,x1,y1,type1,title) {
  splot <- ggplot(data=df1, aes(x = x1, y = y1,shape=type1)) +
    geom_point(color="#00C091",size=3.6)
  
  splot.edited <- splot + 
    scale_y_continuous(name="COVID-19 VE (%)",breaks=seq(0, 100, 25)) + scale_x_continuous(name = "Genetic distance",breaks=seq(0, 3, 1)) +
    theme_classic() + labs(title = title) +
    theme(plot.title = element_text(size=12, face="bold")) +
    theme(legend.position = "none") +
    expand_limits(y = c(0,100),x=c(0,3))
  print(splot.edited)
  ggsave(filename = paste(title,".png",sep=""), plot = splot.edited,width = 4.28, height =3)
}

df1 <- subset(df,developer=="BNT")
sp(df1,df1$RBD,df1$VE,df1$type,"Pfizer–BioNTech: BNT162b2")
df1 <- subset(df,developer=="Janssen")
sp(df1,df1$RBD,df1$VE,df1$type,"Janssen: Ad26.COV2.S")
df1 <- subset(df,developer=="Moderna" )
sp(df1,df1$RBD,df1$VE,df1$type,"Moderna: mRNA-1273")
df1 <- subset(df,developer=="Novavax")
sp(df1,df1$RBD,df1$VE,df1$type,"Novavax: NVX-CoV2373")
df1 <- subset(df,developer=="Oxford")
sp(df1,df1$RBD,df1$VE,df1$type,"Oxford–AstraZeneca: AZD1222")
df1 <- subset(df,developer=="Sinovac")
sp(df1,df1$RBD,df1$VE,df1$type,"Sinovac: CoronaVac")

df1 <- subset(df,developer=="CanSino")
sp.nse(df1,df1$RBD,df1$VE,df1$type,"CanSino: Convidecia")
df1 <- subset(df,developer=="Gamaleya")
sp.nse(df1,df1$RBD,df1$VE,df1$type,"Gamaleya: Sputnik V")
df1 <- subset(df,developer=="Bharat")
sp.nse(df1,df1$RBD,df1$VE,df1$type,"Bharat: Covaxin")
df1 <- subset(df,developer=="Sinopharm")
sp.nse(df1,df1$RBD,df1$VE,df1$type,"Sinopharm: BBIBP-CorV")


library(lmerTest)
p_model <- lmer(VE ~ RBD+(RBD|platforms) , data=df) # age + TimeAfterSecondDose
df1 <- subset(df,developer=="BNT"|"Janssen"|"Moderna"|"Novavax"|"Oxford"|"Sinovac")
d_model <- lmer(VE ~ RBD+(RBD|developer), data=df1) # age + TimeAfterSecondDose

cols_psp <- c("#00BFC4","#F8766D","#FFB90F","#C77CFF") 
p.sp <- function(df1,x1,y1,group1,type1,title) {
  splot <- ggplot() + geom_point(data=df1, aes(x = x1, y = y1,colour=group1,shape=type1),size=3.6) 
  
  splot.edited <- splot + 
    scale_y_continuous(name="COVID-19 VE (%)") + scale_x_continuous(name = "Genetic distance") + 
    scale_colour_manual(values = cols_psp) + 
    theme_classic() + labs(title = title) + 
    theme(plot.title = element_text(size=12, face="bold")) + 
    theme(legend.position = "none") + expand_limits(y = 0)
  print(splot.edited)
  ggsave(filename = paste(title,".png",sep=""), plot = splot.edited,width = 4.28, height =3)
}

p.sp(df,df$ORF1a,df$VE,df$platforms,df$type,"(a) ORF1a")
p.sp(df,df$ORF1b,df$VE,df$platforms,df$type,"(b) ORF1b")
p.sp(df,df$ORF3a,df$VE,df$platforms,df$type,"(c) ORF3a")
p.sp(df,df$E,df$VE,df$platforms,df$type,"(d) E protein")
p.sp(df,df$M,df$VE,df$platforms,df$type,"(e) M protein")
p.sp(df,df$ORF6,df$VE,df$platforms,df$type,"(f) ORF6")
p.sp(df,df$ORF7a,df$VE,df$platforms,df$type,"(g) ORF7a")
p.sp(df,df$ORF7b,df$VE,df$platforms,df$type,"(h) ORF7b")
p.sp(df,df$ORF8,df$VE,df$platforms,df$type,"(i) ORF8")
p.sp(df,df$N,df$VE,df$platforms,df$type,"(j) N protein")
p.sp(df,df$ORF10,df$VE,df$platforms,df$type,"(k) ORF10")


p.sp(df,df$nsp1,df$VE,df$platforms,df$type,"(a) NSP1") 
p.sp(df,df$nsp2,df$VE,df$platforms,df$type,"(b) NSP2")
p.sp(df,df$nsp3,df$VE,df$platforms,df$type,"(c) NSP3")
p.sp(df,df$nsp4,df$VE,df$platforms,df$type,"(d) NSP4")
p.sp(df,df$nsp5,df$VE,df$platforms,df$type,"(e) NSP5")
p.sp(df,df$nsp6,df$VE,df$platforms,df$type,"(f) NSP6")
p.sp(df,df$nsp7,df$VE,df$platforms,df$type,"(g) NSP7")
p.sp(df,df$nsp8,df$VE,df$platforms,df$type,"(h) NSP8")
p.sp(df,df$nsp9,df$VE,df$platforms,df$type,"(i) NSP9")
p.sp(df,df$nsp10,df$VE,df$platforms,df$type,"(j) NSP10")
p.sp(df,df$nsp11,df$VE,df$platforms,df$type,"(k) NSP11")
p.sp(df,df$rdrp,df$VE,df$platforms,df$type,"(l) RdRp")
p.sp(df,df$nsp13,df$VE,df$platforms,df$type,"(m) NSP13")
p.sp(df,df$nsp14,df$VE,df$platforms,df$type,"(n) NSP14")
p.sp(df,df$nsp15,df$VE,df$platforms,df$type,"(o) NSP15")
p.sp(df,df$nsp16,df$VE,df$platforms,df$type,"(p) NSP16")

