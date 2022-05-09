df <- read.csv("VE and whole genome mismatch.csv", header = TRUE, fill = TRUE)
df$type<- factor(df$type,c("efficacy","effectiveness"))
df$platforms <- factor(df$platforms,c("mRNA","Protein subunit","Inactivated","Viral vector"))

# scatter plot for each vaccine platform
cols_lmsp <- c("#00BFC4","#F8766D","#FFB90F","#C77CFF")  
lm.sp <- function(df1,x1,y1,group1,type1,title) {
  splot <- ggplot() + 
    stat_smooth(data=df1,aes(x = x1, y = y1,colour=group1,fill=group1),method,size=1,alpha=0.2) + 
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
    stat_smooth(aes(x = x1, y = y1),method,color="#00C091",fill="#00C091",size=1,alpha=0.2)
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
