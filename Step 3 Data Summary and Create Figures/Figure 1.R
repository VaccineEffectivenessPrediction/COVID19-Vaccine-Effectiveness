
df <- read.csv("VE and whole genome mismatch.csv", header = TRUE, fill = TRUE)
df$type<- factor(df$type,c("efficacy","effectiveness"))
df$platforms <- factor(df$platforms,c("mRNA","Protein subunit","Inactivated","Viral vector"))
# VE plot
cols <- c("#ffc8b9","#99CCFF")
p <- ggplot(df,aes(x=platforms,y=VE,fill=type)) + 
  geom_violin(width=0.8, color="white",trim = FALSE,position=position_dodge(0.8),scale = "width") + 
  geom_boxplot(width=0.08,outlier.colour = NA,size=0.3,position=position_dodge(0.8)) + 
  stat_summary(fun  = mean,geom="point",position=position_dodge(0.8),colour="white",size=0.7)
vp <- p  + 
  theme_classic() +
  labs(x = "Vaccine platforms", y = "COVID-19 VE (%)") +
  scale_fill_manual(values = alpha(cols,c(1,1))) + labs(title = "(a)") + 
  scale_y_continuous(breaks=seq(20, 100, 20),limits = c(20,100)) +
  theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"))
vp
ggsave(filename = "violin plot for VE.png", plot = vp,width = 4.8, height = 3.4)

# RBD mismatch plot
cols <- c("#FFC0CB","#C0C3FF")
p <- ggplot(df,aes(x=platforms,y=RBD,fill=type)) + 
  geom_violin(width=0.8, color="white",trim = FALSE,position=position_dodge(0.8),scale = "width") + 
  geom_boxplot(width=0.08,outlier.colour = NA,size=0.3,position=position_dodge(0.8)) + 
  stat_summary(fun  = mean,geom="point",position=position_dodge(0.8),colour="white",size=0.7)
vp <- p  + theme_classic() +
  labs(x = "Vaccine platforms", y = "Genetic distance on RBD") +  scale_fill_manual(values = cols) + 
  labs(title = "(b)") +
  scale_y_continuous(breaks=seq(0, 4, 1),limits = c(0,4))  +
  theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"))
vp
ggsave(filename = "violin plot for RBD mismatch.png", plot = vp,width = 4.8, height = 3.4)

# NTD mismatch plot
cols <- c("#FFC0CB","#C0C3FF")
p <- ggplot(df1,aes(x=platforms,y=NTD,fill=type)) + 
  geom_violin(width=0.8, color="white",trim = FALSE,position=position_dodge(0.8),scale = "width") + 
  geom_boxplot(width=0.08,outlier.colour = NA,size=0.3,position=position_dodge(0.8)) + 
  stat_summary(fun  = mean,geom="point",position=position_dodge(0.8),colour="white",size=0.7)
vp <- p  + theme_classic() +
  labs(x = "Vaccine platforms", y = "Genetic distance on NTD") +  scale_fill_manual(values = cols) + 
  labs(title = "(a)") +
  scale_y_continuous(breaks=seq(0, 8, 2),limits = c(0,8))  +
  theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"))
vp
ggsave(filename = "violin plot for NTD mismatch.png", plot = vp,width = 4.8, height = 3.4)

# entire S protein
cols <- c("#FFC0CB","#C0C3FF")
p <- ggplot(df1,aes(x=platforms,y=All,fill=type)) + 
  geom_violin(width=0.8, color="white",trim = FALSE,position=position_dodge(0.8),scale = "width") + 
  geom_boxplot(width=0.08,outlier.colour = NA,size=0.3,position=position_dodge(0.8)) + 
  stat_summary(fun  = mean,geom="point",position=position_dodge(0.8),colour="white",size=0.7)
vp <- p  + theme_classic() +
  labs(x = "Vaccine platforms", y = "Genetic distance on S protein") +  scale_fill_manual(values = cols) + 
  labs(title = "(b)") +
  expand_limits(y = 0) +
  theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"))
vp
ggsave(filename = "violin plot for S protein mismatch.png", plot = vp,width = 4.8, height = 3.4)
