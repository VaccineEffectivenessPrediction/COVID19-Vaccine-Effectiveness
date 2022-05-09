
vm <- read.csv("variants mismatch input.csv", header = TRUE, fill = TRUE)
library(merTools)
predVE <- predictInterval(p_model,vml)
predVE_vr <- cbind(vm,pred_ve)
predVE_ve_df <- predVE_vr[order(predVE_vr$platforms,predVE_vr$RBD),]
write.csv(predVE_ve_df,"predict variant-specific VE.csv",row.names = FALSE)

ve_df <- read.csv("predict variant-specific VE.csv")
ve_df$platforms <- factor(ve_df$platforms,levels = c("mRNA","Protein subunit","Inactivated","Viral vector"))
cols <- c("#F8766D","#FFB90F","#00BFC4","#C77CFF")
plot <- ggplot(data = ve_df, aes(x = variants,y=fit,fill=platforms)) + 
  geom_bar(position="dodge",stat="identity",width = 0.7) +
  geom_hline(yintercept = 50,size=0.5,linetype="dotted",color="grey36") 

plot_edited <- plot + theme_classic() +
  scale_fill_manual(values=cols,name="Platforms") +
  scale_y_continuous(name="Estimated COVID-19 VE (%)",breaks = seq(-20,100,20),limits = c(-20,100),
                     labels = c("<-20",0,20,40,60,80,100)) +
  xlab("Variants") +
  theme(axis.text.x = element_text(vjust = 0.95,hjust = 0.9, angle = 30,colour = "black"),
        legend.title = element_text(size = 9.5), 
        legend.text = element_text(size = 7), legend.key.size = unit(14, "pt")) 

plot_edited 
ggsave(filename = "predict VE against variants.png", plot = plot_edited,width = 7, height = 4.1)


mm_df <- read.csv("California mismatch input.csv")
pred_ve <- predictInterval(p_model,mm_df)
ve_df <- cbind(mm_df,pred_ve)
write.csv(ve_df,file="California VE predictions.csv",row.names = FALSE)

ve_df <- read.csv("California VE predictions.csv")
ve_df$weeks <- as.Date(ve_df$weeks ,format = "%d/%m/%Y")
ob <- read.csv("observed VE.csv")
ob$start <- as.Date(ob$start ,format = "%d/%m/%Y")
ob$end <- as.Date(ob$end,format = "%d/%m/%Y")
ob$platforms <- toupper(ob$platforms)
id_date <- as.Date("26/11/2021",format = "%d/%m/%Y")
plot <- ggplot() + 
  geom_ribbon(data = ve_df,aes(x=weeks,ymin=lwr,ymax=upr,fill=platforms)) + 
  geom_line(data = ve_df, mapping = aes(x = weeks, y = fit,colour=platforms),size=1) +
  geom_rect(data=ob,aes(xmin=start,xmax=end,ymin=VE-1.5,ymax=VE+1.5,fill=platforms),size=0.3) 

ve_plot <- plot + theme_classic() +
  annotate("text",x=id_date+14,y=105,label="    Omicron\npredominance",hjust=0,vjust=1.5,size=2.8,colour="black") +
  geom_hline(yintercept = 50,size=0.5,linetype="dotted",color="grey36") +
  geom_hline(yintercept = 0,size=0.5,linetype="dotted",color="grey36") +
  geom_vline(xintercept =id_date,size=0.5,linetype="dotted",color="grey36") +  
  scale_color_manual(name="Estimated VE",breaks = c("mRNA","Inactivated","Viral vector","Protein subunit"),
                     values=c('#FF6666','#66CC99',"#6666FF","#FFCC33")) +  
  scale_fill_manual(name="Observed VE",breaks = c("mRNA","Inactivated","Viral vector","Protein subunit","MRNA","VIRAL VECTOR","PROTEIN SUBUNIT"),
                    values=alpha(c('#FFCCCC','#CCFF99',"#99CCFF","#FFFFCC","#FF9999","mediumpurple2","#FFFF00"),c(0.6,0.5,0.5,0.6,0.6,0.5,0.6))) +
  scale_y_continuous(name="COVID-19 VE (%)",limits = c(-26,105),breaks = seq(-20,105,20),labels = c("<-20",0,20,40,60,80,100)) + 
  scale_x_date(date_breaks = "4 weeks",name = "") +
  theme(axis.text.x = element_text(vjust = 0.8,hjust = 0.9, angle = 30,colour = "black"),
        legend.position="none")
ve_plot
ggsave(filename = "California time-varying VE.png", plot = ve_plot,width = 6, height = 4)
