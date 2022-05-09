
# validation plot
df <- read.csv("VE and whole genome mismatch.csv", header = TRUE, fill = TRUE)
df$start <- as.Date(df$start,'%d/%m/%Y')
df$end <- as.Date(df$end,'%d/%m/%Y')
df$type<- factor(df$type,c("efficacy","effectiveness"))
nvr_df <- subset(df,variants == "")
p_model <- lmer(VE ~ RBD+(RBD|platforms) , data=df) # age + TimeAfterSecondDose

library(merTools)
vr_df <- read.csv("validation data input.csv", header = TRUE, fill = TRUE)
vr_df$start <- as.Date(vr_df$start,'%d/%m/%Y')
vr_df$end <- as.Date(vr_df$end,'%d/%m/%Y')
pre_ve_slope <- predictInterval(p_model,vr_df)
vr_out <- cbind(vr_df,pre_ve_slope)
write.csv(vr_out,"validation data output.csv",row.names = FALSE)

df <- read.csv("validation data.csv", header = TRUE, fill = TRUE)
df$id <- factor(df$id,levels=rev(ord))

vplot <- ggplot(df, aes(x= id, y = VE, color = VE_type)) + geom_hline(yintercept = 50,size=0.5,color="grey60",linetype="dashed") +
  geom_point(position = position_dodge(0.4),size=3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.2, position = position_dodge(0.4) )+coord_flip()

vplot_edited <- vplot + scale_color_manual(values=alpha(c("#0099CC","#FF6666"),c(1,1)),  
                               breaks=c("observed","estimated"),
                               labels=c("Observed\n(95% confidence interval)","Estimated\n(95% prediction interval)"),name="VE evaluation in\nindependent studies") +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.key.size = unit(1.6, 'lines'), legend.background = element_blank(),
        legend.key = element_blank(), legend.title = element_text(size = 10), 
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.ticks.y=element_blank(),axis.line.x = element_line(size=0.5,colour = "black") ) +
  scale_y_continuous(name="COVID-19 VE (%)",limits = c(-10,110),breaks=seq(-50, 100, 25),expand = c(0, 0)) +
  scale_x_discrete(name="",labels=label)
vplot_edited
ggsave(filename = "(220318) - validation.png", plot = vplot_edited,width = 6.5, height = 6.5)

# calibration plot
splot <- ggplot(data=df, aes(x = pred_ve, y = VE)) + 
  geom_abline(intercept = 0, slope = 1,linetype="dotted",colour="#009999") +
  geom_point(size=3.6,color="#66CCCC") 

splot_edited <- splot + 
  annotate("text",label="CCC: 0.95 (95% CI: 0.88-0.98)", parse=FALSE, x=-Inf,y=Inf,hjust=-0.05,vjust=2,size=3.5) + 
  scale_y_continuous(name="Observed COVID-19 VE (%)",limits = c(0,100),breaks = seq(0,100,25)) +
  scale_x_continuous(name = "Estimated COVID-19 VE (%)",limits = c(0,100),breaks = seq(0,100,25)) + 
  theme_classic() + 
  theme(legend.position = "none") + coord_fixed(ratio = 1)
print(splot_edited)
ggsave(filename = "(220318) - calibration plot.png", plot = splot_edited,width = 3.5, height = 3.5)
