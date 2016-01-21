library("ggplot2")

degust_csv <- read.csv("Bq4hTffh.csv")
ggplot(data= degust_csv, aes(x= avg.expression, y= HM, colour = FDR <0.05)) +
  geom_hline(yintercept = 0, colour ="green", size =0.8)+
  geom_point(size =0.5)+
 
  labs(title= "MA Plot", x= "Average Log Expression", 
       y= "Fold Change in Expression Relative to HM")+
  scale_color_manual(values = c("black", "red"))+

  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.border = element_blank())+
  ggsave("MA_Plot.eps")
