
# Figure 1

plotC<-function(CIN=ration_pa){
  
  mm           <-  melt(CIN)
  colnames(mm) <-  c("pred","age","consumed")
  mm$pred      <-  factor(paste("pred",mm$pred),levels=paste("pred",1:n_a))
  
   ggplot(data=mm,aes(x=age,y=consumed)) +
    geom_line(aes(x=age,y=consumed,col=pred),size=1)+
    #facet_grid(~prey) + 
    theme_minimal() +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
  
  
}