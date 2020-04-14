

# Figure 1

plotB<-function(B_ijIN=B_ij){
    B_ijIN<-data.frame(B_ijIN)
    meta<- melt(t(outer(colnames(B_ijIN), rownames(B_ijIN), FUN=paste)))
    colnames(meta) <- c("prey","age","label")
  
    m            <-  data.frame(melt(B_ijIN)[,-1])
    colnames(m)  <-  "Biomass"
    m            <-  data.frame(m,meta)
    m$prey       <-  factor(paste("prey",m$prey),levels=paste("prey",1:n_i))
    
    m            <-  data.frame(melt(B_ijIN)[,-1])
    colnames(m)  <-  "Biomass"
    m            <-  data.frame(m,meta)
    m$prey       <-  factor(paste("prey",m$prey),levels=paste("prey",1:n_i))
    
    ggplot(data=m,aes(x=age,y=Biomass)) +
      geom_line(aes(x=age,y=Biomass,col=prey),size=1)+
      #facet_grid(~prey) + 
      theme_minimal() +
      scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

}