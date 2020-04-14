

plotM2<-function(M2_list=list(
  Curti   = M2_ij_curti,
  Holsman = M2_ij_holsman,
  Adams   = M2_ij_adams,
  Holsman_AdamsSpaij = M2_ij_holsman_adamsSpaij)){
  
  for(l in 1:length(M2_list)){
    if(l==1){
      mm           <-  melt(M2_list[[l]])
      colnames(mm) <-  c("prey","age","M2")
      mm$prey      <-  factor(paste("prey",mm$prey),levels=paste("prey",1:max(mm$age)))
      mm$lab       <-  names(M2_list)[l]
      out          <-  mm
    }else{
      mm           <-  melt(M2_list[[l]])
      colnames(mm) <-  c("prey","age","M2")
      mm$prey      <-  factor(paste("prey",mm$prey),levels=paste("prey",1:max(mm$age)))
      mm$lab       <-  names(M2_list)[l]
      out          <-  rbind(out,mm)
    }
    out$lab        <-  factor(out$lab,levels=names(M2_list))
  }

  ggplot(data=out,aes(x=age,y=M2)) +
    geom_line(aes(x=age,y=M2,col=prey),size=1)+
    facet_grid(.~lab,scales="free") + 
    theme_minimal() +
    scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
  

}