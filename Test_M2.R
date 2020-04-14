

rm(list=ls())

library(ggplot2)
library(dplyr)
library(reshape)

setwd("/Users/kholsman/GitHub_new/Test_M2")

# ---------------------------------------------
# Load Functions:
# ---------------------------------------------
source("R/Test_M2_FUN.R")


# ---------------------------------------------
# Step 0: Demo controls
# ---------------------------------------------
n_i          <-  3  # 3 prey species
n_p          <-  2  # 2 predator species
n_a          <-  n_j <- 10  # 10 ages

B_other_scalar  <- 1 # value of 1 means other prey = sum(Bij) was about 118 in curti et al (has effect of scaling the M2)
setTo1          <- F # set all prey i,j pref ( and other to be equal)
setCurtiVTo1    <- F
v_curti         <- 3  # set to 1 as per curti.... not CEATTLE
p_pi      <- matrix(0,n_p,n_i)

if(!setCurtiVTo1){
 
  p_pi[1,]  <- c(.618,1.549,3.255)
  p_pi[2,]  <- c(.7,5.4,3.2)
  p_other   <- p_pi*0+v_curti
  
}else{
  p_pi[1,]  <- c(1,7.17,2.1)
  p_pi[2,]  <- c(10,5,7.6)
  
}



# ---------------------------------------------
# Step 1: Simulate some data:
# ---------------------------------------------

source("R/Test_M2_sub1.R")

# ---------------------------------------------
# Step 2: Plot results
# ---------------------------------------------
source("R/loadplots.R")
graphics.off()
dev.new()
dev.new(height=5,width=7)

# Fig 1: Plot M2
# ---------------------------------------------
plotM2(M2_list=list(
  Curti   = M2_ij_curti,
  Holsman = M2_ij_holsman,
  Adams   = M2_ij_adams,
  Holsman_AdamsSpaij = M2_ij_holsman_adamsSpaij,
  Adams_fix = M2_ij_adams_fixed))

# other plots:
plotB()  # plot B_ij
plotC()  # plot ration

# Fig 2: Plot U_paij
# ---------------------------------------------
p<-1; a<-6
apply(S_paij_adams[p,a,,],1,sum)
apply(S_paij[p,a,,],1,sum)

apply(U_paij_adams[p,a,,],1,sum)
apply(U_paij[p,a,,],1,sum)

par(mfrow=c(1,3))

plot(U_paij_curti[p,a,1,],type="b",ylim=c(0,max(U_paij_curti[p,a,,])))
points(U_paij_curti[p,a,2,],type="b",col="red")
points(U_paij_curti[p,a,3,],type="b",col="blue")

plot(U_paij_adams[p,a,1,],type="b",ylim=c(0,max(U_paij_adams[p,a,,])))
points(U_paij_adams[p,a,2,],type="b",col="red")
points(U_paij_adams[p,a,3,],type="b",col="blue")

plot(U_paij_adams_rescaleVother[p,a,1,],type="b",ylim=c(0,max(U_paij_adams_rescaleVother[p,a,,])))
points(U_paij_adams_rescaleVother[p,a,2,],type="b",col="red")
points(U_paij_adams_rescaleVother[p,a,3,],type="b",col="blue")



# Fig 3: Plot S_paij
# ---------------------------------------------

par(mfrow=c(1,2))

plot(S_paij[p,a,1,],type="b",ylim=c(0,max(S_paij[p,a,,])))
points(S_paij[p,a,2,],type="b",col="red")
points(S_paij[p,a,3,],type="b",col="blue")

plot(S_paij_adams[p,a,1,],type="b",ylim=c(0,max(S_paij_adams[p,a,,])))
points(S_paij_adams[p,a,2,],type="b",col="red")
points(S_paij_adams[p,a,3,],type="b",col="blue")


## Now let's think beyon type II:

# first plot MSVPA/type II function:

Both_sim <-  seq(0,1000,10)
nsim     <-  length(Both_sim)

U_sim <- S_sim  <-  array(0, c(nsim,n_p,n_a,n_i,n_j))
M2_sim <-  M2_sim_holsman  <-  array(0, c(nsim,n_i,n_j))

for(itr in 1:nsim){
  U_sim[itr,,,,]    <-   curti_U_FUN(np    = n_p,
                                     na    = n_a,
                                     ni    = n_i,
                                     nj    = n_j,
                                     ppi   = p_pi,
                                     gpaij = g_paij,
                                     v_ot  = v_other,
                                     Bij   = B_ij,
                                     Bot   = Both_sim[itr])
  
  for(p in 1:n_p)
    for(a in 1:n_a)
      S_sim[itr, p,a,,]           <-  S_paijFUN(p=p,a=a,Bij=B_ij,Upaij=U_sim[itr,,,,], Bot=Both_sim[itr],ni=n_i,nj=n_j)
  
    for(i in 1:n_i){
      for(j in 1:n_j){
          M2_sim[itr,i,j]         <- M2_Curti_FUN(  i = i, j = j, np=n_p,na=n_a, Bij=B_ij,
                                                ration = ration_use,Upaij=U_sim[itr,,,,],Bot=Both_sim[itr])
          M2_sim_holsman[itr,i,j] <- M2_Holsman_FUN(i = i, j = j, np=n_p,na=n_a, Bij=B_ij,
                                                    ration = ration_use,Spaij=S_sim[itr,,,,],Bot=Both_sim[itr])
      }
    }
}
    
mmm <- melt(M2_sim); colnames(mmm)<-c("nitr","prey","age","M2")
mmm2 <- melt(M2_sim_holsman); colnames(mmm2)<-c("nitr","prey","age","M2")
dat<-mmm2%>%filter(prey==1,age==2)
plot(Both_sim,dat$M2)

#now lets try it with some simulated stomach data:
# 
# for(p in 1:n_p){
#   for(a in 1:n_a){
#     U_paij[p,a,1,] <- exp(p*a*.025*(1:n_j)-.03*(1:n_j)^2)
#     U_paij[p,a,2,] <- exp(p*a*.015*(1:n_j)-.03*(1:n_j)^2)
#     U_paij[p,a,3,] <- exp(p*a*.045*(1:n_j)-.03*(1:n_j)^2)
#     U_paij_sum0[p,a,,]  <- U_paij[p,a,,]/sum(U_paij[p,a,,])
#   }
# }
# 

# 
# 
# 
# 
# # set preference to 1 for all ages and species
# p_pi[1,]  <- c(1,1,1)
# p_pi[2,]  <- c(1,1,1)
# 
# p_other   <- 1-apply(p_pi,1,sum)
# p_other   <- 1+0*p_pi
# v_other   <- suit_other *p_other
# 
# for(a in 1:n_a){
#   g_paij[1,a,,] <- g_paij[2,a,,] <- 1
#   if(a<n_j) {
#     g_paij[1,a,,][log(W_pa[1,a]/W_ij)<0]<-0
#     g_paij[2,a,,][log(W_pa[2,a]/W_ij)<0]<-0
#   }
# }
# 
# 
# # now set different preferences for each spp
# p_pi[1,]  <- c(.5,.05,.1)
# p_pi[2,]  <- c(.15,.1,.6)
# 
# 
# #v_other   <-v_other*0+1
# 
# # now set different preferences for each age:
# for(a in 1:n_a){
#   g_paij[1,a,,] <- matrix(1-inv.logit((1:n_a)-(a/4)),n_i,n_j,byrow=T)
#   g_paij[2,a,,] <- matrix(1-inv.logit((1:n_a)-(a/2)),n_i,n_j,byrow=T)
#   if(a<n_j) {
#     
#     g_paij[1,a,,][log(W_pa[1,a]/W_ij)<0]<-0
#     g_paij[2,a,,][log(W_pa[2,a]/W_ij)<0]<-0
#   }
# }
# 
# p_other   <- 1-apply(p_pi,1,sum)
# v_other   <- suit_other *p_other
# 
# 
# U_paij_new <-U_paij*0
# 
# for(p in 1:n_p){
#   for(a in 1:n_a){
#     v_scaled_new <-(p_pi[p,]*g_paij[p,a,,])/(sum(p_pi[p,]*g_paij[p,a,,])+v_other[p])
#     for(i in 1:n_i){
#       U_paij_new[p,a,i,]<-(p_pi[p,i]*g_paij[p,a,i,]*B_ij[i,])/(sum(p_pi[p,]*g_paij[p,a,,]*B_ij)+v_other[p]*Bother)
#     }
#  
# }}
# 
# S_paij_new <- S_paij*0
# # Ceattle suitability = scaled suitability from curti:
# for(p in 1:n_p){
#   for(a in 1:n_a){
#     S_paij_new[p,a,,]  <-  S_paijFUN(p=p,a=a,Bij=B_ij,Upaij=U_paij_new, Bot=Bother,ni=n_i,nj=n_j)
#   }
# }
# 
# M2_ij_holsman_new<-M2_ij_curti_new<-M2_ij_adams_new<-M2_ij_curti
# for(i in 1:n_i){
#   for(j in 1:n_j){
#     M2_ij_holsman_new[i,j]     <- M2_Holsman_FUN(    i = i, j = j, np=n_p,na=n_a, Bij=B_ij,
#                                                  ration = ration_use,Spaij=S_paij_new,Bot=Bother)
#     M2_ij_curti_new[i,j]       <- M2_Curti_FUN(i = i, j = j, np=n_p,na=n_a, Bij=B_ij,
#                                                  ration = ration_use,Upaij=U_paij_new,Bot=Bother) 
#     M2_ij_adams_new[i,j]       <- M2_Adams_FUN(i = i, j = j, np=n_p,na=n_a, Bij=B_ij,
#                                                      ration = ration_use,Ppi=p_pi,gpaij=g_paij,Bot=Bother) 
#     
#   }
# }
# 
# 
# par(mfrow=c(1,3))
# 
# plot(M2_ij_adams_new[1,],type="b",ylim=c(0,max(M2_ij_adams_new)))
# points(M2_ij_adams_new[2,],type="b",col="red")
# points(M2_ij_adams_new[3,],type="b",col="blue")
# 
# plot(M2_ij_curti_new[1,],type="b",ylim=c(0,max(M2_ij_curti_new)))
# points(M2_ij_curti_new[2,],type="b",col="red")
# points(M2_ij_curti_new[3,],type="b",col="blue")
# 
# 
# plot(M2_ij_holsman_new[1,],type="b",ylim=c(0,max(M2_ij_holsman_new)))
# points(M2_ij_holsman_new[2,],type="b",col="red")
# points(M2_ij_holsman_new[3,],type="b",col="blue")
# 
# 
# 
# 
# 
# 
# 
# B_ij<-data.frame(B_ij)
# rownames(B_ij)  <-  paste0("spp_",1:n_i)
# colnames(B_ij)  <-  paste0("age_",1:n_j)
# 
# meta<- melt(t(outer(colnames(B_ij), rownames(B_ij), FUN=paste)))
# colnames(meta) <- c("prey","age","label")
# 
# m            <-  data.frame(melt(B_ij)[,-1])
# colnames(m)  <-  "Biomass"
# m            <-  data.frame(m,meta)
# m$prey       <-  factor(paste("prey",m$prey),levels=paste("prey",1:n_i))
# 
# m            <-  data.frame(melt(B_ij)[,-1])
# colnames(m)  <-  "Biomass"
# m            <-  data.frame(m,meta)
# m$prey       <-  factor(paste("prey",m$prey),levels=paste("prey",1:n_i))
# 
# mm           <-  melt(consum_pa)
# colnames(mm) <-  c("pred","age","consumed")
# mm$pred      <-  factor(paste("pred",mm$pred),levels=paste("pred",1:n_a))
# 
