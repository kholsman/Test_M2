

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
  Adams_Votherincl = M2_ij_adams_Votherincl))

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

B_sim <-  c(seq(0,10,.1),seq(10,100,1))
nsim     <-  length(B_sim)

U_sim <- S_sim  <-  array(0, c(nsim,n_p,n_a,n_i,n_j))
M2_sim <-  M2_sim_holsman  <-  M2_sim_adams  <-  M2_sim_Adams_Votherincl  <- array(0, c(nsim,n_i,n_j))

for(itr in 1:nsim){
  BIN <-B_ij
  #BIN <-BIN*B_sim[itr]
  BotherIN <- Bother*B_sim[itr]
  
  U_sim[itr,,,,]    <-   curti_U_FUN(np    = n_p,
                                     na    = n_a,
                                     ni    = n_i,
                                     nj    = n_j,
                                     ppi   = p_pi,
                                     gpaij = g_paij,
                                     v_ot  = v_other,
                                     Bij   = BIN,
                                     Bot   = BotherIN)
  
  for(p in 1:n_p)
    for(a in 1:n_a)
      S_sim[itr, p,a,,]           <-  S_paijFUN(p=p,a=a,Bij=BIN,Upaij=U_sim[itr,,,,], Bot=BotherIN,ni=n_i,nj=n_j)
  
    for(i in 1:n_i){
      for(j in 1:n_j){
          M2_sim[itr,i,j]         <- M2_Curti_FUN(  i = i, j = j, np=n_p,na=n_a, Bij=BIN,
                                                ration = ration_use,Upaij=U_sim[itr,,,,],Bot=BotherIN)
          M2_sim_holsman[itr,i,j] <- M2_Holsman_FUN(i = i, j = j, np=n_p,na=n_a, Bij=BIN,
                                                    ration = ration_use,Spaij=S_sim[itr,,,,],Bot=BotherIN)
          M2_sim_adams[itr,i,j]        <-  M2_Adams_FUN(i = i, j = j, np=n_p,na=n_a, Bij=BIN,
                                                    ration = ration_use,Ppi=p_pi,gpaij=g_paij,Bot=BotherIN)
          M2_sim_Adams_Votherincl[itr,i,j]   <-  M2_Adams_Vother_FUN(i = i, j = j, np=n_p,na=n_a, Bij=BIN,
                                                   ration = ration_use,Ppi=p_pi,gpaij=g_paij,Bot=BotherIN,Vot=v_curti)
          
      }
    }
}
    
mmm <- melt(M2_sim); colnames(mmm)<-c("itr","prey","age","M2"); mmm$B_other<- Bother*B_sim[mmm$itr]
mmm2 <- melt(M2_sim_holsman); colnames(mmm2)<-c("itr","prey","age","M2"); mmm2$B_other<- Bother*B_sim[mmm2$itr]
mmm3 <- melt(M2_sim_adams); colnames(mmm3)<-c("itr","prey","age","M2"); mmm3$B_other<- Bother*B_sim[mmm3$itr]
mmm4 <- melt(M2_sim_Adams_Votherincl); colnames(mmm4)<-c("itr","prey","age","M2"); mmm4$B_other<- Bother*B_sim[mmm4$itr]

mmm$M2_model <-"Curti"
mmm2$M2_model <-"Holsman"
mmm3$M2_model <-"Adams"
mmm4$M2_model <-"Adams_Votherincl"

mmm$prey<-factor(mmm$prey,levels=1:n_i)
mmm2$prey<-factor(mmm2$prey,levels=1:n_i)
mmm3$prey<-factor(mmm3$prey,levels=1:n_i)
mmm4$prey<-factor(mmm4$prey,levels=1:n_i)

mmm$age<-factor(mmm$age,levels=1:n_a)
mmm2$age<-factor(mmm2$age,levels=1:n_a)
mmm3$age<-factor(mmm3$age,levels=1:n_a)
mmm4$age<-factor(mmm4$age,levels=1:n_a)

dat<-mmm2%>%filter(prey==1)

#plot(B_sim,exp(-dat$M2))
#plot(Bother*B_sim,exp(-dat$M2))
coll <- colorRampPalette( c("#00AFBB", "#E7B800", "#FC4E07"))
ggplot(data=dat,aes(x=B_other,y=exp(-dat$M2))) +
  geom_line(aes(x=B_other,y=exp(-dat$M2),col=age),size=1)+
  facet_grid(.~prey,scales="free") + 
  theme_minimal() +
  scale_color_manual(values = coll(n_a))

#plot(B_sim,exp(-dat$M2))
#plot(Bother*B_sim,exp(-dat$M2))
dat<-rbind(mmm,mmm2,mmm3,mmm4)%>%filter(prey==1)
dat$M2_model<-factor(dat$M2_model,levels=c("Curti","Holsman","Adams","Adams_Votherincl"))

coll <- colorRampPalette( c("#00AFBB", "#E7B800", "#FC4E07"))
# survival
ggplot(data=dat,aes(x=B_other/1e6,y=exp(-dat$M2))) +
  geom_line(aes(x=B_other/1e6,y=exp(-dat$M2),col=age),size=1)+
  facet_grid(prey~M2_model,scales="free") + 
  theme_minimal() +
  scale_color_manual(values = coll(n_a))

# M2
if(1 ==10){
ggplot(data=dat,aes(x=B_other/1e6,y=dat$M2)) +
  geom_line(aes(x=B_other/1e6,y=dat$M2,col=age),size=1)+
  facet_grid(prey~M2_model,scales="free") + 
  theme_minimal() +
  scale_color_manual(values = coll(n_a))
}

## n