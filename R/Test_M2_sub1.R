#Test_M2_sub1.R

B_ij           <-  matrix(0, n_i,n_j)
M2_ij_sum0     <-  M2_ij  <-  N_y  <-  W_ij  <-  B_ij
S_paij_sum0    <-  S_paij <- U_paij_sum0  <-  S_paij_curti  <- U_paij    <-  array(0,c(n_p,n_a,n_i,n_j))

g_paij         <- u_paij  <- v_paij  <-  v_scaled_paij  <- U_paij*0
O_paij         <- S_paij_curti*0
O_other_curti  <-  S_other_curti <- S_paij_curti[,,1,1]*0
g_aj           <- matrix(0,n_a,n_j)
m2_ij          <- M2_ij*0
n_pi    <- sigma_pi <-p_pi*0

sigma_pi[1,]       <- c(2.259,1.875,1.433)
sigma_pi[2,]       <- c(1.259,2.979,1.093)
n_pi[1,]           <-c(4.159,4.833,3.996)
n_pi[2,]           <-c(2.159,3.95,2.26)


# simulate some vals:
N_y[1,]      <-  100*exp(-.5*(1:n_j))
N_y[2,]      <-  40*exp(-.25*(1:n_j))
N_y[3,]      <-  50*exp(-.7*(1:n_j))

K            <-  c(.22,.45,1.08)
H            <-  c(16.34,9.3, 5.19)
t0           <-  c(.53,-.16,-1.0)
d            <-  c(-.817,-.375,-0.213)
Winf         <-  (H/K)^(1/(1-d))

for(i in 1:n_i)
  for(j in 1:n_j)
    W_ij[i,j]     <- Winf[i]*(1-exp(-K[i]*(1-d[i])*(j-t0[i])))^(1/(1-d[i]))

B_ij         <-  N_y*W_ij
W_pa         <-  W_ij[1:2,]
B_pa         <-  B_ij[1:2,]
N_pa         <-  N_y[1:2,]
Bother       <-  sum(B_ij)*B_other_scalar

#consum_pa    <-  matrix(0, n_p,n_a)
ration_pa    <-  W_pa*0
#consum_pa[1,]     <-  .10*colSums(B_ij)
#consum_pa[2,]     <-  .05*colSums(B_ij)
ration_pa[1,]     <- 365*N_pa[1,]*(0.119*W_pa[1,]^(1-0.460))
ration_pa[2,]     <- 365*N_pa[2,]*0.041*W_pa[2,]^(1-0.122)


# predator pref for a species:

if(setTo1){
  p_pi[1,]  <- rep(1,1)
  p_pi[2,]  <- rep(1,1)
}

if(setCurtiVTo1){
  if(setTo1){
    p_pi[1,]  <- c(1,1,1)
    p_pi[2,]  <- c(1,1,1)
  }
  p_other<- p_pi*0+v_curti
} 
v_other   <- p_other

# size spectrum lit:
# g_paij = exp(dnorm(ln(W_pa/W_ij)-n_pi)^2,sigma_pi)
# v_other (and p_other ) == 1, NOT the scaled v_other eq 6 top of page 3)
# 1/(1+sumV)


# size sprectrum model:
for(p in 1:n_p)
  for(a in 1:n_a)
    for(i in 1:n_i)
      for(j in 1:n_j)
        g_paij[p,a,i,j] <- exp((-1/(2*sigma_pi[p,i]^2))*(log(W_pa[p,a]/W_ij[i,j]) -n_pi[p,i] )^2 )

# for(a in 1:n_a){
#   g_paij[1,a,,] <- matrix(1-inv.logit((1:n_a)-(a/2)),n_i,n_j,byrow=T)
#   g_paij[2,a,,] <- matrix(1-inv.logit((1:n_a)-(a/2)),n_i,n_j,byrow=T)
#   if(a<n_j) {
#     
#     g_paij[1,a,,][log(W_pa[1,a]/W_ij)<0]<-0
#     g_paij[2,a,,][log(W_pa[2,a]/W_ij)<0]<-0
#   }
# }

if(setTo1) {
  for(a in 1:n_a){
      g_paij[1,a,,] <- g_paij[2,a,,] <- 1
      if(a<n_j) {
        
        g_paij[1,a,,][log(W_pa[1,a]/W_ij)<0]<-0
        g_paij[2,a,,][log(W_pa[2,a]/W_ij)<0]<-0
      }
    }
}

# for(a in 1:n_a){
#   g_paij[1,,,] <-(g_paij[1,,,]/apply(g_paij,1,sum)[1])
#   g_paij[2,,,] <-(g_paij[2,,,]/apply(g_paij,1,sum)[2])
# }



for(a in 1:n_a){
  g_aj[a,] <- 1-inv.logit((1:n_a)-(a/2))
}


v_other_scaled<-v_paij[,,1,1]*0

for(p in 1:n_p){
  for(a in 1:n_a){
    v_paij[p,a,,]        <-  p_pi[p,]*g_paij[p,a,,]
    v_scaled_paij[p,a,,] <-  v_paij[p,a,,]/(sum(v_paij[p,a,,])+(v_other[p]))
    v_other_scaled[p,a]  <-  v_other[p]   /(sum(v_paij[p,a,,])+(v_other[p]))
  }
}

v_other_scaled[p,a]+sum(v_scaled_paij[p,a,,])
# versus
v_other[p]+sum(v_paij[p,a,,])

O_paij_wrong <-O_paij*0
# eq 7 : available prey biomass for each predator p,a
for(p in 1:n_p){
  for(a in 1:n_a){
    O_paij[p,a,,]        <- v_scaled_paij[p,a,,]*B_ij
    O_paij_wrong[p,a,,]  <- v_paij[p,a,,]*B_ij
  }
}

# eq 8 : other prey available:
O_other_pa        <- v_other_scaled*Bother
O_other_wrong_pa  <- (v_other_scaled*0+1)*Bother

# eq 10: total prey availble:
O_pa         <-  O_other_pa + apply(O_paij,1:2,sum)
O_pa_wrong   <-  O_other_wrong_pa + apply(O_paij_wrong,1:2,sum)


# eq 11: total prey availble:
U_paij_curti  <-  U_paij_curti_wrong  <-  U_paij_adams  <-  U_paij_adams_rescaleVother  <- O_paij*0
for(p in 1:n_p){
  for(a in 1:n_a){
    U_paij[p,a,,]              <- O_paij[p,a,,]/O_pa[p,a]
    U_paij_curti[p,a,,]        <- O_paij[p,a,,]/O_pa[p,a]
    U_paij_curti_wrong[p,a,,]  <- ration_pa[p,a]*O_paij_wrong[p,a,,]/O_pa_wrong[p,a]
    U_paij_adams[p,a,,]        <- (p_pi[p,]*g_paij[p,a,,]*B_ij)/(sum(p_pi[p,]*g_paij[p,a,,]*B_ij)+Bother)
    v_tmp <- (p_pi[p,]*g_paij[p,a,,])/(sum(p_pi[p,]*g_paij[p,a,,])+v_curti)
    U_paij_adams_rescaleVother[p,a,,]  <- (p_pi[p,]*g_paij[p,a,,]*B_ij)/(sum(p_pi[p,]*g_paij[p,a,,]*B_ij)+v_curti*Bother)
  }
}

S_paij <-  S_paij_adams  <-  v_scaled_paij*0
# Ceattle suitability = scaled suitability from curti:
for(p in 1:n_p){
  for(a in 1:n_a){
    S_paij[p,a,,]           <-  S_paijFUN(p=p,a=a,Bij=B_ij,Upaij=U_paij_curti, Bot=Bother,ni=n_i,nj=n_j)
    S_paij_curti[p,a,,]     <-  S_paijFUN(p=p,a=a,Bij=B_ij,Upaij=U_paij_curti, Bot=Bother,ni=n_i,nj=n_j)
    S_paij_adams[p,a,,]     <-  S_paijFUN(p=p,a=a,Bij=B_ij,Upaij=U_paij_adams, Bot=Bother,ni=n_i,nj=n_j)
  }
}

v_scaled_paij[p,a,,]
S_paij[p,a,,] 

#ration_use <- consum_pa  # Say that consumption rates by each pred were proportional to total available prey (prey limited):
ration_use <- ration_pa  # ration is a function of predator weight


M2_ij_curti <- M2_ij_curti_wrong <-  M2_ij_adams  <-  M2_ij_holsman  <-  M2_ij_holsman_adamsSpaij  <-  M2_ij_adams_fixed  <-  B_ij*0
for(i in 1:n_i){
  for(j in 1:n_j){
    M2_ij_holsman[i,j]                <- M2_Holsman_FUN(    i = i, j = j, np=n_p,na=n_a, Bij=B_ij,ration = ration_use,Spaij=S_paij,Bot=Bother)
    M2_ij_holsman_adamsSpaij[i,j]     <- M2_Holsman_FUN(    i = i, j = j, np=n_p,na=n_a, Bij=B_ij,ration = ration_use,Spaij=S_paij_adams,Bot=Bother)
    
    M2_ij_curti[i,j]       <- M2_Curti_FUN(i = i, j = j, np=n_p,na=n_a, Bij=B_ij,ration = ration_use,Upaij=U_paij_curti,Bot=Bother)
    M2_ij_curti_wrong[i,j] <- M2_Curti_FUN(i = i, j = j, np=n_p,na=n_a, Bij=B_ij,ration = ration_use,Upaij=U_paij_curti_wrong,Bot=Bother)
    M2_ij_adams[i,j]       <- M2_Adams_FUN       (i = i, j = j, np=n_p,na=n_a, Bij=B_ij,ration = ration_use,Ppi=p_pi,gpaij=g_paij,Bot=Bother) 
    M2_ij_adams_fixed[i,j] <- M2_Adams_Vother_FUN(i = i, j = j, np=n_p,na=n_a, Bij=B_ij,ration = ration_use,Ppi=p_pi,gpaij=g_paij,Bot=Bother,Vot=v_curti) 
     }
}
