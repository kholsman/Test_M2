# Functions for demo


logit<-function(p) log(p/(1-p))
inv.logit <- function(x) exp(x)/(1+exp(x))


# CEATTLE- Holsman:
S_paijFUN <-function(p,a,Bij,Upaij=U_paij, Bot=Bother,ni=n_i,nj=n_j){
  
  sumA <- 0
  for(i in 1:ni)
    for(j in 1:nj)
      sumA <-  sumA + Upaij[p,a,i,j]/Bij[i,j]
    
    sumB <- 0
    for(i in 1:ni)
      for(j in 1:nj)
        sumB <-  sumB + Upaij[p,a,i,j]
    
    x <- x2 <- 0
    for(i in 1:ni){
      for(j in 1:nj){
        x  <- x  + Upaij[p,a,i,j]/Bij[i,j]
        x2 <- x2 + Upaij[p,a,i,j]
      }
    }
    out <- (Upaij[p,a,,]/Bij)/(sumA+((1-sumB)/Bot))
    return(out)
}


#MSVPA function of Juardo-Molina and Holsman:
M2_Holsman_FUN <- function(i = 3, j = 4, np=n_p,na=n_a, Bij=B_ij,ration=ration_pa,Spaij=S_paij,Bot=Bother){
  M2  <-  0
  for(p in 1:np){
    for(a in 1:na){
      M2 <- M2 +(ration[p,a]*Spaij[p,a,i,j])/((sum(Spaij[p,a,,]*Bij))+(Bot*( 1-sum(Spaij[p,a,,]) ) ) )
    }
  }
  #return(list(M2=M2,test=test))
  return(M2)
}

M2_Curti_FUN <- function(i = 3, j = 4, np=n_p,na=n_a, ration=ration_pa,Bij=B_ij,Cpa=consum_pa,Upaij=U_paij_curti,Bot=Bother){
  M2  <- 0
  for(p in 1:np){
    for(a in 1:na){
      M2  <-  M2 + Upaij[p,a,i,j]*ration[p,a]
    }
  }
  M2 <- M2/Bij[i,j]
  #return(list(M2=M2,test=test))
  return(M2)
}

M2_Adams_FUN <- function(i = 3, j = 4, np=n_p,na=n_a, ration=ration_pa,Bij=B_ij,Cpa=consum_pa,Ppi=p_pi,gpaij=g_paij,Bot=Bother){
  M2  <-  test  <- 0
  for(p in 1:np){
    for(a in 1:na){
      M2  <-  M2 + (Ppi[p,i]*gpaij[p,a,i,j]*ration[p,a])/(sum(Ppi[p,]*gpaij[p,a,,]*Bij)+Bot)
    }
  }
  #return(list(M2=M2,test=test))
  return(M2)
}
M2_Adams_Vother_FUN <- function(i = 3, j = 4, np=n_p,na=n_a, ration=ration_pa,Bij=B_ij,Cpa=consum_pa,Ppi=p_pi,gpaij=g_paij,Vot=V_curti,Bot=Bother){
  M2  <-  test  <- 0
  for(p in 1:np){
    for(a in 1:na){
      M2  <-  M2 + (Ppi[p,i]*gpaij[p,a,i,j]*ration[p,a])/(sum(Ppi[p,]*gpaij[p,a,,]*Bij)+(Vot*Bot))
    }
  }
  #return(list(M2=M2,test=test))
  return(M2)
}
# Fit parms functions:
LL_Ufit<-function(par=c(
  logitPpi =rep(-1.3,n_p),
  logitGaj=rep(-1.3,n_a),
  logsigma=.2),
  data=list(i=1,j=2,UIN=U_paij,Bij=B_ij,Bot=Bother)){
  
  
  Ppi   <- inv.logit(par[1:n_p])
  Gaj   <- inv.logit(par[n_p+(1:n_a)])
  logsigma <- par[length(par)]
  Bij   <- data$Bij
  UIN   <- data$UIN
  Bot   <- data$Bot
  i     <- data$i
  j     <- data$j
  
  Upaij <- UIN[,,i,j]
  Uhat  <- UIN[,,i,j]*0
  
  for(p in 1:n_p){
    for(a in 1:n_a){
      Uhat[p,a] <- (ppi[p]*Gaj[a]*Bij[i,j])/((sum(matrix(ppi[p,],n_i,n_a,byrow = F)*matrix( gaj[a,],n_i,n_a,byrow = T)*Bij))+Bot )
    }
  }
  LL <- -sum(dnorm(Upaij-Uhat,exp(logsigma),log=T))
  
  return(LL)
  
}

