EM_miss_var = function(Sigma_est, par_est, G,Y,X){

  R = rep(1,length(X))
  R[which(is.na(X)==TRUE)]=0
  X[which(is.na(X)==TRUE)]=0

  G = as.matrix(G)
  K = ncol(G)
  Omega_est<-solve(Sigma_est) #the estimate of the inverse matrix of Sigma
  beta0_est = par_est[1]
  beta1_est = par_est[2]
  alpha0_est = par_est[3]
  alphav_est = par_est[c(4:(4+K-1))]

  detOmega_est<-det(Omega_est)
  ########some fixed values
  mean_gk = apply(G,2,mean)
  #mean_Y = mean(Y)
  #c_gg = t(G-rep(1,N)%*%t(mean_gk))%*%(G-rep(1,N)%*%t(mean_gk))/N
  #c_ygk = t(Y-mean(Y))%*%(G-rep(1,N)%*%t(mean_gk))/N

  ########Calculate the condition expectations

  av<-alpha0_est+G%*%alphav_est
  bv<-Y-beta0_est

  muv_nume_est<-cbind(-av,bv)%*%Omega_est%*%c(1,-beta1_est)
  muv_deno_est<-as.numeric(t(c(1,-beta1_est))%*%Omega_est%*%c(1,-beta1_est))
  muv_est<--muv_nume_est/muv_deno_est
  sigma_est<-(as.numeric(t(c(1,-beta1_est))%*%Omega_est%*%c(1,-beta1_est)))^{-1/2}

  CEx1v_est<-muv_est
  CEx2v_est<-muv_est^2+sigma_est^2
  CEx3v_est<-muv_est^3+3*muv_est*sigma_est^2
  CEx4v_est<-muv_est^4+6*muv_est^2*sigma_est^2+3*sigma_est^4

  Ex1v_est<-R*X+(1-R)*CEx1v_est
  Ex2v_est<-R*X^2+(1-R)*CEx2v_est
  Ex3v_est<-R*X^3+(1-R)*CEx3v_est
  Ex4v_est<-R*X^4+(1-R)*CEx4v_est


  ########Calculate the D20Q matrix

  D20Q<-matrix(0,6+K,6+K)
  D20Q[1,1]<--N*Omega_est[2,2]^2/(2*detOmega_est^2)
  D20Q[1,2]<-N*Omega_est[1,2]*Omega_est[2,2]/detOmega_est^2
  D20Q[1,3]<--N*Omega_est[1,2]^2/(2*detOmega_est^2)
  D20Q[1,4]<-0
  D20Q[1,5]<-0
  D20Q[1,6]<-sum(Ex1v_est-av)
  D20Q[1,7:(6+K)]<-t(Ex1v_est-av)%*%G
  D20Q[2,2]<--N*(Omega_est[1,1]*Omega_est[2,2]+Omega_est[1,2]^2)/detOmega_est^2
  D20Q[2,3]<-N*Omega_est[1,1]*Omega_est[1,2]/detOmega_est^2
  D20Q[2,4]<-sum(Ex1v_est-av)
  D20Q[2,5]<-sum(Ex2v_est)-sum(av*Ex1v_est)
  D20Q[2,6]<-sum(bv-beta1_est*Ex1v_est)
  D20Q[2,7:(6+K)]<-t(bv-beta1_est*Ex1v_est)%*%G
  D20Q[3,3]<--N*Omega_est[1,1]^2/(2*detOmega_est^2)
  D20Q[3,4]<-sum(bv-beta1_est*Ex1v_est)
  D20Q[3,5]<--beta1_est*sum(Ex2v_est)+sum(bv*Ex1v_est)
  D20Q[3,6]<-0
  D20Q[3,7:(6+K)]<-0
  D20Q[4,4]<--N*Omega_est[2,2]
  D20Q[4,5]<--Omega_est[2,2]*sum(Ex1v_est)
  D20Q[4,6]<--N*Omega_est[1,2]
  D20Q[4,7:(6+K)]<--Omega_est[1,2]*N*mean_gk
  D20Q[5,5]<--Omega_est[2,2]*sum(Ex2v_est)
  D20Q[5,6]<--Omega_est[1,2]*sum(Ex1v_est)
  D20Q[5,7:(6+K)]<--Omega_est[1,2]*(t(Ex1v_est)%*%G)
  D20Q[6,6]<--N*Omega_est[1,1]
  D20Q[6,7:(6+K)]<--Omega_est[1,1]*N*mean_gk

  D20Q<-D20Q+t(D20Q)-diag(diag(D20Q))
  D20Q[7:(6+K),7:(6+K)]<--Omega_est[1,1]*t(G)%*%G

  ########Calculate the H matrix

  Ex4vmEx2vsq_est<-Ex4v_est-Ex2v_est^2
  Ex3vmEx2vEx1v_est<-Ex3v_est-Ex2v_est*Ex1v_est
  Ex2vmEx1vsq_est<-Ex2v_est-Ex1v_est^2
  bvpbeta1av<-bv+beta1_est*av
  Ome22bvmOme21av<-Omega_est[2,2]*bv-Omega_est[2,1]*av
  Ome12bvmOme11av<-Omega_est[1,2]*bv-Omega_est[1,1]*av

  sumEx1_est<-sum(Ex1v_est)
  sumEx2_est<-sum(Ex2v_est)
  sumasq<-sum(av^2)
  sumab<-sum(av*bv)
  sumbsq<-sum(bv^2)
  sumaEx1_est<-sum(av*Ex1v_est)
  sumbEx1_est<-sum(bv*Ex1v_est)

  H<-matrix(0,6+K,6+K)
  H[1,1]<-(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)/4+sum(av^2*Ex2vmEx1vsq_est)+sumaEx1_est^2-sum(av*Ex3vmEx2vEx1v_est)-sumEx2_est*sumaEx1_est
  H[1,1]<-H[1,1]+(N*Omega_est[2,2]/detOmega_est-sumasq)*(sumaEx1_est-sumEx2_est/2)+(N*Omega_est[2,2]/detOmega_est-sumasq)^2/4

  H[1,2]<--beta1_est*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)/2+(sum(bvpbeta1av*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(bvpbeta1av*Ex1v_est))/2
  H[1,2]<-H[1,2]+beta1_est*(sum(av*Ex3vmEx2vEx1v_est)+sumEx2_est*sumaEx1_est)
  H[1,2]<-H[1,2]-(sum(av*bvpbeta1av*Ex2vmEx1vsq_est)+sumaEx1_est*sum(bvpbeta1av*Ex1v_est))
  H[1,2]<-H[1,2]+(N*Omega_est[2,2]/detOmega_est-sumasq)*(beta1_est*sumEx2_est-sum(bvpbeta1av*Ex1v_est))/2
  H[1,2]<-H[1,2]+(N*Omega_est[2,1]/detOmega_est-sumab)*(sumEx2_est/2-sumaEx1_est)
  H[1,2]<-H[1,2]-(N*Omega_est[2,2]/detOmega_est-sumasq)*(N*Omega_est[2,1]/detOmega_est-sumab)/2

  H[1,3]<-beta1_est^2*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)/4
  H[1,3]<-H[1,3]-beta1_est*(sum(bv*Ex3vmEx2vEx1v_est)+sumEx2_est*sumbEx1_est)/2
  H[1,3]<-H[1,3]-beta1_est^2*(sum(av*Ex3vmEx2vEx1v_est)+sumEx2_est*sumaEx1_est)/2
  H[1,3]<-H[1,3]+beta1_est*(sum(av*bv*Ex2vmEx1vsq_est)+sumaEx1_est*sumbEx1_est)
  H[1,3]<-H[1,3]-beta1_est*(N*Omega_est[2,2]/detOmega_est-sumasq)*(beta1_est*sumEx2_est/4-sumbEx1_est/2)
  H[1,3]<-H[1,3]-(N*Omega_est[1,1]/detOmega_est-sumbsq)*(sumEx2_est/4-sumaEx1_est/2)
  H[1,3]<-H[1,3]+(N*Omega_est[2,2]/detOmega_est-sumasq)*(N*Omega_est[1,1]/detOmega_est-sumbsq)/4

  H[1,4]<--(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(Ex3vmEx2vEx1v_est)+sumEx2_est*sumEx1_est)/2
  H[1,4]<-H[1,4]+(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(av*Ex2vmEx1vsq_est)+sumaEx1_est*sumEx1_est)
  H[1,4]<-H[1,4]-sum(Ome22bvmOme21av)*sumEx2_est/2+sum(Ome22bvmOme21av)*sumaEx1_est
  H[1,4]<-H[1,4]+(N*Omega_est[2,2]/detOmega_est-sumasq)*((Omega_est[2,1]-beta1_est*Omega_est[2,2])*sumEx1_est+sum(Ome22bvmOme21av))/2

  H[1,5]<--(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)/2
  H[1,5]<-H[1,5]-(sum(Ome22bvmOme21av*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(Ome22bvmOme21av*Ex1v_est))/2
  H[1,5]<-H[1,5]+(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(av*Ex3vmEx2vEx1v_est)+sumEx2_est*sumaEx1_est)
  H[1,5]<-H[1,5]+sum(av*Ome22bvmOme21av*Ex2vmEx1vsq_est)+sumaEx1_est*sum(Ome22bvmOme21av*Ex1v_est)
  H[1,5]<-H[1,5]+(N*Omega_est[2,2]/detOmega_est-sumasq)*((Omega_est[2,1]-beta1_est*Omega_est[2,2])*sumEx2_est+sum(Ome22bvmOme21av*Ex1v_est))/2

  H[1,6]<--(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(Ex3vmEx2vEx1v_est)+sumEx2_est*sumEx1_est)/2
  H[1,6]<-H[1,6]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(av*Ex2vmEx1vsq_est)+sumaEx1_est*sumEx1_est)
  H[1,6]<-H[1,6]-sum(Ome12bvmOme11av)*sumEx2_est/2+sum(Ome12bvmOme11av)*sumaEx1_est
  H[1,6]<-H[1,6]+(N*Omega_est[2,2]/detOmega_est-sumasq)*((Omega_est[1,1]-beta1_est*Omega_est[1,2])*sumEx1_est+sum(Ome12bvmOme11av))/2

  for(i in 1:K){
    H[1,6+i]<--(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(G[,i]*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(G[,i]*Ex1v_est))/2
    H[1,6+i]<-H[1,6+i]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(av*G[,i]*Ex2vmEx1vsq_est)+sumaEx1_est*sum(G[,i]*Ex1v_est))
    H[1,6+i]<-H[1,6+i]-sum(Ome12bvmOme11av*G[,i])*(sumEx2_est/2-sumaEx1_est)
    H[1,6+i]<-H[1,6+i]+(N*Omega_est[2,2]/detOmega_est-sumasq)*((Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(G[,i]*Ex1v_est)+sum(Ome12bvmOme11av*G[,i]))/2
  }

  H[2,2]<-beta1_est^2*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)
  H[2,2]<-H[2,2]+(sum(bvpbeta1av^2*Ex2vmEx1vsq_est)+(sum(bvpbeta1av*Ex1v_est))^2)
  H[2,2]<-H[2,2]-2*beta1_est*(sum(bvpbeta1av*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(bvpbeta1av*Ex1v_est))
  H[2,2]<-H[2,2]-2*(N*Omega_est[2,1]/detOmega_est-sumab)*(beta1_est*sumEx2_est-sum(bvpbeta1av*Ex1v_est))
  H[2,2]<-H[2,2]+(N*Omega_est[2,1]/detOmega_est-sumab)^2

  H[2,3]<--beta1_est^3*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)/2
  H[2,3]<-H[2,3]+beta1_est^2*(sum(bv*Ex3vmEx2vEx1v_est)+sumEx2_est*sumbEx1_est)
  H[2,3]<-H[2,3]+beta1_est^2*(sum(bvpbeta1av*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(bvpbeta1av*Ex1v_est))/2
  H[2,3]<-H[2,3]-beta1_est*(sum(bv*bvpbeta1av*Ex2vmEx1vsq_est)+sum(bvpbeta1av*Ex1v_est)*sumbEx1_est)
  H[2,3]<-H[2,3]+(N*Omega_est[1,1]/detOmega_est-sumbsq)*(beta1_est*sumEx2_est-sum(bvpbeta1av*Ex1v_est))/2
  H[2,3]<-H[2,3]+beta1_est*(N*Omega_est[2,1]/detOmega_est-sumab)*(beta1_est*sumEx2_est/2-sum(bv*Ex1v_est))
  H[2,3]<-H[2,3]-(N*Omega_est[2,1]/detOmega_est-sumab)*(N*Omega_est[1,1]/detOmega_est-sumbsq)/2

  H[2,4]<-beta1_est*(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(Ex3vmEx2vEx1v_est)+sumEx2_est*sumEx1_est)
  H[2,4]<-H[2,4]-(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(bvpbeta1av*Ex2vmEx1vsq_est)+sumEx1_est*sum(bvpbeta1av*Ex1v_est))
  H[2,4]<-H[2,4]+sum(Ome22bvmOme21av)*(beta1_est*sumEx2_est-sum(bvpbeta1av*Ex1v_est))
  H[2,4]<-H[2,4]-(N*Omega_est[2,1]/detOmega_est-sumab)*((Omega_est[2,1]-beta1_est*Omega_est[2,2])*sumEx1_est+sum(Ome22bvmOme21av))

  H[2,5]<-beta1_est*(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)
  H[2,5]<-H[2,5]+beta1_est*(sum(Ome22bvmOme21av*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(Ome22bvmOme21av*Ex1v_est))
  H[2,5]<-H[2,5]-(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(bvpbeta1av*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(bvpbeta1av*Ex1v_est))
  H[2,5]<-H[2,5]-(sum(bvpbeta1av*Ome22bvmOme21av*Ex2vmEx1vsq_est)+sum(bvpbeta1av*Ex1v_est)*sum(Ome22bvmOme21av*Ex1v_est))
  H[2,5]<-H[2,5]-(N*Omega_est[2,1]/detOmega_est-sumab)*((Omega_est[2,1]-beta1_est*Omega_est[2,2])*sumEx2_est+sum(Ome22bvmOme21av*Ex1v_est))

  H[2,6]<-beta1_est*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(Ex3vmEx2vEx1v_est)+sumEx2_est*sumEx1_est)
  H[2,6]<-H[2,6]-(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(bvpbeta1av*Ex2vmEx1vsq_est)+sumEx1_est*sum(bvpbeta1av*Ex1v_est))
  H[2,6]<-H[2,6]+sum(Ome12bvmOme11av)*(beta1_est*sumEx2_est-sum(bvpbeta1av*Ex1v_est))
  H[2,6]<-H[2,6]-(N*Omega_est[2,1]/detOmega_est-sumab)*((Omega_est[1,1]-beta1_est*Omega_est[1,2])*sumEx1_est+sum(Ome12bvmOme11av))

  for(i in 1:K){
    H[2,6+i]<-beta1_est*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(G[,i]*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(G[,i]*Ex1v_est))
    H[2,6+i]<-H[2,6+i]-(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(bvpbeta1av*G[,i]*Ex2vmEx1vsq_est)+sum(G[,i]*Ex1v_est)*sum(bvpbeta1av*Ex1v_est))
    H[2,6+i]<-H[2,6+i]+sum(Ome12bvmOme11av*G[,i])*(beta1_est*sumEx2_est-sum(bvpbeta1av*Ex1v_est))
    H[2,6+i]<-H[2,6+i]-(N*Omega_est[2,1]/detOmega_est-sumab)*((Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(G[,i]*Ex1v_est)+sum(Ome12bvmOme11av*G[,i]))
  }

  H[3,3]<-beta1_est^4*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)/4
  H[3,3]<-H[3,3]-beta1_est^3*(sum(bv*Ex3vmEx2vEx1v_est)+sumEx2_est*sumbEx1_est)
  H[3,3]<-H[3,3]+beta1_est^2*(sum(bv^2*Ex2vmEx1vsq_est)+sumbEx1_est^2)
  H[3,3]<-H[3,3]-beta1_est*(N*Omega_est[1,1]/detOmega_est-sumbsq)*(beta1_est*sumEx2_est/2-sumbEx1_est)
  H[3,3]<-H[3,3]+(N*Omega_est[1,1]/detOmega_est-sumbsq)^2/4

  H[3,4]<--beta1_est^2*(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(Ex3vmEx2vEx1v_est)+sumEx2_est*sumEx1_est)/2
  H[3,4]<-H[3,4]+beta1_est*(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(bv*Ex2vmEx1vsq_est)+sumbEx1_est*sumEx1_est)
  H[3,4]<-H[3,4]-beta1_est*sum(Ome22bvmOme21av)*(beta1_est*sumEx2_est/2-sumbEx1_est)
  H[3,4]<-H[3,4]+(N*Omega_est[1,1]/detOmega_est-sumbsq)*((Omega_est[2,1]-beta1_est*Omega_est[2,2])*sumEx1_est+sum(Ome22bvmOme21av))/2

  H[3,5]<--beta1_est^2*(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)/2
  H[3,5]<-H[3,5]-beta1_est^2*(sum(Ome22bvmOme21av*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(Ome22bvmOme21av*Ex1v_est))/2
  H[3,5]<-H[3,5]+beta1_est*(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(bv*Ex3vmEx2vEx1v_est)+sumEx2_est*sumbEx1_est)
  H[3,5]<-H[3,5]+beta1_est*(sum(bv*Ome22bvmOme21av*Ex2vmEx1vsq_est)+sumbEx1_est*sum(Ome22bvmOme21av*Ex1v_est))
  H[3,5]<-H[3,5]+(N*Omega_est[1,1]/detOmega_est-sumbsq)*((Omega_est[2,1]-beta1_est*Omega_est[2,2])*sumEx2_est+sum(Ome22bvmOme21av*Ex1v_est))/2

  H[3,6]<--beta1_est^2*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(Ex3vmEx2vEx1v_est)+sumEx2_est*sumEx1_est)/2
  H[3,6]<-H[3,6]+beta1_est*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(bv*Ex2vmEx1vsq_est)+sumbEx1_est*sumEx1_est)
  H[3,6]<-H[3,6]-beta1_est*sum(Ome12bvmOme11av)*(beta1_est*sumEx2_est/2-sumbEx1_est)
  H[3,6]<-H[3,6]+(N*Omega_est[1,1]/detOmega_est-sumbsq)*((Omega_est[1,1]-beta1_est*Omega_est[1,2])*sumEx1_est+sum(Ome12bvmOme11av))/2

  for(i in 1:K){
    H[3,6+i]<--beta1_est^2*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(G[,i]*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(G[,i]*Ex1v_est))/2
    H[3,6+i]<-H[3,6+i]+beta1_est*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(bv*G[,i]*Ex2vmEx1vsq_est)+sumbEx1_est*sum(G[,i]*Ex1v_est))
    H[3,6+i]<-H[3,6+i]-beta1_est*sum(Ome12bvmOme11av*G[,i])*(beta1_est*sumEx2_est/2-sumbEx1_est)
    H[3,6+i]<-H[3,6+i]+(N*Omega_est[1,1]/detOmega_est-sumbsq)*((Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(G[,i]*Ex1v_est)+sum(Ome12bvmOme11av*G[,i]))/2
  }

  H[4,4]<-(Omega_est[2,1]-beta1_est*Omega_est[2,2])^2*(sum(Ex2vmEx1vsq_est)+sumEx1_est^2)
  H[4,4]<-H[4,4]+2*(Omega_est[2,1]-beta1_est*Omega_est[2,2])*sum(Ome22bvmOme21av)*sumEx1_est
  H[4,4]<-H[4,4]+(sum(Ome22bvmOme21av))^2

  H[4,5]<-(Omega_est[2,1]-beta1_est*Omega_est[2,2])^2*(sum(Ex3vmEx2vEx1v_est)+sumEx2_est*sumEx1_est)
  H[4,5]<-H[4,5]+(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(Ome22bvmOme21av*Ex2vmEx1vsq_est)+sumEx1_est*sum(Ome22bvmOme21av*Ex1v_est))
  H[4,5]<-H[4,5]+sum(Ome22bvmOme21av)*((Omega_est[2,1]-beta1_est*Omega_est[2,2])*sumEx2_est+sum(Ome22bvmOme21av*Ex1v_est))

  H[4,6]<-(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(Ex2vmEx1vsq_est)+sumEx1_est^2)
  H[4,6]<-H[4,6]+(Omega_est[2,1]-beta1_est*Omega_est[2,2])*sum(Ome12bvmOme11av)*sumEx1_est
  H[4,6]<-H[4,6]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(Ome22bvmOme21av)*sumEx1_est
  H[4,6]<-H[4,6]+sum(Ome22bvmOme21av)*sum(Ome12bvmOme11av)

  for(i in 1:K){
    H[4,6+i]<-(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(G[,i]*Ex2vmEx1vsq_est)+sumEx1_est*sum(G[,i]*Ex1v_est))
    H[4,6+i]<-H[4,6+i]+(Omega_est[2,1]-beta1_est*Omega_est[2,2])*sum(Ome12bvmOme11av*G[,i])*sumEx1_est
    H[4,6+i]<-H[4,6+i]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(Ome22bvmOme21av)*sum(G[,i]*Ex1v_est)
    H[4,6+i]<-H[4,6+i]+sum(Ome22bvmOme21av)*sum(Ome12bvmOme11av*G[,i])
  }

  H[5,5]<-(Omega_est[2,1]-beta1_est*Omega_est[2,2])^2*(sum(Ex4vmEx2vsq_est)+sumEx2_est^2)
  H[5,5]<-H[5,5]+2*(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(sum(Ome22bvmOme21av*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(Ome22bvmOme21av*Ex1v_est))
  H[5,5]<-H[5,5]+sum(Ome22bvmOme21av^2*Ex2vmEx1vsq_est)+(sum(Ome22bvmOme21av*Ex1v_est))^2

  H[5,6]<-(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(Ex3vmEx2vEx1v_est)+sumEx2_est*sumEx1_est)
  H[5,6]<-H[5,6]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(Ome22bvmOme21av*Ex2vmEx1vsq_est)+sumEx1_est*sum(Ome22bvmOme21av*Ex1v_est))
  H[5,6]<-H[5,6]+(Omega_est[2,1]-beta1_est*Omega_est[2,2])*sum(Ome12bvmOme11av)*sumEx2_est
  H[5,6]<-H[5,6]+sum(Ome12bvmOme11av)*sum(Ome22bvmOme21av*Ex1v_est)

  for(i in 1:K){
    H[5,6+i]<-(Omega_est[2,1]-beta1_est*Omega_est[2,2])*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(G[,i]*Ex3vmEx2vEx1v_est)+sumEx2_est*sum(G[,i]*Ex1v_est))
    H[5,6+i]<-H[5,6+i]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*(sum(Ome22bvmOme21av*G[,i]*Ex2vmEx1vsq_est)+sum(G[,i]*Ex1v_est)*sum(Ome22bvmOme21av*Ex1v_est))
    H[5,6+i]<-H[5,6+i]+(Omega_est[2,1]-beta1_est*Omega_est[2,2])*sum(Ome12bvmOme11av*G[,i])*sumEx2_est
    H[5,6+i]<-H[5,6+i]+sum(Ome12bvmOme11av*G[,i])*sum(Ome22bvmOme21av*Ex1v_est)
  }

  H[6,6]<-(Omega_est[1,1]-beta1_est*Omega_est[1,2])^2*(sum(Ex2vmEx1vsq_est)+sumEx1_est^2)
  H[6,6]<-H[6,6]+2*(Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(Ome12bvmOme11av)*sumEx1_est
  H[6,6]<-H[6,6]+(sum(Ome12bvmOme11av))^2

  for(i in 1:K){
    H[6,6+i]<-(Omega_est[1,1]-beta1_est*Omega_est[1,2])^2*(sum(G[,i]*Ex2vmEx1vsq_est)+sumEx1_est*sum(G[,i]*Ex1v_est))
    H[6,6+i]<-H[6,6+i]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(Ome12bvmOme11av*G[,i])*sumEx1_est
    H[6,6+i]<-H[6,6+i]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(Ome12bvmOme11av)*sum(G[,i]*Ex1v_est)
    H[6,6+i]<-H[6,6+i]+sum(Ome12bvmOme11av)*sum(Ome12bvmOme11av*G[,i])
  }

  for(i in 1:K){
    for(j in i:K){
      H[6+i,6+j]<-(Omega_est[1,1]-beta1_est*Omega_est[1,2])^2*(sum(G[,i]*G[,j]*Ex2vmEx1vsq_est)+sum(G[,i]*Ex1v_est)*sum(G[,j]*Ex1v_est))
      H[6+i,6+j]<-H[6+i,6+j]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(Ome12bvmOme11av*G[,i])*sum(G[,j]*Ex1v_est)
      H[6+i,6+j]<-H[6+i,6+j]+(Omega_est[1,1]-beta1_est*Omega_est[1,2])*sum(Ome12bvmOme11av*G[,j])*sum(G[,i]*Ex1v_est)
      H[6+i,6+j]<-H[6+i,6+j]+sum(Ome12bvmOme11av*G[,i])*sum(Ome12bvmOme11av*G[,j])
    }
  }
  H<-H+t(H)-diag(diag(H))

  ########Calculate the information matrix

  Infmatrix<--D20Q-H

  ########estimate the covariance matrix of par_est

  Covmatrixpara_est<-solve(Infmatrix)

  ########the vector of the variance of parameter

  varparav_est<-diag(Covmatrixpara_est)

  se  = sqrt(varparav_est)

  ########the SE of causal effect beta_1
  #sebeta1<-sqrt(varparav_est[5])
  return(list(cov_matrix = Covmatrixpara_est, se = se))
}
