#' Estimate parameters using the expectation-maximization algorithm
#'
#' \code{EM_miss} Estimate parameters using the expectation-maximization algorithm based on the complete log-likelihood
#'    of observed exposure data in the presence of missing data.
#'
#' @param para_start numeric vectors of starting values for beta and alpha. Beta is the parameter in the linear regression model of
#'    response variable and exposure variable of interest. Alpha is the parameter in the linear regression model of exposure variable
#'    and
#' @param Sigma_start numeric vectors of starting values for Sigma, Sigma is the covariance matrix of the error terms.
#' @param G A numeric matrix of genetic variants
#' @param Y A numeric vector of the outcome variable
#' @param X A numeric vector of the exposure variable
#' @param maxit Maximum number of iterations to solve the estimating equations
#' @param tol Numerical precision
#'
#' @return A list
#' \describe{
#' \item{step}{number of iterations of the estimation algorithm}
#' \item{Sigma}{the final estimates of the covariance matrix.}
#' \item{para}{the final estimates of the parameters organized as (beta and alpha).}
#' }
#' @export
#'
#' @examples
#' library
library(MASS)
EM_miss = function(para_start,Sigma_start,G,Y,X,maxit,tol){
  R = rep(1,length(X))
  R[which(is.na(X)==TRUE)]=0
  X[which(is.na(X)==TRUE)]=0
  G = as.matrix(G)
  K = ncol(G)
  para = para_start
  Sigma = Sigma_start
  #######calculate the expectation
  Expect= function(R,para,Sigma,G,Y,X){
    K = ncol(G)
    beta0 = para[1]
    beta1 = para[2]
    alpha0 = para[3]
    alphav = para[c(4:(4+K-1))]

    CEx_nomi = cbind(-(alpha0+G%*%alphav),Y-beta0)%*%solve(Sigma)%*%c(1,-beta1)
    CEx_domi = as.numeric(t(c(1,-beta1))%*%solve(Sigma)%*%c(1,-beta1))
    CEx2 = (CEx_nomi/CEx_domi)^2 + (CEx_domi)^(-1)
    CEx1 = -(CEx_nomi/CEx_domi)

    Ex1 = R*X+(1-R)*CEx1
    Ex2 = R*X^2+(1-R)*CEx2

    return(cbind(Ex1,Ex2))
  }

  ########some fixed values
  mean_gk = apply(G,2,mean)
  mean_Y = mean(Y)
  c_gg = t(G-rep(1,N)%*%t(mean_gk))%*%(G-rep(1,N)%*%t(mean_gk))/N
  c_ygk = t(Y-mean(Y))%*%(G-rep(1,N)%*%t(mean_gk))/N


  ###########other parameters fixed, just iterate the Sigma matrix
  loglikeEM_para_fixed = function(para,Sigma_start,Y,X,G,R,tol = 1e-8, maxit = 1000){
    K = ncol(G)
    Sigma = Sigma_start
    beta0 = para[1]
    beta1 = para[2]
    alpha0 = para[3]
    alphav = para[c(4:(4+K-1))]
    it = 1
    converged = F
    Sigma_seq = as.vector(Sigma_start)
    while(!converged & it<maxit){

      ############################
      ##########E step############
      ############################
      #######calculate the coefficient matrix of linear equations
      EX = Expect(R,para,Sigma,G,Y,X)
      EX1 = EX[,1]
      EX2 = EX[,2]
      mean_EX1 = mean(EX1)
      mean_EX2 = mean(EX2)
      ##########M step
      Sigma_1 = 1/N* sum((EX2 - (EX1)^2))*c(1,-beta1)%*%t(c(1,-beta1))
      Sigma_2_1 = EX1 - mean_EX1 - (G-rep(1,N)%*%t(mean_gk))%*%alphav
      Sigma_2_2 = Y-mean_Y - as.numeric(beta1)*(EX1 - mean_EX1)
      Sigma_2 = 1/N* (t(cbind(Sigma_2_1, Sigma_2_2))%*% (cbind(Sigma_2_1, Sigma_2_2)))
      Sigma = Sigma_1 + Sigma_2

      Sigma_seq = rbind(Sigma_seq,as.vector(Sigma))

      it = it +1
      if(max(abs(Sigma_seq[it,]-Sigma_seq[it-1,]))<tol)
        converged = TRUE
    }
    return(list(it= it,Sigma = Sigma,Sigma_seq = Sigma_seq))
  }

  ###########Sigma matrix fixed, just iterate the other parameters
  loglikeEM_Sigma_fixed = function(starts,Sigma,Y,X,G,R,tol = 1e-8, maxit = 1000){
    para = starts
    it = 1
    converged = F
    beta0_seq = starts[1]
    beta1_seq = starts[2]
    alpha0_seq = starts[3]
    alphav_seq = starts[4:(4+K-1)]
    while(!converged & it<maxit){

      ############################
      ##########E step############
      ############################
      #######calculate the coefficient matrix of linear equations
      EX = Expect(R,para,Sigma,G,Y,X)
      EX1 = EX[,1]
      EX2 = EX[,2]
      mean_EX1 =  mean(EX1)
      mean_EX2 = mean(EX2)
      c_gEX1 = t(EX1-mean(EX1))%*%(G-rep(1,N)%*%t(mean_gk))/N
      c_yEX1 = t(Y-mean(Y))%*%(EX1-mean(EX1))/N

      CM = matrix(NA,1+K,1+K)
      CM[1,1:K] = solve(Sigma)[2,1]*c_gEX1
      CM[1,K+1] = solve(Sigma)[2,2]*(mean_EX2- (mean_EX1)^2)
      CM[c(2:(K+1)),c(1:K)] = solve(Sigma)[1,1]*c_gg
      CM[c(2:(K+1)),K+1] = solve(Sigma)[1,2]*c_gEX1

      ######the constant vector of the linear equations
      CON = rep(0,K+1)
      CON[1] = solve(Sigma)[2,1]*(mean_EX2 - (mean_EX1)^2)+solve(Sigma)[2,2]*c_yEX1
      CON[2:(K+1)] = as.numeric(solve(Sigma)[1,1]*(c_gEX1)+solve(Sigma)[1,2]*(c_ygk))


      ######################################
      ################M step################
      ######################################

      alphav = as.numeric(solve(CM)[1:K,]%*%CON)
      beta1 = as.numeric(solve(CM)[K+1,]%*%CON)
      beta0 = as.numeric(mean_Y - mean_EX1*beta1)
      alpha0 = as.numeric(mean_EX1 - t(mean_gk)%*%alphav)

      para = c(beta0,beta1,alpha0,alphav)
      alphav_seq = rbind(alphav_seq,alphav)
      beta0_seq = c(beta0_seq,beta0)
      beta1_seq  = c(beta1_seq,beta1)
      alpha0_seq = c(alpha0_seq,alpha0)


      it = it +1
      if(max(max(abs(alphav_seq[it,] - alphav_seq[it-1,])),abs(beta1_seq[it]-beta1_seq[it-1]),abs(beta0_seq[it]-beta0_seq[it-1]),abs(beta1_seq[it]-beta1_seq[it-1]))<tol)
        converged = TRUE
    }
    return(list(it = it,par = para,alphav_seq = alphav_seq,beta0_seq = beta0_seq,beta1_seq = beta1_seq, alpha0_seq=alpha0_seq))
  }

  ##############iterated Sigma first and then iterate other parameters
  Diff = function(x,y) sum((x-y)^2)/sum(x^2+tol)
  #max.step = 1000
  #thres = 1e-8
  diff = tol + 1; step = 0
  #Sigma_start= matrix(c(0.4,0.2,0.2,0.4),2,2)
  #Sigma  = Sigma_start
  #para_start = c(0.6,0.1,rep(0.05,6))
  #para = para_start
  while(diff > tol & step < maxit){
    step = step + 1
    #print(step)

    opt1 = loglikeEM_para_fixed(para,Sigma,Y,X,G,R,tol = 1e-8, maxit = 1000)
    diff1 = Diff(as.vector(opt1$Sigma),as.vector(Sigma))
    Sigma = opt1$Sigma

    opt2 = loglikeEM_Sigma_fixed(para,Sigma,Y,X,G,R,tol = 1e-8, maxit = 1000)
    diff  = max(diff1,Diff(opt2$par,para))
    para = opt2$par

    #  print(list(Sigma = Sigma, par = para))
  }
  return(list(step = step,Sigma = Sigma, para = para))
}
