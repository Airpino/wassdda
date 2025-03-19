library(progress)


WH.mat.prod_discr<-function(T1,T2,transpose1=F,transpose2=F){# Product of distributional matrices
  T1<-tibble(T1)
  T2<-tibble(T2)
  if(transpose1){
    if (nrow(T1)==1){
      T1<-TransposeTib_C(T1)
      tmp<-list()
      tmp[[1]]<-T1
      T1<-tmp
    }else{
      T1<-TransposeTib_C(T1)}
  }
  if(transpose2){
    if (nrow(T2)==1){
      T2<-TransposeTib_C(T2)
      tmp<-list()
      tmp[[1]]<-T2
      T2<-tmp
    }else{
      T2<-TransposeTib_C(T2)}
  }
  
  r1<-length(T1[[1]])
  c1<-length(T1)
  r2<-length(T2[[1]])
  c2<-length(T2)
  res<-matrix(0,r1,c2)
  for(i1 in 1:r1){
    for (i2 in 1:c2){
      for(j in 1:c1){
        #get common cdf
        comm_cdf<-c(0,sort(unique(c(T1[[j]][[i1]]$cdf,T2[[i2]][[j]]$cdf))))
        d1_dom<-extract_quantiles(T1[[j]][[i1]],p=comm_cdf)
        d2_dom<-extract_quantiles(T2[[i2]][[j]],p=comm_cdf)
        res[i1,i2]<-res[i1,i2]+sum(d1_dom*d2_dom*c(0,diff(comm_cdf)))
      }
    }
  }
  return(res)
}

WH.mat.prod_discr_for_boot<-function(T1,T2,transpose1=F,transpose2=F){
  # Product of distributional matrices for bootstrap
  T1<-tibble(T1)
  T2<-tibble(T2)
  if(transpose1){
    if (nrow(T1)==1){
      T1<-TransposeTib_C(T1)
      tmp<-list()
      tmp[[1]]<-T1
      T1<-tmp
    }else{
      T1<-TransposeTib_C(T1)}
  }
  if(transpose2){
    if (nrow(T2)==1){
      T2<-TransposeTib_C(T2)
      tmp<-list()
      tmp[[1]]<-T2
      T2<-tmp
    }else{
      T2<-TransposeTib_C(T2)}
  }
  
  r1<-length(T1[[1]])
  c1<-length(T1)
  r2<-length(T2[[1]])
  c2<-length(T2)
  res<-array(0,dim=c(r1,c2,c1))
  for(i1 in 1:r1){
    for (i2 in 1:c2){
      for(j in 1:c1){
        #get common cdf
        comm_cdf<-c(0,sort(unique(c(T1[[j]][[i1]]$cdf,T2[[i2]][[j]]$cdf))))
        d1_dom<-extract_quantiles(T1[[j]][[i1]],p=comm_cdf)
        d2_dom<-extract_quantiles(T2[[i2]][[j]],p=comm_cdf)
        res[i1,i2,j]<-sum(d1_dom*d2_dom*c(0,diff(comm_cdf)))
      }
    }
  }
  return(res)
}

#regression_discrete
WH.regression.two.components_discr <- function(Tib, #the data table
                                               Yvar, # variable to predict
                                               Xvars, # 
                                               mY,mX){#, simplify = FALSE, qua = 20) {
  selected <- c(1:nrow(Tib))
  
  if(is_tibble(Tib)){Y <- Tib[selected, Yvar]}else{ Y <- tibble(Y=Tib[selected, Yvar])}
  
  
  names(Y)<-names(Tib)[Yvar]
  if(is_tibble(Tib)){
    X <- Tib[selected, Xvars]}else{
      X <- Tib[selected, Xvars]
      if(length(Xvars)>1){
        X <- Tib[selected, Xvars]}else{
          X <- tibble(X=Tib[selected, Xvars])
        }
    }
  
  
  n <- length(selected)
  d <- ncol(X)
  
  # extract means and do multiple regression on means
  MatAver <- cbind(mY,mX) # the matrix of means of the distributions
  colnames(MatAver) <- paste0("AV_", c(colnames(Y), colnames(X)))
  rownames(MatAver) <- rownames(Y)
  MatAver <- as.data.frame(MatAver)
  xnam <- colnames(MatAver)[2:(d + 1)]
  fmla <- as.formula(paste(colnames(MatAver)[1], paste(xnam, collapse = "+"), sep = "~"))
  fit <- lm(fmla, data = MatAver)
  AveCoeff <- fit$coefficient
  names(AveCoeff)[1] <- "(AV_Intercept)"
  
  # center data
  YC<-Y
  XC<-X
  
  for (i in 1:n) {
    YC[[1]][[i]][,1] <- Y[[1]][[i]][,1] - MatAver[i,1]
    for (j in 1:d) {
      XC[[j]][[i]][,1] <- X[[j]][[i]][,1] - MatAver[i,j+1]
    }
  }
  
  
  gammas <- WH.NNLS_discr(XC, YC)
  # print(gammas2)
  # gammas <- WH.NNLS_discr(XC, YC)
  # print(gammas)
  gammas <- as.vector(gammas)
  names(gammas) <- paste0("CEN_", colnames(X))
  return(parameters = c(AveCoeff, gammas))
}
## NNLS for histogram variables -----

WH.NNLS_discr<- function(X, Y) {
  # Non Negative Least Squares for histogram variables modifying
  # Lawson, Charles L.; Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
  tol <- 1e-8
  # Lawson and Hanson NNLS algorithm
  # step 1
  P <- integer(0)
  Z <- c(1:ncol(X))
  gamma <- matrix(0, nrow = ncol(X), ncol = 1)#x
  stp <- 0
  
  XpX<- WH.mat.prod_discr(X, X, transpose1 = TRUE)
  XpY<- WH.mat.prod_discr(X, Y, transpose1 = TRUE)
  w <- XpY - XpX %*% gamma
  while((length(Z)>0)&(max(w)>tol)){
    
    j<-which.max(w)
    P<-sort(unique(c(P,j)))
    # print(P)
    # print(w[,1])
    Z<-Z[Z!=j]#R
    s<-matrix(0,ncol(X),1)
    GG<-solve(XpX[P,P]) %*%XpY[P,]#S^p
    s[P,1]<-GG
    while(min(GG)<=0){
      tmp<-gamma[P,1]/(gamma[P,1]-s[P,1])
      alpha<-min(tmp[GG<=0,1])
      gamma<-gamma+alpha*(s-gamma)
      
      Z<-which(gamma<=0)
      P<-c(1:ncol(X))[-c(1:ncol(X))%in%Z]
      # print(P)
      s<-matrix(0,ncol(X),1)
      GG<-solve(XpX[P,P]) %*%XpY[P,]#S^p
      s[P,1]<-GG
    }
    gamma<-s
    w <- XpY - XpX %*% gamma
  }
  return(gamma)
  
}

### PREDICTING ------



WH.two_comp_predict<-function(Tib,Xvars,mX,coeffs){
  Tib<-as_tibble(Tib)
  mX<-t(t(mX))
  m_coeff<-coeffs[1:(length(Xvars)+1)]
  ce_coeff<-coeffs[(length(Xvars)+2):length(coeffs)]
  n<-nrow(Tib)
  Xdata<-Tib[,Xvars]
  predicted<-list(Pred=list())
  for(i in 1:n){
    tmp_df<-data.frame(X=coeffs[1],freq=1,cdf=1)
    for(j in 1:length(Xvars)){
      tmpX<-Tib[[Xvars[j]]][[i]]
      tmp<-Tib[[Xvars[j]]][[i]][,1]
      tmpc<-Tib[[Xvars[j]]][[i]][,1]-mX[i,j]
      tmpf<-m_coeff[j+1]*mX[i,j]+ce_coeff[j]*tmpc
      tmpX[,1]<-tmpf
      tmp_df<-WH.sum_two_dist_dicr(tmp_df,tmpX)
    }
    predicted[[1]][[i]]<-tmp_df
  }
  return(as_tibble(predicted))
}

### GOF MEASURES -----
WH.GOFs_discr<-function(Obs,Pred){
  Obs<-as_tibble(Obs)
  Pred<-as_tibble(Pred)
  n<-nrow(Obs)
  WassD2 <- matrix(0, nrow = n, ncol = 1)
  
  MO <- mean_bargraph(Obs)
  MOm<-sum(MO[[1]]*MO[["freq"]])
  MOs<-sqrt(sum(MO[[1]]^2*MO[["freq"]])-MOm^2)
  TOTSSQ <- variance_Wass_discrete(Obs,mean_distr = MO)*n

  SSQR <- 0
  RMSE_W <- 0
  NUM_OMEGA <- 0
  DEN_OMEGA <- 0
  for (i in 1:n) {
    WassD2[i, 1] <- Wass1d_discrete_2(Obs[[1]][[i]][[1]],Pred[[1]][[i]][[1]],
                                      Obs[[1]][[i]][["freq"]],Pred[[1]][[i]][["freq"]])^2
    
    SSQR <- SSQR + Wass1d_discrete_2(Pred[[1]][[i]][[1]],MO[[1]],
                                     Pred[[1]][[i]][["freq"]],MO[["freq"]])^2
    
    OBSm <-     sum(Obs[[1]][[i]][[1]]  *Obs[[1]][[i]][["freq"]])
    OBSs <-sqrt(sum(Obs[[1]][[i]][[1]]^2*Obs[[1]][[i]][["freq"]])-OBSm^2)
    
    PREm <-     sum(Pred[[1]][[i]][[1]]  *Pred[[1]][[i]][["freq"]])
    PREs <-sqrt(sum(Pred[[1]][[i]][[1]]^2*Pred[[1]][[i]][["freq"]])-PREm^2)
    
    NUM_OMEGA <- NUM_OMEGA + (PREm - MOm)^2 + (PREs)^2
    DEN_OMEGA <- DEN_OMEGA + (OBSm - MOm)^2 + (OBSs)^2
  }
  #  browser()
  return(indices = list(
    RMSE_W = sqrt(sum(WassD2) / n),
    OMEGA = NUM_OMEGA / DEN_OMEGA,
    PSEUDOR2 = list(
      index = max(0, min(SSQR / TOTSSQ, 1 - sum(WassD2) / TOTSSQ)),
      details = c(
        TotSSQ = TOTSSQ, SSQ.R = SSQR, SSQ.E = sum(WassD2),
        Bias = TOTSSQ - SSQR - sum(WassD2), SSQ.R.rel = SSQR / TOTSSQ,
        SSQ.E.rel = sum(WassD2) / TOTSSQ,
        SSQ.bias.rel = (TOTSSQ - SSQR - sum(WassD2)) / TOTSSQ
      )
    )
  ))
  
}


### Boostrapping -----

# rr<-7
# pl<-Obs[[1]][[rr]] %>% ggplot(aes(x=it,y=cdf))+geom_line()
# pl+geom_line(data=Pred[[1]][[rr]],aes(x=X,y=cdf),color="red")

WH.regression.two.components_discr_boot <- function(Tib, #the data table
                                                    Yvar, # variable to predict
                                                    Xvars, # predictors
                                                    mY,mX,# averages of distributtions
                                                    boot_rep=1000, #boot replicates
                                                    GOFS=T){ #Do you want also GOF measures
  
  selected <- c(1:nrow(Tib))
  
  if(is_tibble(Tib)){Y <- Tib[selected, Yvar]}else{ Y <- tibble(Y=Tib[selected, Yvar])}
  
  
  names(Y)<-names(Tib)[Yvar]
  if(is_tibble(Tib)){
    X <- Tib[selected, Xvars]}else{
      X <- Tib[selected, Xvars]
      if(length(Xvars)>1){
        X <- Tib[selected, Xvars]}else{
          X <- tibble(X=Tib[selected, Xvars])
        }
    }
  
  
  n <- length(selected)
  d <- ncol(X)
  
  #create the selection matrix
  sels_M<-matrix(0,n,boot_rep)
  for(rr in 1:boot_rep){
    sels_M[,rr]<-sample(n,replace = T)
  }
  
  # extract means and do multiple regression on means
  MatAver <- cbind(mY,mX) # the matrix of means of the distributions
  colnames(MatAver) <- paste0("AV_", c(colnames(Y), colnames(X)))
  rownames(MatAver) <- rownames(Y)
  MatAver <- as.data.frame(MatAver)
  xnam <- colnames(MatAver)[2:(d + 1)]
  
  
  fmla <- as.formula(paste(colnames(MatAver)[1], paste(xnam, collapse = "+"), sep = "~"))
  
  
  # center data
  YC<-Y
  XC<-X
  
  for (i in 1:n) {
    YC[[1]][[i]][,1] <- Y[[1]][[i]][,1] - MatAver[i,1]
    for (j in 1:d) {
      XC[[j]][[i]][,1] <- X[[j]][[i]][,1] - MatAver[i,j+1]
    }
  }
  
  #compute parameters
  AV_coeffs<-matrix(0,boot_rep,d+1)
  for(rr in 1:boot_rep){
    fit <- lm(fmla, data = MatAver[sels_M[,rr],])
    AV_coeffs[rr,] <- fit$coefficient
  }
  
  colnames(AV_coeffs) <- c("(AV_Intercept)",paste0("AV_", c( colnames(X))))
  gammas <- WH.NNLS_discr_boot(XC, YC,sels_M)
  colnames(gammas) <- paste0("CEN_", colnames(X))
  parameters = cbind(AV_coeffs, gammas)
  GOFS_M<-matrix(0,boot_rep,3)
  colnames(GOFS_M)<-c("RMSE_W","OMEGA","Pseudo_R2")
  if(GOFS){
    pb <- progress_bar$new(
      format = "  Computing GOFS [:bar] :percent eta: :eta",
      total = boot_rep, clear = FALSE, width= 60)
    print("computing GOFS")
    for(rr in 1:boot_rep){
      pb$tick()
      Pred<-WH.two_comp_predict(Tib[sels_M[,rr],], 
                                Xvars,
                                mX[sels_M[,rr],],
                                parameters[rr,])
      GOF_vals<-WH.GOFs_discr(Tib[sels_M[,rr],Yvar],Pred)
      GOFS_M[rr,]=c(GOF_vals$RMSE_W,GOF_vals$OMEGA,GOF_vals$PSEUDOR2$index)
    }
    parameters=cbind(parameters,GOFS_M)
  }else{
    return(parameters = parameters)}
}

WH.NNLS_discr_boot<- function(X, Y,selected) {#selected is a matrix on n times rep 
  # Non Negative Least Squares for histogram variables modifying
  # Lawson, Charles L.; Hanson, Richard J. (1995). Solving Least Squares Problems. SIAM.
  # X<-X[selected,]
  # Y<-Y[selected,]
  rep<-ncol(selected)
  
  tol <- 1e-8
  
  XpX<- WH.mat.prod_discr_for_boot(X, X, transpose1 = TRUE)
  XpY<- WH.mat.prod_discr_for_boot(X, Y, transpose1 = TRUE)
  gammas<-matrix(0,rep,ncol(X))
  for(rept in 1:rep){
    # Lawson and Hanson NNLS algorithm
    # step 1
    P <- integer(0)
    Z <- c(1:ncol(X))
    gamma <- matrix(0, nrow = ncol(X), ncol = 1)#x
    stp <- 0
    
    XpXtmp<-apply(XpX[,,selected[,rept]], c(1,2), sum)
    #browser()
    XpYtmp<-apply(XpY[,,selected[,rep],drop=F], c(1,2), sum) #NON FUNZIONA!!
    
    w <- XpYtmp - XpXtmp %*% gamma
    while((length(Z)>0)&(max(w)>tol)){
      
      j<-which.max(w)
      P<-sort(unique(c(P,j)))
      # print(P)
      # print(w[,1])
      Z<-Z[Z!=j]#R
      s<-matrix(0,ncol(X),1)
      GG<-solve(XpXtmp[P,P]) %*%XpYtmp[P,]#S^p
      s[P,1]<-GG
      while(min(GG)<=0){
        tmp<-gamma[P,1]/(gamma[P,1]-s[P,1])
        alpha<-min(tmp[GG<=0,1])
        gamma<-gamma+alpha*(s-gamma)
        
        Z<-which(gamma<=0)
        P<-c(1:ncol(X))[-c(1:ncol(X))%in%Z]
        # print(P)
        s<-matrix(0,ncol(X),1)
        GG<-solve(XpXtmp[P,P]) %*%XpYtmp[P,]#S^p
        s[P,1]<-GG
      }
      gamma<-s
      w <- XpYtmp - XpXtmp %*% gamma
    }
    gammas[rept,]<-gamma
  }
  return(gammas)
  
}

