# FUNCTION FOR DISCRETE DISTRIBUTIONS
Rcpp::sourceCpp("C:/Users/anton/Il mio Drive/CROCETTA CUSTOMER/R_code/FUNS.cpp")
library(tidyverse)
library(transport)
library(matrixStats)
library(e1071)
### extract means and standard deviations from a Tibble of Discrete distributions ----
means_stds_NestTibb<-function(MAT,labels=1){
  res<-matrix(0,nrow (MAT),ncol(MAT)-length(labels))
  res2<-matrix(0,nrow (MAT),ncol(MAT)-length(labels))
  for(i in 1:nrow(MAT)){
    c <- 0
   
    for(j in 1:ncol(MAT)){
      if (j==labels){}else{
        c <- c+1 
        tmp_m=sum(MAT[[j]][[i]][[1]]*MAT[[j]][[i]][["freq"]])
        res[i,c]=tmp_m
        res2[i,c]=sqrt(sum((MAT[[j]][[i]][[1]])^2*MAT[[j]][[i]][["freq"]])-tmp_m^2)
      }
    }
  }
  colnames(res)=names(MAT)[-labels]
  colnames(res2)=colnames(res)
  if(!is.na(labels)){
    rownames(res)=MAT[[labels]]
    rownames(res2)=MAT[[labels]]
  }
  return(list(means=res,stds=res2))
}

means_stds_NestTibb_nol<-function(Tib){
  means<-sapply(Tib,function(x){ sapply(x,function(x) sum(x[[1]]*x[["freq"]]))})
  stds<-sapply(Tib,function(x){
    sapply(x,
           function(x){
             sum(x[[1]]^2*x[["freq"]])
           }
    )
  })
  stds<-sqrt(stds-means^2)
  return(list(means=means,stds=stds))    
}

m_sd_skew_NestTibb_nol<-function(Tib){
  means<-sapply(Tib,function(x){ sapply(x,function(x) sum(x[[1]]*x[["freq"]]))})
  stds<-sapply(Tib,function(x){
    sapply(x,
           function(x){
             sum(x[[1]]^2*x[["freq"]])
           }
    )
  })
  stds<-sqrt(stds-means^2)
  skew<-sapply(Tib,function(x){sapply(x,function(x){
    m<-sum(x[[1]]*x[["freq"]])
    s<-sqrt(sum(x[[1]]^2*x[["freq"]])-m^2)
    st<-(x[[1]]-m)/s
    return(sum(st^3*x[["freq"]]))
  })})
  return(list(means=means,stds=stds,skew=skew))    
}

### counts of each tibble-cell-distribution ----
COUNTS_NestTibb<-function(MAT,labels=1){
  res<-matrix(0,nrow (MAT),ncol(MAT)-length(labels))
  for(i in 1:nrow(MAT)){
    c <- 0
    for(j in 1:ncol(MAT)){
      if (j==labels){}else{
        c <- c+1
        res[i,c]=sum(MAT[i,j][[1]]$n)
      }
    }
  }
  colnames(res)=names(MAT)[-labels]
  if(!is.na(labels)){
    rownames(res)=MAT[,labels]
  }
  return(res)
}

# computes the correlation of QQ plot of two discrete distributions ----
corr_QQ_two_discrete_distr<-function(subM1,subM2,
                                     m1,m2,
                                     std1,std2){# Computes the rhoQQ from OT theory
  require(transport)  
  # fake_cdf_1 <- round(subM1[,4],r)*(10^r)
  # fake_freq_1<- c() 
  # fake_cdf_2 <- round(subM2[,4],r)*(10^r)
  # tmp<-transport::wasserstein1d(unlist(subM1[,1]), unlist(subM2[,1]), p = 2, wa = subM1$freq, wb = subM2$freq)
  tmp<-Wass1d_discrete_2(unlist(subM1[,1]), unlist(subM2[,1]),  wa = subM1$freq, wb = subM2$freq)
  val<-tmp^2-(m1-m2)^2-(std1-std2)^2
  val<-1-val/(2*std1*std2)
  return(val)
}

# computes the correlation of QQ plot of two discrete distributions by simulation -----
corr_QQ_two_discrete_distr2<-function(subM1,subM2,
                                      m1,m2,
                                      std1,std2,
                                      r=3){# Computes the rhoQQ via approximation
  
  fake_cdf_1 <- ceiling(unlist(subM1$cdf)*(10^r))
  fake_freq_1<- c(fake_cdf_1[1],diff(fake_cdf_1)) 
  fake_cdf_2 <- ceiling(unlist(subM2$cdf)*(10^r))
  fake_freq_2<- c(fake_cdf_2[1],diff(fake_cdf_2))
  v1<-numeric()
  for(i in 1:nrow(subM1)){
    v1=c(v1,rep(subM1[i,1],fake_freq_1[i]))
  }
  v2<-numeric()
  for(i in 1:nrow(subM2)){
    v2=c(v2,rep(subM2[i,1],fake_freq_2[i]))
  }
  val<-cor(v1,v2)
  return(val)
}

WD_rqq_discr<-function(d1,d2,m1,m2,s1,s2){#script per calcolare rhoQQ
  #get common cdf
  comm_cdf<-c(0,sort(unique(c(d1$cdf,d2$cdf))))
  d1_dom<-extract_quantiles(d1,p=comm_cdf)
  d2_dom<-extract_quantiles(d2,p=comm_cdf)
  w=c(0,diff(comm_cdf))
  
  XY<-sum(d1_dom*d2_dom*w)
#  m1<-sum(d1$B1*d1$freq)
#  m2<-sum(d2$B1*d2$freq)
#  s1<-sqrt(sum(d1$B1^2*d1$freq)-m1^2)
#  s2<-sqrt(sum(d2$B1^2*d2$freq)-m2^2)
  rqq<-(XY-m1*m2)/(s1*s2)
  return(rqq)
  }


# extracts quantiles from a discrete distribution  -------
extract_quantiles<-function(subM,p=0){
  qua<-e1071::qdiscrete(p, subM$freq, values = subM[[1]])
  return(qua)
}

# sum two distributions in the sense of Wasserstein
WH.sum_two_dist_dicr<-function(d1,d2){
  comm_cdf<-c(0,sort(unique(c(d1$cdf,d2$cdf))))
  d1_dom<-extract_quantiles(d1,p=comm_cdf)
  d2_dom<-extract_quantiles(d2,p=comm_cdf)
  w=c(0,diff(comm_cdf))
  dom<-d1_dom+d2_dom
  Ds<-data.frame(X=dom[-1],freq=w[-1],cdf=comm_cdf[-1])
  
  return(Ds)
}

# computes the mean Frechet distribution using 2Wasserstein distance  ------
mean_bargraph <- function(Tibble_VAR,
                          w=1,
                          rounddec=5){
  #browser()
  data <- tibble(Tibble_VAR)
  n<-nrow(data)
  vector_of_cdf <- numeric()
  for(i in 1:n){
    #    vector_of_cdf<-c(vector_of_cdf,ceiling(Tibble_VAR[[i]]$cdf*10^rounddec)/(10^rounddec))
    vector_of_cdf<-c(vector_of_cdf,data[[1]][[i]]$cdf)
  }
  vector_of_cdf<-sort(unique(vector_of_cdf))
  if(length(w)==1){
    w<-rep(1/n,n)
  }else{
    w<-w/sum(w)
  }
  
  mean<-extract_quantiles(data[[1]][[1]],vector_of_cdf)*w[1]
  if(n>1){
    for(i in 2:n){
      mean<-mean+extract_quantiles(data[[1]][[i]],vector_of_cdf)*w[i]
    }
  }
  DF<-data.frame(x=mean,freq=diff(c(0,vector_of_cdf)),cdf=vector_of_cdf)
  return(DF)
}

# a function for checking if it works properly (Can be deleted)
mean_bargraph2 <- function(Tibble_VAR,w=1,rounddec=5){
  
  n<-length(Tibble_VAR)
  vector_of_cdf <- numeric()
  for(i in 1:n){
    vector_of_cdf<-c(vector_of_cdf,ceiling(Tibble_VAR[[i]]$cdf*10^rounddec)/(10^rounddec))
  }
  vector_of_cdf<-sort(unique(vector_of_cdf))
  if(length(w)==1){
    w<-rep(1/n,n)
  }else{
    w<-w/sum(w)
  }
  
  mean<-extract_quantiles(Tibble_VAR[[1]],vector_of_cdf)*w[1]
  DFQ<-matrix(0,length(mean),n)
  DFQ[,1]<-extract_quantiles(Tibble_VAR[[1]],vector_of_cdf)
  if(n>1){
    for(i in 2:n){
      DFQ[,i]<-extract_quantiles(Tibble_VAR[[i]],vector_of_cdf)
      mean<-mean+extract_quantiles(Tibble_VAR[[i]],vector_of_cdf)*w[i]
    }
  }
  DF<-data.frame(x=mean,
                 freq=diff(c(0,vector_of_cdf)),
                 cdf=vector_of_cdf)
  return(list(DFQ=DFQ,MEA=DF,FF=diff(c(0,vector_of_cdf))))
}

# computes the 2Wasserstein Frechet variance of a discrete distributional variable
variance_Wass_discrete<-function(Tibble_VAR,
                                 w=1,
                                 rounddec=5,
                                 given_mean=F, mean_distr=NA){
  Tibble_VAR<-tibble(Tibble_VAR)
  n<-nrow(Tibble_VAR)
  if(length(w)==1){
    w<-rep(1/n,n)
  }else{
    w<-w/sum(w)
  }
  if(!given_mean){
    mean_distr<-mean_bargraph (Tibble_VAR,w=w,rounddec=rounddec)
  }
  
  
  VV<-0
  for(i in 1:n){
    #VV<-VV+((transport::wasserstein1d(unlist(Tibble_VAR[[i]][,1]), mean_distr[,1], p = 2, 
    #                               wa = Tibble_VAR[[i]]$freq, wb = mean_distr$freq))^2)*w[i]
    tmp<-((Wass1d_discrete_2(unlist(Tibble_VAR[[1]][[i]][,1]), mean_distr[,1],  
                               wa = Tibble_VAR[[1]][[i]][["freq"]], wb = mean_distr$freq))^2)
#    print(tmp)
    VV<-VV+tmp*w[i]
  }
 # browser()
  return(VV)
}
# computes the 2Wassertein Frechet covariance between two
# discrete distributional variable2
covariance_Wass_discrete<-function( Tibble_VAR1,
                                    Tibble_VAR2,
                                    w=1,rounddec=5,
                                    given_mean1=F, mean_distr1=NA,
                                    given_mean2=F, mean_distr2=NA){
  Tibble_VAR1<-tibble(Tibble_VAR1)
  Tibble_VAR2<-tibble(Tibble_VAR2)
  n<-nrow(Tibble_VAR1)
  
  if(length(w)==1){
    w<-rep(1/n,n)
  }else{
    w<-w/sum(w)
  }
  if(!given_mean1){
    mean_distr1<-mean_bargraph (Tibble_VAR1,w=w,rounddec=rounddec)
  }
  if(!given_mean2){
    mean_distr2<-mean_bargraph (Tibble_VAR2,w=w,rounddec=rounddec)
  }
  md1 <- sum(mean_distr1[,1]*mean_distr1$freq)
  md2 <- sum(mean_distr2[,1]*mean_distr2$freq)
  sd1 <- sqrt(sum((mean_distr1[,1])^2*mean_distr1$freq)-md1^2)
  sd2 <- sqrt(sum((mean_distr2[,1])^2*mean_distr2$freq)-md2^2)
  MD1_MD2<-  corr_QQ_two_discrete_distr(mean_distr1,
                                        mean_distr2,
                                        md1,md2,sd1,sd2)*sd1*sd2+md1*md2
  CV <- 0
  
  for(i in 1:n){
    mxi<-sum(Tibble_VAR1[[1]][[i]][[1]]*Tibble_VAR1[[1]][[i]][["freq"]])
    myi<-sum(Tibble_VAR2[[1]][[i]][[1]]*Tibble_VAR2[[1]][[i]][["freq"]])
    sxi<-sqrt(sum((Tibble_VAR1[[1]][[i]][[1]])^2*Tibble_VAR1[[1]][[i]][["freq"]])-mxi^2)
    syi<-sqrt(sum((Tibble_VAR2[[1]][[i]][[1]])^2*Tibble_VAR2[[1]][[i]][["freq"]])-myi^2)
    # browser()
    xi_yi <- corr_QQ_two_discrete_distr(Tibble_VAR1[[1]][[i]],
                                        Tibble_VAR2[[1]][[i]],
                                        mxi,myi,sxi,syi)*sxi*syi+mxi*myi
    #print(xi_yi)
    CV<-CV+xi_yi*w[i]
  }
  # print(w)
  # print(CV)
  # print(MD1_MD2)
  # browser()
  CV <- CV - MD1_MD2
  return(CV)
}

cov_mat_Wass_discr<-function(Tib){
  CVmat<-matrix(0,ncol(Tib),ncol(Tib))
  for (i in 1:(ncol(Tib)-1)){
    CVmat[i,i]<-variance_Wass_discrete(Tib[,i])
    for(j in (i+1):ncol(Tib)){
      CVmat[i,j]<-CVmat[j,i]<-covariance_Wass_discrete(Tib[[i]],Tib[[j]])
    }
    
  }
  CVmat[ncol(Tib),ncol(Tib)]<-variance_Wass_discrete(Tib[,ncol(Tib)])
  return(CVmat)
}

## Correlation matrix ----
corr_mat_Wass_discr<-function(Tib){
  cov_mat<-cov_mat_Wass_discr(Tib)
  
  nr<-nrow(cov_mat)
  corr_mat<-diag(rep(1,nr))
  for(i in 1:(nr-1)){
    for(j in (i+1):nr){
      corr_mat[i,j]<-corr_mat[j,i]<-cov_mat[i,j]/sqrt(cov_mat[i,i]*cov_mat[j,j])
    }
  }
  return(corr_mat)
}

# A second version for variance
variance_Wass_discrete2<-function(Tibble_VAR,
                                  w=1,rounddec=5,
                                  given_mean=F, 
                                  mean_distr=NA){#NON CONSIDERO I PESI!! Da aggiustare
  n<-length(Tibble_VAR)
  if(length(w)==1){
    w<-rep(1/n,n)
  }else{
    w<-w/sum(w)
  }
  if(!given_mean){
    mean_distr<-mean_bargraph(Tibble_VAR,w=w,rounddec=rounddec)
  }
  md <- sum(mean_distr[,1]*mean_distr$freq)
  sd <- sqrt(sum((mean_distr[,1])^2*mean_distr$freq)-md^2)
  MDQ<-  sd^2+md^2
  
  VV<-0
  for(i in 1:n){
    mxi<-sum(Tibble_VAR[[i]][,1]*Tibble_VAR[[i]]$freq)
    sxi<-sqrt(sum((Tibble_VAR[[i]][,1])^2*Tibble_VAR[[i]]$freq)-mxi^2)
    # browser()
    xiq <- sxi^2+mxi^2
    
    
    VV<-VV+xiq*w[i]
  }
  VV<-VV-MDQ
  return(VV)
}



## DISTANCE MATRIX BETWEEN ROWS  -------

MAT_DIST_2W_discrete<-function(Tib,labels=c(1:nrow(Tib))){
  n<-nrow(Tib)
  dmat<-matrix(0,n,n)
  rownames(dmat)<-as.character(labels)
  colnames(dmat)<-rownames(dmat)
  for (i1 in 1:(n-1)){
    for (i2 in (i1+1):n){
      for(j in 1:ncol(Tib)){
        
        dmat[i1,i2]<-dmat[i1,i2]+
          #  transport::wasserstein1d
          Wass1d_discrete_2(#unname(
            Tib[i1,j][[1]][,1][[1]],#),
            #unname(
            Tib[i2,j][[1]][,1][[1]],#), 
            #p = 2, 
            wa = Tib[i1,j][[1]]$freq, 
            wb = Tib[i2,j][[1]]$freq)^2
        
      }
      dmat[i2,i1]<-dmat[i1,i2]
    }
    
  }
  
  
  dmat<-sqrt(dmat)
  return(dmat)
}
MAT_DIST_2W_discrete_and_C<-function(Tib,labels=c(1:nrow(Tib))){
  n<-nrow(Tib)
  dmat<-dmatC(Tib)
  
  
  rownames(dmat)<-as.character(labels)
  colnames(dmat)<-rownames(dmat)
  return(dmat)
  
}
### The fastest function for matrix of distances of discrete distributions ----
MAT_DIST_2W_discrete_and_C2<-function(Tib,labels=c(1:nrow(Tib))){
  n<-nrow(Tib)
  dmat<-dmatC2(Tib)
  
  
  rownames(dmat)<-as.character(labels)
  colnames(dmat)<-rownames(dmat)
  return(dmat)
  
}
## the oldes and slower one
MAT_DIST_2W_discrete_ori<-function(Tib,labels=c(1:nrow(Tib))){
  n<-nrow(Tib)
  dmat<-matrix(0,n,n)
  rownames(dmat)<-as.character(labels)
  colnames(dmat)<-rownames(dmat)
  for (i1 in 1:(n-1)){
    for (i2 in (i1+1):n){
      for(j in 1:ncol(Tib)){
        
        dmat[i1,i2]<-dmat[i1,i2]+
          transport::wasserstein1d(#unname(
            Tib[i1,j][[1]][,1][[1]],#),
            #unname(
            Tib[i2,j][[1]][,1][[1]],#), 
            p = 2, 
            wa = Tib[i1,j][[1]]$freq, 
            wb = Tib[i2,j][[1]]$freq)^2
        
      }
      dmat[i2,i1]<-dmat[i1,i2]
    }
    
  }
  
  
  dmat<-sqrt(dmat)
  return(dmat)
}



# PLOTTING FUNCTIONS ----

PLOT_tib_discrete<-function(Tib,labels=1){
  # create DF
  c<-0
  for(j in 1:ncol(Tib)){
    if (j%in%labels){}else{
      if(c==0){
        c<-c+1
        DF<-unnest(Tib[,c(labels,j)],cols = 2) %>% mutate(VAR=names(Tib)[j])
        names(DF)=c("LAB","domain","n","freq","cdf","VAR")
      }else{
        TDF<-unnest(Tib[,c(labels,j)],cols = 2) %>% mutate(VAR=names(Tib)[j])
        names(TDF)=c("LAB","domain","n","freq","cdf","VAR")
        DF<-rbind(DF,TDF)
      }
    }
    
  }
  return(DF)
}


# 1D wassertein ----


conta<-function(cub,cua){
  n=length(cub)
  cua<-c(-Inf,cua,Inf)
  m=length(cua)
  counts=rep(1,m-1)
  for (i in 1:(m-1)){
    for(j in 1:n){
      if((cub[j]<=cua[i+1])&(cub[j]>cua[i])) counts[i]=counts[i]+1
    }
  }
  return(counts)
}
Wass1d_discrete_2<-function (a, b,  
                             wa = NULL, wb = NULL) 
{
  m <- length(a)
  n <- length(b)
  if (is.null(wa)) {
    wa <- rep(1, m)
  }
  else {
    wha <- wa > 0
    wa <- wa[wha]
    a <- a[wha]
    m <- length(a)
  }
  if (is.null(wb)) {
    wb <- rep(1, n)
  }
  else {
    whb <- wb > 0
    wb <- wb[whb]
    b <- b[whb]
    n <- length(b)
  }
  ua <- (wa/sum(wa))[-m]
  ub <- (wb/sum(wb))[-n]
  cua <- c(cumsum(ua))
  cub <- c(cumsum(ub))
  arep<- countC(cub,cua)
  brep<- countC(cua,cub)
  aa <- rep(a, times = arep)
  bb <- rep(b, times = brep)
  uu <- sort(c(cua, cub),method="quick")
  uu0 <- c(0, uu)
  uu1 <- c(uu, 1)
  areap <- sqrt(sum((uu1 - uu0) * (bb - aa)^2))
  return(areap)
}

# 
# microbenchmark::microbenchmark(a=transport::wasserstein1d(unname(unlist(Tib[i1,j][[1]][,1])),
#                                                           unname(unlist(Tib[i2,j][[1]][,1])), 
#                                                           p = 2, 
#                                                           wa = Tib[i1,j][[1]]$freq, 
#                                                           wb = Tib[i2,j][[1]]$freq)^2,
#                                b=Wass1d_discrete_2(unname(unlist(Tib[i1,j][[1]][,1])),
#                                                    unname(unlist(Tib[i2,j][[1]][,1])), 
#                                                    wa = Tib[i1,j][[1]]$freq, 
#                                                    wb = Tib[i2,j][[1]]$freq)^2)

# VISUALIZATION TECNIQUES -------
## Compare two distributions via transport maps -----
plot_compare_two_distr<-function(d1,d2,m1=1,m2=10,interact=F){
  q1<-d1$cdf
  q2<-d2$cdf
  pp<-sort(unique(c(q1,q2)))
  QUA1<-extract_quantiles(d1,p=pp)
  QUA2<-extract_quantiles(d2,p=pp)
  DF<-data.frame(Q1=QUA1,Q2=QUA2,lev=pp)
  p<-ggplot(DF,aes(x=Q1,y=Q2,text=round(pp,3)))+geom_point(color="darkred")+geom_path(color="red")+
    geom_segment(aes(x = m1, y = m1, xend = m2, yend = m2),color="black")
  if(interact) p<-ggplotly(p)
  p
}

Discretize_a_distr<-function(df){
  minv<-floor(min(df[,1]))
  maxv<-ceiling(max(df[,1]))
  vals<-c(minv:maxv)
  dom<-minv
  freq<-sum(df$freq[which((df[,1]>=minv)&(df[,1]<=(vals[2]+vals[1])/2))])
  for(i in 2:(length(vals)-1)){
    dom<-c(dom,vals[i])
    freq<-c(freq,
            sum(df$freq[which((df[,1]>(vals[i-1]+vals[i])/2)&
                                (df[,1]<=(vals[i]+vals[i+1])/2))])
            )
  }
  l<-length(vals)
  dom<-c(dom,vals[l])
  freq<-c(freq,sum(df$freq[which(df[,1]>(vals[l-1]+vals[l])/2)]))
  DF<-data.frame(x=dom,freq=freq,cdf=cumsum(freq))
  return(DF)
}

# Methods of analysis -----

## To do CLustering  ----
### Kmeans ----

WD2_Kmeans<-function(Tib,k, nrep=5,labels=c(1:nrow(Tib))){
  TOL<-1e-16
  n<-nrow(Tib) #the number of rows
  p<-ncol(Tib) #the number of variables
  k<-ceiling(k)
  if (nrep<1) nrep=1
  if(k<2) stop("The number of cluster (k) must be grater than 1")
  BestSolCrit<-Inf
  BAGofSOLS<-list()
  BEST_SOL<-list()
  REP<-0
  while (REP<nrep){
    REP<-REP+1
    print(paste0("Executing REPETITION -> ", REP))
    start_time <- Sys.time()
    #-core kmeans----------------------------
    res<-KM_base(Tib,k,n,p)
    if(res$EMPTY_clu==T){
      REP<-REP-1
      print("New REP because of EMPTY CLUST")
      start_time <- Sys.time()
    }else{
      BAGofSOLS[[REP]]<- res
      if(res$Crit<BestSolCrit){
        BestSolCrit <- res$Crit
        BEST_SOL<- res
      }
      end_time <- Sys.time()
      print(paste0("Elapsed time ", end_time-start_time))
    }
      #-end core kmeans----------------------------
    }
 return(BEST_SOL)  
}

KM_base<-function(Tib,k,n,p,maxiter=100){
  TOL<-1e-16
  # n<-nrow(Tib) #the number of rows
  # p<-ncol(Tib) #the number of variables
  # k<-ceiling(k)
  # 
  # if(k<2) stop("The number of cluster (k) must be grater than 1")
  # BestSolCrit=Inf
  
  # INITIALIZE CENTERS
  ActualCrit=Inf
  # sample k centers
  centers<-Tib[sample(1:n,k),]
  STOP<-F
  ID<-0
  Criterion<-Inf
  iter<-0
  Min_d_to_prot<-0
  EMPTY_CLUST<-F
 # browser()
  while(!STOP){
    iter<-iter+1
    ## AFFECT
    EMPTY_CLUST<-F
    #compute distance
    # Dist_to_prot<-matrix(0,n,k)
    # for(i in 1:n){
    #   for (clu in 1:k){
    #     for (j in 1:p){
    #       Dist_to_prot[i,clu]=Dist_to_prot[i,clu]+
    #         DW2_discr_C(a=Tib[i,j][[1]][[1]],
    #                     b=centers[clu,j][[1]][[1]],
    #                     wa=Tib[i,j][[1]]$freq,
    #                     wb=centers[clu,j][[1]]$freq)^2
    #     }
    #   }
    # }
    
    Dist_to_prot<-DiToCen_C(Tib,
              centers,
              n, p, k)
    #browser()
    #first assignment
    IDclu<-max.col(-Dist_to_prot, ties="first")
    #check empty clusters
    if(length(table(IDclu))<k){#there is an empty cluster
      EMPTY_CLUST<-TRUE
      STOP<-T
      print("Empty cluster found!")
    }
   
    if(!EMPTY_CLUST){
      ## CENTER
      #compute new centers
      for (clu in 1:k){
        for (j in 1:p){
          centers[[j]][[clu]] <- tibble(mean_bargraph(Tib[which(IDclu==clu),j]))
        }
      }
      
      ## CRITERION
      # compute criterion
      Min_d_to_prot<-rowMins(Dist_to_prot)
      Criterion<-sum(Min_d_to_prot)
      
      ## CHECKS
      
      if (((ActualCrit-Criterion)>TOL)&(iter<maxiter)){
        print(paste0("IT ---> ",iter," Crit --->", Criterion))
        ActualCrit <- Criterion
      }else{
        ID<-IDclu
        STOP<-T
      }
    }
  }
  resu<-list(ID=ID,Proto=centers,Crit=Criterion,EMPTY_clu=EMPTY_CLUST)
  return(resu)
}

  ### Hierarchical ----
WD2_HC<-function(Tib,method = "complete",labels=c(1:nrow(Tib))){
  DMAT<-MAT_DIST_2W_discrete_and_C2(Tib,labels)
  Hc<-hclust(as.dist(DMAT),method = method )
  return(Hc)
}


  
  ### Regression  ----
  ### MFA ----

### Wass dist multiple ----

Wass_SQ_discr_MULT<-function(Row1,Row2){
  ncols<-ncol(Row1)
  DD<-DM<-DS<-DSH<-RHO<-M1<-M2<-S1<-S2<-numeric()
  # browser()
  for(i in 1:ncols){
    d2<-(Wass1d_discrete_2(Row1[[i]][[1]][[1]],Row2[[i]][[1]][[1]],Row1[[i]][[1]]$freq,Row2[[i]][[1]]$freq))^2
    m1<-sum(Row1[[i]][[1]][[1]]*Row1[[i]][[1]]$freq/sum(Row1[[i]][[1]]$freq))
    m2<-sum(Row2[[i]][[1]][[1]]*Row2[[i]][[1]]$freq/sum(Row2[[i]][[1]]$freq))
    M1<-c(M1,m1)
    M2<-c(M2,m2)
    
    dm<-(m1-m2)^2
    s1<-sqrt(sum(Row1[[i]][[1]][[1]]^2*Row1[[i]][[1]]$freq/sum(Row1[[i]][[1]]$freq))-m1^2)
    s2<-sqrt(sum(Row2[[i]][[1]][[1]]^2*Row2[[i]][[1]]$freq/sum(Row2[[i]][[1]]$freq))-m2^2)
    S1<-c(S1,s1)
    S2<-c(S2,s2)
    ds<-(s1-s2)^2
    dsh<-(d2-dm-ds)
    rho<-1-dsh/(2*s1*s2)
    DD<-c(DD,d2)
    DM<-c(DM,dm)
    DS<-c(DS,ds)
    DSH<-c(DSH,dsh)
    RHO<-c(RHO,rho)
  }
  details<-data.frame(DD,DM,DS,DSH,RHO,M1,M2,S1,S2)
  return(list(DIST=sqrt(sum(DD)),DIST2=sum(DD),details=details))
}

## transport plots or QQ -----
Transp_plots_MULT<-function(Row1,Row2,q=1000){
  ncols<-ncol(Row1)
  one_q<-ref_q<-IDV<-numeric()
  for(i in 1:ncols){
    
    #browser()
    one<-Row1[[i]][[1]]
    ref<-Row2[[i]][[1]]
    one_q<-c(one_q,extract_quantiles(one,p=c(0:q)/q))
    ref_q<-c(ref_q,extract_quantiles(ref,p=c(0:q)/q))
    IDV<-c(IDV,rep(i,q+1))
    # tmpl[[i]]<-data.frame(one_q=one_q,ref_q=ref_q)
  }
  tmpl<-data.frame(IDV,one_q,ref_q)
  return(tmpl)
}


Cronbach_discr <-function(TibbM){
  cov_m <- cov_mat_Wass_discr(TibbM)
  Bstat <-m_sd_skew_NestTibb_nol(TibbM)
  
  cov_means<-cov(Bstat$means)*(nrow(TibbM)-1)/nrow(TibbM)
  
  cov_distrib <- cov_m -cov_means
  
  p <- ncol(TibbM)
  
  cro_al<-p/(p-1)*(1-  sum(diag(cov_m))/sum(cov_m))
  cro_al_m<-p/(p-1)*(1-  sum(diag(cov_means))/sum(cov_means))
    cro_al_v<-p/(p-1)*(1-  sum(diag(cov_distrib))/sum(cov_distrib))
  
  
  cro_without<-cro_without_m<-cro_without_v<-numeric()
  for (j in 1:p){
    cro_without<-c(cro_without,
                   (p-1)/(p-2)*(1-  sum(diag(cov_m[-j,-j]))/sum(cov_m[-j,-j])))
    cro_without_m<-c(cro_without_m,
                   (p-1)/(p-2)*(1-  sum(diag(cov_means[-j,-j]))/sum(cov_means[-j,-j])))
    cro_without_v<-c(cro_without_v,
                   (p-1)/(p-2)*(1-  sum(diag(cov_distrib[-j,-j]))/sum(cov_distrib[-j,-j])))
  }
  res<-list(cro_al=c(GLOB=cro_al,AVE=cro_al_m,VAR=cro_al_v), 
                     cro_without=cro_without,
                     cro_without_m=cro_without_m,
                     cro_without_v=cro_without_v,
            STATS=Bstat)
  return(res)
}

