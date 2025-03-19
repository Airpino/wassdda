targets_points<-function(highliths,routing,gra,slope){
  ID<-numeric()
  #crea_segmenti
  ini_x<-ini_y<-end_x<-end_y<-numeric()
  for(i in 1:(length(routing)-1)){
    ini_x<-c(ini_x,gra$ch$X[routing[i]])
    ini_y<-c(ini_y,gra$ch$Y[routing[i]])
    end_x<-c(end_x,gra$ch$X[routing[i+1]])
    end_y<-c(end_y,gra$ch$Y[routing[i+1]])
    
  }
  DF_seg<-data.frame(ini_x,end_x,ini_y,end_y) %>% 
    mutate(diff_x=end_x-ini_x,diff_y=end_y-ini_y,
           #diff_x=if_else(abs(diff_x)<1e-8,1e-12,diff_x),
           slo=diff_y/diff_x,con=ini_x*slo+ini_y,
           cx=diff_y,cy=-diff_x,tn=diff_y*ini_x-ini_y*diff_x)
  #trova il punto di intersezione RICOMINCIA DA QUI -------
  #pp<-gra$p
  ini_x<-end_x<-ini_y<-end_y<-IDori<-seg_CH<-tval<-numeric()
  for(i in highliths){
    ID<-c(ID,i)
    #calcola intersezione con CH
    # equazione retta passante per punto
    a<-slope#x coeff
    b<--1#y coeff
    co<-gra$transf$X[i]*slope-gra$transf$Y[i]#const coeff
    
    for (j in 1:nrow(DF_seg)){
      xy<-solve(matrix(c(DF_seg$cx[j],a,DF_seg$cy[j],b),2,2),c(DF_seg$tn[j],co))
      #find t
      
      if(DF_seg$diff_x[j]==0){
        tt<-(xy[2]-DF_seg$ini_y[j])/DF_seg$diff_y[j]
      }else{
        tt<-(xy[1]-DF_seg$ini_x[j])/DF_seg$diff_x[j]
      }
    
      if( tt>=0 && tt<=1 && xy[1]>=gra$transf$X[i] && xy[2]>=gra$transf$Y[i]) {
        #print(paste0("That's the point ===> ",i, " t ==>>  ",t))
#        pp<-pp+geom_segment(aes(x=gra$transf$X[i],y=gra$transf$Y[i],xend=xy[1],yend=xy[2]))
        ini_x<- c(ini_x,gra$transf$X[i])
        end_x<- c(end_x,xy[1])
        ini_y<- c(ini_y,gra$transf$Y[i])
        end_y<- c(end_y,xy[2])
        IDori<- c(IDori,i)
        seg_CH<- c(seg_CH,j)
        tval<-c(tval,tt)
        #geom_segment(aes(x=DF_seg$ini_x[j],y=DF_seg$ini_y[j],xend=DF_seg$end_x[j],yend=DF_seg$end_y[j]))+
        # geom_abline(slope=slope,intercept=0)+
      }
    }
    #ne escono due prendi quello con coordinate maggiori
    #ricorda di conservare i pesi
    
    #calcola lunghezza segmento
  }
  DF_ori_targ<-data.frame(IDori,seg_CH,ini_x,end_x,ini_y,end_y,tval)
  #pp<-pp+coord_fixed()
  # show(pp)
  # browser()
  # ranking efficiency
  return(list(CH_seg=DF_seg,DF_segmenti=DF_ori_targ))
}

plot_CH_map<-function(firstdim,secdim,
                      resu=resQ_expected$PCAout,highliths,alpha=1){
  DF_tmp<-data.frame(X=resu$ind$coord[,firstdim],
                     Y=resu$ind$coord[,secdim]) %>% 
    mutate(colo="black",rowID=c(1:nrow(.)),
           colo=if_else(rowID%in%highliths,"red",colo),
           alp=if_else(rowID%in%highliths,1,alpha)) %>% 
    select(X,Y,colo,alp)
  #browser()
  firstplane2<-DF_tmp %>% mutate(X2=X-min(X),Y2=Y-min(Y))
  firstplane2<-rbind(firstplane2,c(min(DF_tmp$X),min(DF_tmp$Y),NA,1,
                                   0,0))
  firstplane2<-rbind(firstplane2,c(min(DF_tmp$X),max(DF_tmp$Y),NA,1,
                                   0,
                                   max(DF_tmp$Y)-min(DF_tmp$Y)))
  firstplane2<-rbind(firstplane2,c(max(DF_tmp$X),min(DF_tmp$Y),NA,1,
                                   max(DF_tmp$X)-min(DF_tmp$X),0))
  firstplane2<-firstplane2 %>% mutate(colo=if_else(is.na(colo),"black",colo),
                                      size=if_else(colo=="red",1.5,1))
  hull <- firstplane2 %>%slice(chull(X, Y)) #contiene la sequenza dei punti
  
  p <- ggplot(firstplane2, aes(X, Y)) + 
    geom_point(data = hull,aes(X=X,Y=Y)) +
    geom_polygon(data = hull, alpha = 0.3,fill="orange") +
    geom_point(shape = 21,aes(fill=colo,size=size,alpha=alp,
                              text=as.character(c(1:nrow(firstplane2)))))+
    
    scale_size_area(
      max_size = 3)
  #browser()
  p2<-plotly::ggplotly(p)
  return(list(p=p,p2=p2,ch=hull,transf=firstplane2))
} 

plot_CH_map2<-function(firstdim,secdim,
                      resu=resQ_expected$PCAout,highliths,alpha=1){
  DF_tmp<-data.frame(X=resu$ind$coord[,firstdim],
                     Y=resu$ind$coord[,secdim]) %>% 
    mutate(colo="black",rowID=c(1:nrow(.)),
           colo=if_else(rowID%in%highliths,"red",colo),
           alp=if_else(rowID%in%highliths,1,alpha)) %>% 
    select(X,Y,colo,alp)
  #browser()
  firstplane2<-DF_tmp %>% mutate(X2=X-min(X),Y2=Y-min(Y))
  # firstplane2<-rbind(firstplane2,c(min(DF_tmp$X),min(DF_tmp$Y),NA,1,
  #                                  0,0))
  # firstplane2<-rbind(firstplane2,c(min(DF_tmp$X),max(DF_tmp$Y),NA,1,
  #                                  0,
  #                                  max(DF_tmp$Y)-min(DF_tmp$Y)))
  # firstplane2<-rbind(firstplane2,c(max(DF_tmp$X),min(DF_tmp$Y),NA,1,
  #                                  max(DF_tmp$X)-min(DF_tmp$X),0))
  firstplane2<-firstplane2 %>% mutate(colo=if_else(is.na(colo),"black",colo),
                                      size=if_else(colo=="red",1.5,1))
  hull <- firstplane2 %>%slice(chull(X, Y)) #contiene la sequenza dei punti
  
  p <- ggplot(firstplane2, aes(X, Y)) + 
    geom_point(data = hull,aes(X=X,Y=Y)) +
    geom_polygon(data = hull, alpha = 0.3,fill="orange") +
    geom_point(shape = 21,aes(fill=colo,size=size,alpha=alp,
                              text=as.character(c(1:nrow(firstplane2)))))+
    
    scale_size_area(
      max_size = 3)
  #browser()
  p2<-plotly::ggplotly(p)
  return(list(p=p,p2=p2,ch=hull,transf=firstplane2))
} 