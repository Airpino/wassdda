#Viz tools
library(tidyverse)
library(ggdendro)
library(patchwork)
library(ggeasy)

Heat_map_Disc<-function(Tib,
                        labels=c(1:nrow(Tib)),#the row names
                        sortR=T,
                        sortC=T
){
  Tib<-as_tibble(Tib)
  n<-nrow(Tib)
  col<-ncol(Tib)
  #browser()
  if(sortR) {
    hcR<-WD2_HC(Tib)
    orderR<-hcR$order}else{orderR<-c(1:n)}
  if(sortC) {
    #qua devo fare il clustering sulle variabili
    covmat<-cov_mat_Wass_discr(Tib)
    cormat<-covmat/t(t(sqrt(diag(covmat))))%*%t(sqrt(diag(covmat)))
    hcC<-hclust(as.dist(1-cormat))
    orderC<-hcC$order
  }else{orderC<-c(1:ncol(Tib))}
  
  
  
  #createDF
  IDr<-IDc<-ORr<-ORc<-valueD<-freq<-numeric()
  for(i in 1:n){
    for(j in 1:col){
      tmpx<-Tib[[j]][[i]][[1]]
      tmpf<-Tib[[j]][[i]][["cdf"]]
      tmpf<-Tib[[j]][[i]][["freq"]]
      tmpl<-length(tmpx)
      IDr<-c(IDr,rep(i,tmpl))
      IDc<-c(IDc,rep(j,tmpl))
      ORr<-c(ORr,rep(orderR[i],tmpl))
      ORc<-c(ORc,rep(orderC[j],tmpl))
      valueD<-c(valueD,tmpx)
      freq<-c(freq,tmpf)
    }
  }
  DF<-data.frame(IDr=IDr,IDc=IDc,
                 ORr=ORr,ORc=ORc,
                 dom=valueD, freq=freq)
  Ranges<-DF %>% group_by(IDc,ORc) %>% summarize(a=min(dom),b=max(dom))
  DF<-DF %>% mutate(IDc2=IDc+0.1*dom-0.05)## Questo solo per voti da 1 10
  p<-ggplot(DF,aes(x=ORr,y=freq))+geom_bar(stat="identity",position=position_fill(reverse = TRUE),
                                           aes(fill=as.factor(dom)))+
    scale_fill_brewer(palette = "RdYlGn",direction=1)+coord_flip()+
    facet_grid(~ ORc)+theme_void()+theme(panel.spacing.x=unit(-0.2, "lines"),
                                         plot.margin = margin(0,0,0,0,"cm"))#+
  
  
  # Rectangular lines
  ddata1 <- dendro_data(hcR, type = "rectangle")
  p1 <- ggplot(segment(ddata1)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    coord_flip() + 
    scale_y_reverse(expand = c(0, 0))+theme_dendro()+theme(plot.margin = margin(0,0,0,0,"cm"))
  ddata2 <- dendro_data(hcC, type = "rectangle")
  p2 <- ggplot(segment(ddata2)) +scale_y_continuous(expand = c(0, 0))+ 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()+
    theme(plot.margin = margin(0,0,0,0,"cm"))
  
  pfin<-plot_spacer()+p2+p1+p+
    plot_layout(heights = c(1, 5), widths = c(1, 4))
  
  
  #devo creare le tiles?
  return(list(DF=DF,pfin=pfin))  
}

squish_data<-function(Tib,labels=c(1:nrow(Tib))) {
  Tib<-as_tibble(Tib)
  n<-nrow(Tib)
  col<-ncol(Tib)
  colna<-colnames(Tib)
  IDr<-IDc<-valueD<-cdf<-freq<-numeric()
  labID<-labVar<-character()
  for(i in 1:n){
    for(j in 1:col){
      tmpx<-Tib[[j]][[i]][[1]]
      tmpcdf<-Tib[[j]][[i]][["cdf"]]
      tmpf<-Tib[[j]][[i]][["freq"]]
      tmpl<-length(tmpx)
      IDr<-c(IDr,rep(i,tmpl))
      IDc<-c(IDc,rep(j,tmpl))
      labID<-c(labID,rep(as.character(labels[i]),tmpl))
      labVar<-c(labVar,rep(colna[j],tmpl))
      valueD<-c(valueD,tmpx)
      freq<-c(freq,tmpf)
      cdf<-c(cdf,tmpcdf)
    }
  }
  DF<-data.frame(IDr=IDr,IDc=IDc,
                 labID=labID,labVar=labVar,
                 dom=valueD, freq=freq,cdf=cdf)
  return(DF)  
}  

Performance_plot<-function(Tib,selected=1,
                           labels=c(1:nrow(Tib)),TITLE=T,notick=F,
                           skewness_plo=FALSE,
                           alpha=1,
                           bg="transparent",BW=F,
                           polar=T,legend=F){
  DF<-squish_data(Tib[selected,],labels = labels[selected])
  if(skewness_plo){
  bstat<-DF %>% 
    group_by(labVar) %>% 
    summarize(
      m=sum(dom*freq),
      s=sqrt(sum(dom^2*freq)-m^2),
      sk=(sum(((dom-m)/s)^3*freq))
    ) %>% ungroup() %>% mutate(sk=sign(sk)*abs(sk)^(1/3), coord=(sk+3)/6) 
  }
  DF$labVar<-factor(DF$labVar,levels=colnames(Tib))
  #stacked perc bars
  # cc<-c( "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B","#D9EF8B", 
  #                 "#A6D96A" ,"#66BD63", "#1A9850", "#006837")
  if (BW){
    
    cols <- c("0" = "#fefefe", "1" = "#f1f1f1", "2" = "#e1e1e1", "3" = "#d9d9d9",
              "4" = "#dddddd", "5" = "#cccccc", "6" = "#a1a1a1", "7" = "#818181",
              "8" = "#666666", "9" = "#414141", "10" = "#000000")  
  }else{
  cols <- c("0" = "#A50026", "1" = "#A50026", "2" = "#D73027", "3" = "#F46D43",
            "4" = "#FDAE61", "5" = "#FEE08B", "6" = "#D9EF8B", "7" = "#66BD63",
            "8" = "#1A9850", "9" = "#006837", "10" = "#003300")
  }
  p<-ggplot(DF,aes(x=labVar,y=freq))
  if (legend){p<-p+geom_bar(stat="identity",
                            aes(fill=as.factor(dom)),
                            width=1,alpha=alpha)}else
  {  p<-p+geom_bar(stat="identity",
             aes(fill=as.factor(dom)),
             width=1,
             show.legend = F,alpha=alpha)
  }
  
  p<-p+
    scale_fill_manual(
      values = cols)+
    #    scale_fill_brewer(palette = "Spectral",direction=1)+
    #    scale_fill_brewer(palette = "RdYlGn",direction=1)+
    geom_hline(yintercept = 0.5,
               linetype="dashed",
               #color="white",
               #size=1,
               alpha=0.6)
  if(skewness_plo){
    p<-p+geom_point(inherit.aes = F,data=bstat,
                    aes(x=labVar,y=coord),
                    color="black",fill="pink",shape=21,
                    alpha=0.6,show.legend = F)
  }
    
  if(polar) p<-p+coord_polar()+theme_minimal()
  p<-p+theme_minimal()
  if(TITLE) p<-p+
    labs(#title = labels[selected]#,
      #subtitle = "subtitle",
      caption = labels[selected]
    )+
    theme(plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5))
if(notick){
  p<-p+
    theme(axis.text.x=element_blank())
}
  p<-p+
    theme(axis.line=element_blank(),
          #        axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.border=element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(color="transparent",fill = "transparent"),
          panel.background = element_rect(color="transparent",fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = bg, color = NA),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))
  return(p) 
}

Performance_plot_maps<-function(Tib,selected=1,
                                labels=c(1:nrow(Tib)),BW=F){
  DF<-squish_data(Tib[selected,],labels = labels[selected])
  DF$min_long<-(DF$IDc-(min(DF$IDc)))/(diff(range(DF$IDc))+1)*360
  DF$max_long<-(DF$IDc+1-(min(DF$IDc)))/(diff(range(DF$IDc))+1)*360
  DF$diff=DF$cdf-DF$freq
  minr<-65
  maxr<-90
  dd<-0
  DF$min_lat<-DF$diff*(maxr-minr)+minr#+(90-minr-dd)
  DF$max_lat<-DF$cdf*(maxr-minr)+minr
  inrl<-(90-minr-dd)/2+minr
  if (BW){
    
    cols <- c("0" = "#fefefe", "1" = "#f1f1f1", "2" = "#e1e1e1", "3" = "#d9d9d9",
              "4" = "#dddddd", "5" = "#cccccc", "6" = "#a1a1a1", "7" = "#818181",
              "8" = "#666666", "9" = "#414141", "10" = "#000000")  
  }else{
    cols <- c("0" = "#A50026", "1" = "#A50026", "2" = "#D73027", "3" = "#F46D43",
              "4" = "#FDAE61", "5" = "#FEE08B", "6" = "#D9EF8B", "7" = "#66BD63",
              "8" = "#1A9850", "9" = "#006837", "10" = "#003300")
  }
  texts<-unique(DF$labVar)
  DF2<-data.frame(y=rep(min(DF$min_lat-5),length(texts)),
                  x=(c(1:length(texts))-0.5)/length(texts)*360,
                  label=texts)
  
  p<-ggplot(DF)+geom_rect(aes(xmin=min_long,xmax=max_long,
                              ymin=min_lat,ymax=max_lat,fill=as.factor(dom)),
                          show.legend = F)+
    geom_hline(yintercept = inrl,linetype="dashed")+
    geom_text(inherit.aes = F,data=DF2,aes(y=y,
                                           x=x,label=label),size=3)+
    scale_fill_manual(
      values = cols)+
    #scale_fill_brewer(palette = "RdYlGn",direction=1)+
    
    coord_map("ortho",orientation = c(90, 180, 0),
              ylim = c(min(DF$min_lat)-20,90))+theme_void()+
    ggtitle(labels[selected])
  return(p) 
}




# Pmap<-Performance_plot_maps(Fin_tibble[,2:13],selected = 1,labels=labs)
# 
# pp<-Performance_plot_maps(Fin_tibble[,2:13],selected = 1,labels=labs)+
#   Performance_plot_maps(Fin_tibble[,2:13],selected = 126,labels=labs)+
#   Performance_plot_maps(Fin_tibble[,2:13],selected = 142,labels=labs)+
#   Performance_plot_maps(Fin_tibble[,2:13],selected = 4,labels=labs)
# 
# 
# pp2<-  Performance_plot(Fin_tibble[,2:13],selected = 126,labels=labs)+
#   Performance_plot(Fin_tibble[,2:13],selected = 142,labels=labs)+
#   Performance_plot(Fin_tibble[,2:13],selected = 1,labels=labs)+
# Performance_plot(Fin_tibble[,2:13],selected = 4,labels=labs)

#create rectangles



## Come pseudo-ordinare delle distribuzione nel senso di Wasserstein
## Massa media-> media o mediana
## Tutti i quantili di 1 sono >= di 2
## Una buona distribuzione di valutazioni deve avere
## Mediana/media alta
## Varaibilità bassa
## skeness a sinistra

Performance_plot_MUL<-function(Tib,selected=1,
                           labels=c(1:nrow(Tib)),
                           #TITLE=T,
                           notick=T,
                           col_wrap=F,
                           col_labs=NULL,
                           #skewness_plo=FALSE,
                           alpha=1,
                           bg="transparent"){
  
  DF<-squish_data(Tib[selected[1],],labels = labels[selected[1]])
  DF$PLO<-labels[1]
  if(col_wrap)DF$CW<-col_labs[1]
  
  if(length(selected)>1){
    for (j in 2:length(selected)){
  DFtmp<-squish_data(Tib[selected[j],],labels = labels[selected[j]])
  DFtmp$PLO<-labels[j]
  if(col_wrap) DFtmp$CW<-col_labs[j]
  DF<-rbind(DF,DFtmp)
    }
  }
  DF$PLO<-factor(DF$PLO,levels=labels)
  # if(skewness_plo){
  #   bstat<-DF %>% 
  #     group_by(labVar) %>% 
  #     summarize(
  #       m=sum(dom*freq),
  #       s=sqrt(sum(dom^2*freq)-m^2),
  #       sk=(sum(((dom-m)/s)^3*freq))
  #     ) %>% ungroup() %>% mutate(sk=sign(sk)*abs(sk)^(1/3), coord=(sk+3)/6) 
  # }
  DF$labVar<-factor(DF$labVar,levels=colnames(Tib))
  #stacked perc bars
  # cc<-c( "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE08B","#D9EF8B", 
  #                 "#A6D96A" ,"#66BD63", "#1A9850", "#006837")
  cols <- c("0" = "#A50026", "1" = "#A50026", "2" = "#D73027", "3" = "#F46D43",
            "4" = "#FDAE61", "5" = "#FEE08B", "6" = "#D9EF8B", "7" = "#66BD63",
            "8" = "#1A9850", "9" = "#006837", "10" = "#003300")
  p<-ggplot(DF,aes(x=labVar,y=freq))+
    geom_bar(stat="identity",
             aes(fill=as.factor(dom)),
             width=1,
             show.legend = F,alpha=alpha)+
    scale_fill_manual(
      values = cols)+
    geom_hline(yintercept = 0.5,
               linetype="dashed",
               alpha=0.6)
  
  p<-p+coord_polar()+theme_minimal()
  if(col_wrap){
    p<-p+facet_grid(rows = vars(PLO),cols = vars(CW))
  }else{
    p<-p+facet_grid(rows = vars(PLO))
    }
  # if(TITLE) p<-p+
  #   labs(#title = labels[selected]#,
  #     #subtitle = "subtitle",
  #     caption = labels[selected]
  #   )+
  #   theme(plot.title = element_text(hjust = 0.5),
  #         plot.caption = element_text(hjust = 0.5))
  if(notick){
    p<-p+
      theme(axis.text.x=element_blank())
  }
  p<-p+
    theme(axis.line=element_blank(),
          #        axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.border=element_blank(),
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(color="transparent",fill = "transparent"),
          panel.background = element_rect(color="transparent",fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = bg, color = NA),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0))
  return(p) 
}
