# ListOfDist <- vector("list", 6)
# ListOfDist[[1]] <- distributionD(c(1, 2, 3), c(0.2, 0.4, 1))
#  ListOfDist[[2]] <- distributionD(c(7, 8, 10, 15), c(0.1, 0.2, 0.7, 1))
#  ListOfDist[[3]] <- distributionD(c(9, 11, 20), c(0.3, 0.5, 1))
#  ListOfDist[[4]] <- distributionD(c(2, 5, 8), c(0.1, 0.3, 1))
#  ListOfDist[[5]] <- distributionD(c(8, 10, 15), c(0.4, 0.75, 1))
#  ListOfDist[[6]] <- distributionD(c(20, 22, 24), c(0.05, 0.12, 1))
#
#  ## create a MatH object filling it by columns
#  MyMAT <- new("MatD",
#    nrows = 3, ncols = 2, ListOfDist = ListOfDist,
#    names.rows = c("I1", "I2", "I3"), names.cols = c("Var1", "Var2"), by.row = FALSE
#  )

From_DF_to_MathD<-function(DF,nV="V1"){#two column one with names second with values
  require(dplyr)
  DF<-DF  %>% group_by(across(1:2)) %>%
    tally() %>%
    mutate(prob=n/sum(n),cdf=cumsum(prob))
  names_id<-as.vector(unique(DF[,1])[[1]])
  ListOfDist <- vector("list", length(names_id))
  cc<-0
  for(i in names_id){
    cc<-cc+1
    tmp<-DF %>% filter(.data[[names(DF)[1]]]==i)
    ListOfDist[[cc]]<-distributionD(as.vector(tmp[,2][[1]]), as.vector(tmp[,5][[1]]))
  }
   MyMAT <- new("MatD",
      nrows = length(names_id), ncols = 1, ListOfDist = ListOfDist,
      names.rows = names_id, names.cols = c("V1"), by.row = FALSE
    )
   return(MyMAT)
}
