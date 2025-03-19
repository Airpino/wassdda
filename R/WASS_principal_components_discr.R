library(mgcv)
WH.1d.PCA_discr <- function(data, var=1, quantiles = 10, plots = TRUE, 
                            listaxes = c(1:4), axisequal = FALSE){#
                      #, qcut = 1, outl = 0) 
  
  # Build the matrix of Quantiles
  VARS <- ncol(data)
  INDIV <- nrow(data)
  if (missing(var)) {
    var <- 1
    varname <- colnames(data)[1]
    cat(paste("Var is missing, We do a PCA on variable ", varname, "\n"))
  } else {
    if (length(var) > 1) {
      varname <- colnames(data)[var[1]]
      cat(paste("Var has several values, We do a PCA on variable -->", var[1], "\n"))
      var <- var[1]
    }
    else {
      if (var > VARS) {
        stop(paste("The variables are less than ", var))
      }
      else {
        varname <- colnames(data)[var]
        cat(paste("We do a PCA on variable ---> ", varname, "\n"))
      }
    }
  }
  data <- tibble(data[, var])
  # check indiv and remove empty rows
  # toremove <- numeric(0)
  # for (i in 1:INDIV) {
  #   if (length(data@M[i, 1][[1]]@x) < 2) {
  #     toremove <- c(toremove, i)
  #   }
  # }
  # if (length(toremove) > 0) {
  #   data@M <- as.matrix(data@M[-toremove, 1])
  # }
  # colnames(data@M) <- varname
   INDIV <- nrow(data)
  if (INDIV < 2) {
    stop("At least two individuals are necessary, please check data")
  }

  # Create the quantile matrix
  MATQ <- matrix(0, nrow = INDIV, ncol = (quantiles + 1))
  p <- c(0, 1:quantiles) / quantiles
  
    namec <- c("Min")
    for (j in 2:(quantiles)) {
      namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
    }
    namec <- c(namec, "Max")
  
  for (i in 1:INDIV)
  {
         MATQ[i, ] <- extract_quantiles(data[[1]][[i]], p)
      }
  rownames(MATQ) <- rownames(data)
  colnames(MATQ) <- namec
# check for quantile vars with no variation.
  # do a PCA
  
  vmeans <- numeric(0)
  vstds <- numeric(0)

  for (ind in 1:INDIV) {
    tmp_m<-sum(data[[1]][[ind]][[1]]*data[[1]][[ind]][["freq"]])
    vmeans <- c(vmeans,tmp_m)
    tmp_s<-sqrt(sum(data[[1]][[ind]][[1]]^2*data[[1]][[ind]][["freq"]])-tmp_m^2)
    vstds <- c(vstds, tmp_s)
  }
  minM <- min(vmeans)
  maxM <- max(vmeans)
  minS <- min(vstds)
  maxS <- max(vstds)

  MATF <- cbind(MATQ / sqrt((quantiles + 1)), vmeans, vstds)
  #browser()
  res.pca <- FactoMineR::PCA(MATF, quanti.sup = c((ncol(MATF) - 1), ncol(MATF)),
                             scale.unit = FALSE, graph = plots)

  VARW <- 0 # WH.var.covar(data)
  TOTINE <- sum(res.pca$eig[, 1])
  ## plotting PCA results ----------------
  plots<-F
  if (plots) {
    # par(mfrow=c(1,1))
    # lt's define the couple of axes
    dev.new(noRStudioGD = FALSE)
    planes <- matrix(0, 1, 2)
    for (cc in 1:floor(length(listaxes) / 2)) {
      planes <- rbind(planes, c(listaxes[(cc * 2 - 1)], listaxes[(cc * 2)]))
    }
    planes <- as.matrix(planes[-1, ])
    dim(planes) <- c(floor(length(listaxes) / 2), 2)
    if ((length(listaxes) %% 2) == 1) {
      planes <- rbind(planes, c(
        listaxes[(length(listaxes) - 1)],
        listaxes[length(listaxes)]
      ))
    }

    # plot Spanish-fan plot
    for (pl in 1:nrow(planes)) {
      axe1 <- planes[pl, 1]
      axe2 <- planes[pl, 2]

      labX <- paste("Axis ", axe1, " (", format(res.pca$eig[axe1, 2], digits = 2, nsmall = 2), "%)")
      labY <- paste("Axis ", axe2, " (", format(res.pca$eig[axe2, 2], digits = 2, nsmall = 2), "%)")
      CVAR <- res.pca$var$coord[, c(axe1, axe2)]
      plot.new()

      if (axisequal) {
        xrange <- c(min(0, min(cbind(CVAR[, 1], CVAR[, 2]))), max(0, max(cbind(CVAR[, 1], CVAR[, 2]))))
        yrange <- xrange
        plot.window(xrange, yrange)
        title(xlab = labX)
        title(ylab = labY)

        #         plot(type="n",
        #              xlab=labX,ylab=labY)
        segments(xrange[1], 0, xrange[2], 0)
        segments(yrange[1], 0, yrange[2], 0)
      }
      else {
        plot.window(
          c(min(0, CVAR[, 1]), max(0, CVAR[, 1])),
          c(min(0, CVAR[, 2]), max(0, CVAR[, 2]))
        )
        title(xlab = labX)
        title(ylab = labY)
        #       plot(c(min(0,CVAR[,1]),max(0,CVAR[,1])), c(min(0,CVAR[,2]),max(0,CVAR[,2])),type="n",
        #            xlab=labX,ylab=labY)Y)
        segments(min(0, CVAR[, 1]), 0, max(0, CVAR[, 1]), 0)
        segments(0, min(0, CVAR[, 2]), 0, max(0, CVAR[, 2]))
      }
      title("PCA Variable plot (Spanish-fan plot)")
      if ((nrow(CVAR) %% 2) == 0) {
        centr <- nrow(CVAR) / 2
      } else {
        centr <- (nrow(CVAR) - 1) / 2
      }
      centr <- nrow(CVAR) / 2
      cxl <- 0.8
      for (tr in 2:nrow(CVAR)) {
        centrality <- (abs(tr - centr - 1) / (centr - 1))

        # cat(centrality,"\n")
        x <- c(0, CVAR[(tr - 1), 1], CVAR[tr, 1], 0)
        y <- c(0, CVAR[(tr - 1), 2], CVAR[tr, 2], 0)
        red <- 1
        green <- centrality
        # cat(red,green,"\n")
        polygon(x, y, col = rgb(red, green, 0, 0.7), lty = 1, lwd = 1, border = "black")
      }
      text(x = CVAR[1, 1], y = CVAR[1, 2], labels = colnames(MATQ)[1], cex = cxl)
      if (nrow(CVAR) < 8) {
        for (tr in 2:(nrow(CVAR) - 1)) {
          text(x = CVAR[tr, 1], y = CVAR[tr, 2], labels = colnames(MATQ)[tr], cex = cxl)
        }
      }
      else {
        text(
          x = CVAR[ceiling(nrow(CVAR) / 8), 1], y = CVAR[ceiling(nrow(CVAR) / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) / 8)], cex = cxl
        )
        text(
          x = CVAR[ceiling(nrow(CVAR) * 2 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 2 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 2 / 8)], cex = cxl
        )
        text(
          x = CVAR[ceiling(nrow(CVAR) * 3 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 3 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 3 / 8)], cex = cxl
        )
        text(x = CVAR[centr + 1, 1], y = CVAR[centr + 1, 2], labels = colnames(MATQ)[centr + 1], cex = cxl)
        text(
          x = CVAR[ceiling(nrow(CVAR) * 5 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 5 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 5 / 8)], cex = cxl
        )
        text(
          x = CVAR[ceiling(nrow(CVAR) * 6 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 6 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 6 / 8)], cex = cxl
        )
        text(
          x = CVAR[ceiling(nrow(CVAR) * 7 / 8), 1], y = CVAR[ceiling(nrow(CVAR) * 7 / 8), 2],
          labels = colnames(MATQ)[ceiling(nrow(CVAR) * 7 / 8)], cex = cxl
        )
      }
      text(
        x = CVAR[nrow(CVAR), 1], y = CVAR[nrow(CVAR), 2],
        labels = colnames(MATQ)[ncol(MATQ)], cex = cxl
      )
      axis(1, pos = 0, cex.axis = 0.7) # , at=xM,labels=labX)
      axis(2, pos = 0, cex.axis = 0.7) # , at=yM,labels=labY)
      dev.new(noRStudioGD = FALSE)

      # plot individuals

      # layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

      xm <- min(res.pca$ind$coord[, axe1])
      xM <- max(res.pca$ind$coord[, axe1])
      ym <- min(res.pca$ind$coord[, axe2])
      yM <- max(res.pca$ind$coord[, axe2])
      Xlen <- xM - xm
      Ylen <- yM - ym

      # Compute scaling factor for domain
      matX <- matrix(0, 128, INDIV)
      matD <- matrix(0, 128, INDIV)
      for (ind in 1:INDIV) {
        distr <- data@M[ind, 1][[1]]
        # generate 200 random points according to the QF
        rn <- 200
        xn <- c(rep(0, rn))
        random_no <- c(0:rn) / rn

        for (i in 1:rn) {
          xn[i] <- compQ(distr, random_no[i])
        }
        d <- density(xn, n = 128)
        matX[, ind] <- d$x
        matD[, ind] <- d$y
      }
      MinX <- min(matX)
      MaxX <- max(matX)
      RX <- MaxX - MinX
      q95D <- quantile(matD, probs = qcut)
      matD[matD > q95D] <- as.numeric(q95D)
      MaxD <- max(matD)

      meanXRange <- diff(range(apply(matX, 2, FUN = range)))
      xfact <- 0.3
      yfact <- 0.1
      if (axisequal) {
        xm1 <- min(xm, ym)
        xM1 <- max(xM, yM)
        ym1 <- xm1
        yM1 <- xM1
        Xlen <- xM1 - xm1
        Ylen <- yM1 - ym1
        xm <- (xM1 + xm1) / 2 - Xlen / 2 * (1 + 1.5 * xfact)
        xM <- (xM1 + xm1) / 2 + Xlen / 2 * (1 + 1.5 * xfact)
        ym <- (yM1 + ym1) / 2 - Ylen / 2 * (1 + 0.5 * yfact)
        yM <- (yM1 + ym1) / 2 + Ylen / 2 * (1 + 1.5 * yfact)
      } else {
        xm1 <- xm
        xM1 <- xM
        ym1 <- ym
        yM1 <- yM
        xm <- (xM1 + xm1) / 2 - Xlen / 2 * (1 + 1.5 * xfact)
        xM <- (xM1 + xm1) / 2 + Xlen / 2 * (1 + 1.5 * xfact)
        ym <- (yM1 + ym1) / 2 - Ylen / 2 * (1 + 0.5 * yfact)
        yM <- (yM1 + ym1) / 2 + Ylen / 2 * (1 + 1.5 * yfact)
      }
      dev.new(noRStudioGD = FALSE)
      plot.new()
      plot.window(c(xm, xM), c(ym, yM))
      title(main = "PCA plot of distributions", sub = "Smoothed distributions")
      for (ind in 1:INDIV) {
        x <- matX[, ind]
        y <- matD[, ind]
        x <- (x - data@M[ind, 1][[1]]@m) / meanXRange * xfact * Xlen + res.pca$ind$coord[ind, axe1]
        x <- c(x, x[length(x)], x[1])
        y <- c(y, 0, 0) / q95D * yfact * Ylen + res.pca$ind$coord[ind, axe2]
        polygon(x, y, col = rgb(1, 1, 0, 0.6))
      }
      for (ind in 1:INDIV) {
        text(res.pca$ind$coord[ind, axe1], res.pca$ind$coord[ind, axe2],
          label = rownames(MATQ)[ind], pos = 1, cex = 0.7
        )
      }
      axis(1, pos = 0, cex.axis = 0.7) # , at=xM,labels=labX)
      axis(2, pos = 0, cex.axis = 0.7) # , at=yM,labels=labY)
      segments(xm, 0, xM, 0)
      segments(0, ym, 0, yM)

      title(xlab = labX)
      title(ylab = labY)
      dev.new(noRStudioGD = FALSE)
      # plot.new()
      plot.new()
      plot.window(c(xm, xM), c(ym, yM))
      title(main = "Coloured using the mean values")
      trasp <- 0.6

      for (ind in 1:INDIV) {
        x <- matX[, ind]
        y <- matD[, ind]
        x <- (x - data@M[ind, 1][[1]]@m) / meanXRange * xfact * Xlen + res.pca$ind$coord[ind, axe1]
        x <- c(x, x[length(x)], x[1])
        y <- c(y, 0, 0) / q95D * yfact * Ylen + res.pca$ind$coord[ind, axe2]
        red <- (data@M[ind, 1][[1]]@m - minM) / (maxM - minM)
        if (red > 0.5) {
          green <- 1 - red
          blue <- 0
        }
        else {
          green <- 1 - red
          red <- 0
          blue <- 1 - green
        }
        polygon(x, y, col = rgb(red, green, blue, trasp))
      }
      for (ind in 1:INDIV) {
        text(res.pca$ind$coord[ind, axe1], res.pca$ind$coord[ind, axe2],
          label = rownames(MATQ)[ind], pos = 1, cex = 0.7
        )
      }
      axis(1, pos = 0, cex.axis = 0.7) # , at=xM,labels=labX)
      axis(2, pos = 0, cex.axis = 0.7) # , at=yM,labels=labY)
      segments(xm, 0, xM, 0)
      segments(0, ym, 0, yM)

      abline(0, res.pca$quanti.sup$coord[1, 2] / res.pca$quanti.sup$coord[1, 1], lty = 3, lwd = 2, col = "red")
      title(xlab = labX)
      title(ylab = labY)
      dev.new(noRStudioGD = FALSE)
      plot.new()
      plot.window(c(xm, xM), c(ym, yM))
      title(main = "Coloured using std values")
      for (ind in 1:INDIV) {
        x <- matX[, ind]
        y <- matD[, ind]
        x <- (x - data@M[ind, 1][[1]]@s) / meanXRange * xfact * Xlen + res.pca$ind$coord[ind, axe1]
        x <- c(x, x[length(x)], x[1])
        y <- c(y, 0, 0) / q95D * yfact * Ylen + res.pca$ind$coord[ind, axe2]
        red <- (data@M[ind, 1][[1]]@s - minS) / (maxS - minS)
        if (red > 0.5) {
          green <- 1 - red
          blue <- 0
        }
        else {
          green <- 1 - red
          red <- 0
          blue <- 1 - green
        }
        polygon(x, y, col = rgb(red, green, blue, trasp))
      }
      for (ind in 1:INDIV) {
        text(res.pca$ind$coord[ind, axe1], res.pca$ind$coord[ind, axe2],
          label = rownames(MATQ)[ind], pos = 1, cex = 0.7
        )
      }
      axis(1, pos = 0, cex.axis = 0.7) # , at=xM,labels=labX)
      axis(2, pos = 0, cex.axis = 0.7) # , at=yM,labels=labY)
      segments(xm, 0, xM, 0)
      segments(0, ym, 0, yM)

      abline(0, res.pca$quanti.sup$coord[2, 2] / res.pca$quanti.sup$coord[2, 1], lty = 3, lwd = 2, col = "red")

      title(xlab = labX)
      title(ylab = labY)
    }
  }
  res.pca$quantile_matrix<-MATF
  return(list(PCAout = res.pca, WASSVARIANCE = VARW, INERTIA = TOTINE))
}

## MFA PCA for quantiles -------
WH.MultiplePCA_discr <- function(data, list.of.vars, 
                                 quantiles = 10,
                                 labels=c(1:nrow(data)),
                                 ncp_in=max(3,ceiling(length(list.of.vars)))/3,
                                 pp=c(0,1)){#, outl = 0) 
  
  data <- data[, list.of.vars]
  VARS <- ncol(data)
  if(VARS<2) {
    print("Selected variables must be at least 2, EXECUTION HALTED!!")
    return(NULL)}
  
  INDIV <- nrow(data)
  MATQ <- matrix(0, nrow = INDIV, ncol = VARS * (quantiles + 1))
  # if ((outl < 0) || (outl > 0.5)) {
  #   outl <- 0
  # }
  # if ((outl == 0.5)) {
  #   outl <- 0.49
  # }
  # browser()
  # create the names of the variables
  #################################################
  COLnames <- list()
  outl<-0
  if (outl == 0) {
    p <-seq(pp[1],pp[2],length.out=quantiles+1)
#    p <- c(0, 1:quantiles) / quantiles
    for (v in 1:VARS) {
      namec <- list("Min")
      for (j in 2:(quantiles)) {
        namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
      }
      namec <- c(namec, "Max")
      COLnames <- c(COLnames, paste(substr(colnames(data)[v], 1, 4), namec, sep = "."))
      for (i in 1:INDIV)
      {
        # tmp=compQ_vect(data@M[i,v][[1]],p)
        MATQ[i, ((v - 1) * (quantiles + 1)+c(1:(quantiles + 1)))] <-extract_quantiles(data[[v]][[i]], p)
        # for (j in 1:(quantiles + 1)) {
        #   MATQ[i, ((v - 1) * (quantiles + 1) + j)] <- compQ(data@M[i, v][[1]], p[j])
        # }
      }
    }
  }
  else {
    p <- c(0, 1:quantiles) / quantiles * (1 - 2 * outl) + outl
    namec <- list(paste0("Q(", format(p[1], digits = 2, nsmall = 2), ")"))
    for (j in 2:length(p)) {
      namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
    }
    for (v in 1:VARS) {
      namec <- list(paste0("Q(", format(p[1], digits = 2, nsmall = 2), ")"))
      for (j in 2:length(p)) {
        namec <- c(namec, paste0("Q(", format(p[j], digits = 2, nsmall = 2), ")"))
      }
      COLnames <- c(COLnames, paste(substr(colnames(data)[v], 1, 4), namec, sep = "."))
      for (i in 1:INDIV)
      {
        MATQ[i, ((v - 1) * (quantiles + 1)+c(1:(quantiles + 1)))] <-extract_quantiles(data[[v]][[i]], p)
        
        # for (j in 1:(quantiles + 1)) {
        #   MATQ[i, ((v - 1) * (quantiles + 1) + j)] <- compQ(data@M[i, v][[1]], p[j])
        # }
      }
    }
  }
 

  rownames(MATQ) <- labels
  colnames(MATQ) <- COLnames

  # do a MFA
  MATF <- cbind(MATQ / sqrt((quantiles + 1)))
  # check all columns
  #supp<-c()
  for (cc in 1:ncol(MATF)) {
    #supp<-c(supp,sd(MATF[, cc]))
    if (sd(MATF[, cc]) < 1e-8) {
      
      MATF[, cc] <- MATF[, cc] + rnorm(nrow(MATF), 0, 1e-8)
    }
  }
  #print(supp)
  #browser()
  MFA.res <- FactoMineR::MFA(MATF,
    group = c(rep(quantiles + 1, length(list.of.vars))), type = c(rep("c", length(list.of.vars))),
    ncp = ncp_in, 
    name.group = colnames(data),
    graph = TRUE
  )
  MFA.res$quantile_matrix<-MATQ
  return(MFA.res)
}

WH.plot_multiple_Spanish.funs_2<-function (res, axes = c(1, 2), var = 1, LABS = TRUE, multi = TRUE, 
          corplot = TRUE) 
{
  if (multi) {
    labX <- paste("Axis ", axes[1], " (", format(res$eig[axes[1], 
                                                         2], digits = 2, nsmall = 2), "%)")
    labY <- paste("Axis ", axes[2], " (", format(res$eig[axes[2], 
                                                         2], digits = 2, nsmall = 2), "%)")
  }
  else {
    labX <- paste("Axis ", axes[1], " (", format(res$PCAout$eig[axes[1], 
                                                                2], digits = 2, nsmall = 2), "%)")
    labY <- paste("Axis ", axes[2], " (", format(res$PCAout$eig[axes[2], 
                                                                2], digits = 2, nsmall = 2), "%)")
  }
  v1 <- numeric()
  v2 <- numeric()
  tt <- character()
  tt2 <- numeric()
  VA <- numeric()
  x1 <- numeric()
  x2 <- numeric()
  hj <- numeric()
  vj <- numeric()
  if (!multi) {
    var <- 1
  }
  for (j in var) {
    if (multi) {
      els <- which(res$summary.quanti$group == j)
    }
    else {
      els <- c(1:nrow(res$PCAout$var$coord))
    }
    tmp_v1 <- numeric()
    tmp_v2 <- numeric()
    tmp_tt <- character()
    tmp_tt2 <- numeric()
    tmp_VA <- numeric()
    if (multi) {
      if (corplot) {
        # CVAR <- res$global.pca$var$cor[els, c(axes[1], 
        #                                       axes[2])]
         CVAR <- res$quanti.var$cor[els, c(axes[1], 
                                              axes[2])]
      }
      else {
        # CVAR <- res$global.pca$var$coord[els, c(axes[1], 
        #                                         axes[2])]
        CVAR <- res$quanti.var$coord[els, c(axes[1], 
                                          axes[2])]
      }
    }
    else {
      if (corplot) {
        CVAR <- res$PCAout$var$cor[els, c(axes[1], axes[2])]
      }
      else {
        CVAR <- res$PCAout$var$coord[els, c(axes[1], 
                                            axes[2])]
      }
    }
    for (pp in 1:(nrow(CVAR) - 1)) {
      tmp_v1 <- c(tmp_v1, CVAR[pp, 1], 0, CVAR[pp + 1, 
                                               1], CVAR[pp, 1])
      tmp_v2 <- c(tmp_v2, CVAR[pp, 2], 0, CVAR[pp + 1, 
                                               2], CVAR[pp, 2])
      tmp_tt <- c(tmp_tt, paste(as.character(j), rep(as.character(pp), 
                                                     4), sep = "_"))
      tmp_tt2 <- c(tmp_tt2, rep(pp, 4))
      tmp_VA <- c(tmp_VA, rep(j, 4))
    }
    tmp_tt2 <- abs(tmp_tt2 - mean(tmp_tt2))
    tmp_tt2 <- 1 - tmp_tt2/max(tmp_tt2)
    v1 <- c(v1, tmp_v1, 0, 0, 0, 0)
    v2 <- c(v2, tmp_v2, 0, 0, 0, 0)
    tt <- c(tt, tmp_tt, paste(as.character(j), rep(as.character(pp), 
                                                   4), sep = "_"))
    tt2 <- c(tt2, tmp_tt2, 1, 1, 1, 1)
    tmp_VA <- c(tmp_VA, rep(j, 4))
    VA <- c(VA, tmp_VA)
    x1 <- c(x1, CVAR[, 1])
    x2 <- c(x2, CVAR[, 2])
  }
  d <- data.frame(v1 = v1, v2 = v2, tt = tt, tt2 = tt2, VA = VA)
  txt_lab <- data.frame(v1 = x1, v2 = x2, hj = (-sign(x1) + 
                                                  1)/2, vj = (-sign(x2) + 1)/2)
  txt_lab$labs<-rownames(CVAR)
  angle <- seq(-pi, pi, length = 250)
  x_c <- sin(angle)
  y_c <- cos(angle)
  df <- data.frame(x_c = x_c, y_c = y_c)
  l1 <- data.frame(x = c(-1.2, 1.2), y = c(0, 0))
  l2 <- data.frame(x = c(0, 0), y = c(-1.2, 1.2))
  p <- ggplot(d, aes(x = v1, y = v2)) + geom_polygon(aes(group = tt, 
                                                         fill = as.factor(VA), alpha = tt2), colour = "grey30")
  if (LABS == TRUE) {
    p <- p + geom_text(data = txt_lab, aes(x = v1, y = v2, 
                                           label = labs, hjust = hj, vjust = vj))
  }
  if (corplot) {
    p <- p + xlim(-1.2, 1.2) + ylim(-1.2, 1.2) + ggtitle("Correlation plot of variables (Spanish-fan plot)")
  }
  else {
    XL <- range(d$v1)
    XL <- (XL - mean(XL)) * 1.2 + mean(XL)
    YL <- range(d$v2)
    YL <- (YL - mean(YL)) * 1.2 + mean(YL)
    p <- p + xlim(XL) + ylim(YL) + ggtitle("Plot of variables (Spanish-fan plot)")
  }
  p <- p + theme_bw() + xlab(labX) + ylab(labY) + theme(legend.position = "none")
  if (corplot) {
    p <- p + coord_fixed(ratio = 1)
  }
  p <- p + theme(plot.title = element_text(hjust = 0.5))
  if (corplot) {
    p <- p + geom_path(data = df, aes(x = x_c, y = y_c), 
                       inherit.aes = F, linetype = 2)
  }
  p <- p + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
  #print(p)
  return(p)
}


# WH.plot_multiple_Spanish.funs_discr <- function(res,
#                                           axes = c(1, 2), var = 1, LABS = TRUE,
#                                           multi = TRUE, corplot = TRUE) {
#   # require(ggplot2)
#   if (multi) {
#     labX <- paste("Axis ", axes[1], " (", format(res$eig[axes[1], 2],
#       digits = 2, nsmall = 2
#     ), "%)")
#     labY <- paste("Axis ", axes[2], " (", format(res$eig[axes[2], 2],
#       digits = 2, nsmall = 2
#     ), "%)")
#   }
#   else {
#     labX <- paste("Axis ", axes[1], " (", format(res$PCAout$eig[axes[1], 2],
#       digits = 2, nsmall = 2
#     ), "%)")
#     labY <- paste("Axis ", axes[2], " (", format(res$PCAout$eig[axes[2], 2],
#       digits = 2, nsmall = 2
#     ), "%)")
#   }
#   v1 <- numeric()
#   v2 <- numeric()
#   tt <- character()
#   tt2 <- numeric()
#   VA <- numeric()
#   x1 <- numeric()
#   x2 <- numeric()
#   hj <- numeric()
#   vj <- numeric()
#   if (!multi) {
#     var <- 1
#   }
#   for (j in var) {
# 
#     # plot Spanish-fan plot
#     if (multi) {
#       els <- which(res$summary.quanti$group == j)
#     }
#     else {
#       els <- c(1:nrow(res$PCAout$var$coord))
#     }
#     tmp_v1 <- numeric()
#     tmp_v2 <- numeric()
#     tmp_tt <- character()
#     tmp_tt2 <- numeric()
#     tmp_VA <- numeric()
#     if (multi) {
#       if (corplot) {
#         CVAR <- res$global.pca$var$cor[els, c(axes[1], axes[2])]
#       }
#       else {
#         CVAR <- res$global.pca$var$coord[els, c(axes[1], axes[2])]
#       }
#     }
#     else {
#       if (corplot) {
#         CVAR <- res$PCAout$var$cor[els, c(axes[1], axes[2])]
#       }
#       else {
#         CVAR <- res$PCAout$var$coord[els, c(axes[1], axes[2])]
#       }
#     }
# 
#     for (pp in 1:(nrow(CVAR) - 1)) {
#       tmp_v1 <- c(tmp_v1, CVAR[pp, 1], 0, CVAR[pp + 1, 1], CVAR[pp, 1])
#       tmp_v2 <- c(tmp_v2, CVAR[pp, 2], 0, CVAR[pp + 1, 2], CVAR[pp, 2])
#       tmp_tt <- c(tmp_tt, paste(as.character(j), rep(as.character(pp), 4), sep = "_"))
#       tmp_tt2 <- c(tmp_tt2, rep(pp, 4))
#       tmp_VA <- c(tmp_VA, rep(j, 4))
#     }
#     tmp_tt2 <- abs(tmp_tt2 - mean(tmp_tt2))
#     tmp_tt2 <- 1 - tmp_tt2 / max(tmp_tt2)
#     v1 <- c(v1, tmp_v1, 0, 0, 0, 0)
#     v2 <- c(v2, tmp_v2, 0, 0, 0, 0)
#     tt <- c(tt, tmp_tt, paste(as.character(j), rep(as.character(pp), 4), sep = "_"))
#     tt2 <- c(tt2, tmp_tt2, 1, 1, 1, 1)
#     tmp_VA <- c(tmp_VA, rep(j, 4))
#     VA <- c(VA, tmp_VA)
#     x1 <- c(x1, CVAR[, 1])
#     x2 <- c(x2, CVAR[, 2])
#   }
# 
#   d <- data.frame(v1 = v1, v2 = v2, tt = tt, tt2 = tt2, VA = VA)
#   # CVAR=res$global.pca$var$cor[,c(axes[1],axes[2])]
#   txt_lab <- data.frame(v1 = x1, v2 = x2, hj = (-sign(x1) + 1) / 2, vj = (-sign(x2) + 1) / 2)
#   ### the circle of correlation
#   angle <- seq(-pi, pi, length = 250)
#   x_c <- sin(angle)
#   y_c <- cos(angle)
#   df <- data.frame(x_c = x_c, y_c = y_c)
#   ##
#   l1 <- data.frame(x = c(-1.2, 1.2), y = c(0, 0))
#   l2 <- data.frame(x = c(0, 0), y = c(-1.2, 1.2))
# 
# 
#   p <- ggplot(d, aes(x = v1, y = v2)) +
#     geom_polygon(aes(group = tt, fill = as.factor(VA), alpha = tt2), colour = "grey30")
#   if (LABS == TRUE) {
#     p <- p + geom_text(data = txt_lab, aes(
#       x = v1, y = v2,
#       label = rownames(txt_lab),
#       hjust = hj, vjust = vj
#     ))
#   }
#   if (corplot) {
#     p <- p + xlim(-1.2, 1.2) + ylim(-1.2, 1.2) +
#       ggtitle("Correlation plot of variables (Spanish-fan plot)")
#   }
#   else {
#     XL <- range(d$v1)
#     XL <- (XL - mean(XL)) * 1.2 + mean(XL)
#     YL <- range(d$v2)
#     YL <- (YL - mean(YL)) * 1.2 + mean(YL)
#     p <- p + xlim(XL) + ylim(YL) +
#       ggtitle("Plot of variables (Spanish-fan plot)")
#   }
# 
#   p <- p + theme_bw() +
#     xlab(labX) + ylab(labY) +
#     theme(legend.position = "none")
#   if (corplot) {
#     p <- p + coord_fixed(ratio = 1)
#   }
#   p <- p + theme(plot.title = element_text(hjust = 0.5))
#   if (corplot) {
#     p <- p + geom_path(data = df, aes(x = x_c, y = y_c), inherit.aes = F, linetype = 2)
#   }
#   p <- p + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
# 
# 
#   print(p)
# 
#   return(p)
# }
# 
# 
# WH.plot_multiple_indivs_discr <- function(data, res, axes = c(1, 2), indiv = 0, var = 1,
#                                     strx = 0.1, stry = 0.1, HISTO = TRUE, coor = 0,
#                                     stat = "mean") {
#   # require(ggplot2)
# 
#   if (indiv[1] == 0) {
#     el_ind <- c(1:nrow(data))
#   } else {
#     el_ind <- indiv
#   }
#   if (length(coor) == 1) {
#     coord <- res$ind$coord[el_ind, axes]
#   } else {
#     coord <- coor
#   }
#   # browser()
# 
#   if (HISTO) {
#     x <- numeric()
#     y <- numeric()
#     ID_O <- numeric()
#     minx <- numeric()
#     mea <- numeric()
#     rmax <- -Inf
#     for (i in el_ind) {
#       tmp_r <- max(data@M[i, var][[1]]@x) - data@M[i, var][[1]]@x[1] + 1e-100
#       mm <- data@M[i, var][[1]]@m
#       if (tmp_r > rmax) {
#         rmax <- tmp_r
#       }
#       vv <- COMP_POLY_H(data@M[i, var][[1]])
#       x <- c(x, min(vv$x), vv$x, max(vv$x))
#       y <- c(y, 0, vv$y, 0)
#       ID_O <- c(ID_O, rep(i, length(vv$x) + 2))
#       mea <- c(mea, rep(mm, length(vv$x) + 2))
#       minx <- c(minx, rep(min(vv$x), length(vv$x) + 2))
#     }
#   } else {
#     # compute densities
#     pt <- 300
#     x <- numeric()
#     y <- numeric()
#     ID_O <- numeric()
#     minx <- numeric()
#     mea <- numeric()
#     rmax <- -Inf
# 
#     for (i in el_ind) {
#       vals <- COMP_MQ(data@M[i, var][[1]], Mp = c(0:pt) / pt)
#       mm <- data@M[i, var][[1]]@m
#       tmp_r <- max(vals) - min(vals)
#       if (tmp_r > rmax) {
#         rmax <- tmp_r
#       }
#       vv <- density(vals, n = 128, from = min(vals), to = max(vals))
#       x <- c(x, min(vv$x), vv$x, max(vv$x))
#       y <- c(y, 0, vv$y, 0)
#       ID_O <- c(ID_O, rep(i, length(vv$x) + 2))
#       mea <- c(mea, rep(mm, length(vv$x) + 2))
#       minx <- c(minx, rep(min(vals), length(vv$x) + 2))
#     }
# 
#     ### end of compute density
#   }
#   ID_o <- ID_O
#   x_o <- x
#   y_o <- y
#   ALL_OBJ <- data.frame(x_o = x_o, y_o = y_o, ID_o = ID_o, minx = minx, mea = mea)
#   limX <- range(ALL_OBJ$x_o)
# 
#   EXTRY <- quantile(ALL_OBJ$y_o, probs = 0.995)
#   # bisogna traslare ogni oggetto e scalarlo
#   # scaliamo gli oggetti tra 0 e 1
#   ALL_OBJ$x <- (ALL_OBJ$x_o - minx) / (rmax)
#   # tagliamo picchi eccessivi
#   ALL_OBJ$y[ALL_OBJ$y_o > EXTRY] <- EXTRY
#   limY <- quantile(ALL_OBJ$y_o, probs = 0.99)
#   ALL_OBJ$y <- ALL_OBJ$y_o / limY
#   # ogni oggetto va traslato e rimpicciolito di un fattore strech x strech y
#   dim1MAX <- max(coord[, 1])
#   dim1Min <- min(coord[, 1])
#   rX <- (dim1MAX - dim1Min)
#   dim2MAX <- max(coord[, 2])
#   dim2Min <- min(coord[, 2])
#   rY <- (dim2MAX - dim2Min)
#   cc <- 0
#   for (i in el_ind) {
#     cc <- cc + 1
#     xc <- coord[cc, 1]
#     yc <- coord[cc, 2]
#     tmp_x <- ALL_OBJ$x_o[which(ALL_OBJ$ID_o == i)]
#     tmp_y <- ALL_OBJ$y_o[which(ALL_OBJ$ID_o == i)]
#     tmp_x <- tmp_x * strx * rX
#     tmp_y <- tmp_y * stry * rY
#     tmp_x <- (tmp_x - mean(tmp_x)) + xc
#     tmp_y <- tmp_y + yc
#     ALL_OBJ$x_o[which(ALL_OBJ$ID_o == i)] <- tmp_x
#     ALL_OBJ$y_o[which(ALL_OBJ$ID_o == i)] <- tmp_y
#   }
# 
# 
#   title <- paste(
#     "Plot of individuals: distributions for variable",
#     colnames(data@M)[var]
#   )
#   labX <- paste("Axis ", axes[1], " (", format(res$eig[axes[1], 2],
#     digits = 2, nsmall = 2
#   ), "%)")
#   labY <- paste("Axis ", axes[2], " (", format(res$eig[axes[2], 2],
#     digits = 2, nsmall = 2
#   ), "%)")
# 
#   colnames(coord) <- c("x", "y")
#   coord <- as.data.frame(coord)
#   p <- ggplot(data = ALL_OBJ, aes(x = x_o, y = y_o)) +
#     geom_polygon(aes(group = ID_o, fill = as.factor(ID_o), alpha = 0.4), colour = "black") +
#     geom_point(data = coord, aes(x = x, y = y)) +
#     geom_text(
#       data = coord,
#       aes(
#         x = x, y = y, label = rownames(coord),
#         hjust = 0.5, vjust = 1
#       )
#     ) +
#     theme_bw() +
#     xlab(labX) +
#     ylab(labY) +
#     ggtitle(title) +
#     theme(legend.position = "none") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     geom_vline(xintercept = 0) +
#     geom_hline(yintercept = 0)
#   ### add a line for mean snd sd
#   Ms <- get.MatH.stats(data[el_ind, var], stat = stat)$mat
#   # SD=get.MatH.stats(data[el_ind,var],stat="std")$mat
#   # Skew=get.MatH.stats(data[el_ind,var],stat="std")$mat
#   # kurtH()=get.MatH.stats(data[el_ind,var],stat="std")$mat
#   ### print a plot with means
# 
#   coord <- cbind(coord, Ms)
#   colnames(coord)[[3]] <- stat
# 
#   p2 <- ggplot(ALL_OBJ, aes(x = x_o, y = y_o)) +
#     scale_fill_gradient(low = "yellow", high = "red") +
#     geom_polygon(aes(group = ID_o, fill = mea, alpha = 0.4), colour = "black") +
#     geom_point(data = coord, aes(x = x, y = y)) +
#     theme_bw() +
#     geom_text(
#       data = coord,
#       aes(
#         x = x, y = y,
#         label = format(mean, digits = 1, nsmall = 2),
#         hjust = 0.5, vjust = 1
#       )
#     ) +
#     xlab(labX) +
#     ylab(labY) +
#     ggtitle(paste(title, "each label is the ", stat, ".")) +
#     theme(legend.position = "none") +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     geom_vline(xintercept = 0) +
#     geom_hline(yintercept = 0)
# 
#   #+ geom_abline(intercept = 0, slope = cor(coord)[1,3]*sd(coord$x),linetype=2)
#   # browser()
#   print(p)
#   print(p2)
# 
#   return(pp = list(p, p2))
# }
# 
# COMP_MQ <- function(object, Mp) {
#   res <- numeric()
#   n <- length(Mp)
# 
#   for (i in 1:n) {
#     p <- Mp[i]
#     # Check for errors
#     if (p < 0 || p > 1) stop("p must be a value between 0 and 1")
# 
#     if (p <= 0) {
#       q <- object@x[1]
#     } else {
#       if (p >= 1) {
#         q <- object@x[length(object@x)]
#       }
#       else {
#         ini <- max(object@x[object@p <= p])
#         pos1 <- which.max(object@x[object@p <= p])
#         pos2 <- pos1 + 1
#         fin <- object@x[pos2]
#         if (ini == fin) {
#           q <- ini
#         }
#         else {
#           q <- ini + (object@x[pos2] - object@x[pos1]) * (p - object@p[pos1]) / (object@p[pos2] - object@p[pos1])
#         }
#       }
#       res <- c(res, q)
#     }
#   }
#   return(res)
# }
# 
# COMP_POLY_H <- function(object) {
#   # browser()
#   x <- object@x[1]
#   y <- 0
#   if (length(object@x) < 2) {
#     x <- c(x, x, x)
#     y <- c(y, 1, 0)
#   }
#   else {
#     for (i in 2:length(object@x)) {
#       tmp_dens <- (object@p[i] - object@p[i - 1]) / (object@x[i] - object@x[i - 1] + 1e-100)
#       x <- c(x, object@x[i - 1], object@x[i])
#       y <- c(y, tmp_dens, tmp_dens)
#     }
#     x <- c(x, object@x[length(object@x)])
#     y <- c(y, 0)
#   }
# 
#   return(res = data.frame(x = x, y = y))
# }

reconstr_MFA<-function(res,ncp,newco,quantiles){
  
  coord.var <- sweep (as.matrix(res$quanti.var$coord)[,1:ncp, drop=F],
                      1,res$call$col.w, FUN="*")
  hatX <- as.matrix(newco) %*% ## ATTENTION HERE !!!!----
  t(sweep(coord.var,2,sqrt(res$eig[1:ncp,1]),FUN="/")) 
  
  
  ecarttype <- res$separate.analyses[[1]]$call$ecart.type
  moy <- res$separate.analyses[[1]]$call$centre
  for(g in 2:length(res$call$group)){
    ecarttype <- c(ecarttype,res$separate.analyses[[g]]$call$ecart.type)
    moy <-       c(moy,res$separate.analyses[[g]]$call$centre)
  }
  
  hatX <-sweep(hatX,2,ecarttype,FUN = "*")
  hatX <-sweep(hatX,2,res$call$col.w,FUN = "/")
  hatX <-sweep(hatX,2,moy,FUN = "+")
  
  Recon_new<-hatX*sqrt(quantiles+1)
}

reconstr_PCA<-function(res, ncp ,newco,quantiles) 
{
  
  ncp <- min(ncp, ncol(res$ind$coord))
  coord.var <- as.matrix(res$var$coord)[, 1:ncp, drop = F]
  
  hatX <- as.matrix(newco) %*% 
    t(sweep(coord.var, 2, sqrt(res$eig[1:ncp, 
                                       1]), FUN = "/"))
  
  hatX <- sweep(hatX, 2, res$call$ecart.type, FUN = "*")
  hatX <- sweep(hatX, 2, res$call$centre, FUN = "+")
  Recon_new<-hatX*sqrt(quantiles+1)
  return(Recon_new)
}


# this function takes a recontructed sequence of quantiles and adjust it to be monotonic
monotonize_seq<-function(x,quantiles){ # a carpenteer function for 
  vect<-x
  
  pattern<-sapply(c(1:(quantiles+1)),function(x){sum(vect>=vect[x])})-c((quantiles+1):1)
  #print(pattern)
  maxit<-100
  it<-1
  while((min(pattern)<0)&&(it<maxit)){
    #plot(x)
    it<-it+1
    checks<-(which(pattern<0))
    for(i in 1:length(checks)) vect[checks[i]]=vect[checks[i]+1]
    pattern<-sapply(c(1:(quantiles+1)),function(x){sum(vect>=vect[x])})-c((quantiles+1):1)
    #points(vect,col="red")
    #print(pattern)
    #browser()
    
  }
  browser()
  return(vect)
}

monotonize_seq2<-function(x,quantiles){ # a carpenteer function for 
  vect<-x
  
  pattern<-sapply(c(1:(quantiles+1)),function(x){sum(vect>=vect[x])})-c((quantiles+1):1)
  #print(pattern)
  maxit<-100
  it<-1
  while((min(pattern)<0)&&(it<maxit)){
    #plot(x)
    it<-it+1
    checks<-(which(pattern<0))
    for(i in 1:length(checks)) {
      if(checks[i]==1){
        vect[checks[i]]=vect[checks[i]+1]
      }else{
        vect[checks[i]]=mean(c(vect[checks[i]-1],vect[checks[i]+1]))
      }
      
    }
    pattern<-sapply(c(1:(quantiles+1)),function(x){sum(vect>=vect[x])})-c((quantiles+1):1)
  }
  return(vect)
}

monotonize_seq3<-function(x,cdf_e,plot=F){ # a carpenteer function for 
  #solution from https://stats.stackexchange.com/questions/197509/how-to-smooth-data-and-force-monotonicity
  
  df <- data.frame(x=cdf_e, y=x)
  
  ## Set up the size of the basis functions/number of knots
  k <- max(c(3,length(cdf_e)-sum(diff(x)<=0)-5))
  ## This fits the unconstrained model but gets us smoothness parameters that
  ## that we will need later
  unc <- gam(y ~ s(x, k = k, bs = "cr"), data = df)
  
  ## This creates the cubic spline basis functions of `x`
  ## It returns an object containing the penalty matrix for the spline
  ## among other things; see ?smooth.construct for description of each
  ## element in the returned object
  sm <- smoothCon(s(x, k = k, bs = "cr"), df, knots = NULL)[[1]]
  
  ## This gets the constraint matrix and constraint vector that imposes
  ## linear constraints to enforce montonicity on a cubic regression spline
  ## the key thing you need to change is `up`.
  ## `up = TRUE` == increasing function
  ## `up = FALSE` == decreasing function (as per your example)
  ## `xp` is a vector of knot locations that we get back from smoothCon
  F1 <- mono.con(sm$xp, up = TRUE)   # get constraints: up = FALSE == Decreasing constraint!
  
  ## Fill in G, the object pcsl needs to fit; this is just what `pcls` says it needs:
  ## X is the model matrix (of the basis functions)
  ## C is the identifiability constraints - no constraints needed here
  ##   for the single smooth
  ## sp are the smoothness parameters from the unconstrained GAM
  ## p/xp are the knot locations again, but negated for a decreasing function
  ## y is the response data
  ## w are weights and this is fancy code for a vector of 1s of length(y)
  G <- list(X = sm$X, C = matrix(0,0,0), sp = unc$sp,
            p = sm$xp, # note the  here! This is for decreasing fits!
            y = df$y,
            w = df$y*0+1)
  G$Ain <- F1$A    # the monotonicity constraint matrix
  G$bin <- F1$b    # the monotonicity constraint vector, both from mono.con
  G$S <- sm$S     # the penalty matrix for the cubic spline
  G$off <- 0      # location of offsets in the penalty matrix
  
  ## Do the constrained fit 
  p <- pcls(G)  # fit spline (using s.p. from unconstrained fit)
  
  ## predict at 100 locations over range of x - get a smooth line on the plot
  newx <- with(df, data.frame(x = seq(min(x), max(x), length = quantiles+1)))
  
  fv <- Predict.matrix(sm, newx) %*% p
  
  newx <- transform(newx, yhat = fv[,1])
  if(plot){browser()}
  fin<-data.frame(x=newx$yhat,cdf=newx$x,freq=c(0,diff(newx$x)),oldx=x)
  fin<-fin %>% mutate(x=round(x,5)) %>% 
    group_by(x) %>% 
    summarize(cdf=max(cdf)) %>% 
    mutate(freq=c(cdf[1],diff(cdf))) %>% ungroup()
  fin<-as.data.frame(fin[,c(1,3,2)])
  return(fin)
}

# this function approximate a continuous distr to a discrete with a fixed domain
discretize_continuous_distr<-function(x,p,dom){
  quantiles<-length(p)-1
  domx<-x
  #browser()
  if(min(diff(x))<0){ 
    x<-monotonize_seq4(x,p)
    domx<-c(x$x[1]-1e-5,unlist(x$x))
    p<-c(0,unlist(x$cdf))
  }
  if(p[1]>0){
    domx<-c(x$x[1]-1e-5,dom(x))
    p<-c(0,p)
  }
  dom_agg<-(dom[1:(length(dom)-1)]+dom[2:length(dom)])*0.5
  # dom_agg<-sort(unique(c(dom,dom_agg)))
    di<-distributionH(domx,p)
  dom_agg<-sort(unique(c(dom_agg,Inf)))
  pro<-sapply(dom_agg,function(x)compP(di,x))
  pro[length(pro)]=1
#  if(max(di@x)>max(dom_agg)) pro[length(pro)]=1
  df<-as.data.frame(cbind(dom,freq=c(pro[1],diff(pro)),cdf=pro))
  
  df<-df[which(df$freq>0),]
  return(df)
}


monotonize_seq4<-function(x,cdf_e,R2m=0.97,plot=F){ # a carpenteer function for 
  #solution from https://stats.stackexchange.com/questions/197509/how-to-smooth-data-and-force-monotonicity
  maxit=100
  df <- data.frame(x=cdf_e, y=x)
  R2=0
  red=5
  it=0
  while((R2<R2m)&&(it<maxit)){
    it=it+1
  ## Set up the size of the basis functions/number of knots
  k <- max(c(5,length(cdf_e)-red))
 
  unc <- gam(y ~ s(x, k = k, bs = "cr"), data = df)
  
  sm <- smoothCon(s(x, k = k, bs = "cr"), df, knots = NULL)[[1]]
  
  F1 <- mono.con(sm$xp, up = TRUE)   # get constraints: up = FALSE == Decreasing constraint!
  
  G <- list(X = sm$X, C = matrix(0,0,0), sp = unc$sp,
            p = sm$xp, # note the  here! This is for decreasing fits!
            y = df$y,
            w = df$y*0+1)
  G$Ain <- F1$A    # the monotonicity constraint matrix
  G$bin <- F1$b    # the monotonicity constraint vector, both from mono.con
  G$S <- sm$S     # the penalty matrix for the cubic spline
  G$off <- 0      # location of offsets in the penalty matrix
  
  ## Do the constrained fit 
  p <- pcls(G)  # fit spline (using s.p. from unconstrained fit)
  
  ## predict at 100 locations over range of x - get a smooth line on the plot
  newx <- with(df, data.frame(x = seq(min(x), max(x), length = quantiles+1)))
  
  fv <- Predict.matrix(sm, newx) %*% p
  
  newx <- transform(newx, yhat = fv[,1])
  fin<-data.frame(x=newx$yhat,cdf=newx$x,freq=c(0,diff(newx$x)),oldx=x)
  R2<-cor(fin$x,fin$oldx)^2
 # print(R2)
  red=red+1
  }
  if(R2<0.99) print(R2)
   fin<-fin %>% mutate(x=round(x,5)) %>% 
    group_by(x) %>% 
    summarize(cdf=max(cdf)) %>% 
    mutate(freq=c(cdf[1],diff(cdf))) %>% ungroup()
  fin<-as.data.frame(fin[,c(1,3,2)])
  return(fin)
}
