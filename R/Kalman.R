#' Kalman
#'
#' Takes list containing ctsem subject matrices, as well as long form data object, and calculates 
#' predicted and updated latent states, likelihoods, and predicted observations using the Kalman filter.
#'
#' @param kpars list object containing DRIFT,T0VAR,DIFFUSION,CINT,T0MEANS,TDPREDEFFECT,
#' MANIFESTMEANS, LAMBDA, and MANIFESTVAR matrices, with list elements named accordingly. 
#' Such a list is returned by \code{\link{ctStanContinuousPars}}.
#' @param datalong long format data object as used by \code{\link{ctStanFit}}, 
#' but must contain only a single subjects' data and does not need an id column.
#' @param manifestNames String vector of names of manifest varifables to use from datalong.
#' @param latentNames String vector of names of latent variables.
#' @param TDpredNames If model contains time dependent predictors, 
#' string vector of their names in the data.
#' @param imputeMissings Logical. If TRUE, randomly generate any missing observations of 
#' manifest variables according to model.
#' @param idcol Character string giving name of subject identification column in data.
#' @param continuoustime Logical, whether to use a continuous time Kalman filter or discrete time. 
#' Refers only to latent states, observations are always at discrete time points.
#' @param timecol name of time column in datalong. Note that time column must be an ascending sequence
#' of numeric values from row 1 to row n. Ignored if continuoustime=FALSE.
#' @param derrind vector of integers denoting which latent variables are involved in covariance calcs.
#' @param optimize Set to TRUE when using for optimization.
#' @param ukf set to TRUE to use the unscented Kalman filter, only necessary for fitting non-linear models, 
#' currently only for optimizing.
#' @param plotoptim set to TRUE to plot / print optimization steps.
#' @return When optimize=TRUE, returns log likelihood. Else, 
#' returns a list containing matrix objects etaprior, etaupd, etasmooth, y, yprior, 
#' yupd, ysmooth, prederror, time, loglik,  with values for each time point in each row. 
#' eta refers to latent states and y to manifest indicators - y itself is thus just 
#' the input data. 
#' Covariance matrices etapriorcov, etaupdcov, etasmoothcov, ypriorcov, yupdcov, ysmoothcov,  
#' are returned in a row * column * time array. 
#' @examples
#' \donttest{
#' ### get parameter matrices
#'
#' kpars <- ctStanContinuousPars(ctstantestfit)
#' 
#' #construct dummy data
#' datalong <- cbind(0:9, 1, matrix(rnorm(20,2,1),ncol=2))
#' datalong[c(1:3,9:10),3:4]<-NA #missing data to pre/fore cast
#' colnames(datalong) <- c('time', 'id', paste0('Y',1:2))
#' print(datalong)
#' 
#' #obtain Kalman filtered estimates
#' kout <- Kalman(kpars=kpars, datalong=datalong,
#'   manifestNames=paste0('Y',1:nrow(kpars$MANIFESTMEANS)),
#'   latentNames=paste0('eta',1:nrow(kpars$DRIFT)))
#' 
#' #print and plot smoothed estimates (conditional on all states) of indicators.
#' print(kout$ysmooth)
#' matplot(kout$time,kout$ysmooth,type='l')
#' matplot(kout$time,datalong[,3:4],type='p',add=TRUE,pch=1)
#' }
#' @export

Kalman<-function(kpars,datalong,
  manifestNames,latentNames,imputeMissings=FALSE,
  TDpredNames=NULL,
  continuoustime=TRUE,idcol='id',
  timecol='time', derrind='all',optimize=FALSE,ukf=FALSE, plotoptim=FALSE){
  
  datalong=as.matrix(datalong)
  
  nmanifest=length(manifestNames)
  nlatent=length(latentNames)
  ntdpred=length(TDpredNames)
  
  if(any(c(nmanifest,nlatent) < 1)) stop('Length of manifestNames and latentNames must be greater than 0!')
  if(all(derrind=='all')) derrind=1:nlatent
  ndiffusion=length(derrind)
  
  Y<-datalong[,manifestNames,drop=FALSE]
  if(ntdpred > 0) {
    tdpreds<-datalong[,TDpredNames,drop=FALSE]
    tdpreds[is.na(tdpreds)] <- 0
    # if(any(is.na(tdpreds))) stop('missingness in time dependent predictors! Kalman cannot run.')
  }
  
  if(continuoustime) {
    DRIFTHATCH <- (kpars$DRIFT[derrind,derrind] %x% diag(ndiffusion) + 
        diag(ndiffusion) %x% kpars$DRIFT[derrind,derrind,drop=FALSE]) 
    asymDIFFUSION <- matrix(-solve(DRIFTHATCH, c(kpars$DIFFUSION[derrind,derrind,drop=FALSE] + 
        diag(1e-5,nlatent)[derrind,derrind,drop=FALSE])), nrow=ndiffusion)
  }
  
  etaprior<-list()
  etapriorcov<-list()
  etaupd<-list()
  etaupdcov<-list()
  err<-list()
  yprior<-list()
  ypriorcov<-list()
  yupd<-list()
  yupdcov<-list()
  
  
  loglik<-rep(0,nrow(datalong))
  observed<-list()
  
  discreteDRIFT<- list()
  discreteCINT<- list()
  discreteDIFFUSION <- array(0,dim=c(nlatent,nlatent,nrow(datalong)))
  
  t0check <- c(1,as.numeric(datalong[-1,idcol]!=datalong[-nrow(datalong),idcol]))
  
  if(continuoustime) dt<-c(-9,datalong[-1,timecol]-datalong[-nrow(datalong),timecol])
  
  if(!ukf){
    for(rowi in 1:(nrow(datalong))){
      
      if(continuoustime){
        discreteDRIFT[[rowi]] <- as.matrix(Matrix::expm(kpars$DRIFT * dt[rowi]))
        discreteCINT[[rowi]] <- solve(kpars$DRIFT, (discreteDRIFT[[rowi]] - diag(nlatent))) %*% kpars$CINT
        discreteDIFFUSION[,,rowi] <- asymDIFFUSION - (discreteDRIFT[[rowi]][derrind,derrind,drop=FALSE] %*% 
            asymDIFFUSION %*% t(discreteDRIFT[[rowi]][derrind,derrind,drop=FALSE]))
      }
      if(!continuoustime){
        discreteDRIFT[[rowi]] <-kpars$DRIFT
        discreteCINT[[rowi]]<- kpars$CINT
        discreteDIFFUSION[derrind,derrind,rowi] <- kpars$DIFFUSION[derrind,derrind,drop=FALSE]
      }
      
      if(t0check[rowi]==1){
        etaprior[[rowi]]<-kpars$T0MEANS #tdpreds added below
        etapriorcov[[rowi]]<-kpars$T0VAR
      }
      
      
      if(t0check[rowi]==0){
        etaprior[[rowi]] <- discreteCINT[[rowi]]  + discreteDRIFT[[rowi]] %*% etaupd[[rowi-1]]
        etapriorcov[[rowi]] <-  discreteDRIFT[[rowi]] %*% 
          etaupdcov[[rowi-1]] %*% t(discreteDRIFT[[rowi]])  + discreteDIFFUSION[derrind,derrind,rowi] #check transpose
      }
      
      if(ntdpred > 0) etaprior[[rowi]] <- etaprior[[rowi]] + kpars$TDPREDEFFECT %*% t(datalong[rowi,TDpredNames,drop=FALSE])
      
      if(imputeMissings) Y[rowi,] <- 0 #etaprior[[rowi]] + t(chol(etapriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
      
      nafilter<-!is.na(Y[rowi,])
      observed[[rowi]]<-nafilter
      err[[rowi]]<-matrix(NA,nrow=nmanifest) #init prediction errors for this row
      
      # // one step ahead predictive distribution of y
      yprior[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etaprior[[rowi]]
      ypriorcov[[rowi]] <- kpars$LAMBDA %*% etapriorcov[[rowi]] %*% 
        t(kpars$LAMBDA) + kpars$MANIFESTVAR
      
      if(imputeMissings) Y[rowi,] <- yprior[[rowi]] + t(chol(ypriorcov[[rowi]])) %*% rnorm(nmanifest,0,1)
      
      y <- Y[rowi,,drop=FALSE][,nafilter,drop=FALSE]
      
      #if all missing...
      if(all(!nafilter)){
        etaupd[[rowi]] <- etaprior[[rowi]]
        etaupdcov[[rowi]] <- etapriorcov[[rowi]]
        yupd[[rowi]] <- yprior[[rowi]]
        yupdcov[[rowi]] <- ypriorcov[[rowi]]
      }
      
      #if any not missing
      if(any(nafilter)){
        
        # // forecast error
        err[[rowi]][nafilter] <- as.numeric(y - yprior[[rowi]][nafilter,])
        
        # // Kalman gain
        # invypriorcov <- solve(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE])
        K<-matrix(0,nrow=nlatent,ncol=nmanifest)
        
        K[,nafilter] <-  t(solve(t(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE]),
          t(etapriorcov[[rowi]] %*%  t(kpars$LAMBDA[nafilter,,drop=FALSE]) ) ))
        
        # updated distribution 
        etaupd[[rowi]] <- etaprior[[rowi]] + K[,nafilter,drop=FALSE] %*% (err[[rowi]][nafilter,,drop=FALSE])
        
        # etaupdcov[[rowi + 1]] <- etapriorcov[[rowi]]
        etaupdcov[[rowi]] <- (diag(ndiffusion) - K[,nafilter,drop=FALSE] %*% 
            kpars$LAMBDA[nafilter,,drop=FALSE]) %*% etapriorcov[[rowi]]
        
        # // log likelihood
        loglik[rowi] <- - 0.5 * (nrow(kpars$LAMBDA[nafilter,,drop=FALSE]) * log(2 * pi)  + 
            log(det(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE]))    + 
            t(err[[rowi]][nafilter,,drop=FALSE]) %*% 
            solve(ypriorcov[[rowi]][nafilter,nafilter,drop=FALSE],err[[rowi]][nafilter,,drop=FALSE]))
        
        yupd[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etaupd[[rowi]]
        yupdcov[[rowi]] <- kpars$LAMBDA[,,drop=FALSE] %*% etaupdcov[[rowi]] %*% 
          t(kpars$LAMBDA[,,drop=FALSE]) + kpars$MANIFESTVAR
      }
      
    }
  }
  

  if(plotoptim && runif(1) > as.numeric(plotoptim)){
    yprior2 <- array(unlist(yprior),dim=c(length(yprior[[1]]),1,nrow(Y)))
    plot(apply(Y,1,mean),type='l')
    # points(1 / (1+exp(-yprior2[1,2,])),type='b',col='blue')
    points(apply(matrix(yprior2[,1,],nrow=dim(yprior2)[1]),2,mean),type='l',col='red')
  }
  
  #extra output
  if(!optimize){
    
    #smoother
    
    etasmooth<-list()
    etasmoothcov<-list()
    ysmooth<-list()
    ysmoothcov<-list()
    
    for(rowi in nrow(datalong):1){
      if(rowi==nrow(datalong)) {
        etasmooth[[rowi]]<-etaupd[[rowi]]
        etasmoothcov[[rowi]]<-etaupdcov[[rowi]]
      } else{
        smoother<-diag(0,nlatent)
        smoother[derrind,derrind]<- t(solve(t(etapriorcov[[rowi+1]]), t( etaupdcov[[rowi]] %*% 
            t(discreteDRIFT[[rowi+1]][derrind,derrind,drop=FALSE]) ) )) #is the rowi+1 correct?
        
        etasmooth[[rowi]]<-etaupd[[rowi]]+smoother %*% (etasmooth[[rowi+1]] - etaprior[[rowi+1]])
        etasmoothcov[[rowi]]<-etaupdcov[[rowi]] + smoother[derrind,derrind,drop=FALSE] %*% 
          ( etasmoothcov[[rowi+1]] - etapriorcov[[rowi+1]])
      }
      ysmooth[[rowi]] <- kpars$MANIFESTMEANS + kpars$LAMBDA %*% etasmooth[[rowi]]
      ysmoothcov[[rowi]] <- kpars$LAMBDA[,derrind,drop=FALSE] %*% etasmoothcov[[rowi]] %*% 
        t(kpars$LAMBDA[,derrind,drop=FALSE]) + kpars$MANIFESTVAR
    }
    
    timedims=paste0('t',datalong[,timecol])
    
    etaprior<-matrix(unlist(etaprior),byrow=T,ncol=nlatent,
      dimnames=list(timedims,latentNames))
    
    etaupd<-matrix(unlist(etaupd),byrow=T,ncol=nlatent,
      dimnames=list(timedims,latentNames))
    
    etasmooth<-matrix(unlist(etasmooth),byrow=T,ncol=nlatent,
      dimnames=list(timedims,latentNames))
    
    yprior<-matrix(unlist(yprior),byrow=T,ncol=nmanifest,
      dimnames=list(timedims,manifestNames))
    
    yupd<-matrix(unlist(yupd),byrow=T,ncol=nmanifest,
      dimnames=list(timedims,manifestNames))
    
    y<-Y #matrix(datalong[,manifestNames,drop=FALSE],ncol=nmanifest) #
    colnames(y)<-manifestNames
    
    if(ntdpred>0) {
      tdpreds<-datalong[,TDpredNames,drop=FALSE]
      colnames(tdpreds)<-TDpredNames
    }
    
    err<-matrix(unlist(err),byrow=T,ncol=nmanifest,
      dimnames=list(timedims,manifestNames))
    
    ysmooth<-matrix(unlist(ysmooth),byrow=T,ncol=nmanifest,
      dimnames=list(timedims,manifestNames))
    
    
    insertdiffusionzeros <- function(x){
      y <- array(0,c(nlatent,nlatent,length(timedims)))
      y[derrind,derrind,] <- x
      dimnames(y) <- list(latentNames,latentNames,timedims)
      return(y)
    }
    
    etapriorcov <- plyr::alply(etapriorcov,1,function(x) x,.dims=TRUE)
    dimnames(etapriorcov)
    
    
    ypriorcov <- array(unlist(ypriorcov),
      dim=c(nmanifest,nmanifest,length(timedims)),
      dimnames=list(manifestNames,manifestNames,timedims))
    
    yupdcov <- array(unlist(yupdcov),
      dim=c(nmanifest,nmanifest,length(timedims)),
      dimnames=list(manifestNames,manifestNames,timedims))
    
    ysmoothcov <- array(unlist(ysmoothcov),
      dim=c(nmanifest,nmanifest,length(timedims)),
      dimnames=list(manifestNames,manifestNames,timedims))
    
    etapriorcov <- insertdiffusionzeros(array(unlist(etapriorcov),
      dim=c(length(derrind),length(derrind),length(timedims))))
    
    etaupdcov <- insertdiffusionzeros(array(unlist(etaupdcov),
      dim=c(length(derrind),length(derrind),length(timedims))))
    
    etasmoothcov <- insertdiffusionzeros(array(unlist(etasmoothcov),
      dim=c(length(derrind),length(derrind),length(timedims))))
    
    names(loglik) = timedims
    
    out<-list(observed,
      etaprior=etaprior,etapriorcov=etapriorcov,
      etaupd=etaupd,etaupdcov=etaupdcov,
      loglik=loglik, 
      prederror=err,y=y,
      yprior=yprior,ypriorcov=ypriorcov,
      yupd=yupd,yupdcov=yupdcov,
      etasmooth=etasmooth, etasmoothcov=etasmoothcov, 
      ysmooth=ysmooth, ysmoothcov=ysmoothcov)
    
    if(ntdpred > 0) out$tdpreds=tdpreds
    if(continuoustime) out$time=datalong[,timecol]
  }
  
  if(optimize) out <- sum(loglik,na.rm=TRUE)
  
  return(out)
}



