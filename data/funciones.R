##########################################################
# Función para calcular el retardo i de x.
##########################################################

retard<-function(x,i){
	l<-length(x)
if(i>0){
	x1<-lag(x,-i)
	x2<-c(rep(NA,i),x1)
	x2<-x2[1:l]
}
else{
	x1<-lag(x,i)
	x2<-c(x1,rep(NA,abs(i)),x1)	
	x2<-x2[abs(i)+1:l+abs(i)-1]
}
}


##########################################################
# Función QAIC
##########################################################

QAICM <- function(model,type="logLik") {
QAICm <- vector("list",0)

  if(!model$family$family%in%c("poisson","quasipoisson")) {
  	stop("only for poisson/quasipoisson family")
  }

  phi <- summary(model)$dispersion
  if(type=="dev") {
	  QAICm <- deviance(model) + 2*summary(model)$df[3]*phi
  } else {
  	loglik <- sum(dpois( model$y, model$fitted.values, log=TRUE))
  	QAICm <- -2*loglik + 2*(summary(model)$n-summary(model)$residual.df)*phi
  }

  return(QAICm)
}


# https://www.google.com/search?q=impressive-package-for-3d-and-4d-graph-r-software-and-data-visualization&rlz=1C1CHBD_esES887ES887&oq=impressibe+package&aqs=chrome.1.69i57j0i22i30j0i8i13i30.8603j1j7&sourceid=chrome&ie=UTF-8


##########################################################
# Función para calcular el intervalo de confianza del
# punto de mínima mortalidad 
##########################################################

findmin <- function(basis,model=NULL,coef=NULL,vcov=NULL,at=NULL,from=NULL,to=NULL,by=NULL,sim=FALSE,nsim=5000) {
  # CREATE THE BASIS AND EXTRACT COEF-VCOV #
  # CHECK AND DEFINE BASIS 
  if(!any(class(basis)%in%c("crossbasis","onebasis")))
    stop("the first argument must be an object of class 'crossbasis' or 'onebasis'") 
  # INFO
  one <- any(class(basis)%in%c("onebasis"))
  attr <- attributes(basis)
  range <- attr(basis,"range")
  if(is.null(by)) by <- 0.1
  lag <- if(one) c(0,0) else cb=attr(basis,"lag") 
  if(is.null(model)&&(is.null(coef)||is.null(vcov)))
    stop("At least 'model' or 'coef'-'vcov' must be provided")
  name <- deparse(substitute(basis))
  cond <- if(one) paste(name,"[[:print:]]*b[0-9]{1,2}",sep="") else
    paste(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}",sep="") 
  # SET COEF, VCOV CLASS AND LINK
  if(!is.null(model)) {
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE] 
    model.link <- dlnm:::getlink(model,model.class)
  } else model.class <- NA
  # CHECK
  if(length(coef)!=ncol(basis) || length(coef)!=dim(vcov)[1] || any(is.na(coef)) || any(is.na(vcov)))
    stop("model or coef/vcov not consistent with basis")
  # DEFINE at
  at <- dlnm:::mkat(at,from,to,by,range,lag,bylag=1) 
  predvar <- if(is.matrix(at)) rownames(at) else at 
  predlag <- dlnm:::seqlag(lag,by=1)
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON TYPE) 
  type <- if(one) "one" else "cb"
  Xpred <- dlnm:::mkXpred(type,basis,at,predvar,predlag,cen=NULL)
  Xpredall <- 0
  for(i in seq(length(predlag))) {
    ind <- seq(length(predvar))+length(predvar)*(i-1)
    Xpredall <- Xpredall + Xpred[ind,,drop=FALSE] 
  }
  # FIND THE MINIMUM
  pred <- drop(Xpredall%*%coef) 
  ind <- which.min(pred)
  min <- predvar[ind]
  # SIMULATIONS
  if(sim) {
    # SIMULATE COEFFICIENTS
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X) 
    # COMPUTE MINIMUM
    minsim <- apply(coefsim,2,function(coefi) { 
      pred <- drop(Xpredall%*%coefi)
      ind <- which.min(pred) 
      return(predvar[ind])
    })
  }
  #
  res <- if(sim) minsim else min
  #
  return(res)
}

##########################################################
# Funciones para calcular el intervalo de confianza de los
# coeficientes relativos al día de la semana
##########################################################

profile.gnm <- function (fitted, which = ofInterest(mfirst)[15:20], alpha = 0.05,
                         maxsteps = 10, stepsize = NULL, trace = FALSE, ...)
{
  fittedCoef <- parameters(fitted)
  coefNames <- names(fittedCoef)
  p <- length(coefNames)
  if (is.null(which))
    which <- 1:p
  else if (is.numeric(which))
    which <- which
  else if (is.character(which))
    which <- match(which, coefNames)
  summ <- summary(fitted)
  sterr <- summ$coefficients[, "Std. Error"]
  fittedDev <- deviance(fitted)
  disp <- summ$dispersion
  ## use z cutoffs as in confint.profile.gnm
  zmax <- abs(qnorm(alpha/2))
  fittedConstrain <- fitted$constrain
  fittedConstrainTo <- fitted$constrainTo
  auto <- is.null(stepsize)
  if (!auto)
    stepsize[1:2] <- stepsize
  prof <-  as.list(rep(NA, length(which)))
  names(prof) <- coefNames[which]
  which <- which[!is.na(sterr)[which]]
  for (i in which) {
    par <- coefNames[i]
    prof[[par]] <- numeric(2 * maxsteps + 1)
    par.vals <- matrix(nrow = 2 * maxsteps + 1, ncol = p,
                       dimnames = list(NULL, coefNames))
    par.vals[maxsteps + 1,] <- fittedCoef
    asymptote <- c(FALSE, FALSE)
    if (auto) {
      ## set defaults
      sub <- 3 # no. of steps from MLE to zmax*se
      stepsize <- c(zmax/sub * sterr[i], zmax/sub * sterr[i])
      ## estimate quadratic in the region MLE +/- zmax*se
      margin <- zmax * sterr[i]
      updatedDev <- numeric(2)
      for (sgn in c(-1, 1)) {
        val <- fittedCoef[i] + sgn * margin
        updated <-
          suppressWarnings(update(fitted, constrain =
                                    c(fittedConstrain, i),
                                  constrainTo =
                                    c(fittedConstrainTo, val),
                                  trace = FALSE, verbose = FALSE,
                                  start = fittedCoef))
        if (is.null(updated))
          break
        updatedDev[(sgn + 1)/2 + 1] <- deviance(updated)
        prof[[par]][maxsteps + 1 + sgn * sub] <-
          sgn * sqrt((deviance(updated) - fittedDev)/disp)
        par.vals[maxsteps + + 1 + sgn * sub,] <- parameters(updated)
      }
      if (all(updatedDev != 0)) {
        quad <- (sum(updatedDev) - 2 * fittedDev)/(2 * margin^2)
        lin <- (fittedDev - updatedDev[1])/margin +
          quad * (margin - 2 * fittedCoef[i])
        int <- fittedDev - lin * fittedCoef[i] - quad * fittedCoef[i]^2
        ## adjust so roots approx where deviance gives z = zmax
        int.adj <- int - zmax^2 * disp - fittedDev
        for (sgn in c(-1, 1)) {
          dir <- (sgn + 1)/2 + 1
          root <- (-lin + sgn * sqrt(lin^2 - 4 * int.adj * quad))/
            (2 * quad)
          firstApprox <- par.vals[maxsteps + 1 + sgn * sub, i]
          ## if likelihood approx quadratic use default stepsize, else
          if (sgn * (root - firstApprox) > 0) {
            ## not gone out far enough, check for asymptote
            val <- fittedCoef[i] + sgn * 10 * sterr[i]
            updated <-
              suppressWarnings(update(fitted, constrain =
                                        c(fittedConstrain, i),
                                      constrainTo =
                                        c(fittedConstrainTo, val),
                                      trace = FALSE,
                                      verbose = FALSE,
                                      start = fittedCoef))
            if (!is.null(updated) &&
                sqrt((deviance(updated) - fittedDev)/disp) < zmax)
              asymptote[dir] <- TRUE
          }
          ## if root more than one step away from firstApprox, i.e.
          ## less than two steps away from fittedCoef, halve stepsize
          if (abs(sgn * (firstApprox - root)) > stepsize[dir] &&
              !asymptote[dir]) {
            prof[[par]][maxsteps + 1 + sgn * sub] <- 0
            par.vals[maxsteps + 1 + sgn * sub, ] <- NA
            stepsize[dir] <- abs(root - fittedCoef[i])/(maxsteps/2)
          }
        }
      }
    }
    for (sgn in c(-1, 1)) {
      if (trace)
        prattle("\nParameter:", par, c("down", "up")[(sgn + 1)/2 + 1],
                "\n")
      step <- 0
      init <- parameters(fitted)
      while ((step <- step + 1) <= maxsteps) {
        if (step > 2 &&
            abs(prof[[par]][maxsteps + 1 + sgn * (step - 2)]) > zmax)
          break
        if (prof[[par]][maxsteps + 1 + sgn * step] != 0)
          next
        val <- fittedCoef[i] + sgn * step * stepsize[(sgn + 1)/2 + 1]
        updated <-
          suppressWarnings(update(fitted, constrain =
                                    c(fittedConstrain, i),
                                  constrainTo =
                                    c(fittedConstrainTo, val),
                                  trace = FALSE, verbose = FALSE,
                                  start = init))
        if (is.null(updated)) {
          message("Could not complete profile for", par, "\n")
          break
        }
        init <- parameters(updated)
        zz <- (deviance(updated) - fittedDev)/disp
        if (zz > -0.001)
          zz <- max(zz, 0)
        else stop("profiling has found a better solution, ",
                  "so original fit had not converged")
        prof[[par]][maxsteps + 1 + sgn * step] <- sgn * sqrt(zz)
        par.vals[maxsteps + 1 + sgn * step,] <- init
        #print(data.frame(step = step, val = bi, deviance = fm$deviance,
        #zstat = z))
      }
    }
    prof[[par]] <- structure(data.frame(prof[[par]][!is.na(par.vals[,1])]),
                             names = "z")
    prof[[par]]$par.vals <- par.vals[!is.na(par.vals[,1]), , drop = FALSE]
    attr(prof[[par]], "asymptote") <- asymptote
  }
  val <- structure(prof, original.fit = fitted, summary = summ)
  class(val) <- c("profile.gnm", "profile.glm", "profile")
  val
}

confint.gnm.diasemana <- function (object, parm = ofInterest(mfirst)[15:20], level = 0.95,
                                   trace = FALSE, ...)
{
  pnames <- names(coef(object))
  if (is.null(parm))
    parm <- seq(along = pnames)
  else if (is.character(parm))
    parm <- match(parm, pnames, nomatch = 0)
  message("Waiting for profiling to be done...")
  flush.console()
  object <- profile.gnm(object, which = parm, alpha = 1 - level, trace = trace)
  confint(object, level = level, ...)
}
