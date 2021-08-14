

twoRateModel <- function(par, schedule) {
  
  # thse values should be zero at the start of the loop:
  Et <- 0 # previous error: none
  St <- 0 # state of the slow process: aligned
  Ft <- 0 # state of the fast process: aligned
  
  # we'll store what happens on each trial in these vectors:
  slow <- c()
  fast <- c()
  total <- c()
  
  # now we loop through the perturbations in the schedule:
  for (t in c(1:length(schedule))) {
    
    # first we calculate what the model does
    # this happens before we get visual feedback about potential errors
    St <- (par['Rs'] * St) - (par['Ls'] * Et)
    Ft <- (par['Rf'] * Ft) - (par['Lf'] * Et)
    Xt <- St + Ft
    
    # now we calculate what the previous error will be for the next trial:
    if (is.na(schedule[t])) {
      Et <- 0
    } else {
      Et <- Xt + schedule[t]
    }
    
    # at this point we save the states in our vectors:
    slow <- c(slow, St)
    fast <- c(fast, Ft)
    total <- c(total, Xt)
    
  }
  
  # after we loop through all trials, we return the model output:
  return(data.frame(slow,fast,total))
  
}


twoRateMSE <- function(par, schedule, reaches, checkStability=FALSE) {
  
  bigError <- mean(schedule^2, na.rm=TRUE) * 10
  
  # learning and retention rates of the fast and slow process are constrained:
  if (par['Ls'] > par['Lf']) {
    return(bigError)
  }
  if (par['Rs'] < par['Rf']) {
    return(bigError)
  }


  if (checkStability) {
    aa <- ((par['Rf'] - par['Lf']) * (par['Rs'] - par['Ls'])) - (par['Lf'] * par['Ls'])
    if (aa <= 0) {
      return(bigError)
    }
    
    p <- par['Rf'] - par['Lf'] - par['Rs'] + par['Ls']
    q <- p^2 + (4 * par['Lf'] * par['Ls'])
    bb <- ((par['Rf'] - par['Lf'] + par['Rs'] - par['Ls'])  +  sqrt(q))
    if (bb >= 2) {
      return(bigError)
    }
    
  }
  
  return( mean((twoRateModel(par, schedule)$total - reaches)^2, na.rm=TRUE) )
  
}


twoRateFit <- function(schedule, reaches, gridpoints=6, gridfits=6, checkStability=FALSE) {
  
  parvals <- seq(1/gridpoints/2,1-(1/gridpoints/2),1/gridpoints)
  
  searchgrid <- expand.grid('Ls'=parvals,
                            'Lf'=parvals,
                            'Rs'=parvals,
                            'Rf'=parvals)
  # evaluate starting positions:
  MSE <- apply(searchgrid, FUN=twoRateMSE, MARGIN=c(1), schedule=schedule, reaches=reaches, checkStability=checkStability)
  
  optimxInstalled <- require("optimx")
  
  if (optimxInstalled) {
    
    # run optimx on the best starting positions:
    allfits <- do.call("rbind",
                       apply( searchgrid[order(MSE)[1:gridfits],],
                              MARGIN=c(1),
                              FUN=optimx::optimx,
                              fn=twoRateMSE,
                              method='L-BFGS-B',
                              lower=c(0,0,0,0),
                              upper=c(1,1,1,1),
                              schedule=schedule,
                              reaches=reaches,
                              checkStability=checkStability) )
    
    # pick the best fit:
    win <- allfits[order(allfits$value)[1],]
    
    # return the best parameters:
    return(unlist(win[1:4]))
    
  } else {
    
    cat('(consider installing optimx, falling back on optim now)\n')
    
    # use optim with Nelder-Mead after all:
    allfits <- do.call("rbind",
                       apply( searchgrid[order(MSE)[1:gridfits],],
                              MARGIN=c(1),
                              FUN=stats::optim,
                              fn=twoRateMSE,
                              method='Nelder-Mead',
                              schedule=schedule,
                              reaches=reaches,
                              checkStability=checkStability) )
    
    # pick the best fit:
    win <- allfits[order(unlist(data.frame(allfits)[,'value']))[1],]
    
    # return the best parameters:
    return(win$par)
    
  }
  
}



oneRateModel <- function(par, schedule) {
  
  # thse values should be zero at the start of the loop:
  Et <- 0 # previous error: none
  Pt <- 0 # state of slow process: aligned

  # we'll store what happens on each trial in these vectors:
  process <- c()

  # now we loop through the perturbations in the schedule:
  for (t in c(1:length(schedule))) {
    
    # first we calculate what the model does
    # this happens before we get visual feedback about potential errors
    Pt <- (par['R'] * Pt) - (par['L'] * Et)

    # now we calculate what the previous error will be for the next trial:
    if (is.na(schedule[t])) {
      Et <- 0
    } else {
      Et <- Pt + schedule[t]
    }
    
    # at this point we save the states in our vectors:
    process <- c(process, Pt)

  }
  
  # after we loop through all trials, we return the model output:
  return(data.frame(process))
  
}





oneRateMSE <- function(par, schedule, reaches) {
  
  bigError <- mean(schedule^2, na.rm=TRUE) * 10
  
  return( mean((oneRateModel(par, schedule)$process - reaches)^2, na.rm=TRUE) )
  
}



oneRateFit <- function(schedule, reaches, gridpoints=6, gridfits=6) {
  
  parvals <- seq(1/gridpoints/2,1-(1/gridpoints/2),1/gridpoints)
  
  searchgrid <- expand.grid('L'=parvals,
                            'R'=parvals)
  # evaluate starting positions:
  MSE <- apply(searchgrid, FUN=oneRateMSE, MARGIN=c(1), schedule=schedule, reaches=reaches)
  
  optimxInstalled <- require("optimx")
  
  if (optimxInstalled) {
    
    # run optimx on the best starting positions:
    allfits <- do.call("rbind",
                       apply( searchgrid[order(MSE)[1:gridfits],],
                              MARGIN=c(1),
                              FUN=optimx::optimx,
                              fn=oneRateMSE,
                              method='L-BFGS-B',
                              lower=c(0,0),
                              upper=c(1,1),
                              schedule=schedule,
                              reaches=reaches ) )
    
    # pick the best fit:
    win <- allfits[order(allfits$value)[1],]
    
    # return the best parameters:
    return(unlist(win[1:2]))
    
  } else {
    
    cat('(consider installing optimx, falling back on optim now)\n')
    
    # use optim with Nelder-Mead after all:
    allfits <- do.call("rbind",
                       apply( searchgrid[order(MSE)[1:gridfits],],
                              MARGIN=c(1),
                              FUN=stats::optim,
                              fn=oneRateMSE,
                              method='Nelder-Mead',
                              schedule=schedule,
                              reaches=reaches ) )
    
    # pick the best fit:
    win <- allfits[order(unlist(data.frame(allfits)[,'value']))[1],]
    
    # return the best parameters:
    return(win$par)
    
  }
  
}


twoRates <- function(groupdata, estN='ac_one', compOne=FALSE) {
  
  groupnames <- names(groupdata)
  
  if (is.numeric(estN)) {
    
    if (!all(groupnames %in% names(estN))) {
      stop('Not all groups in the data have an N specified.\n')
    } else {
      # all is good
      # we put this in a vector that aclag will create:
      observations <- estN
    }


    # all should be good now...
  } else if (is.character(estN)) {
    
    observations <- c()
    
    # determine number of independent observations for each group:
    for (group in groupnames) {
      
      # use reaches to determine number of independent observations:
      reaches <- groupdata[group]$reaches
      
      observations[group] <- seriesEffectiveSampleSize(reaches, method=estN)
    
    }
    
  } else {
    stop('Either specify "estN" as a numeric vector of N or a character naming an available method.\n')
  }
  
  # now do the actual modelling on all datasets?
  
  
  
  
}



seriesEffectiveSampleSize <- function(series, method='ac_one') {
  
  # https://dsp.stackexchange.com/questions/47560/finding-number-of-independent-samples-using-autocorrelation
  # The "Effective" Number of Independent Observations in an An′ = n * ( (1−ρ)/(1+ρ) )utocorrelated Time Series is a defined statistical term - https://www.jstor.org/stable/2983560?seq=1#page_scan_tab_contents
  # The number of independent observations n' of n observations with a constant variance but having a lag 1 autocorrelation ρ equals
  #
  # n′ = n * ( (1−ρ)/(1+ρ) )
  #
  # Also note this is an approximation valid for large n, [1] (reference provided by Ed V) equation 7 is more accurate for small n.
  # [1] N.F. Zhang, "Calculation of the uncertainty of the mean of autocorrelated measurements", Metrologia 43 (2006) S276-S281.
  
  
  
  # Maybe also relevant:
  # https://www.researchgate.net/post/How_do_you_choose_the_optimal_laglength_in_a_time_series
  
  
  
  if (method == 'ac_one') {
  
    # create empty vector for number of independent observations per group:
    observations <- c()
  
    rho <- acf(series, lag.max=1, plot=FALSE)$acf[2]
    
    return( length(series) *  ((1-rho)/(1+rho)) )
    
  } else if (method == 'ac_lag.10') {
  
    critlag <- which(stats::acf(series, lag.max=length(series)-1, plot=FALSE)$acf < 0.1)
    
    if (length(critlag) == 0) {
      
      # autocorrelation is high throughout: we have 1 observation?
      
      return(1)
      
    } else {
      
      # the sequence of reaches doesn't autocorrelate well throughout,
      # so it can be split into (more or less) independent observations
      # (but we have at least 1)
      
      return( max( ( length(series) / (critlag[1]-2) ), 1) )
      
    }
    
  } else if (method == 'ac_lag95%CI') {
    
    Neff_found <- FALSE
    
    critlag <- 1
    
    Neff <- length(series)
    
    while (!Neff_found) {
      
      lagpoints <- length(series) - critlag
      
      point_one <- series[c(1:lagpoints)]
      point_two <- series[c((critlag+1):length(series))]
      
      lag_cor <- stats::cor(point_one, point_two)
      
      shuffle_cor <- rep(NA, 1000)
      
      for (bootstrap in c(1:1000)) {
        
        shuffle_cor[bootstrap] <- stats::cor(point_one, sample(point_two))
        
      }
      
      upperlimit <- stats::quantile(shuffle_cor, probs=0.95)
      
      if (lag_cor < upperlimit) {
        
        return( length(series) / max((critlag - 1), 1) )
        
      }
      
      critlag <- critlag + 1
      
      # lag can only go up to a certain value, determined by the length of the sequence
      # make sure there are at least 3 data points for acf...
      if (critlag > (length(series) - 2)) {
        
        return( 1 ) # or length(series) - 1?
        
      }
      
    }
    
  }
  
  stop('Unrecognized method for determining effective sample size.\nUse one of: ac_one, ac_lag.10 or ac_lag95%CI\n')
  
}


modelCriteriaMSE <- function(MSE, k, N, n=NA, MSEmean) {
  
  # MSE : our goodness of fit measure in lieu of actual likelihood
  # k   : number of parameters
  # N   : number of independent observations
  # n   : number of observations (number of trials for two-rate models)

  if (length(MSE) != length(k)) {
    stop('Arguments MSE and k need to be of the same length.\n')
  }
  if (any(c(length(MSE), length(k), length(N)) < 1)) {
    stop('All arguments must be at least of length 1.\n')
  }
  if (!is.numeric(MSE) | !is.numeric(k)) {
    stop('MSE and k need to be numeric.\n')
  }
  if (length(N) > 1 | !is.numeric(N)) {
    stop('N has to be a single numeric value.\n')
  }
  if (!is.na(n) && (length(n) > 1 | !is.numeric(n))) {
    stop('n has to be NA, or a single numeric value.\n')
  }
  
  # maximum likelihood:
  # L <- (-(n/2) * log(2*pi)) - ((n/2)*log(MSE)) - (1/(2*MSE)*(MSE*n))
  # but this sometimes results in negative likelihoods
  # the log_e of which causes problems later on 
  
  # n <- N
  
  # without the "constant" that Wikipedia mentions:
  # this is simpler, and I might replace the constant
  # L <- -(n/2) * log(MSE)
  
  # sometimes we now get inf or nan output, 
  # so we replace the constant to avoid this:
  # if (any(L < 1)) {
  #   L <- (L - min(L)) + 1
  # }
  
  #-- AIC --# 
  
  # Thomas calculation:
  # C <- N*(log(2*pi)+1) # what is this for? a penalty for large number of observations?
  # AIC <- (2 * k) + N*log(MSE) + C
  
  AIC <- (N * log(MSE)) + (2 * k)
  # AIC <- (2 * k) - (N * log(L))
  
  #-- AICc --#
  
  if (!is.na(n)) {
    # correction for low N (compared to k):
    AICc <- AIC + ( (2 * k^2) / (n - k - 1) )
  }
  
  #-- BIC --#
  
  #BIC <- log(N)*k - (2 * log(L))
  
  #-- Hannan-Quinn --#
  
  #HQC <- (-2 * L) + (2 * k * log(log(N)))
  
  if (length(MSE) == 1) {
    
    # return(data.frame('AIC'=AIC, 'AICc'=AICc, 'BIC'=BIC, 'HQC'=HQC))
    if (is.na(n)) {
      return(data.frame('AIC'=AIC))
    } else {
      return(data.frame('AIC'=AIC, 'AICc'=AICc))
    }
    
  } else {
    
    AIC.rl  <- relativeLikelihood( AIC  ) # exp( ( min( AIC  ) - AIC  ) / 2 )
    
    if (is.na(n)) {
      return(data.frame('AIC'=AIC,'AIC.rl'=AIC.rl))
    } else {
      AICc.rl <- relativeLikelihood( AICc )
      return(data.frame('AIC'=AIC,'AIC.rl'=AIC.rl, 'AICc'=AICc, 'AICc.rl'=AICc.rl))
    }
    # BIC.rl  <- relativeLikelihood( BIC  )
    # HQC.rl  <- relativeLikelihood( HQC  )
    
    # return(data.frame('AIC'=AIC,'AIC.rl'=AIC.rl, 'AICc'=AICc, 'AICc.rl'=AICc.rl, 'BIC'=BIC, 'BIC.rl'=BIC.rl, 'HQC'=HQC, 'HQC.rl'=HQC.rl))
    
  }
  
}


relativeLikelihood <- function(crit) {
  
  return( exp( ( min( crit  ) - crit  ) / 2 ) )
  
}

modelCriteriaLikelihood <- function(MSE, k, N, n=NA, MSEmean=NA) {
  
  # MSE : our goodness of fit measure in lieu of actual likelihood
  # k   : number of parameters
  # N   : number of independent observations
  # n   : number of observations (number of trials for two-rate models)
  
  if (length(MSE) != length(k)) {
    stop('Arguments MSE and k need to be of the same length.\n')
  }
  if (any(c(length(MSE), length(k), length(N)) < 1)) {
    stop('All arguments must be at least of length 1.\n')
  }
  if (!is.numeric(MSE) | !is.numeric(k)) {
    stop('MSE and k need to be numeric.\n')
  }
  if (length(N) > 1 | !is.numeric(N)) {
    stop('N has to be a single numeric value.\n')
  }
  if (!is.na(n) && (length(n) > 1 | !is.numeric(n))) {
    stop('n has to be NA, or a single numeric value.\n')
  }
  
  # maximum likelihood:
  # L <- (-(n/2) * log(2*pi)) - ((n/2)*log(MSE)) - (1/(2*MSE)*(MSE*n))
  # but this sometimes results in negative likelihoods
  # the log_e of which causes problems later on 
  
  # n <- N
  
  # without the "constant" that Wikipedia mentions:
  # this is simpler, and I might replace the constant
  # print(-(N/2))
  # print(MSE)
  # print(log(MSE))
  
  
  
  L <- -(N/2) * log(MSE) # wikipedia
  #print(L)
  
  # sometimes we now get inf or nan output, 
  # so we replace the constant to avoid this:
  # if (any(L < 1)) {
  #   L <- (L - min(L)) + 1
  # }
  
  # R2 <- 1 - (MSE/MSEmean)
  # print(R2)
  # L <- R2
  
  #-- AIC --# 
  
  # Thomas calculation:
  # C <- N*(log(2*pi)+1) # what is this for? a penalty for large number of observations?
  # AIC <- (2 * k) + N*log(MSE) + C
  
  # AIC <- (N * log(MSE)) + (2 * k) # MSE based
  AIC <- (2 * k) - (N * log(L)) # L based?
  
  #-- AICc --#
  
  if (!is.na(n)) {
    # correction for low N (compared to k):
    AICc <- AIC + ( (2 * k^2) / (n - k - 1) )
  }
  
  #-- BIC --#
  
  BIC <- log(N)*k - (2 * log(L)) # L based
  
  #-- Hannan-Quinn --#
  
  HQC <- (-2 * L) + (2 * k * log(log(N))) # L based
  
  if (length(MSE) == 1) {
    
    # return(data.frame('AIC'=AIC, 'AICc'=AICc, 'BIC'=BIC, 'HQC'=HQC))
    if (is.na(n)) {
      return(data.frame('AIC'=AIC))
    } else {
      return(data.frame('AIC'=AIC, 'AICc'=AICc))
    }
    
  } else {
    
    AIC.rl  <- relativeLikelihood( AIC  ) # exp( ( min( AIC  ) - AIC  ) / 2 )
    BIC.rl  <- relativeLikelihood( BIC  )
    HQC.rl  <- relativeLikelihood( HQC  )
    
    if (is.na(n)) {
      return(data.frame('AIC'=AIC,'AIC.rl'=AIC.rl, 'BIC'=BIC, 'BIC.rl'=BIC.rl, 'HQC'=HQC, 'HQC.rl'=HQC.rl))
    } else {
      AICc.rl <- relativeLikelihood( AICc )
      return(data.frame('AIC'=AIC,'AIC.rl'=AIC.rl, 'AICc'=AICc, 'AICc.rl'=AICc.rl, 'BIC'=BIC, 'BIC.rl'=BIC.rl, 'HQC'=HQC, 'HQC.rl'=HQC.rl))
    }

    # return(data.frame('AIC'=AIC,'AIC.rl'=AIC.rl, 'AICc'=AICc, 'AICc.rl'=AICc.rl, 'BIC'=BIC, 'BIC.rl'=BIC.rl, 'HQC'=HQC, 'HQC.rl'=HQC.rl))
    
  }
  
}

