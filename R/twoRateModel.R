

# Two-Rate Model -----

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


twoRateFit <- function(schedule, reaches, gridpoints=9, gridfits=5, checkStability=FALSE) {
  
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

# One-Rate model -----

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


# model fitting extras -----


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

AIC <- function(MSE, k, N) {

  return( (N * log(MSE)) + (2 * k) )
  
}  

AICc <- function(MSE, k, N) {
  
  AIC <- AIC(MSE, k, N)
  
  return( AIC + ( (2 * k^2) / (n - k - 1) ) )
  
}

relativeLikelihood <- function(crit) {
  
  return( exp( ( min( crit  ) - crit  ) / 2 ) )
  
}

# bootstrapping -----

bootstrapFits <- function(exp=NULL, groups=NULL, iterations=1000) {
  
  info <- getInfo()
  
  if (!is.null(exp) & is.null(groups)) {
    groups <- info$group[which(info$experiment == exp)]
  }
  
  alldata <- getSelectedGroupsData(groups=groups)
  alldata <- parseData(alldata, FUN=baseline)
  alldata <- parseData(alldata, FUN=addtrial)
  
  for (group in groups) {
    
    data <- alldata[[group]]
    
    # bootstrapping is based on number of participants
    participants <- unique(data[['abrupt']]$participant)
    
    # get group schedules
    schedules <- list()
    schedules[['abrupt']]  <- data[['abrupt']]$rotation_deg[which(data[['abrupt']]$participant == participants[1])]
    schedules[['gradual']] <- data[['gradual']]$rotation_deg[which(data[['gradual']]$participant == participants[1])]
    
    # random participant sampling: 
    set.seed(sum(participants))
    bs_participants <- matrix(data=sample(x=c(1:length(participants)), size=length(participants)*iterations, replace=TRUE), nrow=iterations, ncol=length(participants))
    
    reachdevs <- list()
    for (condition in c('gradual','abrupt')) {
      
      # empty matrix, to be used for storing & drawing actual data:
      reachdevs[[condition]] <- matrix(data=NA, nrow=max(data[[condition]]$trial), ncol=length(participants) )
      
      # fill matrix
      for (ppno in c(1:length(participants))) {
        reachdevs[[condition]][,ppno] <- data[[condition]]$reachdeviation_deg[which(data[[condition]]$participant == participants[ppno])]
      }
      
    }
    
    # prepared everything for bootstrapping model fits
    
    fits <- list('abrupt.Ls'  = c(),
                 'abrupt.Lf'  = c(),
                 'abrupt.Rs'  = c(),
                 'abrupt.Rf'  = c(),
                 'gradual.Ls' = c(),
                 'gradual.Lf' = c(),
                 'gradual.Rs' = c(),
                 'gradual.Rf' = c()  )
    
    cat(sprintf('bootstrapping %s\n',group))
    for (bs in c(1:iterations)) {
      
      # participant sample for this iteration:
      ps <- as.vector(bs_participants[bs,])
      
      for (condition in c('abrupt','gradual')) {
        
        bs_data <- reachdevs[[condition]][,ps]
        reaches   <- rowMeans(bs_data, na.rm=TRUE)
        
        fit <- twoRateFit(schedule=schedules[[condition]], 
                          reaches=reaches, 
                          gridpoints=9, 
                          gridfits=5, 
                          checkStability=TRUE)
        
        for (parname in names(fit)) {
          fits[[sprintf('%s.%s',condition,parname)]] <- c(fits[[sprintf('%s.%s',condition,parname)]], as.numeric(fit[parname]))
        }
        
      }
      
      cat(sprintf('finished iteration %d / %d\n',bs,iterations))
      
    } 
    
    write.csv(as.data.frame(fits),sprintf('data/%s/bootstrapped_TwoRateFits.csv',group), row.names = FALSE, quote=FALSE)
    
  }
  
}