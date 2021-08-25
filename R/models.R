

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


twoRateFit <- function(schedule, reaches, gridpoints=14, gridfits=10, checkStability=TRUE) {
  
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


getGroupFits <- function(exp=NULL, groups=NULL) {
  
  info <- getInfo()
  
  if (!is.null(exp) & is.null(groups)) {
    groups <- info$group[which(info$experiment == exp)]
  }
  
  alldata <- getSelectedGroupsData(groups=groups)
  alldata <- parseData(alldata, FUN=baseline)
  alldata <- parseData(alldata, FUN=addtrial)
  
  for (group in groups) {
    
    groupdata <- alldata[[group]]
    participants <- unique(groupdata[['abrupt']]$participant)
    
    fits <- list('participant' = c(),
                 'abrupt.Ls'   = c(),
                 'abrupt.Lf'   = c(),
                 'abrupt.Rs'   = c(),
                 'abrupt.Rf'   = c(),
                 'gradual.Ls'  = c(),
                 'gradual.Lf'  = c(),
                 'gradual.Rs'  = c(),
                 'gradual.Rf'  = c()  )
    
    for (participant in participants) {
      
      for (condition in c('abrupt','gradual')) {
        
        pdf <- groupdata[[condition]][which(groupdata[[condition]]$participant == participant),]
        
        schedule <- as.numeric(pdf$rotation_deg)
        reaches  <- as.numeric(pdf$reachdeviation_deg)
        
        fit <- twoRateFit(schedule=schedule,
                          reaches=reaches)
        
        for (parameter in c('Ls','Lf','Rs','Rf')) {
          fits[[sprintf('%s.%s',condition,parameter)]] <- c(fits[[sprintf('%s.%s',condition,parameter)]], fit[parameter])
        }
        
      }
      
      fits[['participant']] <- c(fits[['participant']], participant)
      print(participant)
    }
    
    for (condition in c('abrupt','gradual')) {
      
      schedule <- groupdata[[condition]]$rotation_deg[which(groupdata[[condition]]$participant == participants[1])]
      
      reaches <- as.numeric( aggregate(reachdeviation_deg ~ trial, data=groupdata[[condition]], FUN=mean)$reachdeviation_deg )
      
      fit <- twoRateFit(schedule=schedule,
                        reaches=reaches)
      
      for (parameter in c('Ls','Lf','Rs','Rf')) {
        fits[[sprintf('%s.%s',condition,parameter)]] <- c(fits[[sprintf('%s.%s',condition,parameter)]], fit[parameter])
      }
      
    }
    
    fits[['participant']] <- c(fits[['participant']], -1)
    
    write.csv(as.data.frame(fits), 
              sprintf('data/%s/twoRateFits.csv',group), 
              row.names = FALSE,
              quote = FALSE)
    
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
      
      cols <- length(participants)
      rows <- max(data[[condition]]$trial, na.rm=TRUE)
      
      # empty matrix, to be used for storing & drawing actual data:
      reachdevs[[condition]] <- matrix(data=NA, nrow=rows, ncol=cols )
      
      # fill matrix
      for (ppno in c(1:length(participants))) {
        pprd <- data[[condition]]$reachdeviation_deg[which(data[[condition]]$participant == participants[ppno])]
        
        if (participants[ppno] == 208 & condition == 'gradual') {pprd <- c(pprd,NA)} # this participant is missing the last clamped trial

        reachdevs[[condition]][,ppno] <- pprd
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
        
        #print(schedules[[condition]])
        
        fit <- twoRateFit(schedule=schedules[[condition]], 
                          reaches=reaches, 
                          gridpoints=14, 
                          gridfits=10, 
                          checkStability=TRUE)
        
        #print(c(fit,'MSE'=twoRateMSE(par=fit,schedule=schedules[[condition]],reaches=reaches,checkStability = FALSE)))
        
        #plot(reaches,col='gray',main=sprintf('MSE: %0.3f',twoRateMSE(par=fit,schedule = schedules[[condition]], reaches=reaches,checkStability = T)))
        #lines(twoRateModel(par=fit,schedule=schedules[[condition]])$total,col='red')
        
        for (parname in names(fit)) {
          fits[[sprintf('%s.%s',condition,parname)]] <- c(fits[[sprintf('%s.%s',condition,parname)]], as.numeric(fit[parname]))
        }
        
      }
      
      cat(sprintf('finished iteration %d / %d\n',bs,iterations))
      
    } 
    
    write.csv(as.data.frame(fits),sprintf('data/%s/bootstrapped_TwoRateFits.csv',group), row.names = FALSE, quote=FALSE)
    
  }
  
}


modelMSEs <- function(exp=NULL, groups=NULL) {

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
    
    
    all_fits <- read.csv(sprintf('data/%s/bootstrapped_TwoRateFits.csv', group), stringsAsFactors = F )
    
    iterations <- dim(all_fits)[1]
    
    # random participant sampling: 
    set.seed(sum(participants))
    bs_participants <- matrix(data=sample(x=c(1:length(participants)), size=length(participants)*iterations, replace=TRUE), nrow=iterations, ncol=length(participants))
    
    
    reachdevs <- list()
    for (condition in c('gradual','abrupt')) {
      
      cols <- length(participants)
      rows <- max(data[[condition]]$trial, na.rm=TRUE)
      
      # empty matrix, to be used for storing & drawing actual data:
      reachdevs[[condition]] <- matrix(data=NA, nrow=rows, ncol=cols )
      
      # fill matrix
      for (ppno in c(1:length(participants))) {
        pprd <- data[[condition]]$reachdeviation_deg[which(data[[condition]]$participant == participants[ppno])]
        
        if (participants[ppno] == 208 & condition == 'gradual') {pprd <- c(pprd,NA)} # this participant is missing the last clamped trial
        
        reachdevs[[condition]][,ppno] <- pprd
      }
      
    } # 
    
    #print(reachdevs)
    
    abfits2abdata.MSE <- c()
    abfits2grdata.MSE <- c()
    grfits2grdata.MSE <- c()
    grfits2abdata.MSE <- c()
    
    for (fit_no in c(1:iterations)) {
      
      # # # # # # # # # # # # 
      #
      #   HERE GET ALL MSEs
      #
      # # # # # # # # # # # #
      
      # participant sample for this iteration:
      ps <- as.vector(bs_participants[fit_no,])
      
      bs_gradual_data <- reachdevs[['gradual']][,ps]
      bs_abrupt_data <- reachdevs[['abrupt']][,ps]
      gradual_reaches   <- rowMeans(bs_gradual_data, na.rm=TRUE)
      abrupt_reaches   <- rowMeans(bs_abrupt_data, na.rm=TRUE)
      
      two_fits <- all_fits[fit_no,]
      fits <- list()
      fits[['abrupt']] <- c()
      fits[['gradual']] <- c()
      for (condition in c('gradual','abrupt')) {
        for (parameter in c('Ls','Lf','Rs','Rf')) {
          fits[[condition]] <- c(fits[[condition]], two_fits[,sprintf('%s.%s',condition,parameter)])
        }
        names(fits[[condition]]) <- c('Ls','Lf','Rs','Rf')
      }
      
      abfits2abdata.MSE <- c(abfits2abdata.MSE,  twoRateMSE(par=fits[['abrupt']],  schedule=schedules[['abrupt']], reaches=abrupt_reaches))
      abfits2grdata.MSE <- c(abfits2grdata.MSE,  twoRateMSE(par=fits[['abrupt']],  schedule=schedules[['gradual']], reaches=gradual_reaches))
      grfits2grdata.MSE <- c(grfits2grdata.MSE, twoRateMSE(par=fits[['gradual']], schedule=schedules[['gradual']], reaches=gradual_reaches))
      grfits2abdata.MSE <- c(grfits2abdata.MSE, twoRateMSE(par=fits[['gradual']], schedule=schedules[['abrupt']], reaches=abrupt_reaches))
      
    }
    
    write.csv(data.frame(abfits2abdata.MSE, grfits2grdata.MSE, abfits2grdata.MSE, grfits2abdata.MSE), 
              sprintf('data/%s/modelMSEs.csv',group), row.names=F, quote=F)
    
    
  } # end loop through groups
    
}

# exponential decay model -----



asymptoticDecayModel <- function(par, schedule) {
  
  # the process and error states are initialized at 0:
  Pt <- 0
  Et <- 0
  
  # the total output is stored here:
  output <- c()
  
  for (t in c(1:length(schedule))) {
    
    Pt <- Pt - (par['lambda'] * Et)
    
    # now we calculate what the previous error will be for the next trial:
    if (is.na(schedule[t])) {
      Et <- 0
    } else {
      Et <- Pt + (schedule[t] * par['N0'])
    }
    #print(Et)
    # at this point we save the process state in our vector:
    output <- c(output, Pt)
    
  }
  
  # done all the trials: return output
  return(data.frame(output))
  
}

asymptoticDecayMSE <- function(par, schedule, signal) {
  
  MSE <- mean((asymptoticDecayModel(par, schedule)$output - signal)^2, na.rm=TRUE)
  
  return( MSE )
  
}


asymptoticDecayFit <- function(schedule, signal, gridpoints=10, gridfits=5) {
  
  # set the search grid:
  parvals <- seq(1/gridpoints/2,1-(1/gridpoints/2),1/gridpoints)
  
  maxAsymptote <- 2*max(abs(signal), na.rm=TRUE)
  
  # define the search grid:
  searchgrid <- expand.grid('lambda' = parvals, 
                            'N0'     = parvals * maxAsymptote)
  
  # evaluate starting positions:
  MSE <- apply(searchgrid, FUN=asymptoticDecayMSE, MARGIN=c(1), schedule=schedule, signal=signal)
  
  # testing if optimx is installed and making it available it so:
  optimxInstalled <- require("optimx")
  
  if (optimxInstalled) {
    
    # run optimx on the best starting positions:
    allfits <- do.call("rbind",
                       apply( data.frame(searchgrid[order(MSE)[1:gridfits],]),
                              MARGIN=c(1),
                              FUN=optimx::optimx,
                              fn=asymptoticDecayMSE,
                              method='L-BFGS-B',
                              lower=c(0,0),
                              upper=c(1,maxAsymptote),
                              schedule=schedule,
                              signal=signal ) )
    
    # pick the best fit:
    win <- allfits[order(allfits$value)[1],]
    
    # return the best parameters:
    return(unlist(win[1:2]))
    
  } else {
    
    cat('(consider installing optimx, falling back on optim now)\n')
    
    # use optim with Nelder-Mead after all:
    allfits <- do.call("rbind",
                       apply( data.frame(searchgrid[order(MSE)[1:gridfits],]),
                              MARGIN=c(1),
                              FUN=optim,
                              fn=asymptoticDecayMSE,
                              method='Nelder-Mead',
                              schedule=schedule,
                              signal=signal ) )
    
    # pick the best fit:
    win <- allfits[order(unlist(data.frame(allfits)[,'value']))[1],]
    
    # return the best parameters:
    return(win$par)
    
  }
  
}


getLearningRates <- function(exp=NULL, groups=NULL) {
  
  info <- getInfo()
  
  if (is.null(groups)) {
    groups <- info$group[which(info$experiment == exp)]
  }
  
  data <- getSelectedGroupsData(groups=groups)
  data <- parseData(data, FUN=baseline)
  data <- parseData(data, FUN=normalize)
  data <- parseData(data, FUN=addtrial)
  data <- parseData(data, FUN=addphase)

  for (group in groups) {
    group_data <- data[[group]]
    abrupt <- group_data[['abrupt']]
    
    gdf <- getVersion(group=group)
    gdf <- subset(gdf, select = -c(version)) 
    gdf$lambda <- NA
    gdf$N0 <- NA
    
    participants <- unique(abrupt$participant)
    
    gdf <- gdf[which(gdf$participant %in% participants),]
    
    for (participant in participants) {
      apd <- as.numeric( abrupt[which(abrupt$participant == participant & abrupt$phase == 2),'reachdeviation_deg'] )
      fit <- asymptoticDecayFit(schedule=rep(-1,length(apd)),
                                signal=apd,
                                gridpoints=15,
                                gridfits=5)
      gdf$lambda[which(gdf$participant == participant)] <- fit['lambda']
      gdf$N0[which(gdf$participant == participant)] <- fit['N0']
    }
    
    write.csv(gdf, 
              file=sprintf('data/%s/abrupt_learningrates.csv',group),
              row.names = FALSE,
              quote = FALSE)
    
  }
  
}


# parameter recovery -----

parameterRecoverySimulation <- function(iterations=1000) {
  
  # we will test if parameters are more recoverable from schedules with actual counter-rotation phases
  # to do this, we simulate data, based on 4 schedules:
  
  schedules <- list()
  schedules[['zero']] <- list()
  schedules[['cntr']] <- list()
  schedules[['zero']][['abrupt']] <- c(rep(0,40),rep(-45,120),                                rep(0,20), rep(NA,40)) 
  schedules[['zero']][['ramped']] <- c(rep(0,40),seq(0,-45,length.out = 61)[1:60],rep(-45,60),rep(0,20), rep(NA,40)) 
  schedules[['cntr']][['abrupt']] <- c(rep(0,40),rep(-45,120),                                rep(45,20),rep(NA,40))
  schedules[['cntr']][['ramped']] <- c(rep(0,40),seq(0,-45,length.out = 61)[1:60],rep(-45,60),rep(45,20),rep(NA,40))
  
  # we will simulate noise, based on data from the rotation magnitude experiment
  # and add the same noise to model-generated data based on all four perturbation schedules
  # and then fit the model to that data
  
  # here we get that noise:
  info   <- getInfo()
  groups <- as.character(info$group[which(info$experiment == 1)])

  alldata <- getSelectedGroupsData(groups=groups)
  alldata <- parseData(alldata, FUN=baseline)
  alldata <- parseData(alldata, FUN=addtrial)
  
  MSEs <- c()
  
  for (group in groups) {
    group_fits <- read.csv(sprintf('data/%s/twoRateFits.csv',group), stringsAsFactors = FALSE)
    group_idx <- which(group_fits$participant == -1)
    
    for (condition in c('abrupt','gradual')) {
      
      # we read the parameters for the group for this condition:
      fit <- c('Ls'=NA,'Lf'=NA,'Rs'=NA,'Rf'=NA)
      for (parameter in names(fit)) {
        column <- sprintf('%s.%s',condition,parameter)
        fit[parameter] <- as.numeric(group_fits[group_idx,column])
      }
      
      # we get the average reaches and the perturbation schedule:
      reaches <- as.numeric(unlist(aggregate(reachdeviation_deg ~ trial, data=alldata[[group]][[condition]], FUN=mean, na.rm=TRUE)$reachdeviation_deg))
      schedule <- alldata[[group]][[condition]]$rotation_deg[which(alldata[[group]][[condition]]$participant == alldata[[group]][[condition]]$participant[1])]
      
      MSE <- twoRateMSE(par=fit,
                        schedule=schedule, 
                        reaches=reaches, 
                        checkStability=FALSE)
      
      MSEs <- c(MSEs, MSE)
      
    }
  }
  
  # the square root of the average MSE should be a reasonable estimate of the 
  # level of variation we can expect in a group (standard deviation)
  std <- sqrt(mean(MSEs))
  
  # we want ground truth parameters in this vector:
  ground_truth_parameters <- c('Ls'=NA,
                               'Lf'=NA,
                               'Rs'=NA,
                               'Rf'=NA)
  
  # we'll use the young45 group's abrupt fit for the group:
  y45fits <- read.csv('data/young45/twoRateFits.csv', 
                      stringsAsFactors = FALSE)
  
  # and put those parameters in our vector:
  group_idx <- which(y45fits$participant == -1)
  for (par in names(ground_truth_parameters)) {
    column <- sprintf('abrupt.%s',par)
    ground_truth_parameters[par] <- as.numeric(unlist(y45fits[group_idx,column]))
  }
  
  # now we want the basic reaches without noise, which will always be the same:
  reaches <- list('zero'=list('abrupt'=c(),'ramped'=c()),'cntr'=list('abrupt'=c(),'ramped'=c()))
  for (counter in c('zero','cntr')) {
    for(condition in c('abrupt','ramped')) {
      reaches[[counter]][[condition]] <- twoRateModel( par=ground_truth_parameters,
                                                       schedule = schedules[[counter]][[condition]])$total
    }
  }
  
  # we set the seed to some fixed, but unknown value:
  set.seed(sum(c(unique(alldata[[groups[1]]][['abrupt']]$participant), 
                 unique(alldata[[groups[2]]][['abrupt']]$participant) )))
  
  # here we'll store the results of the simulation:
  empty <- list('Ls'=rep(NA,iterations),
                'Lf'=rep(NA,iterations),
                'Rs'=rep(NA,iterations),
                'Rf'=rep(NA,iterations))
  recovered_fits <- list('zero'=list('abrupt'=empty,
                                     'ramped'=empty),
                         'cntr'=list('abrupt'=empty,
                                     'ramped'=empty))
  
  # then we loop through the iterations:
  for (bs in c(1:iterations)) {
    
    # first, let's make some noise:
    noise <- rnorm(220, mean=0, sd=std)
    
    for (counter in c('zero','cntr')) {
      for(condition in c('abrupt','ramped')) {
        
        fit <- twoRateFit(schedule = schedules[[counter]][[condition]],
                          reaches  =   reaches[[counter]][[condition]]+noise)
        
        for (par in names(fit)) {
          parval <- fit[par]
          names(parval) <- c()
          recovered_fits[[counter]][[condition]][[par]][bs] <- parval
        }

      }
    }
    cat(sprintf('finished iteration: %d / %d\n',bs,iterations))
  }


  # by repeating this 1000 times for a given set of parameter values
  # we can get confidence intervals for how far off the fitted parameter values are
  # in each of the four perturbation schedules
  
  for (counter in c('zero','cntr')) {
    for(condition in c('abrupt','ramped')) {
      df <- as.data.frame(recovered_fits[[counter]][[condition]])
      #print(recovered_fits[[counter]][[condition]])
      #print(data.frame(recovered_fits[[counter]][[condition]]))
      write.csv(df, 
                file=sprintf('data/young45/recovered_fits_%s_%s.csv',counter,condition), 
                row.names = FALSE, 
                quote = FALSE)
    }
  }
  
  # and we're done with the simulation... now it still needs to be interpreted
  
}