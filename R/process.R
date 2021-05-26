
getReachDeviations <- function(groups=c('young30', 'young60', 'mcebat45', 'older45', 'young45')) {
  
  for (group in groups) {
    
    for (condition in c('abrupt','gradual')) {
      
      cat(sprintf('working on: %s / %s\n', group, condition))
      
      gcrd <- NA
      
      participants <- getGroupParticipants(group)
      
      for (participant in participants) {
        
        prd <- getParticipantReachDeviations(group, participant, condition)
        
        if (is.data.frame(gcrd)) {
          gcrd <- rbind(gcrd, prd)
        } else {
          gcrd <- prd
        }
        
        
      }
      
      write.csv(gcrd, file=sprintf('data/%s/%s_%s.csv',group,group,condition), row.names = FALSE)
      
    }
    
  }
  
}

getParticipantReachDeviations <- function(group, participant, condition) {
  
  df <- readReachFile(group, participant, condition)
  
  trials <- unique(df$trial)
  
  output <- NA
  
  for (trial in trials) {
    
    trialdf <- df[which(df$trial == trial),]
    trialdata <- as.data.frame(getTrialReachDeviation(trialdf))
    
    if (is.data.frame(output)) {
      output <- rbind(output, trialdata)
    } else {
      output <- trialdata
    }

  }
  
  rd_idx <- which(output$rotation_deg != 0)
  if (length(rd_idx) == 0) {
    
    cat(sprintf('WARNING: no non-zero reach deviations for %s/%d/%s\n',group,participant,condition))
    
  } else {
    
    if (output$rotation_deg[rd_idx[1]] > 0) {
      output$rotation_deg       <- -1 * output$rotation_deg
      output$reachdeviation_deg <- -1 * output$reachdeviation_deg
    }
    
  }
  
  output$participant <- participant
  
  return(output)
  
}

readReachFile <- function(group,participant,condition) {
  
  df <- read.csv(sprintf('data/%s/%s_p%03d_%s.csv',group,group,participant,condition))
  
  if ('rotation_angle' %in% names(df)) {
    names(df)[which(names(df) == 'rotation_angle')] <- 'rotation_deg'
  }
  
  # participants from Munich already have NA for error-clamp trials, but...
  
  if ('phase' %in% names(df)) {
    # data from 30/60 group:
    phases <- unique(df$phase)
    EC_phase <- max(phases)
    df$rotation_deg[which(df$phase == EC_phase)] <- NA
  }
  
  if ('feedback' %in% names(df)) {
    # data from Toronto controls on 45 deg paradigm
    df$rotation_deg[which(is.na(df$feedback))] <- NA
  }
  
  return(df)
  
}

# getTrialReachDeviation <- function(trialdf,at='p25') {
#   
#   
#   
# }

getSplinedTrajectory <- function(x, y, t, length.out=length(t), spar=0.01) {
  
  spl.Xt <- smooth.spline(t, x, spar=spar, keep.data=F) 
  spl.Yt <- smooth.spline(t, y, spar=spar, keep.data=F) 
  
  tt <- seq(min(t), max(t), length.out = length.out)
  
  XX <- predict(spl.Xt, tt)$y
  YY <- predict(spl.Yt, tt)$y
  
  return(data.frame('x'=XX, 'y'=YY, 't'=tt))
  
}


getSplinedVelocity <- function(x, y, t, spar=0.01) {
  
  # spline interpolate the X and Y coordinates over time:
  # (separately... no multi-dimensional splining in base R)
  ST <- getSplinedTrajectory(x, y, t, length.out=length(t), spar=spar)
  
  # velocity on spline interpolated data
  V <- sqrt(diff(ST$x)^2 + diff(ST$y)^2) / diff(ST$t)
  # add velocity = 0 for first sample:
  V <- c(0, V)
  
  return(data.frame('velocity'=V, 'time'=ST$t))
  
}

getSplinedMaxVelocity <- function(x, y, t) {
  
  # we get the splined velocity for the trajectory, moderately smoothed
  VT <- getSplinedVelocity(x, y, t, spar=0.20) # what is good smoothing?
  
  # V is the velocity
  V <- VT$V
  
  # do we need to check if there was some weirdly large velocity?
  print(max(V))
  
  # find the (first) peak in velocity?
  # idx <- which.max(V) # this finds the largest peak
  
  idx <- which(diff(V) < 0)
  
  if (length(idx) > 0) {
    
    smvt <- VT$time[idx[1]]
    smvt.idx <- which.min(abs(t-mvt))
    mvt <- t[smvt.idx]
    
  } else {
    
    smvt <- NA
    smvt.idx <- NA
    mvt <- NA
    
  }
  
  # we return the maximum velocity time stamp, and index, in the actual (not splined) data
  return(c('mvt'=mvt,'mvi'=smvt.idx))
  
}


rotateCoordinates <- function(df,angle,origin=c(0,0)) {
  
  df.names <- names(df)
  
  # create rotation matrix to rotate the X,Y coordinates
  th <- (angle/180) * pi
  R <- t(matrix(data=c(cos(th),sin(th),-sin(th),cos(th)),nrow=2,ncol=2))
  
  # put coordinates in a matrix, and subtract origin
  coordinates <- sweep(as.matrix(df), 2, origin)
  
  # rotate the coordinates, add the origin back in
  df <- as.data.frame(sweep(coordinates %*% R, 2, origin*-1))
  
  # restore column names
  names(df) <- df.names
  
  # return the rotated coordinates
  return(df)
  
}


getTrialReachDeviation <- function(trialdf, location='pr0.25', posunit='cm', timeunit='s', device='stylus', holdvelocity=NA, holdduration=NA) {
  
  # location (string) determines where the angle of the reach is determined, it is one of:
  
  # maxvel: maximum velocity (if data has been manually screened)
  # smv: splined maximum velocity ()
  # endpoint: end of the reach (only makes sense after selection)
  # cmX: the last sample before this distance from home, where X is replaced by a numeral (deprecated: use prX)
  # prX: first sample at or beyond a proportion of distance from home to target, given by X (e.g. 'pr0.333333')
  # hold: at a hold point; also specify how long and at what maximum velocity the hold has to be in other arguments
  # smvX: first velocity peak in spline-smoothed trajectory, beyond a percentage distance from home to target given by X (e.g. 'smv0.10')
  
  
  # return a matrix of two numbers:
  reachangle = matrix(data=NA,nrow=1,ncol=6)
  colnames(reachangle) <- c( 'reachdeviation_deg', 
                             'targetangle_deg',
                             'rotation_deg',
                             sprintf('%sx_%s',device,posunit), 
                             sprintf('%sy_%s',device,posunit), 
                             sprintf('time_%s',timeunit) )
  
  # extract the relevant reach information
  x <- trialdf[,sprintf('%sx_%s',device,posunit)]
  y <- trialdf[,sprintf('%sy_%s',device,posunit)]
  t <- trialdf[,sprintf('time_%s',timeunit)]
  t <- t - min(t)
  
  angle <- trialdf[1,'targetangle_deg']
  target <- as.numeric(trialdf[1,c(sprintf('targetx_%s',posunit),sprintf('targety_%s',posunit))])
  
  # always return the target angle?
  reachangle[1,2] <- angle
  rot <- trialdf$rotation_deg[1]
  reachangle[1,3] <- rot
  
  # rotate the trajectory:
  # - avoids problems with atan2 angles around 180 / -180
  # - puts the target at 0, so angular deviation is easy to get
  trajectory <- rotateCoordinates(data.frame(x,y),-1*angle)
  x <- trajectory[,1]
  y <- trajectory[,2]
  
  # now try find the specified location in this reach:
  # if we can't find it, we need to know
  invalidlocation <- TRUE
  
  # maximum velocity, should be a column in the data...
  # this only happens with preprocessing or manual screening
  # use 'smv' if this is not the case...
  if (location == 'maxvel') {
    MV <- trialdf[,'maxvel']
    rown <- which(MV == 1)
    if (length(rown) > 1) {
      rown <- rown[1]
    }
    if (length(rown) == 0) {
      # no maximum velocity defined!
      return(reachangle)
    }
    invalidlocation <- FALSE
  }
  
  if (location == 'smv') {
    
  }
  
  # end point, just the last point in the reach
  if (location == 'endpoint') {
    rown <- length(x)
    invalidlocation <- FALSE
  }
  
  # cutoff in centimers, the last sample before this cutoff distance is reached
  # this assumes that people don't go back, or that there is only one movement from home to target
  if (substring(location,1,2) == 'cm') {
    distance <- as.numeric(substring(location, 3))
    
    # get the distance from home:
    dist <- sqrt(x^2 + y^2)
    
    # if there are no selected samples below 3 cm: return NAs
    if (length(which(dist < distance)) == 0) {
      return(reachangle)
    }
    
    # find the last sample, where dist < 3
    rown <- max(which(dist < distance))
    invalidlocation <- FALSE
  }
  
  # cutoff at a percentage from home to target in whatever unit is used
  if (substring(location,1,2) == 'pr') {
    distance <- as.numeric(substring(location, 3))
    #distance <- distance * sqrt(trialdf$targetx_pix[1]^2 + trialdf$targety_pix[1]^2)
    distance <- distance * sqrt(sum(target^2))
    
    # get the distance from home:
    dist <- sqrt(x^2 + y^2)
    
    # if there are no selected samples above 3 cm: return NAs
    if (length(which(dist > distance)) == 0) {
      return(reachangle)
    }
    
    # find the first sample, where dist > X
    rown <- min(which(dist > distance))
    invalidlocation <- FALSE
  }
  
  # find the first 'hold':
  if (substring(location,1,4) == 'hold') {
    holddistance <- as.numeric(substring(location, 5))
    #holdvelocity
    #holdduration # in timeunit: 's' or 'ms'?
    
    
    
  }
  
  # use smooth-splined trajectory to get angle at *first* velocity peak:
  if (substring(location,1,3) == 'smv') {
    
    # how far does the vleocity peak have to be away from the home position
    # (as percentage of home-target distance)
    if (nchar(location) > 3) {
      distance <- as.numeric(substring(location, 4))
    } else {
      distance <- 0.05
    }
    distance <- distance * sqrt(sum(target^2))
    
    dist <- sqrt(x^2 + y^2)
    
    VT <- getSplinedVelocity(x, y, t, spar=0.20)
    v <- c(0, 0, VT$velocity)
    
    peaks <- which(diff(sign(diff(v))) == -2 & dist > distance)
    if (length(peaks) > 0) {
      rown <- peaks[1]
      invalidlocation <- FALSE
    }
    
  }
  
  
  
  # if we don't have a valid location, we can't calculate an angle
  if (invalidlocation) {
    return(reachangle)
  }
  
  # calculate the angle at that point for the rotated trajectory
  # this is the angular deviation we are looking for
  angulardeviation <- (atan2(y[rown],x[rown]) / pi) * 180
  
  # put the result in the little matrix:
  reachangle[1,1] <- angulardeviation
  reachangle[1,2] <- angle
  reachangle[1,3] <- rot
  reachangle[1,4] <- trialdf[rown,sprintf('%sx_%s',device,posunit)]
  reachangle[1,5] <- trialdf[rown,sprintf('%sy_%s',device,posunit)]
  reachangle[1,6] <- t[rown]
  
  return(reachangle)
  
}
