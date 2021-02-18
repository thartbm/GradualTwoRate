
# here there be functions to get data sets from OSF

library('osfr')

sort30s <- function() {
  
  group <- 'young30'
    
  participants <- c(1:26,28:31)
    
  for (ppno in participants) {
    
    for (condition in c('abrupt','gradual')) {
      
      filename <- sprintf('data/young30/%d/p%02d_%s_selected.txt', ppno, ppno, condition)
      
      df <- read.table(filename)
      
      colnames(df) <- c("phase", "task_name", "trial_type", "trial_num", "terminalfeedback", "rotation_angle", "targetangle_deg", "targetdistance_percmax", "homex_px", "homey_px", "targetx_px", "targety_px", "time_s", "mousex_px", "mousey_px", "cursorx_px", "cursory_px", "trial", "selected", "reachsamples", "NA1", "max_velocity", "NA2")
      
      df <- df[c("phase",  "trial", "rotation_angle", "targetangle_deg", "targetx_px", "targety_px", "time_s", "mousex_px", "mousey_px", "cursorx_px", "cursory_px", "selected", "max_velocity")]
      
      filename <- sprintf('data/young30/young30_p%02d_%s.csv', ppno, condition)
      
      write.csv(df, filename, row.names = F)
      
    }
    
  }
  
}

sort60s <- function() {
  
  group <- 'young60'

  participants <- c(1,2,4:6,8:28,30:32)
  
  for (ppno in participants) {
    
    for (condition in c('abrupt','gradual')) {
      
      filename <- sprintf('data/young60/p%02d_%s_selected.csv', ppno, condition)
      
      df <- read.csv(filename)

      colnames(df) <- c("phase", 
                        "task_name", 
                        "trial_type",
                        "trial_num",
                        "terminalfeedback",
                        "rotation_angle",
                        "targetangle_deg",
                        "targetdistance_percmax",
                        "homex_px",
                        "homey_px",
                        "targetx_px",
                        "targety_px",
                        "time_s",
                        "mousex_px",
                        "mousey_px",
                        "cursorx_px",
                        "cursory_px",
                        "trial",
                        "selected",
                        "max_velocity")
      

      df <- df[c("phase",  "trial", "rotation_angle", "targetangle_deg", "targetx_px", "targety_px", "time_s", "mousex_px", "mousey_px", "cursorx_px", "cursory_px", "selected", "max_velocity")]
      
      filename <- sprintf('data/young60/young60_p%02d_%s.csv', ppno, condition)
      
      write.csv(df, filename, row.names = F)
      
    }
    
  }
  
}

library('R.matlab')

downsampleMatlab2CSV <- function() {
  
  feedback <- rep(1,520)
  feedback[41:80]   <- NA
  feedback[261:300] <- NA
  feedback[481:520] <- NA
  
  digits <- 4
  
  # get rotations from model
  
  if (!dir.exists('data/older45'))  { dir.create('data/older45',recursive=T) }
  if (!dir.exists('data/ataxic45')) { dir.create('data/ataxic45',recursive=T) }
  
  # here is a list of participants:
  participants <- read.csv('/home/marius/Science/StepGradual/Data/participant_list.csv', stringsAsFactors = F)
  
  # we'll down/re-sample all signals to 60 Hz, to be more or less compatible with the other files
  N <- nrow(participants)
  
  for (ppno in c(1:N)) {
    
    output <- NA
    
    participant <- participants$ID[ppno]
    version     <- participants$version[ppno]
    diagnosis   <- participants$diagnosis[ppno]
    
    if (diagnosis == 'control') {
      group <- 'older45'
    } else {
      group <- 'ataxic45'
    }
    
    if (!(version %in% c(2,3,4,5))) {
      next
    }
    print(participant)
    print(version)
    print(group)
    filename <- sprintf('/home/marius/Science/StepGradual/Data/%s/StepGradual_%d__dump.mat',participant,version)
    struct <- readMat(filename)
    
    events <- readMat(sprintf('/home/marius/Science/StepGradual/Data/%s/StepGradual_%d__events.mat',participant,version))
    home <- as.numeric(events[['PrivGlobals']][,,]$CenterPos.cm)
    
    # get distortion from event file (used for fitting)
    #triallist <- unlist(readMat(sprintf('/home/marius/Science/StepGradual/Data/%s/StepGradual_%d__events.mat',participant,version))[['exp.dscr']][[2]][[1]][[2]])
    triallist <- unlist(events[['exp.dscr']][[2]][[1]][[2]])
    trialmat <- matrix(triallist,ncol=4)
    distortion <- c(trialmat[,4])
    distortion[which(trialmat[,2] == 1)] <- NA
    distortion[which(trialmat[,2] == 4)] <- NA
    #print(distortion)
    
    # look at the 'model_data.csv' files to fill in the distortion information: rotations!
    
    data <- as.matrix( struct[[1]][[3]] )

    colnames(data) <- c('trial',
                        'penx_cm',
                        'peny_cm',
                        'penpress',
                        'cursorx_cm',
                        'cursory_cm',
                        'time_s',
                        'state',
                        'targetx_cm',
                        'targety_cm')

    data <- as.data.frame(data)
    data <- data[which(data$trial > 0),]

    for (trial in unique(data$trial)) {
      
      # subset of all data for only this trial
      trialdata <- data[which(data$trial == trial),]

      # data for interpolation:
      dat_time <- trialdata$time_s
      dat_penx <- trialdata$penx_cm - home[1]
      dat_peny <- trialdata$peny_cm - home[2]
      dat_curx <- trialdata$cursorx_cm - home[1]
      dat_cury <- trialdata$cursory_cm - home[2]
      dat_time <- dat_time - dat_time[1]

      # interpolate at 100 Hz:
      time_s <- seq(0,max(dat_time),1/60)

      # mouse x
      smspl <- smooth.spline(x = dat_time,
                             y = dat_penx)
      stylusx_cm <- predict(smspl, x=time_s)$y
      stylusx_cm <- round(stylusx_cm, digits=digits)

      # mouse y
      smspl <- smooth.spline(x = dat_time,
                             y = dat_peny)
      stylusy_cm <- predict(smspl, x=time_s)$y
      stylusy_cm <- round(stylusy_cm, digits=digits)
      
      # cursor x
      smspl <- smooth.spline(x = dat_time,
                             y = dat_curx)
      cursorx_cm <- predict(smspl, x=time_s)$y
      cursorx_cm <- round(cursorx_cm, digits=digits)
      
      # cursor y
      smspl <- smooth.spline(x = dat_time,
                             y = dat_cury)
      cursory_cm <- predict(smspl, x=time_s)$y
      cursory_cm <- round(cursory_cm, digits=digits)
      
      # make new data frame for the trial
      intdat <- data.frame(time_s, stylusx_cm, stylusy_cm, cursorx_cm, cursory_cm)
      intdat$trial <- trial
      intdat$targetx_cm <- round(trialdata$targetx_cm[1] - home[1], digits=digits)
      intdat$targety_cm <- round(trialdata$targety_cm[1] - home[2], digits=digits)
      intdat$feedback   <- feedback[trial]
      intdat$rotation   <- distortion[trial]
      intdat$time_s     <- round(time_s, digits=digits)
      
      intdat$targetangle_deg <- round((atan2(intdat$targety_cm[1], intdat$targetx_cm[1])/pi)*180)


      intdat <- intdat[c('trial','targetangle_deg','feedback','rotation','targetx_cm','targety_cm','time_s','stylusx_cm','stylusy_cm','cursorx_cm','cursory_cm')]

      if (is.data.frame(output)) {
        output <- rbind(output, intdat)
      } else {
        output <- intdat
      }

    }

    write.csv(output, sprintf('data/%s/%s.csv', group, participant), row.names = F)
    
  }
  
}