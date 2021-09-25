
# The functions in this file won't run on your computer.
# They were used to get the data from all experiments in the same,
# more readable format. 
# They are provided here as a form of documentation.

library('osfr')

sort30s <- function() {
  
  group <- 'young30'
  
  participants <- c(1:26,28:31)
    
  for (ppno in participants) {
    
    for (condition in c('abrupt','gradual')) {
      #cat(sprintf('%d - %s\n',ppno, condition))
      filename <- sprintf('data/ORG_data/young30org/%d/p%02d_%s_selected.txt', ppno, ppno, condition)
      # print(filename)
      df <- read.table(filename)
      
      colnames(df) <- c("phase", "task_name", "trial_type", "trial_num", "terminalfeedback", "rotation_angle", "targetangle_deg", "targetdistance_percmax", "homex_px", "homey_px", "targetx_px", "targety_px", "time_s", "mousex_px", "mousey_px", "cursorx_px", "cursory_px", "trial", "selected", "reachsamples", "NA1", "max_velocity", "NA2")
      #print(which(names(df) == 'time_s'))
      
      # monitor resolution: 1680 x 1050
      # screen size: 433.5 x 271 mm
      # tablet size: 311   x 216 mm
      #             0.717    0.797
      
      # --> SOUNDS LIKE 72% (PyVMEC screen resize setting)
      
      ppc <- 38.754789272
      ppc <- ppc / 0.72
      ppc <- ppc / 18.3333333333
      #ppc <- 1
      df$targetx_cm <- round(df$targetx_px / ppc, digits=3)
      df$targety_cm <- round(df$targety_px / ppc, digits=3)
      df$mousex_cm  <- round(df$mousex_px / ppc, digits=3)
      df$mousey_cm  <- round(df$mousey_px / ppc, digits=3)
      df$cursorx_cm <- round(df$cursorx_px / ppc, digits=3)
      df$cursory_cm <- round(df$cursory_px / ppc, digits=3)
      
      # print(df$targetx_px[1])
      # print(df$targetx_cm[1])
      
      df <- df[c("phase",  "trial", "rotation_angle", "targetangle_deg", "targetx_cm", "targety_cm", "time_s", "mousex_cm", "mousey_cm", "cursorx_cm", "cursory_cm", "selected", "max_velocity")]
      
      df$time_s <- round(df$time_s / 1000, digits=3)
      
      filename <- sprintf('data/young30/young30_p%03d_%s.csv', ppno, condition)
      write.csv(standardNames(df), filename, row.names = F)
      
    }
    
  }
  
}

sort60s <- function() {
  
  group <- 'young60'

  participants <- c(1,2,4:6,8:28,30:32)
  
  for (ppno in participants) {
    
    for (condition in c('abrupt','gradual')) {
      
      filename <- sprintf('data/ORG_data/young60org/p%02d_%s_selected.csv', ppno, condition)
      
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
      
      # monitor resolution: 1680 x 1050
      ppc <- 38.754789272
      ppc <- ppc / 0.72
      df$targetx_cm <- round(df$targetx_px / ppc, digits=3)
      df$targety_cm <- round(df$targety_px / ppc, digits=3)
      df$mousex_cm  <- round(df$mousex_px / ppc, digits=3)
      df$mousey_cm  <- round(df$mousey_px / ppc, digits=3)
      df$cursorx_cm <- round(df$cursorx_px / ppc, digits=3)
      df$cursory_cm <- round(df$cursory_px / ppc, digits=3)
      
      df <- df[c("phase",  "trial", "rotation_angle", "targetangle_deg", "targetx_cm", "targety_cm", "time_s", "mousex_cm", "mousey_cm", "cursorx_cm", "cursory_cm", "selected", "max_velocity")]
      
      #df <- df[c("phase",  "trial", "rotation_angle", "targetangle_deg", "targetx_px", "targety_px", "time_s", "mousex_px", "mousey_px", "cursorx_px", "cursory_px", "selected", "max_velocity")]
      
      df$time_s <- round(df$time_s, digits=3)
      
      filename <- sprintf('data/young60/young60_p%03d_%s.csv', ppno+100, condition)
      
      write.csv(standardNames(df), filename, row.names = F)
      
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
  if (!dir.exists('data/mcebat45')) { dir.create('data/mcebat45',recursive=T) }
  
  # here is a list of participants:
  participants <- read.csv('/home/marius/Science/StepGradual/Data/participant_list_2.csv', stringsAsFactors = F)
  
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
      group <- 'mcebat45'
    }
    
    if (!(version %in% c(2,3,4,5))) {
      next
    }
    #if (participant != 'ShSh') { next }
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

      # interpolate at 60 Hz:
      time_s <- seq(0,max(dat_time),1/60)

      # mouse x
      # smspl <- smooth.spline(x = dat_time,
      #                        y = dat_penx)
      # stylusx_cm <- predict(smspl, x=time_s)$y
      stylusx_cm <- approx(x = dat_time,
                           y = dat_penx,
                           xout = time_s)$y
      stylusx_cm <- round(stylusx_cm, digits=digits)

      # mouse y
      # smspl <- smooth.spline(x = dat_time,
      #                        y = dat_peny)
      # stylusy_cm <- predict(smspl, x=time_s)$y
      stylusy_cm <- approx(x = dat_time,
                           y = dat_peny,
                           xout = time_s)$y
      stylusy_cm <- round(stylusy_cm, digits=digits)
      
      # cursor x
      # smspl <- smooth.spline(x = dat_time,
      #                        y = dat_curx)
      # cursorx_cm <- predict(smspl, x=time_s)$y
      cursorx_cm <- approx(x = dat_time,
                           y = dat_curx,
                           xout = time_s)$y
      cursorx_cm <- round(cursorx_cm, digits=digits)
      
      # cursor y
      # smspl <- smooth.spline(x = dat_time,
      #                        y = dat_cury)
      # cursory_cm <- predict(smspl, x=time_s)$y
      cursory_cm <- approx(x = dat_time,
                           y = dat_cury,
                           xout = time_s)$y
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

getTorontoControls <- function() {
  
  sourcedir <- '/home/marius/Science/StepGradual/ctrl/data/combined'
  
  columns <- c('trialnumber',
               'trialtype',
               'distortion_deg',
               'targetdirection_deg',
               'time_s',
               'step',
               'targetx_cm',
               'targety_cm',
               'handx_cm',
               'handy_cm',
               'cursorx_cm',
               'cursory_cm')
  
  for (group in c('older45','young45')) {
  #for (group in c('young45')) {
    
    targetdir <- sprintf('data/%s/',group)
    
    if (!dir.exists(targetdir)) { dir.create(targetdir, recursive=TRUE) }
    
    participants <- list('young45'=c(0:11,13:27),
                         'older45'=c(200:213))[[group]]
    
    print(participants)
    
    for (ppno in participants) {
      
      # version is decided by participant number:
      version <- (ppno %% 4) + 2
      
      filename <- sprintf('%s/P%03dV%d.txt', sourcedir, ppno, version)
      
      df <- read.table(filename, header=TRUE)
      names(df)[which(names(df) == 'time_ms')] <- 'time_s'
      # remove the X11 mouse after all
      df <- df[,columns]
      
      # *****************************
      # time_s should be converted to time_ms !!!!
      # or at least it should actually be seconds...
      # *****************************
      
      #names(df)[which(names(df) == 'time_s')] <- 'time_s'
      #print(df$time_s)
      df$time_s <- round(df$time_s, digits=3)
      #print(df$time_s)
      
      df$targetx_cm <- round(df$targetx_cm, digits=3)
      df$targety_cm <- round(df$targety_cm, digits=3)
      df$handx_cm   <- round(df$handx_cm,   digits=3)
      df$handy_cm   <- round(df$handy_cm,   digits=3)
      df$cursorx_cm <- round(df$cursorx_cm, digits=3)
      df$cursory_cm <- round(df$cursory_cm, digits=3)
      
      #cursortrials <- (unique(df$trialnumber[which(df$trialtype == 'Feedback')]))
      errorclamptrials <- (unique(df$trialnumber[which(df$trialtype == 'Error_Clamp')]))
      df$trialtype <- 1
      df$trialtype <- as.numeric(df$trialtype)
      df$trialtype[which(df$trialnumber %in% errorclamptrials)] <- NA
      
      for (start in c(81,301)) {
        
        trials <- c(start:(start+219))
        partdf <- df[which(df$trialnumber %in% trials),]
        
        if (abs(partdf$distortion_deg[which(partdf$trialnumber == (start+50))[1]]) == 45) {
          condition <- 'abrupt'
        } else {
          condition <- 'gradual'
        }
        
        useppno <- ppno
        
        if (group == 'young45') {
          useppno <- useppno + 500
        }
        
        write.csv(x = standardNames(partdf),
                 file = sprintf('data/%s/%s_p%03d_%s.csv', group, group, useppno, condition),
                 row.names = FALSE)
        #cat(sprintf('data/%s/%s_p%03d_%s.csv\n', group, group, useppno, condition))
        #print(str(partdf))
        #print(df$time_s[1:10])
        #print(names(df))
      }
      
    }
    
  }
  
}

splitMunichData <- function() {
  
  for (group in c('older45', 'mcebat45')) {
    
    print(group)
    
    participants <- read.csv(sprintf('data/%s_demographics_org.csv', group), stringsAsFactors = FALSE)
    
    participants$participant <- NA
    
    ppno_offset <- list('mcebat45'=400, 'older45'=300)[[group]]
    
    for (ppno in c(1:length(participants$ID))) {
      
      participant <- participants$ID[ppno]
      participants$participant[ppno] <- ppno + ppno_offset

      # print(participants$participant[ppno])
      # print(ppno)
      # print(group)
      print(sprintf('data/%s/%s.csv',group,participant))
      df <- read.csv(sprintf('data/%s/%s.csv',group,participant))
      
      part1 <- df[which(df$trial %in% c(81:300)),]
      part1$trial <- part1$trial - 80
      part2 <- df[which(df$trial %in% c(301:520)),]
      part2$trial <- part2$trial - 300
      
      part1 <- standardNames(part1)
      part2 <- standardNames(part2)
      
      if ( abs(df$rotation[which(df$trial == 131)[1]]) == 45 ) {
        # part 1 is abrupt
        write.csv(part1, sprintf('data/%s/%s_p%03d_abrupt.csv', group, group, ppno+ppno_offset), row.names = FALSE)
        write.csv(part2, sprintf('data/%s/%s_p%03d_gradual.csv', group, group, ppno+ppno_offset), row.names = FALSE)
      } else {
        # part 1 is gradual
        write.csv(part1, sprintf('data/%s/%s_p%03d_gradual.csv', group, group, ppno+ppno_offset), row.names = FALSE)
        write.csv(part2, sprintf('data/%s/%s_p%03d_abrupt.csv', group, group, ppno+ppno_offset), row.names = FALSE)
      }
      
      
    }
    
    write.csv(participants, sprintf('data/%s_demographics.csv', group), row.names = FALSE)
    
  }
  
}


standardNames <- function(df) {
  
  # change 'distortion_deg' and 'rotation' into 'rotation_deg'
  if ('distortion_deg' %in% names(df)) {
    names(df)[which(names(df) == 'distortion_deg')] <- 'rotation_deg'
  }
  if ('rotation' %in% names(df)) {
    names(df)[which(names(df) == 'rotation')] <- 'rotation_deg'
  }
  
  if ('targetdirection_deg' %in% names(df)) {
    names(df)[which(names(df) == 'targetdirection_deg')] <- 'targetangle_deg'
  }
  
  # rename 'trialnumber' to 'trial'
  if ('trialnumber' %in% names(df)) {
    names(df)[which(names(df) == 'trialnumber')] <- 'trial'
  }
  
  
  # rename hand and mouse to stylus
  if ('handx_cm' %in% names(df)) {
    names(df)[which(names(df) == 'handx_cm')] <- 'stylusx_cm'
  }
  if ('handy_cm' %in% names(df)) {
    names(df)[which(names(df) == 'handy_cm')] <- 'stylusy_cm'
  }
  if ('mousex_cm' %in% names(df)) {
    names(df)[which(names(df) == 'mousex_cm')] <- 'stylusx_cm'
  }
  if ('mousey_cm' %in% names(df)) {
    names(df)[which(names(df) == 'mousey_cm')] <- 'stylusy_cm'
  }
  
  
  
  if ('trialtype' %in% names(df)) {
    names(df)[which(names(df) == 'trialtype')] <- 'feedback'
  }
  
  return(df)
  
}


downloadOSFdata <- function(groups=c('young30', 'young60', 'mcebat45', 'older45', 'young45'), overwrite=FALSE, removezip=FALSE) {
  
  if (overwrite) {
    conflicts = 'overwrite'
  } else {
    conflicts = 'skip'
  }
  
  osfr::osf_auth(Sys.getenv("OSF_PAT"))
  
  OSFnode <- osfr::osf_retrieve_node("c5ezv")
  
  # get a list of files for the year and semester that is requested:
  files <- osfr::osf_ls_files(OSFnode)
  
  print(files)
  
  for (group in groups) {
    
    # get list of files:
    filenames <- c(sprintf('%s.zip',group))
    
    for (filename in filenames) {
      
      print(filename)
      
      # find which line corresponds to the file:
      idx <- which(files$name == filename)
      
      # check that the file exists on OSF, and is unique:
      # if not: abort
      if (length(idx) > 1) {
        # no unique file found: aborting this one
        if (length(idx) == 0) {
          return(FALSE)
        }
        return(NULL)
      }
      
      # download the file:
      if (!file.exists(sprintf('data/%s',files$name[idx])) | overwrite) {
        osfr::osf_download(x = files[idx,], 
                           path = 'data/', 
                           conflicts=conflicts)
      }
      
      # check if it is a zip file:
      if (grepl('\\.zip$', filename)) {
        
        # then unzip it there:
        unzip(sprintf('data/%s',files$name[idx]), exdir='data/')
        
        # and remove the zip file, if that is wanted:
        if (removezip) {
          file.remove(sprintf('data/%s',files$name[idx]))
        }
        
      }
      
    }
    
  }
  
}
