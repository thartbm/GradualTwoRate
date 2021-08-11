
library('osfr')


# meta functions -----

# this provides scripts with info

getInfo <- function() {
  
  group      <- c(   'young30',    'young60',               'mcebat45',  'older45',    'young45')
  
  label      <- c('younger 30', 'younger 60', 'mild cerebellar ataxia', 'older 45', 'younger 45')
  
  rotation   <- c(          30,           60,                       45,         45,           45)
  
  experiment <- c(           1,            1,                        2,          2,            2)
  
  window     <- c(          45,           45,                       45,         45,           45)
  
  return(data.frame(group, 
                    label, 
                    rotation, 
                    experiment, 
                    window))
  
}

getStyle <- function() {
  
  condition <- c('abrupt', 'gradual')
  
  label     <- c('abrupt', 'ramped')
  
  color_s   <- c('#ff8200ff', '#b400e4ff')
  
  color_t   <- c('#ff82002f', '#b400e42f')
  
  return(data.frame(condition,
                    label,
                    color_s,
                    color_t))
  
}

# get OSF data -----

# by default, this function:
# - downloads data of all 5 groups from OSF
# - but does not overwrite it if there already is data, and
# - does not remove any downloaded zip files
# for now, the OSF repo is private, so you need a personal access token set in the OSF_PAT environment variable

downloadOSFdata <- function(groups=getInfo()$group, overwrite=FALSE, removezip=FALSE) {
  
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



# data handling -----

# this returns parsed data
# using the parse function specified
# for every group and participant present
# 
# data can be a (nested) list with data frames
# or a data frame (parseData is recursive)
#
# functions defined here:
# - block (get average reach deviations in a few specified blocks)
# - baseline (subtract average reach deviations from the second half of the aligned phase, per target)
# - addtrial (adds a column in the data frame with trial numbers)
# - addphase (adds a column in the data frame with phase numbers)

parseData <- function(data, FUN) {
  
  
  # data has to either be a list of data frames
  # or a data frame
  if (is.data.frame(data)) {
    # if it is a data frame, we process it
    return(FUN(data))
  }
  
  if (is.list(data)) {
    # if it is a list, we loop through it's entries:
    
    for (itemno in c(1:length(data))) {
      
      data[[itemno]] <- parseData(data=data[[itemno]], FUN=FUN)
      
    }
    
    return(data)
  }
  
  # this should never happen:
  print('can not parse this data, returning unchanged')
  return(data)
  
}

block <- function(df) {
  
  # we want to loop through participants:
  participants <- unique(df$participant)
  
  # we collect blocked data here:
  blocked <- NA
  
  # what are the blocks?
  ntrials <- length(which(df$participant == participants[1]))
  
  if (ntrials == 164) {
    rotated_idx <- c(123:132) # last ten?
    counter_idx <- c(141:144) # last four?
    clamped_idx <- c(155:164) # last ten?
  }
  if (ntrials == 220) {
    rotated_idx <- c(151:160) # last ten?
    counter_idx <- c(177:180) # last four?
    clamped_idx <- c(211:220) # last ten?
  }
  
  # loop through participants
  for (participant in participants) {
    
    # get the participants data
    pp_idx <- which(df$participant == participant)
    pp_df <- df[pp_idx,]
    
    # get their blocked data...
    rotated <- mean(pp_df$reachdeviation_deg[rotated_idx], na.rm=TRUE)
    counter <- mean(pp_df$reachdeviation_deg[counter_idx], na.rm=TRUE)
    clamped <- mean(pp_df$reachdeviation_deg[clamped_idx], na.rm=TRUE)
    
    pp_block <- data.frame( block=c('rotated', 'counter', 'clamped'),
                            reachdeviation_deg=c(rotated, counter, clamped))
    pp_block$participant <- participant
    
    if (is.data.frame(blocked)) {
      blocked <- rbind(blocked, pp_block)
    } else {
      blocked <- pp_block
    }
    
  }
  
  return(blocked)
  
}


baseline <- function(df) {
  
  # we want to loop through participants and targets:
  participants <- unique(df$participant)
  
  # loop through participants:
  for (participant in participants) {
    
    # get the participants' data
    pp_idx <- which(df$participant == participant) 
    pp_df <- df[pp_idx,]
    
    # targets for the participant:
    targets <- unique(pp_df$targetangle_deg)
    
    if (dim(pp_df)[1] == 164) {
      baseline_idx <- c(17:32)
    }
    if (dim(pp_df)[1] == 220) {
      baseline_idx <- c(21:40)
    }
    
    for (target in targets) {
      
      # which trials correspond to this target
      target_idx <- which(pp_df$targetangle_deg == target)
      
      # this is the average reach deviation for the target during aligned phase:
      baseline <- mean( pp_df$reachdeviation_deg[ intersect(target_idx, baseline_idx) ], na.rm=TRUE )
      
      # subtract from all reaches to the target:
      pp_df$reachdeviation_deg[target_idx] <- pp_df$reachdeviation_deg[target_idx] - baseline
      
    }
    
    # put back in main data frame
    df$reachdeviation_deg[pp_idx] <- pp_df$reachdeviation_deg
    
  }
  
  # return baselined data frame:
  return(df)
  
}

addtrial <- function(df) {
  
  # we want to loop through participants and targets:
  participants <- unique(df$participant)
  
  output <- NA
  
  # loop through participants:
  for (participant in participants) {
    
    # get the participants' data
    pp_idx <- which(df$participant == participant) 
    pp_df <- df[pp_idx,]
    
    pp_df$trial <- c(1:dim(pp_df)[1])
    
    if (is.data.frame(output)) {
      output <- rbind(output, pp_df)
    } else {
      output <- pp_df
    }
    
  }
  
  # return data frame with trial numbers:
  return(output)
  
}


addphase <- function(df) {
  
  # we want to loop through participants and targets:
  participants <- unique(df$participant)
  
  output <- NA
  
  # loop through participants:
  for (participant in participants) {
    
    # get the participants' data
    pp_idx <- which(df$participant == participant) 
    pp_df <- df[pp_idx,]
    
    if (dim(pp_df)[1] == 164) {
      pp_df$phase <- c(rep(1,32), rep(2,100), rep(3,12), rep(4,20))
    }
    if (dim(pp_df)[1] == 220) {
      pp_df$phase <- c(rep(1,40), rep(2,120), rep(3,20), rep(4,40))
    }
    
    if (is.data.frame(output)) {
      output <- rbind(output, pp_df)
    } else {
      output <- pp_df
    }
    
  }
  
  # return data frame with phase added:
  return(output)
  
}

normalize <- function(df) {
  
  rotation <- max(abs(df$rotation_deg), na.rm = TRUE)
  
  df$reachdeviation_deg <- df$reachdeviation_deg / rotation
  
  return(df)
  
}


# statistical functions -----

getConfidenceInterval <- function(data, variance = var(data, na.rm=TRUE), conf.level = 0.95, method='t-distr', resamples=1000, FUN=mean) {
  
  data <- data[which(!is.na(data))]
  
  if (method %in% c('t-distr','t')) {
    
    z = qt((1 - conf.level)/2, df = length(data) - 1, lower.tail = FALSE)
    
    xbar = mean(data)
    sdx = sqrt(variance/length(data))
    
    return(c(xbar - z * sdx, xbar + z * sdx))
    
  }
  
  # add sample z-distribution?
  
  if (method %in% c('bootstrap','b')) {
    
    data <- data[which(is.finite(data))] #need is.finite due to NA values
    
    samplematrix <- matrix(sample(data, size = resamples*length(data), replace = TRUE), nrow = resamples)
    BS <- apply(samplematrix, c(1), FUN=FUN) 
    
    lo <- (1-conf.level)/2.
    hi <- 1 - lo
    
    return(quantile(BS, probs = c(lo,hi)))
    
  }
  
}

